"""Interactive static HTML dashboard generator for Geroscience Compound Atlas + Predictive Bench."""

import html
import json
import logging
import math
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from atlas.evidence_review import apply_evidence_review

logger = logging.getLogger(__name__)


def generate_molecule_svg(smiles: str | None, width: int = 280, height: int = 180) -> str:
    """Generate inline SVG vector image for a chemical SMILES using RDKit."""
    if not smiles or not isinstance(smiles, str) or not smiles.strip():
        return (
            f'<div class="no-svg" style="width:{width}px;height:{height}px;'
            "display:flex;align-items:center;justify-content:center;background:#1e293b;color:#64748b;"
            'border-radius:8px;font-size:12px;">No Structure Available</div>'
        )

    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            return (
                f'<div class="no-svg" style="width:{width}px;height:{height}px;'
                "display:flex;align-items:center;justify-content:center;background:#1e293b;color:#ef4444;"
                'border-radius:8px;font-size:12px;">Unparseable Structure</div>'
            )

        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        opts = drawer.drawOptions()
        opts.clearBackground = False
        opts.addStereoAnnotation = True
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg_text = drawer.GetDrawingText()

        # Strip xml declaration if present for clean inline embedding
        if "<?xml" in svg_text:
            svg_text = svg_text[svg_text.find("<svg") :]
        return svg_text
    except (ValueError, RuntimeError, TypeError) as exc:
        logger.debug("SVG generation failed for %s: %s", smiles, exc)
        return (
            f'<div class="no-svg" style="width:{width}px;height:{height}px;'
            "display:flex;align-items:center;justify-content:center;background:#1e293b;color:#64748b;"
            'border-radius:8px;font-size:12px;">2D Render Error</div>'
        )


def external_validation_notice(results_path: str | Path) -> tuple[str, str]:
    """Summarize the frozen external result without treating it as a biology outcome."""
    fallback = (
        "External model qualification is unestablished",
        "The external evaluation result was missing or invalid when this dashboard was built. "
        "Model probabilities and generated scores remain exploratory outputs, not compound-prioritization evidence.",
    )
    path = Path(results_path)
    try:
        if path.is_symlink() or not path.is_file() or path.stat().st_size > 2_000_000:
            return fallback
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return fallback
    if not isinstance(data, dict) or data.get("status") != "external_target_assigned_mtor_assay_evaluation":
        return fallback
    results = data.get("results")
    if not isinstance(results, dict):
        return fallback

    comparisons = []
    for seed in ("42", "123", "456"):
        row = results.get(seed)
        baseline = row.get("training_prevalence_baseline") if isinstance(row, dict) else None
        models = row.get("models") if isinstance(row, dict) else None
        if not isinstance(baseline, dict) or not isinstance(models, dict):
            return fallback
        baseline_brier = baseline.get("brier")
        n_rows = baseline.get("n")
        logistic = models.get("logistic")
        hgb = models.get("hgb")
        if not isinstance(logistic, dict) or not isinstance(hgb, dict):
            return fallback
        values = (baseline_brier, logistic.get("brier"), hgb.get("brier"), hgb.get("auroc"))
        if (any(isinstance(value, bool) or not isinstance(value, (int, float))
                or not math.isfinite(value) for value in values)
                or isinstance(n_rows, bool) or not isinstance(n_rows, int) or n_rows < 1):
            return fallback
        comparisons.append((n_rows, baseline_brier, logistic["brier"], hgb["brier"], hgb["auroc"]))
    if len({item[0] for item in comparisons}) != 1:
        return fallback

    if all(logistic_brier > baseline and hgb_brier > baseline
           for _, baseline, logistic_brier, hgb_brier, _ in comparisons):
        hgb_aurocs = [item[4] for item in comparisons]
        return (
            "External validation: classifiers are not qualified for compound prioritization",
            f"Across {comparisons[0][0]:,} held-out molecule–assay rows, both classifiers had worse pooled "
            f"Brier error than a training-prevalence baseline in all {len(comparisons)} seeds; "
            f"HistGradientBoosting AUROC ranged from {min(hgb_aurocs):.3f} to {max(hgb_aurocs):.3f}. "
            "This is a source-document and structure holdout within ChEMBL, not independent biological or "
            "rejuvenation validation. Treat displayed probabilities and generated scores as exploratory artifacts.",
        )
    return (
        "External validation is available; model qualification remains unestablished",
        "This source-document and structure holdout within ChEMBL does not establish independent biological "
        "validation or qualify model scores for compound prioritization. Treat displayed probabilities and "
        "generated scores as exploratory artifacts.",
    )


def build_dashboard_html(
    compounds_parquet: str | Path = "artifacts/compounds.parquet",
    edges_parquet: str | Path = "artifacts/evidence_edges.parquet",
    coverage_parquet: str | Path = "artifacts/chembl_coverage.parquet",
    metrics_json: str | Path = "artifacts/metrics.json",
    gen_parquet: str | Path = "artifacts/generated_molecules.parquet",
    gen_metrics_json: str | Path = "artifacts/generator_metrics.json",
    output_html: str | Path = "artifacts/dashboard.html",
    *, site_html: str | Path | None = None,
    external_results_json: str | Path = "studies/chembl_external_2026-10-04/results.json",
) -> Path:
    """Compile processed atlas tables, benchmark metrics, and generated molecules into an interactive standalone HTML dashboard."""
    compounds_df = pd.read_parquet(compounds_parquet)
    edges_df = apply_evidence_review(pd.read_parquet(edges_parquet))
    coverage_df = pd.read_parquet(coverage_parquet)

    with open(metrics_json, "r", encoding="utf-8") as f:
        metrics_data = json.load(f)
    external_notice_title, external_notice_text = external_validation_notice(external_results_json)

    # Join coverage information with compounds
    merged_compounds = compounds_df.merge(
        coverage_df[
            ["inchikey", "has_chembl", "has_binding", "num_activities", "num_binding_assays"]
        ],
        on="inchikey",
        how="left",
    )

    # Map edges to compounds
    edges_by_inchikey: dict[str, list[dict[str, Any]]] = {}
    for _, edge_row in edges_df.iterrows():
        key = str(edge_row.get("inchikey", ""))
        edges_by_inchikey.setdefault(key, []).append(edge_row.to_dict())

    # Build Compound Cards JSON
    cards_data = []
    for _, comp in merged_compounds.iterrows():
        inchikey = str(comp.get("inchikey", ""))
        smiles = comp.get("canonical_smiles")
        comp_edges = edges_by_inchikey.get(inchikey, [])

        svg_content = (
            generate_molecule_svg(smiles) if pd.notna(smiles) else generate_molecule_svg(None)
        )

        cards_data.append(
            {
                "raw_name": comp.get("raw_name"),
                "resolved_name": comp.get("resolved_name"),
                "cid": int(comp["cid"]) if pd.notna(comp.get("cid")) else None,
                "inchikey": inchikey if pd.notna(inchikey) else "N/A",
                "canonical_smiles": smiles if pd.notna(smiles) else "N/A",
                "modality": comp.get("modality", "small_molecule"),
                "resolved": bool(comp.get("resolved", False)),
                "has_chembl": bool(comp.get("has_chembl", False)),
                "has_binding": bool(comp.get("has_binding", False)),
                "num_activities": int(comp.get("num_activities", 0))
                if pd.notna(comp.get("num_activities"))
                else 0,
                "num_binding_assays": int(comp.get("num_binding_assays", 0))
                if pd.notna(comp.get("num_binding_assays"))
                else 0,
                "mol_wt": float(comp["mol_wt"]) if pd.notna(comp.get("mol_wt")) else None,
                "log_p": float(comp["log_p"]) if pd.notna(comp.get("log_p")) else None,
                "tpsa": float(comp["tpsa"]) if pd.notna(comp.get("tpsa")) else None,
                "num_h_donors": int(comp["num_h_donors"])
                if pd.notna(comp.get("num_h_donors"))
                else None,
                "num_h_acceptors": int(comp["num_h_acceptors"])
                if pd.notna(comp.get("num_h_acceptors"))
                else None,
                "qed": float(comp["qed"]) if pd.notna(comp.get("qed")) else None,
                "notes": comp.get("notes", ""),
                "svg": svg_content,
                "edges": comp_edges,
            }
        )

    # Load Phase 5 generated candidates if available
    gen_cards_data = []
    gen_file = Path(gen_parquet)
    if gen_file.exists():
        try:
            gen_df = pd.read_parquet(gen_file)
            for _, grow in gen_df.head(12).iterrows():
                smi = str(grow["smiles"])
                svg = generate_molecule_svg(smi, width=280, height=180)
                gen_cards_data.append(
                    {
                        "smiles": smi,
                        "inchikey": str(grow["inchikey"]),
                        "generation": int(grow.get("generation", 0)),
                        "reward": round(float(grow.get("reward", 0.0)), 3),
                        "mtor_prob": round(float(grow.get("mtor_prob", 0.0)), 3),
                        "qed": round(float(grow.get("qed", 0.0)), 3),
                        "has_pains": bool(grow.get("has_pains", False)),
                        "svg": svg,
                    }
                )
        except (ValueError, OSError, KeyError) as exc:
            logger.warning("Failed to load generated molecules for dashboard: %s", exc)

    # Convert data to JSON for embedding
    cards_json = json.dumps(cards_data)
    top_fps_json = json.dumps(metrics_data.get("top_false_positives", []))
    gen_cards_json = json.dumps(gen_cards_data)

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Geroscience Compound Atlas + Predictive Bench</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&family=JetBrains+Mono:wght@400;500;600&display=swap" rel="stylesheet">
    <style>
        :root {{
            --bg-primary: #0b0f19;
            --bg-card: #151d30;
            --bg-card-hover: #1c263f;
            --border-color: #24304f;
            --text-primary: #f8fafc;
            --text-secondary: #94a3b8;
            --text-muted: #64748b;
            --accent-blue: #3b82f6;
            --accent-blue-hover: #2563eb;
            --accent-green: #10b981;
            --accent-amber: #f59e0b;
            --accent-red: #ef4444;
            --accent-purple: #8b5cf6;
        }}

        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }}

        body {{
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
            background-color: var(--bg-primary);
            color: var(--text-primary);
            line-height: 1.5;
            padding: 0;
        }}

        header {{
            background: linear-gradient(180deg, #111827 0%, #0b0f19 100%);
            border-bottom: 1px solid var(--border-color);
            padding: 2.5rem 2rem 2rem;
        }}

        .container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 0 1.5rem;
        }}

        .header-content {{
            display: flex;
            justify-content: space-between;
            align-items: flex-start;
            flex-wrap: wrap;
            gap: 1.5rem;
        }}

        .title-area h1 {{
            font-size: 2rem;
            font-weight: 800;
            letter-spacing: -0.025em;
            background: linear-gradient(135deg, #60a5fa 0%, #a78bfa 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            margin-bottom: 0.5rem;
        }}

        .title-area p.pitch {{
            color: var(--text-secondary);
            font-size: 1.05rem;
            max-width: 800px;
        }}

        .safety-banner {{
            display: inline-flex;
            align-items: center;
            gap: 0.5rem;
            background: rgba(239, 68, 68, 0.1);
            border: 1px solid rgba(239, 68, 68, 0.3);
            color: #fca5a5;
            padding: 0.5rem 1rem;
            border-radius: 9999px;
            font-size: 0.825rem;
            font-weight: 600;
        }}

        /* Metrics Bar */
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 1rem;
            margin: 2rem 0;
        }}

        .metric-card {{
            background: var(--bg-card);
            border: 1px solid var(--border-color);
            border-radius: 12px;
            padding: 1.25rem;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.2);
            transition: transform 0.2s, border-color 0.2s;
        }}

        .metric-card:hover {{
            transform: translateY(-2px);
            border-color: #3b82f6;
        }}

        .metric-label {{
            font-size: 0.8rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: var(--text-muted);
            font-weight: 600;
            margin-bottom: 0.25rem;
        }}

        .metric-value {{
            font-size: 1.8rem;
            font-weight: 800;
            color: var(--text-primary);
        }}

        .metric-sub {{
            font-size: 0.775rem;
            color: var(--text-secondary);
            margin-top: 0.25rem;
        }}

        /* Navigation Tabs */
        .tabs {{
            display: flex;
            gap: 1rem;
            border-bottom: 1px solid var(--border-color);
            margin-bottom: 2rem;
        }}

        .tab-btn {{
            background: none;
            border: none;
            color: var(--text-secondary);
            font-size: 1rem;
            font-weight: 600;
            padding: 0.75rem 0.5rem;
            cursor: pointer;
            border-bottom: 2px solid transparent;
            transition: all 0.2s;
        }}

        .tab-btn:hover {{
            color: var(--text-primary);
        }}

        .tab-btn.active {{
            color: var(--accent-blue);
            border-bottom-color: var(--accent-blue);
        }}

        /* Filters */
        .filter-panel {{
            background: var(--bg-card);
            border: 1px solid var(--border-color);
            border-radius: 12px;
            padding: 1.25rem;
            margin-bottom: 2rem;
            display: flex;
            flex-wrap: wrap;
            gap: 1rem;
            align-items: center;
        }}

        .search-box {{
            flex: 1;
            min-width: 260px;
            position: relative;
        }}

        .search-box input {{
            width: 100%;
            background: #0b0f19;
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 0.6rem 1rem;
            color: var(--text-primary);
            font-size: 0.9rem;
            outline: none;
            transition: border-color 0.2s;
        }}

        .search-box input:focus {{
            border-color: var(--accent-blue);
        }}

        .filter-select {{
            background: #0b0f19;
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 0.6rem 1rem;
            color: var(--text-primary);
            font-size: 0.9rem;
            outline: none;
            cursor: pointer;
        }}

        /* Compound Card Grid */
        .cards-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(360px, 1fr));
            gap: 1.5rem;
            margin-bottom: 3rem;
        }}

        .compound-card {{
            background: var(--bg-card);
            border: 1px solid var(--border-color);
            border-radius: 14px;
            padding: 1.5rem;
            display: flex;
            flex-direction: column;
            transition: all 0.2s;
            position: relative;
            overflow: hidden;
        }}

        .compound-card:hover {{
            border-color: #3b82f6;
            background: var(--bg-card-hover);
            box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.3);
        }}

        .card-header {{
            display: flex;
            justify-content: space-between;
            align-items: flex-start;
            margin-bottom: 1rem;
        }}

        .card-title {{
            font-size: 1.25rem;
            font-weight: 700;
            color: var(--text-primary);
        }}

        .card-subtitle {{
            font-size: 0.8rem;
            color: var(--text-muted);
            margin-top: 0.1rem;
        }}

        .modality-pill {{
            font-size: 0.7rem;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            padding: 0.25rem 0.6rem;
            border-radius: 9999px;
            background: rgba(59, 130, 246, 0.15);
            color: #93c5fd;
            border: 1px solid rgba(59, 130, 246, 0.3);
        }}

        .modality-pill.peptide {{
            background: rgba(16, 185, 129, 0.15);
            color: #6ee7b7;
            border-color: rgba(16, 185, 129, 0.3);
        }}

        .modality-pill.other {{
            background: rgba(239, 68, 68, 0.15);
            color: #fca5a5;
            border-color: rgba(239, 68, 68, 0.3);
        }}

        .svg-container {{
            background: #0d1322;
            border: 1px solid #1e293b;
            border-radius: 10px;
            height: 180px;
            display: flex;
            align-items: center;
            justify-content: center;
            margin-bottom: 1.2rem;
            overflow: hidden;
        }}

        .svg-container svg {{
            max-width: 100%;
            max-height: 100%;
        }}

        .identifiers-box {{
            background: #0f172a;
            border: 1px solid #1e293b;
            border-radius: 8px;
            padding: 0.75rem;
            margin-bottom: 1rem;
            font-family: 'JetBrains Mono', monospace;
            font-size: 0.75rem;
        }}

        .id-row {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 0.35rem;
        }}

        .id-row:last-child {{
            margin-bottom: 0;
        }}

        .id-label {{
            color: var(--text-muted);
            font-weight: 600;
        }}

        .id-val {{
            color: #cbd5e1;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
            max-width: 200px;
        }}

        .copy-btn {{
            background: none;
            border: none;
            color: var(--accent-blue);
            cursor: pointer;
            font-size: 0.75rem;
            margin-left: 0.4rem;
        }}

        .copy-btn:hover {{
            text-decoration: underline;
        }}

        /* Property Badges */
        .properties-grid {{
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            gap: 0.5rem;
            margin-bottom: 1rem;
            font-size: 0.75rem;
        }}

        .prop-item {{
            background: #0e1526;
            padding: 0.4rem 0.5rem;
            border-radius: 6px;
            border: 1px solid #1b253b;
        }}

        .prop-key {{
            color: var(--text-muted);
            display: block;
            font-size: 0.675rem;
            font-weight: 600;
        }}

        .prop-val {{
            font-weight: 700;
            color: #f1f5f9;
        }}

        /* Evidence Section */
        .evidence-section {{
            margin-top: auto;
            border-top: 1px solid var(--border-color);
            padding-top: 1rem;
        }}

        .evidence-header {{
            font-size: 0.8rem;
            font-weight: 700;
            color: var(--text-secondary);
            text-transform: uppercase;
            letter-spacing: 0.05em;
            margin-bottom: 0.5rem;
            display: flex;
            justify-content: space-between;
        }}

        .edge-pill {{
            background: #0f172a;
            border: 1px solid #1e293b;
            border-radius: 8px;
            padding: 0.5rem 0.75rem;
            margin-bottom: 0.5rem;
            font-size: 0.8rem;
        }}

        .edge-pill:last-child {{
            margin-bottom: 0;
        }}

        .edge-top {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 0.25rem;
        }}

        .hallmark-tag {{
            font-weight: 600;
            color: #93c5fd;
        }}

        .grade-badge {{
            font-size: 0.7rem;
            font-weight: 800;
            padding: 0.15rem 0.5rem;
            border-radius: 4px;
            letter-spacing: 0.05em;
        }}

        .grade-E4 {{ background: #059669; color: #ecfdf5; }}
        .grade-E3 {{ background: #d97706; color: #fffbeb; }}
        .grade-E2 {{ background: #2563eb; color: #eff6ff; }}
        .grade-E1 {{ background: #475569; color: #f8fafc; }}
        .grade-E0 {{ background: #334155; color: #94a3b8; border: 1px dashed #64748b; }}

        .review-status {{
            margin-top: 0.25rem;
            color: #fbbf24;
            font-size: 0.68rem;
            line-height: 1.35;
        }}

        .edge-notes {{
            font-size: 0.75rem;
            color: var(--text-secondary);
            line-height: 1.35;
        }}

        /* Tables */
        .table-panel {{
            background: var(--bg-card);
            border: 1px solid var(--border-color);
            border-radius: 12px;
            padding: 1.5rem;
            margin-bottom: 3rem;
            overflow-x: auto;
        }}

        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 0.9rem;
            text-align: left;
        }}

        th {{
            background: #0f172a;
            color: var(--text-secondary);
            padding: 0.85rem 1rem;
            font-size: 0.775rem;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            border-bottom: 2px solid var(--border-color);
        }}

        td {{
            padding: 0.85rem 1rem;
            border-bottom: 1px solid var(--border-color);
            color: #cbd5e1;
        }}

        tr:hover td {{
            background: var(--bg-card-hover);
        }}

        .highlight-row {{
            background: rgba(59, 130, 246, 0.08);
            font-weight: 600;
        }}

        /* Coverage Holes Box */
        .holes-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
            gap: 1rem;
            margin-top: 1rem;
        }}

        .hole-card {{
            background: #111827;
            border: 1px solid #374151;
            border-left: 4px solid var(--accent-red);
            border-radius: 8px;
            padding: 1rem;
        }}

        .hole-title {{
            font-weight: 700;
            color: #f87171;
            font-size: 1rem;
            margin-bottom: 0.35rem;
        }}

        .hole-desc {{
            font-size: 0.825rem;
            color: var(--text-secondary);
        }}

        /* Hidden Utility */
        .hidden {{
            display: none !important;
        }}

        footer {{
            border-top: 1px solid var(--border-color);
            padding: 2.5rem 1.5rem;
            text-align: center;
            color: var(--text-muted);
            font-size: 0.85rem;
        }}
    </style>
</head>
<body>

<header>
    <div class="container header-content">
        <div class="title-area">
            <h1>Geroscience Compound Atlas + Predictive Bench</h1>
            <p class="pitch">
                A computational atlas of public small molecules, peptides and recorded aging-mechanism claims, alongside a frozen Bemis-Murcko scaffold predictive benchmark.
            </p>
            <p class="pitch" style="margin-top: 0.75rem; font-size: 0.9rem;">
                Personal hobby and learning project, developed with substantial assistance from AI coding tools.
            </p>
            <p style="margin-top:0.75rem;"><a href="research.html" style="color:#93c5fd;">Explore research decisions: oral delivery, organ targeting and stimulant reversal</a></p>
        </div>
        <div>
            <span class="safety-banner">
                <svg width="16" height="16" fill="currentColor" viewBox="0 0 20 20"><path fill-rule="evenodd" d="M8.257 3.099c.765-1.36 2.722-1.36 3.486 0l5.58 9.92c.75 1.334-.213 2.98-1.742 2.98H4.42c-1.53 0-2.493-1.646-1.743-2.98l5.58-9.92zM11 13a1 1 0 11-2 0 1 1 0 012 0zm-1-8a1 1 0 00-1 1v3a1 1 0 002 0V6a1 1 0 00-1-1z" clip-rule="evenodd"></path></svg>
                Research & Educational Artifacts Only — No Medical or Protocol Advice
            </span>
        </div>
    </div>
</header>

<section class="container" role="note" aria-label="External model evaluation" style="margin-top:1.25rem;padding:1rem 1.25rem;background:#251b23;border:1px solid #7f1d3c;border-left:5px solid #fb7185;border-radius:10px;">
    <h2 style="margin:0 0 0.5rem;color:#fda4af;font-size:1rem;">{html.escape(external_notice_title)}</h2>
    <p style="margin:0;color:#e2e8f0;line-height:1.6;">{html.escape(external_notice_text)}</p>
    <a href="https://github.com/dylanstechmann/geroscience-compound-atlas/blob/main/studies/chembl_external_2026-10-04/README.md" style="display:inline-block;margin-top:0.5rem;color:#93c5fd;">Read the evaluation, protocol, and limitations</a>
</section>

<main class="container">
    <section role="note" aria-label="Citation review" style="margin-top:1.25rem;padding:1rem 1.25rem;background:#1e293b;border:1px solid #475569;border-radius:10px;">
        <h2 style="margin:0 0 0.5rem;color:#fbbf24;font-size:1rem;">Claim support remains unreviewed</h2>
        <p style="margin:0;color:#e2e8f0;line-height:1.6;">Citation screens check identifier titles, not whether a paper supports a claim. Current displayed grades are E0 pending source-to-claim review. Historical catalog grades and review labels are retained separately; the E4 filter currently has no eligible claims.</p>
        <a href="https://github.com/dylanstechmann/geroscience-compound-atlas/blob/main/docs/citation-audit-2026-10-07/README.md" style="display:inline-block;margin-top:0.5rem;color:#93c5fd;">Read the dated citation audit and its limits</a>
    </section>
    <!-- Top Summary Metrics -->
    <div class="metrics-grid">
        <div class="metric-card">
            <div class="metric-label">Watchlist Ingested</div>
            <div class="metric-value">20</div>
            <div class="metric-sub">19 Resolved (95%) | 1 Negative Control</div>
        </div>
        <div class="metric-card">
            <div class="metric-label">ChEMBL Binding Coverage</div>
            <div class="metric-value">75.0%</div>
            <div class="metric-sub">15/20 with ≥1 structured binding assay</div>
        </div>
        <div class="metric-card">
            <div class="metric-label">Biological Targets</div>
            <div class="metric-value">209</div>
            <div class="metric-sub">662 bioactivity records mapped</div>
        </div>
        <div class="metric-card">
            <div class="metric-label">Recorded Claim Edges</div>
            <div class="metric-value">{len(edges_df)}</div>
            <div class="metric-sub">Claim support unreviewed; displayed as E0</div>
        </div>
        <div class="metric-card">
            <div class="metric-label">Internal Scaffold AUROC</div>
            <div class="metric-value" style="color: #60a5fa;">0.9735</div>
            <div class="metric-sub">Original frozen split; external holdout underperformed baseline</div>
        </div>
    </div>

    <!-- Navigation Tabs -->
    <div class="tabs">
        <button class="tab-btn active" onclick="switchTab('cards')">Compound Cards (<span id="visible-count">20</span>)</button>
        <button class="tab-btn" onclick="switchTab('holes')">Coverage Holes (3)</button>
        <button class="tab-btn" onclick="switchTab('bench')">Benchmark Bake-Off & Leakage</button>
        <button class="tab-btn" onclick="switchTab('errors')">Top False Positives (10)</button>
        <button class="tab-btn" onclick="switchTab('generator')">Generated Candidates (200)</button>
    </div>

    <!-- TAB 1: COMPOUND CARDS -->
    <div id="tab-cards">
        <!-- Filter Controls -->
        <div class="filter-panel">
            <div class="search-box">
                <input type="text" id="search-input" placeholder="Search by name, InChIKey, SMILES, target, or hallmark..." onkeyup="filterCards()">
            </div>
            <select class="filter-select" id="hallmark-filter" onchange="filterCards()">
                <option value="">All Hallmarks</option>
                <option value="nutrient_sensing">Nutrient Sensing</option>
                <option value="cellular_senescence">Cellular Senescence</option>
                <option value="mitochondrial_dysfunction">Mitochondrial Dysfunction</option>
                <option value="autophagy">Autophagy</option>
                <option value="chronic_inflammation">Chronic Inflammation</option>
                <option value="stem_cell_exhaustion">Stem Cell Exhaustion</option>
                <option value="proteostasis">Proteostasis</option>
                <option value="intercellular_communication">Intercellular Communication</option>
                <option value="ecm_fibrosis">ECM / Fibrosis</option>
                <option value="telomere_attrition">Telomere Attrition</option>
                <option value="dysbiosis">Dysbiosis</option>
                <option value="epigenetic_alteration">Epigenetic Alteration</option>
            </select>
            <select class="filter-select" id="grade-filter" onchange="filterCards()">
                <option value="">All Evidence Grades</option>
                <option value="E4">E4 (Mammalian Lifespan / Human)</option>
                <option value="E3">E3+ (In Vivo Model Organism)</option>
                <option value="E2">E2+ (Target Engagement)</option>
                <option value="E1">E1+ (In Vitro Assay)</option>
            </select>
            <select class="filter-select" id="modality-filter" onchange="filterCards()">
                <option value="">All Modalities</option>
                <option value="small_molecule">Small Molecule</option>
                <option value="peptide">Peptide</option>
                <option value="other">Other</option>
            </select>
            <select class="filter-select" id="chembl-filter" onchange="filterCards()">
                <option value="">All ChEMBL Data</option>
                <option value="binding">Has Binding Assays</option>
                <option value="has_chembl">Any ChEMBL Activity</option>
                <option value="none">No ChEMBL Activity (Holes)</option>
            </select>
            <button class="filter-preset-btn" onclick="resetFilters()" style="padding:0.4rem 0.8rem;background:#334155;color:#f8fafc;border:1px solid #475569;border-radius:6px;cursor:pointer;font-size:0.8rem;">Reset Filters</button>
        </div>
        <div style="display:flex;gap:0.5rem;align-items:center;flex-wrap:wrap;margin-bottom:1rem;font-size:0.8rem;color:#94a3b8;">
            <span>Quick Filters:</span>
            <button onclick="applyPreset('mtor')" style="padding:0.25rem 0.5rem;background:#1e293b;border:1px solid #3b82f6;color:#93c5fd;border-radius:4px;cursor:pointer;font-size:0.75rem;">mTOR Pathway</button>
            <button onclick="applyPreset('senolyt')" style="padding:0.25rem 0.5rem;background:#1e293b;border:1px solid #3b82f6;color:#93c5fd;border-radius:4px;cursor:pointer;font-size:0.75rem;">Senolytics</button>
            <button onclick="applyPreset('ampk')" style="padding:0.25rem 0.5rem;background:#1e293b;border:1px solid #3b82f6;color:#93c5fd;border-radius:4px;cursor:pointer;font-size:0.75rem;">AMPK / Metformin</button>
            <button onclick="applyPreset('peptide')" style="padding:0.25rem 0.5rem;background:#1e293b;border:1px solid #3b82f6;color:#93c5fd;border-radius:4px;cursor:pointer;font-size:0.75rem;">Peptides Only</button>
                <button onclick="applyPreset('E4')" style="padding:0.25rem 0.5rem;background:#1e293b;border:1px solid #3b82f6;color:#93c5fd;border-radius:4px;cursor:pointer;font-size:0.75rem;">E4 Reported Lifespan / Human Outcome</button>
            <span id="filter-status-text" style="font-size:0.8rem;color:#60a5fa;margin-left:auto;">Showing <span id="visible-count">20</span> of 20 compounds</span>
        </div>

        <div class="cards-grid" id="cards-container">
            <!-- Rendered by JS -->
        </div>
    </div>

    <!-- TAB 2: COVERAGE HOLES -->
    <div id="tab-holes" class="hidden">
        <div class="table-panel">
            <h2 style="font-size: 1.3rem; margin-bottom: 0.5rem;">Coverage Holes in Historical Chemical Repositories</h2>
            <p style="color: var(--text-secondary); font-size: 0.95rem; margin-bottom: 1.5rem;">
                Scientific humility requires exposing data absence. ChEMBL is historically optimized for approved oral small molecules and kinase/GPCR campaigns. Peptides, natural extracts, and gray-market compounds represent systematic coverage holes.
            </p>
            <div class="holes-grid">
                <div class="hole-card">
                    <div class="hole-title">BPC-157 (Synthetic Pentadecapeptide)</div>
                    <div class="hole-desc">
                        <strong>Reason for Hole:</strong> Modality exclusion. Pentadecapeptides are excluded from standard small-molecule high-throughput binding screens. Zero structured binding assays exist in ChEMBL.<br>
                        <strong>Scientific Status:</strong> Preclinical tissue-healing literature in rodents (E1); marketed heavily in gray-market wellness communities without mammalian lifespan proof (E0).
                    </div>
                </div>
                <div class="hole-card">
                    <div class="hole-title">Epitalon (Synthetic Tetrapeptide)</div>
                    <div class="hole-desc">
                        <strong>Reason for Hole:</strong> Modality & provenance. Ala-Glu-Asp-Gly peptide with Russian literature citations; no curated target engagement entries in ChEMBL.<br>
                        <strong>Scientific Status:</strong> Telomerase elongation claims in cultured human fibroblasts (E1); widely sold as an anti-aging elixir without replicated ITP lifespan validation (E0).
                    </div>
                </div>
                <div class="hole-card">
                    <div class="hole-title">FakeGeroCompoundX99 (Negative Control)</div>
                    <div class="hole-desc">
                        <strong>Reason for Hole:</strong> Unresolvable test entity. Intentionally inserted into the seed watchlist to verify that our ingestion pipeline detects and reports resolution failures rather than silently hallucinating identifiers.<br>
                        <strong>Scientific Status:</strong> Unresolved (404); correctly flagged as unmapped.
                    </div>
                </div>
            </div>
        </div>
    </div>

    <!-- TAB 3: BENCHMARK BAKE-OFF -->
    <div id="tab-bench" class="hidden">
        <div class="table-panel">
            <h2 style="font-size: 1.3rem; margin-bottom: 0.5rem;">Benchmark Model Comparison — Scaffold Split vs Random Split</h2>
            <p style="color: var(--text-secondary); font-size: 0.95rem; margin-bottom: 1.5rem;">
                Evaluating L2-regularized Logistic Regression (Baseline) vs HistGradientBoosting (Contender) across 3 seeds (42, 123, 456) on 2048-bit Morgan circular fingerprints + standardized RDKit 2D descriptors on ChEMBL mTOR kinase activity (N = 561 molecules).
            </p>
            <table>
                <thead>
                    <tr>
                        <th>Split Strategy</th>
                        <th>Model Pipeline</th>
                        <th>AUROC</th>
                        <th>AUPRC</th>
                        <th>Recall @ 5% FPR</th>
                        <th>Brier Score</th>
                    </tr>
                </thead>
                <tbody>
                    <tr class="highlight-row">
                        <td><strong>Bemis-Murcko Scaffold (Honest)</strong></td>
                        <td><strong>Baseline (Logistic Regression)</strong></td>
                        <td><strong>0.9735 ± 0.0203</strong></td>
                        <td><strong>0.9929 ± 0.0057</strong></td>
                        <td><strong>0.8496 ± 0.1207</strong></td>
                        <td><strong>0.0486 ± 0.0238</strong></td>
                    </tr>
                    <tr>
                        <td>Bemis-Murcko Scaffold (Honest)</td>
                        <td>Contender (HistGradientBoosting)</td>
                        <td>0.9434 ± 0.0294</td>
                        <td>0.9855 ± 0.0095</td>
                        <td>0.7187 ± 0.1550</td>
                        <td>0.0820 ± 0.0392</td>
                    </tr>
                    <tr style="color: #94a3b8; font-style: italic;">
                        <td>Stratified Random (Diagnostic)</td>
                        <td>Baseline (Logistic Regression)</td>
                        <td>0.8912 ± 0.0482</td>
                        <td>0.9780 ± 0.0106</td>
                        <td>0.6528 ± 0.0856</td>
                        <td>0.0922 ± 0.0223</td>
                    </tr>
                    <tr style="color: #94a3b8; font-style: italic;">
                        <td>Stratified Random (Diagnostic)</td>
                        <td>Contender (HistGradientBoosting)</td>
                        <td>0.8766 ± 0.0467</td>
                        <td>0.9748 ± 0.0117</td>
                        <td>0.6667 ± 0.1623</td>
                        <td>0.1159 ± 0.0147</td>
                    </tr>
                </tbody>
            </table>
            <div style="margin-top: 1.5rem; padding: 1rem; background: #0f172a; border-radius: 8px; border: 1px solid #1e293b;">
                <h4 style="color: #93c5fd; margin-bottom: 0.5rem;">Benchmark observations:</h4>
                <ul style="padding-left: 1.25rem; font-size: 0.875rem; color: #cbd5e1; line-height: 1.6;">
                    <li><strong>Linear Baseline Generalization:</strong> Regularized linear models (Logistic Regression) outperform non-linear tree ensembles on sparse circular fingerprint representations ($p = 2057$), avoiding overfitting to dominant scaffold clusters.</li>
                    <li><strong>Zero Scaffold Leakage:</strong> Enforcing Bemis-Murcko clustering ensures test compounds share zero scaffold topology with training compounds (Train &cap; Test = &empty;).</li>
                </ul>
            </div>
        </div>
    </div>

    <!-- TAB 4: FALSE POSITIVES -->
    <div id="tab-errors" class="hidden">
        <div class="table-panel">
            <h2 style="font-size: 1.3rem; margin-bottom: 0.5rem;">Top Highest-Confidence False Positive Predictions</h2>
            <p style="color: var(--text-secondary); font-size: 0.95rem; margin-bottom: 1.5rem;">
                Compounds with inactive ground truth (<span style="font-family: monospace;">pChEMBL &lt; 6.0</span>) assigned the highest internal model scores under the frozen scaffold split; these scores are not externally qualified.
            </p>
            <table id="fps-table">
                <thead>
                    <tr>
                        <th>ChEMBL ID</th>
                        <th>Predicted Prob</th>
                        <th>True pChEMBL</th>
                        <th>Murcko Scaffold</th>
                        <th>Chemical & Pharmacophore Analysis</th>
                    </tr>
                </thead>
                <tbody id="fps-tbody">
                    <!-- Populated by JS -->
                </tbody>
            </table>
        </div>
    <!-- TAB 5: GENERATED CANDIDATES (PHASE 5) -->
    <div id="tab-generator" class="hidden">
        <div class="table-panel">
            <h2 style="font-size: 1.3rem; margin-bottom: 0.5rem;">Constrained Molecular Generator (Genetic Algorithm)</h2>
            <p style="color: var(--text-secondary); font-size: 0.95rem; margin-bottom: 1.5rem;">
                Exploring the behavior of the frozen ChEMBL mTOR surrogate under QED-bias and PAINS constraints. External source-document and structure validation underperformed a prevalence baseline; displayed probabilities and generated scores are unqualified software artifacts.
            </p>
            
            <div class="metrics-grid" style="margin-bottom: 2rem;">
                <div class="metric-card">
                    <div class="metric-label">Generated Structures</div>
                    <div class="metric-value">200</div>
                    <div class="metric-sub">100% Valid &bull; 100% Unique</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Novelty vs Training</div>
                    <div class="metric-value" style="color: #34d399;">99.5%</div>
                    <div class="metric-sub">199 / 200 unobserved InChIKeys</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Internal Diversity</div>
                    <div class="metric-value">0.655</div>
                    <div class="metric-sub">1 - mean pairwise Tanimoto</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Mean surrogate probability</div>
                    <div class="metric-value" style="color: #60a5fa;">0.991</div>
                    <div class="metric-sub">Surrogate active probability</div>
                </div>
            </div>

            <h3 style="font-size: 1.1rem; margin-bottom: 0.75rem; color: #93c5fd;">QED Bias Sensitivity Analysis (&lambda;<sub>QED</sub>)</h3>
            <p style="color: var(--text-secondary); font-size: 0.875rem; margin-bottom: 1rem;">
                Evaluating how the molecular generator behaves when historical small-molecule oral drug-likeness priors are removed vs enforced.
            </p>
            <table>
                <thead>
                    <tr>
                        <th>QED Weight (&lambda;<sub>QED</sub>)</th>
                        <th>Objective Character</th>
                        <th>Mean surrogate probability</th>
                        <th>Mean QED</th>
                        <th>Internal Diversity</th>
                        <th>Novelty Rate</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td>&lambda; = 0.0</td>
                        <td>Pure Target Affinity (Unconstrained)</td>
                        <td>1.000</td>
                        <td>0.059</td>
                        <td>0.563</td>
                        <td>100.0%</td>
                    </tr>
                    <tr class="highlight-row">
                        <td><strong>&lambda; = 0.2</strong></td>
                        <td><strong>Balanced Tradeoff (Primary)</strong></td>
                        <td><strong>0.991</strong></td>
                        <td><strong>0.754</strong></td>
                        <td><strong>0.655</strong></td>
                        <td><strong>99.5%</strong></td>
                    </tr>
                    <tr>
                        <td>&lambda; = 0.5</td>
                        <td>Oral Drug-Likeness Constrained</td>
                        <td>0.985</td>
                        <td>0.635</td>
                        <td>0.670</td>
                        <td>96.5%</td>
                    </tr>
                </tbody>
            </table>

            <div style="margin-top: 1.5rem; padding: 1rem; background: #0f172a; border-radius: 8px; border: 1px solid #1e293b; margin-bottom: 2rem;">
                <h4 style="color: #93c5fd; margin-bottom: 0.5rem;">Chemical Space Takeaway (PCA Embedding):</h4>
                <p style="font-size: 0.875rem; color: #cbd5e1; line-height: 1.6;">
                    Principal Component projection (PC1: 13.3%, PC2: 9.1% variance) demonstrates that GA-evolved candidates explore the mTOR training active chemotype manifold adjacent to benchmark binders, rather than drifting into unphysical junk space. When unconstrained (&lambda;<sub>QED</sub> = 0.0), molecules drift toward high molecular weight macrocyclic configurations (mean QED 0.059), demonstrating the restrictive nature of historical small-molecule heuristics.
                </p>
            </div>

            <div style="display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:1rem;margin-bottom:1rem;">
                <h3 style="font-size: 1.1rem; color: #93c5fd; margin: 0;">Top Generated Candidate Molecules (Inline 2D Vectors)</h3>
                <span id="gen-count-text" style="font-size:0.8rem;color:#60a5fa;">Showing <span id="gen-visible-count">200</span> candidates</span>
            </div>
            <div class="filter-panel" style="margin-bottom: 1.5rem;">
                <div class="search-box">
                    <input type="text" id="gen-search-input" placeholder="Search generated candidates by InChIKey or generation..." onkeyup="filterGenCards()">
                </div>
                <select class="filter-select" id="gen-pains-filter" onchange="filterGenCards()">
                    <option value="">All PAINS Status</option>
                    <option value="pass">Pass Only (Zero Motifs)</option>
                    <option value="flagged">Flagged Only</option>
                </select>
                <select class="filter-select" id="gen-sort-select" onchange="filterGenCards()">
                    <option value="reward">Sort: Composite Reward</option>
                    <option value="mtor">Sort: mTOR surrogate score</option>
                    <option value="qed">Sort: QED Drug-Likeness</option>
                </select>
            </div>
            <div class="cards-grid" id="gen-cards-container">
                <!-- Rendered by JS -->
            </div>
        </div>
    </div>
</main>

<footer>
    <div class="container">
        <p><strong>Geroscience Compound Atlas + Predictive Bench</strong> | Dylan Stechmann</p>
        <p style="margin-top: 0.5rem; font-size: 0.8rem;">Research & educational repository. All outputs are computational artifacts. Refusing all dosing, clinical protocol, or rejuvenation claims.</p>
    </div>
</footer>

<script>
    const compoundsData = {cards_json};
    const topFPsData = {top_fps_json};
    const genCardsData = {gen_cards_json};

    function escapeHtml(value) {{
        return String(value == null ? '' : value).replace(/[&<>"']/g, char => ({{
            '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'
        }})[char]);
    }}

    function safeHttpUrl(value) {{
        try {{
            const url = new URL(String(value));
            return ['http:', 'https:'].includes(url.protocol) ? url.href : '';
        }} catch (_error) {{
            return '';
        }}
    }}

    function edgeReviewLabel(edge) {{
        if (edge.claim_support_status === 'unreviewed') {{
            return 'Claim support unreviewed; curator grade ' + (edge.curator_grade || 'unknown') +
                ' withheld, shown as E0. Citation screen: ' + edge.citation_audit_category +
                ' (' + (edge.citation_audit_date || edge.citation_audit_status) + '), titles only.';
        }}
        if (edge.source_review_status === 'source_reviewed_with_claim_limits') {{
            return 'Source reviewed; claim limits recorded.';
        }}
        if (edge.source_review_status === 'claim_support_unreviewed') {{
            return 'Claim support unreviewed; curator grade ' + (edge.curator_grade || 'unknown') + ' withheld, shown as E0.';
        }}
        if (edge.source_review_status === 'unverified_vendor_or_gray_market_claim') {{
            return 'Unverified vendor/gray-market claim; shown as E0.';
        }}
        return 'Source status unverified; shown as E0.';
    }}

    function renderCards(data) {{
        const container = document.getElementById('cards-container');
        container.innerHTML = '';

        document.getElementById('visible-count').innerText = data.length;

        if (data.length === 0) {{
            container.innerHTML = '<div style="grid-column: 1/-1; text-align: center; padding: 4rem; color: #64748b;">No compounds match your filter criteria.</div>';
            return;
        }}

        data.forEach(comp => {{
            const card = document.createElement('div');
            card.className = 'compound-card';

            const edgesHtml = comp.edges.map(e => {{
                const grade = /^E[0-4]$/.test(String(e.grade || '')) ? String(e.grade) : 'E0';
                const sourceUrl = safeHttpUrl(e.source_url);
                const sourceLabel = e.source_title || e.source_url || '';
                const sourceHtml = sourceUrl
                    ? `<a href="${{escapeHtml(sourceUrl)}}" target="_blank" rel="noopener noreferrer">${{escapeHtml(sourceLabel || sourceUrl)}}</a>`
                    : escapeHtml(sourceLabel);
                const context = [e.study_design, e.species, e.sex].filter(value => value != null && String(value).trim());
                const evidenceRows = [
                    ['Endpoint', e.endpoint], ['Comparator', e.comparator],
                    ['Effect estimate', e.effect_estimate], ['Uncertainty', e.uncertainty],
                    ['Source location', e.source_locator], ['Evidence basis', e.evidence_basis],
                    ['Historical catalog review label', e.source_review_status]
                ].filter(([, value]) => value != null && String(value).trim())
                 .map(([label, value]) => `<div><strong>${{escapeHtml(label)}}:</strong> ${{escapeHtml(value)}}</div>`)
                 .join('');
                return `
                <div class="edge-pill">
                    <div class="edge-top">
                        <span class="hallmark-tag">${{escapeHtml(e.hallmark ? String(e.hallmark).replace(/_/g, ' ') : 'unspecified')}}</span>
                        <span class="grade-badge grade-${{grade}}">${{grade}}</span>
                    </div>
                    <div class="edge-notes">
                        <strong>${{escapeHtml(e.target_symbol || e.target_id || '')}}</strong> ${{escapeHtml(e.relation)}} &bull; Recorded claim: ${{escapeHtml(e.notes || '')}}
                        <div class="review-status">${{escapeHtml(edgeReviewLabel(e))}}</div>
                        ${{(e.citation_audit_identifiers || []).map(item =>
                            `<div>${{escapeHtml(item.identifier_kind)}} ${{escapeHtml(item.identifier)}}: ${{escapeHtml(item.screen)}} — ${{escapeHtml(item.screen_note)}}</div>`
                        ).join('')}}
                        ${{e.document_ids ? `<div style="margin-top: 2px; color: #60a5fa;">${{escapeHtml(e.document_ids)}}</div>` : ''}}
                        ${{sourceHtml ? `<div style="margin-top: 2px;">Source: ${{sourceHtml}}</div>` : ''}}
                        ${{context.length ? `<div style="margin-top: 2px;">Context: ${{context.map(escapeHtml).join(' · ')}}</div>` : ''}}
                        ${{evidenceRows ? `<div style="margin-top: 4px;">${{evidenceRows}}</div>` : ''}}
                    </div>
                </div>
            `;
            }}).join('');

            card.innerHTML = `
                <div class="card-header">
                    <div>
                        <div class="card-title">${{comp.resolved_name || comp.raw_name}}</div>
                        <div class="card-subtitle">${{comp.raw_name !== comp.resolved_name ? 'Ingest: ' + comp.raw_name : ''}}</div>
                    </div>
                    <span class="modality-pill ${{comp.modality}}">${{comp.modality.replace('_', ' ')}}</span>
                </div>

                <div class="svg-container">
                    ${{comp.svg}}
                </div>

                <div class="identifiers-box">
                    <div class="id-row">
                        <span class="id-label">InChIKey:</span>
                        <span class="id-val" title="${{comp.inchikey}}">${{comp.inchikey}}</span>
                        <button class="copy-btn" onclick="navigator.clipboard.writeText('${{comp.inchikey}}')">copy</button>
                    </div>
                    <div class="id-row">
                        <span class="id-label">PubChem:</span>
                        <span class="id-val">${{comp.cid ? `<a href="https://pubchem.ncbi.nlm.nih.gov/compound/${{comp.cid}}" target="_blank" style="color:#60a5fa;text-decoration:none;">CID ${{comp.cid}}</a>` : 'Unresolved'}}</span>
                    </div>
                    <div class="id-row">
                        <span class="id-label">ChEMBL:</span>
                        <span class="id-val">${{comp.has_chembl ? `<span style="color:#34d399;">${{comp.num_activities}} assays (${{comp.num_binding_assays}} binding)</span>` : '<span style="color:#f87171;">None</span>'}}</span>
                    </div>
                </div>

                <div class="properties-grid">
                    <div class="prop-item"><span class="prop-key">Mol Weight</span><span class="prop-val">${{comp.mol_wt || 'N/A'}}</span></div>
                    <div class="prop-item"><span class="prop-key">LogP</span><span class="prop-val">${{comp.log_p !== null ? comp.log_p : 'N/A'}}</span></div>
                    <div class="prop-item"><span class="prop-key">TPSA (Å²)</span><span class="prop-val">${{comp.tpsa || 'N/A'}}</span></div>
                    <div class="prop-item"><span class="prop-key">H-Donors</span><span class="prop-val">${{comp.num_h_donors !== null ? comp.num_h_donors : 'N/A'}}</span></div>
                    <div class="prop-item"><span class="prop-key">H-Acceptors</span><span class="prop-val">${{comp.num_h_acceptors !== null ? comp.num_h_acceptors : 'N/A'}}</span></div>
                    <div class="prop-item" title="QED heuristic: Fitted on historical approved oral small molecules; systematically unkind to peptides/macrocycles."><span class="prop-key">QED*</span><span class="prop-val" style="color: #60a5fa;">${{comp.qed !== null ? comp.qed : 'N/A'}}</span></div>
                </div>

                <div class="evidence-section">
                    <div class="evidence-header">
                        <span>Curated Assertions · unreviewed support is E0</span>
                        <span>${{comp.edges.length}} Edges</span>
                    </div>
                    ${{edgesHtml || '<div style="font-size:0.75rem; color:#64748b;">No curated edges</div>'}}
                </div>
            `;
            container.appendChild(card);
        }});
    }}

    function renderFalsePositives() {{
        const tbody = document.getElementById('fps-tbody');
        tbody.innerHTML = '';

        const explanations = {{
            "CHEMBL1088790": "Morpholino-pyrimidine hinge motif; true affinity is 3.6 uM (barely misses strict 1.0 uM active cutoff).",
            "CHEMBL98350": "Chromen-4-one core resembling LY294002 analog series; true IC50 is 3.16 uM.",
            "CHEMBL1088831": "Sub-micromolar borderline; true affinity is 2.0 uM (pChEMBL 5.70). Intact dual kinase pharmacophore.",
            "CHEMBL591338": "Steric clash penalty in deep mTOR active site pocket drops affinity to 1.7 uM.",
            "CHEMBL1088832": "Fluorinated benzyl analog displaying similar borderline micromolar activity (3.6 uM).",
            "CHEMBL188678": "Bulky dibenzothiophene core shifts angle of morpholine entry into ATP binding pocket (1.7 uM).",
            "CHEMBL435507": "Truncated chromone scaffold retains weak residual binding (4.8 uM)."
        }};

        topFPsData.forEach(fp => {{
            const tr = document.createElement('tr');
            const exp = explanations[fp.molecule_chembl_id] || "Near-threshold borderline affinity; pharmacophore matches active kinase inhibitors.";
            tr.innerHTML = `
                <td><a href="https://www.ebi.ac.uk/chembl/compound_report_card/${{fp.molecule_chembl_id}}" target="_blank" style="color:#60a5fa;text-decoration:none;font-weight:600;">${{fp.molecule_chembl_id}}</a></td>
                <td style="font-weight:700; color:#f87171;">${{fp.predicted_prob.toFixed(4)}}</td>
                <td>${{fp.true_pchembl}}</td>
                <td style="font-family:'JetBrains Mono',monospace;font-size:0.75rem;max-width:240px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;" title="${{fp.murcko_scaffold}}">${{fp.murcko_scaffold}}</td>
                <td style="font-size:0.825rem;color:#cbd5e1;">${{exp}}</td>
            `;
            tbody.appendChild(tr);
        }});
    }}

    function filterCards() {{
        const searchVal = document.getElementById('search-input').value.toLowerCase().trim();
        const hallmarkVal = document.getElementById('hallmark-filter').value.toLowerCase();
        const gradeVal = document.getElementById('grade-filter').value;
        const modalityVal = document.getElementById('modality-filter').value;
        const chemblVal = document.getElementById('chembl-filter') ? document.getElementById('chembl-filter').value : '';

        const gradeRank = {{ "E0": 0, "E1": 1, "E2": 2, "E3": 3, "E4": 4 }};
        const minGrade = gradeVal ? gradeRank[gradeVal] : 0;

        const filtered = compoundsData.filter(comp => {{
            // Multi-field search text matching
            const matchesSearch = !searchVal || 
                (comp.raw_name && comp.raw_name.toLowerCase().includes(searchVal)) ||
                (comp.resolved_name && comp.resolved_name.toLowerCase().includes(searchVal)) ||
                (comp.inchikey && comp.inchikey.toLowerCase().includes(searchVal)) ||
                (comp.canonical_smiles && comp.canonical_smiles.toLowerCase().includes(searchVal)) ||
                (comp.cid && String(comp.cid).includes(searchVal)) ||
                (comp.synonyms && comp.synonyms.some(s => s.toLowerCase().includes(searchVal))) ||
                (comp.edges && comp.edges.some(e => 
                    (e.target_symbol && e.target_symbol.toLowerCase().includes(searchVal)) ||
                    (e.hallmark && e.hallmark.toLowerCase().includes(searchVal)) ||
                    (e.notes && e.notes.toLowerCase().includes(searchVal)) ||
                    (e.document_ids && e.document_ids.toLowerCase().includes(searchVal))
                ));

            // Modality
            const matchesModality = !modalityVal || comp.modality === modalityVal;

            // Hallmark
            const matchesHallmark = !hallmarkVal || comp.edges.some(e => e.hallmark === hallmarkVal);

            // Grade
            const matchesGrade = !gradeVal || comp.edges.some(e => gradeRank[e.grade] >= minGrade);

            // ChEMBL status
            let matchesChembl = true;
            if (chemblVal === 'binding') {{
                matchesChembl = comp.has_binding && comp.num_binding_assays > 0;
            }} else if (chemblVal === 'has_chembl') {{
                matchesChembl = comp.has_chembl;
            }} else if (chemblVal === 'none') {{
                matchesChembl = !comp.has_chembl;
            }}

            return matchesSearch && matchesModality && matchesHallmark && matchesGrade && matchesChembl;
        }});

        renderCards(filtered);
    }}

    function resetFilters() {{
        document.getElementById('search-input').value = '';
        document.getElementById('hallmark-filter').value = '';
        document.getElementById('grade-filter').value = '';
        document.getElementById('modality-filter').value = '';
        if (document.getElementById('chembl-filter')) document.getElementById('chembl-filter').value = '';
        filterCards();
    }}

    function applyPreset(preset) {{
        resetFilters();
        if (preset === 'E4') {{
            document.getElementById('grade-filter').value = 'E4';
        }} else if (preset === 'peptide') {{
            document.getElementById('modality-filter').value = 'peptide';
        }} else if (preset) {{
            document.getElementById('search-input').value = preset;
        }}
        filterCards();
    }}

    function filterGenCards() {{
        const searchInput = document.getElementById('gen-search-input');
        const searchVal = searchInput ? searchInput.value.toLowerCase().trim() : '';
        const painsFilter = document.getElementById('gen-pains-filter');
        const painsVal = painsFilter ? painsFilter.value : '';
        const sortSelect = document.getElementById('gen-sort-select');
        const sortVal = sortSelect ? sortSelect.value : 'reward';

        let filtered = genCardsData.filter(comp => {{
            const matchesSearch = !searchVal ||
                (comp.inchikey && comp.inchikey.toLowerCase().includes(searchVal)) ||
                String(comp.generation).includes(searchVal);
            const matchesPains = !painsVal ||
                (painsVal === 'pass' && !comp.has_pains) ||
                (painsVal === 'flagged' && comp.has_pains);
            return matchesSearch && matchesPains;
        }});

        if (sortVal === 'mtor') {{
            filtered.sort((a, b) => b.mtor_prob - a.mtor_prob);
        }} else if (sortVal === 'qed') {{
            filtered.sort((a, b) => b.qed - a.qed);
        }} else {{
            filtered.sort((a, b) => b.reward - a.reward);
        }}

        const countEl = document.getElementById('gen-visible-count');
        if (countEl) countEl.innerText = filtered.length;
        renderGenCards(filtered);
    }}

    function switchTab(tabId) {{
        document.querySelectorAll('.tab-btn').forEach(btn => btn.classList.remove('active'));
        document.getElementById('tab-cards').classList.add('hidden');
        document.getElementById('tab-holes').classList.add('hidden');
        document.getElementById('tab-bench').classList.add('hidden');
        document.getElementById('tab-errors').classList.add('hidden');
        document.getElementById('tab-generator').classList.add('hidden');

        if (tabId === 'cards') {{
            document.querySelector('.tab-btn:nth-child(1)').classList.add('active');
            document.getElementById('tab-cards').classList.remove('hidden');
        }} else if (tabId === 'holes') {{
            document.querySelector('.tab-btn:nth-child(2)').classList.add('active');
            document.getElementById('tab-holes').classList.remove('hidden');
        }} else if (tabId === 'bench') {{
            document.querySelector('.tab-btn:nth-child(3)').classList.add('active');
            document.getElementById('tab-bench').classList.remove('hidden');
        }} else if (tabId === 'errors') {{
            document.querySelector('.tab-btn:nth-child(4)').classList.add('active');
            document.getElementById('tab-errors').classList.remove('hidden');
        }} else if (tabId === 'generator') {{
            document.querySelector('.tab-btn:nth-child(5)').classList.add('active');
            document.getElementById('tab-generator').classList.remove('hidden');
        }}
    }}

    function renderGenCards(data) {{
        const container = document.getElementById('gen-cards-container');
        if (!container) return;
        container.innerHTML = '';
        if (!data || data.length === 0) {{
            container.innerHTML = '<div style="grid-column: 1/-1; text-align: center; padding: 3rem; color: #64748b;">No generated candidates match criteria.</div>';
            return;
        }}
        data.forEach((comp, idx) => {{
            const card = document.createElement('div');
            card.className = 'compound-card';
            card.innerHTML = `
                <div class="card-header">
                    <div>
                        <div class="card-title">Candidate #${{idx + 1}}</div>
                        <div class="card-subtitle">Evolved via GA (Gen ${{comp.generation}})</div>
                    </div>
                    <span class="modality-pill small_molecule">Generated</span>
                </div>
                <div class="svg-container">
                    ${{comp.svg}}
                </div>
                <div class="identifiers-box">
                    <div class="id-row">
                        <span class="id-label">InChIKey:</span>
                        <span class="id-val" title="${{comp.inchikey}}">${{comp.inchikey}}</span>
                    </div>
                    <div class="id-row">
                        <span class="id-label">PAINS:</span>
                        <span class="id-val">${{comp.has_pains ? '<span style="color:#ef4444;">Flagged</span>' : '<span style="color:#34d399;">Pass (Zero Motif)</span>'}}</span>
                    </div>
                </div>
                <div class="properties-grid">
                    <div class="prop-item"><span class="prop-key">mTOR P(act)</span><span class="prop-val" style="color:#60a5fa;">${{comp.mtor_prob}}</span></div>
                    <div class="prop-item"><span class="prop-key">QED</span><span class="prop-val">${{comp.qed}}</span></div>
                    <div class="prop-item"><span class="prop-key">Composite</span><span class="prop-val" style="color:#34d399;">${{comp.reward}}</span></div>
                </div>
            `;
            container.appendChild(card);
        }});
    }}

    // Initial render
    renderCards(compoundsData);
    renderFalsePositives();
    renderGenCards(genCardsData);
</script>
</body>
</html>
"""

    out_file = Path(output_html)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file, "w", encoding="utf-8") as f:
        f.write(html_content)

    # A custom preview/test output must not implicitly overwrite the hosted page.
    if site_html is not None:
        site_file = Path(site_html)
        site_file.parent.mkdir(parents=True, exist_ok=True)
        site_file.write_text(html_content, encoding="utf-8")
        logger.info("Saved hosted dashboard copy to %s", site_file)

    logger.info("Saved interactive dashboard to %s", out_file)
    return out_file


if __name__ == "__main__":
    build_dashboard_html(site_html="site/index.html")
