"""Single-compound lookup for Geroscience Compound Atlas.

Queries compound metadata, physicochemical properties, curated aging hallmark
evidence edges, and predictive benchmark (mTOR surrogate) scoring.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem

from atlas.features import compute_rdkit_descriptors
from gen.scorer import CompositeSurrogateScorer

logger = logging.getLogger("atlas.lookup")


def lookup_compound(
    query: str,
    *,
    compounds_parquet: str | Path = "artifacts/compounds.parquet",
    edges_parquet: str | Path = "artifacts/evidence_edges.parquet",
    coverage_parquet: str | Path = "artifacts/chembl_coverage.parquet",
    activities_parquet: str | Path = "artifacts/chembl_activities.parquet",
    benchmark_parquet: str | Path = "artifacts/benchmark_dataset.parquet",
    splits_json: str | Path = "artifacts/splits.json",
    scorer: CompositeSurrogateScorer | None = None,
) -> dict[str, Any]:
    """Look up a compound by name, synonym, CID, InChIKey, or SMILES.

    Returns:
        Structured dictionary with compound card, hallmark associations,
        and benchmark predictions.
    """
    comp_file = Path(compounds_parquet)
    if not comp_file.exists():
        raise FileNotFoundError(f"Compounds table not found at {comp_file}")

    compounds_df = pd.read_parquet(comp_file)
    q_str = str(query).strip()
    q_lower = q_str.lower()

    # Search in compounds_df
    matched_row: pd.Series | None = None

    # 1. Exact match by name or resolved_name
    for col in ["name", "resolved_name", "raw_name"]:
        if col in compounds_df.columns:
            m = compounds_df[compounds_df[col].astype(str).str.lower() == q_lower]
            if not m.empty:
                matched_row = m.iloc[0]
                break

    # 2. InChIKey match
    if matched_row is None and "inchikey" in compounds_df.columns:
        m = compounds_df[compounds_df["inchikey"].astype(str).str.lower() == q_lower]
        if not m.empty:
            matched_row = m.iloc[0]

    # 3. CID match
    if matched_row is None and q_str.isdigit() and "cid" in compounds_df.columns:
        m = compounds_df[compounds_df["cid"] == int(q_str)]
        if not m.empty:
            matched_row = m.iloc[0]

    # 4. Canonical SMILES match
    if matched_row is None and "canonical_smiles" in compounds_df.columns:
        m = compounds_df[compounds_df["canonical_smiles"].astype(str) == q_str]
        if not m.empty:
            matched_row = m.iloc[0]

    # 5. Synonyms search
    if matched_row is None and "synonyms" in compounds_df.columns:
        for _, row in compounds_df.iterrows():
            syns = row.get("synonyms")
            if isinstance(syns, (list, np.ndarray)):
                if any(str(s).lower() == q_lower for s in syns):
                    matched_row = row
                    break

    # 6. Fallback: Parse query as SMILES if RDKit recognizes it
    in_atlas = matched_row is not None
    if in_atlas:
        name = str(matched_row.get("resolved_name") or matched_row.get("name") or q_str)
        raw_name = str(matched_row.get("raw_name") or name)
        inchikey = str(matched_row.get("inchikey") or "")
        cid = int(matched_row["cid"]) if pd.notna(matched_row.get("cid")) else None
        smiles = str(matched_row.get("canonical_smiles") or "")
        modality = str(matched_row.get("modality") or "small_molecule")
    else:
        mol_test = Chem.MolFromSmiles(q_str)
        if mol_test is not None:
            smiles = Chem.MolToSmiles(mol_test)
            inchikey = Chem.MolToInchiKey(mol_test)
            name = f"Custom Structure ({inchikey[:14]})"
            raw_name = q_str
            cid = None
            modality = "small_molecule"
        else:
            raise ValueError(f"Compound '{query}' not found in atlas and is not a valid SMILES.")

    # Molecular properties
    descriptors = {}
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol:
        descriptors = compute_rdkit_descriptors(mol)
    elif in_atlas:
        for desc_col in ["mol_wt", "log_p", "tpsa", "num_h_donors", "num_h_acceptors", "qed"]:
            if desc_col in matched_row and pd.notna(matched_row[desc_col]):
                descriptors[desc_col] = float(matched_row[desc_col])

    # Hallmark evidence edges
    edges_file = Path(edges_parquet)
    hallmark_edges = []
    if edges_file.exists():
        edges_df = pd.read_parquet(edges_file)
        if "compound_inchikey" in edges_df.columns:
            matched_edges = edges_df[edges_df["compound_inchikey"] == inchikey]
        elif "inchikey" in edges_df.columns:
            matched_edges = edges_df[edges_df["inchikey"] == inchikey]
        else:
            matched_edges = pd.DataFrame()

        for _, e in matched_edges.iterrows():
            hallmark_edges.append({
                "hallmark": str(e.get("hallmark", "unspecified")),
                "grade": str(e.get("grade", "E0")),
                "target_symbol": str(e.get("target_symbol") or e.get("target_id") or "N/A"),
                "relation": str(e.get("relation", "modulates")),
                "notes": str(e.get("notes", "")),
                "document_ids": str(e.get("document_ids", "")),
            })

    # ChEMBL coverage & measured activities
    chembl_info: dict[str, Any] = {"has_chembl": False, "num_activities": 0, "activities": []}
    cov_file = Path(coverage_parquet)
    if cov_file.exists():
        cov_df = pd.read_parquet(cov_file)
        cov_match = cov_df[cov_df["inchikey"] == inchikey]
        if not cov_match.empty:
            c_row = cov_match.iloc[0]
            chembl_info["has_chembl"] = bool(c_row.get("has_chembl", False))
            chembl_info["num_activities"] = int(c_row.get("num_activities", 0))
            chembl_info["has_binding"] = bool(c_row.get("has_binding", False))
            chembl_info["num_binding_assays"] = int(c_row.get("num_binding_assays", 0))

    acts_file = Path(activities_parquet)
    if acts_file.exists():
        acts_df = pd.read_parquet(acts_file)
        if "inchikey" in acts_df.columns:
            act_match = acts_df[acts_df["inchikey"] == inchikey]
            if not act_match.empty:
                chembl_info["has_chembl"] = True
                chembl_info["num_activities"] = len(act_match)
                for _, a in act_match.head(10).iterrows():
                    chembl_info["activities"].append({
                        "target_name": str(a.get("target_name") or a.get("pref_name") or "Target"),
                        "standard_type": str(a.get("standard_type") or "IC50"),
                        "standard_value": float(a["standard_value"]) if pd.notna(a.get("standard_value")) else None,
                        "standard_units": str(a.get("standard_units") or "nM"),
                        "pchembl_value": float(a["pchembl_value"]) if pd.notna(a.get("pchembl_value")) else None,
                    })

    # Benchmark prediction (mTOR kinase surrogate)
    prediction_info: dict[str, Any] = {
        "available": False,
        "model": "LogisticRegression (Bemis-Murcko scaffold trained baseline)",
        "task": "mTOR kinase inhibition (pChEMBL >= 6.0)",
    }
    if mol:
        try:
            if scorer is None:
                scorer = CompositeSurrogateScorer(
                    benchmark_parquet=benchmark_parquet,
                    splits_json=splits_json,
                )
            score_res = scorer.score_molecule(mol)
            prediction_info["available"] = True
            prediction_info["mtor_prob"] = round(float(score_res["mtor_prob"]), 4)
            prediction_info["prediction"] = "ACTIVE" if score_res["mtor_prob"] >= 0.5 else "INACTIVE"
            prediction_info["qed"] = round(float(score_res["qed"]), 3)
            prediction_info["has_pains"] = bool(score_res["has_pains"])
            prediction_info["pains_status"] = "FLAGGED" if score_res["has_pains"] else "PASS (Zero PAINS motifs)"
            prediction_info["composite_reward"] = round(float(score_res["reward"]), 4)
            prediction_info["disclaimer"] = (
                "Surrogate in-silico prediction from frozen scaffold benchmark. "
                "Not experimental confirmation, medical guidance, or a protocol recommendation."
            )
        except Exception as exc:
            logger.debug("Benchmark prediction failed: %s", exc)
            prediction_info["error"] = str(exc)

    return {
        "query": query,
        "in_atlas": in_atlas,
        "name": name,
        "raw_name": raw_name,
        "inchikey": inchikey,
        "cid": cid,
        "canonical_smiles": smiles,
        "modality": modality,
        "properties": descriptors,
        "hallmark_edges": hallmark_edges,
        "chembl_info": chembl_info,
        "benchmark_prediction": prediction_info,
    }


def format_lookup_report(data: dict[str, Any]) -> str:
    """Format single compound lookup result as an auditable text report."""
    lines = []
    lines.append("=" * 78)
    lines.append(f"  GEROSCIENCE COMPOUND EVIDENCE CARD: {data['name']}")
    lines.append("=" * 78)

    # Identifiers
    lines.append("\n[IDENTIFIERS]")
    lines.append(f"  Common Name     : {data['name']}")
    if data["raw_name"] != data["name"]:
        lines.append(f"  Ingested Name   : {data['raw_name']}")
    lines.append(f"  Modality        : {data['modality']}")
    lines.append(f"  InChIKey        : {data['inchikey']}")
    cid_str = f"CID {data['cid']} (https://pubchem.ncbi.nlm.nih.gov/compound/{data['cid']})" if data.get("cid") else "Unresolved"
    lines.append(f"  PubChem         : {cid_str}")
    lines.append(f"  Canonical SMILES: {data['canonical_smiles']}")
    lines.append(f"  Atlas Membership: {'Verified watchlist entity' if data['in_atlas'] else 'Custom molecular query'}")

    # Physicochemical Properties
    props = data.get("properties", {})
    if props:
        lines.append("\n[PHYSICOCHEMICAL PROPERTIES]")
        if "mol_wt" in props:
            lines.append(f"  Molecular Weight : {props['mol_wt']:.2f} g/mol")
        if "log_p" in props:
            lines.append(f"  Calculated LogP  : {props['log_p']:.2f}")
        if "tpsa" in props:
            lines.append(f"  TPSA             : {props['tpsa']:.1f} Å²")
        if "num_h_donors" in props and "num_h_acceptors" in props:
            lines.append(f"  H-Donors / Accept: {int(props['num_h_donors'])} / {int(props['num_h_acceptors'])}")
        if "qed" in props:
            lines.append(f"  QED Score        : {props['qed']:.3f} (Drug-likeness heuristic)")

    # Hallmark Associations
    edges = data.get("hallmark_edges", [])
    lines.append(f"\n[AGING HALLMARK ASSOCIATIONS ({len(edges)} curated edges)]")
    if edges:
        for idx, e in enumerate(edges, 1):
            h_clean = e['hallmark'].replace('_', ' ').title()
            lines.append(f"  • [{e['grade']}] {h_clean}")
            lines.append(f"    Target    : {e['target_symbol']} ({e['relation']})")
            if e.get("notes"):
                lines.append(f"    Evidence  : {e['notes']}")
            if e.get("document_ids"):
                lines.append(f"    Citations : {e['document_ids']}")
    else:
        lines.append("  (No curated evidence edges in current atlas catalog)")

    # ChEMBL Bioactivities
    c_info = data.get("chembl_info", {})
    lines.append("\n[CHEMBL EXPERIMENTAL BIOACTIVITIES]")
    if c_info.get("has_chembl"):
        lines.append(f"  ChEMBL Coverage : {c_info.get('num_activities', 0)} recorded assays ({c_info.get('num_binding_assays', 0)} binding)")
        acts = c_info.get("activities", [])
        for a in acts[:5]:
            p_val = f", pChEMBL={a['pchembl_value']:.2f}" if a.get("pchembl_value") else ""
            lines.append(f"    - {a['target_name']}: {a['standard_type']} = {a.get('standard_value')} {a['standard_units']}{p_val}")
    else:
        lines.append("  No matching structured ChEMBL bioactivities mapped.")

    # Benchmark Prediction
    pred = data.get("benchmark_prediction", {})
    lines.append("\n[PREDICTIVE BENCHMARK INFERENCE (mTOR Kinase)]")
    if pred.get("available"):
        lines.append(f"  Model Architecture  : {pred['model']}")
        lines.append(f"  Predicted Activity  : {pred['prediction']} (P(Active) = {pred['mtor_prob']:.4f})")
        lines.append(f"  PAINS Filter Status : {pred['pains_status']}")
        lines.append(f"  Composite GA Score  : {pred.get('composite_reward', 'N/A')}")
        lines.append(f"  Cautionary Note     : {pred['disclaimer']}")
    elif "error" in pred:
        lines.append(f"  Prediction Unavailable: {pred['error']}")
    else:
        lines.append("  Structure unavailable for in-silico inference.")

    lines.append("=" * 78)
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Geroscience Compound Atlas: Evidence Card & Benchmark Lookup",
    )
    parser.add_argument("compound", help="Compound name, synonym, CID, InChIKey, or canonical SMILES")
    parser.add_argument("--json", action="store_true", help="Output machine-readable JSON instead of formatted text")
    parser.add_argument("--out", type=str, help="Optional output file destination")
    parser.add_argument("--compounds", default="artifacts/compounds.parquet", help="Path to compounds table")
    parser.add_argument("--edges", default="artifacts/evidence_edges.parquet", help="Path to evidence edges table")
    parser.add_argument("--coverage", default="artifacts/chembl_coverage.parquet", help="Path to ChEMBL coverage table")
    parser.add_argument("--benchmark", default="artifacts/benchmark_dataset.parquet", help="Path to benchmark parquet")
    parser.add_argument("--splits", default="artifacts/splits.json", help="Path to splits json")

    args = parser.parse_args(argv)

    try:
        result = lookup_compound(
            args.compound,
            compounds_parquet=args.compounds,
            edges_parquet=args.edges,
            coverage_parquet=args.coverage,
            benchmark_parquet=args.benchmark,
            splits_json=args.splits,
        )

        if args.json:
            output_str = json.dumps(result, indent=2, default=str)
        else:
            output_str = format_lookup_report(result)

        if args.out:
            out_path = Path(args.out)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(output_str, encoding="utf-8")
            print(f"Report saved to {out_path}")
        else:
            print(output_str)

        return 0
    except Exception as exc:
        print(f"Error looking up compound '{args.compound}': {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
