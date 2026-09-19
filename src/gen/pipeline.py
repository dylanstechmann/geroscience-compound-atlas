"""Pipeline for molecular generation, bias sensitivity analysis, and chemical space PCA."""

import json
import logging
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.decomposition import PCA

from gen.ga import MolecularGA
from gen.scorer import CompositeSurrogateScorer

logger = logging.getLogger(__name__)


def compute_internal_diversity(smiles_list: list[str]) -> float:
    """Compute internal diversity (1 - mean pairwise Tanimoto similarity) over 2048-bit Morgan FPs."""
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = []
    for s in smiles_list:
        m = Chem.MolFromSmiles(s)
        if m:
            fps.append(mfpgen.GetFingerprint(m))

    n = len(fps)
    if n < 2:
        return 0.0

    similarities = []
    for i in range(n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1 :])
        similarities.extend(sims)

    if not similarities:
        return 0.0
    return float(1.0 - np.mean(similarities))


def run_generator_pipeline(
    benchmark_parquet: str | Path = "artifacts/benchmark_dataset.parquet",
    compounds_parquet: str | Path = "artifacts/compounds.parquet",
    splits_json: str | Path = "artifacts/splits.json",
    output_parquet: str | Path = "artifacts/generated_molecules.parquet",
    output_metrics: str | Path = "artifacts/generator_metrics.json",
    output_figure: str | Path = "figures/generator_chemical_space.png",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Execute molecular GA under variable QED bias, compute metrics, and project PCA chemical space."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    logger.info("Initializing CompositeSurrogateScorer...")
    scorer = CompositeSurrogateScorer(benchmark_parquet=benchmark_parquet, splits_json=splits_json)

    bench_df = pd.read_parquet(benchmark_parquet)
    training_inchikeys = set(bench_df["inchikey"].dropna())

    # Select top 25 active molecules as seeds
    actives_df = bench_df[bench_df["active"] == 1].sort_values(by="pchembl_value", ascending=False)
    seed_smiles = actives_df["canonical_smiles"].dropna().head(25).tolist()

    # 1. Evaluate QED Bias Sensitivity across lambda_qed in [0.0, 0.2, 0.5]
    qed_weights = [0.0, 0.2, 0.5]
    sensitivity_results = {}
    primary_df = pd.DataFrame()
    primary_stats = {}

    for w in qed_weights:
        logger.info("Running Genetic Algorithm with lambda_qed = %.1f...", w)
        ga = MolecularGA(
            scorer=scorer,
            pop_size=40,
            n_generations=20,
            mutation_rate=0.7,
            crossover_rate=0.3,
            elite_size=5,
            qed_weight=w,
            seed=42,
        )
        archive_df, stats = ga.run(seed_smiles=seed_smiles, max_archive_size=200)

        mean_mtor = float(archive_df["mtor_prob"].mean()) if not archive_df.empty else 0.0
        mean_qed = float(archive_df["qed"].mean()) if not archive_df.empty else 0.0
        novelty = (
            float(
                sum(k not in training_inchikeys for k in archive_df["inchikey"]) / len(archive_df)
            )
            if not archive_df.empty
            else 0.0
        )
        diversity = (
            compute_internal_diversity(archive_df["smiles"].tolist())
            if not archive_df.empty
            else 0.0
        )

        sensitivity_results[f"qed_weight_{w:.1f}"] = {
            "qed_weight": w,
            "generated_count": len(archive_df),
            "mean_mtor_prob": mean_mtor,
            "mean_qed": mean_qed,
            "novelty_rate": novelty,
            "internal_diversity": diversity,
            "validity_rate": stats["validity_rate"],
        }

        if w == 0.2:
            primary_df = archive_df.copy()
            primary_stats = stats

    # Export generated structures
    Path(output_parquet).parent.mkdir(parents=True, exist_ok=True)
    primary_df.to_parquet(output_parquet, index=False)
    csv_path = Path("data/processed/generated_molecules.csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    primary_df.to_csv(csv_path, index=False)
    logger.info("Exported %d generated molecules to %s", len(primary_df), output_parquet)

    # 2. PCA Chemical Space Projection
    # Prepare datasets:
    # A. Generated molecules
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    gen_fps = []
    gen_smiles_valid = []
    for s in primary_df["smiles"]:
        m = Chem.MolFromSmiles(s)
        if m:
            gen_fps.append(np.array(mfpgen.GetFingerprint(m), dtype=np.float32))
            gen_smiles_valid.append(s)

    # B. ChEMBL mTOR actives
    active_mols = bench_df[bench_df["active"] == 1]
    active_fps_list = []
    for s in active_mols["canonical_smiles"]:
        m = Chem.MolFromSmiles(str(s)) if pd.notna(s) else None
        if m:
            active_fps_list.append(np.array(mfpgen.GetFingerprint(m), dtype=np.float32))
    active_fps = np.array(active_fps_list, dtype=np.float32)

    # C. Curated E3/E4 Atlas Geroprotectors
    atlas_df = pd.read_parquet(compounds_parquet)
    known_gero = atlas_df[atlas_df["resolved"] & atlas_df["canonical_smiles"].notna()].copy()
    gero_fps = []
    gero_names = []
    for _, row in known_gero.iterrows():
        m = Chem.MolFromSmiles(row["canonical_smiles"])
        if m:
            gero_fps.append(np.array(mfpgen.GetFingerprint(m), dtype=np.float32))
            gero_names.append(row["resolved_name"] or row["raw_name"])

    # Combine for PCA fitting
    all_fps = np.vstack([active_fps, np.array(gen_fps), np.array(gero_fps)])
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(all_fps)

    n_act = len(active_fps)
    n_gen = len(gen_fps)
    n_gero = len(gero_fps)

    coords_act = coords[:n_act]
    coords_gen = coords[n_act : n_act + n_gen]
    coords_gero = coords[n_act + n_gen :]

    # Plot PCA Chemical Space
    Path(output_figure).parent.mkdir(parents=True, exist_ok=True)
    plt.style.use(
        "seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default"
    )
    fig, ax = plt.subplots(figsize=(10, 7), dpi=300)

    ax.scatter(
        coords_act[:, 0],
        coords_act[:, 1],
        c="#94a3b8",
        alpha=0.45,
        s=30,
        label=f"ChEMBL mTOR Training Actives (N={n_act})",
        edgecolors="none",
    )
    ax.scatter(
        coords_gen[:, 0],
        coords_gen[:, 1],
        c="#3b82f6",
        alpha=0.75,
        s=45,
        label=f"GA Generated Candidates (N={n_gen}, $\\lambda_{{QED}}=0.2$)",
        edgecolors="#1d4ed8",
        linewidth=0.5,
    )
    ax.scatter(
        coords_gero[:, 0],
        coords_gero[:, 1],
        c="#ef4444",
        alpha=0.95,
        s=80,
        marker="^",
        label=f"Known Atlas Geroprotective Chemotypes (N={n_gero})",
        edgecolors="#7f1d1d",
        linewidth=1.0,
        zorder=5,
    )

    # Annotate landmark geroprotectors
    for i, name in enumerate(gero_names):
        if name in {
            "rapamycin",
            "metformin",
            "dasatinib",
            "canagliflozin",
            "acarbose",
            "navitoclax",
        }:
            ax.annotate(
                name.capitalize(),
                (coords_gero[i, 0], coords_gero[i, 1]),
                xytext=(6, 4),
                textcoords="offset points",
                fontsize=8,
                fontweight="bold",
                color="#991b1b",
                bbox={
                    "boxstyle": "round,pad=0.2",
                    "fc": "white",
                    "ec": "#ef4444",
                    "alpha": 0.8,
                    "lw": 0.8,
                },
            )

    var1 = pca.explained_variance_ratio_[0] * 100
    var2 = pca.explained_variance_ratio_[1] * 100
    ax.set_title(
        "Chemical Space Exploration — Generated Candidates vs mTOR Actives vs Atlas Geroprotectors",
        fontsize=12,
        fontweight="bold",
        pad=12,
    )
    ax.set_xlabel(f"Principal Component 1 ({var1:.1f}% explained variance)", fontsize=10)
    ax.set_ylabel(f"Principal Component 2 ({var2:.1f}% explained variance)", fontsize=10)
    ax.legend(frameon=True, loc="upper right", framealpha=0.9, fontsize=9)
    plt.tight_layout()
    fig.savefig(output_figure, dpi=300)
    plt.close(fig)
    logger.info("Saved PCA chemical space visualization to %s", output_figure)

    # 3. Assemble Metrics Summary
    primary_novelty = (
        float(sum(k not in training_inchikeys for k in primary_df["inchikey"]) / len(primary_df))
        if not primary_df.empty
        else 0.0
    )
    primary_diversity = (
        compute_internal_diversity(primary_df["smiles"].tolist()) if not primary_df.empty else 0.0
    )

    metrics_payload = {
        "primary_run": {
            "qed_weight": 0.2,
            "num_generated": len(primary_df),
            "validity_rate": primary_stats.get("validity_rate", 0.0),
            "uniqueness_rate": 1.0,  # Archive deduplicates by canonical InChIKey
            "internal_diversity": primary_diversity,
            "novelty_vs_training": primary_novelty,
            "mean_mtor_prob": float(primary_df["mtor_prob"].mean()),
            "mean_qed": float(primary_df["qed"].mean()),
            "pains_pass_rate": float((~primary_df["has_pains"]).mean()),
        },
        "qed_bias_sensitivity": sensitivity_results,
        "pca_variance_explained": {
            "pc1": float(pca.explained_variance_ratio_[0]),
            "pc2": float(pca.explained_variance_ratio_[1]),
        },
    }

    with open(output_metrics, "w", encoding="utf-8") as f:
        json.dump(metrics_payload, f, indent=2)
    logger.info("Exported generator metrics to %s", output_metrics)

    return primary_df, metrics_payload


if __name__ == "__main__":
    run_generator_pipeline()
