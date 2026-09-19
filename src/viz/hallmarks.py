"""Visualizations for ChEMBL assay coverage and hallmark evidence grade distributions."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_chembl_and_hallmark_coverage(
    coverage_df: pd.DataFrame,
    edges_df: pd.DataFrame,
    output_path: str | Path = "figures/chembl_and_hallmark_coverage.png",
) -> Path:
    """Generate 2-panel figure showing ChEMBL assay coverage and hallmark evidence grades.

    Args:
        coverage_df: DataFrame with compound-level ChEMBL coverage columns (has_binding, has_chembl, etc.).
        edges_df: DataFrame with curated evidence edges (hallmark, grade, etc.).
        output_path: Destination path for figure PNG.

    Returns:
        Path to generated figure.
    """
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    plt.style.use(
        "seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default"
    )
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5), dpi=300)

    # Panel A: ChEMBL Assay Coverage Categories
    total_n = len(coverage_df)
    has_binding = coverage_df["has_binding"].sum()
    functional_only = (coverage_df["has_chembl"] & ~coverage_df["has_binding"]).sum()
    no_chembl = (~coverage_df["has_chembl"]).sum()

    cats = ["≥1 Binding Assay", "Functional / Other Only", "No ChEMBL Activity"]
    vals = [has_binding, functional_only, no_chembl]
    colors = ["#16a34a", "#3b82f6", "#dc2626"]

    bars = ax1.bar(cats, vals, color=colors, width=0.55, edgecolor="black", linewidth=0.8)
    ax1.set_title(
        f"Watchlist ChEMBL Assay Coverage ({has_binding / total_n * 100:.1f}% with Binding)",
        fontsize=11,
        fontweight="bold",
        pad=12,
    )
    ax1.set_ylabel("Number of Compounds", fontsize=10)
    ax1.set_ylim(0, max(vals) * 1.25 if vals else 10)

    for bar in bars:
        height = bar.get_height()
        ax1.annotate(
            f"{height} ({height / total_n * 100:.1f}%)" if total_n > 0 else f"{height}",
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="semibold",
        )

    # Panel B: Hallmark Evidence Grades (Stacked Bar)
    hallmark_order = edges_df["hallmark"].value_counts().index.tolist()
    grades = ["E0", "E1", "E2", "E3", "E4"]
    grade_palette = {
        "E0": "#9ca3af",  # Gray (vendor copy / unstructured)
        "E1": "#93c5fd",  # Light blue (single in vitro)
        "E2": "#3b82f6",  # Blue (target engagement / multi assay)
        "E3": "#f59e0b",  # Amber (in vivo model organism)
        "E4": "#10b981",  # Green (mammalian lifespan / human trial)
    }

    # Group counts
    crosstab = pd.crosstab(edges_df["hallmark"], edges_df["grade"])
    # Reindex to ensure all grades and hallmarks are present
    for g in grades:
        if g not in crosstab.columns:
            crosstab[g] = 0
    crosstab = crosstab.reindex(index=hallmark_order, columns=grades, fill_value=0)

    # Format y-axis labels
    clean_labels = [h.replace("_", " ").title() for h in hallmark_order]

    bottom = pd.Series(0, index=hallmark_order)
    for g in grades:
        ax2.barh(
            clean_labels,
            crosstab[g],
            left=bottom,
            label=f"{g}",
            color=grade_palette[g],
            edgecolor="black",
            linewidth=0.5,
        )
        bottom += crosstab[g]

    ax2.set_title(
        f"Curated Evidence Grades by Hallmark (N = {len(edges_df)} Edges)",
        fontsize=11,
        fontweight="bold",
        pad=12,
    )
    ax2.set_xlabel("Number of Evidence Edges", fontsize=10)
    ax2.legend(title="Evidence Grade", loc="lower right", frameon=True)
    ax2.invert_yaxis()

    fig.suptitle(
        "Geroscience Compound Atlas — Phase 2 ChEMBL Join & Evidence Grading",
        fontsize=13,
        fontweight="bold",
        y=1.02,
    )
    plt.tight_layout()
    plt.savefig(out_file, bbox_inches="tight")
    plt.close(fig)

    return out_file
