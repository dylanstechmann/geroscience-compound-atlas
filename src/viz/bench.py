"""Visualization functions for benchmark curves, scaffold cluster sizes, and leakage gaps."""

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve, roc_curve


def plot_benchmark_curves(
    results_records: list[dict[str, Any]],
    output_path: str | Path = "figures/benchmark_roc_pr_curves.png",
) -> Path:
    """Generate 2-panel figure showing ROC and PR curves for Scaffold vs Random splits across models.

    Args:
        results_records: List of dictionaries containing split_type, model_name, y_true, y_prob, etc.
        output_path: Filepath for output PNG.

    Returns:
        Path to generated figure.
    """
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    plt.style.use(
        "seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default"
    )
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5), dpi=300)

    # Styling colors
    curve_styles = {
        ("scaffold", "baseline"): {
            "color": "#3b82f6",
            "linestyle": "-",
            "label": "Scaffold: Baseline (Logistic Reg)",
        },
        ("scaffold", "contender"): {
            "color": "#1d4ed8",
            "linestyle": "--",
            "label": "Scaffold: Contender (HistGradBoost)",
        },
        ("random", "baseline"): {
            "color": "#f97316",
            "linestyle": ":",
            "label": "Random Leakage: Baseline",
        },
        ("random", "contender"): {
            "color": "#dc2626",
            "linestyle": "-.",
            "label": "Random Leakage: Contender",
        },
    }

    # Reference baseline lines
    ax1.plot([0, 1], [0, 1], "k--", alpha=0.4, label="Chance (AUC = 0.50)")

    for rec in results_records:
        split = rec["split_type"]
        model = rec["model_name"]
        y_true = rec["y_true"]
        y_prob = rec["y_prob"]
        auroc = rec["auroc"]
        auprc = rec["auprc"]

        style = curve_styles.get(
            (split, model), {"color": "black", "linestyle": "-", "label": f"{split}-{model}"}
        )

        # ROC Curve
        fpr, tpr, _ = roc_curve(y_true, y_prob)
        ax1.plot(
            fpr,
            tpr,
            color=style["color"],
            linestyle=style["linestyle"],
            linewidth=1.8,
            label=f"{style['label']} (AUROC = {auroc:.3f})",
        )

        # Precision-Recall Curve
        prec, rec_pts, _ = precision_recall_curve(y_true, y_prob)
        ax2.plot(
            rec_pts,
            prec,
            color=style["color"],
            linestyle=style["linestyle"],
            linewidth=1.8,
            label=f"{style['label']} (AUPRC = {auprc:.3f})",
        )

    # Base prevalence line for PR curve
    if results_records:
        prevalence = float(np.mean(results_records[0]["y_true"]))
        ax2.axhline(
            prevalence, color="k", linestyle="--", alpha=0.4, label=f"Prevalence ({prevalence:.2f})"
        )

    ax1.set_title("Receiver Operating Characteristic (ROC)", fontsize=11, fontweight="bold")
    ax1.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=10)
    ax1.set_ylabel("True Positive Rate (Sensitivity)", fontsize=10)
    ax1.set_xlim([0.0, 1.0])
    ax1.set_ylim([0.0, 1.05])
    ax1.legend(loc="lower right", fontsize=8, frameon=True)

    ax2.set_title("Precision-Recall (PR) Curve", fontsize=11, fontweight="bold")
    ax2.set_xlabel("Recall", fontsize=10)
    ax2.set_ylabel("Precision (PPV)", fontsize=10)
    ax2.set_xlim([0.0, 1.0])
    ax2.set_ylim([0.0, 1.05])
    ax2.legend(loc="lower left", fontsize=8, frameon=True)

    fig.suptitle(
        "Geroscience Benchmark v1 — Scaffold Split vs Random Leakage Diagnostic",
        fontsize=13,
        fontweight="bold",
        y=1.02,
    )
    plt.tight_layout()
    plt.savefig(out_file, bbox_inches="tight")
    plt.close(fig)

    return out_file


def plot_scaffold_size_distribution(
    df: pd.DataFrame,
    output_path: str | Path = "figures/scaffold_size_distribution.png",
    scaffold_col: str = "murcko_scaffold",
) -> Path:
    """Generate histogram of Bemis-Murcko scaffold cluster sizes."""
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    counts = df[scaffold_col].value_counts().values

    plt.style.use(
        "seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default"
    )
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=300)

    # Histogram of cluster sizes
    bins = np.arange(1, min(counts.max() + 2, 30)) - 0.5
    ax.hist(counts, bins=bins, color="#2563eb", edgecolor="black", linewidth=0.8, rwidth=0.8)

    ax.set_title(
        f"Bemis-Murcko Scaffold Cluster Size Distribution (N = {len(df)} molecules, {df[scaffold_col].nunique()} scaffolds)",
        fontsize=11,
        fontweight="bold",
        pad=12,
    )
    ax.set_xlabel("Number of Molecules Sharing Scaffold", fontsize=10)
    ax.set_ylabel("Number of Scaffolds", fontsize=10)
    ax.set_yscale("log")

    singletons = int(np.sum(counts == 1))
    total_scaffolds = len(counts)
    ax.annotate(
        f"Singletons: {singletons}/{total_scaffolds} ({singletons / total_scaffolds * 100:.1f}%)",
        xy=(0.65, 0.85),
        xycoords="axes fraction",
        fontsize=10,
        fontweight="semibold",
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "white", "edgecolor": "#9ca3af"},
    )

    plt.tight_layout()
    plt.savefig(out_file, bbox_inches="tight")
    plt.close(fig)

    return out_file
