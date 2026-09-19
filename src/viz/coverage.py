"""Visualizations for identifier resolution coverage, failure diagnostics, and modalities."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_resolution_coverage(
    df: pd.DataFrame, output_path: str | Path = "figures/resolution_coverage.png"
) -> Path:
    """Generate diagnostic coverage figure showing resolved vs failed names and modality distribution.

    Args:
        df: DataFrame containing resolution pipeline outputs.
        output_path: Path where output PNG figure will be saved.

    Returns:
        Path to the saved figure.
    """
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    # Styling settings
    plt.style.use(
        "seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default"
    )
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5), dpi=300)

    # 1. Resolution Status Breakdown
    res_counts = df["resolved"].value_counts()
    resolved_n = res_counts.get(True, 0)
    failed_n = res_counts.get(False, 0)
    total_n = len(df)
    res_pct = (resolved_n / total_n * 100) if total_n > 0 else 0

    bar_labels = ["Resolved", "Unresolved / Failed"]
    bar_values = [resolved_n, failed_n]
    bar_colors = ["#2563eb", "#dc2626"]

    bars = ax1.bar(
        bar_labels, bar_values, color=bar_colors, width=0.5, edgecolor="black", linewidth=0.8
    )
    ax1.set_title(
        f"PubChem Name Resolution ({res_pct:.1f}% Resolved)", fontsize=12, fontweight="bold", pad=12
    )
    ax1.set_ylabel("Number of Compounds", fontsize=10)
    ax1.set_ylim(0, max(bar_values) * 1.25 if bar_values else 10)

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

    # 2. Modality Breakdown
    mod_counts = df["modality"].value_counts()
    mod_labels = [m.replace("_", " ").title() for m in mod_counts.index]
    mod_colors = ["#0284c7", "#10b981", "#f59e0b", "#6b7280"][: len(mod_counts)]

    _wedges, _texts, autotexts = ax2.pie(
        mod_counts.values,
        labels=mod_labels,
        autopct="%1.1f%%",
        startangle=140,
        colors=mod_colors,
        wedgeprops={"edgecolor": "white", "linewidth": 1.5},
    )
    for at in autotexts:
        at.set_color("white")
        at.set_fontsize(9)
        at.set_weight("bold")

    ax2.set_title("Modality Distribution", fontsize=12, fontweight="bold", pad=12)

    fig.suptitle(
        "Geroscience Compound Atlas — Phase 1 Coverage Diagnostic",
        fontsize=14,
        fontweight="bold",
        y=1.02,
    )
    plt.tight_layout()
    plt.savefig(out_file, bbox_inches="tight")
    plt.close(fig)

    return out_file
