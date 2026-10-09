"""Render the explicitly synthetic sensitivity analysis with matplotlib."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from analysis_models import free_guest_fraction

ROOT = Path(__file__).resolve().parent


def main():
    with (ROOT / "synthetic_exposure_sensitivity.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    labels = ["Reference", "10x oral F", "4x half-life", "4x half-life,\nmore binding /\nless cell access", "10x oral F,\nless activity"]
    values = [float(row["exposure_proxy_relative_to_reference"]) for row in rows]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.7), layout="constrained")
    bars = axes[0].bar(range(len(values)), values, color=["#64748b", "#0d9488", "#2563eb", "#dc2626", "#d97706"])
    axes[0].bar_label(bars, labels=[f"{value:g}x" for value in values], padding=4)
    axes[0].set_xticks(range(len(values)), labels, fontsize=9)
    axes[0].set_ylim(0, 12)
    axes[0].set_ylabel("Effective-exposure proxy relative to reference")
    axes[0].set_title("Longer plasma residence can lose useful exposure", loc="left", fontsize=11)
    axes[0].axhline(1, color="#64748b", linewidth=1, linestyle="--")
    hosts = [i/100 for i in range(201)]
    for kd, color in ((1, "#64748b"), (0.01, "#2563eb"), (0.000001, "#0d9488")):
        axes[1].plot(hosts, [free_guest_fraction(host, kd) for host in hosts],
                     label=f"Kd / guest total = {kd:g}", color=color, linewidth=2)
    axes[1].plot(hosts, [max(0, 1-host) for host in hosts], color="#dc2626", linestyle="--", label="Infinite-affinity capacity limit")
    axes[1].set_xlabel("Binding-site total / guest total (dimensionless)")
    axes[1].set_ylabel("Equilibrium free guest fraction")
    axes[1].set_title("Tighter binding cannot replace binding-site capacity", loc="left", fontsize=11)
    axes[1].set_ylim(-0.03, 1.05)
    axes[1].legend(fontsize=8, loc="upper right")
    fig.suptitle("Synthetic examples: delivery and sequestration tradeoffs", fontsize=15, fontweight="bold")
    fig.supxlabel("Illustrative mathematics only. No human PK, caffeine-binding affinity, dosing or sleep outcome is predicted.", fontsize=9)
    fig.savefig(ROOT / "synthetic_sensitivity.png", dpi=170, facecolor="white")
    plt.close(fig)
    print(ROOT / "synthetic_sensitivity.png")


if __name__ == "__main__":
    main()
