"""Standalone scientific figure of synthetic selectivity counterexamples."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from selectivity_model import ROOT, Tissue, aggregate_responses, matched_target


def main():
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    figure, axes = plt.subplots(1, 2, figsize=(12, 5.6), constrained_layout=True)
    target = Tissue("target", 1, 1, 1)
    demand = [index / 100 for index in range(1, 96)]
    for ratio, color in ((1, "#9c3d2e"), (0.1, "#217a80"), (0.01, "#4159a5")):
        other = Tissue("other", ratio, 1, 1)
        response = [100 * matched_target(value, target, other)["other_tissue_inhibition"]
                    for value in demand]
        axes[0].plot([100 * value for value in demand], response, color=color,
                     label=f"Assumed activity ratio = {ratio:g}", linewidth=2.4)
    axes[0].set(xlabel="Assumed target pathway inhibition (%)",
                ylabel="Other-tissue pathway inhibition (%)", ylim=(0, 100), xlim=(0, 100),
                title="Selectivity must hold at matched target response")
    axes[0].legend(loc="upper left", frameon=False)
    axes[0].grid(alpha=0.15)

    result = aggregate_responses(0.5, target,
        [(Tissue("basal", 1, 1, 1), 0.1), (Tissue("adluminal", 0.01, 1, 1), 0.9)])
    values = [100 * row["other_tissue_inhibition"] for row in result["compartments"]]
    values.append(100 * result["weighted_mean_inhibition"])
    axes[1].bar(["Basal proxy", "Adluminal proxy", "Weighted mean"], values,
                color=["#9c3d2e", "#4159a5", "#73808f"], width=0.65)
    for index, value in enumerate(values):
        axes[1].text(index, value + 1.3, f"{value:.2f}%", ha="center")
    axes[1].set(ylabel="Synthetic pathway inhibition (%)", ylim=(0, 65),
                title="An average can hide a compartment's response")
    axes[1].text(0.5, 0.87, "Assumed weights: 10% basal / 90% adluminal\n"
                 "Target response matched at 50%", transform=axes[1].transAxes,
                 ha="center", va="top", fontsize=9)
    axes[1].grid(axis="y", alpha=0.15)
    figure.suptitle("Synthetic tissue selectivity: no measured rapalog or human parameters", fontsize=13)
    figure.supxlabel("Hill-1 equilibrium examples; no prediction of fertility, aging benefit, oral PK or safe exposure.",
                     fontsize=9)
    figure.savefig(ROOT / "synthetic_selectivity.png", dpi=170)
    plt.close(figure)


if __name__ == "__main__":
    main()
