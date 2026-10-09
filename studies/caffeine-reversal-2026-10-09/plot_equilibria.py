"""Plot source-equilibrium reconstructions at arbitrary assay concentrations."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent


def main():
    with (ROOT / "illustrative_equilibria.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    labels = {
        "receptor2_published_NMR": "Receptor 2: full NMR species",
        "receptor6_published_NMR": "Receptor 6: full NMR species",
        "receptor2_11_species_only": "Receptor 2: 1:1 species only (incomplete)",
        "incorrect_BC500_as_Kd": "BC500 used as Kd (incorrect)",
    }
    colors = ("#176b87", "#6b8e23", "#777777", "#b44444")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.8), layout="constrained")
    for (model_id, label), color in zip(labels.items(), colors, strict=True):
        selected = [r for r in rows if r["model_id"] == model_id and float(r["illustrative_guest_total_M"]) == 1e-5]
        axes[0].plot(
            [float(r["host_to_guest_total_ratio"]) for r in selected],
            [100 * float(r["guest_not_host_captured_fraction"]) for r in selected],
            label=label, color=color, marker="o", linestyle="--" if "incorrect" in model_id or "only" in model_id else "-",
        )
        selected = [r for r in rows if r["model_id"] == model_id and float(r["host_to_guest_total_ratio"]) == 1]
        axes[1].plot(
            [1e6 * float(r["illustrative_guest_total_M"]) for r in selected],
            [100 * float(r["guest_not_host_captured_fraction"]) for r in selected],
            color=color, marker="o", linestyle="--" if "incorrect" in model_id or "only" in model_id else "-",
        )
    axes[0].set_title("Chosen guest total: 10 µM")
    axes[0].set_xlabel("Total host / total guest (molar ratio)")
    axes[0].legend(loc="upper right", fontsize=8.5)
    axes[1].set_title("Equal total host and guest")
    axes[1].set_xlabel("Chosen guest total in aqueous model (µM)")
    axes[1].set_xscale("log")
    axes[1].set_xticks([1, 10, 100], labels=["1", "10", "100"])
    for ax in axes:
        ax.set_ylabel("Guest not captured by host (%)")
        ax.set_ylim(0, 105)
        ax.grid(alpha=0.2)
    fig.suptitle(
        "Caffeine binding: published water equilibria ≠ an oral antidote\n"
        "NMR central constants; illustrative totals; no plasma, kinetics or sleep model",
        fontsize=13,
    )
    fig.savefig(ROOT / "aqueous_equilibria.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
