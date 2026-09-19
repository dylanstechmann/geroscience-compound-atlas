"""Notebook-free training and evaluation CLI for Geroscience Benchmark v1."""

import argparse
import json
import logging
from pathlib import Path
from typing import Any

import numpy as np

from bench.data import prepare_benchmark_dataset
from bench.models import (
    extract_top_false_positives,
    prepare_feature_matrices,
    train_and_eval_baseline,
    train_and_eval_contender,
)
from bench.splits import generate_and_save_splits
from viz.bench import plot_benchmark_curves, plot_scaffold_size_distribution

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("bench.train")


def run_benchmark_pipeline(
    dataset_parquet: str | Path = "artifacts/benchmark_dataset.parquet",
    splits_json: str | Path = "artifacts/splits.json",
    metrics_json: str | Path = "artifacts/metrics.json",
    roc_figure: str | Path = "figures/benchmark_roc_pr_curves.png",
    scaffold_figure: str | Path = "figures/scaffold_size_distribution.png",
    seeds: list[int] = (42, 123, 456),
    pchembl_threshold: float = 6.0,
) -> dict[str, Any]:
    """Execute complete multi-seed benchmark bake-off comparing Baseline vs Contender under Scaffold vs Random splits."""
    data_file = Path(dataset_parquet)
    if not data_file.exists():
        logger.info("Dataset parquet not found at %s. Generating from ChEMBL...", data_file)
        df = prepare_benchmark_dataset(
            output_parquet=data_file, pchembl_threshold=pchembl_threshold
        )
    else:
        logger.info("Loading existing benchmark dataset from %s...", data_file)
        df = prepare_benchmark_dataset(
            output_parquet=data_file, pchembl_threshold=pchembl_threshold
        )

    # Plot scaffold cluster distribution
    scaffold_path = plot_scaffold_size_distribution(df, output_path=scaffold_figure)
    logger.info("Saved scaffold cluster distribution plot to %s", scaffold_path)

    # Generate and save frozen splits
    logger.info("Generating Bemis-Murcko and random splits across seeds %s...", seeds)
    splits = generate_and_save_splits(df, seeds=seeds, output_path=splits_json)

    # Structure to hold metric results
    results_by_setting: dict[str, dict[str, list[dict[str, float]]]] = {
        "scaffold": {"baseline": [], "contender": []},
        "random": {"baseline": [], "contender": []},
    }

    # For plotting representative curves (seed 42)
    curve_records = []
    representative_false_positives: list[dict[str, Any]] = []

    for split_type in ("scaffold", "random"):
        for seed in seeds:
            split_dict = splits[split_type][str(seed)]
            train_idx = np.array(split_dict["train"])
            val_idx = np.array(split_dict["val"])
            test_idx = np.array(split_dict["test"])

            X_train, y_train, _X_val, _y_val, X_test, y_test = prepare_feature_matrices(
                df, train_idx, val_idx, test_idx
            )

            # 1. Baseline: Logistic Regression
            base_metrics, base_probs = train_and_eval_baseline(
                X_train, y_train, X_test, y_test, seed=seed
            )
            results_by_setting[split_type]["baseline"].append(base_metrics)

            # 2. Contender: HistGradientBoosting
            cont_metrics, cont_probs = train_and_eval_contender(
                X_train, y_train, X_test, y_test, seed=seed
            )
            results_by_setting[split_type]["contender"].append(cont_metrics)

            # Capture seed 42 curves for publication figure
            if seed == 42:
                curve_records.append(
                    {
                        "split_type": split_type,
                        "model_name": "baseline",
                        "y_true": y_test,
                        "y_prob": base_probs,
                        "auroc": base_metrics["auroc"],
                        "auprc": base_metrics["auprc"],
                    }
                )
                curve_records.append(
                    {
                        "split_type": split_type,
                        "model_name": "contender",
                        "y_true": y_test,
                        "y_prob": cont_probs,
                        "auroc": cont_metrics["auroc"],
                        "auprc": cont_metrics["auprc"],
                    }
                )

                if split_type == "scaffold":
                    representative_false_positives = extract_top_false_positives(
                        df, test_idx, y_test, cont_probs, top_n=10
                    )

    # Plot ROC & PR curves
    curves_path = plot_benchmark_curves(curve_records, output_path=roc_figure)
    logger.info("Saved ROC and PR curves to %s", curves_path)

    # Aggregate metrics across seeds
    aggregated_summary: dict[str, Any] = {}
    metric_keys = ["auroc", "auprc", "recall_at_fpr_0.05", "brier_score"]

    for split_type in ("scaffold", "random"):
        aggregated_summary[split_type] = {}
        for model in ("baseline", "contender"):
            raw_list = results_by_setting[split_type][model]
            model_summary = {}
            for k in metric_keys:
                vals = [m[k] for m in raw_list]
                model_summary[f"{k}_mean"] = round(float(np.mean(vals)), 4)
                model_summary[f"{k}_std"] = round(float(np.std(vals)), 4)
            aggregated_summary[split_type][model] = model_summary

    # Calculate headline leakage gaps
    leakage_gaps = {
        "baseline_auroc_gap": round(
            aggregated_summary["random"]["baseline"]["auroc_mean"]
            - aggregated_summary["scaffold"]["baseline"]["auroc_mean"],
            4,
        ),
        "contender_auroc_gap": round(
            aggregated_summary["random"]["contender"]["auroc_mean"]
            - aggregated_summary["scaffold"]["contender"]["auroc_mean"],
            4,
        ),
        "baseline_auprc_gap": round(
            aggregated_summary["random"]["baseline"]["auprc_mean"]
            - aggregated_summary["scaffold"]["baseline"]["auprc_mean"],
            4,
        ),
        "contender_auprc_gap": round(
            aggregated_summary["random"]["contender"]["auprc_mean"]
            - aggregated_summary["scaffold"]["contender"]["auprc_mean"],
            4,
        ),
    }

    metrics_payload = {
        "task": "chembl_mtor_activity_classification",
        "label_definition": "active = 1 if pchembl_value >= 6.0 else 0",
        "dataset_size": len(df),
        "num_active": int(df["active"].sum()),
        "num_inactive": int(len(df) - df["active"].sum()),
        "num_scaffolds": int(df["murcko_scaffold"].nunique()),
        "seeds": list(seeds),
        "summary": aggregated_summary,
        "leakage_gaps": leakage_gaps,
        "raw_results": results_by_setting,
        "top_false_positives": representative_false_positives,
    }

    out_json = Path(metrics_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(metrics_payload, f, indent=2)

    logger.info("Saved complete metrics artifact to %s", out_json)

    # Console summary table
    logger.info("=" * 70)
    logger.info("BENCHMARK v1 RESULTS SUMMARY (Mean ± Std over 3 seeds):")
    logger.info("-" * 70)
    logger.info("SCAFFOLD SPLIT (Honest Generalization):")
    logger.info(
        "  Baseline (LogisticReg): AUROC = %.4f ± %.4f | AUPRC = %.4f ± %.4f",
        aggregated_summary["scaffold"]["baseline"]["auroc_mean"],
        aggregated_summary["scaffold"]["baseline"]["auroc_std"],
        aggregated_summary["scaffold"]["baseline"]["auprc_mean"],
        aggregated_summary["scaffold"]["baseline"]["auprc_std"],
    )
    logger.info(
        "  Contender (HistGradB):  AUROC = %.4f ± %.4f | AUPRC = %.4f ± %.4f",
        aggregated_summary["scaffold"]["contender"]["auroc_mean"],
        aggregated_summary["scaffold"]["contender"]["auroc_std"],
        aggregated_summary["scaffold"]["contender"]["auprc_mean"],
        aggregated_summary["scaffold"]["contender"]["auprc_std"],
    )
    logger.info("-" * 70)
    logger.info("RANDOM SPLIT (Chemical Cousin Leakage):")
    logger.info(
        "  Baseline (LogisticReg): AUROC = %.4f ± %.4f | AUPRC = %.4f ± %.4f",
        aggregated_summary["random"]["baseline"]["auroc_mean"],
        aggregated_summary["random"]["baseline"]["auroc_std"],
        aggregated_summary["random"]["baseline"]["auprc_mean"],
        aggregated_summary["random"]["baseline"]["auprc_std"],
    )
    logger.info(
        "  Contender (HistGradB):  AUROC = %.4f ± %.4f | AUPRC = %.4f ± %.4f",
        aggregated_summary["random"]["contender"]["auroc_mean"],
        aggregated_summary["random"]["contender"]["auroc_std"],
        aggregated_summary["random"]["contender"]["auprc_mean"],
        aggregated_summary["random"]["contender"]["auprc_std"],
    )
    logger.info("-" * 70)
    logger.info("LEAKAGE GAP (Random AUROC - Scaffold AUROC):")
    logger.info("  Baseline Gap:  +%.4f AUROC", leakage_gaps["baseline_auroc_gap"])
    logger.info("  Contender Gap: +%.4f AUROC", leakage_gaps["contender_auroc_gap"])
    logger.info("=" * 70)

    return metrics_payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Execute Geroscience Benchmark Bake-Off")
    parser.add_argument("--parquet", default="artifacts/benchmark_dataset.parquet")
    parser.add_argument("--splits", default="artifacts/splits.json")
    parser.add_argument("--metrics", default="artifacts/metrics.json")
    parser.add_argument("--roc-fig", default="figures/benchmark_roc_pr_curves.png")
    parser.add_argument("--scaffold-fig", default="figures/scaffold_size_distribution.png")
    args = parser.parse_args()

    run_benchmark_pipeline(
        dataset_parquet=args.parquet,
        splits_json=args.splits,
        metrics_json=args.metrics,
        roc_figure=args.roc_fig,
        scaffold_figure=args.scaffold_fig,
    )
