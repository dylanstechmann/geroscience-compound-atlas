#!/usr/bin/env python3
"""Continuous pChEMBL check on the existing, frozen mTOR scaffold splits."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

from atlas.features import compute_morgan_fingerprint
from bench.models import prepare_feature_matrices


def checked_scaffold_indices(
    df: pd.DataFrame, split: dict[str, list[int]]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fail if a stored split is incomplete, duplicated, or scaffold-leaking."""
    indices = tuple(np.asarray(split[name], dtype=int) for name in ("train", "val", "test"))
    combined = np.concatenate(indices)
    if len(combined) != len(df) or set(combined.tolist()) != set(range(len(df))):
        raise ValueError("Frozen split must partition every dataset row exactly once")
    scaffold_sets = [set(df.iloc[idx]["murcko_scaffold"]) for idx in indices]
    if any(scaffold_sets[i] & scaffold_sets[j] for i, j in ((0, 1), (0, 2), (1, 2))):
        raise ValueError("Frozen split has overlapping Murcko scaffolds")
    return indices


def regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    return {
        "mae": round(float(mean_absolute_error(y_true, y_pred)), 4),
        "rmse": round(float(np.sqrt(mean_squared_error(y_true, y_pred))), 4),
        "r2": round(float(r2_score(y_true, y_pred)), 4),
    }


def run_regression(
    dataset_parquet: str | Path = "artifacts/benchmark_dataset.parquet",
    splits_json: str | Path = "artifacts/splits.json",
    metrics_json: str | Path = "artifacts/regression_metrics.json",
) -> dict:
    """Fit fixed Ridge and train-mean baselines; evaluate untouched scaffold tests."""
    df = pd.read_parquet(dataset_parquet)
    required = {"canonical_smiles", "murcko_scaffold", "pchembl_value", "active"}
    if not required.issubset(df.columns):
        raise ValueError(f"Dataset is missing columns: {sorted(required - set(df.columns))}")
    y = pd.to_numeric(df["pchembl_value"], errors="raise").to_numpy(dtype=float)
    if not np.isfinite(y).all():
        raise ValueError("pChEMBL measurements must be finite")

    fingerprints = [compute_morgan_fingerprint(s) for s in df["canonical_smiles"]]
    if any(fp is None for fp in fingerprints):
        raise ValueError("Dataset contains an invalid SMILES")
    df["fingerprint"] = fingerprints

    frozen = json.loads(Path(splits_json).read_text(encoding="utf-8"))
    if "scaffold" not in frozen or not frozen["scaffold"]:
        raise ValueError("Frozen scaffold splits are required")

    results = {}
    for seed, split in sorted(frozen["scaffold"].items(), key=lambda item: int(item[0])):
        train_idx, val_idx, test_idx = checked_scaffold_indices(df, split)
        X_train, _, _, _, X_test, _ = prepare_feature_matrices(df, train_idx, val_idx, test_idx)
        model = Ridge(alpha=1.0)
        model.fit(X_train, y[train_idx])
        baseline = np.full(len(test_idx), y[train_idx].mean())
        results[seed] = {
            "n_train": len(train_idx),
            "n_val": len(val_idx),
            "n_test": len(test_idx),
            "train_pchembl_mean": round(float(y[train_idx].mean()), 4),
            "test_pchembl_mean": round(float(y[test_idx].mean()), 4),
            "ridge": regression_metrics(y[test_idx], model.predict(X_test)),
            "train_mean": regression_metrics(y[test_idx], baseline),
        }

    summary = {
        model: {
            metric: round(float(np.mean([result[model][metric] for result in results.values()])), 4)
            for metric in ("mae", "rmse", "r2")
        }
        for model in ("ridge", "train_mean")
    }
    payload = {
        "task": "chembl_mtor_pchembl_regression",
        "split": "frozen_bemis_murcko_scaffold",
        "model": "Ridge(alpha=1.0); Morgan radius 2, 2048 bits plus train-scaled descriptors",
        "note": "Computational prediction of recorded ChEMBL values, not measured activity for generated cards or evidence of rejuvenation.",
        "n_molecules": len(df),
        "per_seed": results,
        "mean_over_seeds": summary,
    }
    output = Path(metrics_json)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", default="artifacts/benchmark_dataset.parquet")
    parser.add_argument("--splits", default="artifacts/splits.json")
    parser.add_argument("--metrics", default="artifacts/regression_metrics.json")
    args = parser.parse_args()
    print(
        json.dumps(
            run_regression(args.dataset, args.splits, args.metrics)["mean_over_seeds"], indent=2
        )
    )
