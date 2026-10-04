"""Reliability diagnostics on the existing frozen mTOR scaffold splits."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import platform
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from rdkit import Chem, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold

from atlas.features import compute_morgan_fingerprint
from bench.config import BenchConfig, load_bench_config
from bench.data import activity_labels
from bench.models import (
    prepare_feature_matrices,
    train_and_eval_baseline,
    train_and_eval_contender,
)
from bench.neighborhood import (
    maximum_training_tanimoto,
    neighborhood_calibration_diagnostics,
    validated_fingerprints,
    validated_row_indexes,
)


def _validate_scaffold_partition(frame: pd.DataFrame, split: dict) -> dict[str, np.ndarray]:
    if not isinstance(split, dict) or any(key not in split for key in ("train", "val", "test")):
        raise ValueError("each frozen scaffold split needs train, val and test indexes")
    parts = {
        key: validated_row_indexes(split[key], len(frame), key) for key in ("train", "val", "test")
    }
    all_rows = np.concatenate(list(parts.values()))
    if len(all_rows) != len(frame) or len(np.unique(all_rows)) != len(frame):
        raise ValueError(
            "frozen partitions must be disjoint and cover every dataset row exactly once"
        )
    scaffold_sets = {key: set(frame.iloc[rows]["murcko_scaffold"]) for key, rows in parts.items()}
    for left, right in (("train", "val"), ("train", "test"), ("val", "test")):
        if scaffold_sets[left].intersection(scaffold_sets[right]):
            raise ValueError(f"frozen split has scaffold leakage between {left} and {right}")
    return parts


def evaluate_frozen_scaffold_neighborhoods(
    frame: pd.DataFrame,
    split_document: dict,
    *,
    seeds: tuple[int, ...] | None = None,
    config: BenchConfig | None = None,
) -> dict:
    """Re-fit fixed model settings on supplied scaffold training partitions.

    Only training rows fit descriptor imputation, scaling and models. Test rows
    evaluate descriptive reliability; no threshold or calibrator is selected.
    """
    config = config or BenchConfig()
    if (
        config.model.baseline != "logistic_regression"
        or config.model.contender != "hist_gradient_boosting"
    ):
        raise ValueError("this diagnostic supports only the frozen logistic and HGB model families")
    if any(column not in frame for column in ("fingerprint", "active", "murcko_scaffold")):
        raise ValueError("benchmark frame requires fingerprint, active and murcko_scaffold columns")
    if frame["murcko_scaffold"].isna().any() or any(
        not isinstance(value, str) for value in frame["murcko_scaffold"]
    ):
        raise ValueError("scaffold annotations must be strings (empty string denotes acyclic)")
    labels = frame["active"].to_numpy()
    if not np.isin(labels, [0, 1]).all():
        raise ValueError("active labels must be exactly zero or one")
    if "pchembl_value" in frame and not np.array_equal(
        labels, activity_labels(frame["pchembl_value"], config.thresholds.pchembl_active)
    ):
        raise ValueError("active labels do not match the configured pChEMBL cutoff")
    scaffold_splits = split_document.get("scaffold") if isinstance(split_document, dict) else None
    if not isinstance(scaffold_splits, dict) or not scaffold_splits:
        raise ValueError("frozen scaffold splits are required")
    if any(
        not isinstance(seed, str) or not seed.isdigit() or str(int(seed)) != seed
        for seed in scaffold_splits
    ):
        raise ValueError("frozen split seed keys must be canonical nonnegative integers")
    available_seeds = tuple(sorted(int(seed) for seed in scaffold_splits))
    selected_seeds = available_seeds if seeds is None else tuple(seeds)
    if (
        not selected_seeds
        or len(selected_seeds) != len(set(selected_seeds))
        or any(
            isinstance(seed, bool)
            or not isinstance(seed, (int, np.integer))
            or str(seed) not in scaffold_splits
            for seed in selected_seeds
        )
    ):
        raise ValueError(
            "requested unique integer seeds must exist in the frozen scaffold split file"
        )

    fingerprint_matrix = validated_fingerprints(np.asarray(frame["fingerprint"].tolist()))
    if len(fingerprint_matrix) != len(frame):
        raise ValueError("fingerprints must have one fixed-width vector per dataset row")
    by_seed = {}
    for seed in selected_seeds:
        parts = _validate_scaffold_partition(frame, scaffold_splits[str(seed)])
        train_idx, val_idx, test_idx = (parts[key] for key in ("train", "val", "test"))
        if len(np.unique(labels[train_idx])) != 2 or len(np.unique(labels[test_idx])) != 2:
            raise ValueError("train and test partitions must each contain both activity classes")
        X_train, y_train, _X_val, _y_val, X_test, y_test = prepare_feature_matrices(
            frame,
            train_idx,
            val_idx,
            test_idx,
            include_descriptors=config.features.include_descriptors,
        )
        nearest = maximum_training_tanimoto(fingerprint_matrix, train_idx, test_idx)
        split_results = {}
        for model_name, fit_predict, settings in (
            ("baseline", train_and_eval_baseline, {}),
            ("contender", train_and_eval_contender, asdict(config.model.hgb_params)),
        ):
            metrics, probabilities = fit_predict(
                X_train, y_train, X_test, y_test, seed=seed, **settings
            )
            split_results[model_name] = {
                "metrics": metrics,
                "reliability_by_nearest_training_similarity": neighborhood_calibration_diagnostics(
                    y_test, probabilities, nearest
                ),
                "test_predictions": [
                    {
                        "dataset_row": int(row),
                        "active": int(label),
                        "probability": float(probability),
                        "nearest_training_tanimoto": float(similarity),
                    }
                    for row, label, probability, similarity in zip(
                        test_idx, y_test, probabilities, nearest
                    )
                ],
            }
        by_seed[str(seed)] = {
            "n_train": len(train_idx),
            "n_validation": len(val_idx),
            "n_test": len(test_idx),
            "models": split_results,
        }
    return {
        "schema_version": 1,
        "task": "chembl_mtor_activity_classification",
        "split": "existing_frozen_bemis_murcko_scaffold_partitions",
        "status": "internal_split_diagnostic_not_external_validation_or_calibration",
        "fit_calibrator": False,
        "decision_threshold_changed": False,
        "split_file_unchanged": True,
        "pchembl_threshold": config.thresholds.pchembl_active,
        "model_settings": {
            "baseline": "LogisticRegression(C=1, max_iter=1000, solver=lbfgs)",
            "contender": asdict(config.model.hgb_params),
            "include_descriptors": config.features.include_descriptors,
        },
        "fingerprint_width": int(fingerprint_matrix.shape[1]),
        "similarity_metric": "maximum binary-fingerprint Tanimoto to the corresponding training partition",
        "interpretation": "Describes probability error by chemical-neighborhood proximity on fixed internal test predictions; does not establish external calibration or applicability domain.",
        "per_seed": by_seed,
    }


def run(
    dataset_path: Path, splits_path: Path, output_path: Path, *, config_path: Path | None = None
) -> dict:
    # A diagnostic must never overwrite one of its frozen input artifacts.
    resolved_config = config_path or Path(__file__).resolve().parents[2] / "configs" / "bench.yaml"
    protected = (
        dataset_path,
        splits_path,
        dataset_path.with_suffix(".fingerprints.npy"),
        resolved_config,
    )
    if output_path.resolve() in {path.resolve() for path in protected}:
        raise ValueError("report output must not overwrite a frozen input or configuration")
    # Snapshot bytes once so hashes refer to exactly the data read and scored.
    dataset_bytes = dataset_path.read_bytes()
    split_bytes = splits_path.read_bytes()
    frame = pd.read_parquet(io.BytesIO(dataset_bytes))
    if "canonical_smiles" not in frame or "murcko_scaffold" not in frame:
        raise ValueError("frozen benchmark requires structure and scaffold columns")
    config = load_bench_config(config_path)
    if (
        config.features.fingerprint != "morgan"
        or config.features.radius != 2
        or config.features.n_bits != 2048
    ):
        raise ValueError(
            "this diagnostic requires the frozen Morgan radius-2 2048-bit feature schema"
        )
    reconstructed = []
    for row, (smiles, saved_scaffold) in enumerate(
        zip(frame["canonical_smiles"], frame["murcko_scaffold"])
    ):
        mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
        if mol is None:
            raise ValueError(f"benchmark has invalid structure at row {row}")
        if MurckoScaffold.MurckoScaffoldSmiles(mol=mol) != saved_scaffold:
            raise ValueError(f"benchmark scaffold does not match its structure at row {row}")
        reconstructed.append(compute_morgan_fingerprint(mol, radius=2, n_bits=2048))
    expected_fingerprints = np.asarray(reconstructed, dtype=np.uint8)
    fingerprint_path = dataset_path.with_suffix(".fingerprints.npy")
    fingerprint_bytes = None
    if fingerprint_path.is_file():
        fingerprint_bytes = fingerprint_path.read_bytes()
        fingerprints = validated_fingerprints(
            np.load(io.BytesIO(fingerprint_bytes), allow_pickle=False)
        )
        if fingerprints.shape != expected_fingerprints.shape or not np.array_equal(
            fingerprints, expected_fingerprints
        ):
            raise ValueError(
                "fingerprint sidecar does not match the benchmark structures and row order"
            )
        fingerprint_source = "sidecar_verified_against_frozen_smiles"
    else:
        fingerprints = expected_fingerprints
        fingerprint_source = "reconstructed_from_frozen_smiles_no_cache_written"
    frame = frame.copy()
    frame["fingerprint"] = list(fingerprints)
    report = evaluate_frozen_scaffold_neighborhoods(frame, json.loads(split_bytes), config=config)
    report["inputs"] = {
        "benchmark_parquet": str(dataset_path),
        "benchmark_sha256": hashlib.sha256(dataset_bytes).hexdigest(),
        "split_json": str(splits_path),
        "split_sha256": hashlib.sha256(split_bytes).hexdigest(),
        "fingerprint_source": fingerprint_source,
        "fingerprint_sidecar": str(fingerprint_path) if fingerprint_bytes is not None else None,
        "fingerprint_sidecar_sha256": hashlib.sha256(fingerprint_bytes).hexdigest()
        if fingerprint_bytes is not None
        else None,
        "ordered_binary_fingerprints_sha256": hashlib.sha256(
            expected_fingerprints.tobytes()
        ).hexdigest(),
        "configuration": asdict(config),
    }
    report["runtime"] = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scikit_learn": sklearn.__version__,
        "rdkit": rdBase.rdkitVersion,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, default=Path("artifacts/benchmark_dataset.parquet"))
    parser.add_argument("--splits", type=Path, default=Path("artifacts/splits.json"))
    parser.add_argument("--out", type=Path, default=Path("artifacts/neighborhood_reliability.json"))
    parser.add_argument("--config", type=Path, default=None)
    args = parser.parse_args(argv)
    result = run(args.dataset, args.splits, args.out, config_path=args.config)
    print(
        json.dumps(
            {"out": str(args.out), "seeds": list(result["per_seed"]), "fit_calibrator": False},
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
