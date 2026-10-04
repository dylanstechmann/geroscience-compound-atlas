"""Internal-split diagnostics stratified by chemical-neighborhood similarity."""

from __future__ import annotations

import numpy as np


def validated_row_indexes(indexes, n_rows: int, name: str, *, nonempty: bool = True) -> np.ndarray:
    """Validate positional indexes before conversion; never truncate floats."""
    values = np.asarray(indexes)
    if (
        values.ndim != 1
        or (nonempty and not len(values))
        or (len(values) and values.dtype.kind not in "iu")
        or len(values) != len(set(values.tolist()))
    ):
        raise ValueError(f"{name} must contain unique integer row indexes")
    if len(values) and (np.any(values < 0) or np.any(values >= n_rows)):
        raise ValueError(f"{name} fingerprint index is out of range")
    return values.astype(np.int64)


def validated_fingerprints(fingerprints: np.ndarray) -> np.ndarray:
    values = np.asarray(fingerprints)
    if (
        values.ndim != 2
        or not values.shape[0]
        or not values.shape[1]
        or values.dtype.kind not in "biuf"
        or not np.isfinite(values).all()
        or not np.isin(values, [0, 1]).all()
    ):
        raise ValueError("fingerprints must be a nonempty finite binary matrix")
    return values


def maximum_training_tanimoto(fingerprints: np.ndarray, train_idx, test_idx) -> np.ndarray:
    """Maximum bit-vector Tanimoto to train rows, with no test-reference leakage.

    Empty/empty vectors use similarity zero. Matrix products avoid allocating a
    test-by-train-by-bit cube, which becomes excessive for 2048-bit fingerprints.
    """
    values = validated_fingerprints(fingerprints)
    if np.asarray(train_idx).size == 0 or np.asarray(test_idx).size == 0:
        raise ValueError("fingerprints require nonempty train and test indexes")
    train_idx = validated_row_indexes(train_idx, len(values), "train")
    test_idx = validated_row_indexes(test_idx, len(values), "test")
    if np.intersect1d(train_idx, test_idx).size:
        raise ValueError("train and test indexes must be disjoint")
    train = values[train_idx].astype(np.int64)
    test = values[test_idx].astype(np.int64)
    intersections = test @ train.T
    unions = test.sum(axis=1)[:, None] + train.sum(axis=1)[None, :] - intersections
    similarities = np.divide(
        intersections, unions, out=np.zeros_like(intersections, dtype=float), where=unions > 0
    )
    return similarities.max(axis=1)


def _wilson_interval(positive: int, total: int) -> list[float]:
    """95% Wilson interval for the descriptive observed active fraction."""
    z = 1.959963984540054
    rate = positive / total
    denominator = 1 + z * z / total
    center = (rate + z * z / (2 * total)) / denominator
    half = z * np.sqrt(rate * (1 - rate) / total + z * z / (4 * total * total)) / denominator
    return [float(max(0, center - half)), float(min(1, center + half))]


def neighborhood_calibration_diagnostics(
    y_true: np.ndarray, y_probability: np.ndarray, nearest_training_tanimoto: np.ndarray
) -> dict:
    """Report binwise probability/observed-rate gaps by training-neighbor band.

    This summarizes supplied frozen-split predictions. It does not fit a
    calibrator, choose a decision threshold, or establish external validity.
    """
    y = np.asarray(y_true, dtype=float)
    probability = np.asarray(y_probability, dtype=float)
    similarity = np.asarray(nearest_training_tanimoto, dtype=float)
    if (
        y.ndim != 1
        or probability.ndim != 1
        or similarity.ndim != 1
        or not len(y)
        or len(y) != len(probability)
        or len(y) != len(similarity)
        or not np.isfinite(y).all()
        or not np.isfinite(probability).all()
        or not np.isfinite(similarity).all()
        or not np.isin(y, [0, 1]).all()
        or np.any(probability < 0)
        or np.any(probability > 1)
        or np.any(similarity < 0)
        or np.any(similarity > 1)
    ):
        raise ValueError(
            "labels, probabilities and similarities must be aligned finite vectors in range"
        )
    bands = []
    edges = (
        (0.0, 0.25, "[0.00, 0.25]"),
        (0.25, 0.50, "(0.25, 0.50]"),
        (0.50, 0.75, "(0.50, 0.75]"),
        (0.75, 1.0, "(0.75, 1.00]"),
    )
    for low, high, label in edges:
        mask = ((similarity >= low) if low == 0 else (similarity > low)) & (similarity <= high)
        indices = np.flatnonzero(mask)
        if not len(indices):
            bands.append(
                {
                    "similarity_band": label,
                    "n": 0,
                    "mean_probability": None,
                    "observed_active_rate": None,
                    "absolute_calibration_gap": None,
                    "n_active": 0,
                    "n_inactive": 0,
                    "observed_active_rate_ci95": None,
                    "brier_score": None,
                }
            )
            continue
        mean_probability = float(probability[indices].mean())
        observed_rate = float(y[indices].mean())
        bands.append(
            {
                "similarity_band": label,
                "n": len(indices),
                "mean_probability": mean_probability,
                "observed_active_rate": observed_rate,
                "n_active": int(y[indices].sum()),
                "n_inactive": int(len(indices) - y[indices].sum()),
                "observed_active_rate_ci95": _wilson_interval(int(y[indices].sum()), len(indices)),
                "absolute_calibration_gap": abs(mean_probability - observed_rate),
                "brier_score": float(np.mean((probability[indices] - y[indices]) ** 2)),
            }
        )
    return {
        "status": "internal_frozen_split_diagnostic_not_external_calibration",
        "similarity_metric": "maximum binary-fingerprint Tanimoto to split training partition",
        "prediction_calibration_fit": False,
        "interval_method": "95% Wilson binomial interval; descriptive only, does not account for within-scaffold dependence",
        "bands": bands,
    }
