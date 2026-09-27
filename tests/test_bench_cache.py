"""Offline regressions for threshold changes and cached benchmark measurements."""

import json
from unittest.mock import Mock

import numpy as np
import pandas as pd
import pytest

from bench import data, train
from bench.config import load_bench_config
from bench.models import DESCRIPTOR_COLS


@pytest.fixture
def offline_pipeline(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fetch = Mock(side_effect=AssertionError("Tests must not fetch live ChEMBL data"))
    monkeypatch.setattr(data, "fetch_chembl_mtor_records", fetch)
    monkeypatch.setattr(train, "plot_scaffold_size_distribution", lambda *a, **kw: "unused.png")
    monkeypatch.setattr(train, "plot_benchmark_curves", lambda *a, **kw: "unused.png")

    fitted_labels = []

    def fit(_x_train, y_train, _x_test, y_test, **kwargs):
        fitted_labels.append((y_train.copy(), y_test.copy()))
        metrics = {key: 0.5 for key in ("auroc", "auprc", "recall_at_fpr_0.05", "brier_score")}
        return metrics, np.full(len(y_test), 0.5)

    monkeypatch.setattr(train, "train_and_eval_baseline", fit)
    monkeypatch.setattr(train, "train_and_eval_contender", fit)
    return fetch, fitted_labels


def write_cache(tmp_path):
    values = np.tile([5.5, 6.0, 6.5, 7.0, 7.5], 12)
    frame = pd.DataFrame({
        "pchembl_value": values,
        "active": (values >= 6.0).astype(int),
        "murcko_scaffold": [f"scaffold-{i}" for i in range(len(values))],
        **{name: np.arange(len(values), dtype=float) for name in DESCRIPTOR_COLS},
    })
    path = tmp_path / "benchmark.parquet"
    frame.to_parquet(path, index=False)
    fingerprint_path = path.with_suffix(".fingerprints.npy")
    np.save(fingerprint_path, np.tile([0, 1, 0, 1], (len(values), 1)))
    return path, frame, fingerprint_path


def test_cached_threshold_changes_relabel_before_splits_and_fit(tmp_path, offline_pipeline):
    fetch, fitted_labels = offline_pipeline
    path, original, fingerprint_path = write_cache(tmp_path)
    original_fingerprints = fingerprint_path.read_bytes()
    config_path = tmp_path / "bench.yaml"
    config_path.write_text("benchmark:\n  thresholds:\n    pchembl_active: 7.0\n")
    threshold = load_bench_config(config_path).thresholds.pchembl_active

    # Exercise raising the threshold, lowering it, and reusing an unchanged cache.
    for cutoff in (threshold, 6.0, 6.0):
        fitted_labels.clear()
        result = train.run_benchmark_pipeline(
            dataset_parquet=path, seeds=[42], pchembl_threshold=cutoff
        )
        expected = (original["pchembl_value"] >= cutoff).astype(int)
        saved = pd.read_parquet(path)
        np.testing.assert_array_equal(saved["active"], expected)
        np.testing.assert_array_equal(saved["pchembl_value"], original["pchembl_value"])
        exported = pd.read_csv("data/processed/benchmark_dataset.csv")
        np.testing.assert_array_equal(exported["active"], expected)
        assert result["num_active"] == int(expected.sum())
        assert result["pchembl_threshold"] == cutoff
        assert result["label_definition"] == f"active = 1 if pchembl_value >= {cutoff} else 0"
        assert json.loads((tmp_path / "artifacts/metrics.json").read_text()) == result

        splits = json.loads((tmp_path / "artifacts/splits.json").read_text())
        assert len(fitted_labels) == 4
        for split_type, fits in zip(("scaffold", "random"), (fitted_labels[:2], fitted_labels[2:])):
            indices = splits[split_type]["42"]
            for train_labels, test_labels in fits:
                np.testing.assert_array_equal(train_labels, expected.iloc[indices["train"]])
                np.testing.assert_array_equal(test_labels, expected.iloc[indices["test"]])
    assert fingerprint_path.read_bytes() == original_fingerprints
    fetch.assert_not_called()


def test_cache_hit_recreates_csv_export_even_when_labels_already_match(tmp_path, offline_pipeline):
    fetch, _ = offline_pipeline
    path, original, _ = write_cache(tmp_path)
    csv_path = tmp_path / "data/processed/benchmark_dataset.csv"
    assert not csv_path.exists()
    train.run_benchmark_pipeline(dataset_parquet=path, seeds=[42], pchembl_threshold=6.0)
    np.testing.assert_array_equal(pd.read_csv(csv_path)["active"], original["active"])
    fetch.assert_not_called()


def test_prepared_measurements_keep_precision_across_cache_roundtrip(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(data, "fetch_chembl_mtor_records", lambda **kwargs: [
        {"canonical_smiles": "CCO", "pchembl_value": "5.9996", "molecule_chembl_id": "A"},
        {"canonical_smiles": "CCO", "pchembl_value": "6.0", "molecule_chembl_id": "A"},
        {"canonical_smiles": "CCN", "pchembl_value": "6.0", "molecule_chembl_id": "B"},
    ])
    path = tmp_path / "benchmark.parquet"
    fresh = data.prepare_benchmark_dataset(output_parquet=path, pchembl_threshold=6.0)
    cached = pd.read_parquet(path)
    below_cutoff = cached.loc[cached["molecule_chembl_id"] == "A"].iloc[0]
    assert below_cutoff["pchembl_value"] == pytest.approx(5.9998)
    assert below_cutoff["active"] == 0
    assert cached.loc[cached["molecule_chembl_id"] == "B", "active"].item() == 1
    np.testing.assert_array_equal(data.activity_labels(cached["pchembl_value"], 6.0), fresh["active"])


@pytest.mark.parametrize("threshold", [float("nan"), float("inf"), -float("inf"), None, True])
def test_invalid_threshold_fails_before_fetch_or_cache_mutation(tmp_path, offline_pipeline, threshold):
    fetch, _ = offline_pipeline
    path, _, _ = write_cache(tmp_path)
    original = path.read_bytes()
    with pytest.raises(ValueError, match="pchembl_threshold must be a finite number"):
        train.run_benchmark_pipeline(dataset_parquet=path, pchembl_threshold=threshold)
    with pytest.raises(ValueError, match="pchembl_threshold must be a finite number"):
        data.prepare_benchmark_dataset(output_parquet=path, pchembl_threshold=threshold)
    assert path.read_bytes() == original
    fetch.assert_not_called()


@pytest.mark.parametrize("bad_value", [float("nan"), float("inf"), -float("inf")])
def test_invalid_cached_measurements_fail_without_overwriting_cache(
    tmp_path, offline_pipeline, bad_value
):
    fetch, fitted_labels = offline_pipeline
    path, frame, _ = write_cache(tmp_path)
    frame.loc[0, "pchembl_value"] = bad_value
    frame.to_parquet(path, index=False)
    original = path.read_bytes()
    with pytest.raises(ValueError, match="measurements must be finite numbers"):
        train.run_benchmark_pipeline(dataset_parquet=path)
    assert path.read_bytes() == original
    assert fitted_labels == []
    fetch.assert_not_called()


def test_cached_measurements_are_required_for_relabeling(tmp_path, offline_pipeline):
    fetch, _ = offline_pipeline
    path, frame, _ = write_cache(tmp_path)
    frame.drop(columns=["pchembl_value"]).to_parquet(path, index=False)
    original = path.read_bytes()
    with pytest.raises(ValueError, match="missing pchembl_value"):
        train.run_benchmark_pipeline(dataset_parquet=path)
    assert path.read_bytes() == original
    fetch.assert_not_called()
