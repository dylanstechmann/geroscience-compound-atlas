#!/usr/bin/env python3
"""Guards for evaluating continuous activity on the frozen scaffold partition."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from bench.regression import checked_scaffold_indices, run_regression


def test_split_guard_rejects_duplicate_rows_and_shared_scaffolds():
    df = pd.DataFrame({"murcko_scaffold": ["A", "A", "B", "C"]})
    valid = {"train": [0, 1], "val": [2], "test": [3]}
    checked_scaffold_indices(df, valid)

    with pytest.raises(ValueError, match="partition"):
        checked_scaffold_indices(df, {"train": [0, 0], "val": [2], "test": [3]})
    with pytest.raises(ValueError, match="overlapping"):
        checked_scaffold_indices(df, {"train": [0], "val": [1, 2], "test": [3]})


def test_regression_uses_frozen_training_mean(tmp_path):
    dataset = "artifacts/benchmark_dataset.parquet"
    splits = "artifacts/splits.json"
    output = tmp_path / "regression.json"
    result = run_regression(dataset, splits, output)
    df = pd.read_parquet(dataset)
    frozen = json.loads(Path(splits).read_text(encoding="utf-8"))["scaffold"]

    assert result["n_molecules"] == len(df)
    assert set(result["per_seed"]) == set(frozen)
    assert output.exists()
    for seed, partition in frozen.items():
        row = result["per_seed"][seed]
        train_mean = float(np.mean(df.iloc[partition["train"]]["pchembl_value"]))
        test_values = df.iloc[partition["test"]]["pchembl_value"].to_numpy()
        assert row["train_pchembl_mean"] == round(train_mean, 4)
        assert row["train_mean"]["mae"] == round(float(np.mean(abs(test_values - train_mean))), 4)
        assert row["n_test"] == len(partition["test"])
