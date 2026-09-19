"""Splitting methodologies: Bemis-Murcko scaffold disjoint splits and random leakage controls."""

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedShuffleSplit


def bemis_murcko_split(
    df: pd.DataFrame,
    train_ratio: float = 0.80,
    val_ratio: float = 0.10,
    test_ratio: float = 0.10,
    seed: int = 42,
    scaffold_col: str = "murcko_scaffold",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Partition dataset into train, validation, and test sets strictly grouped by Murcko scaffold.

    Guarantees:
        scaffolds(train) ∩ scaffolds(val) = ∅
        scaffolds(train) ∩ scaffolds(test) = ∅
        scaffolds(val) ∩ scaffolds(test) = ∅

    Args:
        df: DataFrame containing chemical records and a scaffold column.
        train_ratio: Fraction allocated to train (default 0.80).
        val_ratio: Fraction allocated to val (default 0.10).
        test_ratio: Fraction allocated to test (default 0.10).
        seed: Random seed for deterministic shuffling of scaffold groups.
        scaffold_col: Name of column containing Murcko scaffold SMILES.

    Returns:
        Tuple of (train_indices, val_indices, test_indices) as numpy arrays.
    """
    total_samples = len(df)
    train_target = round(train_ratio * total_samples)
    val_target = round(val_ratio * total_samples)

    # Group row indices by scaffold
    scaffold_to_indices: dict[str, list[int]] = {}
    for idx, scaffold in enumerate(df[scaffold_col]):
        scaffold_to_indices.setdefault(scaffold, []).append(idx)

    # Sort scaffolds deterministically then shuffle with seed
    rng = np.random.RandomState(seed)
    scaffolds = sorted(scaffold_to_indices.keys())
    rng.shuffle(scaffolds)

    train_indices: list[int] = []
    val_indices: list[int] = []
    test_indices: list[int] = []

    for scaffold in scaffolds:
        group = scaffold_to_indices[scaffold]
        if len(train_indices) + len(group) <= train_target:
            train_indices.extend(group)
        elif len(val_indices) + len(group) <= val_target:
            val_indices.extend(group)
        else:
            test_indices.extend(group)

    # In rare cases where test set is empty due to group sizes, balance from largest set
    if not test_indices and val_indices:
        test_indices = val_indices
        val_indices = []

    train_idx = np.array(train_indices, dtype=int)
    val_idx = np.array(val_indices, dtype=int)
    test_idx = np.array(test_indices, dtype=int)

    # Verify disjointness
    train_scaffolds = set(df.iloc[train_idx][scaffold_col])
    val_scaffolds = set(df.iloc[val_idx][scaffold_col])
    test_scaffolds = set(df.iloc[test_idx][scaffold_col])

    assert not (train_scaffolds & test_scaffolds), (
        "Scaffold leakage detected between Train and Test!"
    )
    if val_scaffolds:
        assert not (train_scaffolds & val_scaffolds), (
            "Scaffold leakage detected between Train and Val!"
        )
        assert not (val_scaffolds & test_scaffolds), (
            "Scaffold leakage detected between Val and Test!"
        )

    return train_idx, val_idx, test_idx


def random_split(
    df: pd.DataFrame,
    train_ratio: float = 0.80,
    val_ratio: float = 0.10,
    test_ratio: float = 0.10,
    seed: int = 42,
    target_col: str = "active",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Standard stratified random train/val/test split used strictly as a leakage diagnostic.

    Args:
        df: Input DataFrame.
        train_ratio: Fraction for training.
        val_ratio: Fraction for validation.
        test_ratio: Fraction for testing.
        seed: Random seed.
        target_col: Stratification label column.

    Returns:
        Tuple of (train_indices, val_indices, test_indices).
    """
    y = df[target_col].to_numpy()
    test_val_ratio = val_ratio + test_ratio

    sss1 = StratifiedShuffleSplit(n_splits=1, test_size=test_val_ratio, random_state=seed)
    train_idx, temp_idx = next(sss1.split(df, y))

    relative_test_ratio = test_ratio / test_val_ratio
    sss2 = StratifiedShuffleSplit(n_splits=1, test_size=relative_test_ratio, random_state=seed)
    val_sub_idx, test_sub_idx = next(sss2.split(temp_idx, y[temp_idx]))

    val_idx = temp_idx[val_sub_idx]
    test_idx = temp_idx[test_sub_idx]

    return train_idx, val_idx, test_idx


def generate_and_save_splits(
    df: pd.DataFrame,
    seeds: list[int] = (42, 123, 456),
    output_path: str | Path = "artifacts/splits.json",
) -> dict[str, dict[str, dict[str, list[int]]]]:
    """Generate and save both scaffold and random splits across multiple seeds."""
    all_splits: dict[str, Any] = {"scaffold": {}, "random": {}}

    for seed in seeds:
        sc_train, sc_val, sc_test = bemis_murcko_split(df, seed=seed)
        all_splits["scaffold"][str(seed)] = {
            "train": sc_train.tolist(),
            "val": sc_val.tolist(),
            "test": sc_test.tolist(),
        }

        rd_train, rd_val, rd_test = random_split(df, seed=seed)
        all_splits["random"][str(seed)] = {
            "train": rd_train.tolist(),
            "val": rd_val.tolist(),
            "test": rd_test.tolist(),
        }

    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(all_splits, f, indent=2)

    return all_splits
