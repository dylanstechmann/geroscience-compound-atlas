"""Tests for Bemis-Murcko scaffold splitting, disjointness, and reproducibility."""

import numpy as np
import pandas as pd
from rdkit.Chem.Scaffolds import MurckoScaffold

from bench.splits import bemis_murcko_split, random_split


def test_bemis_murcko_scaffold_disjointness():
    """Requirement: Train, Validation, and Test sets must have zero intersecting Murcko scaffolds."""
    # Create test dataframe with diverse molecules and repeating scaffolds
    test_smiles = [
        "c1ccccc1CC",  # benzene scaffold
        "c1ccccc1CCC",  # benzene scaffold
        "c1ccccc1CCCC",  # benzene scaffold
        "c1ccc2ccccc2c1C",  # naphthalene scaffold
        "c1ccc2ccccc2c1CC",  # naphthalene scaffold
        "c1ccncc1C",  # pyridine scaffold
        "c1ccncc1CC",  # pyridine scaffold
        "c1ccncc1CCC",  # pyridine scaffold
        "C1CCCCC1O",  # cyclohexane scaffold
        "C1CCCCC1CO",  # cyclohexane scaffold
        "CCCC",  # acyclic (empty scaffold)
        "CCCCC",  # acyclic (empty scaffold)
    ]

    scaffolds = [MurckoScaffold.MurckoScaffoldSmiles(s) for s in test_smiles]
    df = pd.DataFrame(
        {
            "canonical_smiles": test_smiles,
            "murcko_scaffold": scaffolds,
            "active": [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
        }
    )

    train_idx, val_idx, test_idx = bemis_murcko_split(
        df, train_ratio=0.6, val_ratio=0.2, test_ratio=0.2, seed=42
    )

    # 1. Zero index overlap
    all_indices = set(train_idx) | set(val_idx) | set(test_idx)
    assert len(all_indices) == len(train_idx) + len(val_idx) + len(test_idx)
    assert len(all_indices) == len(df)

    # 2. Strict scaffold disjointness
    train_scaffolds = set(df.iloc[train_idx]["murcko_scaffold"])
    test_scaffolds = set(df.iloc[test_idx]["murcko_scaffold"])
    val_scaffolds = set(df.iloc[val_idx]["murcko_scaffold"])

    assert not (train_scaffolds & test_scaffolds), "Scaffold leakage between train and test!"
    if val_scaffolds:
        assert not (train_scaffolds & val_scaffolds), "Scaffold leakage between train and val!"
        assert not (val_scaffolds & test_scaffolds), "Scaffold leakage between val and test!"


def test_split_reproducibility():
    """Verify deterministic split generation given the same random seed."""
    df = pd.DataFrame(
        {
            "murcko_scaffold": ["scaff_A"] * 10 + ["scaff_B"] * 10 + ["scaff_C"] * 10,
            "active": [1] * 15 + [0] * 15,
        }
    )

    train_1, _, test_1 = bemis_murcko_split(df, seed=123)
    train_2, _, test_2 = bemis_murcko_split(df, seed=123)

    assert np.array_equal(train_1, train_2)
    assert np.array_equal(test_1, test_2)


def test_random_split_stratification():
    """Verify that random split preserves label distribution."""
    df = pd.DataFrame(
        {
            "murcko_scaffold": ["s1"] * 50 + ["s2"] * 50,
            "active": [1] * 40 + [0] * 60,
        }
    )

    train_idx, val_idx, test_idx = random_split(
        df, train_ratio=0.8, val_ratio=0.1, test_ratio=0.1, seed=42
    )
    assert len(train_idx) == 80
    assert len(val_idx) == 10
    assert len(test_idx) == 10

    train_prevalence = df.iloc[train_idx]["active"].mean()
    test_prevalence = df.iloc[test_idx]["active"].mean()
    assert abs(train_prevalence - test_prevalence) <= 0.05
