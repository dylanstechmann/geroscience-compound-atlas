"""Tests for benchmark feature preparation, model evaluation, and false positive diagnostics."""

import numpy as np
import pandas as pd

from bench.models import (
    compute_metrics,
    extract_top_false_positives,
    prepare_feature_matrices,
    train_and_eval_baseline,
    train_and_eval_contender,
)


def test_compute_metrics():
    """Verify classification metric calculation."""
    y_true = np.array([1, 1, 0, 0, 1, 0])
    y_prob = np.array([0.9, 0.8, 0.2, 0.1, 0.7, 0.3])

    metrics = compute_metrics(y_true, y_prob)

    assert "auroc" in metrics
    assert "auprc" in metrics
    assert "recall_at_fpr_0.05" in metrics
    assert "brier_score" in metrics

    assert metrics["auroc"] == 1.0  # Perfectly separated
    assert 0.0 <= metrics["brier_score"] <= 1.0


def test_prepare_feature_matrices_and_models():
    """Verify feature matrix building and model training on synthetic data."""
    n_samples = 40
    n_bits = 2048

    rng = np.random.RandomState(42)
    synthetic_fps = rng.randint(0, 2, size=(n_samples, n_bits)).tolist()

    df = pd.DataFrame(
        {
            "fingerprint": synthetic_fps,
            "mol_wt": rng.uniform(200, 600, size=n_samples),
            "log_p": rng.uniform(-1, 5, size=n_samples),
            "tpsa": rng.uniform(20, 140, size=n_samples),
            "num_h_donors": rng.randint(0, 5, size=n_samples),
            "num_h_acceptors": rng.randint(1, 10, size=n_samples),
            "num_rotatable_bonds": rng.randint(0, 10, size=n_samples),
            "ring_count": rng.randint(1, 5, size=n_samples),
            "fraction_csp3": rng.uniform(0.1, 0.9, size=n_samples),
            "qed": rng.uniform(0.2, 0.9, size=n_samples),
            "active": rng.randint(0, 2, size=n_samples),
            "canonical_smiles": ["CC"] * n_samples,
            "murcko_scaffold": ["c1ccccc1"] * n_samples,
            "inchikey": [f"KEY_{i}" for i in range(n_samples)],
        }
    )

    train_idx = np.arange(0, 28)
    val_idx = np.arange(28, 34)
    test_idx = np.arange(34, 40)

    X_train, y_train, _X_val, _y_val, X_test, y_test = prepare_feature_matrices(
        df, train_idx, val_idx, test_idx
    )

    # Feature dimensionality: 2048 Morgan bits + 9 RDKit descriptors = 2057
    assert X_train.shape == (28, 2057)
    assert X_test.shape == (6, 2057)

    # Train baseline
    base_metrics, base_probs = train_and_eval_baseline(X_train, y_train, X_test, y_test, seed=42)
    assert "auroc" in base_metrics
    assert len(base_probs) == 6

    # Train contender
    cont_metrics, cont_probs = train_and_eval_contender(X_train, y_train, X_test, y_test, seed=42)
    assert "auroc" in cont_metrics
    assert len(cont_probs) == 6


def test_extract_top_false_positives():
    """Verify false positive ranking utility."""
    df = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
            "inchikey": ["K1", "K2", "K3"],
            "canonical_smiles": ["C1", "C2", "C3"],
            "murcko_scaffold": ["S1", "S2", "S3"],
            "pchembl_value": [5.2, 4.8, 5.5],
        }
    )

    test_idx = np.array([0, 1, 2])
    y_test = np.array([0, 0, 1])  # First two are inactive
    y_prob = np.array([0.95, 0.75, 0.85])

    fps = extract_top_false_positives(df, test_idx, y_test, y_prob, top_n=2)

    assert len(fps) == 2
    # Highest false positive confidence first
    assert fps[0]["molecule_chembl_id"] == "CHEMBL1"
    assert fps[0]["predicted_prob"] == 0.95
    assert fps[1]["molecule_chembl_id"] == "CHEMBL2"
    assert fps[1]["predicted_prob"] == 0.75
