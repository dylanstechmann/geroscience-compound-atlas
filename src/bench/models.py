"""Model definitions, feature pipeline, evaluation metrics, and error diagnostics."""

from typing import Any

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler

DESCRIPTOR_COLS = [
    "mol_wt",
    "log_p",
    "tpsa",
    "num_h_donors",
    "num_h_acceptors",
    "num_rotatable_bonds",
    "ring_count",
    "fraction_csp3",
    "qed",
]


def prepare_feature_matrices(
    df: pd.DataFrame,
    train_idx: np.ndarray,
    val_idx: np.ndarray,
    test_idx: np.ndarray,
    include_descriptors: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Assemble standardized feature matrices combining Morgan fingerprints and scaled RDKit descriptors.

    Scales descriptors using training set parameters only to strictly prevent feature leakage.
    """
    fps = np.array(df["fingerprint"].tolist(), dtype=np.float32)

    if include_descriptors:
        desc_raw = df[DESCRIPTOR_COLS].to_numpy(dtype=np.float32)
        # Impute any missing values with median of training set
        train_desc = desc_raw[train_idx]
        col_medians = np.nanmedian(train_desc, axis=0)

        inds = np.where(np.isnan(desc_raw))
        desc_raw[inds] = np.take(col_medians, inds[1])

        scaler = StandardScaler()
        desc_train = scaler.fit_transform(desc_raw[train_idx])
        desc_val = scaler.transform(desc_raw[val_idx])
        desc_test = scaler.transform(desc_raw[test_idx])

        X_train = np.hstack([fps[train_idx], desc_train])
        X_val = np.hstack([fps[val_idx], desc_val])
        X_test = np.hstack([fps[test_idx], desc_test])
    else:
        X_train = fps[train_idx]
        X_val = fps[val_idx]
        X_test = fps[test_idx]

    y_train = df["active"].to_numpy()[train_idx]
    y_val = df["active"].to_numpy()[val_idx]
    y_test = df["active"].to_numpy()[test_idx]

    return X_train, y_train, X_val, y_val, X_test, y_test


def compute_metrics(y_true: np.ndarray, y_prob: np.ndarray) -> dict[str, float]:
    """Compute comprehensive classification and calibration metrics."""
    auroc = float(roc_auc_score(y_true, y_prob))
    auprc = float(average_precision_score(y_true, y_prob))
    brier = float(brier_score_loss(y_true, y_prob))

    # Recall at <= 5% FPR
    fpr, tpr, _ = roc_curve(y_true, y_prob)
    valid_tpr = tpr[fpr <= 0.05]
    recall_at_fpr05 = float(valid_tpr[-1]) if len(valid_tpr) > 0 else 0.0

    return {
        "auroc": round(auroc, 4),
        "auprc": round(auprc, 4),
        "recall_at_fpr_0.05": round(recall_at_fpr05, 4),
        "brier_score": round(brier, 4),
    }


def train_and_eval_baseline(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    y_test: np.ndarray,
    seed: int = 42,
) -> tuple[dict[str, float], np.ndarray]:
    """Train L2 Logistic Regression baseline and return metrics + test predictions."""
    model = LogisticRegression(C=1.0, max_iter=1000, random_state=seed, solver="lbfgs")
    model.fit(X_train, y_train)
    y_prob = model.predict_proba(X_test)[:, 1]
    metrics = compute_metrics(y_test, y_prob)
    return metrics, y_prob


def train_and_eval_contender(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    y_test: np.ndarray,
    seed: int = 42,
) -> tuple[dict[str, float], np.ndarray]:
    """Train HistGradientBoostingClassifier contender and return metrics + test predictions."""
    model = HistGradientBoostingClassifier(
        max_iter=200, learning_rate=0.05, min_samples_leaf=10, random_state=seed
    )
    model.fit(X_train, y_train)
    y_prob = model.predict_proba(X_test)[:, 1]
    metrics = compute_metrics(y_test, y_prob)
    return metrics, y_prob


def extract_top_false_positives(
    df: pd.DataFrame,
    test_idx: np.ndarray,
    y_test: np.ndarray,
    y_prob: np.ndarray,
    top_n: int = 10,
) -> list[dict[str, Any]]:
    """Identify top N false positives with highest confidence."""
    test_df = df.iloc[test_idx].copy()
    test_df["predicted_prob"] = y_prob
    test_df["true_label"] = y_test

    # False positives: true_label == 0, sorted descending by predicted_prob
    fps = test_df[test_df["true_label"] == 0].sort_values(by="predicted_prob", ascending=False)
    selected = fps.head(top_n)

    mistakes = []
    for _, row in selected.iterrows():
        mistakes.append(
            {
                "molecule_chembl_id": row.get("molecule_chembl_id"),
                "inchikey": row.get("inchikey"),
                "canonical_smiles": row.get("canonical_smiles"),
                "murcko_scaffold": row.get("murcko_scaffold"),
                "predicted_prob": round(float(row["predicted_prob"]), 4),
                "true_pchembl": row.get("pchembl_value"),
                "true_label": int(row["true_label"]),
            }
        )
    return mistakes
