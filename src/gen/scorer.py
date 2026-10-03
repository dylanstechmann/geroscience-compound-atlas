"""Composite surrogate scorer with QED bias control and PAINS penalty."""

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import QED, FilterCatalog, rdFingerprintGenerator
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from atlas.features import compute_morgan_fingerprint, compute_rdkit_descriptors
from bench.config import load_bench_config
from bench.data import activity_labels
from bench.models import DESCRIPTOR_COLS


class CompositeSurrogateScorer:
    """Surrogate scoring engine evaluating predicted mTOR activity, QED bias, and PAINS."""

    def __init__(
        self,
        benchmark_parquet: str | Path = "artifacts/benchmark_dataset.parquet",
        splits_json: str | Path = "artifacts/splits.json",
        seed: int = 42,
    ):
        self.benchmark_df = pd.read_parquet(benchmark_parquet)
        config = load_bench_config()
        self.pchembl_threshold = float(config.thresholds.pchembl_active)
        self.label_definition = (
            f"active = 1 if pchembl_value >= {self.pchembl_threshold:g} else 0"
        )
        if "pchembl_value" not in self.benchmark_df or "active" not in self.benchmark_df:
            raise ValueError("Benchmark requires pchembl_value and active columns")
        expected_labels = activity_labels(
            self.benchmark_df["pchembl_value"], self.pchembl_threshold
        ).to_numpy(dtype=int)
        actual_labels = pd.to_numeric(self.benchmark_df["active"], errors="raise").to_numpy(dtype=int)
        if not np.array_equal(actual_labels, expected_labels):
            raise ValueError(
                "Benchmark active labels do not match configured threshold: "
                f"{self.label_definition}"
            )

        with open(splits_json, "r", encoding="utf-8") as f:
            splits_data = json.load(f)

        # Retrieve scaffold train indices
        seed_key = str(seed)
        train_indices = splits_data["scaffold"][seed_key]["train"]
        train_idx = np.array(train_indices, dtype=int)
        self.split_seed = seed
        # Generation must use the same positional partition as model fitting.
        self.training_df = self.benchmark_df.iloc[train_idx].copy()

        # Assemble training matrices
        fp_sidecar = Path(benchmark_parquet).with_suffix(".fingerprints.npy")
        if fp_sidecar.exists():
            fps = np.load(fp_sidecar)
        else:
            mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
            fps_list = []
            for smi in self.benchmark_df["canonical_smiles"]:
                m = Chem.MolFromSmiles(str(smi)) if pd.notna(smi) else None
                if m:
                    fps_list.append(np.array(mfpgen.GetFingerprint(m), dtype=np.float32))
                else:
                    fps_list.append(np.zeros(2048, dtype=np.float32))
            fps = np.array(fps_list, dtype=np.float32)
        if fps.ndim != 2 or fps.shape[0] != len(self.benchmark_df) or fps.shape[1] != 2048:
            raise ValueError("Fingerprint matrix must have one 2048-bit row per benchmark compound")
        self.training_fingerprints = fps[train_idx].astype(bool)
        desc_raw = self.benchmark_df[DESCRIPTOR_COLS].to_numpy(dtype=np.float32)

        # Fit scaler on training descriptors only
        self.scaler = StandardScaler()
        desc_train_scaled = self.scaler.fit_transform(desc_raw[train_idx])
        X_train = np.hstack([fps[train_idx], desc_train_scaled])
        y_train = self.benchmark_df["active"].to_numpy()[train_idx]

        # Train frozen baseline model
        self.model = LogisticRegression(C=1.0, max_iter=1000, random_state=seed, solver="lbfgs")
        self.model.fit(X_train, y_train)

        # Load RDKit PAINS filter catalog
        pains_params = FilterCatalog.FilterCatalogParams()
        pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        self.pains_catalog = FilterCatalog.FilterCatalog(pains_params)

    def has_pains(self, mol: Chem.Mol) -> bool:
        """Check if a molecule matches PAINS sub-structures."""
        try:
            return bool(self.pains_catalog.HasMatch(mol))
        except (ValueError, RuntimeError):
            return False

    def score_molecule(
        self,
        mol: Chem.Mol | None,
        qed_weight: float = 0.2,
        mtor_weight: float = 1.0,
        pains_penalty: float = 0.5,
    ) -> dict[str, Any]:
        """Compute multi-objective reward for a single candidate molecule."""
        if mol is None:
            return {
                "reward": 0.0,
                "mtor_prob": 0.0,
                "qed": 0.0,
                "has_pains": False,
                "valid": False,
                "nearest_training_tanimoto": None,
                "applicability_domain_status": "not_validated",
            }

        try:
            # 1. Compute 2048-bit Morgan Fingerprint
            fp = compute_morgan_fingerprint(mol, n_bits=2048, radius=2)
            fp_arr = np.array([float(bit) for bit in fp], dtype=np.float32).reshape(1, -1)
            query_bits = fp_arr.astype(bool)
            intersections = np.logical_and(self.training_fingerprints, query_bits).sum(axis=1)
            unions = np.logical_or(self.training_fingerprints, query_bits).sum(axis=1)
            similarities = np.divide(
                intersections,
                unions,
                out=np.zeros_like(intersections, dtype=np.float64),
                where=unions > 0,
            )
            nearest_training_tanimoto = float(similarities.max())

            # 2. Compute 9 RDKit 2D Descriptors
            desc_dict = compute_rdkit_descriptors(mol)
            desc_vals = np.array(
                [desc_dict[col] for col in DESCRIPTOR_COLS], dtype=np.float32
            ).reshape(1, -1)

            # 3. Standardize descriptors
            desc_scaled = self.scaler.transform(desc_vals)
            X = np.hstack([fp_arr, desc_scaled])

            # 4. Predict probability of mTOR activity
            mtor_prob = float(self.model.predict_proba(X)[0, 1])

            # 5. Compute QED score (bias control)
            qed_score = float(QED.qed(mol))

            # 6. Check PAINS filter
            is_pains = self.has_pains(mol)
            penalty = pains_penalty if is_pains else 0.0

            # 7. Composite Reward
            reward = (mtor_weight * mtor_prob) + (qed_weight * qed_score) - penalty
            reward = max(0.0, reward)

            return {
                "reward": reward,
                "mtor_prob": mtor_prob,
                "qed": qed_score,
                "has_pains": is_pains,
                "valid": True,
                "nearest_training_tanimoto": nearest_training_tanimoto,
                "applicability_domain_status": "not_validated",
            }
        except (ValueError, RuntimeError, KeyError):
            return {
                "reward": 0.0,
                "mtor_prob": 0.0,
                "qed": 0.0,
                "has_pains": False,
                "valid": False,
                "nearest_training_tanimoto": None,
                "applicability_domain_status": "not_validated",
            }
