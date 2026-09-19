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

        with open(splits_json, "r", encoding="utf-8") as f:
            splits_data = json.load(f)

        # Retrieve scaffold train indices
        seed_key = str(seed)
        train_indices = splits_data["scaffold"][seed_key]["train"]
        train_idx = np.array(train_indices, dtype=int)

        # Assemble training matrices
        mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        fps_list = []
        for smi in self.benchmark_df["canonical_smiles"]:
            m = Chem.MolFromSmiles(str(smi)) if pd.notna(smi) else None
            if m:
                fps_list.append(np.array(mfpgen.GetFingerprint(m), dtype=np.float32))
            else:
                fps_list.append(np.zeros(2048, dtype=np.float32))
        fps = np.array(fps_list, dtype=np.float32)
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
            }

        try:
            # 1. Compute 2048-bit Morgan Fingerprint
            fp = compute_morgan_fingerprint(mol, n_bits=2048, radius=2)
            fp_arr = np.array([float(bit) for bit in fp], dtype=np.float32).reshape(1, -1)

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
            }
        except (ValueError, RuntimeError, KeyError):
            return {
                "reward": 0.0,
                "mtor_prob": 0.0,
                "qed": 0.0,
                "has_pains": False,
                "valid": False,
            }
