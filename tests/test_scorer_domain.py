"""Checks that the surrogate label definition and novelty diagnostic stay honest."""

import pandas as pd
import pytest
from rdkit import Chem

from gen.scorer import CompositeSurrogateScorer


def test_scorer_uses_configured_threshold_and_nearest_training_similarity():
    scorer = CompositeSurrogateScorer()
    assert scorer.label_definition == "active = 1 if pchembl_value >= 6 else 0"
    smiles = scorer.training_df.iloc[0]["canonical_smiles"]
    result = scorer.score_molecule(Chem.MolFromSmiles(smiles))
    assert result["valid"] is True
    assert result["nearest_training_tanimoto"] == pytest.approx(1.0)
    assert result["applicability_domain_status"] == "not_validated"


def test_scorer_rejects_labels_that_disagree_with_configured_threshold(monkeypatch):
    read_parquet = pd.read_parquet

    def mismatched_labels(path, *args, **kwargs):
        frame = read_parquet(path, *args, **kwargs)
        if str(path).endswith("benchmark_dataset.parquet"):
            frame = frame.copy()
            frame.loc[frame.index[0], "active"] = 1 - int(frame.loc[frame.index[0], "active"])
        return frame

    monkeypatch.setattr(pd, "read_parquet", mismatched_labels)
    with pytest.raises(ValueError, match="do not match configured threshold"):
        CompositeSurrogateScorer()
