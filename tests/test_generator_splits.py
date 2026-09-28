"""Offline regressions for held-out data leaking into molecular generation."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem

from atlas.features import compute_rdkit_descriptors
from bench.models import DESCRIPTOR_COLS
from gen import hypothesis, pipeline
from gen.scorer import CompositeSurrogateScorer

HYPOTHESIS_CONFIG = Path(__file__).resolve().parents[1] / "configs/hypothesis.yaml"


@pytest.fixture
def frozen_data(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Held-out actives deliberately outrank the training active by pChEMBL.
    smiles = [
        "Cc1nc(Nc2ncc(s2)C(=O)Nc2c(C)cccc2Cl)cc(n1)N1CCN(CCO)CC1",
        "CCCCN",
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "CCCO",
    ]
    frame = pd.DataFrame(
        [
            {
                "canonical_smiles": smi,
                "inchikey": Chem.MolToInchiKey(Chem.MolFromSmiles(smi)),
                **compute_rdkit_descriptors(smi),
            }
            for smi in smiles
        ],
        index=[101, 205, 309, 413],
    )
    frame["active"] = [1, 0, 1, 1]
    frame["pchembl_value"] = [12.0, 4.0, 11.0, 7.0]
    benchmark = tmp_path / "benchmark.parquet"
    frame.to_parquet(benchmark)
    splits = tmp_path / "splits.json"
    splits.write_text(
        json.dumps(
            {
                "scaffold": {
                    "42": {"train": [3, 1], "val": [0], "test": [2]},
                    "123": {"train": [1, 2], "val": [3], "test": [0]},
                }
            }
        )
    )
    return frame, benchmark, splits


@pytest.mark.parametrize("seed,positions", [(42, [3, 1]), (123, [1, 2])])
def test_scorer_exposes_the_same_training_partition_it_fits(frozen_data, seed, positions):
    frame, benchmark, splits = frozen_data
    scorer = CompositeSurrogateScorer(benchmark, splits, seed=seed)
    expected = frame.iloc[positions]
    pd.testing.assert_frame_equal(scorer.training_df, expected)
    assert scorer.split_seed == seed
    np.testing.assert_allclose(
        scorer.scaler.mean_,
        expected[DESCRIPTOR_COLS].to_numpy(dtype=np.float32).mean(axis=0),
        rtol=1e-6,
    )


def stub_archive(frame, monkeypatch, module):
    """Return held-out structures independently of seeds to audit reference scope."""
    calls = []
    held_out = frame.iloc[[0, 2]]
    archive = pd.DataFrame(
        {
            "smiles": held_out["canonical_smiles"].tolist(),
            "inchikey": held_out["inchikey"].tolist(),
            "generation": [0, 0],
            "mtor_prob": [0.9, 0.8],
            "qed": held_out["qed"].tolist(),
            "has_pains": [False, False],
            "reward": [1.0, 0.9],
        }
    )

    class StubGA:
        def __init__(self, **kwargs):
            self.qed_weight = kwargs["qed_weight"]

        def run(self, seed_smiles, max_archive_size):
            calls.append((self.qed_weight, seed_smiles))
            return archive.copy(), {"validity_rate": 1.0, "archive_size": len(archive)}

    monkeypatch.setattr(module, "MolecularGA", StubGA)
    return calls


def test_hypothesis_uses_only_training_seeds_and_active_references(frozen_data, monkeypatch):
    frame, benchmark, splits = frozen_data
    calls = stub_archive(frame, monkeypatch, hypothesis)
    accepted, metrics = hypothesis.run_hypothesis_pipeline(
        config_path=HYPOTHESIS_CONFIG,
        benchmark_parquet=benchmark,
        splits_json=splits,
    )
    assert calls == [(0.25, ["CCCO"])]
    # A held-out exact match must not be mistaken for a training clone.
    assert set(accepted["inchikey"]) == set(frame.iloc[[0, 2]]["inchikey"])
    refs = [hypothesis.fingerprint(Chem.MolFromSmiles("CCCO"))]
    for row in accepted.itertuples():
        assert row.nearest_tanimoto == pytest.approx(
            hypothesis.nearest_tanimoto(Chem.MolFromSmiles(row.smiles), refs)
        )
    assert metrics["reject_counts"] == {}
    assert metrics["training_reference"] == {
        "split": "scaffold",
        "seed": 42,
        "rows": 2,
        "active_rows": 1,
    }


def test_exact_training_check_includes_inactive_rows(frozen_data):
    frame, _, _ = frozen_data
    fps, keys = hypothesis.training_active_fingerprints(frame.iloc[[3, 1]])
    assert len(fps) == 1
    assert keys == set(frame.iloc[[3, 1]]["inchikey"])


def test_play_mode_uses_training_seeds_and_novelty_for_every_qed_weight(
    frozen_data, monkeypatch, tmp_path
):
    frame, benchmark, splits = frozen_data
    calls = stub_archive(frame, monkeypatch, pipeline)
    atlas = tmp_path / "compounds.parquet"
    pd.DataFrame(
        {
            "resolved": [True],
            "canonical_smiles": ["CCNCC"],
            "resolved_name": ["fixture"],
            "raw_name": ["fixture"],
        }
    ).to_parquet(atlas)
    _, metrics = pipeline.run_generator_pipeline(
        benchmark_parquet=benchmark,
        splits_json=splits,
        compounds_parquet=atlas,
    )
    assert calls == [(0.0, ["CCCO"]), (0.2, ["CCCO"]), (0.5, ["CCCO"])]
    assert metrics["primary_run"]["novelty_vs_training"] == 1.0
    for result in metrics["qed_bias_sensitivity"].values():
        assert result["novelty_rate"] == 1.0


def test_held_out_changes_do_not_change_hypothesis_outputs(frozen_data, monkeypatch, tmp_path):
    frame, benchmark, splits = frozen_data
    calls = stub_archive(frame, monkeypatch, hypothesis)
    kwargs = {
        "config_path": HYPOTHESIS_CONFIG,
        "benchmark_parquet": benchmark,
        "splits_json": splits,
    }
    original, original_metrics = hypothesis.run_hypothesis_pipeline(
        **kwargs, output_dir=tmp_path / "original"
    )
    # Alter every held-out field used for seed ranking, novelty, similarity or scoring.
    changed = frame.copy()
    for label in changed.index[[0, 2]]:
        changed.loc[label, "canonical_smiles"] = "CCCCCCCC"
        changed.loc[label, "inchikey"] = Chem.MolToInchiKey(Chem.MolFromSmiles("CCCCCCCC"))
        changed.loc[label, "active"] = 0
        changed.loc[label, "pchembl_value"] = 2.0
        for key, value in compute_rdkit_descriptors("CCCCCCCC").items():
            changed.loc[label, key] = value
    changed.to_parquet(benchmark)
    repeated, repeated_metrics = hypothesis.run_hypothesis_pipeline(
        **kwargs, output_dir=tmp_path / "changed"
    )
    pd.testing.assert_frame_equal(original, repeated)
    assert original_metrics == repeated_metrics
    assert calls == [(0.25, ["CCCO"]), (0.25, ["CCCO"])]
