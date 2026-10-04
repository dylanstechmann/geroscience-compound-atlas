import numpy as np
import pandas as pd
import pytest

from bench.neighborhood import maximum_training_tanimoto, neighborhood_calibration_diagnostics
from bench.neighborhood_report import evaluate_frozen_scaffold_neighborhoods


def test_maximum_training_tanimoto_uses_training_rows_only():
    fingerprints = np.array(
        [
            [1, 1, 0, 0],  # train A
            [0, 0, 1, 1],  # train B
            [1, 1, 1, 0],  # test: 2/3 with A, 1/4 with B
        ]
    )
    np.testing.assert_allclose(maximum_training_tanimoto(fingerprints, [0, 1], [2]), [2 / 3])
    with pytest.raises(ValueError, match="nonempty train and test"):
        maximum_training_tanimoto(fingerprints, [], [2])


def test_neighborhood_calibration_reports_empty_bands_and_does_not_fit():
    result = neighborhood_calibration_diagnostics(
        np.array([0, 1, 1]), np.array([0.1, 0.7, 0.8]), np.array([0.1, 0.6, 0.9])
    )
    assert result["prediction_calibration_fit"] is False
    assert result["status"] == "internal_frozen_split_diagnostic_not_external_calibration"
    assert result["bands"][0]["absolute_calibration_gap"] == pytest.approx(0.1)
    assert result["bands"][1]["n"] == 0
    assert result["bands"][2]["n"] == 1
    assert result["bands"][3]["brier_score"] == pytest.approx(0.04)
    with pytest.raises(ValueError, match="aligned finite vectors"):
        neighborhood_calibration_diagnostics([0, 1], [0.3], [0.2, 0.4])


def test_frozen_split_evaluation_retains_splits_thresholds_and_reports_per_model_bands():
    rng = np.random.default_rng(4)
    fingerprints = rng.integers(0, 2, size=(60, 16)).tolist()
    labels = np.array([0, 1] * 30)
    frame = pd.DataFrame(
        {
            "fingerprint": fingerprints,
            "active": labels,
            "murcko_scaffold": [f"scaffold-{i}" for i in range(60)],
            "mol_wt": rng.uniform(200, 500, 60),
            "log_p": rng.uniform(-1, 5, 60),
            "tpsa": rng.uniform(20, 140, 60),
            "num_h_donors": rng.integers(0, 5, 60),
            "num_h_acceptors": rng.integers(1, 10, 60),
            "num_rotatable_bonds": rng.integers(0, 10, 60),
            "ring_count": rng.integers(1, 5, 60),
            "fraction_csp3": rng.uniform(0.1, 0.9, 60),
            "qed": rng.uniform(0.2, 0.9, 60),
        }
    )
    split = {"train": list(range(40)), "val": list(range(40, 50)), "test": list(range(50, 60))}
    original = {"scaffold": {"42": {key: list(value) for key, value in split.items()}}}
    result = evaluate_frozen_scaffold_neighborhoods(frame, original)
    assert result["split_file_unchanged"] is True
    assert result["fit_calibrator"] is False
    assert result["decision_threshold_changed"] is False
    assert result["per_seed"]["42"]["n_test"] == 10
    for model in ("baseline", "contender"):
        diagnostics = result["per_seed"]["42"]["models"][model][
            "reliability_by_nearest_training_similarity"
        ]
        assert diagnostics["prediction_calibration_fit"] is False
        assert sum(band["n"] for band in diagnostics["bands"]) == 10
    assert original["scaffold"]["42"] == split


@pytest.mark.parametrize(
    "train,test",
    [
        ([0.5], [2]),
        ([True], [2]),
        ([[0]], [2]),
        ([0, 0], [2]),
        ([0, 1], [1, 2]),
        ([0], [2, 2]),
        ([-1], [2]),
        ([0], [3]),
    ],
)
def test_neighborhood_rejects_invalid_or_leaking_indexes(train, test):
    with pytest.raises(ValueError):
        maximum_training_tanimoto(np.array([[1, 0], [0, 1], [1, 1]]), train, test)


@pytest.mark.parametrize(
    "values",
    [
        [[1, np.nan], [0, 1]],
        [[1, np.inf], [0, 1]],
        [[1, 2], [0, 1]],
        [[1, -1], [0, 1]],
        [[], []],
        ["bits", "bits"],
    ],
)
def test_neighborhood_rejects_invalid_fingerprints(values):
    with pytest.raises(ValueError, match="binary matrix"):
        maximum_training_tanimoto(values, [0], [1])


def test_tanimoto_matches_rdkit_without_bit_count_overflow():
    from rdkit import DataStructs
    from rdkit.DataStructs.cDataStructs import ExplicitBitVect

    values = np.ones((3, 2048), dtype=np.uint8)
    values[1, :1024] = 0
    values[2, :256] = 0

    def bitvector(row):
        vector = ExplicitBitVect(len(row))
        for bit in np.flatnonzero(row):
            vector.SetBit(int(bit))
        return vector

    reference = [bitvector(row) for row in values]
    expected = max(DataStructs.TanimotoSimilarity(reference[2], row) for row in reference[:2])
    np.testing.assert_allclose(maximum_training_tanimoto(values, [0, 1], [2]), [expected])


def test_observed_rate_intervals_show_uncertainty_and_empty_bands():
    result = neighborhood_calibration_diagnostics([1], [0.99], [0.2])
    band = result["bands"][0]
    assert band["n_active"] == 1
    assert band["n_inactive"] == 0
    assert band["observed_active_rate_ci95"][0] < 0.3
    assert band["observed_active_rate_ci95"][1] == pytest.approx(1)
    assert result["bands"][1]["observed_active_rate_ci95"] is None


def test_frozen_split_refuses_overlap_incomplete_fractional_and_scaffold_leakage():
    from bench.neighborhood_report import _validate_scaffold_partition

    frame = pd.DataFrame({"murcko_scaffold": ["A", "A", "B", "C", "D", "E"]})
    valid = {"train": [0, 1], "val": [2, 3], "test": [4, 5]}
    result = _validate_scaffold_partition(frame, valid)
    np.testing.assert_array_equal(result["train"], [0, 1])
    for invalid in (
        {"train": [0, 1], "val": [1, 2], "test": [3, 4, 5]},
        {"train": [0, 1], "val": [2, 3], "test": [4]},
        {"train": [0.5, 1], "val": [2, 3], "test": [4, 5]},
        {"train": [0, 2], "val": [1, 3], "test": [4, 5]},
    ):
        with pytest.raises(ValueError):
            _validate_scaffold_partition(frame, invalid)


def _structure_fixture():
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    from atlas.features import compute_morgan_fingerprint

    smiles = ["c1ccccc1", "c1ccncc1", "CCO", "CCN"]
    mols = [Chem.MolFromSmiles(smile) for smile in smiles]
    frame = pd.DataFrame(
        {
            "canonical_smiles": smiles,
            "murcko_scaffold": [MurckoScaffold.MurckoScaffoldSmiles(mol=mol) for mol in mols],
            "active": [0, 1, 0, 1],
        }
    )
    fingerprints = np.asarray([compute_morgan_fingerprint(mol) for mol in mols], dtype=np.uint8)
    return frame, fingerprints


def test_report_reconstructs_missing_fingerprints_and_hashes_exact_scored_snapshots(
    tmp_path, monkeypatch
):
    import hashlib
    import json

    import bench.neighborhood_report as report

    frame, fingerprints = _structure_fixture()
    dataset = tmp_path / "dataset.parquet"
    splits = tmp_path / "splits.json"
    frame.to_parquet(dataset)
    splits.write_text(json.dumps({"scaffold": {"42": {}}}))
    dataset_bytes, split_bytes = dataset.read_bytes(), splits.read_bytes()

    def evaluate(scored_frame, split_document, **kwargs):
        np.testing.assert_array_equal(
            np.asarray(scored_frame["fingerprint"].tolist()), fingerprints
        )
        dataset.write_bytes(b"changed after parse")
        splits.write_text("{}")
        return {"per_seed": {}}

    monkeypatch.setattr(report, "evaluate_frozen_scaffold_neighborhoods", evaluate)
    result = report.run(dataset, splits, tmp_path / "report.json")
    assert result["inputs"]["benchmark_sha256"] == hashlib.sha256(dataset_bytes).hexdigest()
    assert result["inputs"]["split_sha256"] == hashlib.sha256(split_bytes).hexdigest()
    assert (
        result["inputs"]["fingerprint_source"]
        == "reconstructed_from_frozen_smiles_no_cache_written"
    )
    assert not dataset.with_suffix(".fingerprints.npy").exists()


def test_report_rejects_same_shape_sidecar_with_wrong_row_order_and_false_scaffold(tmp_path):
    import bench.neighborhood_report as report

    frame, fingerprints = _structure_fixture()
    dataset = tmp_path / "dataset.parquet"
    splits = tmp_path / "splits.json"
    frame.to_parquet(dataset)
    splits.write_text("{}")
    np.save(dataset.with_suffix(".fingerprints.npy"), fingerprints[::-1])
    with pytest.raises(ValueError, match="row order"):
        report.run(dataset, splits, tmp_path / "report.json")
    frame.loc[0, "murcko_scaffold"] = "C1CCCCC1"
    frame.to_parquet(dataset)
    with pytest.raises(ValueError, match="scaffold does not match"):
        report.run(dataset, splits, tmp_path / "report.json")


def test_report_cannot_overwrite_frozen_input_or_configuration(tmp_path):
    import bench.neighborhood_report as report

    dataset = tmp_path / "dataset.parquet"
    splits = tmp_path / "splits.json"
    config = tmp_path / "config.yaml"
    for output in (dataset, splits, config, dataset.with_suffix(".fingerprints.npy")):
        with pytest.raises(ValueError, match="must not overwrite"):
            report.run(dataset, splits, output, config_path=config)
