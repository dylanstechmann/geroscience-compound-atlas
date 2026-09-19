"""Tests for RDKit molecular descriptor computation and Morgan fingerprints."""

from atlas.features import (
    compute_morgan_fingerprint,
    compute_rdkit_descriptors,
    featurize_compound_record,
)


def test_compute_rdkit_descriptors_valid_molecule():
    """Verify descriptor values on a well-characterized reference molecule (Metformin)."""
    metformin_smiles = "CN(C)C(=N)NC(=N)N"
    desc = compute_rdkit_descriptors(metformin_smiles)

    assert desc["mol_wt"] is not None
    assert 125.0 < desc["mol_wt"] < 135.0
    assert desc["num_h_donors"] == 4
    assert desc["num_h_acceptors"] == 2
    assert desc["qed"] is not None
    assert 0.0 <= desc["qed"] <= 1.0


def test_compute_rdkit_descriptors_invalid_molecule():
    """Verify that unparseable input returns None for all descriptor fields without throwing exceptions."""
    desc = compute_rdkit_descriptors("INVALID_SMILES")
    for key, val in desc.items():
        assert val is None, f"Expected None for {key}, got {val}"


def test_compute_morgan_fingerprint():
    """Verify 2048-bit Morgan fingerprint generation."""
    smiles = "CCO"  # Ethanol
    fp = compute_morgan_fingerprint(smiles, radius=2, n_bits=2048)

    assert fp is not None
    assert len(fp) == 2048
    assert sum(fp) > 0  # Should have set bits
    assert all(b in (0, 1) for b in fp)

    # Invalid smiles should return None
    assert compute_morgan_fingerprint("INVALID") is None


def test_featurize_compound_record():
    """Verify dictionary record featurization helper."""
    record = {
        "raw_name": "Quercetin",
        "canonical_smiles": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
        "cid": 5280343,
    }
    featurized = featurize_compound_record(record)
    assert featurized["mol_wt"] > 300.0
    assert featurized["num_h_donors"] == 5
    assert featurized["ring_count"] == 3
