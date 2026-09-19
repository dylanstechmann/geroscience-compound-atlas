"""Tests for chemical structure normalization and SMILES canonicalization."""

import pytest

from atlas.normalize import canonicalize_smiles, smiles_to_inchikey


def test_canonicalize_equivalent_smiles_to_same_inchikey():
    """Requirement: Two equivalent SMILES representations must canonicalize to one identical InChIKey."""
    # Test case 1: Ethanol (different atom orderings)
    smiles_ethanol_1 = "CCO"
    smiles_ethanol_2 = "OCC"

    can_1, key_1 = canonicalize_smiles(smiles_ethanol_1)
    can_2, key_2 = canonicalize_smiles(smiles_ethanol_2)

    assert key_1 == key_2, f"Expected identical InChIKey, got {key_1} vs {key_2}"
    assert can_1 == can_2
    assert key_1 == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"

    # Test case 2: Benzene (aromatic vs Kekulé form)
    smiles_benzene_1 = "c1ccccc1"
    smiles_benzene_2 = "C1=CC=CC=C1"

    can_b1, key_b1 = canonicalize_smiles(smiles_benzene_1)
    can_b2, key_b2 = canonicalize_smiles(smiles_benzene_2)

    assert key_b1 == key_b2, f"Expected identical InChIKey, got {key_b1} vs {key_b2}"
    assert can_b1 == can_b2
    assert key_b1 == "UHOVQNZJYSORNB-UHFFFAOYSA-N"

    # Test case 3: Metformin (tautomer / syntax representations)
    smiles_metformin_1 = "CN(C)C(=N)NC(=N)N"
    smiles_metformin_2 = "NC(=N)NC(=N)N(C)C"

    _, key_m1 = canonicalize_smiles(smiles_metformin_1)
    _, key_m2 = canonicalize_smiles(smiles_metformin_2)

    assert key_m1 == key_m2


def test_smiles_to_inchikey_helper():
    """Verify smiles_to_inchikey helper produces correct InChIKey."""
    key = smiles_to_inchikey("CCO")
    assert key == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"


def test_invalid_smiles_raises_value_error():
    """Unparseable or nonsensical SMILES strings must raise ValueError."""
    invalid_examples = [
        "NOT_A_SMILES",
        "C1CC",  # unclosed ring
        "Xx[invalid]",
    ]
    for bad_smiles in invalid_examples:
        with pytest.raises(ValueError):
            canonicalize_smiles(bad_smiles)


def test_empty_smiles_raises_value_error():
    """Empty or whitespace SMILES input must raise ValueError."""
    with pytest.raises(ValueError):
        canonicalize_smiles("")

    with pytest.raises(ValueError):
        canonicalize_smiles("   ")
