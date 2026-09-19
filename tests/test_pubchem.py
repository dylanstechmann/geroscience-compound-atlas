"""Tests for PubChem client caching, resolution logic, and rate limiting."""

import json

from atlas.pubchem import PubChemClient


def test_pubchem_cache_hit(tmp_path):
    """Verify that cached entries are returned without network activity."""
    cache_file = tmp_path / "test_cache.json"
    fake_entry = {
        "raw_name": "Rapamycin",
        "cid": 5284616,
        "title": "Sirolimus",
        "inchikey": "ZZSNKZQZMQGXPY-UHFFFAOYSA-N",
        "canonical_smiles": "CC1CCC2CC(=O)C...",
        "resolved": True,
    }
    with open(cache_file, "w", encoding="utf-8") as f:
        json.dump({"rapamycin": fake_entry}, f)

    client = PubChemClient(cache_path=cache_file)
    res = client.resolve_name("Rapamycin")

    assert res is not None
    assert res["cid"] == 5284616
    assert res["title"] == "Sirolimus"
    assert res["resolved"] is True


def test_pubchem_cache_negative_result(tmp_path):
    """Verify that negative resolution results (not found) are cached."""
    cache_file = tmp_path / "test_cache.json"
    fake_negative = {
        "raw_name": "NonExistentCompound123",
        "resolved": False,
        "error": "PubChem compound not found (404)",
    }
    with open(cache_file, "w", encoding="utf-8") as f:
        json.dump({"nonexistentcompound123": fake_negative}, f)

    client = PubChemClient(cache_path=cache_file)
    res = client.resolve_name("NonExistentCompound123")

    assert res is not None
    assert res["resolved"] is False
    assert "error" in res


def test_pubchem_batch_resolve(tmp_path):
    """Verify batch resolution with cache."""
    cache_file = tmp_path / "test_cache.json"
    with open(cache_file, "w", encoding="utf-8") as f:
        json.dump(
            {
                "comp_a": {"raw_name": "comp_a", "cid": 1, "resolved": True},
                "comp_b": {"raw_name": "comp_b", "cid": 2, "resolved": True},
            },
            f,
        )

    client = PubChemClient(cache_path=cache_file)
    results = client.batch_resolve(["comp_a", "comp_b"])
    assert len(results) == 2
    assert results["comp_a"]["cid"] == 1
    assert results["comp_b"]["cid"] == 2
