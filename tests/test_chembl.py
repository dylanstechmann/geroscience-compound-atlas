"""Tests for ChEMBL REST client caching, queries, and rate-limiting."""

import json

from atlas.chembl import ChEMBLClient


def test_chembl_cache_hit_molecule(tmp_path):
    """Verify ChEMBLClient returns cached molecule mapping without network requests."""
    cache_file = tmp_path / "test_chembl_cache.json"
    cache_data = {
        "molecules": {"XZWYZXLIPXDOLR-UHFFFAOYSA-N": "CHEMBL1431"},
        "activities": {},
        "targets": {},
    }
    with open(cache_file, "w", encoding="utf-8") as f:
        json.dump(cache_data, f)

    client = ChEMBLClient(cache_path=cache_file)
    chembl_id = client.get_chembl_id_by_inchikey("XZWYZXLIPXDOLR-UHFFFAOYSA-N")

    assert chembl_id == "CHEMBL1431"


def test_chembl_cache_hit_activities(tmp_path):
    """Verify ChEMBLClient returns cached activity rows."""
    cache_file = tmp_path / "test_chembl_cache.json"
    cache_data = {
        "molecules": {},
        "activities": {
            "CHEMBL1431": [
                {
                    "activity_id": 12345,
                    "assay_chembl_id": "CHEMBL779287",
                    "assay_type": "B",
                    "target_chembl_id": "CHEMBL376",
                    "standard_type": "IC50",
                    "standard_value": "100.0",
                }
            ]
        },
        "targets": {},
    }
    with open(cache_file, "w", encoding="utf-8") as f:
        json.dump(cache_data, f)

    client = ChEMBLClient(cache_path=cache_file)
    acts = client.get_compound_activities("CHEMBL1431")

    assert len(acts) == 1
    assert acts[0]["assay_type"] == "B"
    assert acts[0]["standard_type"] == "IC50"


def test_chembl_cache_hit_targets(tmp_path):
    """Verify ChEMBLClient returns cached target details."""
    cache_file = tmp_path / "test_chembl_cache.json"
    cache_data = {
        "molecules": {},
        "activities": {},
        "targets": {
            "CHEMBL376": {
                "target_chembl_id": "CHEMBL376",
                "pref_name": "AMP-activated protein kinase catalytic subunit alpha-1",
                "target_type": "SINGLE PROTEIN",
                "organism": "Homo sapiens",
                "uniprot_accession": "Q13131",
            }
        },
    }
    with open(cache_file, "w", encoding="utf-8") as f:
        json.dump(cache_data, f)

    client = ChEMBLClient(cache_path=cache_file)
    target = client.get_target_details("CHEMBL376")

    assert target is not None
    assert target["pref_name"].startswith("AMP-activated")
    assert target["uniprot_accession"] == "Q13131"
