import hashlib
import json

import httpx
import numpy as np
import pandas as pd
import pytest

from atlas.normalize import canonicalize_smiles
from bench.external_study import (
    document_bootstrap,
    exclusion_reason,
    fetch_collection,
    load_snapshot,
    metrics,
    qualify,
)


def assay():
    return {
        "assay_chembl_id": "A1",
        "document_chembl_id": "D1",
        "target_chembl_id": "CHEMBL2842",
        "assay_type": "B",
        "confidence_score": 9,
        "relationship_type": "D",
        "assay_organism": "Homo sapiens",
        "variant_sequence": None,
        "description": "mTOR inhibition",
    }


def activity():
    return {
        "activity_id": 1,
        "assay_chembl_id": "A1",
        "document_chembl_id": "D1",
        "target_chembl_id": "CHEMBL2842",
        "assay_type": "B",
        "standard_type": "IC50",
        "standard_relation": "=",
        "standard_units": "nM",
        "standard_flag": 1,
        "potential_duplicate": 0,
        "pchembl_value": "6",
        "standard_value": "1000",
        "canonical_smiles": "CCO",
        "molecule_chembl_id": "M1",
        "parent_molecule_chembl_id": "M1",
    }


@pytest.mark.parametrize(
    "field,value",
    [
        ("confidence_score", 8),
        ("relationship_type", "H"),
        ("assay_organism", None),
        ("assay_type", "F"),
        ("variant_sequence", {"mutation": "A1V"}),
        ("description", "FKBP12 binding"),
        ("description", "cell viability"),
    ],
)
def test_assay_scope_rejects(field, value):
    changed = assay()
    changed[field] = value
    assert exclusion_reason(activity(), changed) is not None


@pytest.mark.parametrize(
    "field,value",
    [
        ("standard_relation", ">"),
        ("standard_units", "uM"),
        ("standard_flag", 0),
        ("potential_duplicate", 1),
        ("data_validity_comment", "Outside typical range"),
        ("pchembl_value", "NaN"),
        ("standard_value", "-1"),
        ("pchembl_value", "8"),
        ("document_chembl_id", "other"),
    ],
)
def test_activity_scope_rejects(field, value):
    changed = activity()
    changed[field] = value
    assert exclusion_reason(changed, assay()) is not None


def test_study_and_connectivity_holdouts_and_median_assay_aggregation():
    original = [
        {
            "document_chembl_id": "OLD",
            "molecule_chembl_id": "OLDM",
            "parent_molecule_chembl_id": None,
        }
    ]
    key = canonicalize_smiles("CCN")[1]
    frozen = pd.DataFrame({"inchikey": [key]})
    one = activity()
    two = {**one, "activity_id": 2, "standard_value": "10", "pchembl_value": "8"}
    overlap = {**one, "activity_id": 3, "canonical_smiles": "CCN"}
    old_document = {**one, "activity_id": 4, "assay_chembl_id": "A2", "document_chembl_id": "OLD"}
    old_molecule = {**one, "activity_id": 5, "molecule_chembl_id": "OLDM"}
    external, counts = qualify(
        [one, two, overlap, old_document, old_molecule],
        [assay(), {**assay(), "assay_chembl_id": "A2", "document_chembl_id": "OLD"}],
        original,
        frozen,
    )
    assert len(external) == 1
    assert external.iloc[0]["pchembl_value"] == 7
    assert external.iloc[0]["activity_ids"] == [1, 2]
    assert (
        counts["original_document"]
        == counts["original_connectivity"]
        == counts["original_molecule"]
        == 1
    )


def test_separate_assays_are_not_averaged():
    second = {
        **activity(),
        "activity_id": 2,
        "assay_chembl_id": "A2",
        "standard_value": "10",
        "pchembl_value": "8",
    }
    external, _ = qualify(
        [activity(), second],
        [assay(), {**assay(), "assay_chembl_id": "A2"}],
        [],
        pd.DataFrame({"inchikey": []}),
    )
    assert sorted(external.pchembl_value) == [6, 8]


def test_one_class_metrics_are_undefined_not_invented():
    result = metrics(np.ones(3), np.full(3, 0.5))
    assert result["auroc"] is None and result["average_precision"] is None
    assert result["brier"] == 0.25


def test_cluster_uncertainty_uses_documents_and_is_deterministic():
    y = np.array([1, 1, 0])
    probabilities = np.array([0.0, 0.0, 0.0])
    baseline = np.full(3, 0.5)
    documents = np.array(["A", "A", "B"])
    result = document_bootstrap(y, probabilities, baseline, documents)
    assert result["difference"] == pytest.approx(5 / 12)
    assert result["secondary_equal_document"]["difference"] == 0.25
    assert result == document_bootstrap(y, probabilities, baseline, documents)
    assert document_bootstrap(y, probabilities, baseline, np.array(["A", "A", "A"])) is None


def test_complete_download_and_hash_verified_reload(tmp_path):
    def handler(request):
        offset = int(request.url.params.get("offset", "0"))
        return httpx.Response(
            200,
            json={
                "activities": [{"activity_id": offset + 1}],
                "page_meta": {
                    "offset": offset,
                    "total_count": 2,
                    "next": "/chembl/api/data/activity.json?offset=1" if offset == 0 else None,
                },
            },
        )

    with httpx.Client(transport=httpx.MockTransport(handler)) as client:
        rows, receipts = fetch_collection(client, "activity", "activities", {}, tmp_path)
    assert len(rows) == 2 and len(receipts) == 2
    (tmp_path / "status.json").write_text("{}")
    receipt = {
        "activity_count": 2,
        "assay_count": 0,
        "pages": receipts,
        "status_sha256": hashlib.sha256(b"{}").hexdigest(),
    }
    (tmp_path / "receipt.json").write_text(json.dumps(receipt))
    assert len(load_snapshot(tmp_path)[0]) == 2
    (tmp_path / receipts[0]["path"]).write_text("changed")
    with pytest.raises(ValueError, match="hash mismatch"):
        load_snapshot(tmp_path)


@pytest.mark.parametrize("mode", ["truncated", "duplicate", "count_changed", "foreign_next"])
def test_partial_duplicate_or_changed_pagination_fails(tmp_path, mode):
    def handler(request):
        offset = int(request.url.params.get("offset", "0"))
        next_page = "/chembl/api/data/activity.json?offset=1" if offset == 0 else None
        if mode == "foreign_next":
            next_page = "https://elsewhere.invalid/page"
        if mode == "truncated":
            next_page = None
        return httpx.Response(
            200,
            json={
                "activities": [{"activity_id": 1 if mode == "duplicate" else offset + 1}],
                "page_meta": {
                    "offset": offset,
                    "total_count": 3 if mode == "count_changed" and offset else 2,
                    "next": next_page,
                },
            },
        )

    with httpx.Client(transport=httpx.MockTransport(handler)) as client, pytest.raises(ValueError):
        fetch_collection(client, "activity", "activities", {}, tmp_path)


def test_fkbp_independent_is_not_fkbp_binding():
    changed = assay()
    changed["description"] = (
        "Inhibition of FKBP12-independent human recombinant mTOR with His6-S6K1 by DELFIA"
    )
    assert exclusion_reason(activity(), changed) is None
