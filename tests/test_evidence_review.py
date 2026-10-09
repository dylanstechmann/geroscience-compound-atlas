"""Claim-bound title screens must never become verified biological evidence."""

import importlib.util
import json
import re
from collections import Counter

import pandas as pd
import pytest

from atlas.chembl_join import build_curated_edges_frame
from atlas.evidence_review import MANIFEST, apply_evidence_review, load_citation_audit
from atlas.lookup import format_lookup_report, lookup_compound
from viz.dashboard import build_dashboard_html


def test_packaged_manifest_reproduces_from_receipts():
    spec = importlib.util.spec_from_file_location("review_builder", "tools/build_evidence_review.py")
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    assert json.loads(MANIFEST.read_text(encoding="utf-8")) == builder.build_manifest()


def test_all_current_claims_match_without_promoting_title_matches():
    curated = pd.read_csv("data/curated/curated_evidence.csv")
    reviewed = build_curated_edges_frame(curated)
    assert len(reviewed) == 57
    assert set(reviewed["citation_audit_status"]) == {"matched"}
    assert set(reviewed["citation_audit_date"]) == {"2026-10-07"}
    assert set(reviewed["grade"]) == {"E0"}
    assert set(reviewed["claim_support_status"]) == {"unreviewed"}
    assert Counter(reviewed["citation_audit_category"]) == {
        "title_unrelated": 28, "not_found": 3, "chembl_not_found": 2,
        "needs_reading": 4, "on_topic_label_differs": 3, "on_topic": 15, "no_identifier": 2,
    }
    # Both identifiers stay separate for a multi-citation row.
    combination = reviewed[reviewed["citation_audit_identifiers"].map(len) == 2]
    assert len(combination) == 4


@pytest.mark.parametrize("field", ["document_ids", "notes", "source_title", "endpoint", "inchikey"])
def test_changed_claim_cannot_inherit_audit_by_line_or_shared_identifier(field):
    edges = build_curated_edges_frame(pd.read_csv("data/curated/curated_evidence.csv"))
    edges.loc[edges.index[0], field] = "changed claim"
    edges.loc[edges.index[0], "grade"] = "E4"
    reviewed = apply_evidence_review(edges)
    assert reviewed.iloc[0]["citation_audit_status"] == "unmatched"
    assert reviewed.iloc[0]["grade"] == "E0"
    assert reviewed.iloc[1]["citation_audit_status"] == "matched"


def test_row_reordering_and_reapplication_preserve_claim_identity():
    edges = build_curated_edges_frame(pd.read_csv("data/curated/curated_evidence.csv"))
    reordered = edges.iloc[::-1].copy()
    pd.testing.assert_frame_equal(apply_evidence_review(reordered), reordered)


@pytest.mark.parametrize("contents", [None, "broken json", '{"schema_version":99}'])
def test_missing_or_invalid_audit_withholds_legacy_grades(tmp_path, contents):
    manifest = tmp_path / "audit.json"
    if contents is not None:
        manifest.write_text(contents, encoding="utf-8")
    legacy = pd.DataFrame([{"inchikey": "legacy", "grade": "E4",
                            "source_review_status": "source_reviewed_with_claim_limits"}])
    reviewed = apply_evidence_review(legacy, manifest_path=manifest)
    assert reviewed.iloc[0]["grade"] == "E0"
    assert reviewed.iloc[0]["curator_grade"] == "E4"
    assert reviewed.iloc[0]["citation_audit_status"] == "unavailable"


@pytest.mark.parametrize("tamper", ["claim_support_status", "duplicate", "category"])
def test_invalid_manifest_cannot_promote_a_claim(tmp_path, tamper):
    data = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if tamper == "duplicate":
        data["records"].append(data["records"][0])
    else:
        data["records"][0][tamper] = "verified"
    path = tmp_path / "audit.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    assert load_citation_audit(path)[1] == "unavailable"


def test_lookup_overlays_old_parquet_and_shows_citation_problem():
    # The cached source still contains six historical reviewed rows.
    cached = pd.read_parquet("artifacts/evidence_edges.parquet")
    assert "E4" in set(cached["grade"])
    result = lookup_compound("Rapamycin", edges_parquet="artifacts/evidence_edges.parquet")
    assert all(edge["grade"] == "E0" for edge in result["hallmark_edges"])
    assert any(edge["citation_audit_category"] == "title_unrelated" for edge in result["hallmark_edges"])
    report = format_lookup_report(result)
    assert "titles only" in report
    assert "historical assertion" in report
    assert "Title is about smoking status" in report
    assert "Review    : Source reviewed" not in report


def test_dashboard_filters_use_withheld_grades_and_show_identifier_findings(tmp_path):
    output = build_dashboard_html(output_html=tmp_path / "review.html")
    html = output.read_text(encoding="utf-8")
    match = re.search(r"const compoundsData = (.*?);\s*\n", html)
    assert match, "Locate the rendered compound-card payload"
    edges = [edge for card in json.loads(match[1]) for edge in card["edges"]]
    assert len(edges) == 57
    assert all(edge["grade"] == "E0" for edge in edges)
    assert sum(edge["catalog_grade"] == "E4" for edge in edges) == 4
    assert "Citation screen:" in html
    assert "escapeHtml(item.screen_note)" in html
    assert "Historical catalog review label" in html
