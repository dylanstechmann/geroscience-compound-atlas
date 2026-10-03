"""Unit tests for single-compound lookup CLI and report formatting."""

import json
from pathlib import Path
import pytest

from atlas.lookup import lookup_compound, format_lookup_report, main as lookup_main
from atlas.cli import main as cli_main


def test_lookup_watchlist_compound_rapamycin():
    res = lookup_compound("Rapamycin")
    assert res["in_atlas"] is True
    assert "Sirolimus" in res["name"] or "Rapamycin" in res["name"]
    assert res["inchikey"] == "QFJCIRLUMZQUOT-HPLJOQBZSA-N"
    assert res["cid"] == 5284616
    assert len(res["hallmark_edges"]) > 0
    assert any("nutrient_sensing" in e["hallmark"] for e in res["hallmark_edges"])
    assert res["benchmark_prediction"]["available"] is True
    assert res["benchmark_prediction"]["prediction"] == "ABSTAIN"
    assert res["benchmark_prediction"]["prediction_status"] == "unvalidated_surrogate"
    assert res["benchmark_prediction"]["applicability_domain_status"] == "not_validated"
    assert res["benchmark_prediction"]["mtor_prob"] > 0.9


def test_lookup_watchlist_compound_case_insensitive_and_cid():
    res_lower = lookup_compound("metformin")
    assert res_lower["in_atlas"] is True
    assert res_lower["cid"] == 4091

    res_cid = lookup_compound("4091")
    assert res_cid["in_atlas"] is True
    assert res_cid["inchikey"] == res_lower["inchikey"]
    assert res_cid["benchmark_prediction"]["prediction"] == "ABSTAIN"


def test_lookup_by_inchikey():
    res = lookup_compound("QFJCIRLUMZQUOT-HPLJOQBZSA-N")
    assert res["in_atlas"] is True
    assert "Sirolimus" in res["name"] or "Rapamycin" in res["name"]


def test_lookup_custom_smiles():
    # Ethanol
    res = lookup_compound("CCO")
    assert res["in_atlas"] is False
    assert "CCO" in res["canonical_smiles"]
    assert res["benchmark_prediction"]["available"] is True
    assert res["benchmark_prediction"]["prediction"] == "ABSTAIN"


def test_lookup_invalid_compound_raises():
    with pytest.raises(ValueError, match="not found in atlas"):
        lookup_compound("NOT_A_VALID_COMPOUND_XYZ_12345")


def test_format_lookup_report():
    res = lookup_compound("Rapamycin")
    report = format_lookup_report(res)
    assert "GEROSCIENCE COMPOUND EVIDENCE CARD" in report
    assert "[IDENTIFIERS]" in report
    assert "[PHYSICOCHEMICAL PROPERTIES]" in report
    assert "[AGING HALLMARK ASSOCIATIONS" in report
    assert "[CHEMBL EXPERIMENTAL BIOACTIVITIES]" in report
    assert "[PREDICTIVE BENCHMARK INFERENCE (mTOR Kinase)]" in report
    assert "Uncalibrated surrogate score" in report
    assert "ABSTAIN" in report
    assert "does not establish binding" in report


def test_lookup_distinguishes_source_reviewed_edges_from_unreviewed_claims():
    result = lookup_compound("Rapamycin")
    unreviewed = [edge for edge in result["hallmark_edges"] if edge["source_review_status"] == "claim_support_unreviewed"]
    assert unreviewed
    assert all(edge["grade"] == "E0" for edge in unreviewed)
    assert all(edge["curator_grade"] in {"E1", "E2", "E3"} for edge in unreviewed)
    report = format_lookup_report(result)
    assert "Claim support unreviewed" in report
    assert "withheld and displayed as E0" in report


def test_invalid_surrogate_score_is_not_reported_as_available():
    class InvalidScorer:
        label_definition = "active = 1 if pChEMBL >= configured threshold"

        def score_molecule(self, mol):
            return {"valid": False, "mtor_prob": 0.0, "qed": 0.0, "has_pains": False, "reward": 0.0}

    result = lookup_compound("CCO", scorer=InvalidScorer())
    assert result["benchmark_prediction"]["available"] is False
    assert result["benchmark_prediction"]["prediction"] == "ABSTAIN"
    assert result["benchmark_prediction"]["prediction_status"] == "invalid_structure_or_score"


def test_lookup_cli_main_stdout_and_json(capsys):
    ret = lookup_main(["Metformin", "--json"])
    assert ret == 0
    captured = capsys.readouterr()
    parsed = json.loads(captured.out)
    assert parsed["cid"] == 4091
    assert parsed["benchmark_prediction"]["available"] is True


def test_atlas_cli_main_lookup_subcommand(tmp_path: Path):
    out_file = tmp_path / "card.txt"
    ret = cli_main(["lookup", "Rapamycin", "--out", str(out_file)])
    assert ret == 0
    assert out_file.exists()
    content = out_file.read_text(encoding="utf-8")
    assert "GEROSCIENCE COMPOUND EVIDENCE CARD" in content
