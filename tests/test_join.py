"""Tests for curated evidence edges, hallmark validation, and graph integrity."""

from pathlib import Path

import pandas as pd

from atlas.graph import EvidenceGraph
from atlas.models import EvidenceEdge, EvidenceGrade, HallmarkSlug
from atlas.chembl_join import build_curated_edges_frame, run_phase2_pipeline, validate_curated_evidence


def test_curated_evidence_file_integrity():
    """Verify that all rows in data/curated/curated_evidence.csv are valid EvidenceEdges."""
    curated_path = Path("data/curated/curated_evidence.csv")
    assert curated_path.exists(), "Curated evidence CSV must exist."

    df = pd.read_csv(curated_path)
    assert len(df) >= 50, f"Specification §9 requires 50-150 curated edges; got {len(df)}"
    compounds = pd.read_parquet("artifacts/compounds.parquet")
    validate_curated_evidence(df, compounds)

    graph = EvidenceGraph()
    valid_grades = {e.value for e in EvidenceGrade}
    valid_hallmarks = {h.value for h in HallmarkSlug}

    for idx, row in df.iterrows():
        # Grade must be valid
        grade_str = str(row["grade"]).strip().upper()
        assert grade_str in valid_grades, f"Row {idx} has invalid grade: {grade_str}"

        # Hallmark must be in controlled vocabulary or None
        hallmark_str = str(row["hallmark"]).strip().lower() if pd.notna(row["hallmark"]) else None
        if hallmark_str and hallmark_str != "none":
            assert hallmark_str in valid_hallmarks, (
                f"Row {idx} has unknown hallmark slug: {hallmark_str}"
            )
            h_enum = HallmarkSlug(hallmark_str)
        else:
            h_enum = None

        doc_raw = str(row.get("document_ids", "")).strip()
        docs = (
            [d.strip() for d in doc_raw.split(";") if d.strip()]
            if doc_raw and doc_raw != "None"
            else []
        )

        # Instantiation through Pydantic model must succeed
        edge = EvidenceEdge(
            compound_inchikey=str(row.get("inchikey", "")),
            relation=str(row["relation"]),
            grade=EvidenceGrade(grade_str),
            target_id=str(row.get("target_id", "")) if pd.notna(row.get("target_id")) else None,
            hallmark=h_enum,
            provenance=str(row.get("provenance", "")),
            document_ids=docs,
        )
        graph.add_edge(edge)

    assert len(graph.edges) == len(df)


def test_offline_curated_edge_materialization_keeps_provenance_and_combination_scope():
    curated = pd.read_csv("data/curated/curated_evidence.csv")
    compounds = pd.read_parquet("artifacts/compounds.parquet")
    validate_curated_evidence(curated, compounds)
    edges = build_curated_edges_frame(curated)
    assert len(edges) == len(curated)
    assert {"source_url", "source_title", "species", "endpoint", "effect_estimate", "source_review_status"}.issubset(edges)
    missing_sources = edges[edges["source_url"].fillna("").str.strip().eq("")]
    assert set(missing_sources["evidence_basis"]) == {"unverified_vendor_claim"}
    assert set(missing_sources["source_review_status"]) == {"unverified_vendor_or_gray_market_claim"}
    combo = edges[edges["provenance"] == "D+Q mouse clearance"]
    assert set(combo["intervention_components"]) == {"Dasatinib + quercetin"}
    assert set(combo["attribution_scope"]) == {"combination"}


def test_unreviewed_claim_grades_are_withheld_until_source_review():
    curated = pd.read_csv("data/curated/curated_evidence.csv")
    compounds = pd.read_parquet("artifacts/compounds.parquet")
    validate_curated_evidence(curated, compounds)
    edges = build_curated_edges_frame(curated)

    unreviewed = edges[edges["source_review_status"] == "claim_support_unreviewed"]
    reviewed = edges[edges["source_review_status"] == "source_reviewed_with_claim_limits"]
    vendor = edges[edges["source_review_status"] == "unverified_vendor_or_gray_market_claim"]
    assert len(unreviewed) == 49
    assert len(reviewed) == 6
    assert len(vendor) == 2
    assert set(unreviewed["grade"]) == {"E0"}
    assert set(vendor["grade"]) == {"E0"}
    assert set(unreviewed["curator_grade"]) == {"E1", "E2", "E3"}
    assert (reviewed["grade"] == reviewed["curator_grade"]).all()
    assert set(reviewed[reviewed["grade"] == "E4"]["evidence_basis"]) == {"mammalian_lifespan"}
    assert "curator_grade" in curated.columns
    assert set(curated.loc[curated["source_review_status"] == "claim_support_unreviewed", "grade"]) == {"E0"}

    bad = curated.copy()
    row_index = bad.index[bad["source_review_status"] == "claim_support_unreviewed"][0]
    bad.loc[row_index, "grade"] = bad.loc[row_index, "curator_grade"]
    try:
        validate_curated_evidence(bad, compounds)
    except ValueError as exc:
        assert "conflicts with review status" in str(exc)
    else:
        raise AssertionError("Unreviewed curator grades must not be promoted in the active CSV")


def test_high_grade_curated_interventions():
    """Verify that E4 lifespan interventions are properly registered."""
    df = pd.read_csv("data/curated/curated_evidence.csv")
    e4_rows = df[df["grade"] == "E4"]

    assert len(e4_rows) >= 3, (
        "Expected multiple E4 validated interventions (e.g. Rapamycin, Acarbose, Metformin, 17a-estradiol)"
    )
    e4_compounds = set(e4_rows["compound_name"])
    assert "Rapamycin" in e4_compounds
    assert "17alpha-estradiol" in e4_compounds
    assert "Acarbose" in e4_compounds
    assert "Canagliflozin" in e4_compounds
    assert "Metformin" not in e4_compounds


def test_compound_ids_and_inchikeys_match_exact_snapshot_names():
    curated = pd.read_csv("data/curated/curated_evidence.csv")
    compounds = pd.read_parquet("artifacts/compounds.parquet").set_index("raw_name")
    assert curated["compound_id"].notna().all()
    for _, row in curated.iterrows():
        resolved = compounds.loc[row["compound_name"]]
        assert row["compound_id"] == f"pubchem:{int(resolved['cid'])}"
        assert row["inchikey"] == resolved["inchikey"]


def test_target_id_symbol_conflicts_fail_closed():
    curated = pd.read_csv("data/curated/curated_evidence.csv")
    compounds = pd.read_parquet("artifacts/compounds.parquet")
    bad = curated.copy()
    bad.loc[bad.index[0], "target_symbol"] = "PIK3CA"
    try:
        validate_curated_evidence(bad, compounds)
    except ValueError as exc:
        assert "target conflict" in str(exc)
    else:
        raise AssertionError("A CHEMBL2842 / PIK3CA contradiction must be rejected")


def test_dangling_compound_and_e4_observational_claim_fail_closed():
    curated = pd.read_csv("data/curated/curated_evidence.csv")
    compounds = pd.read_parquet("artifacts/compounds.parquet")
    bad_identity = curated.copy()
    bad_identity.loc[bad_identity.index[0], "compound_id"] = "pubchem:999999999"
    try:
        validate_curated_evidence(bad_identity, compounds)
    except ValueError as exc:
        assert "unresolved compound_id" in str(exc)
    else:
        raise AssertionError("A dangling compound edge must be rejected")

    bad_grade = curated.copy()
    metformin = bad_grade.index[
        (bad_grade["compound_name"] == "Metformin")
        & (bad_grade["provenance"] == "Human epidemiological & clinical")
    ][0]
    bad_grade.loc[metformin, "grade"] = "E4"
    try:
        validate_curated_evidence(bad_grade, compounds)
    except ValueError as exc:
        assert "E4 without" in str(exc)
    else:
        raise AssertionError("Observational association must not pass as E4")


def test_combination_evidence_is_labeled_as_combination():
    curated = pd.read_csv("data/curated/curated_evidence.csv")
    combo = curated[curated["provenance"] == "D+Q mouse clearance"]
    assert set(combo["compound_name"]) == {"Dasatinib", "Quercetin"}
    assert set(combo["attribution_scope"]) == {"combination"}
    assert set(combo["intervention_components"]) == {"Dasatinib + quercetin"}


def test_pipeline_rejects_identity_conflict_before_client_construction(tmp_path, monkeypatch):
    curated = pd.read_csv("data/curated/curated_evidence.csv")
    curated.loc[curated.index[0], "target_symbol"] = "PIK3CA"
    bad_path = tmp_path / "bad_curated.csv"
    curated.to_csv(bad_path, index=False)

    def forbidden_client(*args, **kwargs):
        raise AssertionError("Network client must not be initialized before curated data validation")

    monkeypatch.setattr("atlas.chembl_join.ChEMBLClient", forbidden_client)
    try:
        run_phase2_pipeline(
            compounds_parquet="artifacts/compounds.parquet",
            curated_edges_csv=bad_path,
            cache_path=tmp_path / "cache.json",
            activities_out=tmp_path / "activities.parquet",
            targets_out=tmp_path / "targets.parquet",
            edges_out=tmp_path / "edges.parquet",
            figure_out=tmp_path / "figure.png",
        )
    except ValueError as exc:
        assert "target conflict" in str(exc)
    else:
        raise AssertionError("The pipeline must reject a target contradiction")
