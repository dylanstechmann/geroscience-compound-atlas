"""Tests for curated evidence edges, hallmark validation, and graph integrity."""

from pathlib import Path

import pandas as pd

from atlas.graph import EvidenceGraph
from atlas.models import EvidenceEdge, EvidenceGrade, HallmarkSlug


def test_curated_evidence_file_integrity():
    """Verify that all rows in data/curated/curated_evidence.csv are valid EvidenceEdges."""
    curated_path = Path("data/curated/curated_evidence.csv")
    assert curated_path.exists(), "Curated evidence CSV must exist."

    df = pd.read_csv(curated_path)
    assert len(df) >= 50, f"Specification §9 requires 50-150 curated edges; got {len(df)}"

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
