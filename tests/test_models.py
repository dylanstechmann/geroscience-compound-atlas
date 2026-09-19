"""Tests for atlas models, entities, and evidence grade validation."""

import pytest
from pydantic import ValidationError

from atlas.graph import EvidenceGraph
from atlas.models import (
    Assay,
    Compound,
    Document,
    EvidenceEdge,
    EvidenceGrade,
    HallmarkSlug,
    Modality,
    Target,
)


def test_refuse_evidence_edge_without_grade():
    """Requirement: EvidenceEdge must refuse instantiation without an explicit EvidenceGrade."""
    # Omitted grade argument
    with pytest.raises(ValidationError) as exc_info:
        EvidenceEdge(
            compound_inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            relation="inhibits",
        )
    errors = exc_info.value.errors()
    assert any(err["loc"] == ("grade",) for err in errors)

    # Explicit None grade
    with pytest.raises(ValidationError) as exc_info:
        EvidenceEdge(
            compound_inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            relation="inhibits",
            grade=None,  # type: ignore[arg-type]
        )
    errors = exc_info.value.errors()
    assert any(err["loc"] == ("grade",) for err in errors)

    # Invalid string grade
    with pytest.raises(ValidationError) as exc_info:
        EvidenceEdge(
            compound_inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            relation="inhibits",
            grade="E99",  # type: ignore[arg-type]
        )
    errors = exc_info.value.errors()
    assert any(err["loc"] == ("grade",) for err in errors)


def test_valid_evidence_edge():
    """Verify successful EvidenceEdge creation with valid grade and attributes."""
    edge = EvidenceEdge(
        compound_inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        relation="inhibits",
        grade=EvidenceGrade.E2,
        target_id="CHEMBL2842",
        hallmark=HallmarkSlug.nutrient_sensing,
        provenance="ChEMBL assay CHEMBL829471",
        document_ids=["PMID:12345678"],
    )
    assert edge.grade == EvidenceGrade.E2
    assert edge.compound_inchikey == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    assert edge.hallmark == HallmarkSlug.nutrient_sensing
    assert edge.grade.value == "E2"


def test_compound_model_validation():
    """Verify Compound entity validation and InChIKey format checking."""
    # Valid compound
    compound = Compound(
        name="Rapamycin",
        synonyms=["Sirolimus", "AY-22989"],
        cid=5284616,
        inchikey="ZZSNKZQZMQGXPY-UHFFFAOYSA-N",
        canonical_smiles="CC1CCC2CC(=O)C...",
        modality=Modality.small_molecule,
    )
    assert compound.name == "Rapamycin"
    assert compound.cid == 5284616
    assert compound.modality == Modality.small_molecule

    # Invalid InChIKey (wrong character structure)
    with pytest.raises(ValidationError):
        Compound(
            name="InvalidMol",
            inchikey="NOT-A-VALID-INCHIKEY",
            canonical_smiles="CC",
        )


def test_target_assay_document_models():
    """Verify basic instantiation of supporting schema models."""
    target = Target(
        target_id="CHEMBL2842", gene_symbol="MTOR", protein_family="PI3K-related kinase"
    )
    assert target.gene_symbol == "MTOR"

    assay = Assay(
        assay_id="CHEMBL829471", assay_type="binding", organism="Homo sapiens", standard_type="IC50"
    )
    assert assay.standard_type == "IC50"

    doc = Document(pmid="20083675", year=2009, species="Mus musculus", endpoint_class="lifespan")
    assert doc.endpoint_class == "lifespan"


def test_evidence_graph_filtering():
    """Verify EvidenceGraph filtering by evidence grade."""
    graph = EvidenceGraph()
    inchikey = "ZZSNKZQZMQGXPY-UHFFFAOYSA-N"

    compound = Compound(
        name="Rapamycin",
        inchikey=inchikey,
        canonical_smiles="CC1CCC2...",
    )
    graph.add_compound(compound)

    e1_edge = EvidenceEdge(
        compound_inchikey=inchikey,
        relation="inhibits",
        grade=EvidenceGrade.E1,
        hallmark=HallmarkSlug.nutrient_sensing,
    )
    e4_edge = EvidenceEdge(
        compound_inchikey=inchikey,
        relation="tested_in",
        grade=EvidenceGrade.E4,
        hallmark=HallmarkSlug.nutrient_sensing,
    )

    graph.add_edge(e1_edge)
    graph.add_edge(e4_edge)

    # All edges
    all_edges = graph.get_edges_for_compound(inchikey)
    assert len(all_edges) == 2

    # Filtered by min_grade E3
    high_grade_edges = graph.get_edges_for_compound(inchikey, min_grade=EvidenceGrade.E3)
    assert len(high_grade_edges) == 1
    assert high_grade_edges[0].grade == EvidenceGrade.E4
