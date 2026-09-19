"""Pydantic data schemas for Geroscience Compound Atlas entities and evidence edges."""

from enum import Enum

from pydantic import BaseModel, Field, field_validator


class EvidenceGrade(str, Enum):
    """Graded evidence hierarchy for geroscience compound mechanisms (E0 to E4).

    E0: Name only / vendor copy / no structured assay.
    E1: In vitro assay, single paper or ChEMBL row.
    E2: Multiple consistent assays or clear target engagement.
    E3: In vivo functional endpoint (not just a clock) in a model organism.
    E4: Mammalian lifespan / healthspan or registered human outcome trial.
    """

    E0 = "E0"
    E1 = "E1"
    E2 = "E2"
    E3 = "E3"
    E4 = "E4"


class Modality(str, Enum):
    """Compound structural modality."""

    small_molecule = "small_molecule"
    peptide = "peptide"
    antibody = "antibody"
    other = "other"


class HallmarkSlug(str, Enum):
    """Controlled vocabulary for aging hallmarks."""

    genomic_instability = "genomic_instability"
    epigenetic_alteration = "epigenetic_alteration"
    telomere_attrition = "telomere_attrition"
    proteostasis = "proteostasis"
    autophagy = "autophagy"
    nutrient_sensing = "nutrient_sensing"
    mitochondrial_dysfunction = "mitochondrial_dysfunction"
    cellular_senescence = "cellular_senescence"
    stem_cell_exhaustion = "stem_cell_exhaustion"
    intercellular_communication = "intercellular_communication"
    chronic_inflammation = "chronic_inflammation"
    ecm_fibrosis = "ecm_fibrosis"
    dysbiosis = "dysbiosis"


class Compound(BaseModel):
    """Chemical compound entity with canonical identifiers."""

    name: str = Field(..., description="Primary common or systematic name")
    synonyms: list[str] = Field(default_factory=list, description="Alternative names / synonyms")
    cid: int | None = Field(default=None, description="PubChem Compound ID (CID)")
    inchikey: str = Field(..., description="Standard 27-character InChIKey")
    canonical_smiles: str = Field(..., description="Canonical SMILES string")
    inchi: str | None = Field(default=None, description="Standard InChI string")
    modality: Modality = Field(default=Modality.small_molecule, description="Molecule modality")

    @field_validator("inchikey")
    @classmethod
    def validate_inchikey_format(cls, v: str) -> str:
        v_stripped = v.strip().upper()
        # InChIKey is 27 characters: 14 letters, hyphen, 10 letters, hyphen, 1 letter
        parts = v_stripped.split("-")
        if len(parts) != 3 or len(parts[0]) != 14 or len(parts[1]) != 10 or len(parts[2]) != 1:
            raise ValueError(f"Invalid InChIKey format: {v}")
        return v_stripped


class Target(BaseModel):
    """Biological target entity (protein / gene)."""

    target_id: str = Field(..., description="UniProt Accession or ChEMBL Target ID")
    gene_symbol: str | None = Field(default=None, description="Approved gene symbol (e.g. MTOR)")
    protein_family: str | None = Field(default=None, description="Protein family classification")


class Assay(BaseModel):
    """Structured experimental assay record."""

    assay_id: str = Field(..., description="ChEMBL Assay ID or internal identifier")
    assay_type: str | None = Field(
        default=None, description="Type: binding, functional, ADMET, phenotypic"
    )
    organism: str | None = Field(default=None, description="Test organism / species")
    standard_type: str | None = Field(
        default=None, description="Standard activity type (IC50, Ki, EC50, etc.)"
    )


class Document(BaseModel):
    """Bibliographic source reference."""

    pmid: str | None = Field(default=None, description="PubMed ID")
    doi: str | None = Field(default=None, description="Digital Object Identifier")
    year: int | None = Field(default=None, description="Publication year")
    species: str | None = Field(default=None, description="Target species investigated")
    endpoint_class: str | None = Field(default=None, description="Endpoint classification")


class EvidenceEdge(BaseModel):
    """Directed relation linking a compound to a target or hallmark with mandatory evidence grade."""

    compound_inchikey: str = Field(..., description="InChIKey of the subject compound")
    relation: str = Field(
        ...,
        description="Relation predicate: binds, inhibits, agonizes, degrades, phenotypic, claimed_for_hallmark, tested_in",
    )
    grade: EvidenceGrade = Field(
        ...,
        description="Mandatory evidence grade E0-E4. Edges without an explicit grade are rejected.",
    )
    target_id: str | None = Field(default=None, description="Target identifier if applicable")
    hallmark: HallmarkSlug | None = Field(default=None, description="Hallmark slug if applicable")
    provenance: str | None = Field(default=None, description="Data provenance / source")
    document_ids: list[str] = Field(default_factory=list, description="Referenced PMIDs or DOIs")
