"""Exact comparison graphs, not optimized medicines or human oral predictions."""

import json
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import Crippen, Descriptors, Draw, Lipinski, rdDepictor, rdMolDescriptors

ROOT = Path(__file__).resolve().parent
DISCOVERY = "https://pubmed.ncbi.nlm.nih.gov/37421886/"
METABOLISM = "https://pubmed.ncbi.nlm.nih.gov/41588687/"


@dataclass(frozen=True)
class Specification:
    identifier: str
    smiles: str
    status: str
    question: str
    source_url: str | None = None


SPECIFICATIONS = (
    Specification("SLU_PP_915", "OB(O)c1cccc(-c2ccc(C(=O)Nc3ccccc3F)s2)c1",
                  "published_parent_10s", "Intact parent reference", DISCOVERY),
    Specification("published_10q", "OB(O)c1cccc(-c2ccc(C(=O)Nc3cccc(F)c3)s2)c1",
                  "published_ERR_comparator", "Aniline fluorine position versus isoform response", DISCOVERY),
    Specification("published_10r", "OB(O)c1ccc(-c2ccc(C(=O)Nc3cccc(F)c3)s2)cc1",
                  "published_ERR_comparator", "Boronic-acid position versus isoform response", DISCOVERY),
    Specification("product_M1", "O=C(O)c1ccc(-c2cccc(O)c2)s1",
                  "published_reference_confirmed_transformation_product",
                  "Combined hydrolysis and boronic-acid-to-phenol conversion; activity unknown here", METABOLISM),
    Specification("product_M3", "OB(O)c1cccc(-c2ccc(C(=O)O)s2)c1",
                  "published_reference_confirmed_transformation_product",
                  "Amide hydrolysis product; activity unknown here", METABOLISM),
    Specification("product_M4", "Oc1cccc(-c2ccc(C(=O)Nc3ccccc3F)s2)c1",
                  "published_reference_confirmed_transformation_product",
                  "Boronic-acid-to-phenol product; activity unknown here", METABOLISM),
    Specification("hypothesis_aniline_2_4_difluoro", "OB(O)c1cccc(-c2ccc(C(=O)Nc3ccc(F)cc3F)s2)c1",
                  "defined_hypothesis_unvalidated_for_oral_improvement",
                  "Added para fluorine: test turnover and retained response; M6 carbon is not mapped"),
    Specification("hypothesis_aniline_2_5_difluoro", "OB(O)c1cccc(-c2ccc(C(=O)Nc3cc(F)ccc3F)s2)c1",
                  "defined_hypothesis_unvalidated_for_oral_improvement",
                  "Added fluorine at an alternative site; distinguish positional effects and metabolic switching"),
    Specification("hypothesis_amide_N_methyl", "OB(O)c1cccc(-c2ccc(C(=O)N(C)c3ccccc3F)s2)c1",
                  "defined_hypothesis_unvalidated_for_oral_improvement",
                  "Amide sterics and donor loss versus hydrolysis and ERR response; not assumed more stable"),
    Specification("hypothesis_boronic_pinacol_ester", "CC1(C)OB(c2cccc(-c3ccc(C(=O)Nc4ccccc4F)s3)c2)OC1(C)C",
                  "defined_precursor_hypothesis_unvalidated_for_oral_improvement",
                  "Masked boronic acid versus intact-parent release, competing breakdown and retained activity"),
)


def build_molecule(spec):
    molecule = Chem.MolFromSmiles(spec.smiles)
    if molecule is None or len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError("Specification must yield one valid covalent graph")
    return molecule


def record(spec):
    mol = build_molecule(spec)
    return {
        "id": spec.identifier,
        "evidence_status": spec.status,
        "comparison_question": spec.question,
        "source_url": spec.source_url,
        "canonical_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
        "standard_inchi_key": Chem.MolToInchiKey(mol),
        "neutral_formula": rdMolDescriptors.CalcMolFormula(mol),
        "neutral_average_mass_Da": round(Descriptors.MolWt(mol), 3),
        "neutral_exact_mass_Da": round(rdMolDescriptors.CalcExactMolWt(mol), 8),
        "neutral_graph_formal_charge": Chem.GetFormalCharge(mol),
        "neutral_graph_descriptors": {
            "rdkit_tpsa_A2": round(rdMolDescriptors.CalcTPSA(mol), 3),
            "rdkit_tpsa_include_S_P_A2": round(rdMolDescriptors.CalcTPSA(mol, includeSandP=True), 3),
            "crippen_logP": round(Crippen.MolLogP(mol), 3),
            "HBD": Lipinski.NumHDonors(mol),
            "HBA": Lipinski.NumHAcceptors(mol),
            "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        },
        "better_oral_exposure": "not_established",
        "human_bioavailability": "not_established_in_this_study",
        "novelty": "not_assessed",
    }


def main():
    records = [record(s) for s in SPECIFICATIONS]
    data = {
        "status": "structure_definitions_not_a_validated_drug_design_result",
        "rdkit_version": rdBase.rdkitVersion,
        "representation": "neutral_covalent_graphs; salts and pH-dependent speciation excluded",
        "limits": [
            "No synthesis, biological testing, oral-exposure improvement or novelty claim for hypotheses.",
            "Transformation-product identity does not establish inactivity, safety or human abundance.",
            "LogP and topological descriptors do not measure logD, solubility, permeability or binding.",
            "Both RDKit default TPSA and sulfur-inclusive TPSA are recorded; descriptor conventions differ.",
            "SDF coordinates and images are 2D depictions, not resolved or simulated bound poses.",
        ],
        "records": records,
    }
    (ROOT / "structural_comparators.json").write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    writer = Chem.SDWriter(str(ROOT / "structural_comparators.sdf"))
    try:
        for s in SPECIFICATIONS:
            mol = build_molecule(s)
            rdDepictor.Compute2DCoords(mol)
            mol.SetProp("_Name", s.identifier)
            mol.SetProp("evidence_status", s.status)
            mol.SetProp("limits", "Graph definition only; no oral improvement or human-use validation")
            writer.write(mol)
    finally:
        writer.close()
    # Normalize RDKit's cosmetic data-header spaces; preserve fixed-width
    # molfile atom/bond lines and all other SDF representation details.
    sdf_path = ROOT / "structural_comparators.sdf"
    lines = sdf_path.read_text().splitlines()
    sdf_path.write_text("\n".join(line.rstrip() if line.startswith("> ") else line for line in lines) + "\n")
    groups = [
        (SPECIFICATIONS[:6], "references.png", ["SLU-PP-915 (10s)", "Published 10q", "Published 10r", "Product M1", "Product M3", "Product M4"]),
        (SPECIFICATIONS[6:], "hypotheses.png", ["2,4-difluoro hypothesis", "2,5-difluoro hypothesis", "Amide N-methyl hypothesis", "Pinacol ester precursor hypothesis"]),
    ]
    for specs, filename, labels in groups:
        mols = [build_molecule(s) for s in specs]
        for mol in mols:
            rdDepictor.Compute2DCoords(mol)
        Draw.MolsToGridImage(mols, molsPerRow=2, subImgSize=(620, 320), legends=labels).save(ROOT / filename)
    print(json.dumps({"graphs": len(records), "rdkit_version": rdBase.rdkitVersion, "validated_improved_drugs": 0}))


if __name__ == "__main__":
    main()
