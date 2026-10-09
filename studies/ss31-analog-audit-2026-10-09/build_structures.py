"""Exact graph bookkeeping for published comparators and untested hypotheses.

Neutral covalent structures exclude salts. RDKit descriptors are not oral PK,
biological potency, human-use suitability, cost estimates or optimized scores.
"""

import json
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdDepictor, rdMolDescriptors

ROOT = Path(__file__).resolve().parent
SAR_URL = "https://elifesciences.org/articles/75531"
PATENT_URL = "https://patents.google.com/patent/US20240108740A1/en"
SIDE_CHAINS = {
    "Arg": "CCCNC(=N)N",
    "Dmt": "Cc1c(C)cc(O)cc1C",
    "Tyr": "Cc1ccc(O)cc1",
    "Phe": "Cc1ccccc1",
    "Trp": "Cc1c[nH]c2ccccc12",
    "Lys": "CCCCN",
}


@dataclass(frozen=True)
class Specification:
    identifier: str
    residues: tuple[tuple[str, str], ...]
    status: str
    question: str
    acetylated: bool = False
    n_methyl_positions: tuple[int, ...] = ()


SPECIFICATIONS = (
    Specification("SS31", (("Arg", "D"), ("Dmt", "L"), ("Lys", "L"), ("Phe", "L")),
                  "published_reference", "Reference for every matched comparison"),
    Specification("SPN4", (("Arg", "D"), ("Tyr", "L"), ("Lys", "L"), ("Phe", "L")),
                  "published_cell_and_membrane_comparator",
                  "Specialty-residue substitution with endpoint-specific measured tradeoffs"),
    Specification("SS20", (("Phe", "L"), ("Arg", "D"), ("Phe", "L"), ("Lys", "L")),
                  "published_cell_and_membrane_comparator",
                  "Different register and aromatic chemistry; not an SS31 substitute"),
    Specification("SPN10", (("Trp", "L"), ("Arg", "L"), ("Trp", "L"), ("Lys", "L")),
                  "published_cell_and_membrane_comparator",
                  "Published functional comparator; human oral PK unknown here"),
    Specification("SS31_terminal_D_Phe", (("Arg", "D"), ("Dmt", "L"), ("Lys", "L"), ("Phe", "D")),
                  "research_hypothesis_unvalidated_here",
                  "C-terminal stereochemistry versus retained function and parent clearance"),
    Specification("SS31_N_acetyl", (("Arg", "D"), ("Dmt", "L"), ("Lys", "L"), ("Phe", "L")),
                  "research_hypothesis_unvalidated_here",
                  "Terminal cap versus nominal charge and mitochondrial activity", True),
    Specification("SPN10_D_Arg", (("Trp", "L"), ("Arg", "D"), ("Trp", "L"), ("Lys", "L")),
                  "research_hypothesis_unvalidated_here",
                  "One internal stereocenter inversion; proteolysis and activity both unknown"),
    Specification("SPN10_terminal_D_Lys", (("Trp", "L"), ("Arg", "L"), ("Trp", "L"), ("Lys", "D")),
                  "research_hypothesis_unvalidated_here",
                  "One terminal stereocenter inversion; degradation and activity both unknown"),
    Specification("SS31_Phe_backbone_N_methyl",
                  (("Arg", "D"), ("Dmt", "L"), ("Lys", "L"), ("Phe", "L")),
                  "patent_described_hypothesis_unvalidated_here",
                  "One backbone donor removed while preserving nominal +3 charge; oral PK unknown",
                  n_methyl_positions=(4,)),
)


def build_molecule(specification):
    if not specification.residues:
        raise ValueError("At least one residue is required")
    positions = specification.n_methyl_positions
    if (len(set(positions)) != len(positions)
            or any(type(position) is not int or position < 2
                   or position > len(specification.residues) for position in positions)):
        raise ValueError("N-methyl positions must be unique internal backbone nitrogens, numbered 2..n")
    parts = ["CC(=O)N" if specification.acetylated else "N"]
    for position, (name, configuration) in enumerate(specification.residues, start=1):
        if name not in SIDE_CHAINS or configuration not in ("L", "D"):
            raise ValueError("Only the explicitly supported residues and L/D configurations are valid")
        if position > 1:
            parts.append("N(C)" if position in positions else "N")
        alpha = "[C@@H]" if configuration == "L" else "[C@H]"
        parts.append(f"{alpha}({SIDE_CHAINS[name]})C(=O)")
    parts.append("N")
    molecule = Chem.MolFromSmiles("".join(parts))
    if molecule is None:
        raise ValueError("Could not construct a valid molecular graph")
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    observed = [label for _, label in Chem.FindMolChiralCenters(molecule)]
    # All supported side chains yield L=S and D=R at their alpha carbons.
    expected = ["S" if configuration == "L" else "R"
                for _, configuration in specification.residues]
    if observed != expected:
        raise ValueError(f"Stereochemical identity mismatch: {observed} != {expected}")
    return molecule


def record(specification):
    molecule = build_molecule(specification)
    sequence = "-".join(
        ("NalphaMe-" if position in specification.n_methyl_positions else "")
        + f"{configuration}-{name}"
        for position, (name, configuration) in enumerate(specification.residues, start=1))
    sequence = ("Ac-" if specification.acetylated else "H-") + sequence + "-NH2"
    return {
        "id": specification.identifier,
        "sequence_specification": sequence,
        "evidence_status": specification.status,
        "research_question": specification.question,
        "source_url": (PATENT_URL if specification.status.startswith("patent") else
                       SAR_URL if specification.status.startswith("published") else None),
        "isomeric_smiles": Chem.MolToSmiles(molecule, isomericSmiles=True),
        "standard_inchi_key": Chem.MolToInchiKey(molecule),
        "neutral_formula": rdMolDescriptors.CalcMolFormula(molecule),
        "neutral_average_mass_da": round(Descriptors.MolWt(molecule), 3),
        "neutral_graph_formal_charge": Chem.GetFormalCharge(molecule),
        "nominal_charge_near_physiological_ph": 2 if specification.acetylated else 3,
        "alpha_cip_in_sequence_order": [label for _, label in Chem.FindMolChiralCenters(molecule)],
        "neutral_graph_descriptors": {
            "tpsa_angstrom_squared": round(rdMolDescriptors.CalcTPSA(molecule), 3),
            "crippen_logp": round(Crippen.MolLogP(molecule), 3),
            "h_bond_donors": Lipinski.NumHDonors(molecule),
            "h_bond_acceptors": Lipinski.NumHAcceptors(molecule),
            "rotatable_bonds_rdkit": Lipinski.NumRotatableBonds(molecule),
        },
        "human_oral_bioavailability": "not_established_in_this_review",
        "improved_human_half_life": "not_established",
        "manufacturing_cost": "not_measured",
    }


def main():
    records = [record(specification) for specification in SPECIFICATIONS]
    output = {
        "status": "exact_graph_bookkeeping_not_a_drug_optimization_result",
        "rdkit_version": rdBase.rdkitVersion,
        "representation": "neutral_covalent_molecules_no_salts_or_physiological_protonation_model",
        "limits": [
            "Published membrane and cell results do not demonstrate human oral exposure.",
            "Hypotheses have no novelty, synthesis, biological activity or PK claim.",
            "Neutral descriptors are not logD, permeability or mitochondrial-access measurements.",
            "A valid molecular graph is not an analytical certificate for a physical sample.",
            "SDF coordinates are a 2D depiction, not a modeled bioactive conformation.",
        ],
        "records": records,
    }
    (ROOT / "structural_comparators.json").write_text(
        json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    writer = Chem.SDWriter(str(ROOT / "structural_comparators.sdf"))
    try:
        for specification in SPECIFICATIONS:
            molecule = build_molecule(specification)
            rdDepictor.Compute2DCoords(molecule)
            molecule.SetProp("_Name", specification.identifier)
            molecule.SetProp("evidence_status", specification.status)
            molecule.SetProp("limits", "Graph specification only; no oral PK or biological validation")
            writer.write(molecule)
    finally:
        writer.close()
    print(json.dumps({"graphs": len(records), "rdkit_version": rdBase.rdkitVersion,
                      "biological_validation_of_hypotheses": "none"}))


if __name__ == "__main__":
    main()
