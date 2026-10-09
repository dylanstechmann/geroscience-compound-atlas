"""Terminal-modification matrices; no oral or pharmacological optimization."""

import json
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, Draw, Lipinski, rdDepictor, rdMolDescriptors

ROOT = Path(__file__).resolve().parent


@dataclass(frozen=True)
class Specification:
    identifier: str
    sequence: str
    acetylated: bool = False
    amidated: bool = False
    status: str = "defined_terminal_comparison_unvalidated"


SPECIFICATIONS = (
    Specification("MOTS_c", "MRWQEMGYIFYPRKLR", status="published_parent_identity_reference"),
    Specification("MOTS_c_N_acetyl", "MRWQEMGYIFYPRKLR", True),
    Specification("MOTS_c_C_amide", "MRWQEMGYIFYPRKLR", False, True),
    Specification("MOTS_c_both_caps", "MRWQEMGYIFYPRKLR", True, True),
    Specification("MOTS_c_K14Q", "MRWQEMGYIFYPRQLR", status="published_sequence_variant_defined_free_acid_graph"),
    Specification("Epitalon", "AEDG", status="published_parent_identity_reference"),
    Specification("N_acetyl_Epitalon", "AEDG", True, status="public_identity_reference_biological_improvement_unvalidated"),
    Specification("Epitalon_C_amide", "AEDG", False, True),
    Specification("Epitalon_both_caps", "AEDG", True, True),
)


def terminal_atom(mol, residue_number, name):
    matches = [a.GetIdx() for a in mol.GetAtoms()
               if a.GetPDBResidueInfo() is not None
               and a.GetPDBResidueInfo().GetResidueNumber() == residue_number
               and a.GetPDBResidueInfo().GetName().strip() == name]
    if len(matches) != 1:
        raise ValueError("Terminal atom identity must be unique")
    return matches[0]


def build_molecule(spec):
    if not spec.sequence or any(letter not in "ARNDCQEGHILKMFPSTWYV" for letter in spec.sequence):
        raise ValueError("Use an explicitly defined standard all-L peptide sequence")
    mol = Chem.MolFromSequence(spec.sequence)
    if mol is None:
        raise ValueError("Sequence construction failed")
    editable = Chem.RWMol(mol)
    if spec.acetylated:
        terminal_n = terminal_atom(mol, 1, "N")
        carbonyl = editable.AddAtom(Chem.Atom(6))
        oxygen = editable.AddAtom(Chem.Atom(8))
        methyl = editable.AddAtom(Chem.Atom(6))
        editable.AddBond(terminal_n, carbonyl, Chem.BondType.SINGLE)
        editable.AddBond(carbonyl, oxygen, Chem.BondType.DOUBLE)
        editable.AddBond(carbonyl, methyl, Chem.BondType.SINGLE)
    if spec.amidated:
        # Residue-aware OXT selection avoids amidating Asp/Glu side chains.
        terminal_o = terminal_atom(mol, len(spec.sequence), "OXT")
        editable.GetAtomWithIdx(terminal_o).SetAtomicNum(7)
    result = editable.GetMol()
    Chem.SanitizeMol(result)
    Chem.AssignStereochemistry(result, cleanIt=True, force=True)
    return result


def nominal_charge(spec):
    """Ordinary R/K/E/D and terminal states; no compound-specific pKa model."""
    return (spec.sequence.count("R") + spec.sequence.count("K")
            - spec.sequence.count("E") - spec.sequence.count("D")
            + int(not spec.acetylated) - int(not spec.amidated))


def record(spec):
    mol = build_molecule(spec)
    return {
        "id": spec.identifier,
        "sequence_specification": ("Ac-" if spec.acetylated else "H-") + spec.sequence + ("-NH2" if spec.amidated else "-OH"),
        "evidence_status": spec.status,
        "canonical_isomeric_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
        "standard_inchi_key": Chem.MolToInchiKey(mol),
        "neutral_formula": rdMolDescriptors.CalcMolFormula(mol),
        "neutral_average_mass_Da": round(Descriptors.MolWt(mol), 3),
        "neutral_monoisotopic_mass_Da": round(rdMolDescriptors.CalcExactMolWt(mol), 8),
        "neutral_graph_formal_charge": Chem.GetFormalCharge(mol),
        "nominal_charge_near_physiological_pH": nominal_charge(spec),
        "stereocenter_count": len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
        "neutral_graph_descriptors": {
            "TPSA_A2": round(rdMolDescriptors.CalcTPSA(mol), 3),
            "HBD": Lipinski.NumHDonors(mol),
            "HBA": Lipinski.NumHAcceptors(mol),
        },
        "improved_human_oral_exposure": "not_established",
        "novelty": "not_assessed",
    }


def main():
    data = {
        "status": "exact_graph_and_charge_bookkeeping_not_drug_validation",
        "rdkit_version": rdBase.rdkitVersion,
        "representation": "all-L neutral covalent peptides; no salts, pH-speciation or conformational ensemble",
        "limits": [
            "Acetylation is a covalent cap, not an acetate counterion.",
            "Nominal charge is a residue-state convention, not a measured pKa or permeability.",
            "MOTS-c K14Q free-acid graph is defined here; original study preparation terminal identity is not independently verified.",
            "Only graph identity, not efficacy, is verified against PubChem references.",
            "SDF coordinates and images are two-dimensional depictions, not structures at a target.",
            "Caps preserve residue stereochemistry but need independent functional and exposure tests.",
        ],
        "records": [record(s) for s in SPECIFICATIONS],
    }
    (ROOT / "structural_comparators.json").write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    sdf_path = ROOT / "structural_comparators.sdf"
    writer = Chem.SDWriter(str(sdf_path))
    try:
        for s in SPECIFICATIONS:
            mol = build_molecule(s)
            rdDepictor.Compute2DCoords(mol)
            mol.SetProp("_Name", s.identifier)
            mol.SetProp("evidence_status", s.status)
            mol.SetProp("limits", "Graph identity only; no demonstrated oral improvement")
            writer.write(mol)
    finally:
        writer.close()
    lines = sdf_path.read_text().splitlines()
    sdf_path.write_text("\n".join(line.rstrip() if line.startswith("> ") else line for line in lines) + "\n")
    # The small Epitalon matrix is readable in 2D; long MOTS-c graphs stay in SDF.
    epis = SPECIFICATIONS[5:]
    labels = ["H-AEDG-OH (parent)", "Ac-AEDG-OH", "H-AEDG-NH2", "Ac-AEDG-NH2"]
    Draw.MolsToGridImage([build_molecule(s) for s in epis], molsPerRow=2, subImgSize=(620, 320), legends=labels).save(ROOT / "epitalon_caps.png")
    print(json.dumps({"graphs": len(SPECIFICATIONS), "rdkit_version": rdBase.rdkitVersion, "validated_oral_improvements": 0}))


if __name__ == "__main__":
    main()
