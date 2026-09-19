"""Molecular featurization, RDKit 2D physico-chemical descriptors, and Morgan fingerprints."""

from typing import Any

from rdkit import Chem
from rdkit.Chem import QED, Descriptors, Lipinski
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

# Morgan fingerprint generator (radius 2, 2048 bits = ECFP4)
_MORGAN_GEN_R2_2048 = GetMorganGenerator(radius=2, fpSize=2048)


def compute_rdkit_descriptors(mol_or_smiles: Chem.Mol | str) -> dict[str, float | None]:
    """Compute standard 2D physico-chemical descriptors for a molecule.

    Scientific Note (§4):
        QED and Lipinski rule-of-five descriptors represent historical priors fitted
        on approved small-molecule sets. They are systematically unkind to natural
        products, macrocycles, and peptides.

    Args:
        mol_or_smiles: RDKit Mol object or canonical SMILES string.

    Returns:
        Dictionary of descriptor values (mol_wt, log_p, tpsa, hbd, hba, rotb, rings, fsp3, qed).
    """
    if isinstance(mol_or_smiles, str):
        mol = Chem.MolFromSmiles(mol_or_smiles)
    else:
        mol = mol_or_smiles

    if mol is None:
        return {
            "mol_wt": None,
            "log_p": None,
            "tpsa": None,
            "num_h_donors": None,
            "num_h_acceptors": None,
            "num_rotatable_bonds": None,
            "ring_count": None,
            "fraction_csp3": None,
            "qed": None,
        }

    try:
        mol_wt = float(Descriptors.MolWt(mol))
        log_p = float(Descriptors.MolLogP(mol))
        tpsa = float(Descriptors.TPSA(mol))
        hbd = int(Lipinski.NumHDonors(mol))
        hba = int(Lipinski.NumHAcceptors(mol))
        rotb = int(Lipinski.NumRotatableBonds(mol))
        rings = int(Lipinski.RingCount(mol))
        fsp3 = float(Lipinski.FractionCSP3(mol))
        qed_val = float(QED.qed(mol))
    except (ValueError, TypeError, ZeroDivisionError, RuntimeError):
        mol_wt = float(Descriptors.MolWt(mol)) if mol else None
        log_p = None
        tpsa = None
        hbd = None
        hba = None
        rotb = None
        rings = None
        fsp3 = None
        qed_val = None

    return {
        "mol_wt": round(mol_wt, 3) if mol_wt is not None else None,
        "log_p": round(log_p, 3) if log_p is not None else None,
        "tpsa": round(tpsa, 2) if tpsa is not None else None,
        "num_h_donors": hbd,
        "num_h_acceptors": hba,
        "num_rotatable_bonds": rotb,
        "ring_count": rings,
        "fraction_csp3": round(fsp3, 3) if fsp3 is not None else None,
        "qed": round(qed_val, 4) if qed_val is not None else None,
    }


def compute_morgan_fingerprint(
    mol_or_smiles: Chem.Mol | str, radius: int = 2, n_bits: int = 2048
) -> list[int] | None:
    """Generate Morgan (ECFP) bit vector fingerprint.

    Args:
        mol_or_smiles: RDKit Mol object or SMILES.
        radius: Radius of circular neighborhood (default 2 -> ECFP4).
        n_bits: Length of bit vector (default 2048).

    Returns:
        List of binary integers (0 or 1), or None if molecule cannot be parsed.
    """
    if isinstance(mol_or_smiles, str):
        mol = Chem.MolFromSmiles(mol_or_smiles)
    else:
        mol = mol_or_smiles

    if mol is None:
        return None

    if radius == 2 and n_bits == 2048:
        fp = _MORGAN_GEN_R2_2048.GetFingerprint(mol)
    else:
        gen = GetMorganGenerator(radius=radius, fpSize=n_bits)
        fp = gen.GetFingerprint(mol)

    return [int(b) for b in fp.ToList()]


def featurize_compound_record(record: dict[str, Any]) -> dict[str, Any]:
    """Featurize a compound record with canonical SMILES and RDKit descriptors."""
    output = dict(record)
    smiles = record.get("canonical_smiles")
    if smiles:
        descriptors = compute_rdkit_descriptors(smiles)
        output.update(descriptors)
    else:
        output.update(compute_rdkit_descriptors(None))
    return output
