"""Chemical mutation and crossover operators for molecular generation."""

import copy
import logging
import random

from rdkit import Chem

logger = logging.getLogger(__name__)

# Common medicinal chemistry building blocks for fragment attachment
COMMON_FRAGMENTS = [
    "[CH3]",  # Methyl
    "[OH]",  # Hydroxyl
    "[NH2]",  # Amine
    "[F]",  # Fluoro
    "[Cl]",  # Chloro
    "C(=O)O",  # Carboxyl
    "C(=O)N",  # Amide
    "C1CC1",  # Cyclopropyl
    "C1CCOCC1",  # Morpholine
    "N1CCOCC1",  # Morpholinyl
]


def sanitize_molecule(mol: Chem.Mol | None) -> Chem.Mol | None:
    """Sanitize and validate an RDKit Mol object, returning None if invalid."""
    if mol is None:
        return None
    try:
        mol_copy = copy.deepcopy(mol)
        Chem.SanitizeMol(mol_copy)
        # Ensure it has at least 4 heavy atoms and no disconnects
        frags = Chem.GetMolFrags(mol_copy, asMols=True)
        if len(frags) != 1:
            # Pick largest fragment
            mol_copy = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        if mol_copy.GetNumHeavyAtoms() < 4:
            return None
        # Check that it converts to valid SMILES and back
        smi = Chem.MolToSmiles(mol_copy)
        if not smi:
            return None
        valid_mol = Chem.MolFromSmiles(smi)
        return valid_mol
    except (ValueError, RuntimeError, TypeError):
        return None


def mutate_atom_type(mol: Chem.Mol, rng: random.Random | None = None) -> Chem.Mol | None:
    """Randomly mutate the element of an allowable atom (e.g. C -> N, O, S, F)."""
    rng = rng or random
    mol_copy = copy.deepcopy(mol)
    atoms = [a for a in mol_copy.GetAtoms() if a.GetAtomicNum() in {6, 7, 8, 16}]  # C, N, O, S
    if not atoms:
        return None

    target_atom = rng.choice(atoms)
    current_z = target_atom.GetAtomicNum()
    candidates = [6, 7, 8, 9, 16, 17]  # C, N, O, F, S, Cl
    candidates = [z for z in candidates if z != current_z]
    new_z = rng.choice(candidates)

    target_atom.SetAtomicNum(new_z)
    return sanitize_molecule(mol_copy)


def mutate_add_atom_or_fragment(mol: Chem.Mol, rng: random.Random | None = None) -> Chem.Mol | None:
    """Attach a small functional group to a hydrogen-bearing heavy atom."""
    rng = rng or random
    mol_copy = copy.deepcopy(mol)
    eligible_atoms = [
        a.GetIdx()
        for a in mol_copy.GetAtoms()
        if a.GetTotalNumHs() > 0 and a.GetAtomicNum() in {6, 7, 8}
    ]
    if not eligible_atoms:
        return None

    attach_idx = rng.choice(eligible_atoms)
    frag_smi = rng.choice(COMMON_FRAGMENTS)
    frag = Chem.MolFromSmiles(frag_smi)
    if frag is None:
        return None

    # Combine molecules
    combined = Chem.CombineMols(mol_copy, frag)
    rw_mol = Chem.RWMol(combined)

    # Connect attach_idx in mol to first atom in fragment
    frag_start_idx = mol_copy.GetNumAtoms()
    rw_mol.AddBond(attach_idx, frag_start_idx, Chem.BondType.SINGLE)

    return sanitize_molecule(rw_mol.GetMol())


def mutate_remove_atom(mol: Chem.Mol, rng: random.Random | None = None) -> Chem.Mol | None:
    """Remove a peripheral/terminal atom if the molecule is large enough."""
    rng = rng or random
    if mol.GetNumHeavyAtoms() <= 5:
        return None

    mol_copy = copy.deepcopy(mol)
    terminal_atoms = [
        a.GetIdx() for a in mol_copy.GetAtoms() if a.GetDegree() == 1 and not a.IsInRing()
    ]
    if not terminal_atoms:
        return None

    del_idx = rng.choice(terminal_atoms)
    rw_mol = Chem.RWMol(mol_copy)
    rw_mol.RemoveAtom(del_idx)
    return sanitize_molecule(rw_mol.GetMol())


def mutate_bond_order(mol: Chem.Mol, rng: random.Random | None = None) -> Chem.Mol | None:
    """Modify acyclic bond order between single and double where chemically valid."""
    rng = rng or random
    mol_copy = copy.deepcopy(mol)
    bonds = [
        b
        for b in mol_copy.GetBonds()
        if not b.IsInRing() and b.GetBondType() in {Chem.BondType.SINGLE, Chem.BondType.DOUBLE}
    ]
    if not bonds:
        return None

    bond = rng.choice(bonds)
    if bond.GetBondType() == Chem.BondType.SINGLE:
        bond.SetBondType(Chem.BondType.DOUBLE)
    else:
        bond.SetBondType(Chem.BondType.SINGLE)

    return sanitize_molecule(mol_copy)


def crossover(
    parent1: Chem.Mol, parent2: Chem.Mol, rng: random.Random | None = None
) -> Chem.Mol | None:
    """Single-cut crossover between two parent molecules at acyclic single bonds."""
    rng = rng or random
    try:
        bonds1 = [
            b.GetIdx()
            for b in parent1.GetBonds()
            if not b.IsInRing() and b.GetBondType() == Chem.BondType.SINGLE
        ]
        bonds2 = [
            b.GetIdx()
            for b in parent2.GetBonds()
            if not b.IsInRing() and b.GetBondType() == Chem.BondType.SINGLE
        ]

        if not bonds1 or not bonds2:
            return None

        # Fragment parent 1
        b1 = rng.choice(bonds1)
        broken1 = Chem.FragmentOnBonds(parent1, [b1], dummyLabels=[(1, 1)])
        frags1 = Chem.GetMolFrags(broken1, asMols=True)
        if len(frags1) != 2:
            return None

        # Fragment parent 2
        b2 = rng.choice(bonds2)
        broken2 = Chem.FragmentOnBonds(parent2, [b2], dummyLabels=[(2, 2)])
        frags2 = Chem.GetMolFrags(broken2, asMols=True)
        if len(frags2) != 2:
            return None

        # Pick one fragment from parent 1 and one from parent 2
        f1 = rng.choice(frags1)
        f2 = rng.choice(frags2)

        # Replace dummy atoms with single bond
        combined = Chem.CombineMols(f1, f2)
        rw_mol = Chem.RWMol(combined)

        dummy_indices = [a.GetIdx() for a in rw_mol.GetAtoms() if a.GetAtomicNum() == 0]
        if len(dummy_indices) != 2:
            return None

        # Find neighbors of dummy atoms
        n1 = next(n.GetIdx() for n in rw_mol.GetAtomWithIdx(dummy_indices[0]).GetNeighbors())
        n2 = next(n.GetIdx() for n in rw_mol.GetAtomWithIdx(dummy_indices[1]).GetNeighbors())

        # Add single bond between neighbors
        rw_mol.AddBond(n1, n2, Chem.BondType.SINGLE)

        # Remove dummy atoms (remove higher index first to maintain valid indices)
        for idx in sorted(dummy_indices, reverse=True):
            rw_mol.RemoveAtom(idx)

        return sanitize_molecule(rw_mol.GetMol())
    except (ValueError, RuntimeError, IndexError):
        return None


def apply_random_mutation(mol: Chem.Mol, rng: random.Random | None = None) -> Chem.Mol:
    """Apply one or more random mutation operators, falling back to original if all fail."""
    rng = rng or random
    mutators = [
        mutate_atom_type,
        mutate_add_atom_or_fragment,
        mutate_remove_atom,
        mutate_bond_order,
    ]
    rng.shuffle(mutators)
    for mutator in mutators:
        mutant = mutator(mol, rng=rng)
        if mutant is not None and mutant.GetNumHeavyAtoms() >= 4:
            return mutant
    return mol
