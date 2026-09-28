"""Hard property gates for hypothesis-mode generation.

These reject common surrogate-spam shapes. They are not a proof of stability,
synthesizability, or biological function.
"""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

DEFAULT_ALLOWED_ATOMIC_NUMS = frozenset({1, 6, 7, 8, 9, 15, 16, 17, 35})


@dataclass(frozen=True)
class HypothesisGates:
    min_heavy_atoms: int = 12
    max_heavy_atoms: int = 38
    min_mol_wt: float = 180.0
    max_mol_wt: float = 520.0
    min_rings: int = 1
    max_rings: int = 5
    max_ring_size: int = 8
    max_rotatable_bonds: int = 10
    min_tpsa: float = 20.0
    max_tpsa: float = 140.0
    min_logp: float = -1.0
    max_logp: float = 5.0
    min_fraction_csp3: float = 0.15
    max_aliphatic_carbon_chain: int = 6
    allowed_atomic_nums: frozenset[int] = DEFAULT_ALLOWED_ATOMIC_NUMS
    reject_pains: bool = True
    require_single_fragment: bool = True

    @classmethod
    def from_mapping(cls, raw: dict) -> "HypothesisGates":
        allowed = raw.get("allowed_atomic_nums")
        return cls(
            min_heavy_atoms=int(raw.get("min_heavy_atoms", 12)),
            max_heavy_atoms=int(raw.get("max_heavy_atoms", 38)),
            min_mol_wt=float(raw.get("min_mol_wt", 180)),
            max_mol_wt=float(raw.get("max_mol_wt", 520)),
            min_rings=int(raw.get("min_rings", 1)),
            max_rings=int(raw.get("max_rings", 5)),
            max_ring_size=int(raw.get("max_ring_size", 8)),
            max_rotatable_bonds=int(raw.get("max_rotatable_bonds", 10)),
            min_tpsa=float(raw.get("min_tpsa", 20)),
            max_tpsa=float(raw.get("max_tpsa", 140)),
            min_logp=float(raw.get("min_logp", -1.0)),
            max_logp=float(raw.get("max_logp", 5.0)),
            min_fraction_csp3=float(raw.get("min_fraction_csp3", 0.15)),
            max_aliphatic_carbon_chain=int(raw.get("max_aliphatic_carbon_chain", 6)),
            allowed_atomic_nums=frozenset(int(z) for z in allowed)
            if allowed is not None
            else DEFAULT_ALLOWED_ATOMIC_NUMS,
            reject_pains=bool(raw.get("reject_pains", True)),
            require_single_fragment=bool(raw.get("require_single_fragment", True)),
        )


def longest_aliphatic_carbon_chain(mol: Chem.Mol) -> int:
    """Longest path of non-aromatic, non-ring carbons. Catches greasy tails."""
    carbon_idx = [
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and not a.GetIsAromatic() and not a.IsInRing()
    ]
    if not carbon_idx:
        return 0
    carbon_set = set(carbon_idx)
    neighbors = {
        i: [n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() in carbon_set]
        for i in carbon_set
    }

    best = 1
    for start in carbon_set:
        stack = [(start, {start})]
        while stack:
            node, seen = stack.pop()
            best = max(best, len(seen))
            for nxt in neighbors[node]:
                if nxt not in seen:
                    stack.append((nxt, seen | {nxt}))
    return best


def largest_ring_size(mol: Chem.Mol) -> int:
    ri = mol.GetRingInfo()
    sizes = ri.AtomRings()
    return max((len(r) for r in sizes), default=0)


def evaluate_gates(
    mol: Chem.Mol | None,
    gates: HypothesisGates,
    has_pains: bool | None = None,
) -> tuple[bool, list[str]]:
    """Return (accepted, reject_reasons)."""
    reasons: list[str] = []
    if mol is None:
        return False, ["unparseable"]

    if gates.require_single_fragment:
        if len(Chem.GetMolFrags(mol, asMols=False)) != 1:
            reasons.append("disconnected")

    heavy = mol.GetNumHeavyAtoms()
    if heavy < gates.min_heavy_atoms:
        reasons.append("too_few_heavy_atoms")
    if heavy > gates.max_heavy_atoms:
        reasons.append("too_many_heavy_atoms")

    if any(a.GetAtomicNum() not in gates.allowed_atomic_nums for a in mol.GetAtoms()):
        reasons.append("disallowed_element")

    mw = float(Descriptors.MolWt(mol))
    if mw < gates.min_mol_wt:
        reasons.append("mol_wt_low")
    if mw > gates.max_mol_wt:
        reasons.append("mol_wt_high")

    rings = int(Lipinski.RingCount(mol))
    if rings < gates.min_rings:
        reasons.append("too_few_rings")
    if rings > gates.max_rings:
        reasons.append("too_many_rings")

    if largest_ring_size(mol) > gates.max_ring_size:
        reasons.append("ring_too_large")

    rotb = int(Lipinski.NumRotatableBonds(mol))
    if rotb > gates.max_rotatable_bonds:
        reasons.append("too_flexible")

    tpsa = float(Descriptors.TPSA(mol))
    if tpsa < gates.min_tpsa:
        reasons.append("tpsa_low")
    if tpsa > gates.max_tpsa:
        reasons.append("tpsa_high")

    logp = float(Descriptors.MolLogP(mol))
    if logp < gates.min_logp:
        reasons.append("logp_low")
    if logp > gates.max_logp:
        reasons.append("logp_high")

    fsp3 = float(Lipinski.FractionCSP3(mol))
    if fsp3 < gates.min_fraction_csp3:
        reasons.append("too_flat")

    if longest_aliphatic_carbon_chain(mol) > gates.max_aliphatic_carbon_chain:
        reasons.append("greasy_tail")

    if gates.reject_pains and has_pains:
        reasons.append("pains")

    return len(reasons) == 0, reasons
