"""Gate and anti-clone tests that do not require the frozen mTOR model."""

from rdkit import Chem

from gen.filters import HypothesisGates, evaluate_gates, longest_aliphatic_carbon_chain
from gen.hypothesis import clone_penalty, nearest_tanimoto


def _fp_list(*smiles: str):
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    assert all(m is not None for m in mols)
    from gen.hypothesis import fingerprint

    return [fingerprint(m) for m in mols]


def test_benzene_rejected_as_too_small_and_flat():
    mol = Chem.MolFromSmiles("c1ccccc1")
    ok, reasons = evaluate_gates(mol, HypothesisGates(), has_pains=False)
    assert ok is False
    assert "too_few_heavy_atoms" in reasons or "mol_wt_low" in reasons


def test_greasy_tail_rejected():
    mol = Chem.MolFromSmiles("CCCCCCCCCCCCN")
    assert mol is not None
    assert longest_aliphatic_carbon_chain(mol) > 6
    ok, reasons = evaluate_gates(mol, HypothesisGates(), has_pains=False)
    assert ok is False
    assert "greasy_tail" in reasons or "too_few_rings" in reasons


def test_reasonable_druglike_passes_size_gates():
    # dasatinib-like enough to pass size/element gates; PAINS left off here
    mol = Chem.MolFromSmiles("Cc1nc(Nc2ncc(s2)C(=O)Nc2c(C)cccc2Cl)cc(n1)N1CCN(CCO)CC1")
    assert mol is not None
    ok, reasons = evaluate_gates(mol, HypothesisGates(), has_pains=False)
    assert ok is True, reasons


def test_disallowed_element():
    mol = Chem.MolFromSmiles("c1ccccc1[Si](C)(C)C")
    ok, reasons = evaluate_gates(mol, HypothesisGates(), has_pains=False)
    assert ok is False
    assert "disallowed_element" in reasons


def test_pains_hard_reject():
    rhodanine = Chem.MolFromSmiles("O=C1NC(=S)SC1=Cc1ccccc1")
    ok, reasons = evaluate_gates(rhodanine, HypothesisGates(), has_pains=True)
    assert ok is False
    assert "pains" in reasons


def test_clone_penalty_and_self_tanimoto():
    smiles = "Cc1ccccc1"
    mol = Chem.MolFromSmiles(smiles)
    fps = _fp_list(smiles)
    assert nearest_tanimoto(mol, fps) == 1.0
    assert clone_penalty(0.30, start=0.40, weight=1.25) == 0.0
    assert clone_penalty(1.0, start=0.40, weight=1.25) > 0.5
