"""Unit tests for molecular genetic algorithm, mutations, and surrogate scorer."""

import random

import pandas as pd
from rdkit import Chem

from gen.ga import MolecularGA
from gen.mutations import (
    apply_random_mutation,
    crossover,
    mutate_add_atom_or_fragment,
    mutate_atom_type,
    mutate_remove_atom,
    sanitize_molecule,
)
from gen.scorer import CompositeSurrogateScorer


def test_ga_seed_repeats_independently_of_global_random_state():
    class FixedScorer:
        def score_molecule(self, mol, qed_weight):
            return {
                "valid": True,
                "reward": 1.0,
                "mtor_prob": 0.5,
                "qed": 0.5,
                "has_pains": False,
            }

    def run():
        ga = MolecularGA(
            scorer=FixedScorer(),
            pop_size=8,
            n_generations=3,
            mutation_rate=0.8,
            crossover_rate=0.5,
            elite_size=2,
            seed=42,
        )
        return ga.run(["c1ccccc1CCO", "c1ccncc1CCCO"], max_archive_size=30)

    saved_state = random.getstate()
    try:
        random.seed(1)
        before_first = random.getstate()
        first, first_stats = run()
        assert random.getstate() == before_first

        random.seed(999)
        before_second = random.getstate()
        second, second_stats = run()
        assert random.getstate() == before_second
    finally:
        random.setstate(saved_state)

    assert first_stats == second_stats
    pd.testing.assert_frame_equal(first, second)


def test_sanitize_molecule():
    valid = sanitize_molecule(Chem.MolFromSmiles("CCCO"))
    assert valid is not None
    assert valid.GetNumHeavyAtoms() == 4

    assert sanitize_molecule(None) is None
    # Molecule with < 4 heavy atoms should be rejected
    assert sanitize_molecule(Chem.MolFromSmiles("CC")) is None
    assert sanitize_molecule(Chem.MolFromSmiles("CCO")) is None


def test_mutation_operators():
    mol = Chem.MolFromSmiles("c1ccccc1CCO")
    assert mol is not None

    mutant = mutate_atom_type(mol)
    if mutant is not None:
        assert isinstance(mutant, Chem.Mol)

    frag_mutant = mutate_add_atom_or_fragment(mol)
    if frag_mutant is not None:
        assert frag_mutant.GetNumHeavyAtoms() >= mol.GetNumHeavyAtoms()

    del_mutant = mutate_remove_atom(mol)
    if del_mutant is not None:
        assert del_mutant.GetNumHeavyAtoms() < mol.GetNumHeavyAtoms()

    random_mutant = apply_random_mutation(mol)
    assert isinstance(random_mutant, Chem.Mol)
    assert random_mutant.GetNumHeavyAtoms() >= 4


def test_crossover():
    p1 = Chem.MolFromSmiles("c1ccccc1CCC")
    p2 = Chem.MolFromSmiles("c1ccncc1CCCO")
    assert p1 is not None and p2 is not None

    child = crossover(p1, p2)
    if child is not None:
        assert isinstance(child, Chem.Mol)
        assert child.GetNumHeavyAtoms() >= 4


def test_pains_filter():
    benzene = Chem.MolFromSmiles("c1ccccc1")
    rhodanine = Chem.MolFromSmiles("O=C1NC(=S)SC1=Cc1ccccc1")
    assert benzene is not None and rhodanine is not None

    scorer = CompositeSurrogateScorer()
    assert not scorer.has_pains(benzene)
    assert scorer.has_pains(rhodanine)


def test_surrogate_scorer():
    scorer = CompositeSurrogateScorer()

    # Score valid candidate
    mol = Chem.MolFromSmiles("CN(C)C(=N)NC(=N)N")  # Metformin
    assert mol is not None
    score = scorer.score_molecule(mol, qed_weight=0.2)

    assert score["valid"] is True
    assert 0.0 <= score["mtor_prob"] <= 1.0
    assert 0.0 <= score["qed"] <= 1.0
    assert isinstance(score["has_pains"], bool)
    assert score["reward"] >= 0.0

    # Score invalid candidate
    none_score = scorer.score_molecule(None)
    assert none_score["valid"] is False
    assert none_score["reward"] == 0.0
