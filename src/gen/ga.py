"""Population-based Molecular Genetic Algorithm for target-directed generation."""

import logging
import random
from typing import Any

import pandas as pd
from rdkit import Chem

from atlas.normalize import smiles_to_inchikey
from gen.mutations import apply_random_mutation, crossover, sanitize_molecule
from gen.scorer import CompositeSurrogateScorer

logger = logging.getLogger(__name__)


class MolecularGA:
    """Genetic Algorithm evolving molecular structures towards surrogate objective."""

    def __init__(
        self,
        scorer: CompositeSurrogateScorer,
        pop_size: int = 40,
        n_generations: int = 25,
        mutation_rate: float = 0.7,
        crossover_rate: float = 0.3,
        elite_size: int = 5,
        qed_weight: float = 0.2,
        seed: int = 42,
    ):
        self.scorer = scorer
        self.pop_size = pop_size
        self.n_generations = n_generations
        self.mutation_rate = mutation_rate
        self.crossover_rate = crossover_rate
        self.elite_size = elite_size
        self.qed_weight = qed_weight
        self.rng = random.Random(seed)

    def initialize_population(self, seed_smiles: list[str]) -> list[Chem.Mol]:
        """Initialize starting population from seed molecules."""
        valid_seeds = []
        for smi in seed_smiles:
            mol = Chem.MolFromSmiles(smi)
            clean = sanitize_molecule(mol)
            if clean is not None:
                valid_seeds.append(clean)

        if not valid_seeds:
            raise ValueError("No valid seed SMILES provided for GA initialization.")

        # Sample or replicate to reach population size
        population = [self.rng.choice(valid_seeds) for _ in range(self.pop_size)]
        return population

    def tournament_selection(
        self, scored_pop: list[tuple[Chem.Mol, dict[str, Any]]], k: int = 3
    ) -> Chem.Mol:
        """Select parent molecule via tournament selection."""
        candidates = self.rng.sample(scored_pop, min(k, len(scored_pop)))
        best = max(candidates, key=lambda item: item[1]["reward"])
        return best[0]

    def run(
        self, seed_smiles: list[str], max_archive_size: int = 200
    ) -> tuple[pd.DataFrame, dict[str, Any]]:
        """Execute evolutionary loop across generations and return archive of unique molecules."""
        population = self.initialize_population(seed_smiles)

        # Archive maps inchikey -> dict of molecular properties and scores
        archive: dict[str, dict[str, Any]] = {}
        total_proposals = 0
        valid_proposals = 0

        for gen in range(self.n_generations):
            # Score current population
            scored_pop: list[tuple[Chem.Mol, dict[str, Any]]] = []
            for mol in population:
                score_dict = self.scorer.score_molecule(mol, qed_weight=self.qed_weight)
                scored_pop.append((mol, score_dict))

                smi = Chem.MolToSmiles(mol)
                try:
                    ikey = smiles_to_inchikey(smi)
                    if ikey not in archive and score_dict["valid"]:
                        archive[ikey] = {
                            "smiles": smi,
                            "inchikey": ikey,
                            "generation": gen,
                            "reward": score_dict["reward"],
                            "mtor_prob": score_dict["mtor_prob"],
                            "qed": score_dict["qed"],
                            "has_pains": score_dict["has_pains"],
                        }
                except ValueError:
                    pass

            # Sort descending by reward
            scored_pop.sort(key=lambda item: item[1]["reward"], reverse=True)

            # Elitism: retain top individuals
            next_generation: list[Chem.Mol] = [item[0] for item in scored_pop[: self.elite_size]]

            # Breed new offspring
            while len(next_generation) < self.pop_size:
                total_proposals += 1
                parent1 = self.tournament_selection(scored_pop)

                if self.rng.random() < self.crossover_rate:
                    parent2 = self.tournament_selection(scored_pop)
                    child = crossover(parent1, parent2)
                else:
                    child = parent1

                if child is None:
                    child = parent1

                if self.rng.random() < self.mutation_rate:
                    child = apply_random_mutation(child)

                child = sanitize_molecule(child)
                if child is not None:
                    valid_proposals += 1
                    next_generation.append(child)
                else:
                    # If mutation was invalid, re-insert parent1
                    next_generation.append(parent1)

            population = next_generation

        # Assemble archive into DataFrame
        records = list(archive.values())
        archive_df = pd.DataFrame(records)
        if not archive_df.empty:
            archive_df = archive_df.sort_values(by="reward", ascending=False).reset_index(drop=True)
            if len(archive_df) > max_archive_size:
                archive_df = archive_df.iloc[:max_archive_size]

        validity_rate = (valid_proposals / total_proposals) if total_proposals > 0 else 0.0
        stats = {
            "total_proposals": total_proposals,
            "valid_proposals": valid_proposals,
            "validity_rate": validity_rate,
            "unique_molecules_generated": len(archive),
            "archive_size": len(archive_df),
        }
        return archive_df, stats
