"""Synthetic tissue-selectivity counterexamples, not rapalog PK or fertility.

Linear tissue partitioning, independent Hill-1 inhibition and an arbitrary
shared exposure scale. No compound-specific parameters or human thresholds.
"""

import json
import math
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent


@dataclass(frozen=True)
class Tissue:
    name: str
    total_partition: float
    free_fraction: float
    ic50: float

    def validate(self):
        if not self.name:
            raise ValueError("Tissue name is required")
        values = (self.total_partition, self.free_fraction, self.ic50)
        if any(not math.isfinite(value) for value in values):
            raise ValueError("All tissue parameters must be finite")
        if self.total_partition < 0 or not 0 <= self.free_fraction <= 1 or self.ic50 <= 0:
            raise ValueError("Partition must be nonnegative, free fraction in [0,1], and IC50 positive")

    def unbound_partition(self):
        self.validate()
        return self.total_partition * self.free_fraction


def inhibition(unbound_index, ic50):
    if not math.isfinite(unbound_index) or unbound_index < 0:
        raise ValueError("Unbound index must be finite and nonnegative")
    if not math.isfinite(ic50) or ic50 <= 0:
        raise ValueError("IC50 must be finite and positive")
    return 0.0 if unbound_index == 0 else 1 / (1 + ic50 / unbound_index)


def odds(fraction):
    if not math.isfinite(fraction) or not 0 < fraction < 1:
        raise ValueError("Illustrative response fractions must be strictly between 0 and 1")
    return fraction / (1 - fraction)


def matched_target(desired_inhibition, target, other):
    """Compare other-tissue inhibition after matching target pathway inhibition.

    rho = (Kp_other*fu_other/IC50_other)/(Kp_target*fu_target/IC50_target).
    A common plasma binding factor cancels. Fractions are synthetic, not an
    efficacy endpoint or acceptable fertility-risk threshold.
    """
    target_odds = odds(desired_inhibition)
    target_factor = target.unbound_partition()
    other_factor = other.unbound_partition()
    if target_factor == 0:
        raise ValueError("Target has no unbound exposure; matching is impossible")
    ratio = (other_factor / target_factor) * (target.ic50 / other.ic50)
    if not math.isfinite(ratio):
        raise ValueError("Combined sensitivity ratio exceeds numeric range")
    return {
        "target_inhibition": desired_inhibition,
        "other_tissue_inhibition": inhibition(ratio * target_odds, 1),
        "relative_other_to_target_normalized_activity": ratio,
        "total_partition_ratio_other_to_target": other.total_partition / target.total_partition,
        "unbound_partition_ratio_other_to_target": other_factor / target_factor,
    }


def largest_ratio_for_bound(target_inhibition, other_inhibition_bound):
    """Pure Hill-1 inequality; these bounds are not biological safety limits."""
    return odds(other_inhibition_bound) / odds(target_inhibition)


def aggregate_responses(target_inhibition, target, weighted_compartments):
    compartments = tuple(weighted_compartments)
    if not compartments:
        raise ValueError("At least one compartment is required")
    weights = [weight for _, weight in compartments]
    if any(not math.isfinite(weight) or weight < 0 for weight in weights):
        raise ValueError("Weights must be finite and nonnegative")
    if not math.isclose(math.fsum(weights), 1, rel_tol=0, abs_tol=1e-12):
        raise ValueError("Illustrative weights must sum to one")
    records = [{"name": tissue.name, "assumed_weight": weight,
                **matched_target(target_inhibition, target, tissue)}
               for tissue, weight in compartments]
    return {
        "compartments": records,
        "weighted_mean_inhibition": math.fsum(row["assumed_weight"]
                                              * row["other_tissue_inhibition"] for row in records),
        "maximum_compartment_inhibition": max(row["other_tissue_inhibition"] for row in records),
    }


def scenarios():
    return [
        ("reference", Tissue("target", 1, 1, 1), Tissue("testis_proxy", 1, 1, 1)),
        ("less_exposure_everywhere", Tissue("target", 0.1, 1, 1), Tissue("testis_proxy", 0.1, 1, 1)),
        ("more_potent_everywhere", Tissue("target", 1, 1, 0.1), Tissue("testis_proxy", 1, 1, 0.1)),
        ("lower_testis_unbound_partition", Tissue("target", 1, 1, 1), Tissue("testis_proxy", 0.1, 1, 1)),
        ("lower_testis_pathway_sensitivity", Tissue("target", 1, 1, 1), Tissue("testis_proxy", 1, 1, 10)),
        ("lower_total_testis_but_same_unbound", Tissue("target", 1, 0.1, 1), Tissue("testis_proxy", 0.1, 1, 1)),
    ]


def main():
    records = [{"scenario": name, "assumed_target": target.__dict__,
                "assumed_other": other.__dict__, **matched_target(0.5, target, other)}
               for name, target, other in scenarios()]
    thresholds = [{"assumed_target_inhibition": target, "assumed_other_bound": other,
                   "largest_relative_other_to_target_activity": largest_ratio_for_bound(target, other)}
                  for target in (0.2, 0.5, 0.8, 0.9) for other in (0.05, 0.1, 0.2)]
    compartments = aggregate_responses(0.5, Tissue("target", 1, 1, 1),
        [(Tissue("basal_proxy", 1, 1, 1), 0.1), (Tissue("adluminal_proxy", 0.01, 1, 1), 0.9)])
    output = {
        "status": "synthetic_counterexamples_no_rapalog_or_fertility_prediction",
        "limits": [
            "No measured compound PK, tissue partition, free fractions or pathway potency.",
            "Target pathway inhibition is not aging benefit or functional efficacy.",
            "All response bounds and compartment weights are arbitrary mathematical examples.",
            "No mTORC2, FKBP12 kinetics, prodrug activation, clearance or reproductive model.",
            "The input scale supplies no human concentration, dose or treatment instruction.",
        ],
        "matched_target_scenarios": records,
        "synthetic_response_bound_grid": thresholds,
        "synthetic_aggregation_counterexample": compartments,
    }
    (ROOT / "synthetic_selectivity.json").write_text(
        json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"matched_target_scenarios": len(records),
                      "synthetic_response_pairs": len(thresholds), "biological_validation": "none"}))


if __name__ == "__main__":
    main()
