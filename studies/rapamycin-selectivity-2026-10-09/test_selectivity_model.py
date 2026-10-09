"""Independent mathematical identities and counterexamples; no biological tests."""

import math
import unittest

from selectivity_model import (
    Tissue,
    aggregate_responses,
    inhibition,
    largest_ratio_for_bound,
    matched_target,
    scenarios,
)


class SelectivityChecks(unittest.TestCase):
    def test_hill_one_limits_monotonicity_and_mass_action_identity(self):
        self.assertEqual(inhibition(0, 1), 0)
        self.assertEqual(inhibition(1, 1), 0.5)
        self.assertGreater(inhibition(2, 1), inhibition(1, 1))
        self.assertLess(inhibition(1, 2), inhibition(1, 1))
        for exposure in (0.001, 0.1, 1, 10, 1000):
            response = inhibition(exposure, 1)
            self.assertAlmostEqual(response / (1 - response), exposure, delta=1e-9)

    def test_matched_target_has_independent_expected_off_tissue_response(self):
        target, other = Tissue("target", 1, 1, 1), Tissue("other", 0.1, 1, 1)
        self.assertAlmostEqual(matched_target(0.5, target, other)["other_tissue_inhibition"], 1 / 11)
        self.assertAlmostEqual(matched_target(0.8, target, other)["other_tissue_inhibition"], 2 / 7)

    def test_global_exposure_or_potency_change_does_not_create_selectivity(self):
        for _, target, other in scenarios()[:3]:
            self.assertEqual(matched_target(0.5, target, other)["other_tissue_inhibition"], 0.5)

    def test_total_partition_reduction_can_hide_same_free_exposure(self):
        _, target, other = scenarios()[-1]
        result = matched_target(0.5, target, other)
        self.assertEqual(result["total_partition_ratio_other_to_target"], 0.1)
        self.assertEqual(result["unbound_partition_ratio_other_to_target"], 1)
        self.assertEqual(result["other_tissue_inhibition"], 0.5)

    def test_bound_formula_meets_boundary_and_responds_to_target_demand(self):
        self.assertAlmostEqual(largest_ratio_for_bound(0.5, 0.1), 1 / 9)
        self.assertAlmostEqual(largest_ratio_for_bound(0.8, 0.1), 1 / 36)
        ratio = largest_ratio_for_bound(0.8, 0.1)
        result = matched_target(0.8, Tissue("target", 1, 1, 1), Tissue("other", ratio, 1, 1))
        self.assertAlmostEqual(result["other_tissue_inhibition"], 0.1)

    def test_aggregation_can_hide_critical_compartment_response(self):
        result = aggregate_responses(0.5, Tissue("target", 1, 1, 1),
                                    [(Tissue("basal", 1, 1, 1), 0.1),
                                     (Tissue("adluminal", 0.01, 1, 1), 0.9)])
        self.assertAlmostEqual(result["weighted_mean_inhibition"], 0.05 + 0.9 / 101)
        self.assertLess(result["weighted_mean_inhibition"], 0.1)
        self.assertEqual(result["maximum_compartment_inhibition"], 0.5)

    def test_zero_off_tissue_exposure_and_impossible_target(self):
        target, other = Tissue("target", 1, 1, 1), Tissue("other", 0, 1, 1)
        self.assertEqual(matched_target(0.5, target, other)["other_tissue_inhibition"], 0)
        with self.assertRaises(ValueError):
            matched_target(0.5, Tissue("target", 1, 0, 1), other)

    def test_invalid_input_rejection(self):
        for value in (-1, math.nan, math.inf):
            with self.assertRaises(ValueError):
                inhibition(value, 1)
        for value in (0, -1, math.nan, math.inf):
            with self.assertRaises(ValueError):
                inhibition(1, value)
        for value in (0, 1, -1, math.nan, math.inf):
            with self.assertRaises(ValueError):
                largest_ratio_for_bound(value, 0.1)
            with self.assertRaises(ValueError):
                largest_ratio_for_bound(0.5, value)
        for tissue in (Tissue("", 1, 1, 1), Tissue("bad", -1, 1, 1),
                       Tissue("bad", 1, 1.1, 1), Tissue("bad", 1, 1, 0),
                       Tissue("bad", math.inf, 1, 1)):
            with self.assertRaises(ValueError):
                tissue.unbound_partition()
        for weights in ((), ((Tissue("a", 1, 1, 1), 0.5),),
                        ((Tissue("a", 1, 1, 1), -1),),
                        ((Tissue("a", 1, 1, 1), math.nan),)):
            with self.assertRaises(ValueError):
                aggregate_responses(0.5, Tissue("target", 1, 1, 1), weights)


if __name__ == "__main__":
    unittest.main()
