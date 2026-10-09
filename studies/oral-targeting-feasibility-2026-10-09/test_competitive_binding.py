"""Independent equilibrium identities and physical limits, not biology tests."""

import math
import unittest

from analysis_models import free_guest_fraction
from competitive_binding import Guest, equilibrate, scenarios


class CompetitiveBindingTests(unittest.TestCase):
    def test_single_guest_matches_analytic_solution(self):
        for host in (0, 0.25, 1, 2):
            for kd in (0.000001, 0.01, 1):
                result = equilibrate(host, [Guest("target", 1, kd)])
                self.assertAlmostEqual(result["guests"][0]["free_fraction"],
                                       free_guest_fraction(host, kd), places=12)

    def test_mass_conservation_capacity_and_mass_action(self):
        for _, host, guests in scenarios():
            result = equilibrate(host, guests)
            bound = sum(row["bound"] for row in result["guests"])
            self.assertLessEqual(bound, host + 1e-12)
            self.assertAlmostEqual(result["free_host"] + bound, host, places=12)
            for row in result["guests"]:
                self.assertAlmostEqual(row["free"] + row["bound"], row["total"], places=12)
                self.assertAlmostEqual(row["free"] * result["free_host"],
                                       row["assumed_kd"] * row["bound"], places=12)

    def test_equal_affinity_guests_have_equal_free_fractions(self):
        rows = equilibrate(1, [Guest("a", 1, 0.01), Guest("b", 3, 0.01)])["guests"]
        self.assertAlmostEqual(rows[0]["free_fraction"], rows[1]["free_fraction"])

    def test_extreme_affinity_preserves_small_free_concentration(self):
        result = equilibrate(1, [Guest("target", 1, 1e-200)])
        # For equal unit totals: h_free^2 + Kd*h_free - Kd = 0,
        # so h_free approaches sqrt(Kd) as Kd tends to zero.
        self.assertTrue(math.isclose(result["free_host"], 1e-100, rel_tol=1e-12))
        self.assertTrue(math.isclose(result["guests"][0]["free"], 1e-100, rel_tol=1e-12))

    def test_competitor_reduces_target_capture(self):
        target = Guest("target", 1, 0.01)
        baseline = equilibrate(1, [target])["guests"][0]["free_fraction"]
        competed = equilibrate(1, [target, Guest("competitor", 20, 0.001)])
        self.assertGreater(competed["guests"][0]["free_fraction"], baseline)

    def test_more_host_improves_capture_under_competition(self):
        guests = [Guest("target", 1, 0.01), Guest("competitor", 20, 0.001)]
        lower = equilibrate(1, guests)["guests"][0]["free_fraction"]
        higher = equilibrate(2, guests)["guests"][0]["free_fraction"]
        self.assertLess(higher, lower)

    def test_empty_and_absent_guests(self):
        self.assertEqual(equilibrate(2, [])["free_host"], 2)
        row = equilibrate(2, [Guest("absent", 0, 0.01)])["guests"][0]
        self.assertEqual(row["free"], 0)
        self.assertEqual(row["bound"], 0)
        self.assertIsNone(row["free_fraction"])

    def test_invalid_inputs(self):
        for host in (-1, math.inf, math.nan):
            with self.assertRaises(ValueError):
                equilibrate(host, [])
        for total, kd in ((-1, 1), (math.inf, 1), (math.nan, 1),
                          (1, 0), (1, -1), (1, math.inf), (1, math.nan)):
            with self.assertRaises(ValueError):
                equilibrate(1, [Guest("target", total, kd)])
        for guests in ([Guest("", 1, 1)], [Guest("a", 1, 1), Guest("a", 2, 2)]):
            with self.assertRaises(ValueError):
                equilibrate(1, guests)
        with self.assertRaises(ValueError):
            equilibrate(1e308, [Guest("target", 1e308, 1)])


if __name__ == "__main__":
    unittest.main()
