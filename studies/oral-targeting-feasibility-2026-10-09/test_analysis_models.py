"""Physical limits and independent reference checks for the illustrative model."""

import math
import unittest

from analysis_models import free_guest_fraction, mass, oral_exposure_proxy


class AnalysisChecks(unittest.TestCase):
    def test_reference_peptide_mass_matches_label_neutral_mass(self):
        # FDA hydrochloride formula C32H49N9O5*3HCl, label mass 749.2 Da.
        neutral = mass({"C": 32, "H": 49, "N": 9, "O": 5})
        self.assertAlmostEqual(neutral + 3 * (1.008 + 35.45), 749.2, delta=0.1)

    def test_host_cannot_bind_more_than_its_capacity(self):
        for host in (0.01, 0.25, 0.5, 1, 2):
            for kd in (0, 1e-12, 1e-6, 0.1, 1, 100):
                free = free_guest_fraction(host, kd)
                self.assertGreaterEqual(free, max(0, 1-host)-1e-12)
                self.assertLessEqual(free, 1)
                bound = 1-free
                if kd > 0:
                    self.assertAlmostEqual(bound * kd, free * (host-bound), delta=1e-10)

    def test_stronger_affinity_does_not_remove_capacity_limit(self):
        self.assertEqual(free_guest_fraction(0.5, 0), 0.5)
        self.assertAlmostEqual(free_guest_fraction(0.5, 1e-12), 0.5, places=10)

    def test_affinity_and_host_capacity_monotonicity(self):
        self.assertLess(free_guest_fraction(2, 0.1), free_guest_fraction(1, 0.1))
        self.assertLess(free_guest_fraction(1, 0.01), free_guest_fraction(1, 0.1))

    def test_high_affinity_equimolar_limit_is_stable(self):
        self.assertGreater(free_guest_fraction(1, 1e-20), 0)
        self.assertEqual(free_guest_fraction(0, 1), 1)

    def test_longer_total_exposure_can_reduce_effective_proxy(self):
        baseline = oral_exposure_proxy(0.01, 1, 1, 1, 1)
        changed = oral_exposure_proxy(0.01, 4, 0.05, 0.5, 1)
        self.assertLess(changed, baseline)
        self.assertEqual(oral_exposure_proxy(0.01, 4, 1, 1, 0), 0)

    def test_invalid_parameters_fail_explicitly(self):
        for value in (-1, math.nan, math.inf):
            with self.assertRaises(ValueError):
                free_guest_fraction(value, 1)
            with self.assertRaises(ValueError):
                oral_exposure_proxy(value, 1, 1, 1, 1)
        with self.assertRaises(ValueError):
            oral_exposure_proxy(1.1, 1, 1, 1, 1)


if __name__ == "__main__":
    unittest.main()
