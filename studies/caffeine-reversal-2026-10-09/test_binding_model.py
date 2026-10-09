"""Independent equilibrium checks; these do not validate a treatment."""

import math
import unittest

from binding_model import Species, load_models, solve


class BindingModelChecks(unittest.TestCase):
    def test_single_species_matches_analytic_reference(self):
        kd = 4.2e-6
        for rt, gt in ((1e-6, 1e-5), (1e-5, 1e-5), (1e-4, 1e-6)):
            total = rt + gt + kd
            bound = 2 * rt * gt / (total + math.sqrt(total**2 - 4 * rt * gt))
            result = solve(rt, gt, [Species("RG", 1, 1, -math.log10(kd))])
            self.assertAlmostEqual(result["host_captured_guest_M"] / bound, 1, places=12)

    def test_independent_host_dimer_reference(self):
        beta, rt = 10**4.68, 2e-5
        free = 2 * rt / (1 + math.sqrt(1 + 8 * beta * rt))
        result = solve(rt, 0, [Species("R2", 2, 0, 4.68)])
        self.assertAlmostEqual(result["free_host_M"] / free, 1, places=12)
        self.assertIsNone(result["guest_not_host_captured_fraction"])

    def test_no_host_does_not_call_guest_dimers_captured(self):
        result = solve(0, 1e-4, [Species("G2", 0, 2, 0.78)])
        self.assertLess(result["free_guest_monomer_M"], 1e-4)
        self.assertEqual(result["host_captured_guest_M"], 0)
        self.assertEqual(result["guest_not_host_captured_fraction"], 1)
        self.assertAlmostEqual(result["reconstructed_guest_M"] / 1e-4, 1, places=12)

    def test_constructed_state_recovers_free_species(self):
        # Independently form totals from prescribed free concentrations.
        r, g = 8e-6, 3e-6
        values = {
            "R2_dimer": 10**4.68 * r**2,
            "R2_tetramer": 10**12.1 * r**4,
            "caffeine_dimer": 10**0.78 * g**2,
            "R2_C11": 10**4.74 * r * g,
            "R2_C21": 10**10.9 * r**2 * g,
        }
        rt = r + 2 * values["R2_dimer"] + 4 * values["R2_tetramer"] + values["R2_C11"] + 2 * values["R2_C21"]
        gt = g + 2 * values["caffeine_dimer"] + values["R2_C11"] + values["R2_C21"]
        result = solve(rt, gt, load_models()["receptor2_published_NMR"])
        self.assertAlmostEqual(result["free_host_M"] / r, 1, places=12)
        self.assertAlmostEqual(result["free_guest_monomer_M"] / g, 1, places=12)
        for name, expected in values.items():
            self.assertAlmostEqual(result["species_M"][name] / expected, 1, places=12)

    def test_two_guests_per_host_accounting(self):
        r, g, beta = 2e-6, 7e-6, 10**7.1
        rg2 = beta * r * g**2
        result = solve(r + rg2, g + 2 * rg2, [Species("RG2", 1, 2, 7.1)])
        self.assertAlmostEqual(result["free_guest_monomer_M"] / g, 1, places=12)
        self.assertAlmostEqual(result["host_captured_guest_M"] / (2 * rg2), 1, places=12)

    def test_mass_balances_and_nonnegative_grid(self):
        for species in load_models().values():
            for rt, gt in ((0, 0), (1e-6, 0), (0, 1e-6), (1e-9, 1e-4), (1e-3, 1e-8), (1e-4, 1e-4)):
                result = solve(rt, gt, species)
                self.assertLessEqual(abs(result["reconstructed_host_M"] - rt), max(rt, 1e-10) * 1e-12)
                self.assertLessEqual(abs(result["reconstructed_guest_M"] - gt), max(gt, 1e-10) * 1e-12)
                self.assertTrue(all(x >= 0 for x in result["species_M"].values()))
                if gt:
                    self.assertTrue(-1e-14 <= result["guest_not_host_captured_fraction"] <= 1 + 1e-14)

    def test_increasing_host_reduces_uncaptured_guest(self):
        for species in load_models().values():
            fractions = [solve(r * 1e-5, 1e-5, species)["guest_not_host_captured_fraction"] for r in (0, 0.5, 1, 4, 16)]
            self.assertEqual(fractions, sorted(fractions, reverse=True))

    def test_bc500_is_not_the_single_species_kd(self):
        models = load_models()
        full = solve(1e-5, 1e-5, models["receptor2_published_NMR"])
        wrong = solve(1e-5, 1e-5, models["incorrect_BC500_as_Kd"])
        self.assertGreater(abs(full["guest_not_host_captured_fraction"] - wrong["guest_not_host_captured_fraction"]), 0.05)

    def test_reject_invalid_inputs_and_species(self):
        for bad in (-1, float("nan"), float("inf"), 1):
            with self.assertRaises(ValueError):
                solve(bad, 1e-5, [])
        for species in (Species("bad", 0, 0, 1), Species("bad", -1, 2, 1), Species("bad", 1.5, 1, 1), Species("bad", 1, 1, float("nan"))):
            with self.assertRaises(ValueError):
                solve(1e-5, 1e-5, [species])
        with self.assertRaises(ValueError):
            solve(1e-5, 1e-5, [Species("same", 1, 1, 1), Species("same", 2, 1, 2)])


if __name__ == "__main__":
    unittest.main()
