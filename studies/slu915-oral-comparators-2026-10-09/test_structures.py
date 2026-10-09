"""Independent identity, positional-isomer and transformation accounting."""

import json
import unittest

from build_structures import ROOT, SPECIFICATIONS, build_molecule, record
from rdkit import Chem


class StructureChecks(unittest.TestCase):
    def setUp(self):
        self.specs = {s.identifier: s for s in SPECIFICATIONS}

    def test_parent_matches_independently_retrieved_identity(self):
        ref = json.loads((ROOT / "parent_identity_reference.json").read_text())["retrieved_properties"]
        parent = record(self.specs["SLU_PP_915"])
        self.assertEqual(parent["neutral_formula"], ref["MolecularFormula"])
        self.assertEqual(parent["standard_inchi_key"], ref["InChIKey"])
        self.assertAlmostEqual(parent["neutral_exact_mass_Da"], float(ref["ExactMass"]), places=6)
        self.assertAlmostEqual(parent["neutral_average_mass_Da"], float(ref["MolecularWeight"]), delta=0.06)

    def test_source_products_match_neutral_formulas_from_ion_table(self):
        # Primary Table 2 reports [M-H]-; add one H for each neutral graph.
        expected = {"product_M1": "C11H8O3S", "product_M3": "C11H9BO4S", "product_M4": "C17H12FNO2S"}
        for key, formula in expected.items():
            self.assertEqual(record(self.specs[key])["neutral_formula"], formula)

    def test_published_isomers_have_same_formula_different_identity(self):
        records = [record(self.specs[k]) for k in ("SLU_PP_915", "published_10q", "published_10r")]
        self.assertEqual(len({r["neutral_formula"] for r in records}), 1)
        self.assertEqual(len({r["standard_inchi_key"] for r in records}), 3)

    def test_aniline_substitution_positions(self):
        for key, expected_n_distances, expected_f_distance in (
            ("hypothesis_aniline_2_4_difluoro", [3, 5], 4),
            ("hypothesis_aniline_2_5_difluoro", [3, 4], 5),
        ):
            mol = build_molecule(self.specs[key])
            nitrogen = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 7)
            fluorines = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 9]
            # Atom-bond distances include N/F bonds as well as ring bonds.
            distances = sorted(len(Chem.GetShortestPath(mol, nitrogen, f)) - 1 for f in fluorines)
            self.assertEqual(distances, expected_n_distances)
            self.assertEqual(len(Chem.GetShortestPath(mol, *fluorines)) - 1, expected_f_distance)
            self.assertEqual(record(self.specs[key])["neutral_formula"], "C17H12BF2NO3S")

    def test_published_names_match_substitution_distances(self):
        for key, nf_distance, b_thiophene_distance in (
            ("SLU_PP_915", 3, 4), ("published_10q", 4, 4), ("published_10r", 4, 5),
        ):
            mol = build_molecule(self.specs[key])
            idx = {a.GetAtomicNum(): a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in (5, 7, 9, 16)}
            self.assertEqual(len(Chem.GetShortestPath(mol, idx[7], idx[9])) - 1, nf_distance)
            thiophene = next(r for r in mol.GetRingInfo().AtomRings() if idx[16] in r)
            distance = min(len(Chem.GetShortestPath(mol, idx[5], c)) - 1 for c in thiophene)
            self.assertEqual(distance, b_thiophene_distance)

    def test_n_methyl_changes_correct_nitrogen_and_one_donor(self):
        mol = build_molecule(self.specs["hypothesis_amide_N_methyl"])
        self.assertTrue(mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[NX3]([CH3])c")))
        parent = record(self.specs["SLU_PP_915"])
        variant = record(self.specs["hypothesis_amide_N_methyl"])
        self.assertEqual(variant["neutral_formula"], "C18H15BFNO3S")
        self.assertEqual(variant["neutral_graph_descriptors"]["HBD"], parent["neutral_graph_descriptors"]["HBD"] - 1)

    def test_pinacol_ester_masks_two_boron_oh_groups(self):
        mol = build_molecule(self.specs["hypothesis_boronic_pinacol_ester"])
        self.assertFalse(mol.HasSubstructMatch(Chem.MolFromSmarts("B([OX2H])[OX2H]")))
        self.assertTrue(mol.HasSubstructMatch(Chem.MolFromSmarts("B1O[C]([CH3])([CH3])[C]([CH3])([CH3])O1")))
        variant = record(self.specs["hypothesis_boronic_pinacol_ester"])
        self.assertEqual(variant["neutral_formula"], "C23H23BFNO3S")
        self.assertEqual(variant["neutral_graph_descriptors"]["HBD"], 1)

    def test_all_graphs_roundtrip_sdf_and_smiles(self):
        sdf = {m.GetProp("_Name"): m for m in Chem.SDMolSupplier(str(ROOT / "structural_comparators.sdf")) if m is not None}
        self.assertEqual(set(sdf), set(self.specs))
        for spec in SPECIFICATIONS:
            rec = record(spec)
            self.assertEqual(Chem.MolToInchiKey(sdf[spec.identifier]), rec["standard_inchi_key"])
            self.assertEqual(Chem.MolToInchiKey(Chem.MolFromSmiles(rec["canonical_smiles"])), rec["standard_inchi_key"])
            self.assertEqual(Chem.GetFormalCharge(sdf[spec.identifier]), 0)


if __name__ == "__main__":
    unittest.main()
