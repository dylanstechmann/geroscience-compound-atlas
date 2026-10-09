"""Meaningful identity, stereochemistry and site-specific cap accounting."""

import json
import unittest
from dataclasses import replace

from build_structures import ROOT, SPECIFICATIONS, build_molecule, nominal_charge, record
from rdkit import Chem


class StructureChecks(unittest.TestCase):
    def setUp(self):
        self.specs = {s.identifier: s for s in SPECIFICATIONS}

    def test_three_structures_match_official_identity_references(self):
        refs = json.loads((ROOT / "identity_references.json").read_text())["properties"]
        for key, ref in refs.items():
            rec = record(self.specs[key])
            self.assertEqual(rec["neutral_formula"], ref["MolecularFormula"])
            self.assertEqual(rec["standard_inchi_key"], ref["InChIKey"])
            self.assertAlmostEqual(rec["neutral_monoisotopic_mass_Da"], float(ref["MonoisotopicMass"]), places=6)

    def test_all_l_stereochemistry_including_isoleucine_beta(self):
        for s in SPECIFICATIONS:
            expected_count = 16 if s.sequence.startswith("MRW") else 3
            centers = Chem.FindMolChiralCenters(build_molecule(s), includeUnassigned=True)
            self.assertEqual(len(centers), expected_count)
            self.assertTrue(all(label == "S" for _, label in centers))

    def test_terminal_caps_change_expected_formulas(self):
        expected = {
            "MOTS_c_N_acetyl": "C103H154N28O23S2",
            "MOTS_c_C_amide": "C101H153N29O21S2",
            "MOTS_c_both_caps": "C103H155N29O22S2",
            "MOTS_c_K14Q": "C100H148N28O23S2",
            "Epitalon_C_amide": "C14H23N5O8",
            "Epitalon_both_caps": "C16H25N5O9",
        }
        for key, formula in expected.items():
            self.assertEqual(record(self.specs[key])["neutral_formula"], formula)

    def test_amidation_preserves_side_chain_acids(self):
        acid = Chem.MolFromSmarts("C(=O)[OX2H]")
        for key, expected_count in (("Epitalon", 3), ("N_acetyl_Epitalon", 3), ("Epitalon_C_amide", 2), ("Epitalon_both_caps", 2), ("MOTS_c", 2), ("MOTS_c_C_amide", 1)):
            self.assertEqual(len(build_molecule(self.specs[key]).GetSubstructMatches(acid)), expected_count)

    def test_acetylation_is_terminal_and_not_lysine(self):
        mol = build_molecule(self.specs["MOTS_c_N_acetyl"])
        self.assertTrue(mol.HasSubstructMatch(Chem.MolFromSmarts("[CH3]C(=O)N[C@@H](CCSC)C(=O)")))
        lysines = [a for a in mol.GetAtoms() if a.GetPDBResidueInfo() is not None and a.GetPDBResidueInfo().GetResidueNumber() == 14 and a.GetPDBResidueInfo().GetName().strip() == "NZ"]
        self.assertEqual(len(lysines), 1)
        self.assertEqual(lysines[0].GetTotalNumHs(), 2)
        self.assertEqual(lysines[0].GetDegree(), 1)

    def test_charge_matrix_and_k14q_change(self):
        expected = {"MOTS_c": 3, "MOTS_c_N_acetyl": 2, "MOTS_c_C_amide": 4, "MOTS_c_both_caps": 3, "MOTS_c_K14Q": 2, "Epitalon": -2, "N_acetyl_Epitalon": -3, "Epitalon_C_amide": -1, "Epitalon_both_caps": -2}
        for key, charge in expected.items():
            self.assertEqual(nominal_charge(self.specs[key]), charge)
            self.assertEqual(record(self.specs[key])["neutral_graph_formal_charge"], 0)

    def test_all_graphs_roundtrip_sdf_and_smiles(self):
        sdf = {m.GetProp("_Name"): m for m in Chem.SDMolSupplier(str(ROOT / "structural_comparators.sdf")) if m is not None}
        self.assertEqual(set(sdf), set(self.specs))
        for s in SPECIFICATIONS:
            rec = record(s)
            self.assertEqual(Chem.MolToInchiKey(sdf[s.identifier]), rec["standard_inchi_key"])
            self.assertEqual(Chem.MolToInchiKey(Chem.MolFromSmiles(rec["canonical_isomeric_smiles"])), rec["standard_inchi_key"])

    def test_invalid_sequence_fails_explicitly(self):
        for sequence in ("", "aedg", "AEDX"):
            with self.assertRaises(ValueError):
                build_molecule(replace(self.specs["Epitalon"], sequence=sequence))


if __name__ == "__main__":
    unittest.main()
