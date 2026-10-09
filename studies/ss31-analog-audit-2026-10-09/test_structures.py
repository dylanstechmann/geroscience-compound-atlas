"""Identity and representation checks; no biological or oral-PK validation."""

import unittest
from dataclasses import replace
from pathlib import Path

from build_structures import SPECIFICATIONS, build_molecule, record
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


class StructuralChecks(unittest.TestCase):
    def setUp(self):
        self.specs = {spec.identifier: spec for spec in SPECIFICATIONS}

    def test_parent_matches_independent_pubchem_identity(self):
        # PubChem PUG REST CID 11764719; retrieval archived in source_receipts.json.
        parent = record(self.specs["SS31"])
        self.assertEqual(parent["neutral_formula"], "C32H49N9O5")
        self.assertEqual(parent["standard_inchi_key"], "SFVLTCAESLKEHH-WKAQUBQDSA-N")
        self.assertAlmostEqual(parent["neutral_average_mass_da"], 639.8, delta=0.02)

    def test_four_correctly_assigned_residue_configurations(self):
        references = {"SS31": ["R", "S", "S", "S"],
                      "SPN4": ["R", "S", "S", "S"],
                      "SS20": ["S", "R", "S", "S"],
                      "SPN10": ["S", "S", "S", "S"],
                      "SS31_terminal_D_Phe": ["R", "S", "S", "R"],
                      "SPN10_D_Arg": ["S", "R", "S", "S"],
                      "SPN10_terminal_D_Lys": ["S", "S", "S", "R"]}
        for identifier, expected in references.items():
            self.assertEqual(record(self.specs[identifier])["alpha_cip_in_sequence_order"], expected)

    def test_stereoisomers_preserve_formula_but_change_identity(self):
        for parent_id, variant_id in (("SS31", "SS31_terminal_D_Phe"),
                                      ("SPN10", "SPN10_D_Arg"),
                                      ("SPN10", "SPN10_terminal_D_Lys")):
            parent, variant = (record(self.specs[key]) for key in (parent_id, variant_id))
            self.assertEqual(parent["neutral_formula"], variant["neutral_formula"])
            self.assertEqual(parent["neutral_average_mass_da"], variant["neutral_average_mass_da"])
            self.assertEqual(parent["neutral_graph_descriptors"], variant["neutral_graph_descriptors"])
            self.assertNotEqual(parent["standard_inchi_key"], variant["standard_inchi_key"])

    def test_acetyl_cap_changes_composition_and_nominal_charge(self):
        parent, capped = (record(self.specs[key]) for key in ("SS31", "SS31_N_acetyl"))
        self.assertEqual(capped["neutral_formula"], "C34H51N9O6")
        self.assertAlmostEqual(capped["neutral_average_mass_da"] - parent["neutral_average_mass_da"],
                               42.037, delta=0.001)
        self.assertEqual(capped["nominal_charge_near_physiological_ph"], 2)
        self.assertEqual(parent["nominal_charge_near_physiological_ph"], 3)

    def test_tryptophan_reference_formula_from_residue_balance(self):
        # 2 Trp + Arg + Lys - 3 H2O; amidation replaces OH with NH2.
        self.assertEqual(record(self.specs["SPN10"])["neutral_formula"], "C34H47N11O4")

    def test_backbone_methylation_preserves_charge_and_removes_one_donor(self):
        parent, changed = (record(self.specs[key]) for key in ("SS31", "SS31_Phe_backbone_N_methyl"))
        self.assertEqual(changed["neutral_formula"], "C33H51N9O5")
        self.assertEqual(changed["nominal_charge_near_physiological_ph"], 3)
        self.assertEqual(changed["neutral_graph_descriptors"]["h_bond_donors"],
                         parent["neutral_graph_descriptors"]["h_bond_donors"] - 1)
        self.assertEqual(changed["alpha_cip_in_sequence_order"], parent["alpha_cip_in_sequence_order"])
        self.assertAlmostEqual(changed["neutral_average_mass_da"] - parent["neutral_average_mass_da"],
                               14.027, delta=0.001)
        molecule = build_molecule(self.specs["SS31_Phe_backbone_N_methyl"])
        tertiary_amide = Chem.MolFromSmarts("[CX3](=O)[NX3]([CH3])[C@@H](Cc1ccccc1)C(=O)N")
        self.assertEqual(len(molecule.GetSubstructMatches(tertiary_amide, useChirality=True)), 1)

    def test_isomeric_smiles_round_trip_and_no_disconnected_salt(self):
        for specification in SPECIFICATIONS:
            row = record(specification)
            molecule = Chem.MolFromSmiles(row["isomeric_smiles"])
            self.assertEqual(Chem.MolToInchiKey(molecule), row["standard_inchi_key"])
            self.assertEqual(len(Chem.GetMolFrags(molecule)), 1)
            self.assertEqual(Chem.GetFormalCharge(molecule), 0)

    def test_sdf_round_trip_preserves_all_identities(self):
        rows = list(Chem.SDMolSupplier(str(Path(__file__).with_name("structural_comparators.sdf"))))
        self.assertEqual(len(rows), len(SPECIFICATIONS))
        for molecule, specification in zip(rows, SPECIFICATIONS, strict=True):
            self.assertIsNotNone(molecule)
            expected = record(specification)
            self.assertEqual(Chem.MolToInchiKey(molecule), expected["standard_inchi_key"])
            self.assertEqual(rdMolDescriptors.CalcMolFormula(molecule), expected["neutral_formula"])
            self.assertEqual(molecule.GetProp("_Name"), specification.identifier)

    def test_invalid_or_unsupported_specifications_fail(self):
        original = self.specs["SS31"]
        for residues in ((), (("Cys", "L"),), (("Arg", "unknown"),)):
            with self.assertRaises(ValueError):
                build_molecule(replace(original, residues=residues))
        for positions in ((1,), (5,), (4, 4), (True,), (2.5,)):
            with self.assertRaises(ValueError):
                build_molecule(replace(original, n_methyl_positions=positions))


if __name__ == "__main__":
    unittest.main()
