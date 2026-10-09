"""Integration and stale-source failures, not repeated chemical/clinical tests."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from build_catalog import REPO, ROOT, build, existing_file, normalize_records, validate_programs


class CatalogChecks(unittest.TestCase):
    def setUp(self):
        self.manifest = json.loads((ROOT / "input_manifest.json").read_text())
        self.programs = json.loads((ROOT / "programs.json").read_text())

    def test_current_archive_inventory_and_all_graph_links(self):
        catalog = build(self.manifest, self.programs)
        self.assertEqual(len(catalog["records"]), 28)
        self.assertEqual(len(catalog["programs"]), 10)
        self.assertEqual(catalog["role_counts"], {"published_biological_reference": 10, "confirmed_transformation_product": 3, "identity_reference_only": 1, "unvalidated_comparison": 14})
        graph_ids = {r["id"] for r in catalog["records"]}
        self.assertEqual({g for p in catalog["programs"] for g in p["linked_graphs"]}, graph_ids)

    def test_original_evidence_and_qualification_are_preserved(self):
        catalog = build(self.manifest, self.programs)
        for r in catalog["records"]:
            raw = json.loads((REPO / r["source_path"]).read_text())
            source = next(x for x in raw["records"] if x["id"] == r["archive_id"])
            self.assertEqual(r["archive_evidence_status"], source["evidence_status"])
            self.assertEqual(r["archive_qualification"], source[r["qualification_field"]])
            self.assertEqual(r["standard_inchi_key"], source["standard_inchi_key"])

    def test_changed_snapshot_fails_before_any_output(self):
        manifest = copy.deepcopy(self.manifest)
        manifest["datasets"][0]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "Snapshot changed"):
            build(manifest, self.programs)

    def test_new_graph_does_not_inherit_a_role(self):
        dataset = self.manifest["datasets"][0]
        document = json.loads((REPO / dataset["path"]).read_text())
        document["records"][0]["id"] = "new_graph"
        with self.assertRaisesRegex(ValueError, "explicit known role"):
            normalize_records(document, dataset, self.manifest["role_definitions"])

    def test_duplicate_record_or_stale_role_map_fails(self):
        dataset = copy.deepcopy(self.manifest["datasets"][0])
        document = json.loads((REPO / dataset["path"]).read_text())
        duplicate = copy.deepcopy(document)
        duplicate["records"].append(copy.deepcopy(duplicate["records"][0]))
        with self.assertRaisesRegex(ValueError, "safe and unique"):
            normalize_records(duplicate, dataset, self.manifest["role_definitions"])
        dataset["roles"]["removed_graph"] = "unvalidated_comparison"
        with self.assertRaisesRegex(ValueError, "match exactly"):
            normalize_records(document, dataset, self.manifest["role_definitions"])

    def test_program_requires_real_graph_and_both_branches(self):
        catalog = build(self.manifest, self.programs)
        ids = {r["id"] for r in catalog["records"]}
        programs = copy.deepcopy(catalog["programs"])
        programs[0]["linked_graphs"] = ["unknown:graph"]
        with self.assertRaisesRegex(ValueError, "unique known"):
            validate_programs(programs, ids, REPO)
        programs = copy.deepcopy(catalog["programs"])
        programs[0]["if_not_supported"] = ""
        with self.assertRaisesRegex(ValueError, "both decision branches"):
            validate_programs(programs, ids, REPO)

    def test_local_sources_cannot_escape_repository(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "repo"
            root.mkdir()
            outside = Path(directory) / "outside.json"
            outside.write_text("{}")
            with self.assertRaises(ValueError):
                existing_file(root, "../outside.json")
            with self.assertRaises(ValueError):
                existing_file(root, str(outside.resolve()))

    def test_nonfinite_mass_and_missing_qualification_fail(self):
        dataset = self.manifest["datasets"][0]
        original = json.loads((REPO / dataset["path"]).read_text())
        for value in (float("nan"), float("inf"), -1, True):
            doc = copy.deepcopy(original)
            doc["records"][0][dataset["average_mass_field"]] = value
            with self.assertRaisesRegex(ValueError, "molecular mass"):
                normalize_records(doc, dataset, self.manifest["role_definitions"])
        doc = copy.deepcopy(original)
        del doc["records"][0][dataset["qualification_field"]]
        with self.assertRaisesRegex(ValueError, "qualification fields"):
            normalize_records(doc, dataset, self.manifest["role_definitions"])


if __name__ == "__main__":
    unittest.main()
