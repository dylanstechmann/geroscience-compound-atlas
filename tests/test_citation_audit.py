"""Checks for the 2026-10-07 identifier audit of data/curated/curated_evidence.csv.

The title screen is a judgement, so these tests check its bookkeeping rather than the judgement:
every identifier is screened, the receipts agree with one another, every row appears in exactly one
decision, and the committed outputs are exactly what the builder produces from the current table.
A new or edited row therefore fails here until it has been screened.
"""

import csv
import importlib.util
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "docs" / "citation-audit-2026-10-07"


def load_builder():
    spec = importlib.util.spec_from_file_location("build_triage", AUDIT / "build_triage.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CitationAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.builder = load_builder()
        cls.rows, cls.summary = cls.builder.build()

    def test_every_identifier_is_screened_and_the_counts_add_up(self):
        # build() raises on an unscreened or stale PMID; here the bookkeeping must also add up.
        table_rows = self.summary["inputs"]["table"]["rows"]
        self.assertEqual(sum(self.summary["rows_by_category"].values()), table_rows)
        self.assertEqual(sum(self.summary["identifiers_by_screen"].values()), self.summary["unique_identifiers"])
        chembl_occurrences = sum(1 for row in self.rows if row["identifier_kind"] == "CHEMBL")
        no_identifier_rows = self.summary["rows_by_category"]["no_identifier"]
        self.assertEqual(len(self.rows), self.summary["pmid_occurrences"] + chembl_occurrences + no_identifier_rows)

    def test_every_row_appears_in_exactly_one_decision(self):
        lines = [line for decision in self.summary["decisions"] for line in decision["row_lines"]]
        self.assertEqual(len(lines), len(set(lines)))
        self.assertEqual(sorted(lines), list(range(2, 2 + self.summary["inputs"]["table"]["rows"])))

    def test_the_not_found_identifiers_have_second_sources_and_controls(self):
        second = json.loads((AUDIT / "second_sources.json").read_text(encoding="utf-8"))
        not_found = {row["identifier"] for row in self.rows if row["refcheck_status"] == "NOT_FOUND"}
        self.assertEqual(not_found, set(second["not_found_pmids"]))
        self.assertEqual(len(not_found), 3)
        for pmid in not_found:
            self.assertEqual(second["europe_pmc"]["results"][pmid]["hit_count"], 0, pmid)
        self.assertEqual(second["europe_pmc"]["positive_control"]["hit_count"], 1)
        self.assertEqual(second["chembl_document"]["positive_control"]["http_status"], 200)
        chembl_ids = {row["identifier"] for row in self.rows if row["identifier_kind"] == "CHEMBL"}
        for chembl_id in chembl_ids:
            self.assertEqual(second["chembl_document"]["results"][chembl_id]["http_status"], 404, chembl_id)

    def test_committed_outputs_match_a_rebuild(self):
        committed = json.loads((AUDIT / "summary.json").read_text(encoding="utf-8"))
        self.assertEqual(committed, json.loads(json.dumps(self.summary, ensure_ascii=False)))
        with (AUDIT / "triage.csv").open(encoding="utf-8", newline="") as handle:
            on_disk = list(csv.DictReader(handle))
        self.assertEqual(on_disk, self.rows)

    def test_the_limits_that_stop_this_audit_being_read_as_claim_support_stay_on_record(self):
        limits = " ".join(self.summary["limits"])
        self.assertIn("does not show that the paper supports", limits)
        self.assertIn("The curated table was not edited", limits)
        self.assertIn("claim support still needs", self.builder.ACTIONS[self.builder.ON_TOPIC])


if __name__ == "__main__":
    unittest.main()
