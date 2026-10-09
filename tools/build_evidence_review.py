"""Rebuild the packaged, claim-bound display audit from the existing receipts.

Run from the repository root after installing the project. No network calls,
curated-table edits, replacement identifiers or claim reviews are performed.
"""

import csv
import importlib.util
import json
from pathlib import Path

from atlas.evidence_review import MANIFEST, claim_digest

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "docs" / "citation-audit-2026-10-07"


def build_manifest():
    spec = importlib.util.spec_from_file_location("citation_triage_builder", AUDIT / "build_triage.py")
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    triage, summary = builder.build()
    if summary != json.loads((AUDIT / "summary.json").read_text(encoding="utf-8")):
        raise ValueError("Citation summary is stale; rebuild and review the audit first")
    with (AUDIT / "triage.csv").open(encoding="utf-8", newline="") as handle:
        if triage != list(csv.DictReader(handle)):
            raise ValueError("Citation triage is stale; rebuild and review the audit first")
    with builder.TABLE.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    records = []
    for line, row in enumerate(rows, 2):
        screens = [item for item in triage if int(item["line"]) == line]
        records.append({
            "claim_sha256": claim_digest(row),
            "source_row_line": line,
            "category": screens[0]["row_category"],
            "claim_support_status": "unreviewed",
            "identifiers": [{
                field: item[field] for field in (
                    "identifier_kind", "identifier", "registry_title", "screen", "screen_note",
                )
            } for item in screens],
        })
    return {
        "schema_version": 1, "scope": "identifier_titles_only",
        "audit_date": summary["audit_date"], "inputs": summary["inputs"],
        "limits": summary["limits"], "records": records,
    }


if __name__ == "__main__":
    manifest = build_manifest()
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {len(manifest['records'])} claim-bound title screens to {MANIFEST}")
