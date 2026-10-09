"""Apply the dated identifier screen without promoting it to claim verification.

The packaged manifest matches the entire materialized claim, rather than a CSV
line number or PMID alone. Missing, changed and legacy records remain unreviewed.
"""

from __future__ import annotations

import hashlib
import json
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd

MANIFEST = Path(__file__).with_name("data") / "citation_audit.json"
CLAIM_FIELDS = (
    "compound_name", "compound_id", "inchikey", "target_id", "target_symbol",
    "hallmark", "relation", "provenance", "document_ids", "notes", "source_url",
    "source_title", "species", "sex", "study_design", "endpoint", "comparator",
    "effect_estimate", "uncertainty", "source_locator", "evidence_basis",
    "source_review_status", "intervention_components", "attribution_scope",
)
CATEGORIES = {
    "title_unrelated", "not_found", "chembl_not_found", "needs_reading",
    "on_topic_label_differs", "on_topic", "no_identifier",
}


def text_cell(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def claim_digest(row: Any) -> str:
    values = {field: text_cell(row.get(field)) for field in CLAIM_FIELDS}
    if not values["inchikey"]:
        values["inchikey"] = text_cell(row.get("compound_inchikey"))
    payload = json.dumps(values, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def load_citation_audit(path: str | Path = MANIFEST) -> tuple[dict, str]:
    """Load a versioned title screen; never accept a claim-verification assertion."""
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        if data["schema_version"] != 1 or data["scope"] != "identifier_titles_only":
            raise ValueError("Unsupported citation audit schema or scope")
        if not isinstance(data["audit_date"], str) or not data["audit_date"]:
            raise ValueError("Missing audit date")
        date.fromisoformat(data["audit_date"])
        records = data["records"]
        if not isinstance(records, list) or not records:
            raise ValueError("Missing audit records")
        by_digest = {}
        for record in records:
            digest = record["claim_sha256"]
            if (not isinstance(digest, str) or len(digest) != 64
                    or any(char not in "0123456789abcdef" for char in digest)
                    or digest in by_digest or record["category"] not in CATEGORIES
                    or record["claim_support_status"] != "unreviewed"
                    or not isinstance(record["identifiers"], list) or not record["identifiers"]):
                raise ValueError("Invalid or duplicate audit record")
            for identifier in record["identifiers"]:
                if any(not isinstance(identifier[field], str) for field in (
                    "identifier_kind", "identifier", "registry_title", "screen", "screen_note",
                )):
                    raise ValueError("Invalid identifier screen")
            by_digest[digest] = record
        return {"date": data["audit_date"], "records": by_digest}, "loaded"
    except (OSError, ValueError, TypeError, KeyError):
        return {"date": "", "records": {}}, "unavailable"


def apply_evidence_review(
    edges: pd.DataFrame, *, manifest_path: str | Path = MANIFEST,
) -> pd.DataFrame:
    """Withhold displayed grades while retaining the source table's assertions.

The current manifest verifies no claims. Higher grades cannot be enabled by
editing a parquet status or by resolving a citation to an on-topic title.
"""
    result = edges.copy()
    audit, status = load_citation_audit(manifest_path)
    overlays = []
    for _, row in result.iterrows():
        record = audit["records"].get(claim_digest(row))
        overlays.append({
            "catalog_grade": text_cell(row.get("catalog_grade")) or text_cell(row.get("grade")) or "E0",
            "curator_grade": text_cell(row.get("curator_grade")) or text_cell(row.get("grade")) or "E0",
            "grade": "E0",
            "claim_support_status": "unreviewed",
            "citation_audit_status": "matched" if record else ("unmatched" if status == "loaded" else status),
            "citation_audit_date": audit["date"] if record else "",
            "citation_audit_category": record["category"] if record else "not_audited",
            "citation_audit_identifiers": record["identifiers"] if record else [],
        })
    # Positional assignment also supports frames whose index labels repeat.
    for field in (
        "catalog_grade", "curator_grade", "grade", "claim_support_status",
        "citation_audit_status", "citation_audit_date", "citation_audit_category",
        "citation_audit_identifiers",
    ):
        result[field] = [overlay[field] for overlay in overlays]
    return result
