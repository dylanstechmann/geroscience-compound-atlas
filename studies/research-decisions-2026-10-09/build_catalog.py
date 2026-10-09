"""Read-only, hash-pinned integration of existing reviews and molecule graphs.

Roles are explicit authored mappings, not inferred evidence grades. This
index copies identity and qualifications; it does not infer drug performance.
"""

import csv
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path
from urllib.parse import urlparse

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent.parent


def existing_file(repo, relative):
    """Validate a repository-local source without allowing traversal/symlink escape."""
    if not isinstance(relative, str) or not relative or Path(relative).is_absolute():
        raise ValueError("Source must be a relative repository path")
    candidate = (repo / relative).resolve()
    if not candidate.is_relative_to(repo.resolve()) or not candidate.is_file():
        raise ValueError(f"Missing or outside-repository source: {relative}")
    return candidate


def normalize_records(document, dataset, allowed_roles):
    records = document.get("records")
    if not isinstance(records, list) or not records:
        raise ValueError("Each archive must contain a nonempty record list")
    rows, seen = [], set()
    for r in records:
        key = r.get("id")
        if not isinstance(key, str) or not re.fullmatch(r"[A-Za-z0-9_]+", key) or key in seen:
            raise ValueError("Archive record IDs must be safe and unique")
        seen.add(key)
        role = dataset["roles"].get(key)
        if role not in allowed_roles:
            raise ValueError(f"An explicit known role is required for {key}")
        mass = r.get(dataset["average_mass_field"])
        if type(mass) not in (int, float) or not math.isfinite(mass) or mass <= 0:
            raise ValueError(f"Invalid average molecular mass for {key}")
        texts = {
            "isomeric_smiles": r.get(dataset["smiles_field"]),
            "neutral_formula": r.get("neutral_formula"),
            "standard_inchi_key": r.get("standard_inchi_key"),
            "archive_evidence_status": r.get("evidence_status"),
            "archive_qualification": r.get(dataset["qualification_field"]),
        }
        if any(not isinstance(value, str) or not value.strip() for value in texts.values()):
            raise ValueError(f"Incomplete identity/qualification fields for {key}")
        if not re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", texts["standard_inchi_key"]):
            raise ValueError(f"Invalid InChIKey representation for {key}")
        rows.append({
            "id": dataset["id"] + ":" + key,
            "archive_id": key,
            "dataset_id": dataset["id"],
            "graph_role": role,
            "neutral_average_mass_Da": mass,
            **texts,
            "qualification_field": dataset["qualification_field"],
            "source_path": dataset["path"],
            "source_sha256": dataset["sha256"],
            "review_path": dataset["review_path"],
        })
    if set(dataset["roles"]) != seen:
        raise ValueError("Role map and archive IDs must match exactly")
    return rows


def validate_programs(programs, graph_ids, repo):
    seen = set()
    needed = ("name", "objective", "starting_point", "first_question", "if_supported", "if_not_supported", "invalid_shortcut")
    for p in programs:
        key = p.get("id")
        if not isinstance(key, str) or not re.fullmatch(r"[a-z][a-z0-9_]+", key) or key in seen:
            raise ValueError("Program IDs must be safe and unique")
        seen.add(key)
        if any(not isinstance(p.get(k), str) or not p[k].strip() for k in needed):
            raise ValueError("Programs need explicit objectives and both decision branches")
        existing_file(repo, p["review_path"])
        url = urlparse(p["primary_url"])
        if url.scheme != "https" or not url.netloc or any(c in p["primary_url"] for c in "\n\r()"):
            raise ValueError("Primary reference must be an HTTPS source URL")
        links = p.get("linked_graphs")
        if not isinstance(links, list) or len(links) != len(set(links)) or not set(links) <= graph_ids:
            raise ValueError("Program graphs must reference unique known catalog entries")
        requirements = p.get("requirements")
        if not isinstance(requirements, list) or not requirements or any(not isinstance(r, str) or not r.strip() for r in requirements):
            raise ValueError("Explicit measurement requirements are required")
    return programs


def build(manifest, program_document, repo=REPO):
    if manifest.get("schema_version") != 1 or program_document.get("schema_version") != 1:
        raise ValueError("Unsupported integration schema")
    datasets = manifest.get("datasets", [])
    if not datasets or len({d["id"] for d in datasets}) != len(datasets):
        raise ValueError("Archive dataset IDs must be nonempty and unique")
    rows, provenance = [], []
    for dataset in datasets:
        source = existing_file(repo, dataset["path"])
        existing_file(repo, dataset["review_path"])
        raw = source.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if digest != dataset["sha256"]:
            raise ValueError(f"Snapshot changed; explicit re-review needed: {dataset['path']}")
        document = json.loads(raw)
        rows.extend(normalize_records(document, dataset, manifest["role_definitions"]))
        provenance.append({"path": dataset["path"], "sha256": digest, "archive_limits": document.get("limits", [])})
    graph_ids = {r["id"] for r in rows}
    if len(graph_ids) != len(rows):
        raise ValueError("Global graph IDs must be unique")
    programs = validate_programs(program_document["programs"], graph_ids, repo)
    if not programs:
        raise ValueError("At least one explicit program is required")
    return {
        "schema_version": 1,
        "status": "review_navigation_and_decision_branches_not_drug_or_evidence_ranking",
        "snapshot_date": "2026-10-09",
        "limits": [
            "No new literature search, human validation or owner source review is implied.",
            "Roles are explicit reviewed mappings, not atlas evidence grades or efficacy rankings.",
            "Archive qualifications retain their original field names and scope.",
            "Descriptors, exact/monoisotopic masses and nominal charges are not pooled across conventions.",
            "No graph or program is selected for human administration.",
        ],
        "role_definitions": manifest["role_definitions"],
        "role_counts": dict(sorted(Counter(r["graph_role"] for r in rows).items())),
        "input_provenance": provenance,
        "records": rows,
        "programs": programs,
    }


def markdown_link(path):
    return "../../" + path


def write_outputs(catalog):
    (ROOT / "research_catalog.json").write_text(json.dumps(catalog, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    with (ROOT / "molecule_index.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(catalog["records"][0]))
        writer.writeheader()
        writer.writerows(catalog["records"])
    lines = ["# Research program decisions", "", "Generated from explicit agent-authored decisions and hash-pinned graph archives.",
             "This is review navigation, not a clinical/evidence ranking or a completed experimental program.", ""]
    for p in catalog["programs"]:
        lines.extend([f"## {p['name']}", "", f"**Objective:** {p['objective']}", "",
                      (f"**Starting evidence:** {p['starting_point']} [Primary reference]({p['primary_url']}); "
                       f"[review and limits]({markdown_link(p['review_path'])})."), "",
                      f"**First question:** {p['first_question']}", "",
                      f"**If resolved in support:** {p['if_supported']}", "",
                      f"**If unresolved or not supported:** {p['if_not_supported']}", "",
                      "**Required measurements:** " + "; ".join(p["requirements"]) + ".", "",
                      f"**Unsupported shortcut:** {p['invalid_shortcut']}", "",
                      "**Defined graph comparisons:** " + (", ".join(f"`{g}`" for g in p["linked_graphs"]) or "No exact graph panel for this program in this snapshot.") , ""])
    (ROOT / "PROGRAMS.md").write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def main():
    manifest = json.loads((ROOT / "input_manifest.json").read_text())
    programs = json.loads((ROOT / "programs.json").read_text())
    catalog = build(manifest, programs)
    write_outputs(catalog)
    print(json.dumps({"graphs": len(catalog["records"]), "programs": len(catalog["programs"]), "roles": catalog["role_counts"]}))


if __name__ == "__main__":
    main()
