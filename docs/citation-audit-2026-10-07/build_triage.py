#!/usr/bin/env python3
"""Identifier triage for data/curated/curated_evidence.csv, audit of 2026-10-07.

Reads the audited table, the refcheck receipt (refcheck_results.json) and the second-source
receipt (second_sources.json), applies the title screen below, and writes triage.csv and
summary.json beside this file. Standard library only. It never edits the curated table and it
makes no network calls.

The screen is a reading of each registry title against the row's compound and recorded claim.
It is not a reading of the paper, and a screen is not a claim-support verdict. Changing SCREEN
needs a written reason in the same commit. tests/test_citation_audit.py fails until the committed
outputs match a rebuild, and build() fails if an identifier has no screen or a receipt does not
cover the table, so a new or edited row cannot slip through unscreened.

Usage: python3 -B docs/citation-audit-2026-10-07/build_triage.py
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
TABLE = REPO / "data" / "curated" / "curated_evidence.csv"
REFCHECK = HERE / "refcheck_results.json"
SECOND = HERE / "second_sources.json"
AUDIT_DATE = "2026-10-07"
RULES_VERSION = 1

ON_TOPIC = "on_topic"
LABEL = "on_topic_label_differs"
READING = "needs_reading"
UNRELATED = "title_unrelated"
NOT_FOUND = "not_found"
CHEMBL_NOT_FOUND = "chembl_not_found"
NO_ID = "no_identifier"
SCREENS = (ON_TOPIC, LABEL, READING, UNRELATED, NOT_FOUND, CHEMBL_NOT_FOUND, NO_ID)
# A row takes the most severe category among its identifiers, in this order.
ROW_PRIORITY = (UNRELATED, NOT_FOUND, CHEMBL_NOT_FOUND, READING, LABEL, ON_TOPIC, NO_ID)

ACTIONS = {
    UNRELATED: (
        "Treat the citation as unsupported for this row until the intended paper is found. Search by "
        "compound and claim, verify any replacement identifier with refcheck, read the paper, and make "
        "the change only with the owner's approval. Do not guess identifiers."),
    NOT_FOUND: (
        "The identifier is absent from PubMed and Europe PMC. Search for the intended paper or mark the "
        "claim unsupported. Keep the row visibly unresolved meanwhile."),
    CHEMBL_NOT_FOUND: (
        "The ChEMBL document ID returns HTTP 404 while a known document returns 200. Confirm in the ChEMBL "
        "web interface, or replace the ID after a search and the owner's approval."),
    READING: (
        "The title does not show the recorded claim. Read the abstract or full text. Keep the row unverified "
        "unless the paper supports the claim."),
    LABEL: (
        "Decide which is wrong, the identifier or the recorded title, and record the decision. Correct the "
        "table only with the owner's approval."),
    ON_TOPIC: (
        "No identifier action. The title fits the claim's subject only; claim support still needs a reading "
        "of the paper before the row is treated as supported."),
    NO_ID: "Nothing to verify. Keep the unverified vendor or gray-market status visible.",
}

NOT_FOUND_NOTE = (
    "No PubMed article record (efetch returned none) and no Europe PMC hit. Europe PMC's positive control "
    "(PMID 19587680) returned one hit.")

# PMID -> (screen, note). Notes quote or paraphrase the registry title as retrieved on AUDIT_DATE; they
# do not describe the paper's contents, which were not read.
SCREEN = {
    # Title names the compound and fits the recorded claim's subject (title level only).
    "19587680": (ON_TOPIC, "Title names rapamycin and lifespan in genetically heterogeneous mice. The recorded title matches."),
    "30279143": (ON_TOPIC, "Title names fisetin as a senotherapeutic that extends health and lifespan. Covers both rows that cite it (senescent-cell burden and SASP); the SASP detail is not in the title."),
    "10839993": (ON_TOPIC, "Title names metformin and inhibition of complex 1 of the mitochondrial respiratory chain, which fits the mild complex I row. The hepatocyte detail is not in the title."),
    "27841876": (ON_TOPIC, "Title names spermidine with cardioprotection and lifespan extension. The macroautophagy mechanism in the row is not in the title."),
    "17086191": (ON_TOPIC, "Title names resveratrol and survival of mice on a high-calorie diet, which fits the high-fat-diet row. The ITP lifespan note in the row is not in the title."),
    "17112576": (ON_TOPIC, "Title names resveratrol with mitochondrial function, SIRT1 and PGC-1alpha, which fits the biogenesis row."),
    "27127236": (ON_TOPIC, "Title names NAD+ repletion with mitochondrial and stem-cell function and lifespan in mice. Fits both nicotinamide riboside rows; the muscle-specific detail is not in the title."),
    "29599478": (ON_TOPIC, "Title names nicotinamide riboside raising NAD(+) in healthy middle-aged and older adults, which fits the human-pilot row. The blood-pressure result is not in the title."),
    "24245565": (ON_TOPIC, "Title names acarbose, 17-alpha-estradiol and NDGA with male-preferential lifespan extension. Fits the sex-difference rows; the recorded title for line 27 agrees."),
    "31195972": (ON_TOPIC, "Title names acarbose with gut-microbiome and fermentation changes and enhanced longevity. The recorded title agrees apart from one article ('the')."),
    "30728281": (ON_TOPIC, "Title names acarbose and the murine gut microbiome, which fits the microbiome row. Line 42 records one title (for 31195972), so this identifier has no recorded title of its own; refcheck's MISMATCH reflects that, not a disagreement."),
    "27400265": (ON_TOPIC, "Title names urolithin A with mitophagy and lifespan in C. elegans and muscle function in rodents. Fits both urolithin A rows."),
    "32990681": (ON_TOPIC, "Title names canagliflozin and male-only lifespan extension in genetically heterogeneous mice. The recorded title matches."),
    "21030672": (ON_TOPIC, "Title names BPC 157 and tendon healing, which fits the tendon row. The row itself says no verified lifespan data exist."),
    # On topic, but the recorded label or title disagrees with this identifier.
    "33788371": (LABEL, "Title is on-topic (17-alpha-estradiol, male lifespan in UM-HET3 mice). The recorded label 'late-start 17alpha-estradiol intervention' is a paraphrase, not this title."),
    "30688027": (LABEL, "Title is on-topic (acarbose, lifespan in aging HET3 mice). The recorded label 'genetically heterogeneous mice' paraphrases it."),
    "25041462": (LABEL, "Title is on-topic (mortality after metformin or sulphonylurea monotherapy in people with type 2 diabetes). The row's recorded title does not describe this paper."),
    "28802803": (LABEL, "Title is on-topic (metformin and all-cause mortality; systematic review and meta-analysis). The row's recorded title is not this paper's title and may describe another paper."),
    # The title does not show the recorded claim.
    "24360282": (READING, "Title is about declining NAD(+) and nuclear-mitochondrial communication during aging. It does not name nicotinamide mononucleotide or the one-week reversal claim."),
    "15545992": (READING, "Title is a general article on type 2 diabetes pathophysiology. It does not name berberine; the row's claim concerns db/db mice."),
    "19149749": (READING, "Title is a general article on chronic inflammation and oxidative stress in age-related disease and cancer. It does not name curcumin."),
    "14501183": (READING, "Title names Epitalon and aging, lifespan and tumor incidence in female Swiss-derived SHR mice. The row's claim (telomerase elongation in a cultured human fibroblast line) is not in the title."),
    # The title is about an unrelated subject.
    "20542383": (UNRELATED, "Title is about smoking status and motivation in heroin-dependent patients. Not rapamycin or autophagy."),
    "26140592": (UNRELATED, "Title is about Argonaute and nucleic-acid guide binding. Not rapamycin or senescence."),
    "19864835": (UNRELATED, "Title is about iron nutrition in infancy. Not rapamycin or hematopoietic stem cells."),
    "28826136": (UNRELATED, "Title is about PCB concentrations in soil along an urban-rural gradient in Shanghai. Not dasatinib, quercetin or senescent cells. Cited by two rows (7 and 11)."),
    "30953051": (UNRELATED, "Title is about left-ventricular remodeling and masked hypertension. Not dasatinib or SASP."),
    "25754376": (UNRELATED, "Title is about synthesizing oxazolidine-2,4-diones (chemistry). Not dasatinib or senescent preadipocytes."),
    "27038166": (UNRELATED, "Title is about antiepileptic drugs and prostate cancer risk. Not quercetin, NF-kB or SASP."),
    "16443825": (UNRELATED, "Title names a mouse palmitoyl-CoA desaturase. Not quercetin or SIRT1."),
    "29704257": (UNRELATED, "Title is about FGFR inhibition in cholangiocarcinoma. Not fisetin or neuronal autophagy."),
    "23708518": (UNRELATED, "Title is about the USP15 deubiquitylase and REST. Not metformin or inflammatory markers."),
    "32130836": (UNRELATED, "Title is about sarcopenia assessment after stroke. Not metformin or epigenetic clocks."),
    "33005886": (UNRELATED, "Title is about RNA-seq library preparation for C. elegans. Not spermidine or eIF5A."),
    "30472097": (UNRELATED, "Title is about VANGL2 stability and integrin alphav. Not spermidine or B-cell defects."),
    "25439502": (UNRELATED, "Title is about intentions to eat low-glycemic-index foods. Not spermidine or EP300."),
    "30139981": (UNRELATED, "Title is about ferromagnetism in a topological semimetal. Not 17-alpha-estradiol or inflammaging."),
    "26711582": (UNRELATED, "Title is about Pacearchaeota and Woesearchaeota in lake surface waters. Not navitoclax or senescent cells. Cited by two rows (30 and 32)."),
    "18698246": (UNRELATED, "Title is about photodynamic therapy for skin lesions in transplant recipients. Not navitoclax or BCL-xL/BCL-2 binding."),
    "20543254": (UNRELATED, "Title is 'Integrated health systems.' Not resveratrol or autophagy."),
    "27797883": (UNRELATED, "Title is about bypass grafting and bilateral brachiocephalic vein occlusion. Not nicotinamide mononucleotide."),
    "31190014": (UNRELATED, "Title is about the isotopic constraint on tropospheric ozone. Not urolithin A or muscle mitochondria."),
    "22086884": (UNRELATED, "Title is about lead, calcium uptake and renal cell carcinoma genetics. Not berberine or AMPK."),
    "25732152": (UNRELATED, "Title is about pain and the bed nucleus of the stria terminalis. Not berberine or hepatocyte autophagy."),
    "28599427": (UNRELATED, "Title is about periostin and PDK1/Akt/mTOR in bladder cancer invasion. Not curcumin or SASP."),
    "16551579": (UNRELATED, "Title names nordihydroguaiaretic acid (NDGA), not curcumin. The row's claim concerns fibrotic stellate-cell gene expression."),
    "26113653": (UNRELATED, "Title is about Aspergillus PCR for clinical use. Not trametinib or lifespan."),
    "29352236": (UNRELATED, "Title is about neural-crest induction and JNK signaling. Not trametinib, MEK or oncogene-induced senescence."),
    # No record in PubMed or Europe PMC.
    "21303889": (NOT_FOUND, NOT_FOUND_NOTE),
    "24095454": (NOT_FOUND, NOT_FOUND_NOTE),
    "26760636": (NOT_FOUND, NOT_FOUND_NOTE),
}

PMID_PART = re.compile(r"PMID:(\d+)")
CHEMBL_PART = re.compile(r"CHEMBL\d+")
TRIAGE_FIELDS = [
    "line", "compound_name", "hallmark", "relation", "grade", "source_review_status", "claim",
    "identifier_kind", "identifier", "refcheck_status", "registry_title", "recorded_source_title",
    "recorded_title_check", "screen", "row_category", "screen_note",
]
LIMITS = [
    "Titles only. The screen reads registry titles as retrieved on the audit date; no abstract or full text was read.",
    "A resolved identifier shows that a record exists and how it is titled. It does not show that the paper supports the row's claim.",
    "Registry records and service responses are dated 2026-10-07 and can change; the receipts here are the record of what was seen.",
    "The curated table was not edited. No replacement identifiers were proposed or applied.",
    "'Not found' means absent from PubMed and Europe PMC on the audit date (ChEMBL: HTTP 404 from its document endpoint, with a positive control passing).",
    "The screen was written by an AI coding assistant. A person with the papers should check each category before any change.",
]


def parse_identifiers(document_ids: str) -> list[tuple[str, str]]:
    """(kind, identifier) pairs in the order written in the table."""
    found = []
    for piece in document_ids.split(";"):
        part = piece.strip()
        if not part:
            continue
        match = PMID_PART.fullmatch(part)
        if match:
            found.append(("PMID", match.group(1)))
        elif CHEMBL_PART.fullmatch(part):
            found.append(("CHEMBL", part))
        else:
            raise ValueError(f"unrecognized identifier in document_ids: {part!r}")
    return found


def normalized_sha256(path: Path) -> str:
    """Hash with line endings normalized, so a CRLF working copy and an LF checkout agree."""
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def build():
    """Return (triage rows as strings, summary dict). Reads only; writes nothing."""
    with TABLE.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        records = list(reader)
    for record in records:
        if len(record) != len(header):
            raise ValueError("a table row does not have one field per column")
        if any("\n" in field or "\r" in field for field in record):
            raise ValueError("a multi-line record would break the line numbers this audit reports")

    registry = {item["ident"]: item for item in json.loads(REFCHECK.read_text(encoding="utf-8"))
                if item["kind"] == "PMID"}
    second = json.loads(SECOND.read_text(encoding="utf-8"))
    chembl = second["chembl_document"]["results"]
    control = second["chembl_document"]["positive_control"]
    if control["http_status"] != 200:
        raise ValueError("the ChEMBL positive control did not return 200, so the 404 results are not interpretable")

    table_pmids = {ident for record in records for kind, ident in
                   parse_identifiers(dict(zip(header, record))["document_ids"]) if kind == "PMID"}
    if set(registry) != table_pmids:
        raise ValueError("refcheck_results.json does not cover exactly the PMIDs in the table; re-run refcheck")
    if set(SCREEN) != table_pmids:
        raise ValueError("SCREEN and the table's PMIDs differ; screen every new identifier and drop stale ones")
    refcheck_not_found = {pmid for pmid, item in registry.items() if item["status"] == "NOT_FOUND"}
    if refcheck_not_found != set(second["not_found_pmids"]):
        raise ValueError("refcheck and second_sources.json disagree about which PMIDs are NOT_FOUND")
    for pmid in refcheck_not_found:
        if second["europe_pmc"]["results"][pmid]["hit_count"] != 0:
            raise ValueError(f"Europe PMC returned hits for {pmid}, which refcheck reported as NOT_FOUND")

    triage: list[dict] = []
    category_by_line: dict[int, str] = {}
    unique: dict[tuple[str, str], str] = {}
    for position, record in enumerate(records):
        line = position + 2  # header is line 1; the table has no multi-line records
        row = dict(zip(header, record))
        base = {
            "line": str(line), "compound_name": row["compound_name"], "hallmark": row["hallmark"],
            "relation": row["relation"], "grade": row["grade"],
            "source_review_status": row["source_review_status"], "claim": row["notes"],
        }
        identifiers = parse_identifiers(row["document_ids"])
        if not identifiers:
            triage.append({**base, "identifier_kind": "", "identifier": "", "refcheck_status": "",
                           "registry_title": "", "recorded_source_title": row["source_title"],
                           "recorded_title_check": "", "screen": NO_ID, "row_category": NO_ID,
                           "screen_note": f"No identifier in the row; source_review_status is {row['source_review_status']}."})
            category_by_line[line] = NO_ID
            continue
        row_rows, screens = [], []
        for kind, ident in identifiers:
            if kind == "PMID":
                screen, note = SCREEN[ident]
                item = registry[ident]
                status, title, check = item["status"], item.get("title") or "", item.get("title_match") or ""
            else:
                result = chembl.get(ident)
                if result is None:
                    raise ValueError(f"line {line}: ChEMBL ID {ident} was not checked")
                if result["http_status"] != 404:
                    raise ValueError(f"line {line}: ChEMBL {ident} returned {result['http_status']}; re-screen it before use")
                screen = CHEMBL_NOT_FOUND
                note = (f"ChEMBL document API returned HTTP {result['http_status']}; the positive control "
                        f"{control['id']} returned {control['http_status']}.")
                status, title, check = f"NOT_FOUND (ChEMBL HTTP {result['http_status']})", "", ""
            screens.append(screen)
            unique[(kind, ident)] = screen
            row_rows.append({**base, "identifier_kind": kind, "identifier": ident, "refcheck_status": status,
                             "registry_title": title, "recorded_source_title": row["source_title"],
                             "recorded_title_check": check, "screen": screen, "screen_note": note})
        category = next(item for item in ROW_PRIORITY if item in screens)
        for item in row_rows:
            item["row_category"] = category
        triage.extend(row_rows)
        category_by_line[line] = category

    rows_by_category = Counter(category_by_line.values())
    identifiers_by_screen = Counter(unique.values())
    decisions = []
    for category in ROW_PRIORITY:
        decisions.append({
            "category": category,
            "row_lines": sorted(line for line, item in category_by_line.items() if item == category),
            "identifiers": sorted(ident for (kind, ident), screen in unique.items() if screen == category),
            "action": ACTIONS[category],
        })
    unique_pmids = [ident for kind, ident in unique if kind == "PMID"]
    summary = {
        "schema_version": 1,
        "audit_date": AUDIT_DATE,
        "screen_rules_version": RULES_VERSION,
        "scope": ("Every row of data/curated/curated_evidence.csv. Identifiers are read from document_ids "
                  "(PubMed IDs and ChEMBL document IDs). The repository's prose documents were scanned with the "
                  "same tool and contain no citation list."),
        "inputs": {
            "table": {"path": TABLE.relative_to(REPO).as_posix(), "rows": len(records),
                      "sha256_lf": normalized_sha256(TABLE)},
            "refcheck": {"path": REFCHECK.relative_to(REPO).as_posix(), "identifiers": len(registry),
                         "sha256_lf": normalized_sha256(REFCHECK)},
            "second_sources": {"path": SECOND.relative_to(REPO).as_posix(),
                               "sha256_lf": normalized_sha256(SECOND)},
        },
        "pmid_occurrences": sum(1 for item in triage if item["identifier_kind"] == "PMID"),
        "unique_pmids": len(unique_pmids),
        "unique_chembl_ids": sum(1 for kind, _ in unique if kind == "CHEMBL"),
        "unique_identifiers": len(unique),
        "rows_by_category": {category: rows_by_category.get(category, 0) for category in ROW_PRIORITY},
        "identifiers_by_screen": {screen: identifiers_by_screen.get(screen, 0) for screen in SCREENS},
        "refcheck_status_counts_pmids": dict(sorted(Counter(registry[ident]["status"] for ident in unique_pmids).items())),
        "decisions": decisions,
        "limits": LIMITS,
    }
    return triage, summary


def main() -> int:
    triage, summary = build()
    with (HERE / "triage.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TRIAGE_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(triage)
    (HERE / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"rows_by_category": summary["rows_by_category"],
                      "identifiers_by_screen": summary["identifiers_by_screen"],
                      "refcheck_status_counts_pmids": summary["refcheck_status_counts_pmids"]},
                     indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
