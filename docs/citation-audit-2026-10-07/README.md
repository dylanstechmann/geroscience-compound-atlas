# Citation identifier audit: curated evidence table (2026-10-07)

This is a personal hobby and learning project, developed with substantial assistance from AI coding tools.
This audit was run by an AI coding assistant. Its screens are a reading of registry titles, and a person with
the papers needs to review them before anything in the curated table changes.

## What was asked, and what was not

The audit asks one question of every identifier in `data/curated/curated_evidence.csv`: does it resolve to a
record, and does that record's title fit the compound and the claim recorded on the row?

It does **not** check whether any paper supports its row's claim. No row has been verified as supporting its
claim. The curated table was not edited, and no replacement identifiers were proposed or applied.

## Result

Line numbers below are line numbers in `data/curated/curated_evidence.csv` (the header is line 1).
A row takes the most severe category among its identifiers, in the order listed.

| Row category | Rows | Identifiers | What it means |
|---|---:|---:|---|
| `title_unrelated` | 28 | 26 PMIDs | The registry title is about a different subject, so the citation does not support the row. |
| `not_found` | 3 | 3 PMIDs | No PubMed or Europe PMC record. |
| `chembl_not_found` | 2 | 2 ChEMBL document IDs | ChEMBL's document endpoint returns HTTP 404. |
| `needs_reading` | 4 | 4 PMIDs | The title does not show the recorded claim; the paper has to be read. |
| `on_topic_label_differs` | 3 | 4 PMIDs | The title fits, but the recorded label or title does not match this identifier. |
| `on_topic` | 15 | 14 PMIDs | The title names the compound and fits the claim's subject. Claim support is still unread. |
| `no_identifier` | 2 | none | Vendor or gray-market claims with nothing to verify. |
| **Total** | **57** | **53 identifiers** | 57 PMID occurrences resolve to 51 unique PMIDs. |

Read together, 37 of the 57 rows (`title_unrelated`, `not_found`, `chembl_not_found` and `needs_reading`) cite an
identifier that does not resolve, resolves to an unrelated paper, or resolves to a title that does not show the
claim. The other 18 rows fit at title level only.

## Decisions for the owner

1. **Unrelated titles (lines 3, 4, 6, 7, 8, 9, 11, 12, 13, 17, 21, 22, 24, 25, 26, 28, 30, 31, 32, 35, 36, 45, 47, 48, 50, 51, 52, 53).**
   Treat these citations as unsupported until a replacement is found, verified and read. Several rows share one
   identifier whose title is unrelated: PMID 28826136 (lines 7 and 11) and PMID 26711582 (lines 30 and 32).
2. **Not found (lines 5, 20, 29).** The three PMIDs are absent from PubMed and Europe PMC. Find the intended
   papers or mark the claims unsupported.
3. **ChEMBL document IDs not found (lines 10, 14).** Dasatinib and quercetin rows. Confirm in the ChEMBL web
   interface or replace after a search.
4. **Needs reading (lines 37, 46, 49, 57).** The titles are general or name neither the compound nor the claim.
   Read the papers before deciding.
5. **Label or title disagreements (lines 18, 27, 41).** Decide whether the identifier or the recorded label is
   wrong. Line 18 is one of the six rows with the `source_reviewed_with_claim_limits` status (lines 2, 18, 27,
   41, 42, 54). Its recorded title matches neither of its two identifiers, so that status did not guarantee the
   identifiers. Line 42 records one title for two identifiers; PMID 30728281 has no recorded title of its own.
6. **No identifier (lines 56, 58).** BPC-157 and Epitalon vendor claims. Nothing to verify; keep the status visible.

Any replacement identifier needs the same steps: find the paper by search, verify the new ID with refcheck, read
the paper against the claim, and get the owner's approval before the table changes.

## Method

- **Identifiers.** PubMed IDs and ChEMBL document IDs in the `document_ids` column of all 57 rows.
  The repository's prose documents were scanned with the same tool and contain no citation list.
- **refcheck** (the citation-check skill's `refcheck.py`, run 2026-10-07 with `--records … --all --json`).
  Resolves each PMID against PubMed and compares the registry title with the row's recorded title where one
  exists. Receipt: `refcheck_results.json`. Result: 48 resolved, 3 NOT_FOUND.
- **Second sources for the three NOT_FOUND PMIDs.** PubMed efetch returned no article records. Europe PMC
  returned no hits for any of them. Europe PMC's positive control, PMID 19587680, returned one hit.
  Receipt: `second_sources.json`.
- **ChEMBL document API** for the two ChEMBL IDs, which refcheck does not cover. Both returned HTTP 404. A known
  document, CHEMBL1132865, returned 200. Receipt: `second_sources.json`.
- **Screen.** Each PMID's registry title was read against the row's compound and recorded claim. The reasons are
  written in `SCREEN` in `build_triage.py`. This is a title-level screen: no abstract or full text was read.

## Files

- `triage.csv`: one line per row and identifier, with the registry title, the recorded title, the screen and the
  reason.
- `summary.json`: counts, decisions by row line, input hashes and limits.
- `build_triage.py`: rebuilds both files from the curated table and the two receipts. Standard library only; it
  makes no network calls and never writes to the curated table.
- `refcheck_results.json` and `second_sources.json`: the receipts.

## Reproduce

```bash
python3 -B docs/citation-audit-2026-10-07/build_triage.py   # rebuilds triage.csv and summary.json
python3 -m pytest tests/test_citation_audit.py -q           # coverage, receipts, and outputs against a rebuild
```

Re-running refcheck against the live services is optional. The receipts are what this audit relies on, and
registry records can change. If the curated table changes, `build_triage.py` fails until every new identifier
has a screen, and the test fails until the committed outputs are rebuilt and reviewed.

## Limits

- The screen is about titles. It says whether a title fits the compound and claim, not whether a paper supports
  a sentence.
- Everything is dated 2026-10-07. Services and registry records can change.
- "Not found" means absent from PubMed and Europe PMC on the audit date, and for ChEMBL, a 404 from its document
  endpoint while a known document resolved.
- The screen was written by an AI coding assistant. Review each category against the papers before any change.
