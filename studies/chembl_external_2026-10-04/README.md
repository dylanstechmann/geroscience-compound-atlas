# Prespecified external target-assigned mTOR assay evaluation — 2026-10-04

The frozen mTOR classifiers do **not** transfer reliably enough for compound prioritization on this study/structure holdout. Across all three original scaffold training splits, both classifiers have worse pooled Brier error than a constant training-prevalence predictor. This is a negative surrogate validation result, not evidence about aging or treatment.

## Results

ChEMBL 37 supplied 5,025 IC50 activities and 1,064 assay records through complete, hash-verified API pagination. Fixed eligibility rules retained 2,910 molecule–assay measurements, 2,731 distinct connectivity structures, 194 assays and 170 source documents. All documents, molecule IDs/parent IDs and connectivity structures overlapping the original frozen dataset were excluded. The 2,461 positive / 449 negative rows are imbalanced; average precision must be read against the prevalence baseline of 0.8457.

| Original training seed | Constant baseline Brier | Logistic Brier | Logistic AUROC | HGB Brier | HGB AUROC |
|---|---:|---:|---:|---:|---:|
| 42 | 0.1314 | 0.1614 | 0.5475 | 0.1946 | 0.4775 |
| 123 | 0.1306 | 0.1553 | 0.5544 | 0.1667 | 0.4668 |
| 456 | 0.1305 | 0.1761 | 0.5504 | 0.1748 | 0.4766 |

Lower Brier is better. Primary uncertainty resamples whole source documents for the pooled measured-row Brier difference; the separately labeled equal-document macro analysis gives each document equal weight. These are different estimands. Exact differences and both intervals are retained for every seed in `results.json`. These intervals condition on fixed training models; documents can share laboratories and chemical series. The median maximum training Tanimoto is approximately 0.26. All eligible documents are reported, including one-class documents whose AUROC/AP are undefined (`null`), rather than selecting favorable studies.

## Training-source problem

Under current ChEMBL assay metadata, only 209 of the original 600 raw activity records satisfy this study's stricter direct-human single-protein target-assignment contract. First-exclusion counts: 370 confidence-score exclusions, 14 potential duplicates and 7 assay-description scope exclusions. For example, activity 873240 (assay CHEMBL677281; document CHEMBL1132865) has the description “The inhibitory activity by using FK506 binding protein 12 SPA binding assay” despite its mTOR target annotation. A target ID alone does not identify a direct mTOR activity assay. Even B/confidence9/D target-assignment metadata do not guarantee biochemical measurement: accepted assays also include cellular 4EBP1, AKT and S6 phosphorylation endpoints. The cohort is a mixed-context target-assigned IC50 transfer test. ChEMBL format labels comprise 158 single-protein, 29 cell-based, 6 generic and 1 tissue-based assays (2,621 / 246 / 42 / 1 measured rows). These labels themselves are curation annotations, not proof of the actual format. Format counts and every accepted assay description are recorded in `results.json`. This audit joins the original raw records to the newly retrieved assay metadata; it does not reconstruct the database version of the old download.

The original dataset, thresholds, scaffold splits, model parameters, generator and hypothesis gates remain frozen. Repairing the training dataset is a separate prospective benchmark: it needs assay qualification, source-aware split design and a new untouched evaluation set. Tuning against these external outcomes would consume this holdout. Generated molecule scores remain unvalidated surrogate artifacts.

## Reporting corrections

Independent review corrected the initial “biochemical” scope label without changing the frozen eligibility or outcomes. The original plan remains intact. The primary bootstrap was corrected to its prespecified pooled estimand; equal-document macro uncertainty is secondary. Review also identified a genuine FKBP12-independent mTOR assay falsely rejected by the keyword importer. A regression test and metadata-only correction admitted four additional measurements from one new document/assay, with no threshold, split or parameter change. The original 2,906-row results remain in `pre_correction_results.json`. Other explicit FKBP mentions remain conservatively excluded rather than attempting to infer mechanisms from descriptions. The final result still has worse pooled Brier than prevalence in every model/seed.

## Reproduce

Use the workspace dev container and this repo's Linux virtual environment (or `python -m pip install -e '.[dev]'`). From the repo root:

```sh
python -m bench.external_study fetch
python -m bench.external_study evaluate
pytest -q
```

`fetch` saves actual API responses under ignored `data/raw/chembl_external_2026-10-04/` and fails if the snapshot already has a receipt. Partial, changed-count, repeated or incomplete pagination cannot produce a valid receipt. `evaluate` verifies source bytes and all frozen input hashes before fitting. It uses only original training rows to fit descriptor imputation, scaling and the fixed classifiers. Original validation data are transformed by the existing helper, but never select parameters or thresholds. Day/date, document identity and measured labels do not enter model features.

The committed `plan.json` was written before external predictions. `results.json` contains complete receipts, software versions, exclusion counts, every document result, assay macro errors and fixed-seed document bootstrap intervals. `qualified_rows.json` preserves measurement IDs and within-assay median aggregation; `predictions.csv` records every evaluated row. The source download is about 11 MB and is intentionally outside git. The database can change; exact reproduction requires the hashed local snapshot, not just another live fetch.

## Scope and sources

This is an external **source-document and chemical-structure** holdout within the same curated database. It is not an independent database, prospective experiment, senolysis assay, human validation or rejuvenation study. Enzyme preparation, ATP concentration and assay formats differ. Repeated structures across external assays are retained with assay/document grouping; the original training-source heterogeneity also confounds transfer.

- [ChEMBL API pagination and resource documentation](https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services)
- [ChEMBL assay confidence, activity types and pChEMBL definitions](https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/chembl-data-questions)
- [mTOR target CHEMBL2842](https://www.ebi.ac.uk/chembl/explore/target/CHEMBL2842)
- [ChEMBL data licensing](https://chembl.gitbook.io/chembl-interface-documentation/about)

ChEMBL-derived `qualified_rows.json`, `predictions.csv`, `results.json` and this study's data summaries are attributed to ChEMBL/EMBL-EBI and redistributed under [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/). They are derived by the documented eligibility, median aggregation and prediction steps; the repository's MIT software license does not replace the data license. Raw response URLs and hashes identify all source records; document IDs link to `https://www.ebi.ac.uk/chembl/explore/document/<ID>`.
