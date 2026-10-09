# Geroscience Compound Atlas + Predictive Bench

## About this project
This atlas maps public compounds to aging-related mechanisms with graded evidence, compares predictive models on frozen chemical splits, and explores constrained molecule generation.

This is a personal hobby and learning project, developed with substantial assistance from AI coding tools.

## Non-Goals (Hard Constraints)
This repository does **not** implement or provide:
- Personal supplement / peptide / SARM / hormone recommenders
- Dose, cycle, source, or stacking advice
- "Upload labs → protocol" workflows
- A single scalar "anti-aging score" sold as biological age
- Scraping gray-market vendor catalogs as ground truth
- Claims that epigenetic-clock decrease equates to rejuvenation
- Training on private medical records

All outputs are strictly research and educational computational biology artifacts.

## Live Interactive Dashboard
Explore the compound cards, evidence grades, benchmark bake-off, and generated candidate gallery:  
👉 **[https://dylanstechmann.github.io/geroscience-compound-atlas/](https://dylanstechmann.github.io/geroscience-compound-atlas/)**

**Evidence-table status (2026-10-07):** an identifier audit found that 37 of the 57 curated evidence rows cite an identifier that does not resolve, resolves to an unrelated paper, or resolves to a title that does not show the recorded claim. No row is yet verified as supporting its claim. The counts, the row-by-row screen and the owner decision list are in [docs/citation-audit-2026-10-07](docs/citation-audit-2026-10-07/README.md).

**Evidence displays (2026-10-09):** lookups, newly materialized edges and the
dashboard apply that dated audit to each matching claim. All current displayed
grades are E0 pending source-to-claim review; historical catalog and curator
grades remain visible separately. The four historical E4 rows therefore do not
pass the current E4 filter. Each card shows the citation category and identifier
findings. A resolved or on-topic title does not verify claim support. Changed,
legacy or unmatched records remain unreviewed, including when the audit is
unavailable. This also applies when reading an older cached parquet.

The versioned display manifest is packaged with the Python library. Rebuild it
from the existing audit receipts with `python -B tools/build_evidence_review.py`;
the builder refuses stale triage or summary files. No source identifiers,
recorded claims or grades in the curated CSV are changed by that command.

## Oral delivery and tissue-targeting feasibility

The [integrated research decision guide](studies/research-decisions-2026-10-09/README.md)
connects 28 defined molecular comparisons to ten programs, with explicit
evidence roles, preserved source qualifications and decision branches.
Hash-pinned graph inputs make changed or unmapped records fail a rebuild.
The guide organizes the existing reviews; it assigns no efficacy ranking.

The [October 9 research analysis](studies/oral-targeting-feasibility-2026-10-09/README.md)
examines organ-specific aging endpoints, oral SLU-PP-915 and peptide delivery,
SS-31 duration/manufacturing hypotheses, caffeine sequestration and
reproductive-sparing mTOR targeting. It includes primary-source links,
registry/search receipts, reference and unvalidated hypothesis specifications and
synthetic sensitivity calculations. It establishes no improved human-use drug,
oral formulation or aging-reversal effect.

The [SS-31 analog audit](studies/ss31-analog-audit-2026-10-09/README.md) adds
published comparators, nine exact molecular definitions in JSON/SDF, and 2D
depictions. It corrects the earlier Tyr comparator to published SPN4 evidence
and distinguishes a patent's oral-delivery predictions from measured PK.

The [rapamycin selectivity review](studies/rapamycin-selectivity-2026-10-09/README.md)
separates complex selectivity, organ distribution and reproductive outcomes,
with primary-source access receipts and synthetic examples at matched target
pathway response. It establishes no fertility-sparing oral rapalog.

The [caffeine reversal audit](studies/caffeine-reversal-2026-10-09/README.md)
adds measured aqueous binders, DNA/RNA aptamer and catalytic-conversion
evidence, source-access receipts and a reconstruction of published binding
equilibria. It distinguishes aggregate affinity from a one-site Kd and
defines the missing oral exposure, metabolite and sleep evidence.

The [SLU-PP-915 comparator study](studies/slu915-oral-comparators-2026-10-09/README.md)
adds ten exact molecular graphs, published ERR response comparisons and
four unvalidated structural hypotheses. It retains microsomal lower bounds
and separates parent, precursors and confirmed transformation products.

The [peptide-delivery transfer audit](studies/peptide-delivery-transfer-2026-10-09/README.md)
adds MOTS-c CK2/tissue-response evidence, the K14Q persistence counterexample,
Semax/P021 modification precedents and nine exact terminal-comparison graphs.
It also records a native-name MOTS-c trial separately from CB4211 and oral claims.

## Architecture

The project is structured into four core modules:
1. **Atlas (`src/atlas`)**: Ingests compound names, resolves identifiers against PubChem, featurizes them (RDKit descriptors, Morgan fingerprints), and joins them to ChEMBL bioactivities to build a curated evidence graph linking compounds to aging hallmarks.
2. **Bench (`src/bench`)**: Freezes a labeled dataset, generates rigorous Murcko scaffold splits, and trains baseline/contender models (Logistic Regression vs HistGradientBoosting) for a selected benchmark task (e.g., mTOR kinase activity).
3. **Gen (`src/gen`)**: Runs a Genetic Algorithm to generate novel candidate molecules optimizing for the predicted benchmark score, balanced by QED (drug-likeness bias) and filtered by PAINS constraints.
4. **Viz (`src/viz`)**: Generates coverage plots, feature distributions, and the standalone static HTML dashboard summarizing all artifacts.

See the full [SPEC.md](SPEC.md) for the detailed design document and [REPORT.md](REPORT.md) for the experimental write-up and metrics.

## Generation modes

| Command | Mode | What it keeps |
|---|---|---|
| `make generate` | Play | QED weight sweep + PAINS penalty. High-score motif copies can survive. Gallery in `artifacts/generated_molecules.parquet`. |
| `make hypothesis` | Hypothesis | Hard property gates, PAINS reject, Tanimoto cap vs training actives. Dated cards in `hypotheses/YYYY-MM-DD-mtor-hypothesis/`. |

Hypothesis cards are in-silico artifacts with a written way to be wrong. They are not a stack.
The original [2026-09-28 hypothesis archive](hypotheses/2026-09-28-mtor-hypothesis/README.md) used full-benchmark generation seeds. The [corrected training-only archive](hypotheses/training-only/2026-09-28-mtor-hypothesis/README.md) uses only the frozen scaffold training partition. Both archives contain surrogate scores, not measured mTOR activity or validated leads. See [REPORT.md](REPORT.md) for the split correction and rejection counts.

## Quickstart
```bash
# 1. Setup virtual environment and install dependencies
make setup

# 2. Run unit test suite
make test

# 3. Run the full end-to-end data pipeline (Atlas -> Bench -> Gen -> Viz)
make data
make train
make generate
make hypothesis
make report
```

To view the generated dashboard locally:
```bash
# MacOS
open artifacts/dashboard.html

# Windows
start artifacts/dashboard.html

# Linux
xdg-open artifacts/dashboard.html
```

## Single-Compound Evidence & Benchmark Lookup CLI

Look up any watchlist compound (or arbitrary SMILES string) from the command line to inspect identifiers, physicochemical descriptors, curated aging hallmark evidence edges (E0–E4), ChEMBL bioactivities, and predictive benchmark inference:

```bash
# Look up by compound common or ingested name
geroatlas lookup rapamycin
# or via python module
python -m atlas lookup metformin

# Output formatted JSON instead of human-readable report
geroatlas lookup dasatinib --json

# Run inference on an arbitrary candidate SMILES string
geroatlas lookup "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1" --out candidate_card.txt
```

### Dashboard Search & Filter Capabilities
The interactive dashboard (`site/index.html` and `artifacts/dashboard.html`) includes:
- **Multi-Field Real-Time Search**: Matches compound names, synonyms, PubChem CID, SMILES, InChIKey, target symbols, hallmark categories, citations (PMIDs), and mechanistic notes.
- **ChEMBL Assay Filters**: Filter by verified binding assays, general bioactivity records, or historical coverage holes.
- **Quick Preset Buttons**: 1-click filtering for `mTOR Pathway`, `Senolytics`, `AMPK / Metformin`, `Peptides Only`, and `E4 Reported Lifespan / Human Outcome`. The E4 preset currently has no eligible claims pending source-to-claim review.
- **Candidate Gallery Explorer**: Filter and sort the 200 GA-generated candidate molecules by PAINS status (pass/flagged) and sort by composite reward, mTOR probability, or QED drug-likeness.

## Continuous mTOR activity check

For a continuous check on the existing mTOR data, run `make regression` after
installing dependencies. It reads `artifacts/splits.json` without regenerating
the scaffold partitions and writes `artifacts/regression_metrics.json`. Its
pChEMBL predictions describe the frozen ChEMBL table; they are not measured
activity for generated hypothesis cards.

## Changing the benchmark activity threshold

Set `benchmark.thresholds.pchembl_active` in `configs/bench.yaml`, or pass a
custom YAML file with `python -m bench.train --config path/to/bench.yaml`.
The cutoff is inclusive: measurements equal to it are labeled active.

Training recalculates labels from cached `pchembl_value` measurements before
splitting or fitting. On cached runs, it refreshes the benchmark parquet and
`data/processed/benchmark_dataset.csv`; it reuses the fingerprint sidecar and
does not contact ChEMBL when both cached files exist. Metrics record the numeric
`pchembl_threshold` and the matching `label_definition`. Regenerate downstream
reports and generator outputs after changing the threshold.

New datasets preserve the full precision of averaged activity measurements.
Older caches can only be relabeled at their stored precision; rebuild them from
the raw activity records if measurements near the cutoff were rounded. Invalid
thresholds or non-finite cached measurements fail explicitly instead of silently
changing the labels.

## Frozen chemical-neighborhood reliability

`make neighborhood PYTHON=.venv/bin/python` audits both classifiers on the
existing scaffold splits, including error by nearest training-fingerprint
similarity, class counts, uncertainty summaries, row predictions, and exact
input hashes. It verifies split isolation and fingerprint/structure alignment,
and reconstructs missing fingerprints offline without rewriting the benchmark.
See [the diagnostic and its limits](docs/neighborhood-reliability.md) and
[the October 4 report](artifacts/neighborhood_reliability_2026-10-04.json).
The existing test folds contain very few chemically remote compounds, so this
does not establish reliability for novel optimized structures.


## Public-data external validation (2026-10-04)

The fixed classifiers were tested against 2,910 qualified molecule–assay rows from 170 source documents and 2,731 connectivity structures excluded from the original dataset. Both models had worse pooled Brier error than the training-prevalence baseline in all three seeds; HGB AUROC was 0.467–0.478. Current assay metadata also shows that only 209/600 original raw records qualify as direct human single-protein target-assignment records under the new strict contract. The frozen benchmark and generator are unchanged; their scores are not qualified for compound prioritization. See the [prespecified study, all-source results and reproduction steps](studies/chembl_external_2026-10-04/README.md).

The interactive dashboard repeats this qualification result beside its internal
split metrics and generated-molecule gallery so surrogate probabilities are
not mistaken for externally validated activity or evidence of rejuvenation.
