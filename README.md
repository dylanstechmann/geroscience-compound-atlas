# Geroscience Compound Atlas + Predictive Bench

## Pitch
> Build a reproducible atlas that maps public compounds onto aging-related mechanisms with graded evidence, then run one honest predictive bake-off (and optionally a constrained generator) so the repo proves I can handle chemical data, splits, baselines, and failure cases — not that I invented a youth pill.

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
