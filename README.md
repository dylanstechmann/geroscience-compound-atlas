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
