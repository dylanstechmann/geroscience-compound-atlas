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

## Quickstart
```bash
make setup
make test
```
