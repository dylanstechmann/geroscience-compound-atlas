# Frozen-split chemical-neighborhood reliability diagnostic

Run `make neighborhood PYTHON=.venv/bin/python` or:

```bash
python -m bench.neighborhood_report --out artifacts/neighborhood_reliability.json
```

The report re-fits the configured logistic and histogram-gradient-boosting
models using only each saved Bemis–Murcko scaffold training partition. It reports
held-out error in fixed bands of maximum Morgan radius-2, 2048-bit Tanimoto
similarity to that partition's training molecules. Descriptor median imputation
and standardization use training rows only. No calibrator, activity cutoff,
model settings, or frozen partitions are tuned using the test results.

Input checks reject fractional/duplicate/out-of-range indexes, overlapping or
incomplete partitions, scaffold leakage, nonbinary fingerprints, and activity
labels inconsistent with the configured pChEMBL cutoff. Fingerprints are
reconstructed from each frozen SMILES offline if the sidecar is absent. A
present sidecar must exactly match reconstructed structure fingerprints and
row order. Scaffold annotations are also verified against the structures.

JSON records immutable dataset/split hashes, the ordered fingerprint hash,
configuration/runtime versions and row-level test predictions. Empty bands
remain explicit. Each occupied band includes positive/negative counts, mean
probability, observed active fraction, Brier score, and a 95% Wilson interval
for the observed fraction. Those intervals are descriptive binomial intervals;
they do not account for dependence within scaffolds, uncertainty of a learned
model, or repeated molecules across the three splits. Do not pool splits as
independent external cohorts.

## October 4, 2026 run

The dated report is
[`neighborhood_reliability_2026-10-04.json`](../artifacts/neighborhood_reliability_2026-10-04.json).
It reuses 561 compounds and frozen seeds 42, 123 and 456, each with 56 test
rows. Only **4, 3 and 0** test rows, respectively, have nearest-training
Tanimoto at most 0.50. All three splits therefore provide very little evidence
for chemically remote predictions.

In seed 42's (0.25, 0.50] band, all four compounds are inactive; baseline and
contender mean probabilities are 0.230 and 0.319, with Brier errors 0.061 and
0.178. These four compounds alone cannot define an applicability cutoff or a
calibration correction. A small mean probability/observed-fraction gap can
also hide opposing individual errors; use the row predictions and Brier
scores.

The analysis does **not** validate an applicability domain or extend the
model to optimized novel molecules. A predicted mTOR probability remains a
model output rather than measured binding or a geroscience effect. The atlas
continues to abstain from interpreting surrogate scores as compound evidence.
