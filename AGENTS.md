# Agent instructions — geroscience-compound-atlas

Work only in this repository. This is a public evidence atlas plus a frozen
ChEMBL mTOR fingerprint bake-off plus an optional generator. It is not a
formulary.

Keep the public framing as a personal hobby and learning project and preserve
the disclosure of substantial AI coding assistance. Do not add employment
pitches, claims of solo authorship, or claims that the repo proves the owner's
skills unless the user explicitly requests that framing.

## Do not

- Add doses, cycles, vendors, stacks, or “take this.”
- Treat a methylation clock, SenMayo score, or logistic probability as rejuvenation.
- Change `configs/bench.yaml` thresholds or the scaffold split without saying so in `REPORT.md`.
- Train the generator on the test split.
- Paste API keys into git or chat.

## First commands

```bash
python -m pip install -e ".[dev]"
pytest tests/test_hypothesis_filters.py tests/test_splits.py -q
```

## Improve, in this order

1. Run `make hypothesis` if `artifacts/benchmark_dataset.parquet` and `artifacts/splits.json` exist. Read `hypotheses/*/metrics.json` reject counts. If almost everything dies on `too_close_to_training_active`, that is a result, not a bug to quietly loosen.
2. Do not retune `configs/hypothesis.yaml` in the same commit as a code change. If you change a gate, say which spam shape you were killing.
3. Optional next science change, one only: add pChEMBL regression (or a probability calibration plot) on the **same frozen scaffold split** so the 1 µM cliff stops being the only label. Do not add a second target family in the same commit.
4. Leave play mode (`make generate`) able to emit high-score junk. Hypothesis mode is the constrained game. Do not delete the QED=0 sensitivity run.
5. If you touch the dashboard, label generated cards as surrogate artifacts.

## Done when

Tests you ran still pass, `REPORT.md` states what changed, and no sentence tells a person to ingest a molecule.
