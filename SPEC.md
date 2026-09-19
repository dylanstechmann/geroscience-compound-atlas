# Geroscience Compound Atlas + Predictive Bench
## Master specification for humans and AI agents

**Owner:** Dylan Stechmann  
**Stack context:** MS CS / AI (FAU), current course = Reinforcement Learning, spring = Deep Learning. Solo build. Portfolio + job applications ~Feb 2027.  
**Domain:** computational geroscience — public compounds, public assays, named aging mechanisms.  
**Not:** a protocol, a stack recommender, a “biological age” consumer app, or medical advice.

Paste this entire document at the start of implementation chats. If an agent proposes features that contradict **Non-goals** or **Safety**, refuse and point back here.

---

## 1. One-sentence pitch

Build a reproducible atlas that maps public compounds onto aging-related mechanisms with graded evidence, then run one honest predictive bake-off (and optionally a constrained generator) so the repo proves I can handle chemical data, splits, baselines, and failure cases — not that I invented a youth pill.

---

## 2. Problem this exists to solve

Aging biology is multi-causal (genomic instability, epigenetic drift, senescence/SASP, nutrient sensing, inflammation, ECM/fibrosis, stem-cell decline, proteostasis). Popular discussion collapses this into “fix the methylome” or “take the vial.”

Drug-like work needs:

1. A **mechanism tag** that is more precise than “anti-aging.”
2. An **evidence grade** (mouse lifespan ≠ clock shift ≠ blog).
3. A **predictive task** with a split that does not leak chemical cousins.
4. Explicit humility: QED and “approved drug” are biased priors, not truth.

This project produces (1)–(3) as artifacts. It does not tell a person what to ingest.

---

## 3. Non-goals (hard)

Do **not** implement:

- Personal supplement / peptide / SARM / hormone recommenders
- Dose, cycle, source, or stacking advice
- “Upload labs → protocol”
- A single scalar “anti-aging score” sold as biological age
- Scraping gray-market vendor catalogs as ground truth
- Claims that epigenetic-clock decrease = rejuvenation
- Training on private medical records

If a user later wants a watchlist of *names they already collected* (PubChem resolve only), that is **identifier resolution**, not recommendation.

---

## 4. Safety and scientific posture

- All outputs are **research / educational**.
- Prefer peer-reviewed assays and registered endpoints over marketing copy.
- Separate **sequence mutations** from **epigenetic marks** from **cell-state (senescence)** from **signaling (mTOR, IL-11, etc.)**. Do not merge them into one “age knob.”
- When using drug-likeness (QED, Lipinski), label it: *prior fitted on historical approved sets; systematically unkind to peptides, natural products, and new modalities.*
- No generation of scheduled-substance synthesis routes. If a generated structure is a controlled scaffold, drop it from public galleries or hash-only.

---

## 5. Success criteria (what “heavy duty and done” means)

A stranger can, in one evening:

1. Clone the repo, run `make data` + `make train` + `make report` (or equivalent).
2. Open a dashboard or static HTML/notebook that shows:
   - compound cards with SMILES, InChIKey, mechanism tags, evidence grade
   - coverage holes (names with no ChEMBL activity)
   - one model comparison table
3. Read `REPORT.md` that states the question, split, metric, what failed, and what we are *not* claiming.

Resume test: three bullets, each with a number (e.g. “N compounds resolved, scaffold-split AUROC X vs Morgan+HGB baseline Y, Z% ChEMBL coverage on watchlist”).

---

## 6. Architecture (three modules, one repo)

```
gerosci-atlas/
  README.md
  REPORT.md                 # living lab notebook; final narrative
  SPEC.md                   # this file, or a short pointer to it
  pyproject.toml
  Makefile / justfile
  configs/                  # all knobs; no magic numbers in code
  data/
    raw/                    # gitignored dumps
    interim/
    processed/
    watchlist/              # optional user name list → CIDs (not advice)
  src/
    atlas/                  # ingest, normalize, evidence graph
    bench/                  # splits, models, metrics
    gen/                    # optional RL / constrained generation
    viz/
  notebooks/                # exploration only; logic lives in src/
  tests/
  figures/
  artifacts/                # frozen tables for the report
```

### Module A — Atlas (knowledge graph / tables)

**Entities**

- `Compound`: name, synonyms, CID, InChIKey, canonical SMILES, InChI, modality (`small_molecule` | `peptide` | `antibody` | `other`)
- `Target`: UniProt / ChEMBL target id, gene symbol, protein family
- `Hallmark`: controlled vocabulary (see §7)
- `Assay`: ChEMBL assay id, type (binding, functional, ADMET, phenotypic), organism, standard_type (IC50, Ki, EC50, …)
- `Document`: PMID / DOI, year, species, endpoint class
- `EvidenceEdge`: compound —[relation]→ target or hallmark, with grade and provenance

**Relations (examples)**

- `binds`, `inhibits`, `agonizes`, `degrades`, `phenotypic`
- `claimed_for_hallmark` (weak; from reviews)
- `tested_in` (species + endpoint)

**Evidence grade (required on every edge)**

| Grade | Meaning |
|------|---------|
| E0 | Name only / vendor copy / no structured assay |
| E1 | In vitro assay, single paper or ChEMBL row |
| E2 | Multiple consistent assays or clear target engagement |
| E3 | In vivo functional endpoint (not just a clock) in a model organism |
| E4 | Mammalian lifespan / healthspan or registered human outcome trial |

A senolytic claim without a senescent-vs-proliferating selectivity number is at most E2.

### Module B — Bench (the job-market core)

Pick **exactly one primary task** for v1. Do not train five models on five fantasies.

**Recommended primary task (default): senolytic-like selectivity**

- Positive: compounds reported to preferentially kill or disable senescent cells (literature + any public screens you can legally use).
- Negative: cytotoxic-to-everything chemotypes and close decoys.
- Features v1: RDKit descriptors + Morgan fingerprints.
- Model v1: HistGradientBoosting or logistic regression (baseline first).
- Split: **Murcko scaffold split** (or Bemis-Murcko). Report random-split numbers only as a *leakage diagnostic*, never as the headline.
- Metrics: AUROC, AUPRC, recall at 5% FPR, calibration. Show 10 false positives.

**Alternative primary task: single target-family activity**

If senolytic labels are too sparse or too dirty, switch to a dense ChEMBL family adjacent to aging signaling, for example:

- PI3K/AKT/mTOR
- BCL-2 family
- a kinase panel with enough rows

Same split discipline. Classification at a declared nM threshold *or* pIC50 regression. Headline metric on scaffold split.

**Spring upgrade (Deep Learning class):** replace or add a graph encoder (Chemprop, GIN, or a small MPNN). Same splits frozen. The point is the encoder change, not a new question.

### Module C — Generator (optional, RL-class hook)

Only after A and B produce a frozen scorer.

- Representation: SMILES (REINVENT-style) **or** a simple graph GA. Start with a genetic algorithm over SMILES/SELFIES if RL is unstable.
- Reward (example): `valid * unique * (QED_weight * qed + sel_weight * bench_score - pains_penalty)`.
- Expose QED weight as a **bias control**, not a moral good.
- Compare generated cloud vs watchlist in UMAP/PCA of fingerprints. Question: do we wander into the same property island as known geroprotective chemotypes, or into junk space?
- Cap compute. Publish 200 structures max in the gallery, with validity/uniqueness/internal-diversity.

If RL training is a mess, ship the GA and write why PPO/REINFORCE failed. That is a valid REPORT section.

---

## 7. Controlled hallmark vocabulary

Use short stable slugs. Map papers onto these; do not invent twenty near-synonyms.

- `genomic_instability`
- `epigenetic_alteration`
- `telomere_attrition`
- `proteostasis`
- `autophagy`
- `nutrient_sensing`
- `mitochondrial_dysfunction`
- `cellular_senescence`
- `stem_cell_exhaustion`
- `intercellular_communication`
- `chronic_inflammation`
- `ecm_fibrosis`
- `dysbiosis` (optional; easy to fake — skip if no data)

A compound may have multiple tags. Each tag needs its own evidence grade.

**Working scientific priors for this repo (do not “fix” these in code comments without a citation):**

- Epigenetic clocks are correlates; 2025 work ties much clock signal to mutation-adjacent methylome remodeling. Resetting marks ≠ spell-checking the genome.
- Senescence + SASP have stronger *functional* intervention evidence in mice than “mutation load is the whole aging program.”
- Nutrient sensing (mTOR) still has the cleanest small-molecule lifespan literature in mammals.
- IL-11 / related inflammatory nodes are a newer, high-interest axis (mouse lifespan with late-life blockade). Treat as E3-class biology, not as a supplement aisle.

---

## 8. Data sources (v1)

Use only what you can download and cite.

| Source | Use |
|--------|-----|
| PubChem | Name → CID, SMILES, InChIKey, synonyms |
| ChEMBL (SQL dump or REST, then cache) | assays, targets, activities, documents |
| RDKit | descriptors, fingerprints, Murcko scaffolds, QED, SA if available |
| PubMed / open reviews | mechanism tags + endpoint class; store PMID |
| Optional: CellAge, Senescence-associated gene lists, DrugAge / Geroprotectors DBs | seed labels — treat as noisy |
| Optional: ZINC250k / GuacaMol subset | unlabeled pool for generation |

License: keep a `DATA_LICENSES.md`. Do not republish entire ChEMBL dumps in GitHub LFS if ToS frowns; publish scripts + frozen *derived* tables of modest size.

**Watchlist file** (optional, private or public): the user’s long chemical name list is **only** for identifier resolution and coverage plots. Schema: `raw_name, resolved_cid, smiles, modality, notes`. No “take this.”

---

## 9. Implementation phases

Agents: do not skip ahead. Each phase ends with a commit + a paragraph in REPORT.md.

### Phase 0 — Repo skeleton (1–2 evenings)

- pyproject with `rdkit-pypi` or conda pin, pandas, numpy, scikit-learn, pyyaml, httpx
- Makefile targets: `setup`, `test`, `lint`
- Empty schemas as pydantic models or typed dicts
- README with the pitch and non-goals

### Phase 1 — Resolve and featurize (week 1–2)

- Watchlist + a seed of well-known gerontology names (rapamycin, dasatinib, quercetin, fisetin, metformin, spermidine, 17α-estradiol, navitoclax, etc.) resolved to InChIKey
- RDKit descriptor table
- Coverage plot: resolved vs failed names; small molecule vs peptide
- **Exit:** `artifacts/compounds.parquet` + figure

### Phase 2 — ChEMBL join (week 2–4)

- Activities joined on InChIKey / ChEMBL molregno
- Target + assay tables
- Hallmark tags for a *hand-curated* subset first (50–150 edges). Quality over auto-NLP at this stage.
- Evidence grades assigned by hand for the curated set; document rules in `configs/grading.yaml`
- **Exit:** graph or relational tables + “% of watchlist with ≥1 binding assay”

### Phase 3 — Bench v1 (week 4–8)  ★ portfolio critical path

- Freeze label definition in config
- Scaffold split, three seeds
- Baseline model + one slightly stronger model
- Plots: ROC, PR, scaffold-size histogram, top mistakes
- **Exit:** `artifacts/metrics.json` + notebook-free CLI

### Phase 4 — Report and viz (week 8–10)

- Static site or Streamlit *read-only* on processed files
- Compound card page
- REPORT.md narrative (see template §12)

### Phase 5 — Optional generator (after Phase 3 is frozen)

- GA first, RL second
- Diversity vs predicted score Pareto
- Side-by-side UMAP: known E3/E4 compounds vs generated

### Phase 6 — Spring DL swap

- Frozen split files from Phase 3
- Graph model vs fingerprint baseline
- One ablation (no atom features / no edges / random features)

---

## 10. Engineering standards

- Pin versions. Record seeds.
- No hidden downloads in `train.py` without a cache dir.
- Unit tests: SMILES canonicalize, scaffold assignment, split disjointness, grade enum
- CI optional; at least `pytest tests/test_splits.py`
- Large raw files gitignored
- Config-driven thresholds (`pchembl >= 6.0`, etc.)

---

## 11. Suggested stack

- Python 3.11+
- RDKit, pandas, pyarrow
- scikit-learn (v1 models)
- Later: PyTorch + PyG or Chemprop
- Optional viz: Streamlit or a generated Plotly HTML
- Optional RL: REINVENT leftovers are heavy; a 300-line REINFORCE-on-SELFIES is enough to *learn*, not to beat Nature

Do not start with a microservices fantasy.

---

## 12. REPORT.md template (fill this; do not delete sections)

```markdown
# Report — Geroscience Compound Atlas + Bench

## Question
What we asked, in one paragraph.

## What this is not

## Data
Sources, dates downloaded, N compounds, N assays, N labels.

## Labels and leakage
Definition. Scaffold split. Random-split AUROC vs scaffold AUROC (leakage gap).

## Models
Baseline. Contender. Compute budget.

## Results
Table. Plots. Ten worst errors and a guess why.

## Hallmark coverage
Which mechanisms have chemistry vs which are empty slogans in this atlas.

## Generator (if any)
Validity, uniqueness, novelty vs training InChIKeys. QED weight sensitivity.

## Claims we refuse
No human dosing. No “reversed aging.” No clock-as-outcome.

## Next
Spring graph encoder / more hand-curated edges / drop generator.
```

---

## 13. Resume bullets (fill numbers when real)

Template:

- Built a geroscience compound atlas resolving N names to PubChem/ChEMBL and grading mechanism evidence from in-vitro assays to in-vivo endpoints.
- Trained fingerprint baselines under a Murcko scaffold split for [task]; reported AUROC/AUPRC vs random-split leakage and analyzed false positives.
- (Optional) Constrained generator (GA or RL) optimized for predicted [task] score with explicit QED-bias control; compared generated chemical space to curated geroprotective compounds.

---

## 14. Prompts for later AI sessions

**Session type: implement Phase K**

> Read SPEC.md. Implement only Phase K. Do not add features from later phases. Show the exit artifact. If data access fails, implement a mocked fixture of 20 compounds so tests pass, and document the blocker.

**Session type: review**

> Read SPEC.md §3 and §5. Review the current repo against non-goals and success criteria. List drift. Do not rewrite style.

**Session type: science check**

> Challenge any place we equated epigenetic clock movement with functional rejuvenation, or treated QED as goodness.

**Session type: job-market cut**

> Cut the repo to what a hiring manager will run in 20 minutes. Propose deletions.

---

## 15. Default decisions (so agents stop asking)

| Fork | Default |
|------|---------|
| Primary bench task | Senolytic-like selectivity if ≥80 usable labeled mols after cleanup; else ChEMBL mTOR/PI3K pChEMBL regression/class |
| Split | Murcko scaffold, 80/10/10 or 5-fold by scaffold groups |
| v1 model | Morgan 2048-r2 + HistGradientBoosting |
| Graph model | Deferred to spring |
| Generator | Off until metrics.json exists |
| UI | Static HTML + parquet; Streamlit only if time |
| Watchlist | Optional; never the training label source by itself |
| Peptides on watchlist | Resolve if possible; exclude from fingerprint models or featurize separately — do not pretend Morgan understands BPC-157 |

---

## 16. Glossary (keep this in the spec so later chats stay readable)

- **SMILES / SELFIES:** text encodings of molecules. SELFIES is more valid-by-construction.
- **InChIKey:** hashed standard identifier; join key of choice.
- **ChEMBL / PubChem:** public chemical + assay databases.
- **pChEMBL:** −log10 of standard activity in molar; ~6 means 1 µM.
- **QED:** drug-likeness heuristic from historical oral drugs.
- **PAINS:** motifs that fool assays.
- **Murcko scaffold:** ring-core used to stop train/test cousins.
- **SASP:** inflammatory secretions of senescent cells.
- **Somatic mutation vs epigenetic mark:** sequence change vs regulation change. Coupled, not identical.
- **OSK / OSKM:** Yamanaka reprogramming factors. Out of scope for v1 implementation.

---

## 17. First command block for an implementation agent

```text
Create the repo skeleton exactly as §6. Add pydantic schemas for Compound, EvidenceEdge, and EvidenceGrade. Add tests that (1) canonicalize two equivalent SMILES to one InChIKey and (2) refuse an EvidenceEdge without a grade. Do not download ChEMBL yet. Write a 20-line README that quotes the pitch and non-goals.
```

When that merges, start Phase 1 with a CSV of seed names (known gerontology + any resolved watchlist names) and PubChem identifier exchange.

---

*End of specification. Update dates and N’s in REPORT.md, not by silently changing the question.*
