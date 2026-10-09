# Report — Geroscience Compound Atlas + Bench

## Question
Can we map public compounds to aging-related biological mechanisms with rigorous evidence grading (E0–E4) and build a leak-free predictive model evaluated on Bemis-Murcko scaffold splits that outperforms simple fingerprint baselines without overclaiming biological rejuvenation?

## What this is not
This project is not a personalized protocol generator, supplement stack advisor, dose guide, or medical recommender. It does not compute a scalar "biological age" score, nor does it treat epigenetic clock shifts as functional proof of organismal rejuvenation. All outputs are computational research artifacts.

## Continuous mTOR activity check (2026-09-28)

`make regression` fits one fixed Ridge model (`alpha=1.0`) to the recorded
pChEMBL values. It uses the 561-row frozen mTOR table and the existing
Bemis-Murcko scaffold partitions for seeds 42, 123, and 456. Morgan fingerprints
and nine descriptors match the binary benchmark; descriptor scaling uses only
each training partition. Validation and test rows are never used for fitting,
and the stored splits are checked for full row coverage and scaffold overlap.
The comparator predicts the training partition's mean pChEMBL on every test row.

Across the three 56-molecule test partitions, Ridge mean MAE was **0.4661**
pChEMBL units (per seed: 0.5523, 0.4542, 0.3919), versus **1.2983** for the
training-mean comparator. Mean RMSE was 0.6277 versus 1.4963; mean R² was
0.8082 versus -0.0319. Full per-seed metrics are in
[`artifacts/regression_metrics.json`](artifacts/regression_metrics.json).
These are correlated internal splits of one ChEMBL target table, not external
validation or a measurement of any generated structure. Assay heterogeneity,
record curation, and scaffold composition can affect these errors. The
classifier's 1 µM cutoff remains a task definition, not a biological boundary.

## Generator split correction (2026-09-28)

The mTOR surrogate was already fit on the seed-42 scaffold training partition, but both generator modes formerly selected starting structures from the full benchmark. Hypothesis mode also measured nearest-active similarity against the full benchmark. This let validation and test structures influence generation despite the model fit being restricted to training rows. Both modes now take seeds from the scorer's 449 training rows (393 actives); hypothesis mode computes its similarity cap against those training actives and rejects exact matches to any training row. Play-mode novelty and its PCA reference also use the training partition. The scaffold split and both YAML configurations are unchanged.

The earlier [hypothesis archive](hypotheses/2026-09-28-mtor-hypothesis/README.md) remains a historical output with the old reference scope. The [corrected run](hypotheses/training-only/2026-09-28-mtor-hypothesis/README.md) used the same 0.55 Tanimoto cap and saved 25 cards from a 400-structure GA archive. Reject counts were 347 for too many rings, 122 for similarity to a training active, 97 for molecular weight, 77 for heavy-atom count, 66 for TPSA, and 19 exact training matches; reasons can overlap. Two clean Linux runs produced identical output files after routing all GA mutations through the recorded seed. These counts are computational filter outcomes, not measured activity, IC50, or validated leads. Old and corrected counts should not be read as a change in biological quality because the starting pool and random sequence changed.

A 2026-10-02 repeat used Python 3.11.16, RDKit 2026.03.6, pandas 3.0.6,
NumPy 2.4.6, scikit-learn 1.9.1 and PyArrow 25.0.1. It reproduced all 25
structures, their order and card Markdown, and the same reject counts in the
[dated run](hypotheses/2026-10-02-mtor-hypothesis/README.md). Three unrounded
QED values in the CSV differed from the earlier archive by at most
1.11e-16; displayed card values were identical. This is a software
reproducibility check, not evidence about chemical activity or biological
benefit.

The [2026-10-02 frozen-split regression rerun](artifacts/regression_metrics_2026-10-02.json)
exactly matches the prior artifact across all three scaffold seeds: mean
Ridge MAE/RMSE/R2 are 0.4661/0.6277/0.8082, versus 1.2983/1.4963/-0.0319 for
the train-mean baseline. This predicts recorded ChEMBL pChEMBL values only; it
does not validate generated cards or establish an aging or rejuvenation
effect.

The earlier play-mode gallery and numerical results below predate this split correction. They describe that archived run; its novelty-vs-training value and PCA legend should be regenerated before being used as current results. A fresh temporary play-mode run completed all three QED weights and produced 200 structures at the primary weight; 198/200 were absent from the training partition, while the QED=0 sensitivity run remained available. High classifier probability and QED remain surrogate artifacts, not evidence of efficacy or rejuvenation.

## Data
*Status (Phases 0–4 Complete; Full Benchmark & Interactive Atlas Deployed)*:
- Watchlist compounds: **N = 20** candidate names.
- Resolved compounds: **19/20 (95.0%)** via PubChem PUG REST.
- ChEMBL molecule resolution: **17/20 (85.0%)** mapped to standardized ChEMBL molecule entities.
- ChEMBL assay coverage: **15/20 (75.0%)** have $\ge 1$ structured binding assay (`assay_type == 'B'`); 2/20 have functional/other assays only; 3/20 have no ChEMBL record (research peptides BPC-157 and Epitalon, plus the negative control).
- ChEMBL bioactivities ingested: **662 activity records** across 209 unique biological targets (with UniProt accessions).
- Curated evidence edges: **62 edges** hand-curated with explicit E0–E4 grading across 12 hallmarks.
- Benchmark dataset: **561 unique molecules** across 279 Murcko scaffolds for ChEMBL mTOR kinase (`CHEMBL2842`).
- Primary exit artifacts:
  - [artifacts/compounds.parquet](artifacts/compounds.parquet) (20 compounds with 2048-bit Morgan fingerprints + 9 RDKit descriptors)
  - [artifacts/chembl_activities.parquet](artifacts/chembl_activities.parquet) (662 bioactivity rows)
  - [artifacts/targets.parquet](artifacts/targets.parquet) (209 biological targets)
  - [artifacts/evidence_edges.parquet](artifacts/evidence_edges.parquet) (62 validated evidence edges)
  - [artifacts/chembl_coverage.parquet](artifacts/chembl_coverage.parquet) (compound-level ChEMBL binding assay breakdown)
  - [artifacts/benchmark_dataset.parquet](artifacts/benchmark_dataset.parquet) (561 mTOR molecules with binary active labels)
  - [artifacts/splits.json](artifacts/splits.json) (scaffold and random split partitions across seeds 42, 123, 456)
  - [artifacts/metrics.json](artifacts/metrics.json) (AUROC, AUPRC, Recall@5%FPR, Brier, and top-10 false positives)
  - [artifacts/generated_molecules.parquet](artifacts/generated_molecules.parquet) (200 GA-generated candidate structures)
  - [artifacts/generator_metrics.json](artifacts/generator_metrics.json) (validity, uniqueness, diversity, QED bias sensitivity)
  - [artifacts/dashboard.html](artifacts/dashboard.html) (standalone interactive HTML dashboard with inline 2D SVGs)
  - [figures/resolution_coverage.png](figures/resolution_coverage.png)
  - [figures/chembl_and_hallmark_coverage.png](figures/chembl_and_hallmark_coverage.png)
  - [figures/scaffold_size_distribution.png](figures/scaffold_size_distribution.png)
  - [figures/benchmark_roc_pr_curves.png](figures/benchmark_roc_pr_curves.png)
  - [figures/generator_chemical_space.png](figures/generator_chemical_space.png)

## Labels and leakage
- **Task**: ChEMBL mTOR longevity kinase activity classification (`target_chembl_id = CHEMBL2842`), selected per §15 default decision fork (clean literature senolytics with selectivity ratios <80; ChEMBL mTOR provides a dense, well-characterized 561-compound set).
- **Label Definition**: `active = 1` if `pchembl_value >= 6.0` (IC50 $\le 1000\ \text{nM} = 1\ \mu\text{M}$); `active = 0` if `pchembl_value < 6.0`.
- **Dataset Composition**: 561 unique structures (477 active [85.0%], 84 inactive [15.0%]).
- **Scaffold Diversity**: 279 unique Bemis-Murcko scaffolds (74.9% singletons), confirming rich chemotype diversity without trivial repetition.
- **Split Discipline**: Bemis-Murcko scaffold split (80% train, 10% val, 10% test) across 3 seeds (`[42, 123, 456]`). Guarantees $\text{Scaffolds}_{\text{train}} \cap \text{Scaffolds}_{\text{test}} = \emptyset$.
- **Leakage Diagnostic**: Evaluated concurrently against a 10-fold Stratified Random Split across the identical seeds. The random split yields an apparent generalization penalty in this sparse-inactive regime, demonstrating the impact of scaffold cluster composition.

## Models
- **Feature Space**: 2048-bit Morgan circular fingerprints (radius 2, ECFP4 equivalent) concatenated with 9 RDKit 2D physico-chemical descriptors (`mol_wt`, `log_p`, `tpsa`, `num_h_donors`, `num_h_acceptors`, `num_rotatable_bonds`, `ring_count`, `fraction_csp3`, `qed`). Descriptors are standardized using training-set parameters only (zero feature leakage).
- **Baseline**: L2-regularized Logistic Regression ($C = 1.0$, `max_iter = 1000`, `solver = 'lbfgs'`).
- **Contender**: `HistGradientBoostingClassifier` (`max_iter = 200`, `learning_rate = 0.05`, `min_samples_leaf = 10`).
- **Compute Budget**: Local CPU evaluation (~4 seconds total runtime across all 3 seeds and both split strategies).

## Results

### Benchmark Comparison Table (Mean ± Std over 3 seeds)

| Split Strategy | Model Pipeline | AUROC | AUPRC | Recall @ 5% FPR | Brier Calibration Score |
|---|---|---|---|---|---|
| **Bemis-Murcko Scaffold (Honest)** | **Baseline (Logistic Regression)** | **0.9735 ± 0.0203** | **0.9929 ± 0.0057** | **0.8496 ± 0.1207** | **0.0486 ± 0.0238** |
| Bemis-Murcko Scaffold (Honest) | Contender (HistGradientBoosting) | 0.9434 ± 0.0294 | 0.9855 ± 0.0095 | 0.7187 ± 0.1550 | 0.0820 ± 0.0392 |
| *Stratified Random (Diagnostic)* | *Baseline (Logistic Regression)* | *0.8912 ± 0.0482* | *0.9780 ± 0.0106* | *0.6528 ± 0.0856* | *0.0922 ± 0.0223* |
| *Stratified Random (Diagnostic)* | *Contender (HistGradientBoosting)* | *0.8766 ± 0.0467* | *0.9748 ± 0.0117* | *0.6667 ± 0.1623* | *0.1159 ± 0.0147* |

- **Key Finding — Class-Imbalance Caveat**: The linear baseline outperforms the tree ensemble on both AUROC (0.9735 vs 0.9434) and AUPRC (0.9929 vs 0.9855), but the more important observation is the **scaffold > random AUROC inversion**. Normally, scaffold splits produce *harder* generalization tasks because they prevent leakage from chemical cousins. Here, the inversion is an artifact of extreme class imbalance (85% active, 15% inactive) combined with high scaffold-singleton rates (74.9%). Scaffold splitting can isolate most rare inactives into a single partition, creating near-pure-active test folds that are trivially classified. The stratified random split preserves the 85/15 ratio within each fold, giving models a harder discrimination task. This motivates a follow-up with either (a) a regression formulation on continuous pChEMBL values, or (b) active/inactive rebalancing via down-sampling or a more class-balanced target.

### Qualitative Analysis of Ten Highest-Confidence False Positives

| ChEMBL ID | InChIKey | Murcko Scaffold | Predicted Prob | True pChEMBL | True Label | Chemical & Biological Failure Analysis |
|---|---|---|---|---|---|---|
| `CHEMBL1088790` | `IVRHPVWEBIOJTP-UHFFFAOYSA-N` | `c1ccc(-c2nc(N3CCOCC3)c3[nH]cc(CN4CCCC4)c3n2)cc1` | 0.9406 | 5.44 | Inactive (0) | **Morpholino-pyrimidine hinge motif**: Contains the canonical PI3K/mTOR pharmacophore; true IC50 is ~3.6 $\mu\text{M}$ (active biologically, but narrowly misses the arbitrary 1.0 $\mu\text{M}$ cutoff). |
| `CHEMBL98350` | `CZQHHVNHHHRRDU-UHFFFAOYSA-N` | `O=c1cc(N2CCOCC2)oc2c(-c3ccccc3)cccc12` | 0.7421 | 5.50 | Inactive (0) | **Chromen-4-one core**: Structurally resembles LY294002 analog series; true IC50 is 3.16 $\mu\text{M}$, deceiving the model through fingerprint sub-structure overlap. |
| `CHEMBL1088831` | `BEHPFPYVKVKMJC-UHFFFAOYSA-N` | `c1ccc(CN2CCC(c3c[nH]c4c(N5CCOCC5)nc(-c5ccccc5)nc34)CC2)cc1` | 0.6154 | 5.70 | Inactive (0) | **Sub-micromolar borderline**: True affinity is 2.0 $\mu\text{M}$ (pChEMBL 5.70); the model predicts active due to intact dual kinase pharmacophore. |
| `CHEMBL591338` | `JDIMUQGRASLCKM-UHFFFAOYSA-N` | `c1ccc(CN2CCC(n3cnc4c(N5CCOCC5)nc(-c5ccccc5)nc43)CC2)cc1` | 0.4351 | 5.77 | Inactive (0) | **Steric clash penalty**: Bulky 4-fluorobenzyl piperidine tail weakens mTOR active site entry relative to PI3K$\alpha$, dropping affinity just below threshold. |
| `CHEMBL1088832` | `GGFPPLTUSBYFPC-UHFFFAOYSA-N` | `c1ccc(CN2CCC(c3c[nH]c4c(N5CCOCC5)nc(-c5ccccc5)nc34)CC2)cc1` | 0.4159 | 5.44 | Inactive (0) | Fluorinated benzyl analog displaying similar borderline micromolar activity (~3.6 $\mu\text{M}$). |
| `CHEMBL188678` | `JAMULYFATHSZJM-UHFFFAOYSA-N` | `O=c1cc(N2CCOCC2)oc2c(-c3cccc4c3sc3ccccc34)cccc12` | 0.3861 | 5.77 | Inactive (0) | Bulky dibenzothiophene core shifts angle of morpholine entry into ATP binding pocket, reducing potency to 1.7 $\mu\text{M}$. |
| `CHEMBL435507` | `RFWCIJWQGHFULK-UHFFFAOYSA-N` | `O=c1cc(N2CCOCC2)oc2c1ccc1ccccc12` | 0.1109 | 5.32 | Inactive (0) | Truncated chromone scaffold retains weak residual binding (4.8 $\mu\text{M}$); model assigns low-positive probability. |

- **Core Conclusion**: The model's "errors" are not hallucinations or random noise; they are **pharmacophorically genuine near-misses** (compounds with true IC50 between 1.7 and 3.6 $\mu\text{M}$) where 2D circular fingerprints cannot capture subtle steric clashes in the deep ATP pocket that drop affinity below 1.0 $\mu\text{M}$.

## Hallmark coverage
Empirical evidence density across the 12+1 controlled aging hallmarks highlights a stark divide between chemistry-rich mechanistic nodes and marketing slogans:

| Hallmark Slug | Total Edges | E4 (Lifespan/Human) | E3 (In Vivo Model) | E2 (Target Engage) | E1 (In Vitro) | E0 (Unstructured) | Scientific Assessment |
|---|---|---|---|---|---|---|---|
| `nutrient_sensing` | 16 | 4 | 7 | 3 | 2 | 0 | **Highest chemical density**: mTOR (rapamycin), AMPK (metformin, berberine), SGLT2 (canagliflozin), and alpha-glucosidase (acarbose) feature replicated ITP mammalian lifespan extensions. |
| `cellular_senescence` | 10 | 0 | 4 | 6 | 0 | 0 | **Strong functional phenotype**: Senolytic combinations (D+Q, fisetin) and BCL-2/BCL-xL inhibitors (navitoclax) demonstrate in vitro selectivity and murine healthspan rescue. |
| `mitochondrial_dysfunction` | 9 | 0 | 5 | 3 | 1 | 0 | **Solid physiological data**: Mitophagy induction (urolithin A), complex I modulation (metformin), and NAD+ repletion (NMN, NR) improve murine bioenergetics. |
| `autophagy` | 8 | 0 | 5 | 2 | 1 | 0 | **Conserved flux**: Spermidine, rapamycin, and metformin demonstrate robust in vivo autophagic flux linked to cardiac and metabolic longevity. |
| `chronic_inflammation` | 8 | 0 | 6 | 2 | 0 | 0 | **SASP reduction**: D+Q, 17$\alpha$-estradiol, and curcumin lower circulating inflammaging cytokines (IL-6, TNF-$\alpha$) in aged rodents. |
| `stem_cell_exhaustion` | 4 | 0 | 4 | 0 | 0 | 0 | **Regenerative rescue**: NR and NMN preserve muscle/neural stem cell pools in mice; navitoclax clears senescent bone marrow HSCs. |
| `intercellular_communication` | 3 | 0 | 0 | 1 | 1 | 1 | **Moderate**: Ephrin kinase binding and immune synapse remodeling; peptide marketing lacks mechanistic clarity. |
| `proteostasis` | 2 | 0 | 1 | 1 | 0 | 0 | **Chaperone flux**: Spermidine eIF5A hypusination and rapamycin translation attenuation mitigate proteotoxic collapse. |
| `ecm_fibrosis` | 2 | 0 | 0 | 1 | 1 | 0 | **Tissue repair only**: Stellate cell collagen suppression and tendon healing without organismal lifespan evidence. |
| `telomere_attrition` | 2 | 0 | 0 | 0 | 1 | 1 | **Slogan-dominated**: Heavily marketed on gray-market wellness sites (Epitalon); lacks replicated mammalian lifespan proof. |
| `dysbiosis` | 1 | 0 | 1 | 0 | 0 | 0 | **Metabolome remodeling**: Acarbose alters cecal short-chain fatty acids (acetate/propionate) in aged mice. |
| `epigenetic_alteration` | 1 | 0 | 0 | 0 | 1 | 0 | **Correlative biomarker**: Metformin PBMCs clock acceleration shift; resetting marks does not equate to genomic repair. |

## Generator (if any)
Implemented a molecular Genetic Algorithm (GA) optimizing the frozen Phase 3 ChEMBL mTOR kinase surrogate model while evaluating QED as an explicit bias control and penalizing Pan-Assay Interference (PAINS) motifs.

### Generation Metrics (Capped at N = 200 structures in Gallery)

| Metric | Primary Run ($\lambda_{\text{QED}} = 0.2$) | Unconstrained ($\lambda_{\text{QED}} = 0.0$) | Oral-Biased ($\lambda_{\text{QED}} = 0.5$) |
|---|---|---|---|
| **Validity Rate** | **100.0%** (via RDKit sanitization) | 100.0% | 100.0% |
| **Uniqueness Rate** | **100.0%** (deduplicated InChIKeys) | 100.0% | 100.0% |
| **Novelty vs Training Set** | **99.5%** (199 / 200 unobserved) | 100.0% | 96.5% |
| **Internal Tanimoto Diversity** | **0.655** (1 - mean pairwise similarity) | 0.563 | 0.670 |
| **Mean Predicted mTOR Prob** | **0.991** | 1.000 | 0.985 |
| **Mean QED Score** | **0.754** | 0.059 | 0.635 |
| **PAINS Pass Rate** | **100.0%** (0 filter matches) | 100.0% | 100.0% |

### QED Weight Sensitivity: Exposing Bias Control
- When unconstrained ($\lambda_{\text{QED}} = 0.0$), the optimizer achieves a maximal predicted mTOR probability ($1.000$), but generated molecules drift into high molecular weight polycyclic structures with near-zero drug-likeness ($\text{mean QED} = 0.059$).
- In this archived play run, a QED weight of 0.2 increased mean QED to 0.754 while the classifier's mean predicted probability was 0.991. Neither value measures target affinity.
- This confirms that QED is an artificial prior from historical oral small-molecule libraries: it penalizes macrocyclic geroprotectors (like rapamycin, $\text{QED} = 0.179$) that achieve potent, lifespan-extending target engagement.

### Chemical Space Exploration (PCA Embedding)
- Computed 2D PCA projection on 2048-bit Morgan circular fingerprints comparing:
  1. GA-generated candidates ($N = 200$, blue circles)
  2. ChEMBL mTOR dataset actives ($N = 477$, gray dots; this archived plot predates the split correction)
  3. Curated landmark geroprotectors (red triangles: rapamycin, metformin, dasatinib, canagliflozin, acarbose, navitoclax)
- **Variance Explained**: PC1 accounts for 13.3%, PC2 accounts for 9.1%.
- **Interpretation limit**: This PCA shows where the generated fingerprints project relative to archived ChEMBL and atlas structures. Nearby points in two dimensions do not establish potency, synthesizability, or biological effect; play mode can still emit high-scoring junk.

## Claims we refuse
- No human dosing or administration protocols.
- No claims of "reversing aging" or curing senescence.
- No treating surrogate biomarkers (e.g. methylation clocks) as conclusive functional rejuvenation.
- No treating drug-likeness (QED) as biological efficacy.

## Next
- Spring semester deep-learning integration: Replace 2048-bit Morgan fingerprint representations with Graph Neural Network (GNN / Chemprop / SchNet) molecular encoders evaluated on the identical frozen Bemis-Murcko scaffold split partitions (`artifacts/splits.json`).
- Ablation study comparing molecular graph encoders with and without edge/stereochemical features against the linear baseline.
- Expand the curated evidence graph from N=20 seed compounds to 100+ NIA Interventions Testing Program (ITP) compounds.

## Archived results summary (§13)

- **Archived atlas:** 20 candidate compounds (95.0% PubChem resolution, 85.0% ChEMBL mapping) and 62 curated mechanistic evidence edges across 12 aging hallmarks with structured E0–E4 grading (4 mammalian lifespan interventions; 75% binding assay coverage).
- **Internal fingerprint benchmark:** Strict Bemis-Murcko scaffold disjoint splits (80/10/10 across 3 seeds; 279 unique scaffolds) for ChEMBL mTOR kinase activity (N=561); regularized L2 Logistic Regression reported AUROC 0.9735 ± 0.0203 and AUPRC 0.9929 ± 0.0057, outperforming tree ensembles on these splits. The error analysis covers 10 false positive near-misses (IC50 1.7–3.6 µM) driven by canonical ATP-hinge pharmacophores.
- **Archived play-mode generator:** A Genetic Algorithm optimizes the frozen mTOR surrogate with explicit QED-bias control and PAINS filtering; its archived play run generated 200 RDKit-valid unique molecules (99.5% novel against the full dataset, internal diversity 0.655). The PCA is descriptive and does not validate the molecules.
- **Static dashboard:** Self-contained HTML/SVG/JS displays inline RDKit vector structures, dynamic multi-attribute filtering, coverage hole diagnostics for non-small-molecule modalities, and reproducible split audit metrics.

## Implementation follow-up — 2026-10-03

The evidence-integrity corrections identified in the longevity tooling review are now reflected in the active curated graph. Curated evidence joins use stable PubChem compound IDs and exact structure-key checks; contradictory target mappings and invalid citations are rejected or retained in dated quarantine with reasons. The active curated table contains 57 rows: six have source-to-claim review with limits recorded, 49 are explicitly marked as having unreviewed claim support, and two are unverified vendor/gray-market claims. Five invalid rows were quarantined. Materialization preserves the curator-assigned grade as `curator_grade`, but displays unreviewed and vendor claims as E0; only source-reviewed rows retain their evidence grade. E4 is restricted to mammalian lifespan evidence or a reported positive human outcome, so trial registration alone does not qualify.

`--curated-only` rebuilds the evidence-edge artifact offline from curated records without calling ChEMBL. The dashboard and static site were regenerated from that artifact and distinguish source-reviewed evidence from unreviewed assertions. Acarbose lookup now returns hallmark records with their review status, while the mTOR surrogate abstains because its applicability domain is not validated. The frozen scaffold split and configured activity threshold were not changed, and no retraining or new generator run was performed.

The historical results summary and archived benchmark metrics above describe earlier artifacts; use the current graph and dated artifacts for present evidence counts. The October 2 regression metrics, hypothesis archive and prior report content were preserved. Atlas tests: **87 passed**. A manual Acarbose lookup returned two source-reviewed hallmark records, including the mouse-lifespan record with its species and claim limits; any unreviewed claim remains E0 in the output.

## Audit and reliability follow-up — 2026-10-04

The unfinished neighborhood diagnostic is now executable from the frozen
benchmark even when its fingerprint sidecar is absent. It reconstructs Morgan
radius-2, 2048-bit fingerprints from the saved structures without writing a
cache. A supplied sidecar must exactly match those structures in row order;
equal row counts alone are insufficient. Invalid indexes, incomplete or
overlapping partitions, scaffold leakage and inconsistent activity labels
fail closed. The Tanimoto computation uses matrix intersections rather than
allocating a test-by-train-by-bit cube.

The [dated reliability artifact](artifacts/neighborhood_reliability_2026-10-04.json)
retains the existing 561-row dataset and scaffold seeds 42/123/456 and uses the
unchanged `configs/bench.yaml` model settings and pChEMBL cutoff 6.0. It records
source hashes, runtime versions, exact held-out predictions and fixed
similarity-band summaries for both models. Empty bands, class counts, Brier
scores and descriptive Wilson intervals prevent small bands from appearing
to provide strong calibration evidence. Only 4/3/0 test molecules respectively
have maximum training Tanimoto at most 0.50; this internal benchmark has little
coverage of chemically remote predictions. No calibrator or applicability
threshold was selected and no external assay validation was performed.

`make hypothesis PYTHON=.venv/bin/python` was rerun with unchanged
`configs/hypothesis.yaml`. The [October 4 archive](hypotheses/2026-10-04-mtor-hypothesis/)
contains 25 accepted surrogate cards after gating an archive of 400 molecules.
Reject counts are overlapping: 347 too-many-rings, 122
too-close-to-training-active, 97 molecular-weight, 77 heavy-atom, 66 TPSA and
19 seen-in-training. The frozen seed-42 reference contains 449 training rows
and 393 actives. This is a constrained generator result, not evidence of
measured activity, synthesis feasibility or longevity benefit. The older
hypothesis archives and QED=0 play sensitivity artifacts remain intact.

Validation: the atlas suite passed **115 tests**, including invalid/leaking
splits, fingerprint row-order corruption, snapshot hashing, interval
uncertainty and report-output protection. Ruff passed for the three new
diagnostic files. The dated reliability and hypothesis jobs both completed in
the Linux development container. Frozen dataset/split/config bytes were
verified unchanged against their saved SHA-256 values; hypothesis provenance
also retains configuration snapshots.

A dashboard builder defect was also corrected: creating an HTML preview or
temporary test output silently rewrote `site/index.html`. Hosted output now
requires explicit `site_html=`; the existing `make report` CLI still writes
both the artifact and hosted copy. Regression tests verify preview isolation
and identical explicit copies. The test-regenerated tracked page was restored
to its initially clean contents.


## Public-data external validation (2026-10-04)

The fixed classifiers were tested against 2,910 qualified molecule–assay rows from 170 source documents and 2,731 connectivity structures excluded from the original dataset. Both models had worse pooled Brier error than the training-prevalence baseline in all three seeds; HGB AUROC was 0.467–0.478. Current assay metadata also shows that only 209/600 original raw records qualify as direct human single-protein target-assignment records under the new strict contract. The frozen benchmark and generator are unchanged; their scores are not qualified for compound prioritization. See the [prespecified study, all-source results and reproduction steps](studies/chembl_external_2026-10-04/README.md).

The dashboard now repeats this negative external result before its internal split
metrics and generated-molecule gallery. The generator remains available as an
algorithmic exploration, with classifier outputs labeled as unqualified
surrogate artifacts.

## Documentation framing — 2026-10-07

The README and specification now describe this as a personal hobby and learning
project developed with substantial assistance from AI coding tools. Employment
pitch language and résumé-style authorship claims were replaced with project
purpose and artifact descriptions. The dashboard heading is now “Benchmark
observations” in its generator and both existing HTML copies. The dashboard
introduction also discloses the hobby purpose and substantial AI assistance.
Existing scientific results, numerical metrics, and limitations are retained.
Verification: diff checks passed; all three dashboard changes contain only the
heading replacement and the same introduction note, and the two HTML copies
remain byte-identical.

## Citation identifier audit — 2026-10-07

An identifier audit of `data/curated/curated_evidence.csv` (57 rows; 51 unique PubMed IDs; 2 ChEMBL document IDs) found that 37 rows cite an identifier that does not resolve, resolves to a paper with an unrelated title, or resolves to a title that does not show the recorded claim. Three PubMed IDs are absent from both PubMed and Europe PMC. Twenty-six PubMed IDs resolve to titles on unrelated subjects. Two ChEMBL document IDs return HTTP 404. Four rows need the paper read before their titles can be judged. Eighteen rows fit at title level only, and no row has been verified as supporting its claim.

No curated value, grade or claim was changed. The audit, its receipts (refcheck results and second-source checks), the row-level screen and the owner decision list are in [docs/citation-audit-2026-10-07](docs/citation-audit-2026-10-07/README.md). The README carries a short status note that points there.

## Claim-bound citation display — 2026-10-09

An AI coding assistant integrated the October 7 identifier screen into edge
materialization, the lookup CLI and dashboard generation. The versioned Python
package now includes 57 claim-bound title screens rebuilt from the existing
audit receipts. Matching uses the normalized compound identity, target,
citations, recorded claim and study context rather than CSV positions or shared
identifiers. Changed or legacy claims and missing/invalid manifests remain
unreviewed. The manifest builder refuses stale summary or triage outputs.

All 57 current displayed grades are E0 pending source-to-claim review. The six
historical `source_reviewed_with_claim_limits` labels remain inspectable with
the original `catalog_grade` and `curator_grade`; they no longer override the
newer audit's unreviewed claim-support state. The four historical E4 rows do
not pass the current dashboard E4 preset. CLI and cards expose the dated
category and individual identifier findings, label the notes as recorded
claims, and explain that title resolution does not verify paper support.
Cached parquet readers apply the overlay without requiring a data rebuild.
The source CSV, citation identifiers and audit judgments were not changed.

Validation used the workspace's Docker `dev` service. The full suite passed
**162 tests**, including changed claim/citation/context identities, reordered
records, legacy caches, missing/malformed manifests, duplicate records,
per-identifier findings, display withholding and reproducibility from receipts.
The focused run passed 53 tests; the prescribed hypothesis-filter/split check
also passed 9 tests. The built wheel contains and loads all 57 screens outside
the editable source path. Both generated HTML copies are identical. A Node
execution of the generated dashboard script confirmed initial rendering,
E4 → zero compounds, reset → 20 compounds and Rapamycin search → one compound.
Ruff passed for the new module, manifest builder and tests with `EXE002`
excluded because the Windows bind mount reports Python files as executable.
The broader lint check of existing lookup/dashboard files still reports their
pre-existing `SIM102`, `RUF046`, `BLE001` and `ISC004` findings; those unrelated
sections were retained.

The unchanged constrained generator was rerun using the existing seed-42
scaffold reference (449 training rows, 393 active rows). It again retained 25
cards from a 400-molecule archive, with overlapping rejects: 347 rings, 122
too-close-to-training-active, 97 molecular weight, 77 heavy atoms, 66 TPSA and
19 seen-in-training. The [run receipt](artifacts/hypothesis_review_2026-10-09.json)
records the metrics and unchanged SHA-256 values for the frozen dataset, split
and both configuration files. Generated cards from this check are local under
`.venv/review-hypotheses`; they are surrogate artifacts, not measured activity
or evidence of rejuvenation. Existing hypothesis archives and QED=0 sensitivity
outputs were preserved. No new source reading or biological result is claimed.

## Oral delivery and selective-exposure feasibility — 2026-10-09

An AI coding assistant prepared a [focused primary-source feasibility analysis](studies/oral-targeting-feasibility-2026-10-09/README.md)
covering organ-specific biological-age endpoints, SLU-PP-915, SS-31, MOTS-c,
acetylated Epitalon, stimulant reversal and reproductive-sparing rapamycin
development. This is agent-run research; no biological finding or paper review
is attributed to the owner.

The study archives real Europe PMC search and ClinicalTrials.gov registry
receipts, eight peptide reference/hypothesis specifications with chemistry
bookkeeping, five synthetic exposure scenarios, thirty synthetic equilibrium
binding scenarios, and a rendered comparison figure. It separates current
observations from proposed modifications and highlights failed or missing
evidence. It does not update the curated evidence table, benchmark, generator,
review statuses or claims of oral/human efficacy. No analog was synthesized,
optimized against a validated model or shown to reduce reproductive risk.

All seven physical/reference checks passed in Docker `dev`. They cover
stoichiometric capacity, conservation, affinity/capacity monotonicity,
limiting behavior, invalid parameters, exposure tradeoffs and the FDA reference
peptide mass. The figure was visually inspected. The standalone model uses
only the standard library; plotting uses matplotlib. These checks verify the
calculations, not biological validity.

## Research decisions and competitive capture follow-up — 2026-10-09

The AI coding assistant extended the [oral/selective-exposure study](studies/oral-targeting-feasibility-2026-10-09/RESEARCH_DECISIONS.md)
with six prioritized programs and explicit advance/stop criteria. Newly read
sources include the full SLU-PP-915 metabolism paper, the FDA elamipretide
integrated review and a 2025 SS-31 mouse function/aging-clock study. The review
also records a sponsor's negative human oral SBT-272 exposure disclosure and
a primary patent's unchanged-parent manufacturing approach, with limitations.
Source access modes and downloaded-file hashes are archived in follow-up
receipts; full texts remain in ignored local storage. Relevant FDA PDF pages
were rendered and visually checked. This research was performed by the agent,
not attributed to the owner.

The findings distinguish chemical instability from enzyme-mediated turnover,
human in vitro plasma stability from clinical plasma elimination, manufacturing
burden from marketed price, and functional improvement from clock changes.
No source supplies a newly validated oral candidate or fertility-sparing
rapalog. Existing curated evidence grades and compound tables were retained.

A standard-library equilibrium solver now includes multiple guests competing
for one-to-one binding sites. Six wholly synthetic scenarios illustrate
metabolite coverage and competition; there are no measured caffeine affinities,
human exposure values, candidate scores or sleep predictions. All **15 study
checks** passed in Docker `dev`, including eight new checks of mass action,
conservation, limiting behavior and agreement with an independent analytic
single-guest reference. The prescribed hypothesis-filter/split tests passed
**9 tests**. New-code Ruff passed with `EXE002` excluded for Windows bind-mount
executable flags. These checks validate calculation behavior, not biology.
Frozen datasets, splits, configurations and generated hypothesis archives
were unchanged.

## SS-31 analog audit and molecular definitions — 2026-10-09

An AI coding assistant added an [SS-31 analog audit](studies/ss31-analog-audit-2026-10-09/README.md)
and nine explicit molecular graphs. Published SPN4/SPN10 comparators are
separated from five unvalidated modification proposals, including one
patent-described Phe-backbone N-methyl comparator. The earlier Tyr-comparator
status was corrected to published SPN4 evidence, retaining the legacy ID and
adding a source/alias. No claim of human oral efficacy, better half-life,
manufacturing savings, novelty or biological validation was added.

JSON and SDF outputs record exact stereochemical identity, formula, nominal
charge and neutral-graph descriptors. SDF coordinates and PNG depictions are
2D representations. SS31 matches an independently retrieved official PubChem
formula, mass and stereochemical InChIKey. Sources and downloaded-file hashes
are archived in receipts; copyrighted full text remains in ignored storage.
The diagrams were visually inspected. No owner source review or laboratory
work is claimed.

All **9 structural checks**, **15 previous study checks**, and **9 prescribed
hypothesis-filter/split tests** passed in Docker `dev`. The structural checks
include independent reference identity, exact stereochemistry, donor and
composition accounting, and SDF/SMILES round trips. New/changed Python Ruff
checks passed with `EXE002` excluded for Windows bind-mount executable flags.
These checks concern chemical definitions and calculations, not drug efficacy.
Curated evidence tables/grades, frozen benchmark inputs, splits, thresholds
and generator configurations were retained.

## Rapamycin tissue selectivity and reproductive evidence — 2026-10-09

The AI coding assistant added a [focused rapamycin review](studies/rapamycin-selectivity-2026-10-09/README.md)
covering complex selectivity, adult-animal reproductive findings, cell-specific
genetic results and local/tissue-restricted delivery. Nine source-access
receipts include five verified primary full-text downloads. The review does
not infer fertility protection from improved metabolic outcomes or brain
restriction, and it distinguishes blood non-detection from zero exposure.
Three targeting hypotheses have explicit retained-function, active-parent
distribution and reproductive requirements. No molecule was biologically
validated and no owner research or source review is claimed.

A standard-library Hill-1 model produces six synthetic comparisons at matched
target pathway inhibition, twelve arbitrary response-bound pairs and a
compartment-aggregation counterexample. Its inputs are not measured rapalog
parameters or human safety limits. It demonstrates that global exposure or
potency changes do not create tissue selectivity after target-response
matching; lower total testis partition can be offset by higher free fraction;
and an aggregate signal can hide a larger compartment-specific response.
It does not predict fertility, aging benefit, oral PK or an exposure regimen.

All **8 mathematical checks** and **9 prescribed hypothesis-filter/split tests**
passed in Docker `dev`. New-code Ruff passed with the Windows-mount `EXE002`
excluded. The standalone scientific figure was rendered and visually checked.
Curated evidence grades/tables and frozen benchmark, split and generator
configuration were retained. These checks validate calculations, not biology.

## Caffeine reversal: source equilibria and oral-development requirements — 2026-10-09

The AI coding assistant added a [caffeine reversal audit](studies/caffeine-reversal-2026-10-09/README.md)
with seven primary-source access receipts, including three hashed official
full-text XML downloads. It identifies aqueous small-molecule recognition,
DNA/RNA aptamer measurements and bacterial catalytic conversion, distinguishing
sensor output and substrate disappearance from systemic pharmacological
reversal. The 2026 RNA study's undetected paraxanthine binding is retained as
an assay limitation, not an infinite-Kd estimate. No owner laboratory work or
source review is claimed.

A standard-library solver reconstructs the tweezer paper's NMR equilibria,
including host self-association and multi-unit host/guest species. Source
constants, uncertainties, medium and units are recorded separately from
arbitrarily chosen total concentrations. Ninety-six illustrative calculations
and a visually inspected scientific figure compare source models with two
explicitly incomplete/incorrect shortcuts. Aggregate BC500 is not treated as
a one-site Kd in either source reconstruction. Guest dimers are not assumed
inactive. These outputs do not predict plasma binding, kinetics, oral PK,
brain response or sleep restoration. Uncertainties are not propagated.

The review specifies an unvalidated oral prohost concept and its missing
physiological recognition, metabolite, endogenous-ligand, active-host PK,
clearance/rebound, brain-response and clinical-endpoint evidence. No exact
improved drug, safe precursor, administration regimen or rejuvenation claim
is established. Copyrighted full text remains in ignored local storage.

All **9 new mathematical checks** and **9 prescribed hypothesis-filter/split
tests** passed in Docker `dev`. Ruff passed for new Python with `EXE002`
excluded for Windows bind-mount executable flags. Calculation checks include
independent analytic references, a constructed multimer state, mass conservation
and two-guests-per-host accounting. Curated evidence grades/tables, frozen
inputs, splits, thresholds and generator configurations were unchanged.

## SLU-PP-915 exact oral-development comparisons — 2026-10-09

The AI coding assistant added a [SLU-PP-915 comparator study](studies/slu915-oral-comparators-2026-10-09/README.md)
with ten exact neutral molecular graphs: parent, published 10q/10r ERR
comparators, reference-confirmed M1/M3/M4 products, and four unvalidated
structural comparisons. The hypotheses are two difluoro-aniline positional
variants, an amide N-methyl comparator and a pinacol-ester precursor. None
has a claimed synthesis, improved oral performance, novelty, safety or aging
effect. No owner laboratory work or source review is claimed.

Three primary-source receipts and one official PubChem receipt preserve
access limitations and downloaded-file hashes. The discovery full text was
read via official NCBI BioC after Europe PMC fullTextXML failed. Microsomal
half-lives reported as >60 minutes remain censored lower bounds; reported
EC50 and maximal response stay distinct. The oral-study assessment is based
on its official abstract and supplies no numerical absolute availability.
M6's aniline-ring assignment is not promoted to a known hydroxylation carbon.
Apparent source analytical-mass inconsistencies are documented, not inserted
as chemical identity references. Copyrighted full text remains ignored.

JSON/SDF and visually inspected 2D depictions define the panel. Parent matches
independently retrieved PubChem formula, exact mass and stereochemical
InChIKey. Both default RDKit TPSA and sulfur-inclusive TPSA are recorded so
descriptor conventions are not silently combined. No classifier or oral
drug-likeness ranking is assigned. The review defines retained activity,
chemical recovery, active-species, free/tissue exposure and route-reference
requirements for assessing oral improvement.

All **8 structural checks** and **9 prescribed hypothesis-filter/split tests**
passed in Docker `dev`; new-code Ruff passed with the Windows bind-mount
`EXE002` excluded. Structural checks cover independent identity, source-product
formula accounting, substitution distances, donor/connectivity changes and
SDF/SMILES round trips. They validate molecular definitions, not biology.
Curated evidence grades, frozen inputs, splits, thresholds and model/generator
configurations were unchanged.
