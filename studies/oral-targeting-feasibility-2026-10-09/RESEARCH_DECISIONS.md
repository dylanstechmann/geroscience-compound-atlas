# Research decisions after the initial feasibility review

Prepared 2026-10-09 by an AI coding assistant for this personal hobby and
learning project. This continuation adds primary-source review and a
reproducible competition calculation. No compound was synthesized or tested
biologically. Priorities reflect how quickly a comparison could resolve
uncertainty; they are not efficacy, commercial-value or safety rankings.

## Findings that change the decisions

**SLU-PP-915: check chemical stability before optimizing metabolism.** The full
2026 paper reports amide hydrolysis, boronic-acid conversion to a phenol and
hydroxylation. M1, M3 and M4 also appeared in enzyme blanks; their identities
were supported by reference compounds. Aniline-ring hydroxylation was assigned
for M6, but a dominant human clearance pathway was not established.
[Full metabolism paper, In Vitro Metabolic Transformation and Table 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC12835572/).
The design inference is to distinguish chemical breakdown from enzyme-dependent
turnover before selecting a substitution. Neither an enzyme inhibitor nor a
particular replacement group is justified by these observations.

**SS-31: human plasma stability is already much longer than systemic residence.**
FDA's integrated review reports in vitro human-plasma half-lives of 26.1 and
30.8 hours in its two test conditions. These estimates came from a six-hour
incubation study, so they are extrapolations rather than day-long observations.
The visually checked human IV Table 103 reports mean plasma elimination
half-lives of 3.23-3.93 hours across cohorts, with six active participants per
cohort. These are different experiments and cannot be divided to apportion
clearance mechanisms. They show why resistance to plasma proteases alone
is an inadequate duration goal. [FDA integrated review, printed pages 181-182
and 186-187](https://www.accessdata.fda.gov/drugsatfda_docs/nda/2025/215244Orig1s000IntegratedR.pdf).
The approved label also documents urinary elimination and C-terminal
degradation. [Label, section 12.3](https://www.accessdata.fda.gov/drugsatfda_docs/label/2025/215244s000lbl.pdf).

**A second-generation comparator encountered the human oral-exposure barrier.**
Stealth's 2021 annual filing states that orally administered SBT-272 did not
reach the desired drug exposure in its completed human phase-I study, despite
earlier animal oral-availability findings. This is a sponsor disclosure, not a
peer-reviewed dataset or a quantitative oral-bioavailability estimate. It does
not establish that every formulation would fail. [Sponsor filing, SBT-272
Safety](https://www.sec.gov/Archives/edgar/data/1696396/000095017022005474/mito-20211231.htm).
The inference is that improved animal exposure and a mitochondria-targeting
peptidomimetic label cannot establish successful human oral delivery.

**Cheaper SS-31 may be a process problem before it is a molecule problem.**
Patent EP3160984B1 describes a solution-phase route with crystallization that
the applicants say avoids final HPLC purification and freeze-drying. This is a
concrete process-development comparator using the same covalent peptide,
not proof of present manufacturing cost, independent process reproducibility
or lower marketed price. [Primary patent, summary and process description](https://patents.google.com/patent/EP3160984B1/en).
Compare material utilization, isolated yield, purification burden, impurity
control and stability. Changing Dmt to Tyr additionally requires biological
equivalence to be re-established.

**Functional benefit can occur without an improved aging-clock result.**
The 2025 peer-reviewed mouse study measured function and molecular clocks after
eight weeks of SS-31. Some functional outcomes improved; cardiac epigenetic
age did not, and transcriptomic-age improvement was absent in most groups.
Young male hearts were an exception for the transcriptomic clocks, so "no
clock changed anywhere" would overstate the finding. Molecular analyses used
four to five animals per group, limiting sensitivity. This is not human age
reversal or proof that clocks are useless. [Primary study, Results and
Future Directions and Limitations](https://onlinelibrary.wiley.com/doi/10.1111/acel.70026).

A separate pilot in aged female mice found preserved exercise tolerance with
intermittent SS-31 treatment, while several other outcomes did not reproduce
benefits reported with continuous delivery. It did not establish a human
schedule or compare a long-lived analog against parent at matched exposure.
[Intermittent-treatment experiment](https://pmc.ncbi.nlm.nih.gov/articles/PMC10651577/).
Measure persistence of function after exposure falls rather than assume longer
plasma residence is always necessary. That is a proposed endpoint, not a
washout result established by this intermittent-treatment pilot.

## Proposed order of work

| Program | First useful comparison | Evidence needed to advance | Result that stops or redirects the idea |
|---|---|---|---|
| SLU-PP-915 oral parent | Parent loss in formulation and enzyme-free controls versus active metabolic systems; then parent formulations | Reproducible intact-parent recovery and unbound exposure, preserved ERR response, identified products | Apparent metabolism is chemical breakdown; increased total signal is a product or artifact |
| SS-31 duration and cost | Parent versus terminal D-Phe comparator for function and degradation; separately compare parent manufacturing processes | Useful intracellular exposure or durable function, acceptable impurities; measured process economics | Longer circulation with less mitochondrial access; unchanged clearance; cheaper chemistry with lost activity |
| Caffeine capture | Parent plus active metabolites and physiological competitors | Independent binding measurements, residual activity, intact-host systemic availability and complex disposition | Strong binding only in buffer; active metabolites remain; competitors consume capacity |
| Reproductive-sparing rapamycin | One declared target organ versus testis at comparable target-organ pathway response | Better functional selectivity with retained benefit, plus independent reproductive endpoints | Less exposure everywhere; favorable total-tissue ratio with unchanged testicular inhibition |
| MOTS-c oral development | Establish active species and parent functional reference first | Exact identity, reproducible function and measured parent exposure | Stable analog or fragment has a different response; endogenous observations replace treatment evidence |
| Acetylated Epitalon / Adamax | Verify identity and comparative activity/PK evidence | Confirmed structure, meaningful functional assay and parent-versus-variant comparison | Salt is mistaken for a covalent cap; telomere length or vendor claims are the sole efficacy basis |

SLU-PP-915 and SS-31 have defined parents and concrete early comparisons.
Caffeine capture remains a discovery program. Rapamycin targeting needs a
chosen organ before "better selectivity" has an interpretable meaning.
MOTS-c and acetylated-peptide work need stronger starting assays or identity
evidence. This order is a research-management judgment, not a recommendation
to use a compound.

## What better means

For **oral delivery**, separate intact-parent absolute bioavailability from
plasma concentration, absorption rate and terminal half-life. Distinguish
parent from active products. Use appropriate route references, chemical
recovery and uncertainty; a larger peak alone is insufficient. Stability
and permeability calculations are screening measurements, not human exposure.

For **SS-31 duration**, separate plasma residence, free-parent exposure,
intracellular mitochondrial access and persistence of function. The D-Phe
comparator changes one stereocenter while preserving formula and nominal
charge. It tests whether C-terminal degradation can be reduced while retaining
function; it is not an optimized drug. Reversible albumin attachment must
establish parent release and cellular access. Conjugate signal cannot be
counted as active parent.

For **cost**, compare unchanged-parent processes first. Fewer expensive
purification steps could reduce manufacturing burden without the biological
revalidation required for a new sequence. This does not quantify savings.
Low systemic oral availability can increase material requirements, so route
convenience and production cost need separate accounting.

For **rapamycin selectivity**, testicular volume, sperm function and endocrine
effects are distinct endpoints. Measure testicular mTORC1 and mTORC2 responses
separately from reproductive outcomes. Retained target-organ function with
less testicular suppression would support distribution selectivity, but still
requires direct reproductive evaluation. The initial review contains the human
reproductive and mouse mechanistic sources.

For **biological age**, retain parallel outcomes: organ function, individual
molecular/clock measurements, adverse effects and persistence. A null clock
result does not erase demonstrated function, and a favorable clock does not
establish function or organism-wide reversal. Do not collapse these into a new
composite rejuvenation score.

## Competitive binding calculation

The new [standard-library model](competitive_binding.py) extends the original
single-guest calculation to guests competing for one class of one-to-one
binding sites. Every total and dissociation constant uses the same arbitrary
scale. None is a measurement of caffeine, metabolites or a candidate host.

```text
C_i,free = C_i,total * Kd_i / (Kd_i + H_free)
C_i,bound = C_i,total * H_free / (Kd_i + H_free)
H_total = H_free + sum(C_i,bound)
```

The solver finds the unique free-host concentration between zero and total
host. Near saturation it preserves small free concentrations rather than
subtracting two rounded, nearly equal occupied-site totals.

| Purely synthetic scenario | Target remaining free | Other guest remaining free |
|---|---:|---:|
| Target alone | 9.51% | None |
| Equal amount of an equally well-bound second guest | 50.49% | 50.49% |
| Equal amount of a poorly bound second guest | 9.95% | 99.10% |
| Abundant, tighter competitor | 99.48% | 95.03% |
| Abundant, weaker competitor | 36.28% | 98.27% |
| Twice the host capacity with the tighter competitor | 98.91% | 90.06% |

Exact assumptions are in [the six-scenario JSON](synthetic_competitive_binding.json).
These counterexamples show why metabolite coverage and competitor selectivity
matter. Capturing an active metabolite consumes sites; failing to capture it
can leave activity outside the modeled target. A weaker competitor can matter
when abundant. Increasing capacity alone may do little in an unfavorable mix.

There are no brain compartments, rate constants, renal elimination, receptor
responses, circadian phase, oral absorption or sleep architecture. The model
cannot estimate a human amount, required affinity, speed of reversal or
benefit. "Active_metabolite" is a hypothetical role, not a potency assignment.

## Retrieval and validation

[Follow-up receipts](followup_source_receipts.json) record access modes,
download hashes and limitations. Full SLU metabolism and 2025 aging-study XML
were read through the official Europe PMC service. The FDA PDF has 378 pages;
relevant complete pages were extracted, and PDF pages 251 and 256 (printed
182 and 187) were rendered and visually inspected. The targeted FDA search
did not supply an oral elamipretide result; oral aspirin coadministration was
a different experiment. This bounded absence does not prove that the entire
development archive lacks oral data.

Discovery and intermittent-study XML endpoints returned HTTP 500. Discovery
claims remain as cited in the initial analysis. The intermittent finding uses
indexed primary-paper text; the SBT-272 statement uses the indexed sponsor
filing. No numerical human oral PK dataset was retrieved. Downloads and PDF
renders remain in ignored local storage; receipts do not redistribute full
text. No biological work or source review is attributed to the owner.

Reproduce in Docker `dev`:

```bash
cd geroscience-compound-atlas/studies/oral-targeting-feasibility-2026-10-09
python3 competitive_binding.py
python3 -m unittest test_analysis_models test_competitive_binding -v
```

All **15 study checks** passed. Eight new checks cover an independent analytic
reference, conservation, capacity, mass action, equal-affinity behavior,
competition, capacity monotonicity, empty cases, extreme-affinity stability
and invalid inputs. The prescribed hypothesis-filter/split check passed
**9 tests**. New-code Ruff passed with `EXE002` excluded for executable Windows
bind-mount behavior. These validate calculations, not drug designs.

No curated grades, target tables, frozen splits, benchmark thresholds or
generator settings changed. The output supports research decisions, without
claiming a proven oral drug or a fertility-sparing rapalog.
