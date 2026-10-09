# Integrated research decisions and molecular catalog

Prepared 2026-10-09 by an AI coding assistant for a personal hobby and learning
project. This snapshot connects the earlier source reviews to **28 defined
molecular graphs and 10 research programs**. It adds navigation, explicit
decision branches and provenance checks; it does not add experimental results
or a new literature search. No owner laboratory work or source review is implied.

Start with [PROGRAMS.md](PROGRAMS.md) to see the unresolved question, the next
step if supported, the alternative if unsupported, and the measurements needed
for each program. The [CSV index](molecule_index.csv) connects every graph to its
original identity, evidence qualification and review. The
[JSON catalog](research_catalog.json) retains both the programs and provenance.

## What the catalog contains

| Graph archive | Total | Published biological references | Confirmed transformation products | Identity reference only | Unvalidated comparisons |
|---|---:|---:|---:|---:|---:|
| [SS-31](../ss31-analog-audit-2026-10-09/README.md) | 9 | 4 | 0 | 0 | 5 |
| [SLU-PP-915](../slu915-oral-comparators-2026-10-09/README.md) | 10 | 3 | 3 | 0 | 4 |
| [MOTS-c / Epitalon](../peptide-delivery-transfer-2026-10-09/README.md) | 9 | 3 | 0 | 1 | 5 |
| Total | **28** | **10** | **3** | **1** | **14** |

These are roles for a molecular comparison, not evidence grades. A published
reference may have a weaker or undesirable endpoint. A confirmed chemical
product has no assumed activity or abundance in humans. A public identity
reference establishes a graph, not an improved treatment. The 14 unvalidated
comparisons remain hypotheses.

The catalog preserves the original qualification text and field name for each
record. Human oral bioavailability, better oral exposure and improved human
oral exposure are distinct source fields, not interchangeable positive results.
SMILES, neutral formula, standard InChIKey and neutral average mass are copied
from the graph archives. Exact/monoisotopic mass, charge and surface-area
conventions remain in their original archives rather than being combined into
a misleading property ranking.

No new caffeine binder or rapamycin analog graph is invented here. Those
programs have source reviews and unresolved design requirements, rather than
an exact candidate panel. Biological-age assessment and reversal of other
stimulants likewise start with a specified endpoint or ligand.

## Questions that decide whether a modification is useful

The following is a synthesis of the linked reviews, not a validated ranking
of experimental programs or a clinical recommendation.

| Program | Resolve first | Why it changes the next step |
|---|---|---|
| Biological age | Desired organ and functional benefit | A molecular-age shift cannot define the missing benefit. |
| Oral SLU-PP-915 | Cause of intact-parent loss | Chemical instability, metabolic loss, dissolution and poor analytical recovery imply different comparisons. |
| Oral SS-31 | Intact exposure plus retained mitochondrial function | Membrane affinity and carrier signal do not establish useful oral delivery. |
| SS-31 duration | Clearance, cell access or persistence of function | Increased stability may leave the actual duration limit unchanged. |
| SS-31 cost | Measured process burden for unchanged parent | A new sequence adds a biological-equivalence question to manufacturing economics. |
| Oral MOTS-c | Retained function and trafficking after a defined change | The reviewed K14Q results show that persistence can accompany diminished function. |
| Oral acetylated Epitalon | Exact cap identity and parent-relative functional effect | Public structure identity cannot substitute for a useful effect. |
| Caffeine reversal | Physiological recognition and active-metabolite coverage | Water binding cannot establish sustained systemic reversal or sleep improvement. |
| Other stimulant reversal | Exact ligand, mechanism and compartment | Recognition and reversal requirements depend on the stimulant. |
| Reproductive-sparing rapamycin | Retained target-organ benefit at matched response | Testicular exposure and reproductive outcomes must then be assessed directly. |

Each row links to a full decision branch in [PROGRAMS.md](PROGRAMS.md), where
primary references accompany the starting evidence. For example, the MOTS-c
review separates an ELISA-based clearance calculation from intact-species
human pharmacokinetics; the rapamycin review separates complex selectivity from
testicular protection. These distinctions must survive any later candidate
comparison.

## Reproducibility and stale-input handling

[input_manifest.json](input_manifest.json) pins the SHA-256 hashes of the three
molecular JSON archives and maps every archive ID to an explicit role.
[programs.json](programs.json) is the authored decision source. The builder
checks the archives and links before writing outputs. A changed archive,
unmapped graph, stale role map, invalid mass, missing qualification or unknown
program graph fails the build. Updating a hash requires an explicit review;
it does not automatically validate a new or changed graph.

Only the three graph archives are hash-pinned. The review documents and external
URLs are navigation references and are not frozen copies of scientific sources.
Their source-access limits are recorded in the linked studies. Full copyrighted
texts remain in ignored local caches.

From this directory, using the repository's review runtime:

```sh
../../.venv/atlas-review-runtime/bin/python build_catalog.py
../../.venv/atlas-review-runtime/bin/python -m unittest test_catalog -v
../../.venv/atlas-review-runtime/bin/ruff check . --ignore EXE002
```

The eight integration checks cover counts and coverage, verbatim qualification
and identity transfer, altered hashes, new/stale/duplicate records, invalid
program branches, repository path containment and malformed required fields.
They establish catalog consistency, not chemistry or biological efficacy;
the earlier studies contain the structural checks.

No atlas evidence grades, trained predictions, frozen chemical splits or
generator settings are changed. Nothing in this snapshot establishes an
improved human-use molecule, oral formulation, rejuvenation effect, stimulant
antidote or fertility-sparing rapalog.
