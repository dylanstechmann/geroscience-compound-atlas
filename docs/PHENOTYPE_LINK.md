# Citing a phenotype table

When an atlas export is discussed next to a phenotype result, record the
`phenotype_manifest_sha256` from that study's `atlas_link.json`.

The hash identifies bytes. It is not a pChEMBL value, an evidence grade, or a
reason to change a compound card.

The frozen synthetic check that writes this file is
[regen-workbench/studies/frozen_cohort](https://github.com/dylanstechmann/regen-workbench/tree/main/studies/frozen_cohort).
Age is aliased with donor on that fixture. Do not treat its estimand as a
biological age effect.
