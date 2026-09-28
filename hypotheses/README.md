# Hypotheses

Dated run folders from `python -m gen.hypothesis`.

Each folder is a **research artifact**, not a shopping list.

A card must include SMILES, InChIKey, the frozen surrogate score, nearest
Tanimoto to a training active, and a short list of ways the card can be wrong.
It must not include a dose, vendor, stack, or "take this."

Play / monster galleries stay in `artifacts/generated_molecules.parquet` from
`python -m gen.pipeline`.
