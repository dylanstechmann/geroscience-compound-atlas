# Hypothesis run 2026-09-28

Research artifact only. Not a drug, supplement, dose, or protocol. A high mTOR surrogate score is not an IC50. QED is a historical oral-drug prior.

Propose valid small molecules that the frozen ChEMBL mTOR logistic
surrogate scores as active, after rejecting property-space spam and
near-duplicates of the training actives.

Accepted cards: 25
GA archive before gates: 400
Anti-clone Tanimoto cap: 0.55

Reject counts:
- too_close_to_training_active: 207
- too_many_rings: 155
- mol_wt_high: 85
- too_many_heavy_atoms: 71
- tpsa_high: 58
- seen_in_training: 21
- logp_high: 13
- too_flat: 4

## What this run actually says

The cap worked. 207 archived graphs were still too close to a training active, and 21 were the training molecules themselves. The 25 cards that passed have nearest Morgan Tanimoto between about 0.29 and 0.42, and surrogate probabilities from about 0.93 to 0.99.

That high probability is the optimizer doing its job on a weak classifier. It is not an IC50, not a lead list, and not something to make. A chemist can still reject a card for a weird sulfur/nitrogen graph that RDKit sanitized.
