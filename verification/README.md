# Verification scope

Run from the repository root:

    python verification/check.py

Python 3.10+ with NumPy and SciPy, and Lean 4.30.0, are required. The scripts are
bounded finite checks, with per-command timeouts in the runner.

## Included checks

- **houghton_labels.py:** exact integer/rational bookkeeping of the 104-dimensional
  gadget: 33 symbols, 35 triangles, 174 entries, 26 distinct group labels,
  singleton anchors and total squared coefficient weight 104. A length-L ray
  word is represented by its action up to height L and its eventual translations.
  The reason this representation is complete is mathematical, not a Lean theorem.
- **gadgetcheck.py:** numerical unitary/Kraus/entropy identities for the S3, Z4,
  and Z2 x Z2 finite-group gadgets. These examples do not prove the general
  group-gadget theorem.
- **tagcheck.py:** numerical verification of the tagging construction and its
  one-half Choi-distance identity in small stand-in dimensions. It does not
  directly construct the 208-dimensional hard channel's infinite factor.
- **accuracycheck.py:** Choi distance, support leakage and product-state distances
  for the rotating-dephasing example. The all-observer adaptive equality follows
  from the manuscript's common-processor proof, not finite numerical sampling.
- **Lean Resources.lean:** exact dimension arithmetic, telescoping integer
  replacement schedules, and a conditional integer budget inequality. The integer
  schedule applies to dyadic rank; the real-logarithmic general schedule and
  entropy budgets are not formalized.

The Lean theorems contain no proof placeholders, native proof-evaluation oracle,
or added mathematical axioms. Their standard logical dependencies are printed
when run.

## External and unformalized dependencies

The companion [group paper](https://github.com/Apsiape/houghton-sofic-profile)
supplies the Houghton far-commutator certificates, reweighting and group-profile
arguments. Its verification notes state the exact formalization boundary.

This repository does **not** formally verify the Haar concentration theorem,
weighted extraction constants, conditional entropy arguments, adaptive process
norms, the complete reservoir compiler, or the cited operator-algebra results.
Closed-device achievability remains a cited external theorem. These analytic
dependencies are proved or cited in the manuscript; automated review and these
checks do not replace their proofs.
