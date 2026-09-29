# Verification scope

Run from the repository root:

    python verification/check.py

Python 3.10+ with NumPy and SciPy, and Lean 4.30.0, are required. The scripts are
bounded finite checks, with per-command timeouts in the runner.

## Included checks

- **houghton_labels.py:** exact integer/rational bookkeeping of the 104-dimensional
  gadget: 33 symbols, 35 triangles, 174 entries, 26 distinct group labels,
  singleton anchors and total squared coefficient weight 104. A length-L ray
  word is represented by its action up to height L+1 and its eventual translations,
  sampled at height L+2. A regression covers the last-letter-alpha boundary case.
  The reason this representation is complete is mathematical, not a Lean theorem.
- **tracecheck.py:** the trace channels of the Houghton face. It rebuilds the
  104-dimensional gadget and its 26 labels and computes the moment matrices of the
  regular trace and of the fermionic traces with one to four flavours. In exact
  rational arithmetic it checks rank 26 and the moment f = tau(gamma) = 4^-k. In
  floating point it checks positive semidefiniteness, the full moment spectrum
  {1 (x12), 1 - 4^-k (x5), 1 + 4^-k, 1 + 4^-k (1 +- sqrt 17)/2 (x4 each)}, the
  inequalities f/104 <= ||J - J_H||_op <= d_J, and the normalized-Choi distance to
  Phi_H (0.0944 for one flavour). It does not prove the face characterization, the
  case f = 1 or the expense bound, which are proved in the manuscript.
- **gadgetcheck.py:** numerical unitary/Kraus/entropy identities for the S3, Z4,
  and Z2 x Z2 finite-group gadgets, and their general spectrum completions:
  the normalized Kraus Gram matrix, trace preservation, unitality, invariant
  corner and clock-label orthogonality. The small Kraus Gram matrix avoids
  allocating a dense Choi matrix. These examples do not prove the general
  group-gadget or spectrum-completion theorems, the approximation-profile
  transfer, or the adaptive tag-recycling argument.
- **flat_choi_check.py:** exact ray-action label and rational slot bookkeeping
  for the dimension-1664 construction: 32 equal coefficient weights, Choi
  eigenvalues 1/32, 1804 slots, translation differences and error allocation.
  This does not prove the cross-block holonomy identity, adaptive tag reduction,
  or reservoir compiler; those arguments are given in the manuscript.
- **tagcheck.py:** numerical verification of the tagging construction and its
  one-half Choi-distance identity in small stand-in dimensions. It does not
  directly construct the 208-dimensional hard channel's infinite factor.
- **cornercheck.py:** deterministic finite matrix tests of corner compression,
  the two-unitary polar repair and its single doubled-bath realization, including
  singular compressions and mixtures with cancelling perturbations. It checks
  trace preservation, unitality, leakage/completion trace equality and the
  dimension-scaled Choi-distance bound. It does not prove the general profile,
  closure, purity-monotonicity or adaptive-transfer theorems.
- **accuracycheck.py:** Choi distance, support leakage and product-state distances
  for the rotating-dephasing example. The all-observer adaptive equality follows
  from the manuscript's common-processor proof, not finite numerical sampling.
  It also checks the flat-Choi SWAP example's rank-four support, distance t/2,
  and quadratic entropy-production expansion at three parameters. The universal
  obstruction to linear trace-distance repair is the written unitality argument.
- **scalarcheck.py:** exact rational upper bounds for the smoothed Choi example
  at k=16, M=2120 and error 10^-9, and the exact integer bracket that forces the
  bath-qubit ceiling 135767. The analytic Choi inequalities and pi < 22/7 are inputs.
  Integer comparisons certify the 107-exponent margins c_H > .0016 and
  c_a > .00017; floating-point tests check only the illustrative rounded displays.
- **reservoir_dimension.py:** a sufficient active-bath dimension selector with
  integer/rational acceptance tests for Appendix A's dimension and padding
  inequalities. Supply an integer upper bound on the analytic buffer and a
  maximum export size. This chooses a dimension, not the Haar-existence mixers.
  For example:

      python verification/reservoir_dimension.py --m 3 --r 4 --n 2 --export-qubits 2 --buffer-bits 3 --epsilon 1/10

- **adaptive_framecheck.py:** a physical two-round SWAP-bath experiment with
  GHZ-entangled inputs/reference and coherent output-to-next-input feedback,
  compared independently with the history-Gram expansion. Exact and perturbed
  mixers check the Gram increments for survivors and first-exit children,
  deletion bounds, and a negative control for the Gram transpose. Artificial
  pruning tests the linear algebra; it does not validate the reflected-cost
  tail, Haar existence, or all adaptive observers.
- **smoothed_modelcheck.py:** independently reconstructs the one-particle link
  matrices from the manuscript, using a least-squares curl solver rather than
  the spanning-tree construction. At M=6,7,8,12 it checks all five relators and
  52 word holonomies (41 literal gadget words plus 11 canonical representatives
  regenerated by exact ray actions and BFS), canonical four-coordinate support,
  and determinant moments. Second-order bounds are tested only for even M.
  No private code or saved model matrices are loaded. This is a finite numerical
  regression, not an all-M proof or an enumeration of the Fock-space matrices.
- **Lean Resources.lean:** exact dimension arithmetic, telescoping integer
  replacement schedules, a conditional integer budget inequality, and the cleared-
  denominator scalar margins in the weighted-extraction ledger. The integer
  schedule applies to dyadic rank; the real-logarithmic general schedule and
  entropy budgets are not formalized. It also checks the integer-scaled budget
  rearrangement for the fixed-error extension, conditional on its analytic
  premises, the numerical example's bath-dimension bracket, omission of the final
  replacement in the integer schedule, and the four-tag-qubit budget rearrangement.

The Lean theorems contain no proof placeholders, native proof-evaluation oracle,
or added mathematical axioms. Their standard logical dependencies are printed
when run.

## External and unformalized dependencies

The companion [group paper, version 1.1.0](https://github.com/Apsiape/houghton-sofic-profile/tree/v1.1.0)
supplies the Houghton far-commutator certificates, reweighting and group-profile
arguments. Its verification notes state the exact formalization boundary.

This repository does **not** formally verify the Haar concentration theorem,
weighted matrix-norm and entropy inequalities (only their listed scalar arithmetic is kernel-checked), conditional entropy arguments, adaptive process
norms, the complete reservoir compiler, or the cited operator-algebra results.
Closed-device achievability remains a cited external theorem. These analytic
dependencies are proved or cited in the manuscript; automated review and these
checks do not replace their proofs.
