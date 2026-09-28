# Causal quantum-channel simulation: memory beyond entropy

Seth Douglas and Nidhal Mghirbi — September 2026.

[Read the manuscript](paper.pdf) · [TeX source](paper.tex) ·
[Verification scope](verification/README.md)

Two channels with matched Choi rank, maximal complementary entropy, flat optimal
complementary spectrum and asymptotic purity cost can nevertheless have constant
versus polynomial memory requirements at the same optimal exchange rate.
The paper studies this through finite-approximation profiles, sequential support
tests, reservoir constructions and group presentations.

The matched-invariant separation is Theorem A, stated in full at the start.
Both the rank-rate compiler and the entropy-rate closure compiler (in its stated
rolling regime) use zero initial and imported purity under the stated
round-dependent control contract. Error always refers to the complete adaptive
experiment, not merely one-call marginals. Gaps between upper and lower exponents
remain explicit.

For any pure probe with a flat output spectrum, memory bounds total purity at
its rank exchange rate for every fixed error below one. This broader converse
does not change the small-error restriction in the Houghton separation.

For full-Choi-rank targets at fixed positive exchange slack, the compiler reads
the finite-factor profile at accuracy of order epsilon/sqrt(n), with zero
initial and imported purity in its rolling branch. This refinement does not
apply to the proper-support Houghton target without additional hypotheses.

This is a **private release candidate**, not an announced publication.
Licenses and the joint public-release decision remain pending; see [RIGHTS.md](RIGHTS.md).

## Reproduce

Requires Python 3.10+, NumPy and SciPy for numerical checks, a TeX distribution with
the packages listed in paper.tex, and Lean 4.30.0 for the scalar formal checks.

    python -m pip install -r verification/requirements.txt
    python build.py
    python verification/check.py

Checks cover specified finite identities and arithmetic, not the full analytical
proofs. See the verification notes before interpreting a passing result.
Generated TeX files stay in build/; paper.pdf is the canonical output.
The bibliography is embedded in paper.tex.

## Dependencies

The group-theoretic companion,
[Houghton's group H_3 has superpolynomial sofic profile](https://github.com/Apsiape/houghton-sofic-profile/tree/0873417a8b8b5f3a6d0f7dfcf4dd255d37922ed9),
is pinned to version `0873417` in the manuscript. It
contains the group-word certificates and their Lean checker. They are not
duplicated here. The manuscript separately cites the published closed-device
achievability theorem and other external results. Companion links require access
while the repositories remain private.
