# Causal quantum-channel simulation: memory beyond entropy

Seth Douglas and Nidhal Mghirbi — September 2026.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23027626.svg)](https://doi.org/10.5281/zenodo.23027626)

[Read the manuscript](paper.pdf) · [TeX source](paper.tex) ·
[Verification scope](verification/README.md)

Two explicit channels on dimension 1664 have the same complete normalized Choi
spectrum (32 eigenvalues 1/32), maximal complementary entropy five and zero
asymptotic purity cost. At the same optimal exchange rate, one uses five memory
qubits exactly; the other needs polynomially many, even at fixed error at most
1/16 and without a separate purity restriction. The latter also has a
zero-purity service using O(n^(1/3) log^(4/3) n) memory at fixed error. The lower
and upper exponents do not match. Theorem A states the precise separation.

Version 1.1 locates the expensive channels (Theorem B): they lie on the relative
boundary of the closure of finite tracial factorizations; any fixed depolarizing
admixture makes them cheap at the rank rate; on a face containing the Houghton
channel a single moment of a trace on the group decides expense; and one free
fermion per site already forces polynomial memory.

Every output is released before the next input arrives. Error refers to the
complete adaptive experiment with quantum references, not merely one-call
marginals. All retained registers count as memory, and exchange counts
replacement slots rather than both directions of traffic.

This repository contains the preprint, its source, and scoped verification
artifacts. The manuscript is the authority for exact statements and hypotheses;
the automated checks do not constitute full formal verification or peer review.

## License and citation

The manuscript and repository content are available under **CC BY 4.0**.
The software and machine-readable certificates are additionally available under
the **MIT License**, at your option; see [RIGHTS.md](RIGHTS.md) for the scope.

Please cite Seth Douglas and Nidhal Mghirbi, *Causal quantum-channel simulation:
memory beyond entropy* (2026), version 1.1.0,
[doi:10.5281/zenodo.23027626](https://doi.org/10.5281/zenodo.23027626) (all versions).
[CITATION.cff](CITATION.cff) provides
machine-readable citation metadata. Tagged releases archive the corresponding
manuscript and verification artifacts together.

The DOI above is the all-versions DOI, which identifies the evolving work. The v1.0.0 archive is
[doi:10.5281/zenodo.23027627](https://doi.org/10.5281/zenodo.23027627).

## Reading guide

The main article occupies pages 1–38, including the abstract and contents.
The complete 74-page manuscript includes all technical appendices and references;
no separate unpublished proof supplement is needed.

| Question | Where to read |
|---|---|
| What is separated, and under which resource contract? | Theorem A; Sections 1–2 |
| What sets the exchange and purity boundary? | Section 3 |
| How is the hard channel constructed, bounded below and served? | Sections 4–6 |
| Why do the complete Choi spectra agree? | Section 7 |
| What extends to exact factors and other group channels? | Sections 8–9 |
| Where do the expensive channels lie? | Theorem B; Section 10 |
| What remains unresolved? | Section 11 |

Appendices A–D contain the compiler, weighted extraction and refined converse
proofs. Appendices E–K retain alternative services, purity examples, the general
group correspondence, the Slofstra application, the smaller matched-invariant
pair, general spectrum completion and adaptive-accuracy comparisons.

The general exchange–purity rate region uses the cited closed-device
achievability theorem; the explicit separation does not. Companion group results
are version-pinned below. The manuscript distinguishes these dependencies from
its own proofs and from the narrower scope of the executable checks.

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
[Houghton's group H_3 has superpolynomial sofic profile](https://github.com/Apsiape/houghton-sofic-profile/tree/v1.1.0),
is cited by the statement numbers of its version 1.1.0 in the manuscript. It
contains the group-word certificates and their Lean checker. They are not
duplicated here. The manuscript separately cites the closed-device preprint's
achievability theorem and other external results. The pinned group source remains
the mathematical dependency even when its release metadata is updated.
The companion's all-versions DOI is
[doi:10.5281/zenodo.23027624](https://doi.org/10.5281/zenodo.23027624).
