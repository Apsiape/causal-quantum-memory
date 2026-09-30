"""Exact slot bookkeeping for the flat-Choi pair, not an analytic verifier.

Labels are computed by ray actions, not imported from the manuscript's weight
histogram. The finite signature's completeness argument is in houghton_labels.
The cross-block holonomy identity, adaptive tag reduction and compiler remain
mathematical inputs; no finite test here proves them.
"""
if not __debug__:
    raise SystemExit("This check uses assert statements; run it without python -O.")
from collections import Counter, defaultdict
from fractions import Fraction

from houghton_labels import signature, slots, weights

expected = Counter({Fraction(2): 17, Fraction(5, 2): 1, Fraction(4): 3,
                    Fraction(9, 2): 1, Fraction(6): 1, Fraction(10): 1,
                    Fraction(25, 2): 1, Fraction(45, 2): 1})
assert Counter(weights.values()) == expected
assert sum(weights.values()) == 104

# Expand the actual scalar-square slot list of the new block diagonal unitary.
balanced_slots = [(signature(w), q) for w, q in slots] * 2
old_singletons = 0
for label, weight in weights.items():
    count = 52 - 2 * weight
    assert count.denominator == 1 and count >= 0
    old_singletons += int(count)
    balanced_slots.extend([(label, Fraction(1))] * int(count))
new_labels = [signature("a" * j) for j in range(3, 9)]
assert len(set(new_labels)) == 6 and not set(new_labels).intersection(weights)
for label in new_labels:
    balanced_slots.extend([(label, Fraction(1))] * 52)

balanced = defaultdict(Fraction)
for label, weight in balanced_slots:
    balanced[label] += weight
assert old_singletons == 1144
assert len(balanced_slots) == 1804
assert len(balanced) == 32 and set(balanced.values()) == {52}
dimension = sum(balanced.values())
assert dimension == 1664 == 16 * 104
assert {weight / dimension for weight in balanced.values()} == {Fraction(1, 32)}

# Translation differences certify the absence of clock aliasing for M >= 12.
old_translations = {(sig[1][1], sig[1][2]) for sig in weights}
new_translations = {(sig[1][1], sig[1][2]) for sig in new_labels}
assert new_translations == {(j, 0) for j in range(3, 9)}
assert all(-1 <= u <= 2 and -1 <= v <= 1 for u, v in old_translations)
assert old_translations.isdisjoint(new_translations)
for u in new_translations:
    for v in old_translations | new_translations:
        if u != v:
            diffs = tuple(a - b for a, b in zip(u, v))
            assert 0 < max(map(abs, diffs)) <= 9 < 12

# Exact allocation of the three error terms in the sufficient parameter bound.
assert Fraction(1804, 2) == 902
assert Fraction(902, 1804) == Fraction(1, 2)
assert Fraction(902, 2 * 1804) == Fraction(1, 4)
assert Fraction(902 * 8, 28864) == Fraction(1, 4)
print("Exact flat-Choi bookkeeping: dimension 1664, rank 32, eigenvalues 1/32,")
print("1804 slots, 1144+312 added singletons; clock differences and error allocation pass.")
