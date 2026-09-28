"""Exact finite bookkeeping of all entry labels of the 104-dimensional gadget.

A length-L word acts by a fixed translation on each ray above L+1. Before
each letter, a start at height L+2 is at height at least 3: translations move
by at most one, and alpha only acts at heights 1 and 2. Thus values at heights
1..L+1 and the three eventual translations determine the action exactly.
The proof of that observation and the entropy/anchor implications remain in the
manuscript; this program verifies the integer/rational bookkeeping.
"""
from fractions import Fraction
from collections import defaultdict

INV = dict(zip("aAbBxX", "AaBbXx"))
def inverse(w): return "".join(INV[c] for c in reversed(w))
def comm(u, v): return u + v + inverse(u) + inverse(v)
REL = ["xx", "xaxA" * 3, comm("x", "aaxAA"), "abABX", "axAbXB"]

def step(point, s):
    ray, j = point
    if s in "xX":
        return (1, 3-j) if ray == 1 and j <= 2 else point
    destination = 2 if s in "aA" else 3
    if s in "ab":
        if ray == 1: return (destination, 1) if j == 1 else (1, j-1)
        if ray == destination: return ray, j+1
    else:
        if ray == destination: return (1, 1) if j == 1 else (ray, j-1)
        if ray == 1: return 1, j+1
    return point

def act(point, word):
    for s in word: point = step(point, s)
    return point

symbols = list("aAbBxX")
for w in REL:
    symbols += [w[:i] for i in range(2, len(w))]
triangles = []
for w in REL:
    for i in range(2, len(w)+1):
        triangles.append((w[:i-1], w[i-1], w[:i] if i < len(w) else ""))
triangles += [(s, INV[s], "") for s in "abx"]
slots = [(w, Fraction(1)) for w in [""] + symbols]
for x, y, z in triangles:
    slots += [(w, Fraction(1, 2)) for w in ["", y, x, x+y]]
L = max(map(len, [w for w, _ in slots] + REL))

def signature(w):
    assert len(w) <= L
    base = tuple(act((ray, j), w) for ray in (1,2,3) for j in range(1, L+2))
    tail = tuple(act((ray, L+2), w)[1]-(L+2) for ray in (1,2,3))
    assert all(act((ray, L+2), w)[0] == ray for ray in (1,2,3))
    return base, tail

# The older height-L cutoff is not valid for general words: the final alpha
# can still act at height 2. This regression protects the extra boundary layer.
for length in range(1, 25):
    word = "a" * (length - 1) + "x"
    assert act((1, length+1), word) == (1, 1)
    assert act((1, length+2), word) == (1, 3)
    assert 3 - (length+2) == -(length-1)

identity = signature("")
assert all(signature(w) == identity for w in REL)
assert all(signature(x+y) == signature(z) for x,y,z in triangles)
assert len(symbols) == 33 and len(triangles) == 35 and len(slots) == 174
anchors = {signature(w) for w in [""] + symbols}
weights = defaultdict(Fraction)
for w, q in slots: weights[signature(w)] += q
assert set(weights) == anchors and len(weights) == 26
assert sum(weights.values()) == 104
print("Exact: 33 symbols, 35 triangles, 174 entries, 26 distinct labels; all labels anchored.")
print("Exact: coefficient-square sum 104; normalized Choi eigenvalues sum to 1.")
print("Squared coefficient weights:", sorted(weights.values()))
