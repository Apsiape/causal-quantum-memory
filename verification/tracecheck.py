"""Trace channels of the 104-dimensional Houghton gadget (Section "The Houghton face").

Reproduces the gadget exactly (33 symbols, 35 triangles, d = 104, right actions) and its 26 labels. For the regular
trace and the fermionic traces tau_F^k (k = 1, 2, 3, 4 flavours) it computes the 26 x 26 moment matrix
T(g,h) = tau(g h^-1) in exact rational arithmetic, its exact rank, the moment f = tau(gamma), its least eigenvalue,
and the normalized-Choi distance to Phi_H from J(Phi_tau) = V D_w^(1/2) T D_w^(1/2) V^* with w_g = ||K_g||_F^2 / 104.
For finitary g, tau_F(g) = det((I + P_g)/2), the product over cycles of length l of (1 - (-1)^l)/2^l; tau_F = 0 off
the finitary subgroup. Exits 1 if an exact value differs from the manuscript.
"""
from fractions import Fraction
import sys
import numpy as np

W, FAR = 40, 90


def act(letter, p):
    i, j = p
    if letter in "xX":
        return (1, 2) if p == (1, 1) else (1, 1) if p == (1, 2) else p
    ray = 2 if letter in "aA" else 3
    if letter.islower():
        if i == 1:
            return (ray, 1) if j == 1 else (1, j - 1)
        return (ray, j + 1) if i == ray else p
    if i == ray:
        return (1, 1) if j == 1 else (ray, j - 1)
    return (1, j + 1) if i == 1 else p


def image(word, p):
    for letter in word:
        p = act(letter, p)
    return p


POINTS = [(i, j) for i in (1, 2, 3) for j in range(1, W + 1)]


def key(w):
    return tuple(image(w, p) for p in POINTS)


def inv(w):
    return "".join(c.swapcase() for c in reversed(w))


def cycles(word):
    """Cycle lengths of a finitary permutation, or None if the word has nonzero translation."""
    for r in (1, 2, 3):
        if image(word, (r, FAR)) != (r, FAR):
            return None
    seen, out = set(), []
    for p in POINTS:
        if p in seen:
            continue
        q, n = p, 0
        while True:
            seen.add(q)
            q = image(word, q)
            n += 1
            if q == p:
                break
        out.append(n)
    return out


FAILS = []
rel = {1: "xx", 2: "xaxA" * 3, 3: "x" + "aaxAA" + "X" + "aaXAA", 4: "abABX", 5: "axAbXB"}
e_key = key("")
if not all(key(w) == e_key for w in rel.values()):
    FAILS.append("relators")

symbols = [(s, s) for s in "aAbBxX"] + [(f"p{i},{l}", rel[i][:l]) for i in range(1, 6) for l in range(2, len(rel[i]))]
word = dict(symbols)
triangles = [((rel[i][0] if k == 2 else f"p{i},{k-1}"), rel[i][k - 1]) for i in range(1, 6) for k in range(2, len(rel[i]) + 1)]
triangles += [("a", "A"), ("b", "B"), ("x", "X")]
d = 1 + len(symbols) + 2 * len(triangles)
if (len(symbols), len(triangles), d) != (33, 35, 104):
    FAILS.append("gadget size")

s2 = 2 ** -0.5
entries = [(1.0, "")] + [(1.0, w) for _, w in symbols]
for x, y in triangles:
    entries += [(s2, ""), (s2, word[y]), (s2, word[x]), (-s2, word[x] + word[y])]
elements, weight = {}, {}
for c, w in entries:
    k = key(w)
    elements.setdefault(k, w)
    weight[k] = weight.get(k, 0.0) + c * c
keys = list(elements)
print(f"gadget: {len(symbols)} symbols, {len(triangles)} triangles, d = {d}, {len(keys)} labels, "
      f"total weight {sum(weight.values()):.6f}")
if len(keys) != 26 or abs(sum(weight.values()) - 104) > 1e-9:
    FAILS.append("labels")

pair_cycles = [[cycles(elements[g] + inv(elements[h])) for h in keys] for g in keys]


def moments(char):
    return [[Fraction(0) if cy is None else char(cy) for cy in row] for row in pair_cycles]


def fermion(k):
    def c(cy):
        out = Fraction(1)
        for l in cy:
            out *= Fraction(1 - (-1) ** l, 2 ** l)
        return out ** k
    return c


def regular(cy):
    return Fraction(1) if all(l == 1 for l in cy) else Fraction(0)


def rank(T):
    M = [row[:] for row in T]
    rk, n = 0, len(M)
    for col in range(n):
        piv = next((i for i in range(rk, n) if M[i][col] != 0), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(n):
            if i != rk and M[i][col] != 0:
                f = M[i][col] / M[rk][col]
                M[i] = [u - f * v for u, v in zip(M[i], M[rk])]
        rk += 1
    return rk


wv = np.array([weight[k] for k in keys]) / d
Nh = np.diag(np.sqrt(wv))


def choi(T):
    return Nh @ np.array([[float(x) for x in row] for row in T]) @ Nh


J_reg = choi(moments(regular))
ig, ie = keys.index(key("xaxA")), keys.index(e_key)

print("\ntrace          rank  tau(gamma)  least eigenvalue  d_J to Phi_H")
expected = {"regular": Fraction(0), "fermion k=1": Fraction(1, 4), "fermion k=2": Fraction(1, 16),
            "fermion k=3": Fraction(1, 64), "fermion k=4": Fraction(1, 256)}
for name, ch in [("regular", regular)] + [(f"fermion k={k}", fermion(k)) for k in (1, 2, 3, 4)]:
    T = moments(ch)
    rk, f = rank(T), T[ig][ie]
    mn = float(np.linalg.eigvalsh(np.array([[float(x) for x in row] for row in T])).min())
    dist = 0.5 * float(np.abs(np.linalg.eigvalsh(choi(T) - J_reg)).sum())
    print(f"{name:14s} {rk:4d}  {str(f):>10s}  {mn:16.4f}  {dist:.4f}")
    if rk != 26 or f != expected[name] or mn < -1e-9:
        FAILS.append(name)
    opn = float(np.abs(np.linalg.eigvalsh(choi(T) - J_reg)).max())
    if not (float(f) / 104 <= opn + 1e-12 and opn <= dist + 1e-12):
        FAILS.append(f"f/104 <= ||J - J_reg||_op <= d_J fails for {name}")
    if name.startswith("fermion"):
        k = int(name[-1])
        x = 4.0 ** -k
        want = sorted([1.0] * 12 + [1 - x] * 5 + [1 + x] + [1 + x * (1 + 17 ** 0.5) / 2] * 4 + [1 + x * (1 - 17 ** 0.5) / 2] * 4)
        got = sorted(np.linalg.eigvalsh(np.array([[float(y) for y in row] for row in T])))
        if max(abs(a - b) for a, b in zip(got, want)) > 1e-9:
            FAILS.append(f"moment spectrum of {name}")
    if name == "fermion k=1" and not 0.0935 < dist < 0.0955:
        FAILS.append("distance of the one-flavour channel")

print("\nall checks passed" if not FAILS else f"\nFAILED: {FAILS}")
sys.exit(1 if FAILS else 0)
