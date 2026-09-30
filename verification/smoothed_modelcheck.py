"""Rebuild finite one-particle models and check the gadget holonomies.

Uses manuscript formulas, with a least-squares curl solver instead of the
spanning-tree construction. No saved model matrices or private-repo imports.
Checks M=6,7,8,12, all five relators and 52 literal/canonical words at every clock.
Second-order moment bounds are tested only for even M. Floating-point regression,
not a proof for all M or a construction of the exponentially large Fock matrices.
"""
if not __debug__:
    raise SystemExit("This check uses assert statements; run it without python -O.")
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import math
import numpy as np
from houghton_labels import REL, slots, signature, inverse, act

TOL = 3e-9
WORDS = sorted({word for word, _ in slots})
EXCEPTIONAL = {"abA", "abAB", "abABX"}
assert len(WORDS) == 41


def cycle(L, r):
    # Right-action reverse r-cycle, represented as a usual column matrix.
    image = [r-1] + list(range(r-1)) + list(range(r, L))
    return np.eye(L, dtype=complex)[:, image]


def swap(L, i, j):
    p = np.eye(L, dtype=complex)
    p[[i, j]] = p[[j, i]]
    return p


def face_edges(M, q):
    i, j = q
    return [(('a', i, j), 1), (('b', (i+1) % M, j), 1),
            (('a', i, (j+1) % M), -1), (('b', i, j), -1)]


def regions(M):
    out = []
    for x in range(4, M+3):  # Zero-based t=x+1.
        faces = {(M-1, M-1)} | {(M-1, j) for j in range(1, M-1)}
        faces |= {(i, j) for i in range(M-2) for j in range(M-1) if i+j <= x+M-7}
        def label(i, j, x=x):
            return (x+M-1-j-(i if i < M-1 else 0), x+2*M-1-j)
        out.append(((M-1, M-1), faces, label, (-1)**x))
    faces = {(0, M-2)} | {(i, j) for i in range(M) for j in range(M) if i+j <= M-3}
    out.append(((0, M-2), faces, lambda i, j: (M+2-i-j, 2*M+2-j), (-1)**(M-1)))
    return out


def curl_solution(M, source, faces, sign):
    rows = sorted(faces)
    owners = {}
    for q in rows:
        for edge, direction in face_edges(M, q):
            owners.setdefault(edge, []).append((q, direction))
    columns = sorted(e for e, qs in owners.items() if len(qs) == 2)
    row_index = {q: i for i, q in enumerate(rows)}
    curl = np.zeros((len(rows), len(columns)))
    for j, edge in enumerate(columns):
        for q, direction in owners[edge]:
            curl[row_index[q], j] = direction
    demand = np.full(len(rows), sign*math.pi/len(rows))
    demand[row_index[source]] -= sign*math.pi
    potentials, _, rank, _ = np.linalg.lstsq(curl, demand, rcond=None)
    assert rank == len(rows)-1
    assert np.max(np.abs(curl @ potentials-demand)) < TOL
    return dict(zip(columns, potentials))


def build(M):
    L, ell = 4*M+4, 2*M+2
    eye, tau = np.eye(L, dtype=complex), swap(L, 0, 1)
    top = [cycle(L, 2*ell-M+1)]
    for i in range(M-1):
        top.append(cycle(L, ell-i) @ top[-1] @ tau @ cycle(L, ell-i-(M-1)).conj().T)
    links = {}
    for i in range(M):
        for j in range(M):
            links['a', i, j] = cycle(L, ell-i-j)
            links['b', i, j] = top[i].copy() if j == M-1 else cycle(L, 2*ell-j)
    links['a', M-1, M-1] = links['a', M-1, M-1] @ swap(L, 0, M+3)
    links['b', 0, M-2] = links['b', 0, M-2] @ swap(L, 4, M+4)
    flux = {q: [] for q in ((i, j) for i in range(M) for j in range(M))}
    region_data = regions(M)
    assert sum(s for _, _, _, s in region_data) == M % 2
    for source, faces, label, sign in region_data:
        for q in faces:
            pair = label(*q)
            assert min(pair) >= 4 and max(pair) < L and pair[0] != pair[1]
            flux[q].append((pair, sign*math.pi/len(faces)))
        for edge, value in curl_solution(M, source, faces, sign).items():
            pair = label(edge[1], edge[2])
            v = np.zeros(L)
            v[pair[0]], v[pair[1]] = 1/np.sqrt(2), -1/np.sqrt(2)
            rotation = eye + np.expm1(1j*value)*np.outer(v, v)
            links[edge] = links[edge] @ rotation
    assert all(np.max(np.abs(U.conj().T @ U-eye)) < TOL for U in links.values())
    return links, tau, flux


def holonomy(M, links, tau, word):
    L = tau.shape[0]
    matrices, endpoints = [], []
    for i0 in range(M):
        for j0 in range(M):
            i, j = i0, j0
            H = np.eye(L, dtype=complex)
            for letter in word:
                if letter == 'a':
                    U = links['a', i, j]; i = (i+1) % M
                elif letter == 'A':
                    i = (i-1) % M; U = links['a', i, j].conj().T
                elif letter == 'b':
                    U = links['b', i, j]; j = (j+1) % M
                elif letter == 'B':
                    j = (j-1) % M; U = links['b', i, j].conj().T
                else:
                    U = tau
                H = U @ H
            matrices.append(H)
            endpoints.append((i, j))
    return endpoints, np.array(matrices)


# Reconstruct normal forms by ray actions and BFS in S4; no supplied label table.
permutation_words = {tuple(range(4)): ""}
frontier = [tuple(range(4))]
for p in frontier:
    for j, generator in enumerate(("x", "axA", "aaxAA")):
        transposition = list(range(4))
        transposition[j], transposition[j+1] = j+1, j
        q = tuple(transposition[x] for x in p)
        if q not in permutation_words:
            permutation_words[q] = permutation_words[p] + generator
            frontier.append(q)
representatives, canonical = {}, {}
for word in WORDS:
    u, v = word.count('a')-word.count('A'), word.count('b')-word.count('B')
    translation = ('a'*u if u >= 0 else 'A'*(-u)) + ('b'*v if v >= 0 else 'B'*(-v))
    hword = word + inverse(translation)
    p = tuple(act((1, j+1), hword)[1]-1 for j in range(4))
    assert p in permutation_words
    for ray in (1, 2, 3):
        for height in range(1, len(hword)+3):
            if ray != 1 or height > 4:
                assert act((ray, height), hword) == (ray, height)
    canonical[word] = permutation_words[p] + translation
    representatives.setdefault(signature(word), canonical[word])
assert len(representatives) == 26
ALL_WORDS = sorted(set(WORDS) | set(canonical.values()))
assert len(ALL_WORDS) == 52  # 41 literal words plus 11 independently selected representatives.

for M in (6, 7, 8, 12):
    links, tau, flux = build(M)
    L, eye = len(tau), np.eye(len(tau), dtype=complex)
    data = {word: holonomy(M, links, tau, word) for word in ALL_WORDS}
    start = [(i, j) for i in range(M) for j in range(M)]
    defect = data['abABX'][1]
    for word in REL:
        end, matrices = data[word]
        assert end == start
        if word != 'abABX':
            assert np.max(np.abs(matrices-eye)) < TOL, (M, word)
    predicted_det = []
    for at, q in enumerate(start):
        predicted = eye.copy()
        used = set()
        for pair, angle in flux[q]:
            assert used.isdisjoint(pair)
            used.update(pair)
            v = np.zeros(L)
            v[pair[0]], v[pair[1]] = 1/np.sqrt(2), -1/np.sqrt(2)
            predicted += np.expm1(1j*angle)*np.outer(v, v)
        assert np.max(np.abs(defect[at]-predicted)) < TOL
        determinant = np.prod([(1+np.exp(1j*f))/2 for _, f in flux[q]])
        assert abs(np.linalg.det((eye+defect[at])/2)-determinant) < TOL
        predicted_det.append(determinant)
    for word in WORDS:
        rep = representatives[signature(word)]
        assert data[word][0] == data[rep][0]
        expected = data[rep][1] @ defect if word in EXCEPTIONAL else data[rep][1]
        assert np.max(np.abs(data[word][1]-expected)) < TOL, (M, word, rep)
    # All same-clock canonical comparisons act on the four gadget coordinates.
    reps = list(representatives.values())
    for i, a in enumerate(reps):
        for b in reps[i:]:
            if data[a][0] != data[b][0]:
                continue
            W = data[a][1].conj().transpose(0, 2, 1) @ data[b][1]
            assert np.max(np.abs(W[:, 4:, 4:]-np.eye(L-4))) < TOL
            assert np.max(np.abs(W[:, :4, 4:])) < TOL
            assert np.max(np.abs(W[:, 4:, :4])) < TOL
            assert np.max(np.abs(W[:, :4, :4]-W[0, :4, :4])) < TOL
    for k in (1, 2, 3):
        z = np.mean(np.array(predicted_det)**k)
        if M % 2 == 0:
            minimum = (M-1)*(M-2)/2+1
            bound = k*math.pi**2/(8*M*minimum) + k*k*math.pi**2/(2*minimum**2)
            assert abs(1-z) <= bound+TOL
    print(f"Smoothed matrices M={M}: five relators, 52 word holonomies, canonical coordinate support and determinant moments pass.")
    del data

print("Second-order moment checks restricted to even M; all-M geometric proofs remain manuscript dependencies.")
