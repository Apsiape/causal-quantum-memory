"""Finite numerical checks of invariant-corner compression and polar repair.

Tests single factors and mixtures (including cancelling perturbations and
singular corner contractions). Not a proof of profile or adaptive theorems.
"""
import numpy as np
from scipy.linalg import block_diag, expm

rng = np.random.default_rng(608)


def unitary(n):
    z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / np.abs(np.diag(r)))


def choi(u, d, m):
    blocks = u.reshape(d, m, d, m)
    vectors = blocks.transpose(1, 3, 0, 2).reshape(m*m, d*d).T
    return vectors @ vectors.conj().T / (d*m)


def trace_norm(x):
    return np.abs(np.linalg.eigvalsh((x + x.conj().T) / 2)).sum()


def repair(u, d, m):
    a = u[:d*m, :d*m]
    left, singular, right = np.linalg.svd(a)
    assert singular.max() <= 1 + 1e-12
    c = (left * np.sqrt(np.maximum(0, 1 - singular**2))) @ right
    plus, minus = a + 1j*c, a - 1j*c
    for w in (plus, minus):
        assert np.allclose(w.conj().T @ w, np.eye(d*m), atol=1e-11)
    repaired = (choi(plus, d, m) + choi(minus, d, m)) / 2
    assert np.allclose(repaired, choi(a, d, m) + choi(c, d, m), atol=1e-11)
    # Same channel from one flat 2m-dimensional bath; keep system index first.
    doubled = np.zeros((d, 2*m, d, 2*m), complex)
    doubled[:, :m, :, :m] = plus.reshape(d, m, d, m)
    doubled[:, m:, :, m:] = minus.reshape(d, m, d, m)
    assert np.allclose(choi(doubled.reshape(2*d*m, 2*d*m), d, 2*m), repaired, atol=1e-11)
    return repaired, choi(a, d, m), choi(c, d, m)


count = 0
for D, d, m in ((3, 2, 1), (4, 2, 2), (5, 3, 2), (4, 2, 3)):
    target_u = block_diag(unitary(d*m), unitary((D-d)*m))
    target = choi(target_u, D, m)
    target_corner = choi(target_u[:d*m, :d*m], d, m)
    x = rng.normal(size=(D*m, D*m)) + 1j*rng.normal(size=(D*m, D*m))
    h = (x + x.conj().T) / (2*np.linalg.norm(x))
    p = expm(0.2j*h) @ target_u
    n = expm(-0.2j*h) @ target_u
    cases = [[(1.0, target_u)], [(1.0, p)], [(0.5, p), (0.5, n)],
             [(0.3, unitary(D*m)), (0.7, unitary(D*m))]]
    if D == 2*d:
        swap = np.kron(np.array([[0, 1], [1, 0]]), np.eye(d*m))
        cases.append([(1.0, swap)])  # zero corner contraction
    for components in cases:
        full = sum(w * choi(u, D, m) for w, u in components)
        repaired = np.zeros((d*d, d*d), complex)
        compressed = repaired.copy()
        added = repaired.copy()
        for weight, u in components:
            jr, ja, jc = repair(u, d, m)
            repaired += weight * jr
            compressed += weight * ja
            added += weight * jc
        indices = [out*D + ref for out in range(D) for ref in range(d)]
        chi = (D/d) * full[np.ix_(indices, indices)]
        t = chi[:d*d, :d*d]
        leakage = np.trace(chi[d*d:, d*d:]).real
        assert np.allclose(t, compressed, atol=1e-11)
        assert abs(np.trace(added).real - leakage) < 1e-11
        assert np.allclose(np.trace(repaired.reshape(d, d, d, d), axis1=0, axis2=2), np.eye(d)/d, atol=1e-11)
        assert np.allclose(np.trace(repaired.reshape(d, d, d, d), axis1=1, axis2=3), np.eye(d)/d, atol=1e-11)
        embedded_target = np.zeros_like(chi, dtype=complex)
        embedded_target[:d*d, :d*d] = target_corner
        error = trace_norm(repaired - target_corner) / 2
        input_error = trace_norm(chi - embedded_target) / 2
        full_error = trace_norm(full - target) / 2
        assert error <= input_error + 1e-10
        assert input_error <= (D/d)*full_error + 1e-10
        count += 1

print(f"Invariant-corner polar repair: {count} deterministic finite matrix cases passed.")
print("Numerical regression only; general mixture, closure, purity and adaptive implications are manuscript proofs.")
