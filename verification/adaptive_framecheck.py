"""Two-round, reference-entangled coherent-feedback regression.

Independently compare a physical purified bath evolution with its history Gram
expansion, then check survivors AND first-exit children. Artificial pruning tests
the linear algebra, not the manuscript's reflected-cost tail or Haar existence.
Finite floating-point tests do not prove the all-observer compiler theorem.
"""
import numpy as np
from scipy.linalg import expm

TOL = 2e-11
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1, -1]).astype(complex)
PAULI = (I, X, Y, Z)
SWAP = np.eye(4)[[0, 2, 1, 3]]
CNOT = np.eye(4)[[0, 1, 3, 2]]
phases = np.exp(1j*np.array([0.11, 0.37, -0.23, 0.53]))
B = [phase*np.kron(p, I) for p, phase in zip(PAULI, phases)]
K = [p/(2*phase) for p, phase in zip(PAULI, phases)]


def apply(state, op, axes, qubits):
    axes = list(axes)
    order = axes + [i for i in range(qubits) if i not in axes]
    data = state.reshape([2]*qubits).transpose(order).reshape(2**len(axes), -1)
    data = (op @ data).reshape([2]*qubits).transpose(np.argsort(order))
    return data.reshape(-1)


def gram(frame):
    return frame.conj().T @ frame


def opnorm(x):
    return np.linalg.norm(x, 2)


rng = np.random.default_rng(91374)
feedback = np.kron(expm(0.31j*Y), expm(0.17j*X)) @ CNOT
psi = np.zeros(8, complex)  # Two inputs and a reference, initially GHZ-entangled.
psi[[0, 7]] = 1/np.sqrt(2)
flat_purifier = np.eye(4).reshape(-1)/2
epsilon, n = 0.1, 2
for strength in (0.0, 1e-5, 1e-4):
    raw = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))
    H = (raw + raw.conj().T)/2
    V = expm(1j*strength*H) @ SWAP
    # Physical collision: SWAP with the first bath qubit in each round.
    physical = np.kron(psi, flat_purifier)
    physical = apply(physical, SWAP, [0, 3], 7)
    physical = apply(physical, feedback, [0, 1], 7)
    physical = apply(physical, V, [3, 4], 7)
    physical = apply(physical, SWAP, [1, 3], 7)
    rho_physical = physical.reshape(8, 16) @ physical.reshape(8, 16).conj().T

    first = [apply(psi, k, [0], 3) for k in K]
    frames1 = np.column_stack([b.reshape(-1)/2 for b in B])
    histories = [(i, j) for i in range(4) for j in range(4)]
    coefficients, bath = [], []
    for i, j in histories:
        coefficients.append(apply(apply(first[i], feedback, [0, 1], 3), K[j], [1], 3))
        bath.append((B[j] @ V @ B[i]).reshape(-1)/2)
    T, F = np.column_stack(coefficients), np.column_stack(bath)
    rho_gram = T @ gram(F).T @ T.conj().T
    assert np.max(np.abs(rho_physical-rho_gram)) < TOL
    assert abs(np.trace(T @ T.conj().T)-1) < TOL
    assert np.max(np.abs(physical.reshape(8, 16) - T @ F.T)) < TOL

    # First exit at time 1: i=3. At time 2: j=2 with surviving i!=3.
    children2 = [h for h, (i, _) in enumerate(histories) if i != 3]
    good = [h for h, (i, j) in enumerate(histories) if i != 3 and j != 2]
    errors = []
    for time, frame in ((1, frames1), (1, frames1[:, :3]),
                        (2, F[:, children2]), (2, F[:, good])):
        error = opnorm(gram(frame)-np.eye(frame.shape[1]))
        assert error <= time*epsilon/(8*n) + TOL
        errors.append(error)
    ideal_exit = sum(np.vdot(coefficients[h], coefficients[h]).real
                     for h in range(16) if h not in good)
    deleted1 = [h for h, (i, _) in enumerate(histories) if i == 3]
    deleted2 = [h for h, (i, j) in enumerate(histories) if i != 3 and j == 2]
    p1 = np.vdot(first[3], first[3]).real
    p2 = sum(np.vdot(coefficients[h], coefficients[h]).real for h in deleted2)
    assert abs(p1+p2-ideal_exit) < TOL
    norms = []
    for subset, p in ((deleted1, p1), (deleted2, p2)):
        removed = T[:, subset] @ F[:, subset].T
        norms.append(np.linalg.norm(removed))
        assert norms[-1] <= np.sqrt((1+epsilon/8)*p) + TOL
    good_vector = T[:, good] @ F[:, good].T
    assert np.linalg.norm(physical.reshape(8, 16)-good_vector) <= sum(norms)+TOL
    assert sum(norms) <= np.sqrt((1+epsilon/8)*n*ideal_exit)+TOL
    # A transposition error is observable in the perturbed complex example.
    if strength == 1e-4:
        wrong = T @ gram(F) @ T.conj().T
        wrong_error = np.max(np.abs(wrong-rho_physical))
        assert wrong_error > TOL, wrong_error
    print(f"Entangled feedback strength={strength:g}: max frame error={max(errors):.3e}; physical/Gram and first-exit checks pass.")

print("Finite two-round checks only; no assertion of the general compiler or its tail bound.")
