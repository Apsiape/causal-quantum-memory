# Numerical check of the matched-invariant tag construction, with small stand-in dimensions.
if not __debug__:
    raise SystemExit("This check uses assert statements; run it without python -O.")
import numpy as np
from scipy.stats import unitary_group
rng = np.random.default_rng(1)

def choi(chan, d):
    # normalized Choi state (chan (x) id)(|Omega><Omega|)
    J = np.zeros((d*d, d*d), complex)
    for i in range(d):
        for j in range(d):
            E = np.zeros((d, d), complex); E[i, j] = 1
            J += np.kron(chan(E), E)
    return J / d

def tn(X): return np.abs(np.linalg.eigvalsh((X + X.conj().T)/2)).sum()

def factor_channel(U, d, m):
    # X -> (id (x) tau_m)[U (X (x) I_m) U^*]
    def ch(X):
        Y = U @ np.kron(X, np.eye(m)) @ U.conj().T / m
        return np.einsum('iaja->ij', Y.reshape(d, m, d, m))
    return ch

D, m, k = 6, 3, 2          # system dim (stand-in for 104), bath of Psi, number of Z powers (stand-in for 6)
ZD = np.diag(np.exp(2j*np.pi*np.arange(D)/D))
Theta = lambda X: sum(np.linalg.matrix_power(ZD, a) @ X @ np.linalg.matrix_power(ZD, a).conj().T for a in range(k)) / k
# target stand-in Phi_H: another exact factor with bath 2
U0 = unitary_group.rvs(D*2, random_state=3); PhiH = factor_channel(U0, D, 2)
U = unitary_group.rvs(D*m, random_state=4); Psi = factor_channel(U, D, m)

def hat(ch):
    def f(X):
        X = X.reshape(2, D, 2, D)
        out = np.zeros((2, D, 2, D), complex)
        out[0, :, 0, :] = ch(X[0, :, 0, :]); out[1, :, 1, :] = Theta(X[1, :, 1, :])
        return out.reshape(2*D, 2*D)
    return f

# Build W on tag(2) (x) sys(D) (x) bathM(m) (x) bathK(k) (x) bathQ(2), per the proof
I = np.eye
P0 = np.diag([1, 0]); P1 = np.diag([0, 1])
V = sum(np.kron(np.linalg.matrix_power(ZD, a), np.outer(I(k)[a], I(k)[a])) for a in range(k))  # sys (x) bathK
# U acts on sys (x) bathM ; embed into sys (x) bathM (x) bathK (x) bathQ
UU = np.kron(np.kron(U, I(k)), I(2))
# V acts on sys (x) bathK ; need sys (x) bathM (x) bathK: permute
Vfull = np.zeros((D*m*k, D*m*k), complex)
Vr = V.reshape(D, k, D, k)
for s in range(D):
    for b in range(k):
        for s2 in range(D):
            for b2 in range(k):
                for mm in range(m):
                    Vfull[(s*m+mm)*k+b, (s2*m+mm)*k+b2] = Vr[s, b, s2, b2]
VV = np.kron(Vfull, I(2))
W1 = np.kron(P0, UU) + np.kron(P1, VV)
Zq = np.diag([1, -1])
CZ = np.kron(P0, I(D*m*k*2)) + np.kron(P1, np.kron(I(D*m*k), Zq))
W = CZ @ W1
assert np.allclose(W @ W.conj().T, I(2*D*m*k*2))
PsiHatBuilt = factor_channel(W, 2*D, m*k*2)
PsiHatFormula = hat(Psi)
Jb = choi(PsiHatBuilt, 2*D); Jf = choi(PsiHatFormula, 2*D)
assert np.allclose(Jb, Jf)
print("built factorization == block formula:", np.allclose(Jb, Jf))
dJ_hat = tn(Jf - choi(hat(PhiH), 2*D)) / 2
dJ = tn(choi(Psi, D) - choi(PhiH, D)) / 2
assert np.isclose(dJ_hat, dJ/2)
print("d_J(PsiHat, PhiHat) =", dJ_hat, " half d_J(Psi, PhiH) =", dJ/2)
