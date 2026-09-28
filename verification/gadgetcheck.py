# Finite examples for the gadget-entropy proposition: Choi rank r = #distinct elements carried,
# flat complementary output at rho_* = (1/r) sum_g |i_g><i_g|, H = log r.
import numpy as np, itertools
def compose(p, q):  # p then q (right action): x -> q[p[x]]
    return tuple(q[i] for i in p)
def inv(p):
    r = [0]*len(p)
    for i, j in enumerate(p): r[j] = i
    return tuple(r)
def closure(gens):
    idp = tuple(range(len(gens[0]))); G = {idp}; frontier = [idp]
    while frontier:
        new = []
        for g in frontier:
            for s in gens:
                h = compose(g, s)
                if h not in G: G.add(h); new.append(h)
        frontier = new
    return sorted(G)
def run(gens, relators, name):
    # gens: dict name->perm ; relators: list of lists of (name, +/-1)
    elems = closure(list(gens.values())); N = len(elems); idx = {g: i for i, g in enumerate(elems)}
    e = tuple(range(len(elems[0])))
    def val(word):
        g = e
        for (s, sg) in word: g = compose(g, gens[s] if sg == 1 else inv(gens[s]))
        return g
    def lam(g):  # right regular rep: basis x -> x g^{-1}  (any faithful unitary rep with trace delta works)
        M = np.zeros((N, N))
        for x in elems: M[idx[compose(x, inv(g))], idx[x]] = 1
        return M
    signed = [((s, 1),) for s in gens] + [((s, -1),) for s in gens]
    symbols = list(signed)
    for r in relators:
        for l in range(2, len(r)): symbols.append(tuple(r[:l]))
    triangles = []
    for r in relators:
        L = len(r)
        for l in range(2, L+1):
            x = tuple(r[:l-1]); y = (r[l-1],); z = tuple(r[:l]) if l < L else ()
            triangles.append((x, y, z))
    for s in gens: triangles.append((((s, 1),), ((s, -1),), ()))
    d = 1 + len(symbols) + 2*len(triangles)
    p, k = len(gens), len(relators)
    assert d == 1 + 4*p + 3*sum(len(r) for r in relators) - 4*k
    # slots: list of (row, col, coeff, element)
    slots = [(0, 0, 1.0, e)]
    for i, w in enumerate(symbols): slots.append((1+i, 1+i, 1.0, val(w)))
    off = 1 + len(symbols); c = 2**-0.5
    for t, (x, y, z) in enumerate(triangles):
        a, b = off + 2*t, off + 2*t + 1
        gx, gy = val(x), val(y)
        slots += [(a, a, c, e), (a, b, c, gy), (b, a, c, gx), (b, b, -c, compose(gx, gy))]
    U = np.zeros((d*N, d*N))
    for (i, j, co, g) in slots: U[i*N:(i+1)*N, j*N:(j+1)*N] += co*lam(g)
    assert np.allclose(U @ U.T, np.eye(d*N)), "gadget not unitary"
    carried = sorted(set(g for (_, _, _, g) in slots)); r = len(carried)
    A = {g: np.zeros((d, d)) for g in carried}
    for (i, j, co, g) in slots: A[g][i, j] += co
    # channel via Kraus A_g (valid since tau(lam(g)lam(h)^*) = delta)
    # cross-check with the tracial formula on a random input
    rho = np.random.default_rng(0).normal(size=(d, d)); rho = rho @ rho.T; rho /= np.trace(rho)
    Y = U @ np.kron(rho, np.eye(N)) @ U.T / N
    out1 = np.einsum('iaja->ij', Y.reshape(d, N, d, N)); out2 = sum(A[g] @ rho @ A[g].T for g in carried)
    assert np.allclose(out1, out2)
    Kr = np.array([A[g].reshape(-1) for g in carried]); rankK = np.linalg.matrix_rank(Kr)
    # anchors: symbol/identity diagonal slot with coefficient 1 carrying g
    anchor = {}
    for (i, j, co, g) in slots[:1+len(symbols)]:
        anchor.setdefault(g, i)
    assert set(anchor) == set(carried), "some carried element has no anchor"
    rs = np.zeros((d, d))
    for g in carried: rs[anchor[g], anchor[g]] = 1/r
    comp = np.array([[np.trace(A[g] @ rs @ A[h].T) for h in carried] for g in carried])
    ev = np.linalg.eigvalsh(comp); S = -sum(x*np.log2(x) for x in ev if x > 1e-14)
    assert rankK == r and np.allclose(comp, np.eye(r)/r)
    print(f"{name}: |G|={N}, d={d}, distinct carried r={r}, Kraus rank={rankK}, comp flat={np.allclose(comp, np.eye(r)/r)}, S={S:.6f}, log r={np.log2(r):.6f}")
# S3 = <s,t | s^2, t^2, (st)^3>
s = (1, 0, 2); t = (0, 2, 1)
run({'s': s, 't': t}, [[('s',1),('s',1)], [('t',1),('t',1)], [('s',1),('t',1)]*3], "S3")
# Z4 = <s | s^4>
run({'s': (1, 2, 3, 0)}, [[('s',1)]*4], "Z4")
# Z2xZ2 = <s,t | s^2, t^2, s t s^-1 t^-1>, with prefix coincidences
run({'s': (1, 0, 2, 3), 't': (0, 1, 3, 2)}, [[('s',1),('s',1)], [('t',1),('t',1)], [('s',1),('t',1),('s',-1),('t',-1)]], "Z2xZ2")
