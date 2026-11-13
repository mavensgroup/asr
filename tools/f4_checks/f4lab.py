"""F4 sandbox: the energy-dependent chain as a bilinear Krylov-Galerkin approximant of the real pencil (A, M).

K(z) = z M - A,  M = T^2 > 0 (site diagonal), A = T H T real symmetric, T = 1 on the start site.
The code's chain (asr_edc::edc_green) is the complex-symmetric Lanczos chain of H(z) = z - K(z) started on e1,
G_m(z) = [ (z - T_m(z))^-1 ]_11 (+ tail).
"""
import numpy as np


def make_pencil(L=5, seed=0, tvals=(1.0, 0.6, 1.5), conc=(0.4, 0.3, 0.3), eps=(0.0, 0.08, -0.10), hop=1.0 / 12.0,
                norb=1):
    """Simple cubic L^3 cluster (open boundary), nn hopping, species-dependent site energy and T.  Site 0 = centre."""
    rng = np.random.default_rng(seed)
    idx = {}
    sites = []
    c = L // 2
    order = sorted(((abs(i - c) + abs(j - c) + abs(k - c), i, j, k) for i in range(L) for j in range(L) for k in range(L)))
    for n, (_, i, j, k) in enumerate(order):
        idx[(i, j, k)] = n
        sites.append((i, j, k))
    N = len(sites)
    H = np.zeros((N, N))
    for (i, j, k), n in idx.items():
        for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            nb = (i + d[0], j + d[1], k + d[2])
            if nb in idx:
                H[n, idx[nb]] = H[idx[nb], n] = -hop
    sp = rng.choice(len(tvals), size=N, p=conc)
    sp[0] = 0
    T = np.array([tvals[s] for s in sp])
    H[np.arange(N), np.arange(N)] = [eps[s] for s in sp]
    A = T[:, None] * H * T[None, :]
    M = np.diag(T ** 2)
    return A, M, H, T


def cs_lanczos(K, m, start=0, reorth=True):
    """Complex-symmetric (bilinear) Lanczos of K from e_start: returns alpha, beta (K-chain), Q (N x (m+1))."""
    N = K.shape[0]
    Q = np.zeros((N, m + 1), dtype=complex)
    al = np.zeros(m, dtype=complex)
    be = np.zeros(m, dtype=complex)
    Q[start, 0] = 1.0
    for k in range(m):
        w = K @ Q[:, k]
        al[k] = Q[:, k] @ w
        w = w - al[k] * Q[:, k]
        if k > 0:
            w = w - be[k - 1] * Q[:, k - 1]
        if reorth:
            for _ in range(2):
                w = w - Q[:, :k + 1] @ (Q[:, :k + 1].T @ w)
        be[k] = np.sqrt(w @ w)
        Q[:, k + 1] = w / be[k]
    return al, be, Q


def tri(al, be, m):
    Tm = np.diag(al[:m]).astype(complex)
    for k in range(m - 1):
        Tm[k, k + 1] = Tm[k + 1, k] = be[k]
    return Tm


def g_chain(al, be, m, t=0.0):
    """[ (Ktilde_m - t e_m e_m^T)^-1 ]_11 for the K-chain (alpha, beta); t = tail self-energy on the last level."""
    g = 1.0 / (al[m - 1] - t)
    for k in range(m - 2, -1, -1):
        g = 1.0 / (al[k] - be[k] ** 2 * g)
    return g


def hull_tail_K(al, be, m, z, real_edges=False):
    """The code's 'hull' tail (cf_species) written for the K-chain.  The code works with T_m = z - Ktilde_m and spectral
    variable z: edges from the eigenvalues of T_m with extreme real part; tail = bb^2 g."""
    Tm = z * np.eye(m) - tri(al, be, m)        # H-chain: a_k = z - alpha_k, b_k = -beta_k (sign irrelevant)
    ev = np.linalg.eigvals(Tm)
    lo = ev[np.argmin(ev.real)]
    hi = ev[np.argmax(ev.real)]
    if real_edges:
        lo, hi = lo.real, hi.real
    aa = 0.5 * (lo + hi)
    bb = 0.25 * (hi - lo)
    d = z - aa
    s = np.sqrt(d * d - 4 * bb * bb + 0j)
    g1 = (d - s) / (2 * bb * bb)
    g2 = (d + s) / (2 * bb * bb)
    return bb * bb * (g1 if abs(bb * g1) <= abs(bb * g2) else g2)


def g_exact(A, M, z, start=0):
    N = A.shape[0]
    e = np.zeros(N)
    e[start] = 1.0
    return np.linalg.solve(z * M - A, e)[start]
