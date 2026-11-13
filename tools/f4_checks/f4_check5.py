import numpy as np, sys
from scipy.optimize import nnls
from f4lab import *
np.set_printoptions(linewidth=200, precision=6)
A, M, H, T = make_pencil(L=7, seed=2, tvals=(1.0, 0.25, 3.0), eps=(0.0, 0.15, -0.2))
N = A.shape[0]; e1 = np.zeros(N); e1[0] = 1

# ---- (d) causal completion: Galerkin pencil on the real word space reproduces the m-step chain at every z
print("(d) word-space Galerkin pencil against the full pencil: m-step chain values at arbitrary z")
for m in (3, 4, 5):
    words = [e1]; level = [e1]
    for _ in range(m - 1):
        level = [X @ v for v in level for X in (A, M)]
        words += level
    Wm = np.array(words).T
    U, s, _ = np.linalg.svd(Wm, full_matrices=False)
    Qw = U[:, s > 1e-10 * s[0]]
    Aw = Qw.T @ A @ Qw; Mw = Qw.T @ M @ Qw; sw = Qw.T @ e1
    # rotate so that the start vector is the first basis vector
    Hh = np.eye(len(sw)) - 2 * np.outer(sw - np.eye(len(sw))[0], sw - np.eye(len(sw))[0]) / max((sw - np.eye(len(sw))[0]) @ (sw - np.eye(len(sw))[0]), 1e-300)
    Aw = Hh @ Aw @ Hh; Mw = Hh @ Mw @ Hh
    d = 0.0
    for z in (0.1 + 0.05j, -0.4 + 0.3j, 0.7 + 0.001j):
        al, be, _ = cs_lanczos(z * M - A, m); al2, be2, _ = cs_lanczos(z * Mw - Aw, m)
        d = max(d, abs(g_chain(al, be, m) - g_chain(al2, be2, m)))
    print("   m = %d: dim W = %3d (N = %d), eig(M_W) > 0: %s, max |G_m(full) - G_m(W)| = %.1e" % (m, Qw.shape[1], N, np.linalg.eigvalsh(Mw).min() > 0, d))

# ---- (e) Pick matrix of the contour data and the positive-measure (NNLS) projection
print("\n(e) Pick matrix and positive-measure projection of contour data (Fermi contour, kT 0.001, 16 poles, 32 + 16 + 8 + 4 nodes)")
def gl(n): return np.polynomial.legendre.leggauss(n)
def fermi_contour(ebot, mu, kt, npole, n_main, n_step):
    gam = 2 * npole * np.pi * kt; zs = []; ws = []
    def f(e):
        return 1.0 / (np.exp(np.clip(((e - mu) / kt).real, -600, 600) + 1j * ((e - mu) / kt).imag) + 1.0)
    x, w = gl(4)
    for xi, wi in zip(x, w):
        zz = ebot + 0.5j * gam * (1 + xi); zs.append(zz); ws.append(0.5j * gam * wi * f(zz))
    for lo, hi, n in ((ebot, mu - 10 * kt, n_main), (mu - 10 * kt, mu + 10 * kt, n_step), (mu + 10 * kt, mu + 40 * kt, 8)):
        x, w = gl(n)
        for xi, wi in zip(x, w):
            e = 0.5 * (lo + hi) + 0.5 * (hi - lo) * xi; zs.append(e + 1j * gam); ws.append(0.5 * (hi - lo) * wi * f(e + 0j))
    for j in range(npole):
        zs.append(mu + 1j * (2 * j + 1) * np.pi * kt); ws.append(-2j * np.pi * kt)
    return np.array(zs), np.array(ws), gam
def pick(zs, G):
    f = -G
    return (f[:, None] - f.conj()[None, :]) / (zs[:, None] - zs.conj()[None, :])
wgen, xgen = np.linalg.eigh(np.diag(1 / np.sqrt(np.diag(M))) @ A @ np.diag(1 / np.sqrt(np.diag(M))))
ebot = wgen.min() - 0.3
grid = np.linspace(wgen.min() - 0.2, wgen.max() + 0.2, 6001)
m = 24
print("   mu      |  max err  hull / tail-free / projected   | N-N_exact  hull / tail-free / proj  | Pick lam_min/lam_max  exact / hull / tail-free | lower bound on max err (hull)")
for mu in (-0.30, -0.10, 0.05, 0.20, 0.35):
    zs, ws, gam = fermi_contour(ebot, mu, 0.001, 16, 32, 16)
    Gx = np.array([g_exact(A, M, z) for z in zs]); G0 = np.zeros_like(Gx); Gh = np.zeros_like(Gx)
    for i, z in enumerate(zs):
        al, be, _ = cs_lanczos(z * M - A, m)
        G0[i] = g_chain(al, be, m); Gh[i] = g_chain(al, be, m, hull_tail_K(al, be, m, z))
    def lam(G):
        ev = np.linalg.eigvalsh(pick(zs, G)); return ev.min() / ev.max()
    # projection: positive measure on the grid, least squares on the node values (real and imaginary parts)
    Cm = 1.0 / (zs[:, None] - grid[None, :])
    X = np.vstack([Cm.real, Cm.imag]); rhs = np.concatenate([Gh.real, Gh.imag])
    wts, rn = nnls(X, rhs, maxiter=20000)
    Gp = Cm @ wts
    Nx, N0, Nh, Np = [-(ws * g).sum().imag / np.pi for g in (Gx, G0, Gh, Gp)]
    P = pick(zs, Gh); lmin = np.linalg.eigvalsh(P).min()
    S = np.sqrt((1.0 / np.abs(zs[:, None] - zs.conj()[None, :]) ** 2).sum())
    lb = max(0.0, -lmin) / (2 * S)
    print("   %+.2f   |  %.2e / %.2e / %.2e      | %+.2e / %+.2e / %+.2e | %+.1e / %+.1e / %+.1e | %.1e   (sum w = %.4f)" % (
        mu, abs(Gh - Gx).max(), abs(G0 - Gx).max(), abs(Gp - Gx).max(), Nh - Nx, N0 - Nx, Np - Nx, lam(Gx), lam(Gh), lam(G0), lb, wts.sum()))
    sys.stdout.flush()
