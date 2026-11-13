"""Checks 1: error identity, disc bound, positivity certificate, where negative -Im G occurs."""
import numpy as np
from f4lab import *

np.set_printoptions(linewidth=200, precision=5)
A, M, H, T = make_pencil(L=7, seed=3)
N = A.shape[0]
Minv = np.diag(1.0 / np.diag(M))
print("N =", N, " M range", np.diag(M).min(), np.diag(M).max())
wgen = np.linalg.eigvalsh(np.linalg.solve(np.sqrt(M), np.linalg.solve(np.sqrt(M), A).T))
print("spectrum of the pencil: %.4f .. %.4f" % (wgen.min(), wgen.max()))

# ---- 1. identity  G - G_m = rho^2 q_{m+1}^T K^-1 q_{m+1}  and the bound  |.| <= |rho|^2 q^H M^-1 q / Im z
print("\n1. error identity and disc bound")
print("   z                 m   |G-Gm|      |identity err|  bound       bound/|G-Gm|  ||q||^2")
for z in (0.05 + 0.1j, -0.2 + 0.1j, 0.05 + 0.02j, 0.05 + 0.00314j):
    K = z * M - A
    al, be, Q = cs_lanczos(K, 40)
    Gx = g_exact(A, M, z)
    for m in (8, 16, 24, 36):
        Gm = g_chain(al, be, m)
        Km = tri(al, be, m)
        rho = np.prod(be[:m]) / np.linalg.det(Km)
        q = Q[:, m]
        gf = q @ np.linalg.solve(K, q)
        ident = abs((Gx - Gm) - rho ** 2 * gf)
        bound = abs(rho) ** 2 * (q.conj() @ Minv @ q).real / z.imag
        print("   %-16s %3d  %.3e   %.1e        %.3e   %8.1f      %.3f" % (z, m, abs(Gx - Gm), ident, bound, bound / abs(Gx - Gm),
                                                                      (q.conj() @ q).real))

# ---- 2. scan: negative -Im G, tail-free and with tails; certificate Y_m = Im(Ktilde_m) > 0
print("\n2. scan over E at fixed Im z: count of -Im G <= 0;  certificate Y_m > 0")
for y in (0.1, 0.02, 0.005, 0.001):
    for m in (12, 24, 36):
        nneg0 = nnegh = nnegr = ncert_fail = ncert_ok_neg = 0
        worst = 0.0
        Es = np.linspace(-0.75, 0.75, 301)
        for E in Es:
            z = E + 1j * y
            K = z * M - A
            al, be, Q = cs_lanczos(K, m)
            Km = tri(al, be, m)
            Y = Km.imag
            cert = np.linalg.eigvalsh(Y).min() > 0
            G0 = g_chain(al, be, m)
            Gh = g_chain(al, be, m, hull_tail_K(al, be, m, z))
            tr = hull_tail_K(al, be, m, z, real_edges=True)
            Gr = g_chain(al, be, m, tr)
            nneg0 += G0.imag >= 0
            nnegh += Gh.imag >= 0
            nnegr += Gr.imag >= 0
            worst = max(worst, Gh.imag)
            ncert_fail += (not cert)
            if cert and (G0.imag >= 0 or (Gr.imag >= 0 and tr.imag <= 0)):
                ncert_ok_neg += 1
        print("   Im z %.3f  m %2d:  -ImG<=0: tail-free %3d, hull %3d, real-edge %3d of %d;  Y_m not PD at %3d;"
              "  certificate violated %d;  worst Im G(hull) %+.3f" % (y, m, nneg0, nnegh, nnegr, len(Es), ncert_fail, ncert_ok_neg, worst))
