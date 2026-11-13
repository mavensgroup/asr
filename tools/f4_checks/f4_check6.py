import numpy as np, sys
from f4lab import *
np.set_printoptions(linewidth=200, precision=6)
A, M, H, T = make_pencil(L=7, seed=2, tvals=(1.0, 0.25, 3.0), eps=(0.0, 0.15, -0.2))
N = A.shape[0]

# ---- (f) real axis: slope of a Ritz value,  kappa' = x^T M_m x + 2 beta_m x_m (c . x),  c = Q_m'^T q_{m+1}
print("(f) slope of the Ritz values of K(E) on K_m(K(E), e1): formula against finite differences; sign = sign of the residue")
m = 12; E0 = 0.13; h = 1e-6
def real_chain(E):
    al, be, Q = cs_lanczos((E + 0j) * M - A, m)
    return al.real, be.real, Q.real
al, be, Q = real_chain(E0); alp, bep, Qp = real_chain(E0 + h); alm, bem, Qm = real_chain(E0 - h)
Tm = tri(al, be, m).real; kap, X = np.linalg.eigh(Tm)
kp = np.linalg.eigvalsh(tri(alp, bep, m).real); km = np.linalg.eigvalsh(tri(alm, bem, m).real)
dQ = (Qp - Qm) / (2 * h)
c = dQ[:, :m].T @ Q[:, m]
Mm = Q[:, :m].T @ M @ Q[:, :m]
print("   Ritz value   slope (FD)   formula      x^T M_m x    rotation term   Ritz residual beta_m |x_m|")
nneg = 0
for i in range(m):
    x = X[:, i]; t1 = x @ Mm @ x; t2 = 2 * be[m - 1] * x[m - 1] * (c @ x)
    print("   %+.5f    %+.5f     %+.5f     %.5f      %+.5f        %.2e" % (kap[i], (kp[i] - km[i]) / (2 * h), t1 + t2, t1, t2, abs(be[m - 1] * x[m - 1])))
# statistics over E and m: how many Ritz values have negative slope, and their residuals
print("   statistics, 200 energies:")
for m in (12, 24):
    neg = tot = 0; resneg = []; respos = []
    for E0 in np.linspace(-0.8, 0.8, 200):
        al, be, Q = real_chain(E0) if m == 12 else (lambda r: (r[0].real, r[1].real, r[2].real))(cs_lanczos((E0 + 0j) * M - A, m))
        r2 = cs_lanczos((E0 + h + 0j) * M - A, m); r3 = cs_lanczos((E0 - h + 0j) * M - A, m)
        k0, X = np.linalg.eigh(tri(al, be, m).real)
        sl = (np.linalg.eigvalsh(tri(r2[0], r2[1], m).real) - np.linalg.eigvalsh(tri(r3[0], r3[1], m).real)) / (2 * h)
        res = np.abs(be[m - 1] * X[m - 1, :])
        tot += m; neg += (sl < 0).sum(); resneg += list(res[sl < 0]); respos += list(res[sl >= 0])
    print("   m = %d: %d of %d Ritz values have negative slope; their residuals: min %.2e median %.2e; residuals of the others: median %.2e; mu_min = %.3f"
          % (m, neg, tot, min(resneg) if resneg else 0, np.median(resneg) if resneg else 0, np.median(respos), np.diag(M).min()))

# ---- (g) pointwise positive alternatives on the same Krylov space, tail-free, at contour-like nodes
print("\n(g) tail-free error |G_approx - G| at nodes: bilinear chain (code), Hermitian Galerkin on V (FOM), Hermitian Galerkin on span(Re V, Im V)")
e1 = np.zeros(N); e1[0] = 1
print("   z                 m   chain       FOM         real span    -ImG>0: chain / FOM / real span")
for z in (0.05 + 0.1j, -0.2 + 0.1j, 0.3 + 0.05j, 0.05 + 0.02j):
    K = z * M - A; Gx = g_exact(A, M, z)
    al, be, Q = cs_lanczos(K, 36)
    for m in (12, 24, 36):
        Gc = g_chain(al, be, m)
        V, _ = np.linalg.qr(Q[:, :m])                       # unitary basis of the complex Krylov space
        Gf = (V.conj().T @ e1) @ np.linalg.solve(V.conj().T @ K @ V, V.conj().T @ e1)
        R = np.hstack([Q[:, :m].real, Q[:, :m].imag]); U, s, _ = np.linalg.svd(R, full_matrices=False); U = U[:, s > 1e-10 * s[0]]
        Gr = (U.T @ e1) @ np.linalg.solve(U.T @ K @ U, U.T @ e1)
        print("   %-16s %3d  %.2e    %.2e    %.2e     %s / %s / %s   (dim %d)" % (z, m, abs(Gc - Gx), abs(Gf - Gx), abs(Gr - Gx), Gc.imag < 0, Gf.imag < 0, Gr.imag < 0, U.shape[1]))
