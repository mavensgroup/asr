import numpy as np, sys, itertools
from fractions import Fraction
from f4lab import *
np.set_printoptions(linewidth=200, precision=6)

# ---- (a) necessity of the Weyl-disc criterion, exact construction of the worst real tail
print("(a) necessity: W not PD  =>  a real tail value gives Im G >= 0")
A, M, H, T = make_pencil(L=7, seed=1, tvals=(1.0, 0.25, 3.0), eps=(0.0, 0.15, -0.2))
nbad = nfound = 0
for y in (0.1, 0.005):
    for E in np.linspace(-0.9, 0.9, 201):
        z = E + 1j * y
        al, be, Q = cs_lanczos(z * M - A, 24)
        Ki = np.linalg.inv(tri(al, be, 24)); g11, g1m, gmm = Ki[0, 0], Ki[0, 23], Ki[23, 23]
        W = -np.array([[g11.imag, g1m.imag], [g1m.imag, gmm.imag]])
        if np.linalg.eigvalsh(W).min() > 0: continue
        nbad += 1
        if gmm.imag >= 0 or g11.imag >= 0:
            nfound += 1; continue                      # t = 0 or t -> infinity side already non-causal / unbounded image
        c = g11 + 1j * g1m ** 2 / (2 * gmm.imag); r = abs(g1m) ** 2 / (2 * abs(gmm.imag))
        Gtop = c + 1j * r
        t = (Gtop - g11) / (g1m ** 2 + gmm * (Gtop - g11))
        ok = abs(t.imag) < 1e-6 * abs(t) and Gtop.imag >= -1e-12
        nfound += ok
print("   W not PD at", nbad, "points; real tail with Im G >= 0 constructed at", nfound)

# ---- (b) minimal example: m = 3.  Search small integer pencils (N = 4) for a non-real pole of G_3
print("\n(b) minimal example, m = 3, N = 4, integer entries")
def G3_poly(A, Mdiag):
    """numerator and denominator of G_3(z) as numpy poly1d, from the Hankel form with moments s_j(z) = e1^T K^j e1."""
    N = len(Mdiag)
    # vectors K^j e1 as polynomial-valued: list of arrays of poly coefficients -> use poly1d objects per component
    z = np.poly1d([1.0, 0.0])
    K = [[z * (Mdiag[i] if i == j else 0.0) - A[i][j] for j in range(N)] for i in range(N)]
    v = [np.poly1d([1.0]) if i == 0 else np.poly1d([0.0]) for i in range(N)]
    vs = [v]
    for _ in range(3):
        v = [sum((K[i][j] * v[j] for j in range(N)), np.poly1d([0.0])) for i in range(N)]
        vs.append(v)
    def s(j):  # s_j = (K^a e1).(K^b e1), a + b = j
        a = j // 2; b = j - a
        return sum((vs[a][i] * vs[b][i] for i in range(N)), np.poly1d([0.0]))
    S = [s(j) for j in range(6)]
    H1 = [[S[i + j + 1] for j in range(3)] for i in range(3)]
    def det3(X):
        return (X[0][0] * (X[1][1] * X[2][2] - X[1][2] * X[2][1]) - X[0][1] * (X[1][0] * X[2][2] - X[1][2] * X[2][0])
                + X[0][2] * (X[1][0] * X[2][1] - X[1][1] * X[2][0]))
    D = det3(H1)
    num = np.poly1d([0.0])
    for i in range(3):                    # G = s^T H1^-1 s = sum_i s_i * det(H1 with column i replaced by s) / det
        Xi = [[(S[r] if cc == i else H1[r][cc]) for cc in range(3)] for r in range(3)]
        num = num + S[i] * det3(Xi)
    return num, D
rng = np.random.default_rng(11)
found = None
for trial in range(20000):
    Mdiag = [1] + list(rng.integers(1, 5, size=3))
    Am = rng.integers(-2, 3, size=(4, 4)); Am = np.triu(Am) + np.triu(Am, 1).T
    num, D = G3_poly(Am.tolist(), Mdiag)
    # cancel common factors numerically: poles = roots of D that are not roots of num
    rD = D.roots
    poles = [r for r in rD if abs(num(r)) > 1e-6 * (1 + abs(np.polyval(np.abs(num.coeffs), abs(r))))]
    cplx = [r for r in poles if abs(r.imag) > 1e-6]
    if cplx:
        # verify against the direct chain at a test point
        z0 = 0.3 + 0.7j
        al, be, Q = cs_lanczos(z0 * np.diag(Mdiag) - Am, 3)
        if abs(g_chain(al, be, 3) - num(z0) / D(z0)) < 1e-9:
            found = (Mdiag, Am, num, D, poles); break
Mdiag, Am, num, D, poles = found
print("   M = diag", [int(x) for x in Mdiag]); print("   A =\n", Am)
print("   poles of G_3:", np.array(sorted(poles, key=lambda r: (r.real, r.imag))))
print("   exact generalized eigenvalues:", np.linalg.eigvalsh(np.diag(1 / np.sqrt(Mdiag)) @ Am @ np.diag(1 / np.sqrt(Mdiag))))
zp = max(poles, key=lambda r: r.imag)
res = num(zp) / D.deriv()(zp)
print("   UHP pole z_p = %s, residue %s" % (zp, res))
zt = zp.real + 0.5j * zp.imag
al, be, Q = cs_lanczos(zt * np.diag(Mdiag) - Am, 3)
print("   -Im G_3 at z = Re z_p + i Im z_p / 2: %.4f   (exact -Im G: %.4f)" % (-g_chain(al, be, 3).imag, -g_exact(Am, np.diag(Mdiag), zt).imag))
np.save("minimal_example.npy", dict(M=Mdiag, A=Am), allow_pickle=True)

# ---- (c) contour formula:  N(mu) = sum_real r f(E_p) + 2 Re sum_{UHP, Im z_p < Gamma} r_p f(z_p)   (tail-free, Fermi contour)
print("\n(c) Fermi-contour charge of the tail-free G_3 against the pole formula")
def gl(n): return np.polynomial.legendre.leggauss(n)
def fermi_contour(ebot, mu, kt, npole, n_main, n_step):
    gam = 2 * npole * np.pi * kt; zs = []; ws = []
    f = lambda e: 1.0 / (np.exp((e - mu) / kt) + 1.0)
    x, w = gl(4)
    for xi, wi in zip(x, w):
        zz = ebot + 0.5j * gam * (1 + xi); zs.append(zz); ws.append(0.5j * gam * wi * f(zz))
    for lo, hi, n in ((ebot, mu - 10 * kt, n_main), (mu - 10 * kt, mu + 10 * kt, n_step), (mu + 10 * kt, mu + 40 * kt, 8)):
        x, w = gl(n)
        for xi, wi in zip(x, w):
            e = 0.5 * (lo + hi) + 0.5 * (hi - lo) * xi; zs.append(e + 1j * gam); ws.append(0.5 * (hi - lo) * wi * f(e))
    for j in range(npole):
        zs.append(mu + 1j * (2 * j + 1) * np.pi * kt); ws.append(-2j * np.pi * kt)
    return np.array(zs), np.array(ws), gam
Gm = lambda z: num(z) / D(z)
allres = [(p, num(p) / D.deriv()(p)) for p in poles]
print("   sum of all residues (weight): %.6f" % sum(r for _, r in allres).real)
kt = 0.02; npole = 12
ebot = min(p.real for p in poles) - 3.0
print("   kT = %.3f, Gamma = %.4f;  UHP pole at Im z_p = %.4f" % (kt, 2 * npole * np.pi * kt, zp.imag))
for mu in np.linspace(zp.real - 6 * kt, zp.real + 6 * kt, 7):
    zs, ws, gam = fermi_contour(ebot, mu, kt, npole, 400, 64)
    Nc = -(ws * Gm(zs)).sum().imag / np.pi
    fz = lambda zz: 1.0 / (np.exp((zz - mu) / kt) + 1.0)
    Nf = 0.0; Nreal = 0.0
    for p, r in allres:
        if abs(p.imag) < 1e-9: Nf += (r * fz(p.real)).real; Nreal += (r * fz(p.real)).real
        elif 0 < p.imag < gam: Nf += 2 * (r * fz(p)).real
    print("   mu = %+.4f: contour N = %+.6f, pole formula = %+.6f, real poles only = %+.6f" % (mu, Nc, Nf, Nreal))
