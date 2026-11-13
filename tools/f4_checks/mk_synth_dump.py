import numpy as np
from f4lab import *
d = np.load("minimal_example.npy", allow_pickle=True).item()
M = np.diag(np.array(d["M"], dtype=float)); A = np.array(d["A"], dtype=float)
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
    return np.array(zs), np.array(ws)
m = 3; kt = 0.02
with open("synth.dat", "w") as u:
    u.write("H %25.16e%6.1f%6d%6d%6d\n" % (1.0, 1.0, 1, 1, m))
    u.write("W %25.16e\n" % 1.0); u.write("O    1.0   0.0   0.0   0.0\n")
    for js in range(-10, 11):
        mu = 1.5551 + js * 0.02
        zs, ws = fermi_contour(-6.0, mu, kt, 12, 32, 16)
        nq = len(zs)
        zs = np.append(zs, mu + 0.01j); ws = np.append(ws, 0.0)
        u.write("E %4d%25.16e%6d\n" % (js, mu, nq))
        for j in range(nq + 1):
            u.write("Z %5d%25.16e%25.16e%25.16e%25.16e\n" % (j + 1, zs[j].real, zs[j].imag, ws[j].real, ws[j].imag))
        for j in range(nq + 1):
            al, be, Q = cs_lanczos(zs[j] * M - A, m)
            a = zs[j] - al; b = be
            u.write("C %5d%5d%5d%5d\n%4d\n" % (j + 1, 1, 1, 1, m))
            for k in range(m):
                u.write("%25.16e%25.16e%25.16e%25.16e\n" % (a[k].real, a[k].imag, b[k].real, b[k].imag))
