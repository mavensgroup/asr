#!/usr/bin/env python3
"""Pole of the SRO screening matrix F of asr_screening against sim_beta (notes/FINDINGS_LOG.md, "B3 closure").

numpy port of asr_screening::screen_lattice (coulomb_k, the bisection for the hardness a, G00, G1).  It prints, per lattice and
beta: a, G00, G1, the pole lambda* = -G00/G1 of F(lambda) = (W00 + W1 lambda) / (G00 + G1 lambda), and F at the Ahat eigenvalues
of the FeNiCr gate triads (0, 7.92, 12), plus the beta at which lambda* reaches z (full clustering of one species).

Usage: python3 -I sro_pole_scan.py [nk]      (nk = BZ mesh per direction; the code uses 28; 16 is faster and close)
"""
import sys
from math import erfc, sqrt, pi
import numpy as np


def prim(lat):
    if lat == 'fcc':
        return np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float).T * 0.5
    return np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]], dtype=float).T * 0.5


def coulomb_k(p, kv):
    """M(k) = sum_{R /= 0} (2/R) cos(k.R) (Ry/e), Ewald, as asr_screening::coulomb_k."""
    om = abs(np.linalg.det(p))
    b = 2 * pi * np.linalg.inv(p).T
    eta = sqrt(pi) / om ** (1 / 3)
    rmax = 7 / eta
    gmax = 14 * eta + np.linalg.norm(b, axis=0).max()
    n = [int(np.ceil(rmax * np.linalg.norm(b[:, i]) / (2 * pi))) + 1 for i in range(3)]
    idx = np.array(np.meshgrid(*[np.arange(-k, k + 1) for k in n], indexing='ij')).reshape(3, -1)
    rv = p @ idx
    r = np.linalg.norm(rv, axis=0)
    keep = (r > 1e-10) & (r <= rmax)
    rv, r = rv[:, keep], r[keep]
    wr = np.array([erfc(eta * x) / x for x in r])
    n = [int(np.ceil(gmax * np.linalg.norm(p[:, i]) / (2 * pi))) + 1 for i in range(3)]
    idx = np.array(np.meshgrid(*[np.arange(-k, k + 1) for k in n], indexing='ij')).reshape(3, -1)
    gv = b @ idx
    gv = gv[:, np.linalg.norm(gv, axis=0) <= gmax]
    out = np.empty(kv.shape[1])
    for j, k in enumerate(kv.T):
        v = k[:, None] + gv
        q2 = (v * v).sum(0)
        s = 4 * pi / om * np.sum(np.exp(-q2 / (4 * eta ** 2)) / q2) + np.sum(np.cos(k @ rv) * wr)
        out[j] = 2 * (s - 2 * eta / sqrt(pi))
    return out


def lattice(lat, alat, nk):
    p = prim(lat) * alat
    b = 2 * pi * np.linalg.inv(p).T
    idx = np.array(np.meshgrid(*[np.arange(-2, 3)] * 3, indexing='ij')).reshape(3, -1)
    v = p @ idx
    r = np.linalg.norm(v, axis=0)
    rmin = r[r > 1e-10].min()
    sh = v[:, np.abs(r - rmin) < 1e-8 * rmin]
    f = (np.array(np.meshgrid(*[np.arange(nk)] * 3, indexing='ij')).reshape(3, -1) + 0.5) / nk
    kv = b @ f
    m = coulomb_k(p, kv)
    c1 = np.cos(kv.T @ sh).mean(1)
    special = {'fcc': [[1, 0, 0], [1, .5, 0], [.5, .5, .5]], 'bcc': [[1, 0, 0], [.5, .5, 0], [.5, .5, .5]]}[lat]
    amin = max(-m.min(), -coulomb_k(p, np.array(special, dtype=float).T * 2 * pi / alat).min())
    return m, c1, rmin, sh.shape[1], amin


def solve(m, c1, r1, amin, beta):
    """hardness a for beta (bisection on log(a - a_min), as the code), then G00, G1, beta_max."""
    def beta_of(a):
        g = np.mean(a / (a + m))
        return a * (g - 1) * r1 / (2 * g)
    lo, hi = np.log(amin * 1e-9), np.log(1e9)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if beta_of(amin + np.exp(mid)) > beta:
            lo = mid
        else:
            hi = mid
    a = amin + np.exp(0.5 * (lo + hi))
    g = a / (a + m)
    return a, g.mean(), (g * c1).mean(), beta_of(amin * (1 + 1e-9))


def main():
    nk = int(sys.argv[1]) if len(sys.argv) > 1 else 28
    for lat, alat in [('fcc', 6.78412), ('bcc', 5.38572)]:
        m, c1, r1, z, amin = lattice(lat, alat, nk)
        print(f'{lat}, alat {alat} Bohr, nk {nk}: z {z}, R1 {r1:.4f}, z*c1 on the mesh in [{z * c1.min():.3f}, {z * c1.max():.3f}]')
        lo, hi = 0.8, 1.3
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            a, g00, g1, _ = solve(m, c1, r1, amin, mid)
            if -g00 / g1 > z:
                lo = mid
            else:
                hi = mid
        print(f'  lambda* reaches z = {z} at beta = {0.5 * (lo + hi):.4f}')
        for beta in [0.5, 0.8, 0.9, 0.95, 1.0, 1.01, 1.02, 1.05, 1.1, 1.2, 1.3]:
            a, g00, g1, bmax = solve(m, c1, r1, amin, beta)
            if beta >= bmax:
                print(f'  beta {beta:4.2f}: above beta_max {bmax:.4f}')
                continue
            w00, w1 = a * (1 - g00), -a * g1
            fl = [(w00 + w1 * lam) / (g00 + g1 * lam) for lam in (0.0, 7.92, float(z))]
            print(f'  beta {beta:4.2f}: a {a:8.4f}  G00 {g00:.4f}  G1 {g1:+.4f}  lambda* {-g00 / g1:8.3f}  '
                  f'F(0) {fl[0]:+7.3f}  F(7.92) {fl[1]:+8.3f}  F(z) {fl[2]:+9.3f} Ry/e')


if __name__ == '__main__':
    main()
