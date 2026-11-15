#!/usr/bin/env python3
"""Numeric check of the S_trunc identity, the positivity argument and the energy condition (docs/FINDINGS_LOG.md, "B3 closure").

Gamma = G00 I + G1 Ahat is the single-site response (dq = Gamma db).  Claims checked here on the code's own lattice functions:
 (1) diag(x) Gamma = <G(k) S_trunc(k)>_BZ on the deviation subspace x.db = 0, S_trunc(k) = diag(x) [(I - 1 x^T) + c1(k) Ahat].
 (2) S_trunc(k) >= 0 for all k  <=>  every eigenvalue of Ahat lies in [-1/c1max, -1/c1min] = [-1, 3] (fcc), [-1, 1] (bcc).
 (3) G00 + G1 lambda > 0 on that interval for every beta < beta_max, so the pole lambda* = -G00/G1 lies outside it.
 (4) Single-site functional Phi(q) = sum_p x_p a/2 (q_p - b_p)^2 + 1/2 q^T diag(x) F q: stationary at q = Gamma db, value
     1/2 db^T diag(x) (W00 + W1 Ahat) db, Hessian diag(x) (a + F) = a diag(x) Gamma^-1.  Minimum <=> Gamma > 0.

Usage: python3 -I sro_strunc_check.py [nk]   (default 16)
"""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))          # python3 -I drops the script directory
from sro_pole_scan import lattice, solve


def random_ahat(rng, n, z, strength):
    """Ahat = z (P - 1 x^T) from a random symmetric pair-probability matrix (rows of P sum to 1, x_p P(p'|p) symmetric)."""
    x = rng.dirichlet(np.ones(n) * 3)
    y = np.outer(x, x) * (1 - strength) + strength * np.diag(x)          # pair probabilities: random + clustered
    noise = rng.uniform(-1, 1, (n, n)); noise = noise + noise.T
    noise -= noise.mean(0)[None, :]; noise -= noise.mean(1)[:, None]      # zero row and column sums keep the marginals
    y = y + 0.1 * (1 - strength) * x.min() ** 2 * noise
    assert (y > 0).all() and np.allclose(y.sum(1), x)
    p = y / x[:, None]
    return x, z * (p - x[None, :])


def main():
    nk = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    rng = np.random.default_rng(1)
    for lat, alat in [('fcc', 6.78412), ('bcc', 5.38572)]:
        m, c1, r1, z, amin = lattice(lat, alat, nk)
        c1min, c1max = {'fcc': (-1 / 3, 1.0), 'bcc': (-1.0, 1.0)}[lat]
        lo, hi = -1 / c1max, -1 / c1min
        print(f'{lat}: z {z}, c1 in [{c1min:.4f}, {c1max:.1f}], realizable first-shell-only range of lambda: [{lo:.0f}, {hi:.0f}]')
        # (1) identity, 3 and 4 components
        a, g00, g1, bmax = solve(m, c1, r1, amin, 1.0)
        g = a / (a + m)
        err = 0.0
        for n in (2, 3, 4):
            for strength in (0.1, 0.5, 0.9):
                x, ah = random_ahat(rng, n, z, strength)
                lhs = np.diag(x) @ (g00 * np.eye(n) + g1 * ah)
                rhs = sum(gk * np.diag(x) @ ((np.eye(n) - np.outer(np.ones(n), x)) + ck * ah) for gk, ck in zip(g, c1)) / g.size
                d = rng.normal(size=(n, 5)); d -= (x @ d)[None, :]        # deviation vectors: x.d = 0
                err = max(err, np.abs(d.T @ (lhs - rhs) @ d).max(), np.abs(lhs - lhs.T).max(), np.abs(x @ ah).max(),
                          np.abs(ah @ np.ones(n)).max())
        print(f'  (1) max error of the identity, of the symmetry of diag(x) Gamma, of x^T Ahat = 0, Ahat 1 = 0: {err:.1e}')
        # (2), (3): pole against beta, and the minimum of G00 + G1 lambda on the realizable interval
        print('  (3) beta   a/a_min   G00      G1       lambda*   min over [lo, hi] of (G00 + G1 lambda)')
        for beta in (0.2, 0.5, 0.8, 1.0, 1.1, 1.2, 1.3, 0.999 * bmax):
            a, g00, g1, _ = solve(m, c1, r1, amin, beta)
            print(f'      {beta:5.3f}  {a / amin:7.3f}  {g00:7.4f}  {g1:+8.4f}  {-g00 / g1:8.3f}   {min(g00 + g1 * lo, g00 + g1 * hi):.4f}')
        # (4) energy functional at beta 1 for a binary mode, lambda below and above the pole
        a, g00, g1, _ = solve(m, c1, r1, amin, 1.0)
        w00, w1 = a * (1 - g00), -a * g1
        print(f'  (4) beta 1.0, lambda* = {-g00 / g1:.3f}: single-site functional per unit x_A x_B db^2')
        print('      lambda   Gamma     F (Ry/e)   a + F     Phi(stationary)  1/2 (W00 + W1 lambda)   kind')
        for lam in (-2.0, 0.0, 3.0, 7.92, 12.0, 12.3, 12.4, 13.0, 16.0):
            gam = g00 + g1 * lam
            f = a * (1 / gam - 1)
            q = gam                                                         # dq per unit db
            phi = 0.5 * a * (q - 1) ** 2 + 0.5 * f * q * q
            kind = 'minimum' if a + f > 0 else 'MAXIMUM (unbounded below)'
            print(f'      {lam:6.2f}  {gam:+8.4f}  {f:+9.3f}  {a + f:+9.3f}  {phi:+10.5f}      {0.5 * (w00 + w1 * lam):+10.5f}        {kind}')


if __name__ == '__main__':
    main()
