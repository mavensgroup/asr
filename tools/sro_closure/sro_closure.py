#!/usr/bin/env python3
"""The closure of the first-shell SRO (docs/FINDINGS_LOG.md, "B3 closure"): a structure factor that is positive by construction.

Krivoglaz-Clapp-Moss / spherical-model form for one SRO mode with first-shell value alpha_1 = lambda / z:
    alpha(k) = N / (1 - t c1(k)),   N = 1 / <1/(1 - t c1)>   (so <alpha(k)> = 1),   alpha_1 = <c1 alpha> = (1 - N) / t.
t runs from the ordering limit (1 - t c1min = 0) to the clustering limit t = 1, where alpha_1 = alpha_c (fcc 0.256).  Above
alpha_c the spherical model condenses at k = 0: alpha(k) = (1 - m) alpha_c(k) + m N_sites delta(k), alpha_1 = (1 - m) alpha_c + m.
G(k = 0) = 0, so the condensate carries no charge:
    Gamma_KCM = <G(k) alpha(k)> = (1 - m) Gamma_c   ->  0 at alpha_1 = 1 for every beta.  F = a (1 / Gamma - 1) has no pole.
Printed against the first-shell model Gamma_1 = G00 + G1 z alpha_1 of asr_screening.

Usage: python3 -I sro_closure.py [nk]     (default 40; offset mesh)
"""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sro_pole_scan import lattice, solve  # noqa: E402


def kcm(g, c1, alpha1, c1min):
    """Gamma_KCM and t (or the condensate fraction m) for the first-shell alpha_1."""
    def a1(t):
        n = 1 / np.mean(1 / (1 - t * c1))
        return (1 - n) / t if abs(t) > 1e-12 else 0.0
    ac = a1(1.0)
    if alpha1 >= ac:
        m = (alpha1 - ac) / (1 - ac)
        return (1 - m) * np.mean(g / (1 - c1)) / np.mean(1 / (1 - c1)), 1.0, m
    lo, hi = 1 / c1min * (1 - 1e-9), 1.0                       # alpha_1 grows with t
    if alpha1 <= a1(lo):
        return np.nan, lo, 0.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if a1(mid) < alpha1:
            lo = mid
        else:
            hi = mid
    t = 0.5 * (lo + hi)
    return np.mean(g / (1 - t * c1)) / np.mean(1 / (1 - t * c1)), t, 0.0


def main():
    nk = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    for lat, alat, betas in (('fcc', 6.78412, (0.8, 1.0, 1.02, 1.2)), ('bcc', 5.38572, (1.0, 1.2))):
        m, c1, r1, z, amin = lattice(lat, alat, nk)
        c1min = {'fcc': -1 / 3, 'bcc': -1.0}[lat]
        print(f'{lat}, nk {nk}: alpha_c = {(1 - 1 / np.mean(1 / (1 - c1))):.4f} (mesh value; exact fcc 0.2563, bcc 0.2822)')
        for beta in betas:
            a, g00, g1, _ = solve(m, c1, r1, amin, beta)
            g = a / (a + m)
            print(f'  beta {beta}: a {a:.4f}  G00 {g00:.4f}  z G1 {z * g1:+.4f}  first-shell pole at alpha_1 = {-g00 / (z * g1):.3f}')
            print('    alpha_1  lambda   Gamma_1   Gamma_KCM   F_1 (Ry/e)   F_KCM (Ry/e)   t or m')
            for al in (-0.30, -0.266, -0.148, -0.05, 0.0, 0.05, 0.10, 0.20, 0.25, 0.34, 0.50, 0.66, 0.80, 0.90, 0.95, 1.0):
                gk, t, mm = kcm(g, c1, al, c1min)
                g1s = g00 + g1 * z * al
                tag = f'm = {mm:.3f}' if mm > 0 else f't = {t:+.3f}'
                fk = a * (1 / gk - 1) if gk > 1e-12 else np.inf
                print(f'    {al:+6.3f}  {z * al:+6.2f}  {g1s:+8.4f}  {gk:+8.4f}   {a * (1 / g1s - 1):+9.3f}   {fk:+9.3f}       {tag}')


if __name__ == '__main__':
    main()
