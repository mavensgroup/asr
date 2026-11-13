"""Where does -Im G <= 0 occur?  Vary disorder strength in T, seed, reorthogonalisation."""
import numpy as np, sys
from f4lab import *
def scan(A, M, y, m, reorth, nE=201):
    n0 = nh = 0; worst0 = worsth = 0.0
    for E in np.linspace(-0.9, 0.9, nE):
        z = E + 1j * y
        al, be, Q = cs_lanczos(z * M - A, m, reorth=reorth)
        G0 = g_chain(al, be, m); Gh = g_chain(al, be, m, hull_tail_K(al, be, m, z))
        n0 += G0.imag >= 0; nh += Gh.imag >= 0
        worst0 = max(worst0, G0.imag); worsth = max(worsth, Gh.imag)
    return n0, nh, worst0, worsth
for tv in ((1.0, 0.6, 1.5), (1.0, 0.4, 2.0), (1.0, 0.25, 3.0)):
    for seed in (1, 2, 3):
        A, M, H, T = make_pencil(L=7, seed=seed, tvals=tv, eps=(0.0, 0.15, -0.2))
        for y in (0.02, 0.005, 0.001):
            for m in (24, 48):
                r = [scan(A, M, y, m, ro) for ro in (True, False)]
                print("T %s seed %d  Im z %.3f m %2d | reorth: neg %3d/%3d worst %+.2f/%+.2f | plain: neg %3d/%3d worst %+.2f/%+.2f"
                      % (tv, seed, y, m, r[0][0], r[0][1], r[0][2], r[0][3], r[1][0], r[1][1], r[1][2], r[1][3]))
        sys.stdout.flush()
