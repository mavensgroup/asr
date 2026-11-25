#!/usr/bin/env python3
"""Exact Gamma of the degenerate alpha_1 = -1/3 ground states of the first-shell fcc model (L = 8 supercell).

Decides the ordering end point of the closure on fcc.  Thermal Monte Carlo freezes near alpha_1 = -0.30 and never reaches the
end segment, so the ground states are built explicitly and enumerated completely.
x = 1/2: s = sigma_j (-1)^i with i = h_x, j = h_y (h = 2 r / a), any signs sigma_j: 2^8 stackings, all enumerated.
x = 1/4: minority on even planes h_z, on (h_x, h_y) = (e, e) or (o, o) chosen per plane pair (L1_2, DO22, ...): 2^4 stackings.
Each family has alpha_1 = -1/3 exactly (checked).  At x = 1/2 the stacking spectrum of random signs is flat, so the ensemble
fills the X-W line (1, q, 0) 2 pi / a evenly: its mean Gamma is the uniform line average of G (equal to 4e-16 here).
Usage: python3 -I sro_fcc_degenerate.py
"""
import os
import sys, itertools, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sro_basis_test as t
_proto = t.proto
def proto(name, coa=None):
    if name == 'fcc':
        return 0.5 * np.array([[0.0, 1, 1], [1, 0, 1], [1, 1, 0]]).T, None, np.zeros((3, 1)), np.array([1])
    return _proto(name, coa)
t.proto = proto
alat = 3.59 / 0.529177210903
L = 8
avec, _, tau0, lab = proto('fcc')
p = avec * alat
for beta in (1.0, 1.2):
    mdl = t.branch_model(avec, tau0, lab, alat, beta=beta, nk=12, ksro=0)
    a = mdl['a']
    sc = t.BasisSupercell('fcc', None, alat, L=L, ksro=0, a=a)
    idx = np.array(np.meshgrid(*[np.arange(L)] * 3, indexing='ij')).reshape(3, -1)   # (i1, i2, i3), cell_id order
    h = np.vstack([idx[1] + idx[2], idx[0] + idx[2], idx[0] + idx[1]])                 # h = 2 r / a
    G = lambda kc: a / (a + t.coulomb_k_basis(p, np.zeros((3, 1)), np.atleast_2d(kc).T * 2 * np.pi / alat)[0, 0, :].real)
    q = 2.0 * np.arange(L) / L                                                          # commensurate stacking wavevectors
    g_line_d = G(np.vstack([np.ones(L), q, np.zeros(L)]).T)
    s = (np.arange(256) + 0.5) / 256 * 2
    g_line_c = G(np.vstack([np.ones_like(s), s, np.zeros_like(s)]).T)
    gX, gW = G([1, 0, 0])[0], G([1, 0.5, 0])[0]
    print(f"\nfcc beta {beta}: a={a:.6f}  G(X)={gX:.6f}  G(W)={gW:.6f}  X-W line average: continuous {g_line_c.mean():.6f}, "
          f"commensurate (L={L}) {g_line_d.mean():.6f}")
    # x = 1/2: all 2^L sign sequences
    gs, als = [], []
    for sig in itertools.product((1, -1), repeat=L):
        sgn = np.array(sig)[h[1] % L] * (-1) ** (h[0] % 2)
        al1, gex = sc.exact_gamma_binary((sgn > 0).astype(np.int8))
        gs.append(gex); als.append(al1)
    gs = np.array(gs)
    print(f"  x=1/2, {len(gs)} stackings: alpha_1 in [{min(als):+.12f}, {max(als):+.12f}]; Gamma_exact mean {gs.mean():.6f} "
          f"(commensurate line average {g_line_d.mean():.6f}, diff {gs.mean() - g_line_d.mean():+.1e}); range [{gs.min():.6f}, {gs.max():.6f}]")
    # x = 1/4: minority on even h_z planes, sublattice (e,e) or (o,o) per plane pair n = h_z / 2 (period L / 2)
    gs4, als4 = [], []
    for tt in itertools.product((0, 1), repeat=L // 2):
        tsel = np.array(tt)[(h[2] // 2) % (L // 2)]
        occ = ((h[2] % 2 == 0) & (h[0] % 2 == tsel) & (h[1] % 2 == tsel)).astype(np.int8)
        al1, gex = sc.exact_gamma_binary(occ)
        gs4.append(gex); als4.append(al1)
    gs4 = np.array(gs4)
    print(f"  x=1/4, {len(gs4)} stackings: mean occupation {occ.mean():.4f}; alpha_1 in [{min(als4):+.12f}, {max(als4):+.12f}]; "
          f"Gamma_exact mean {gs4.mean():.6f}; range [{gs4.min():.6f}, {gs4.max():.6f}] (L1_2 {gs4[0]:.6f})")
    print(f"  former code (largest, W) {gW:.6f}; code now (uniform line average) {g_line_c.mean():.6f}; numpy closure g_ord {mdl['g_ord']:.6f}")
