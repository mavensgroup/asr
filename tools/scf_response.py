#!/usr/bin/env python3
"""Linear response of the SCF from an [asr] ef_scan_dump file (route edc, spin 2, at least 3 scan energies): the
self-consistent change of the moments and charges that follows from a change after one iteration at a fixed potential.
See docs/DERIVATION_scf_response.md, part C.

From the dump, per component p and spin s (1/Ry per atom):
  n_ps    = dq_ps / dE_F                                thermal DOS: response to a uniform shift of the potential
  chi_ps  = (1/pi) Im sum_j w_j sum_L G_pL^2             on-site response: shift on the site alone (Result B1)
  B_pp's  = w_p' (1/pi) Im sum_j w_j sum_L (-G'_pL - G_pL^2) kappa_p'L,   kappa_p' = G_p'^2 / sum_p'' w_p'' G_p''^2
            response of a p site to a shift on all other sites of species p' (single-site vertex weights);
            sum_p' B_pp' = n_p - chi_p holds exactly.
From asr.out of the same potential: the moments m_p and the Stoner parameters I_p = (C_d,down - C_d,up) / m_p.

Model.  A change (dm, dq) of the output moments and charges changes the input potential of spin s = +-1 by
  dv_ps = -s I_p dm_p / 2 + kappa dq_p          (exchange splitting; kappa: hardness of the component charge, Ry/e)
and the output by  dn_ps = -chi_ps dv_ps - sum_p' B_pp's dv_p's + n_ps dmu  (dmu from neutrality).  With the change
(dm0, dq0) at the fixed potential as the source, the self-consistent change solves the linear system
  dm = dm0 + sum_s s dn_s,    dq = dq0 + sum_s dn_s,    sum_p w_p dq_p = 0.

Usage:
  scf_response.py chains.dat asr.out                          response matrices and the gain matrix dm / dm0
  scf_response.py chains.dat asr.out --dm0 a,b,.. --dq0 a,b,..  the predicted self-consistent change
Options: --kappa x (default 1.0), --closure vertex | conc | onsite (default vertex)
"""
import re
import sys
import numpy as np
from f4_dump import read_dump, tail_hull, chain_green


def response(path):
    hdr, scans = read_dump(path)
    ncomp, nsp = hdr['ncomp'], hdr['nsp']
    if nsp != 2:
        sys.exit('scf_response: the dump must be spin polarised')
    w = hdr['w']
    sc = {s['js']: s for s in scans}
    if not all(k in sc for k in (-1, 0, 1)):
        sys.exit('scf_response: the dump needs the scan energies -1, 0, 1 (ef_scan >= 1)')

    def greens(s):
        out = {}
        for (j, o, p, g), (a, b) in s['chains'].items():
            if j <= s['nq'] and len(a) > 0:
                z = s['z'][j]
                out[(j, o, p, g)] = chain_green(a, b, z, tail_hull(a, b, z))
        return out
    g0, gp, gm = greens(sc[0]), greens(sc[1]), greens(sc[-1])
    n = np.zeros((ncomp, 2))
    chi = np.zeros((ncomp, 2))
    bm = np.zeros((ncomp, ncomp, 2))
    for g in (1, 2):
        for o in range(1, 5):
            wo = hdr['worb'][o - 1]
            for j in range(1, sc[0]['nq'] + 1):
                wj = sc[0]['w'][j]
                gs = np.array([g0[(j, o, p, g)] for p in range(1, ncomp + 1)])
                # the main line of the contour is stretched with E_F: use the shift of the node itself.  Nodes that do
                # not move lie far below the bands, where G' = -G^2.
                dz = sc[1]['z'][j] - sc[-1]['z'][j]
                if abs(dz) > 1e-4:
                    dg = np.array([(gp[(j, o, p, g)] - gm[(j, o, p, g)]) / dz for p in range(1, ncomp + 1)])
                else:
                    dg = -gs**2
                kap = gs**2 / np.sum(w * gs**2)
                n[:, g - 1] += -wo * (wj * dg).imag / np.pi
                chi[:, g - 1] += wo * (wj * gs**2).imag / np.pi
                for p in range(ncomp):
                    bm[p, :, g - 1] += wo * w * (wj * (-dg[p] - gs[p]**2) * kap).imag / np.pi
    return w, n, chi, bm


def stoner(path):
    """moments and I = (C_d,down - C_d,up) / m per component from the converged-results block of asr.out"""
    names, ms, cu, cd = [], [], [], []
    for ln in open(path):
        f = ln.split()
        if ln.startswith(' class '):
            names.append(f[2])
        elif 'spin moment in sphere' in ln:
            ms.append(float(f[4]))
        elif len(f) > 3 and f[1] == 'd' and f[0] in ('1', '2') and re.match(r'^-?\d', f[2]):
            (cu if f[0] == '1' else cd).append(float(f[3]))
    k = len(names)
    m = np.array(ms[:k])
    return names, m, (np.array(cd[:k]) - np.array(cu[:k])) / m


def solve(w, n, chi, bm, stoner_i, dm0, dq0, kappa=1.0, closure='vertex', gain=False):
    nc = len(w)
    if closure == 'onsite':
        b = np.zeros((nc, nc, 2))
    elif closure == 'conc':
        b = (n - chi)[:, None, :] * np.asarray(w)[None, :, None]
    else:
        b = bm
    r = b.copy()
    for p in range(nc):
        r[p, p, :] += chi[p, :]
    nn = 2 * nc + 1
    a = np.eye(nn)
    rhs = np.zeros(nn)
    sg = (1.0, -1.0)
    for p in range(nc):
        rhs[p], rhs[nc + p] = dm0[p], dq0[p]
        for s in range(2):
            for row, f in ((p, sg[s]), (nc + p, 1.0)):
                for p2 in range(nc):
                    a[row, p2] -= f * r[p, p2, s] * sg[s] * stoner_i[p2] / 2
                    a[row, nc + p2] += f * r[p, p2, s] * kappa
                a[row, 2 * nc] -= f * n[p, s]
    a[2 * nc, :] = 0.0
    a[2 * nc, nc:2 * nc] = w
    if gain:
        return np.linalg.inv(a)[:nc, :nc]
    x = np.linalg.solve(a, rhs)
    return x[:nc], x[nc:2 * nc]


def main():
    args = sys.argv[1:]
    if len(args) < 2:
        sys.exit(__doc__)

    def opt(name, default=None):
        return args[args.index(name) + 1] if name in args else default
    w, n, chi, bm = response(args[0])
    names, m, si = stoner(args[1])
    kappa = float(opt('--kappa', '1.0'))
    closure = opt('--closure', 'vertex')
    np.set_printoptions(precision=4, suppress=True, linewidth=200)
    print('components      :', ' '.join(names))
    print('moments (mu_B)  :', m)
    print('Stoner I (Ry)   :', si)
    print('n   up / down   :', n[:, 0], '/', n[:, 1])
    print('chi up / down   :', chi[:, 0], '/', chi[:, 1])
    print('B up (row: responding component, column: shifted species):\n', bm[:, :, 0])
    print('B down:\n', bm[:, :, 1])
    print('gain matrix dm / dm0 (closure %s, kappa %.2f):\n' % (closure, kappa),
          solve(w, n, chi, bm, si, [0] * len(w), [0] * len(w), kappa, closure, gain=True).round(2))
    if '--dm0' in args:
        dm0 = [float(x) for x in opt('--dm0').split(',')]
        dq0 = [float(x) for x in opt('--dq0', ','.join(['0'] * len(w))).split(',')]
        dm, dq = solve(w, n, chi, bm, si, dm0, dq0, kappa, closure)
        print('fixed potential : dm0', np.array(dm0), ' dq0', np.array(dq0))
        print('self-consistent : dm ', dm, ' dq ', dq)
        print('                  dm / |m| (%)', 100 * dm / abs(m), '  new moments', m + dm)


if __name__ == '__main__':
    main()
