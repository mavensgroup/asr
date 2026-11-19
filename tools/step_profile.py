#!/usr/bin/env python3
"""Step profile from an [asr] ef_scan_dump file (route edc): the charges and moments of every component that the same chains
give when they are cut after m' <= m steps and closed with the hull tail (the Lanczos coefficients of the m'-step chain are
the first m' of the m-step chain).  E_F is set again for every m' from the scan energies (N = N_elec, interpolation).
Usage: step_profile.py chains.dat [m1,m2,...] [--notail]
See docs/DERIVATION_levels_steps.md."""
import sys
import numpy as np
from f4_dump import read_dump, tail_hull, chain_green

hdr, scans = read_dump(sys.argv[1])
ncomp, nsp, nlev = hdr['ncomp'], hdr['nsp'], hdr['nlev']
ms = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else list(range(4, nlev + 1, 2))
tailfree = '--notail' in sys.argv
evs = np.array([sc['ev'] for sc in scans])
q = np.zeros((len(ms), len(scans), ncomp, nsp))
mused = np.zeros(len(ms))
for isc, sc in enumerate(scans):
    nq = sc['nq']
    for (j, o, p, g), (a, b) in sc['chains'].items():
        if j > nq or len(a) == 0:
            continue
        z = sc['z'][j]; wj = sc['w'][j]
        for im, m in enumerate(ms):
            mm = min(m, len(a))
            aa, bb = a[:mm], b[:mm]
            t = 0.0 if tailfree else tail_hull(aa, bb, z)
            G = chain_green(aa, bb, z, t)
            q[im, isc, p - 1, g - 1] += -hdr['worb'][o - 1] * (wj * G).imag / np.pi
            mused[im] += mm
fs = hdr['fs']; w = hdr['w']
print('# dump %s: %d scan energies, components %d, spins %d, steps %d%s' % (sys.argv[1], len(scans), ncomp, nsp, nlev, ', no tail' if tailfree else ', hull tail'))
print('# m   E_F(m)      ' + ' '.join('q%d      ' % (p + 1) for p in range(ncomp)) + (' '.join('m%d      ' % (p + 1) for p in range(ncomp)) if nsp == 2 else '') + ' DOS(E_F)')
for im, m in enumerate(ms):
    N = fs * np.einsum('p,spg->s', w, q[im])
    d = N - hdr['nelec']
    if len(scans) > 1:
        k = np.argsort(evs)
        ef = np.interp(0.0, d[k], evs[k]) if np.all(np.diff(d[k]) > 0) else evs[np.argmin(abs(d))]
        qe = np.array([[np.interp(ef, evs[k], q[im, k, p, g]) for g in range(nsp)] for p in range(ncomp)])
        dos = np.gradient(N[k], evs[k])[len(k) // 2]
    else:
        ef = evs[0]; qe = q[im, 0]; dos = float('nan')
    ch = fs * qe.sum(axis=1)
    out = '%3d %10.6f  ' % (m, ef) + ' '.join('%8.4f' % x for x in ch)
    if nsp == 2:
        out += '  ' + ' '.join('%8.4f' % (qe[p, 0] - qe[p, 1]) for p in range(ncomp))
    print(out + '  %7.2f' % dos)
