#!/usr/bin/env python3
"""Exact charge-excess model on minimal-tail configurations against the closures (notes/FINDINGS_LOG.md, "B3 closure").

The ASR treats the sites beyond the first shell as uncorrelated, and the SQS cells of gate B3 (tools/sqs.py --alpha) are built the
same way.  This test makes such configurations by inverse Monte Carlo: pair couplings on shells 1..NS are adjusted until the
first-shell alpha has its target and the alphas of shells 2..NS are as close to zero as the lattice allows.  For every
configuration the exact Gamma (all lattice vectors of the periodic cell, as sro_exact_test.py) is compared with the closures at
the measured first-shell alpha: linear (first shell only, the form before the closure), clip (minimal-tail closure) and KCM
(spherical model).  F = a (1/Gamma - 1).

Usage: python3 -I sro_mintail_test.py [L]     (default 16; all cores, NPROC to limit; SRO_BUILD: directory of sro_mc.so)
"""
import ctypes
import itertools
import multiprocessing
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sro_exact_test import Cell, load_mc, ALAT  # noqa: E402
from sro_pole_scan import lattice, solve, coulomb_k, prim  # noqa: E402
from sro_closure import kcm, clip  # noqa: E402

NS = 6                                       # shells with couplings: 1 is the target, 2..NS are driven to zero


def shell_offsets(cell, ns):
    """Grid offsets of the lattice vectors of shells 1..ns and the shell index (0-based) of each."""
    r2s = cell.shell_r2[1:ns + 1]
    off, sh = [], []
    for d in itertools.product(range(-6, 7), repeat=3):
        d = np.array(d)
        if cell.lat == 'fcc' and d.sum() % 2:
            continue
        if cell.lat == 'bcc' and not (d[0] % 2 == d[1] % 2 == d[2] % 2):
            continue
        r2 = int((d ** 2).sum())
        if r2 in r2s:
            off.append(d); sh.append(int(np.where(r2s == r2)[0][0]))
    return np.array(off, dtype=np.int32), np.array(sh, dtype=np.int32)


def alphas(cell, spec, ns):
    x, pr = cell.correlations(spec)
    al = []
    for s in range(1, ns + 1):
        sel = cell.r2sorted == cell.shell_r2[s]
        al.append(1.0 - pr[sel, 0, 1].mean() / x[1])
    return x, pr, np.array(al)


def make_config(job):
    """Inverse Monte Carlo for one target first-shell alpha (binary, x = 0.5)."""
    lat, L, target, seed = job
    cell = Cell(lat, L)
    lib = load_mc()
    lib.mc_shells.restype = ctypes.c_long
    lib.mc_shells.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_long, ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p,
                              ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p, ctypes.c_long, ctypes.c_uint64]
    off, sh = shell_offsets(cell, NS)
    spec = cell.random_config([0.5, 0.5], np.random.default_rng(seed))
    tgt = np.zeros(NS); tgt[0] = target
    K = np.zeros(NS)
    for it in range(60):
        eps = np.zeros((NS, 2, 2))
        eps[:, 0, 1] = eps[:, 1, 0] = K
        eps = np.ascontiguousarray(eps)
        lib.mc_shells(cell.n, spec.ctypes.data, cell.nsite, cell.sites.ctypes.data, len(sh), off.ctypes.data, sh.ctypes.data, 2,
                      eps.ctypes.data, 40, seed + 1000 * it)
        _, _, al = alphas(cell, spec, NS)
        K += 1.5 * (tgt - al) * np.where(np.arange(NS) == 0, 1.0, 0.5)
    return spec


def main():
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    cases = [('fcc', 1.0, (0.05, 0.10, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)),
             ('fcc', 1.2, ()),
             ('bcc', 1.0, (-0.05, -0.10, -0.148, -0.20, -0.266, 0.10, 0.20, 0.40, 0.60, 0.80))]
    jobs = []
    for lat, beta, targets in cases:
        for t in targets:
            jobs.append((lat, L, t, 7 + int(abs(t) * 1000)))
    nproc = min(len(jobs), int(os.environ.get('NPROC', os.cpu_count())))
    with multiprocessing.get_context('fork').Pool(nproc) as pool:
        specs = pool.map(make_config, jobs, chunksize=1)
    configs = {}
    for (lat, L_, t, _), sp in zip(jobs, specs):
        configs.setdefault(lat, []).append((t, sp))
    for lat, beta, _ in cases:
        cell = Cell(lat, L)
        m, c1, r1, z, amin = lattice(lat, ALAT[lat], 28)
        a, g00, g1, _ = solve(m, c1, r1, amin, beta)
        g = a / (a + m)
        c1min = {'fcc': -1.0 / 3.0, 'bcc': -1.0}[lat]
        kord = {'fcc': [1.0, 0.5, 0.0], 'bcc': [1.0, 0.0, 0.0]}[lat]           # ordering point: fcc W, bcc H (lowest M at min c1)
        g_ord = a / (a + coulomb_k(prim(lat) * ALAT[lat], np.array([kord]).T * 2 * np.pi / ALAT[lat])[0])
        print(f'\n=== {lat} binary x = 0.5, L = {L} ({cell.nsite} sites), sim_beta {beta}: a {a:.4f}, linear range of lambda '
              f'[-1, {-1 / c1min:.0f}], pole of the linear form at {-g00 / g1:.2f}')
        print('  target | alpha_1..alpha_6 measured                         | lambda |  Gamma exact  linear    clip     KCM  |'
              '   F exact  linear    clip     KCM (Ry/e)')
        for t, spec in configs[lat]:
            x, pr = cell.correlations(spec)
            out, exact, zg, ps, g00s = cell.gammas(a, x, pr, 1)
            al = np.array([1.0 - pr[cell.r2sorted == cell.shell_r2[s], 0, 1].mean() / x[1] for s in range(1, 7)])
            lam = z * al[0]
            ge = exact[0, 0] - exact[0, 1]
            gl, gc, gk = g00 + g1 * lam, clip(g, c1, lam, z, c1min, g_ord), kcm(g, c1, al[0], c1min)[0]
            f = [a * (1 / v - 1) if v > 1e-9 else np.inf for v in (ge, gl, gc, gk)]
            print(f'  {t:+6.3f} | ' + ' '.join(f'{v:+6.3f}' for v in al) + f' | {lam:+6.2f} | {ge:+7.4f} {gl:+7.4f} {gc:+7.4f} {gk:+7.4f} |'
                  + ''.join(f' {v:+8.3f}' for v in f))


if __name__ == '__main__':
    main()
