#!/usr/bin/env python3
"""Causality analysis of an [asr] ef_scan_dump file (route edc).  See docs/DERIVATION_f4.md.

The dump holds, for every scan energy, the Fermi contour (nodes z, weights w) and every chain (a_k, b_k) of
asr_edc::edc_green.  From them this tool rebuilds the Green function with the code's `hull` tail and without a tail and
reports, per scan energy and in total:

  neg     nodes with -Im G <= 0 (a Herglotz function has none)
  weyl    nodes where the 2 x 2 matrix -Im [[g_11, g_1m], [g_1m, g_mm]] of (z - T_m)^-1 is not positive definite
          (then some tail with Im t <= 0 gives -Im G <= 0: the head of the chain alone does not guarantee causality)
  timag   nodes where the hull tail itself has Im t > 0
  chi     bare on-site charge response of every component, chi_p = (1/pi) Im sum_j w_j sum_{L, spin} G^2 (1/Ry; see
          docs/DERIVATION_scf_response.md, Result B1), at the central scan energy
  N       electron count from the dump (must equal the `ef_scan` line of the run), and its increments along the scan
          (a decreasing N(mu) is a negative thermal DOS)
  vortex  phase vortices of G on the grid (scan energy) x (Matsubara frequency): +1 = zero, -1 = pole of G inside the
          cell.  A Herglotz function has neither in the upper half plane.

Usage: f4_dump.py chains.dat [--verbose]
"""
import sys
import numpy as np


def read_dump(path):
    scans = []
    with open(path) as fh:
        lines = fh.read().split('\n')
    i = 0
    hdr = {}
    cur = None
    while i < len(lines):
        ln = lines[i]
        if not ln.strip():
            i += 1
            continue
        tag = ln[:2]
        f = ln[2:].split()
        if tag == 'H ':
            hdr = dict(nelec=float(f[0]), fs=float(f[1]), ncomp=int(f[2]), nsp=int(f[3]), nlev=int(f[4]))
        elif tag == 'W ':
            hdr['w'] = np.array([float(x) for x in f])
        elif tag == 'O ':
            hdr['worb'] = np.array([float(x) for x in f])
        elif tag == 'E ':
            cur = dict(js=int(f[0]), ev=float(f[1]), nq=int(f[2]), z={}, w={}, chains={})
            scans.append(cur)
        elif tag == 'Z ':
            j = int(f[0])
            cur['z'][j] = complex(float(f[1]), float(f[2]))
            cur['w'][j] = complex(float(f[3]), float(f[4]))
        elif tag == 'C ':
            j, o, p, g = (int(x) for x in f)
            m = int(lines[i + 1])
            ab = np.array([[float(x) for x in lines[i + 2 + k].split()] for k in range(m)]).reshape(m, 4)
            cur['chains'][(j, o, p, g)] = (ab[:, 0] + 1j * ab[:, 1], ab[:, 2] + 1j * ab[:, 3])
            i += 1 + m
        i += 1
    return hdr, scans


def tail_hull(a, b, z):
    """asr_endep::cf_species, ed_tail = 'hull'."""
    m = len(a)
    if m < 2:
        return 0.0
    tm = np.diag(a).astype(complex)
    for k in range(m - 1):
        tm[k, k + 1] = tm[k + 1, k] = b[k]
    ev = np.linalg.eigvals(tm)
    lo = ev[np.argmin(ev.real)]
    hi = ev[np.argmax(ev.real)]
    aa = 0.5 * (lo + hi)
    bb = 0.25 * (hi - lo)
    if abs(bb) < 1e-14:
        return 0.0
    d = z - aa
    s = np.sqrt(d * d - 4 * bb * bb)
    g1 = (d - s) / (2 * bb * bb)
    g2 = (d + s) / (2 * bb * bb)
    return bb * bb * (g1 if abs(bb * g1) <= abs(bb * g2) else g2)


def chain_green(a, b, z, t):
    m = len(a)
    g = 1.0 / (z - a[m - 1] - t)
    for k in range(m - 2, -1, -1):
        g = 1.0 / (z - a[k] - b[k] ** 2 * g)
    return g


def corner(a, b, z):
    """g_11, g_1m, g_mm of (z - T_m)^-1 by the forward and backward fractions."""
    m = len(a)
    up = np.zeros(m, dtype=complex)       # up[k]: [(z - T_{k..m})^-1]_kk
    up[m - 1] = 1.0 / (z - a[m - 1])
    for k in range(m - 2, -1, -1):
        up[k] = 1.0 / (z - a[k] - b[k] ** 2 * up[k + 1])
    dn = 1.0 / (z - a[0])                 # [(z - T_{1..k})^-1]_kk
    for k in range(1, m):
        dn = 1.0 / (z - a[k] - b[k - 1] ** 2 * dn)
    g1m = up[0]
    for k in range(m - 1):
        g1m = g1m * b[k] * up[k + 1]
    return up[0], g1m, dn


def vortices(G):
    """winding number of arg G around every cell of a 2D array G[i, j]."""
    ph = np.angle(G)

    def d(x, y):
        return (y - x + np.pi) % (2 * np.pi) - np.pi
    w = d(ph[:-1, :-1], ph[1:, :-1]) + d(ph[1:, :-1], ph[1:, 1:]) + d(ph[1:, 1:], ph[:-1, 1:]) + d(ph[:-1, 1:], ph[:-1, :-1])
    return np.rint(w / (2 * np.pi)).astype(int)


def main():
    path = sys.argv[1]
    verbose = '--verbose' in sys.argv
    hdr, scans = read_dump(path)
    ncomp, nsp, nlev = hdr['ncomp'], hdr['nsp'], hdr['nlev']
    print('dump %s: %d scan energies, %d components, %d spins, %d levels, N_elec %.6f' % (path, len(scans), ncomp, nsp, nlev, hdr['nelec']))
    part_names = ['vert', 'main', 'step', 'tail', 'poles', 'dos']
    tot = dict(nodes=0, neg_h=0, neg_0=0, weyl=0, timag=0)
    by_part = {k: dict(nodes=0, neg_h=0, neg_0=0, weyl=0) for k in part_names}
    worst = (0.0, None)
    Ns = []
    chis = {}
    grids = {}
    print('  js      E (Ry)     N - N_elec    dN        nodes  neg(hull) neg(no tail) weyl  timag')
    for sc in scans:
        nq = sc['nq']
        zs = sc['z']
        npole = 0                                # Matsubara nodes: the last nodes, all with Re z = scan energy
        for j in range(nq, 0, -1):
            if abs(zs[j].real - sc['ev']) > 1e-12:
                break
            npole += 1
        n_main = nq - 4 - 16 - 8 - npole
        bnd = [4, 4 + n_main, 4 + n_main + 16, 4 + n_main + 24, nq]

        def part(j):
            for ip, bd in enumerate(bnd):
                if j <= bd:
                    return part_names[ip]
            return 'dos'
        loc = dict(nodes=0, neg_h=0, neg_0=0, weyl=0, timag=0)
        N = 0.0
        chi = np.zeros(ncomp)
        for (j, o, p, g), (a, b) in sc['chains'].items():
            z = zs[j]
            m = len(a)
            if m == 0:
                continue
            t = tail_hull(a, b, z)
            Gh = chain_green(a, b, z, t)
            g11, g1m, gmm = corner(a, b, z)
            W = -np.array([[g11.imag, g1m.imag], [g1m.imag, gmm.imag]])
            wbad = np.linalg.eigvalsh(W).min() <= 0
            pt = part(j)
            if j <= nq:
                N += -hdr['fs'] * hdr['w'][p - 1] * hdr['worb'][o - 1] * (sc['w'][j] * Gh).imag / np.pi
                chi[p - 1] += hdr['fs'] * hdr['worb'][o - 1] * (sc['w'][j] * Gh * Gh).imag / np.pi
            for dct in (loc, tot, by_part[pt]):
                dct['nodes'] += 1
                dct['neg_h'] += Gh.imag >= 0
                dct['neg_0'] += g11.imag >= 0
                dct['weyl'] += wbad
            loc['timag'] += t.imag > 0
            tot['timag'] += t.imag > 0
            if Gh.imag > worst[0]:
                worst = (Gh.imag, (sc['js'], j, o, p, g, z))
            if pt == 'poles':
                grids.setdefault((o, p, g), {})[(sc['js'], j - bnd[3])] = (Gh, g11)
        Ns.append(N)
        chis[sc['js']] = chi
        dN = N - Ns[-2] if len(Ns) > 1 else 0.0
        print('  %3d  %11.6f  %+.6e  %+.2e  %6d  %6d  %6d  %6d  %6d' % (sc['js'], sc['ev'], N - hdr['nelec'], dN, loc['nodes'], loc['neg_h'],
                                                                 loc['neg_0'], loc['weyl'], loc['timag']))
    print('total: %d node values; -Im G <= 0: hull %d, no tail %d; Weyl test fails at %d; hull tail with Im t > 0 at %d' % (
        tot['nodes'], tot['neg_h'], tot['neg_0'], tot['weyl'], tot['timag']))
    for k in part_names:
        v = by_part[k]
        if v['nodes']:
            print('  part %-5s: %7d values, neg(hull) %5d, neg(no tail) %5d, weyl %6d' % (k, v['nodes'], v['neg_h'], v['neg_0'], v['weyl']))
    if worst[1] is not None:
        print('worst Im G (hull) = %+.4f at scan %d node %d orbital %d component %d spin %d, z = %s' % ((worst[0],) + worst[1]))
    jmid = min(chis, key=abs)
    print('bare on-site response chi_p (1/Ry) at scan %d: %s;  1/chi_p (Ry/e): %s' % (
        jmid, ' '.join('%.4f' % c for c in chis[jmid]), ' '.join('%.4f' % (1.0 / c) for c in chis[jmid])))
    dNs = np.diff(Ns)
    print('N(mu) increments: min %+.3e, max %+.3e; decreasing steps: %d of %d' % (dNs.min(), dNs.max(), (dNs <= 0).sum(), len(dNs)))
    # vortices on the (scan energy, Matsubara index) grid
    jss = sorted({sc['js'] for sc in scans})
    nv_h = nv_0 = 0
    for key, gd in sorted(grids.items()):
        rows = sorted({r for (_, r) in gd})
        if len(jss) < 2 or len(rows) < 2:
            continue
        Gh = np.array([[gd[(js, r)][0] for r in rows] for js in jss])
        G0 = np.array([[gd[(js, r)][1] for r in rows] for js in jss])
        vh = vortices(Gh)
        v0 = vortices(G0)
        nv_h += np.abs(vh).sum()
        nv_0 += np.abs(v0).sum()
        if verbose or np.abs(vh).sum() or np.abs(v0).sum():
            zs0 = {sc['js']: sc for sc in scans}
            for name, v in (('hull', vh), ('no tail', v0)):
                for (i, r) in zip(*np.nonzero(v)):
                    sc = zs0[jss[i]]
                    print('  vortex %+d (%s) orbital %d component %d spin %d: E in [%.4f, %.4f], Matsubara index %d to %d' % (
                        v[i, r], name, key[0], key[1], key[2], sc['ev'], zs0[jss[i + 1]]['ev'], rows[r], rows[r + 1]))
    print('phase vortices on the Matsubara grid: hull %d, no tail %d (a Herglotz function has none)' % (nv_h, nv_0))


if __name__ == '__main__':
    main()
