#!/usr/bin/env python3
"""Short-range order parameters from the effective pair interactions of asr (the [epi] section, file epi.toml).

Pair model.  A bond between two sites R apart with species c and s costs U(c, s) (asr_epi; band + screened Coulomb).
The energy of a configuration is the sum over the bonds of the listed shells.  Exchanges inside a sublattice keep the
composition, and site terms drop out.

Two solvers.
  Monte Carlo (default): canonical species exchange, Metropolis, kernel epi_mc.c.  Exact for the pair model.  N species,
      any number of shells, any decorated lattice of asr.  Output: Warren-Cowley alpha of every pair type, the energy and
      the heat capacity per site.
  --kcm: Krivoglaz-Clapp-Moss with the spherical-model sum rule (one sublattice only).  Analytic check at high temperature.
      It also gives the instability temperature of the random alloy and its wave vector.

Warren-Cowley convention of asr (asr_sro): P(c on a site, s on a shell site) = x_c y_s (1 - alpha_cs).  alpha < 0: the pair
is favoured.  alpha > 0: the pair avoids contact.

Usage:  sro_alpha.py epi.toml -T 600 800 1000 [--size 12] [--sweeps 4000] [--equil 1000] [--kcm] [--band-only]
The last block of the output is the [[sro]] input of asr for the first shell at each temperature.
Environment: SRO_BUILD (directory for epi_mc.so; default: the system temporary directory).
"""
import argparse
import ctypes
import itertools
import os
import subprocess
import sys
import tempfile

import numpy as np

KB_RY = 8.617333262e-5 / 13.605693123          # Boltzmann constant, Ry / K
HERE = os.path.dirname(os.path.abspath(__file__))


def load_epi(path, band_only=False):
    """Read epi.toml (asr_driver::write_epi) into the dictionary used by Model."""
    import tomllib
    with open(path, 'rb') as f:
        d = tomllib.load(f)
    subs = d['sublattice']
    names = [s['name'] for s in subs]
    pairs = []
    for p in d.get('pair', []):
        u = np.array(p['u_band'], float)
        if not band_only:
            u = u + np.array(p['u_es'], float)
        pairs.append(dict(sub=(names.index(p['sublattices'][0]), names.index(p['sublattices'][1])), shell=p['shell'],
                          vectors=np.array(p['vectors'], int), u=u))
    return dict(lattice=d['lattice'], alat=d['alat'], period=d.get('period', 1),
                sublattices=[dict(name=s['name'], species=list(s['species']), conc=np.array(s['concentrations'], float),
                                  sites=np.array(s['sites'], int)) for s in subs],
                pairs=pairs)


def on_graph(kind, h):
    h = np.asarray(h)
    if kind == 'sc':
        return np.all(h % 2 == 0, axis=-1)
    if kind == 'bcc':
        return (h[..., 0] % 2 == h[..., 1] % 2) & (h[..., 1] % 2 == h[..., 2] % 2)
    if kind == 'fcc':
        return h.sum(axis=-1) % 2 == 0
    raise ValueError('unknown graph lattice ' + kind)


class Model:
    """The pair model on a periodic box of size^3 conventional cubes (coordinates in units of alat / 2)."""

    def __init__(self, epi, size):
        self.epi = epi
        per = epi['period']
        if size % per:
            raise ValueError('size must be a multiple of the decoration period %d' % per)
        self.m = m = 2 * size
        g = np.stack(np.meshgrid(*[np.arange(m)] * 3, indexing='ij'), -1).reshape(-1, 3)
        g = g[on_graph(epi['lattice'], g)]
        self.coord = g
        self.nsite = len(g)
        idx = -np.ones((m, m, m), int)
        idx[g[:, 0], g[:, 1], g[:, 2]] = np.arange(self.nsite)
        self.idx = idx
        mc = 2 * per                                  # decoration cell edge
        tabc = -np.ones((mc, mc, mc), int)
        for k, s in enumerate(epi['sublattices']):
            for h in s['sites']:
                tabc[tuple(np.mod(h, mc))] = k
        self.sub = tabc[g[:, 0] % mc, g[:, 1] % mc, g[:, 2] % mc]
        self.tabc = tabc
        if np.any(self.sub < 0):
            raise ValueError('a graph site belongs to no sublattice')
        self.nsub = len(epi['sublattices'])
        self.off = np.cumsum([0] + [len(s['species']) for s in epi['sublattices']])   # global species index offsets
        self.nsp = int(self.off[-1])
        self.labels = [sp for s in epi['sublattices'] for sp in s['species']]
        # neighbour entries and energy tables: table 2p (owner on the first sublattice of pair p), 2p + 1 (transposed)
        ent_i, ent_j, ent_t = [], [], []
        self.bonds = []                               # per pair type: (i, j) with i on the first sublattice
        self.tables = np.zeros((2 * len(epi['pairs']), self.nsp, self.nsp))
        for p, pr in enumerate(epi['pairs']):
            k1, k2 = pr['sub']
            u = np.array(pr['u'], float)
            vec = pr['vectors']
            if k1 == k2:
                u = 0.5 * (u + u.T)                   # exact U is symmetric; the two halves come from different chains
            o1, o2 = self.off[k1], self.off[k2]
            self.tables[2 * p, o1:o1 + u.shape[0], o2:o2 + u.shape[1]] = u
            self.tables[2 * p + 1, o2:o2 + u.shape[1], o1:o1 + u.shape[0]] = u.T
            bi, bj = [], []
            pos = np.mod(g, mc)
            key = (pos[:, 0] * mc + pos[:, 1]) * mc + pos[:, 2]
            for kk in np.unique(key[self.sub == k1]):
                # the sites of one position of the decoration cell see the vector list, or its inverse (a sublattice
                # with two orientations, e.g. the 8c sites of a Heusler structure)
                i = np.nonzero((self.sub == k1) & (key == kk))[0]
                h = pos[i[0]]
                for sg in (1, -1):
                    q = np.mod(h + sg * vec, mc)
                    if np.all(tabc[q[:, 0], q[:, 1], q[:, 2]] == k2):
                        break
                else:
                    raise ValueError('pair type %d: the bonds of a site do not end on the partner sublattice' % p)
                for v in sg * vec:
                    q = np.mod(g[i] + v, m)
                    bi.append(i); bj.append(idx[q[:, 0], q[:, 1], q[:, 2]])
            bi = np.concatenate(bi); bj = np.concatenate(bj)
            if k1 == k2:                              # every bond must be listed by both of its sites
                a1 = np.sort(bi * self.nsite + bj); a2 = np.sort(bj * self.nsite + bi)
                if not np.array_equal(a1, a2):
                    raise ValueError('pair type %d: the bond list is not symmetric' % p)
            self.bonds.append((bi, bj))
            ent_i.append(bi); ent_j.append(bj); ent_t.append(np.full(len(bi), 2 * p))
            if k1 != k2:                              # the partner's entry of the same bond
                ent_i.append(bj); ent_j.append(bi); ent_t.append(np.full(len(bi), 2 * p + 1))
        ei = np.concatenate(ent_i); order = np.argsort(ei, kind='stable')
        self.nbr = np.ascontiguousarray(np.concatenate(ent_j)[order], np.int64)
        self.tab = np.ascontiguousarray(np.concatenate(ent_t)[order], np.int32)
        self.ptr = np.ascontiguousarray(np.concatenate([[0], np.cumsum(np.bincount(ei, minlength=self.nsite))]), np.int64)
        so = np.argsort(self.sub, kind='stable')
        self.sl = np.ascontiguousarray(so, np.int64)
        self.sp = np.ascontiguousarray(np.concatenate([[0], np.cumsum(np.bincount(self.sub, minlength=self.nsub))]), np.int64)
        self.sub32 = np.ascontiguousarray(self.sub, np.int32)
        self.lib = None

    # ------------------------------------------------------------------ configurations
    def random_config(self, rng):
        """Random occupation with the composition of every sublattice (largest-remainder rounding)."""
        spec = np.zeros(self.nsite, np.int8)
        for k, s in enumerate(self.epi['sublattices']):
            sites = np.nonzero(self.sub == k)[0]
            n = len(sites)
            cnt = np.floor(s['conc'] * n).astype(int)
            rem = s['conc'] * n - cnt
            for t in np.argsort(-rem)[:n - cnt.sum()]:
                cnt[t] += 1
            lab = np.repeat(np.arange(len(cnt)) + self.off[k], cnt)
            rng.shuffle(lab)
            spec[sites] = lab
        return spec

    def concentrations(self, spec):
        """Concentration of every global species on its sublattice."""
        x = np.zeros(self.nsp)
        for k in range(self.nsub):
            sel = spec[self.sub == k]
            x[self.off[k]:self.off[k + 1]] = np.bincount(sel - self.off[k], minlength=self.off[k + 1] - self.off[k]) / len(sel)
        return x

    # ------------------------------------------------------------------ Monte Carlo
    def _kernel(self):
        if self.lib is None:
            so = os.path.join(os.environ.get('SRO_BUILD', tempfile.gettempdir()), 'epi_mc.so')
            src = os.path.join(HERE, 'epi_mc.c')
            if not os.path.exists(so) or os.path.getmtime(so) < os.path.getmtime(src):
                subprocess.run(['gcc', '-O2', '-shared', '-fPIC', '-o', so, src, '-lm'], check=True)
            lib = ctypes.CDLL(so)
            vp = ctypes.c_void_p
            lib.epi_mc.restype = ctypes.c_long
            lib.epi_mc.argtypes = [ctypes.c_long, vp, vp, vp, vp, ctypes.c_int, vp, vp, vp, vp, ctypes.c_long, ctypes.c_uint64]
            lib.epi_energy.restype = ctypes.c_double
            lib.epi_energy.argtypes = [ctypes.c_long, vp, vp, vp, vp, ctypes.c_int, vp]
            self.lib = lib
        return self.lib

    def sweeps(self, spec, kt, nsweep, seed):
        """nsweep Monte Carlo sweeps at kT (Ry), in place.  Returns the acceptance ratio."""
        eps = np.ascontiguousarray(self.tables / kt)
        acc = self._kernel().epi_mc(self.nsite, spec.ctypes.data, self.ptr.ctypes.data, self.nbr.ctypes.data, self.tab.ctypes.data,
                                    self.nsp, eps.ctypes.data, self.sub32.ctypes.data, self.sp.ctypes.data, self.sl.ctypes.data,
                                    nsweep, seed)
        return acc / (nsweep * self.nsite)

    def energy(self, spec):
        """Pair energy per site (Ry)."""
        eps = np.ascontiguousarray(self.tables)
        return self._kernel().epi_energy(self.nsite, spec.ctypes.data, self.ptr.ctypes.data, self.nbr.ctypes.data,
                                         self.tab.ctypes.data, self.nsp, eps.ctypes.data) / self.nsite

    def pair_counts(self, spec):
        """Per pair type the bond counts n(c, s): c on the first sublattice of the type, s on the second."""
        out = []
        for p, pr in enumerate(self.epi['pairs']):
            k1, k2 = pr['sub']
            n1, n2 = self.off[k1 + 1] - self.off[k1], self.off[k2 + 1] - self.off[k2]
            bi, bj = self.bonds[p]
            c = (spec[bi] - self.off[k1]) * n2 + (spec[bj] - self.off[k2])
            out.append(np.bincount(c, minlength=n1 * n2).reshape(n1, n2).astype(float))
        return out

    def alpha_from_counts(self, counts, x):
        """Warren-Cowley tables from accumulated bond counts and the concentrations x (global species)."""
        out = []
        for p, pr in enumerate(self.epi['pairs']):
            k1, k2 = pr['sub']
            pc = counts[p] / counts[p].sum()
            x1, x2 = x[self.off[k1]:self.off[k1 + 1]], x[self.off[k2]:self.off[k2 + 1]]
            out.append(1.0 - pc / np.outer(x1, x2))
        return out

    def run(self, temp, nsweep=4000, nequil=1000, every=5, seed=1, spec=None):
        """Monte Carlo at temperature temp (K).  Returns a dictionary: alpha (list of tables), e (Ry per site), cv (k_B per
        site), acc, x, spec (last configuration)."""
        kt = KB_RY * temp
        rng = np.random.default_rng(seed)
        if spec is None:
            spec = self.random_config(rng)
        self.sweeps(spec, kt, nequil, seed * 7919 + 1)
        counts = [np.zeros_like(c) for c in self.pair_counts(spec)]
        es, acc, nb = [], 0.0, max(1, nsweep // every)
        for b in range(nb):
            acc += self.sweeps(spec, kt, every, seed * 7919 + 2 + b)
            for t, c in zip(counts, self.pair_counts(spec)):
                t += c
            es.append(self.energy(spec))
        es = np.array(es)
        x = self.concentrations(spec)
        return dict(alpha=self.alpha_from_counts(counts, x), e=es.mean(), cv=self.nsite * es.var() / kt**2, acc=acc / nb, x=x,
                    spec=spec, temp=temp)

    # ------------------------------------------------------------------ Krivoglaz-Clapp-Moss, one sublattice
    def kcm(self, temp, nk=24, spherical=True):
        """Pair correlations of the random-phase form on one sublattice.  S(q) = (A + U(q) / kT)^-1 for the independent
        species 2..n (species 1 eliminated); A = C^-1 with C = diag(x) - x x^T (Krivoglaz-Clapp-Moss), or, spherical = True,
        the matrix A that restores the on-site sum rule <S(q)> = C.  Returns alpha (list of tables), ok."""
        if self.nsub != 1:
            raise ValueError('kcm: one sublattice only')
        x = self.epi['sublattices'][0]['conc']
        n = len(x)
        beta = 1.0 / (KB_RY * temp)
        q = (np.stack(np.meshgrid(*[np.arange(nk)] * 3, indexing='ij'), -1).reshape(-1, 3) + 0.5) / nk
        uq = np.zeros((len(q), n - 1, n - 1))
        phase = []
        for pr in self.epi['pairs']:
            u = 0.5 * (pr['u'] + pr['u'].T)[1:, 1:]
            ph = np.cos(2.0 * np.pi * q @ pr['vectors'].T)          # (nq, nvec)
            phase.append(ph)
            uq += ph.sum(1)[:, None, None] * u
        c = np.diag(x[1:]) - np.outer(x[1:], x[1:])
        a = np.linalg.inv(c)
        ok = True
        if spherical:
            for it in range(500):
                mat = a + beta * uq
                if np.linalg.eigvalsh(mat).min() <= 0.0:
                    ok = False
                    break
                mean = np.linalg.inv(mat).mean(0)
                step = np.linalg.inv(c) - np.linalg.inv(mean)      # <S> too large: stiffen A
                a = a + 0.5 * step
                if np.abs(step).max() < 1e-12 * np.abs(a).max():
                    break
            else:
                ok = False
        mat = a + beta * uq
        if np.linalg.eigvalsh(mat).min() <= 0.0:
            ok = False
        s = np.linalg.inv(mat)
        out = []
        for ph in phase:
            sr = (s * ph[:, :1, None]).mean(0)                      # S(R) on the first vector of the orbit
            al = np.zeros((n, n))
            al[1:, 1:] = -sr / np.outer(x[1:], x[1:])
            al[0, 1:] = -(al[1:, 1:] * x[1:, None]).sum(0) / x[0]   # marginals: sum_c x_c alpha_cs = 0
            al[1:, 0] = al[0, 1:]
            al[0, 0] = -(x[1:] * al[0, 1:]).sum() / x[0]
            out.append(al)
        return out, ok

    def instability(self, nk=24):
        """Mean-field instability of the random alloy on one sublattice: the highest temperature (K) at which
        C^-1 + U(q) / kT has a zero eigenvalue, and its wave vector in units of 2 pi / alat.  T <= 0: stable at all T."""
        x = self.epi['sublattices'][0]['conc']
        q = np.stack(np.meshgrid(*[np.arange(2 * nk + 1)] * 3, indexing='ij'), -1).reshape(-1, 3) / (2.0 * nk)   # q in [0, 1]
        uq = np.zeros((len(q), len(x) - 1, len(x) - 1))
        for pr in self.epi['pairs']:
            u = 0.5 * (pr['u'] + pr['u'].T)[1:, 1:]
            uq += np.cos(2.0 * np.pi * q @ pr['vectors'].T).sum(1)[:, None, None] * u
        c = np.diag(x[1:]) - np.outer(x[1:], x[1:])
        w, v = np.linalg.eigh(c)
        ch = (v * np.sqrt(w)) @ v.T
        lam = np.linalg.eigvalsh(-ch @ uq @ ch)[:, -1]
        k = int(np.argmax(lam))
        return lam[k] / KB_RY, 2.0 * q[k]


def fit_marginals(p, x, y, iters=200):
    """Iterative proportional fitting of the joint table p to the marginals x (rows) and y (columns)."""
    p = np.maximum(np.array(p, float), 0.0)
    for _ in range(iters):
        p *= (x / p.sum(1))[:, None]
        p *= (y / p.sum(0))[None, :]
    return p


def sro_block(model, alphas, temp):
    """The [[sro]] input of asr for the first-shell pair types, with the input concentrations (marginals exact)."""
    lines = []
    subs = model.epi['sublattices']
    for pr, al in zip(model.epi['pairs'], alphas):
        if pr['shell'] != 1:
            continue
        k1, k2 = pr['sub']
        x, y = subs[k1]['conc'], subs[k2]['conc']
        p = fit_marginals(np.outer(x, y) * (1.0 - al), x, y)
        if k1 == k2:
            p = fit_marginals(0.5 * (p + p.T), x, y)
        a = 1.0 - p / np.outer(x, y)
        lines.append('[[sro]]   # %.0f K' % temp)
        if len(subs) > 1:
            lines.append('sublattices = ["%s", "%s"]' % (subs[k1]['name'], subs[k2]['name']))
        lines.append('shell = 1')
        ent = []
        for i, j in itertools.product(range(len(x)), range(len(y))):
            if k1 == k2 and j <= i:
                continue
            ent.append('["%s", "%s", %.10f]' % (subs[k1]['species'][i], subs[k2]['species'][j], a[i, j]))
        lines.append('pairs = [' + ', '.join(ent) + ']')
    return '\n'.join(lines)


def print_tables(model, alphas, head):
    subs = model.epi['sublattices']
    for pr, al in zip(model.epi['pairs'], alphas):
        k1, k2 = pr['sub']
        print('  %s  %s - %s  shell %d  (%d sites)' % (head, subs[k1]['name'], subs[k2]['name'], pr['shell'], len(pr['vectors'])))
        print('      ' + ''.join('%10s' % s for s in subs[k2]['species']))
        for i, c in enumerate(subs[k1]['species']):
            print('  %4s' % c + ''.join('%10.5f' % v for v in al[i]))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('epi', help='pair tables written by asr ([epi] file)')
    ap.add_argument('-T', '--temperatures', type=float, nargs='+', required=True, help='temperatures in K')
    ap.add_argument('--size', type=int, default=12, help='box edge in conventional cubes (default 12)')
    ap.add_argument('--sweeps', type=int, default=4000)
    ap.add_argument('--equil', type=int, default=1000)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--band-only', action='store_true', help='leave out the screened Coulomb part')
    ap.add_argument('--kcm', action='store_true', help='also the Krivoglaz-Clapp-Moss / spherical-model result (one sublattice)')
    ap.add_argument('--cool', action='store_true', help='start each temperature from the last configuration (give T descending)')
    a = ap.parse_args(argv)
    epi = load_epi(a.epi, a.band_only)
    if not epi['pairs']:
        sys.exit('no [[pair]] tables in ' + a.epi)
    mdl = Model(epi, a.size)
    print('pair model: %s, %d sites, %d pair types, species %s' % (epi['lattice'], mdl.nsite, len(epi['pairs']), ' '.join(mdl.labels)))
    for pr in epi['pairs']:
        k1, k2 = pr['sub']
        print('  %s - %s shell %d: U (mRy)' % (epi['sublattices'][k1]['name'], epi['sublattices'][k2]['name'], pr['shell']))
        for row in pr['u']:
            print('     ' + ''.join('%10.4f' % (1e3 * v) for v in row))
    if mdl.nsub == 1:
        ts, qs = mdl.instability()
        if ts > 0.0:
            print('mean-field instability of the random alloy: %.0f K at q = (%.3f, %.3f, %.3f) 2 pi / a' % (ts, *qs))
            print('  (q = 0: phase separation; zone boundary: ordering.  Mean field is an upper bound of the transition.)')
        else:
            print('mean-field: the random alloy is stable at every temperature')
    spec, blocks = None, []
    for t in a.temperatures:
        r = mdl.run(t, a.sweeps, a.equil, seed=a.seed, spec=spec if a.cool else None)
        spec = r['spec']
        print('\nT = %.0f K   E = %.5f mRy/site   C_v = %.3f k_B/site   acceptance %.3f' % (t, 1e3 * r['e'], r['cv'], r['acc']))
        print_tables(mdl, r['alpha'], 'Monte Carlo alpha')
        if a.kcm and mdl.nsub == 1:
            al, ok = mdl.kcm(t)
            if ok:
                print_tables(mdl, al, 'spherical model alpha')
            else:
                print('  spherical model: below its instability at this temperature')
        blocks.append(sro_block(mdl, r['alpha'], t))
    print('\n# first-shell [[sro]] input of asr (Monte Carlo)')
    print('\n\n'.join(blocks))


if __name__ == '__main__':
    main()
