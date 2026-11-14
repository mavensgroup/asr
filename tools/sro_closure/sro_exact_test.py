#!/usr/bin/env python3
"""Exact solution of the charge-excess model on supercells against the single-site F (docs/FINDINGS_LOG.md, "B3 closure").

Model: q_i = b_i - (V_i - mu)/a, V = M q (M = 2/R, periodic Ewald).  On a periodic supercell q(k) = G(k) b(k), G = a/(a + M(k)),
G = 0 at reciprocal lattice vectors (neutrality).  Species averages for ANY b:
    dq = Gamma db,   Gamma(p, p') = sum_R G(R) [P_R(p' | p) - x_p'] + G00 x_p'      (convention of the code: Gamma 1 = G00 1)
    V  = a (1 - Gamma) db,   F = a (Gamma^-1 - 1).
"exact": all lattice vectors R of the supercell, P_R measured in the configuration.  "n shells": the sum cut after shell n
(beyond it P_R = x).  n = 1 with the measured first-shell P is the single-site model of asr_screening (code: G00, G1 of the
infinite lattice on the offset mesh nk = 28; "sc": G00, G1 of the supercell itself).
gamma = eigenvalues of Gamma on the deviation subspace x.db = 0; F = a (1/gamma - 1).

The lattice (fcc or bcc) is embedded in a simple-cubic grid of spacing alat/2 with n = 2 L points per direction, so that the
solve and the pair correlations are FFTs.  Configurations: Monte Carlo (sro_mc.c) on a first-shell pair Hamiltonian, or built.

Usage: python3 -I sro_exact_test.py [L] [part ...]     parts: shells binary ternary domains bcc   (default: all; L = 24)
Environment: NSWEEP (default 2000), NPROC (default all cores; the Monte Carlo runs are independent), SRO_BUILD (dir of sro_mc.so).
"""
import ctypes
import itertools
import multiprocessing
import os
import subprocess
import sys
import tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                                    # python3 -I drops the script directory
from sro_pole_scan import coulomb_k, prim, lattice, solve   # noqa: E402
from sro_closure import kcm                                 # noqa: E402

ALAT = {'fcc': 6.78412, 'bcc': 5.38572}                     # Bohr: FeNiCr (3.59 A), CuZn (2.85 A)


def load_mc():
    so = os.path.join(os.environ.get('SRO_BUILD', tempfile.gettempdir()), 'sro_mc.so')
    src = os.path.join(HERE, 'sro_mc.c')
    if not os.path.exists(so) or os.path.getmtime(so) < os.path.getmtime(src):
        subprocess.run(['gcc', '-O2', '-shared', '-fPIC', '-o', so, src, '-lm'], check=True)
    lib = ctypes.CDLL(so)
    lib.mc.restype = ctypes.c_long
    lib.mc.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_long, ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p, ctypes.c_int,
                       ctypes.c_void_p, ctypes.c_long, ctypes.c_uint64, ctypes.c_int]
    return lib


class Cell:
    """Lattice lat in a periodic sc grid of n = 2 L points per direction (spacing alat/2)."""

    def __init__(self, lat, L):
        self.lat, self.L, self.n, self.alat = lat, L, 2 * L, ALAT[lat]
        n = self.n
        ix = np.indices((n, n, n))
        if lat == 'fcc':
            self.mask = (ix.sum(0) % 2 == 0)
            self.nb = np.array([d for d in itertools.product((-1, 0, 1), repeat=3) if sum(abs(c) for c in d) == 2], dtype=np.int32)
        else:
            self.mask = (ix[0] % 2 == ix[1] % 2) & (ix[1] % 2 == ix[2] % 2)
            self.nb = np.array(list(itertools.product((-1, 1), repeat=3)), dtype=np.int32)
        self.z = len(self.nb)
        self.sites = np.flatnonzero(self.mask.ravel()).astype(np.int64)
        self.nsite = self.sites.size
        w = (ix + L) % n - L                                  # wrapped grid vectors
        self.r2 = (w ** 2).sum(0)                             # |R|^2 in units (alat/2)^2
        # M(k) on the grid mesh through the cubic irreducible wedge
        key = np.sort(np.abs(w.reshape(3, -1)), axis=0)
        code = (key[0] * (L + 1) + key[1]) * (L + 1) + key[2]
        ucode, inv = np.unique(code, return_inverse=True)
        ukey = np.array([ucode // (L + 1) ** 2, (ucode // (L + 1)) % (L + 1), ucode % (L + 1)], dtype=float)
        p = prim(lat) * self.alat
        kv = ukey * (4 * np.pi / (n * self.alat))
        ph = (kv.T @ p) / (2 * np.pi)
        recip = np.all(np.abs(ph - np.round(ph)) < 1e-9, axis=1)   # reciprocal lattice vectors: G = 0 there
        m = np.full(ucode.size, np.inf)
        m[~recip] = coulomb_k(p, kv[:, ~recip])
        self.mk = m[inv].reshape(n, n, n)
        self.r1 = np.linalg.norm(self.nb[0]) * self.alat / 2
        # shells of the lattice (relative vectors of two lattice sites are lattice vectors)
        r2l = self.r2[self.mask]
        self.shell_r2 = np.unique(r2l)                        # shell_r2[0] = 0
        self.order = np.argsort(self.r2.ravel()[self.sites], kind='stable')
        self.r2sorted = self.r2.ravel()[self.sites][self.order]

    def green(self, a):
        """G(R) of the supercell (grid array; zero off the lattice)."""
        gk = a / (a + self.mk)
        return np.fft.ifftn(gk).real

    def random_config(self, x, rng):
        cnt = np.rint(np.array(x) * self.nsite).astype(int)
        cnt[0] += self.nsite - cnt.sum()
        s = np.repeat(np.arange(len(x)), cnt).astype(np.int8)
        rng.shuffle(s)
        spec = np.full(self.n ** 3, -1, dtype=np.int8)
        spec[self.sites] = s
        return spec

    def run_mc(self, lib, spec, eps, nsweep, seed, local=0):
        eps = np.ascontiguousarray(eps, dtype=np.float64)
        return lib.mc(self.n, spec.ctypes.data, self.nsite, self.sites.ctypes.data, self.z, self.nb.ctypes.data, eps.shape[0],
                      eps.ctypes.data, nsweep, seed, local)

    def correlations(self, spec):
        """x_p and P_R(p' | p) on the lattice vectors, sorted by |R| (array [nR, nsp, nsp])."""
        nsp = spec.max() + 1
        g = spec.reshape(self.n, self.n, self.n)
        nk = [np.fft.fftn((g == p).astype(float)) for p in range(nsp)]
        x = np.array([(spec == p).sum() for p in range(nsp)]) / self.nsite
        pr = np.empty((self.nsite, nsp, nsp))
        for p in range(nsp):
            for p2 in range(nsp):
                c = np.fft.ifftn(np.conj(nk[p]) * nk[p2]).real / self.nsite       # (1/N) sum_r n_p(r) n_p2(r + R)
                pr[:, p, p2] = c.ravel()[self.sites][self.order] / x[p]
        return x, pr

    def gammas(self, a, x, pr, nshell):
        """Gamma (code convention) cut after shell 1 .. nshell, and exact (all R).  Also z_s G_s and the shell-mean P_s."""
        gr = self.green(a).ravel()[self.sites][self.order]
        nsp = x.size
        dev = pr - x[None, None, :]
        term = gr[:, None, None] * dev                           # G(R) (P_R - x)
        g00 = gr[0]
        base = g00 * np.eye(nsp)
        out, zg, ps = [], [], []
        acc = np.zeros((nsp, nsp))
        for s in range(1, nshell + 1):
            sel = self.r2sorted == self.shell_r2[s]
            acc = acc + term[sel].sum(0)
            out.append(base + acc)
            zg.append(gr[sel].sum())
            ps.append(pr[sel].mean(0))
        exact = base + term[1:].sum(0)
        return out, exact, np.array(zg), np.array(ps), g00


_CELLS = {}


def _mc_job(job):
    """Worker: one Monte Carlo run; returns the configurations after the sweep counts of the schedule."""
    lat, L, x, eps, schedule, seed, local = job
    if (lat, L) not in _CELLS:
        _CELLS[(lat, L)] = Cell(lat, L)
    cell, lib = _CELLS[(lat, L)], load_mc()
    spec = cell.random_config(x, np.random.default_rng(seed))
    out, done = [], 0
    for tot in schedule:
        cell.run_mc(lib, spec, np.array(eps), tot - done, seed + 7919 * tot, local)
        done = tot
        out.append(spec.copy())
    return out


def run_jobs(jobs):
    """Monte Carlo runs in parallel (independent configurations); NPROC from the environment, default all cores."""
    nproc = min(len(jobs), int(os.environ.get('NPROC', os.cpu_count())))
    with multiprocessing.get_context('fork').Pool(nproc) as pool:
        return pool.map(_mc_job, jobs, chunksize=1)


class Closure:
    """KCM / spherical-model closure (sro_closure.py) on the infinite lattice, offset mesh nk: Gamma for a mode alpha_1."""

    def __init__(self, lat, nk=40):
        self.m, self.c1, _, self.z, _ = lattice(lat, ALAT[lat], nk)
        self.c1min = {'fcc': -1 / 3, 'bcc': -1.0}[lat]

    def gamma(self, a, alpha1):
        return kcm(a / (a + self.m), self.c1, alpha1, self.c1min)[0]

    def matrix(self, a, ahat, x):
        """Gamma_KCM for several components: the scalar closure on every eigenmode of Ahat (diag(x) Ahat is symmetric)."""
        sx = np.sqrt(x)
        y = sx[:, None] * ahat / sx[None, :]
        w, v = np.linalg.eigh(0.5 * (y + y.T))
        gm = np.array([self.gamma(a, lam / self.z) for lam in w])
        return (v * gm[None, :]) @ v.T * sx[None, :] / sx[:, None]


CLOS = {}


def dev_eigs(gam, x):
    """Eigenvalues of Gamma on the deviation subspace (diag(x) Gamma is symmetric; the trivial eigenvector is 1)."""
    sx = np.sqrt(x)
    y = sx[:, None] * gam / sx[None, :]
    y = 0.5 * (y + y.T)
    w, v = np.linalg.eigh(y)
    triv = np.argmax(np.abs(v.T @ sx))
    return np.delete(w, triv)


def f_of(a, gam):
    return a * (1.0 / gam - 1.0)


def fmt(v, w=8, d=3):
    return ' '.join(f'{t:+{w}.{d}f}' for t in np.atleast_1d(v))


def code_lattice(lat, betas):
    """a, G00, G1 of the code (infinite lattice, offset mesh nk = 28) for each beta."""
    m, c1, r1, z, amin = lattice(lat, ALAT[lat], 28)
    return {b: solve(m, c1, r1, amin, b)[:3] for b in betas}


def part_shells(cell, lib, code):
    print(f'\n=== shells: z_s G_s of the supercell ({cell.lat}, L = {cell.L}, {cell.nsite} sites) against beta ===')
    x = np.array([0.5, 0.5])
    spec = cell.random_config(x, np.random.default_rng(0))
    x, pr = cell.correlations(spec)
    ns = 12
    print('shell  R^2/(a/2)^2  z_s   ' + '   '.join(f'beta {b:4.2f}' for b in code))
    cols = []
    for b, (a, g00, g1) in code.items():
        _, _, zg, _, g00s = cell.gammas(a, x, pr, ns)
        cols.append((g00s, zg, g00, g1))
    zs = [int((cell.r2sorted == cell.shell_r2[s]).sum()) for s in range(ns + 1)]
    print(f'    0  {0:6d}      {1:3d}   ' + '   '.join(f'{c[0]:+9.4f}' for c in cols) + '    (G00, supercell)')
    for s in range(1, ns + 1):
        print(f'  {s:3d}  {int(cell.shell_r2[s]):6d}      {zs[s]:3d}   ' + '   '.join(f'{c[1][s - 1]:+9.4f}' for c in cols))
    print('  sum over shells >= 2    ' + '   '.join(f'{-c[0] - c[1][0]:+9.4f}' for c in cols) + '    (= -G00 - z G1; sum over all R is 0)')
    print('  code G00 (nk 28)        ' + '   '.join(f'{c[2]:+9.4f}' for c in cols))
    print('  code z G1               ' + '   '.join(f'{cell.z * c[3]:+9.4f}' for c in cols))


def analyse_binary(cell, spec, code, label):
    x, pr = cell.correlations(spec)
    rows = []
    for b, (a, g00, g1) in code.items():
        out, exact, zg, ps, g00s = cell.gammas(a, x, pr, 12)
        al = 1 - ps[:, 0, 1] / x[1]                              # alpha_s
        lam = cell.z * al[0]
        g_ss = g00 + g1 * lam                                    # the code's single-site value at the measured alpha_1
        gn = [dev_eigs(o, x)[0] for o in out]
        ge = dev_eigs(exact, x)[0]
        rows.append((b, a, al, lam, ge, g_ss, gn, CLOS[cell.lat].gamma(a, al[0])))
    b, a, al, lam, ge, g_ss, gn, gk = rows[0]
    s = f'{label:>9s} {al[0]:+6.3f} {al[1]:+6.3f} {al[2]:+6.3f} {al[5]:+6.3f} {al[11]:+6.3f} {lam:+6.2f} |'
    for b, a, al, lam, ge, g_ss, gn, gk in rows:
        s += f' {ge:+7.4f} {g_ss:+7.4f} {gk:+7.4f} {gn[1]:+7.4f} {gn[3]:+7.4f} {gn[11]:+7.4f} | {f_of(a, ge):+8.2f} {f_of(a, g_ss):+8.2f} {f_of(a, gk):+8.2f} |'
    print(s)


def head_binary(code):
    print('        K  alpha1 alpha2 alpha3 alpha6 alph12 lambda |' + ''.join(
        f' beta {b:4.2f}: gamma exact, single-site (code), KCM, 2 / 4 / 12 shells   | F exact, code, KCM (Ry/e)  |' for b in code))


def part_binary(cell, lib, code, nsweep):
    for xb in (0.5, 0.2):
        print(f'\n=== binary {cell.lat} A(1-x)B(x), x = {xb}, L = {cell.L} ({cell.nsite} sites): MC with unlike-pair energy K kT, {nsweep} sweeps ===')
        print('    (fcc Ising critical point K_c = 0.204 at x = 0.5; K > K_c: phase separation, domains coarsen with the sweeps)')
        head_binary(code)
        ks = (0.0, 0.05, 0.10, 0.15, 0.18, 0.20, 0.22, 0.25, 0.30, 0.40, 0.60)
        snaps = (2, 5, 10, 20, 50, 100, 300, 1000, 3000)
        jobs = [(cell.lat, cell.L, [1 - xb, xb], K * (1 - np.eye(2)), (nsweep,), 1234 + int(K * 1000), 0) for K in ks]
        if xb == 0.5:
            jobs.append((cell.lat, cell.L, [1 - xb, xb], 0.40 * (1 - np.eye(2)), snaps, 77, 1))
        res = run_jobs(jobs)
        for K, r in zip(ks, res):
            analyse_binary(cell, r[0], code, f'{K:5.2f}')
        if xb == 0.5:
            print('  coarsening at K = 0.40 (same run, snapshots; local exchange, so domains grow slowly):')
            for tot, spec in zip(snaps, res[-1]):
                analyse_binary(cell, spec, code, f'sw {tot}')


def part_domains(code):
    print('\n=== built configurations, fcc binary x = 0.5: fully segregated slabs (100), thickness L/2 cells each ===')
    head_binary(code)
    for L in (4, 8, 16, 24, 32):
        cell = Cell('fcc', L)
        g = np.full((cell.n,) * 3, -1, dtype=np.int8)
        ix = np.indices((cell.n,) * 3)[0]
        g[cell.mask] = (ix[cell.mask] >= cell.n // 2).astype(np.int8)
        analyse_binary(cell, g.ravel(), code, f'L {L}')
    print('=== built: periodic array of cubic B precipitates of edge D cells in A (x = 1/8), cell L = 2 D ===')
    head_binary(code)
    for D in (2, 4, 8, 12, 16):
        cell = Cell('fcc', 2 * D)
        ix = np.indices((cell.n,) * 3)
        inside = np.all(ix < 2 * D, axis=0)
        g = np.full((cell.n,) * 3, -1, dtype=np.int8)
        g[cell.mask] = inside[cell.mask].astype(np.int8)
        analyse_binary(cell, g.ravel(), code, f'D {D}')


def part_ternary(cell, lib, code, nsweep):
    x0 = [0.66, 0.14, 0.20]
    names = ['Fe', 'Ni', 'Cr']
    dq_rand = np.array([0.0, 0.146, -0.190])                    # random-alloy net charges of the triad table (Ni, Cr)
    dq_rand[0] = -(x0[1] * dq_rand[1] + x0[2] * dq_rand[2]) / x0[0]
    print(f'\n=== ternary fcc Fe66Ni14Cr20, L = {cell.L} ({cell.nsite} sites): MC, pair energy K kT on the unlike pairs of the triad, {nsweep} sweeps ===')
    print('    alpha: measured first shell (FeNi, FeCr, NiCr).  lam: nontrivial eigenvalues of the measured Ahat.')
    print('    F: nontrivial eigenvalues (Ry/e).  dq: species charges for the b that gives the random-alloy charges '
          f'({fmt(dq_rand, 6, 3)}); model-internal.')
    lines = {b: [] for b in code}
    jobs, tags = [], []
    for tri in ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)):
        eps = np.zeros((3, 3))
        eps[0, 1] = eps[1, 0] = tri[0]; eps[0, 2] = eps[2, 0] = tri[1]; eps[1, 2] = eps[2, 1] = tri[2]
        for K in (0.1, 0.2, 0.3, 0.5, 1.0):
            jobs.append((cell.lat, cell.L, x0, K * eps, (nsweep,), 4321 + int(K * 1000), 0))
            tags.append((tri, K))
    res = run_jobs(jobs)
    for (tri, K), r in zip(tags, res):
        if True:
            spec = r[0]
            x, pr = cell.correlations(spec)
            for b, (a, g00, g1) in code.items():
                db = dq_rand / g00
                out, exact, zg, ps, g00s = cell.gammas(a, x, pr, 4)
                p1 = ps[0]
                al = [1 - p1[0, 1] / x[1], 1 - p1[0, 2] / x[2], 1 - p1[1, 2] / x[2]]
                ahat = cell.z * (p1 - x[None, :])
                lam = dev_eigs(ahat, x)
                g_ss = g00 * np.eye(3) + g1 * ahat               # the code's matrix at the measured first-shell P
                g_k = CLOS[cell.lat].matrix(a, ahat, x)
                fe, fs, f2, f4 = (f_of(a, dev_eigs(g, x)) for g in (exact, g_ss, g_k, out[3]))
                tname = ''.join(str(t) for t in tri)
                lines[b].append(f'  {tname}  {K:4.2f} | {al[0]:+6.3f} {al[1]:+6.3f} {al[2]:+6.3f} | {lam[0]:+5.2f} {lam[1]:+6.2f} | {fmt(fe)} | {fmt(fs)} | {fmt(f2)} | {fmt(f4)}'
                                f' | {fmt(exact @ db, 6, 3)} | {fmt(g_ss @ db, 6, 3)} | {fmt(g_k @ db, 6, 3)}')
    for b, (a, g00, g1) in code.items():
        print(f'--- beta {b:4.2f}: a {a:.4f}, G00 {g00:.4f}, G1 {g1:+.4f}, pole lambda* {-g00 / g1:.2f}')
        print('  triad   K | alpha FeNi  FeCr   NiCr | lam          | F exact           | F single-site     | F KCM             | F 4 shells'
              '        | dq exact (Fe Ni Cr)  | dq single-site       | dq KCM')
        print('\n'.join(lines[b]))


def part_bcc(lib, L, nsweep):
    cell = Cell('bcc', L)
    CLOS['bcc'] = Closure('bcc')
    code = code_lattice('bcc', (1.0,))
    print(f'\n=== binary bcc A50B50 (CuZn lattice), L = {cell.L} ({cell.nsite} sites): ordering, K < 0 (B2 transition at K = -0.315) ===')
    head_binary(code)
    ks = (0.0, -0.05, -0.10, -0.15, -0.20, -0.25, -0.30, 0.05, 0.10, 0.15)
    res = run_jobs([(cell.lat, cell.L, [0.5, 0.5], K * (1 - np.eye(2)), (nsweep,), 999 + int(abs(K) * 1000), 0) for K in ks])
    for K, r in zip(ks, res):
        analyse_binary(cell, r[0], code, f'{K:5.2f}')


def main():
    args = sys.argv[1:]
    L = int(args.pop(0)) if args and args[0].isdigit() else 24
    parts = args or ['shells', 'binary', 'ternary', 'domains', 'bcc']
    lib = load_mc()
    nsweep = int(os.environ.get('NSWEEP', 2000))
    code = code_lattice('fcc', (1.0, 1.2))
    CLOS['fcc'] = Closure('fcc')
    for b, (a, g00, g1) in code.items():
        print(f'fcc code lattice functions, beta {b}: a {a:.5f}  G00 {g00:.5f}  G1 {g1:+.5f}  pole lambda* {-g00 / g1:.3f}')
    cell = Cell('fcc', L) if set(parts) & {'shells', 'binary', 'ternary'} else None
    if 'shells' in parts:
        part_shells(cell, lib, code_lattice('fcc', (0.8, 1.0, 1.02, 1.1, 1.2)))
    if 'binary' in parts:
        part_binary(cell, lib, code, nsweep)
    if 'ternary' in parts:
        part_ternary(cell, lib, code, nsweep)
    if 'domains' in parts:
        part_domains(code)
    if 'bcc' in parts:
        part_bcc(lib, L, nsweep)


if __name__ == '__main__':
    main()
