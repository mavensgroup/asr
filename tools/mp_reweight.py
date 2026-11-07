#!/usr/bin/env python3
"""Exact Mookerjee-Prasad (MP) SRO averages from site-resolved supercell charges, by reweighting.

The first-shell MP model of the ASR [[sro]] route: the centre has species c; given c, its first-shell sites are independent with
X_s^(c) = x_s (1 - alpha_cs); every other site is random.  This measure is the random measure times the weight
    w(c; shell) = prod_{j in shell 1} X_{s_j}^(c) / x_{s_j},
so for any site quantity q of the centre
    E_MP[q | c] = E_rand[q w | c] / E_rand[w | c].
This is exact for independent (grand-canonical) occupations.  A cell with fixed composition narrows the distribution of the
neighbour count (hypergeometric), and the reweighted value is then about 6 % low at alpha = -0.148 (128 sites).  So the tool also
prints the binomial average of a quadratic fit q = g(n_unlike), which has no such bias.
Random cells with one potential per site (site_classes = true) therefore give the exact MP average without a new cell per alpha.
The hopping charge of the MP model is E_MP[q | c](alpha) - E_rand[q | c].  At identical (free-atom) potentials this is the quantity
the ASR [[sro]] run computes for component c, up to the ASR truncation.
The reweighting holds every cell at its own E_F, so the model cell is not neutral: dq_A + dq_B != 0 at x = 1/2.  The ASR fixes E_F by
neutrality.  The comparable number therefore adds the first-order E_F response, dq_A - f (dq_A + dq_B) with f = N_A / (N_A + N_B) at E_F.
Pass f as --dos-fraction; a cell of equal concentrations is assumed, and f comes from the k-space DOS at E_F of a random cell.

  python3 -I tools/mp_reweight.py random  --template in.toml --seed 3 --out cell/in.toml   (shuffled species, site_classes = true)
  python3 -I tools/mp_reweight.py siteify --template in.toml --out cell/in.toml            (same occupation, site_classes = true)
  python3 -I tools/mp_reweight.py analyse --runs r1 r2 ... --pair Cu,Zn --alpha -0.148 0.148 -0.266
  (each run directory holds in.toml and run.log with the per-class "class <Sym><n> ... qnet=" lines)
"""
import argparse, itertools, re, sys
import numpy as np


def read_cell(path):
    """alat (unused), cell vectors, site positions and species (cubic units) from a cell TOML written by tools/sqs.py."""
    txt = open(path).read()
    vec = re.search(r'^vectors\s*=\s*(\[.*\])\s*$', txt, re.M)
    if not vec: sys.exit(f'{path}: no vectors line')
    sup = np.array(eval(vec.group(1)), dtype=float)
    pos, spc = [], []
    for blk in txt.split('[[site]]')[1:]:
        p = re.search(r'position\s*=\s*(\[[^\]]*\])', blk); s = re.search(r'species\s*=\s*"([^"]+)"', blk)
        pos.append(np.array(eval(p.group(1)), dtype=float)); spc.append(s.group(1))
    return txt, sup, np.array(pos), spc


def first_shell(sup, pos):
    """List of first-shell neighbour indices (with periodic images, multiplicity kept) of every site."""
    n = len(pos)
    imgs = np.array([a * sup[0] + b * sup[1] + c * sup[2] for a, b, c in itertools.product(range(-2, 3), repeat=3)])
    dmin = np.inf; dl = []
    for i in range(n):
        for j in range(n):
            r = np.linalg.norm(pos[j] - pos[i] + imgs, axis=1)
            dl.append((i, j, r)); m = r[r > 1e-6].min() if np.any(r > 1e-6) else np.inf; dmin = min(dmin, m)
    nb = [[] for _ in range(n)]
    for i, j, r in dl:
        nb[i] += [j] * int(np.sum(np.abs(r - dmin) < 1e-5))
    return nb


def read_qnet(path, nsite):
    q = {}
    for line in open(path):
        m = re.match(r'\s*class\s+([A-Za-z]+)(\d+)\s+S\(Bohr\)=.*qnet=\s*(\S+)', line)
        if m: q[int(m.group(2))] = float(m.group(3))
    if len(q) != nsite: sys.exit(f'{path}: {len(q)} per-site classes, expected {nsite} (site_classes = true?)')
    return np.array([q[i + 1] for i in range(nsite)])


def write_cell(txt, spc, out, name):
    head, *blks = txt.split('[[site]]')
    head = re.sub(r'^site_classes\s*=.*$', 'site_classes = true', head, flags=re.M)
    head = re.sub(r'^name\s*=.*$', f'name = "{name}"', head, count=1, flags=re.M)
    nb = []
    for b, s in zip(blks, spc):
        nb.append(re.sub(r'species\s*=\s*"[^"]+"', f'species = "{s}"', b, count=1))
    import os; os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
    open(out, 'w').write(head + ''.join('[[site]]' + b for b in nb))


def analyse(a):
    sA, sB = a.pair.split(',')
    rows = []                                       # (run, species, q, n_unlike)
    bino_all = {}                                   # (species, alpha) -> binomial hopping charge at fixed E_F
    for d in a.runs:
        txt, sup, pos, spc = read_cell(f'{d}/in.toml')
        nb = first_shell(sup, pos)
        z = len(nb[0])
        q = read_qnet(f'{d}/run.log', len(pos))
        for i, s in enumerate(spc):
            rows.append((d, s, q[i], sum(spc[j] != s for j in nb[i]), z))
    for c in (sA, sB):
        u = sB if c == sA else sA
        R = [r for r in rows if r[1] == c]
        q = np.array([r[2] for r in R]); n = np.array([r[3] for r in R]); z = R[0][4]
        x = {sA: 0.5, sB: 0.5} if a.x is None else {sA: a.x, sB: 1.0 - a.x}
        print(f'centre {c}: {len(q)} sites in {len(a.runs)} cells, z = {z}, <q> = {q.mean():+.5f}, '
              f'<n_unlike> = {n.mean():.3f} (random {z * x[u]:.3f})')
        k, b = np.polyfit(n, q, 1)
        res = q - (k * n + b)
        print(f'  fit q = b + k n_unlike: k = {k:+.5f} e per unlike neighbour, rms {res.std():.4f}; '
              f'per-cell <q> spread {np.std([q[[r[0] == d for r in R]].mean() for d in a.runs]):.5f}')
        for al in a.alpha:
            xu, xc = x[u] * (1 - al), None
            xl = 1.0 - xu                                   # X_like^(c) (binary: diagonal follows)
            w = (xu / x[u]) ** n * (xl / (1.0 - x[u])) ** (z - n)
            mp = np.sum(q * w) / np.sum(w)
            rng = np.random.default_rng(1); bs = []
            for _ in range(2000):                           # bootstrap over cells, then sites (cell-to-cell correlation)
                cs = rng.choice(len(a.runs), len(a.runs))
                idx = np.concatenate([np.where([r[0] == a.runs[ci] for r in R])[0] for ci in cs])
                idx = rng.choice(idx, len(idx))
                bs.append(np.sum(q[idx] * w[idx]) / np.sum(w[idx]) - q[idx].mean())
            dn = z * (xu - x[u])
            from math import comb                           # binomial estimate from a quadratic fit g(n): no fixed-composition bias
            g = np.poly1d(np.polyfit(n, q, 2))
            bino = sum(comb(z, m) * (xu ** m * (1 - xu) ** (z - m) - x[u] ** m * (1 - x[u]) ** (z - m)) * g(m) for m in range(z + 1))
            bino_all[(c, al)] = bino
            print(f'  alpha {al:+.3f}: reweighted {mp - q.mean():+.5f} +- {np.std(bs):.5f} e (Neff {np.sum(w) ** 2 / np.sum(w * w):.0f}); '
                  f'binomial with quadratic g(n) {bino:+.5f}; linear k dn {k * dn:+.5f} (dn = {dn:+.3f})')
    if a.dos_fraction is not None and a.x is None:
        for al in a.alpha:
            qa, qb = bino_all[(sA, al)], bino_all[(sB, al)]
            print(f'neutral alpha {al:+.3f}: {sA} {qa - a.dos_fraction * (qa + qb):+.5f}, cell imbalance at fixed E_F {qa + qb:+.5f} '
                  f'(f = {a.dos_fraction:.3f})')


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('random'); r.add_argument('--template', required=True); r.add_argument('--seed', type=int, required=True)
    r.add_argument('--out', required=True)
    s = sub.add_parser('siteify'); s.add_argument('--template', required=True); s.add_argument('--out', required=True)
    s.add_argument('--name', default=None)
    n = sub.add_parser('analyse'); n.add_argument('--runs', nargs='+', required=True); n.add_argument('--pair', required=True)
    n.add_argument('--alpha', type=float, nargs='+', default=[-0.148]); n.add_argument('--x', type=float, default=None)
    n.add_argument('--dos-fraction', type=float, default=None, help='N_A / (N_A + N_B) at E_F: also print the neutral hopping charge')
    a = ap.parse_args()
    if a.cmd == 'analyse': return analyse(a)
    txt, sup, pos, spc = read_cell(a.template)
    if a.cmd == 'random':
        spc = list(np.random.default_rng(a.seed).permutation(spc)); name = f'random cell seed {a.seed}'
    else:
        name = a.name or 'site classes'
    write_cell(txt, spc, a.out, f'{name} ({len(pos)} sites)')
    print('wrote', a.out)


if __name__ == '__main__':
    main()
