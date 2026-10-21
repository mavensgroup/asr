#!/usr/bin/env python3
"""Special quasirandom structures (SQS) for a random alloy of N >= 2 species on a cubic Bravais lattice (bcc / fcc).

Searches supercells (Hermite normal forms of determinant N of the primitive cell, compact shapes only) and, in each, the
occupation by simulated annealing (swaps of two sites with different species) that best matches the random alloy on the first
`npair` pair shells: the pair probabilities
    P_ss'(k) = (bonds s-s' in shell k) / (bonds in shell k) = x_s x_s'      (all species pairs),
objective sum_k w_k sum_ss' |P_ss'(k) - x_s x_s'|.  For a binary this is exactly |Pi_k - (2x - 1)^2| with the spin correlation
Pi_k = <sigma_i sigma_j>_k (sigma = +1 for A, -1 for B), and the smallest triplet (Pi_3 = (2x - 1)^3) is matched as well; the
binary search is the one of the binary-only version (same random stream, same SQS for the same seed).  Output: one TOML fragment per SQS with
[structure] lattice = "cell" (vectors in units of the cubic lattice constant) and the [[site]] list; add [[species]], [scf],
[kspace] and [calculation] method = "kspace" (an explicit cell runs in k-space only).

  python3 -I tools/sqs.py --lattice bcc --n 16 --x 0.5 --species Fe Co --nsqs 2 --out sqs_dir
  python3 -I tools/sqs.py --lattice fcc --n 40 --species Co Cr Fe Mn Ni --x 0.2 0.2 0.2 0.2 0.2 --out hea   (--x: all
  concentrations, or only x of the first species for a binary)

Partial (sublattice) disorder: the SQS is built on the disordered sublattice only and every site is decorated with the ordered
atoms at fixed offsets (cubic lattice-constant units), e.g. L2_1 Co2(Fe,Cr)Al with the Fe/Cr 4b sublattice as the fcc SQS lattice:
  python3 -I tools/sqs.py --lattice fcc --n 10 --x 0.4 --species Fe Cr --decorate Al:0.5,0.5,0.5 Co:0.25,0.25,0.25 Co:-0.25,-0.25,-0.25 --out d
"""
import argparse, itertools, os, sys
import numpy as np

PRIM = {'bcc': np.array([[-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]]),
        'fcc': np.array([[0., .5, .5], [.5, 0., .5], [.5, .5, 0.]])}


def hnfs(n):
    """All 3x3 lower-triangular Hermite normal forms with determinant n."""
    for a in range(1, n + 1):
        if n % a: continue
        for c in range(1, n // a + 1):
            if (n // a) % c: continue
            f = n // (a * c)
            for b in range(c):
                for d in range(f):
                    for e in range(f):
                        yield np.array([[a, 0, 0], [b, c, 0], [d, e, f]])


def reduce_cell(v):
    """Simple LLL-like reduction of the rows of v (shorter, more orthogonal vectors)."""
    v = v.copy()
    changed = True
    while changed:
        changed = False
        for i, j in itertools.permutations(range(3), 2):
            mu = round(np.dot(v[i], v[j]) / np.dot(v[j], v[j]))
            if mu:
                w = v[i] - mu * v[j]
                if np.dot(w, w) < np.dot(v[i], v[i]) - 1e-9:
                    v[i] = w; changed = True
    return v


def sites_of(sup, prim):
    """Lattice points of the supercell (rows of sup are the supercell vectors, Cartesian)."""
    inv = np.linalg.inv(sup)
    n = int(round(abs(np.linalg.det(sup) / np.linalg.det(prim))))
    rng = range(-n, n + 1)
    pts = []
    for i, j, k in itertools.product(rng, rng, rng):
        r = i * prim[0] + j * prim[1] + k * prim[2]
        f = r @ inv
        if np.all(f > -1e-9) and np.all(f < 1 - 1e-9):
            pts.append(r)
        if len(pts) == n: break
    assert len(pts) == n, (len(pts), n)
    return np.array(pts)


def shells_and_pairs(sup, pos, nshell):
    """Pair lists per shell with periodic images: list over shells of (i, j, multiplicity) arrays."""
    n = len(pos)
    imgs = np.array([a * sup[0] + b * sup[1] + c * sup[2] for a, b, c in itertools.product(range(-3, 4), repeat=3)])
    d = {}
    for i in range(n):
        for j in range(n):
            dv = pos[j] - pos[i] + imgs
            for r in np.round(np.linalg.norm(dv, axis=1), 6):
                if r > 1e-6: d.setdefault((i, j), []).append(r)
    radii = sorted({r for v in d.values() for r in v})[:nshell]
    shells = []
    for r0 in radii:
        m = np.zeros((n, n))
        for (i, j), rs in d.items():
            m[i, j] = sum(1 for r in rs if abs(r - r0) < 1e-5)
        shells.append(m)
    return radii, shells


def triplet_list(sup, pos, radii):
    """Smallest triplet: the triangles of least perimeter with all sides in the first two shells, (R1, R1, R2) on bcc,
    (R1, R1, R1) on fcc -> list of (i, j, k) with images."""
    n = len(pos)
    imgs = [a * sup[0] + b * sup[1] + c * sup[2] for a, b, c in itertools.product(range(-2, 3), repeat=3)]
    nb = []
    for i in range(n):
        lst = []
        for j in range(n):
            for t in imgs:
                v = pos[j] - pos[i] + t
                r = np.linalg.norm(v)
                if abs(r - radii[0]) < 1e-5 or abs(r - radii[1]) < 1e-5: lst.append((j, v, r))
        nb.append(lst)
    cand = []
    for i in range(n):
        for (j, vj, rj), (k, vk, rk) in itertools.combinations(nb[i], 2):
            rjk = np.linalg.norm(vk - vj)
            if min(abs(rjk - radii[0]), abs(rjk - radii[1])) < 1e-5: cand.append((rj + rk + rjk, (i, j, k)))
    pmin = min(p for p, _ in cand)
    return np.array([t for p, t in cand if p < pmin + 1e-5])


def correlations(sig, shells, tri):
    pis = [sig @ m @ sig / m.sum() for m in shells]
    pt = np.mean(sig[tri[:, 0]] * sig[tri[:, 1]] * sig[tri[:, 2]]) if len(tri) else 0.0
    return np.array(pis), pt


def pair_probs(lab, shells, ns):
    """P[k, s, s'] = bond fraction s-s' in shell k (lab: species index per site)."""
    occ = np.eye(ns)[lab]
    return np.array([occ.T @ m @ occ / m.sum() for m in shells])


def anneal(counts, shells, tri, x, w, rng, steps=20000):
    """counts: sites per species; x: concentrations.  Returns (objective, species index per site)."""
    ns = len(counts)
    lab = np.concatenate([np.full(c, s) for s, c in enumerate(counts)]); rng.shuffle(lab)
    xx = np.outer(x, x)
    if ns == 2:
        t2, t3 = (2 * x[0] - 1) ** 2, (2 * x[0] - 1) ** 3
        def obj(l):
            p, t = correlations(1.0 - 2.0 * l, shells, tri)
            return np.sum(w * np.abs(p - t2)) + np.abs(t - t3)
    else:
        def obj(l):
            return np.sum(w * np.abs(pair_probs(l, shells, ns) - xx).sum(axis=(1, 2)))
    cur = obj(lab); best, blab = cur, lab.copy()
    for it in range(steps):
        T = 0.3 * (1 - it / steps) + 1e-4
        if ns == 2:
            i = rng.choice(np.where(lab == 0)[0]); j = rng.choice(np.where(lab == 1)[0])
        else:
            i = rng.integers(len(lab)); j = rng.choice(np.where(lab != lab[i])[0])
        lab[i], lab[j] = lab[j], lab[i]
        new = obj(lab)
        if new <= cur or rng.random() < np.exp(-(new - cur) / T):
            cur = new
            if cur < best - 1e-12: best, blab = cur, lab.copy()
        else:
            lab[i], lab[j] = lab[j], lab[i]
    return best, blab


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--lattice', default='bcc'); ap.add_argument('--n', type=int, default=16)
    ap.add_argument('--x', type=float, nargs='+', default=[0.5]); ap.add_argument('--species', nargs='+', required=True)
    ap.add_argument('--npair', type=int, default=6); ap.add_argument('--ncells', type=int, default=40)
    ap.add_argument('--nsqs', type=int, default=2); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', required=True)
    ap.add_argument('--decorate', nargs='*', default=[], help='symbol:dx,dy,dz atoms added at every SQS site (cubic units)')
    a = ap.parse_args()
    deco = []
    for d in a.decorate:
        sym, off = d.split(':'); deco.append((sym, np.array([float(c) for c in off.split(',')])))
    prim = PRIM[a.lattice]; rng = np.random.default_rng(a.seed)
    ns = len(a.species)
    if ns < 2: sys.exit('at least two species')
    x = np.array(a.x if len(a.x) == ns else ([a.x[0], 1.0 - a.x[0]] if ns == 2 and len(a.x) == 1 else []))
    if len(x) != ns or abs(x.sum() - 1.0) > 1e-9 or np.any(x <= 0): sys.exit('--x: one positive concentration per species, sum 1')
    counts = np.rint(x * a.n).astype(int)
    if np.any(np.abs(counts - x * a.n) > 1e-9): sys.exit('x_s * n must be integers')
    w = 1.0 / np.arange(1, a.npair + 1)                      # nearer shells weigh more
    cands = []
    for h in hnfs(a.n):
        sup = reduce_cell(h @ prim)
        lens = np.linalg.norm(sup, axis=1)
        cands.append((lens.max() / lens.min(), sup))
    cands.sort(key=lambda c: c[0])
    seen, results = [], []
    for aspect, sup in cands:
        key = tuple(np.round(sorted(np.linalg.norm(sup, axis=1)), 4))
        if key in seen: continue
        seen.append(key)
        if len(seen) > a.ncells: break
        pos = sites_of(sup, prim)
        radii, shells = shells_and_pairs(sup, pos, a.npair)
        tri = triplet_list(sup, pos, radii)
        best, lab = anneal(counts, shells, tri, x, w, rng)
        dp = np.abs(pair_probs(lab, shells, ns) - np.outer(x, x)).max(axis=(1, 2))   # worst pair deviation per shell
        pt = correlations(1.0 - 2.0 * lab, shells, tri) if ns == 2 else None
        results.append((best, sup, pos, lab, dp, pt, aspect))
        print(f'cell {len(seen):3d} aspect {aspect:5.2f}  objective {best:.4f}  max |P - x x| per shell {np.round(dp, 3)}', flush=True)
    results.sort(key=lambda r: r[0])
    os.makedirs(a.out, exist_ok=True)
    for q, (best, sup, pos, lab, dp, pt, aspect) in enumerate(results[:a.nsqs], 1):
        fn = os.path.join(a.out, f'sqs{a.n}_{q}.toml')
        with open(fn, 'w') as f:
            f.write(f'# SQS {a.lattice} N = {a.n}, ' + ' '.join(f'{s}{c:g}' for s, c in zip(a.species, x)) +
                    f': objective {best:.4f}, max |P_ss\' - x_s x_s\'| per shell (1..{a.npair}) {np.round(dp, 4).tolist()}\n')
            if ns == 2:
                p, t = pt
                f.write(f'# pair correlations {np.round(p, 4).tolist()}, triplet {t:+.4f} '
                        f'(random: {(2 * x[0] - 1) ** 2:.3f} / {(2 * x[0] - 1) ** 3:+.3f})\n')
            if deco: f.write('# decorated: ' + ' '.join(a.decorate) + '\n')
            f.write('# fragment: add [[species]], [scf], [kspace] and [calculation] method = "kspace"\n')
            f.write('[structure]\nlattice = "cell"\nsite_classes = true\nvectors = [' +
                    ', '.join('[' + ', '.join(f'{c:.6f}' for c in v) + ']' for v in sup) + ']\n\n')
            for r, s in zip(pos, lab):
                f.write('[[site]]\nposition = [' + ', '.join(f'{c:.6f}' for c in r) + f']\nspecies = "{a.species[s]}"\n\n')
            for sym, off in deco:
                for r in pos:
                    f.write('[[site]]\nposition = [' + ', '.join(f'{c:.6f}' for c in r + off) + f']\nspecies = "{sym}"\n\n')
        print('wrote', fn, f'objective {best:.4f}')


if __name__ == '__main__':
    main()
