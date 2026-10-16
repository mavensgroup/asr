#!/usr/bin/env python3
"""Special quasirandom structures (SQS) for a binary A_x B_(1-x) on a cubic Bravais lattice (bcc / fcc).

Searches supercells (Hermite normal forms of determinant N of the primitive cell, compact shapes only) and, in each, the
occupation by simulated annealing (pair swaps) that best matches the random-alloy correlation functions
    Pi_k = <sigma_i sigma_j>_{shell k} = (2x - 1)^2      (sigma = +1 for A, -1 for B)
for the first `npair` pair shells and the smallest triplet (Pi_3 = (2x - 1)^3).  Output: one TOML fragment per SQS with
[structure] lattice = "cell" (vectors in units of the cubic lattice constant) and the [[site]] list for `asr` kind = "lmto".

  python3 -I tools/sqs.py --lattice bcc --n 16 --x 0.5 --species Fe Co --nsqs 2 --out sqs_dir
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
    """Smallest triplet: triangles with sides (R1, R1, R2) (bcc) or (R1, R1, R1) (fcc) -> list of (i, j, k) with images."""
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
    tri = []
    sides = sorted([radii[0], radii[0], radii[1]])
    for i in range(n):
        for (j, vj, rj), (k, vk, rk) in itertools.combinations(nb[i], 2):
            s = sorted([rj, rk, np.linalg.norm(vk - vj)])
            if np.allclose(s, sides, atol=1e-5): tri.append((i, j, k))
    return np.array(tri)


def correlations(sig, shells, tri):
    pis = [sig @ m @ sig / m.sum() for m in shells]
    pt = np.mean(sig[tri[:, 0]] * sig[tri[:, 1]] * sig[tri[:, 2]]) if len(tri) else 0.0
    return np.array(pis), pt


def anneal(n_a, shells, tri, target2, target3, w, rng, steps=20000):
    n = shells[0].shape[0]
    sig = np.array([1.0] * n_a + [-1.0] * (n - n_a)); rng.shuffle(sig)
    def obj(s):
        p, t = correlations(s, shells, tri)
        return np.sum(w * np.abs(p - target2)) + np.abs(t - target3)
    cur = obj(sig); best, bsig = cur, sig.copy()
    for it in range(steps):
        T = 0.3 * (1 - it / steps) + 1e-4
        i = rng.choice(np.where(sig > 0)[0]); j = rng.choice(np.where(sig < 0)[0])
        sig[i], sig[j] = sig[j], sig[i]
        new = obj(sig)
        if new <= cur or rng.random() < np.exp(-(new - cur) / T):
            cur = new
            if cur < best - 1e-12: best, bsig = cur, sig.copy()
        else:
            sig[i], sig[j] = sig[j], sig[i]
    return best, bsig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--lattice', default='bcc'); ap.add_argument('--n', type=int, default=16)
    ap.add_argument('--x', type=float, default=0.5); ap.add_argument('--species', nargs=2, required=True)
    ap.add_argument('--npair', type=int, default=6); ap.add_argument('--ncells', type=int, default=40)
    ap.add_argument('--nsqs', type=int, default=2); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    prim = PRIM[a.lattice]; rng = np.random.default_rng(a.seed)
    n_a = int(round(a.x * a.n))
    if abs(n_a - a.x * a.n) > 1e-9: sys.exit('x * n must be an integer')
    t2, t3 = (2 * a.x - 1) ** 2, (2 * a.x - 1) ** 3
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
        best, sig = anneal(n_a, shells, tri, t2, t3, w, rng)
        p, t = correlations(sig, shells, tri)
        results.append((best, sup, pos, sig, p, t, aspect))
        print(f'cell {len(seen):3d} aspect {aspect:5.2f}  objective {best:.4f}  pairs {np.round(p, 3)}  triplet {t:+.3f}', flush=True)
    results.sort(key=lambda r: r[0])
    os.makedirs(a.out, exist_ok=True)
    for q, (best, sup, pos, sig, p, t, aspect) in enumerate(results[:a.nsqs], 1):
        fn = os.path.join(a.out, f'sqs{a.n}_{q}.toml')
        with open(fn, 'w') as f:
            f.write(f'# SQS {a.lattice} N = {a.n}, x({a.species[0]}) = {a.x}: objective {best:.4f}, pair correlations (shells 1..{a.npair}) '
                    f'{np.round(p, 4).tolist()}, triplet {t:+.4f} (random: {t2:.3f} / {t3:+.3f})\n')
            f.write('[structure]\nlattice = "cell"\nsite_classes = true\nvectors = [' +
                    ', '.join('[' + ', '.join(f'{c:.6f}' for c in v) + ']' for v in sup) + ']\n\n')
            for r, s in zip(pos, sig):
                f.write('[[site]]\nposition = [' + ', '.join(f'{c:.6f}' for c in r) + f']\nspecies = "{a.species[0] if s > 0 else a.species[1]}"\n\n')
        print('wrote', fn, f'objective {best:.4f}')


if __name__ == '__main__':
    main()
