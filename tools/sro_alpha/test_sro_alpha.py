#!/usr/bin/env python3
"""Self-test of sro_alpha.py (about a minute).  Run: python3 -I tools/sro_alpha/test_sro_alpha.py

 (1) bcc binary, first-shell interaction V: the exchange model is the Ising model with |J| = V / 4, k T_c = 1.5887 V
     (bcc Ising, k T_c / J = 6.3547).  The Monte Carlo is ordered below T_c and disordered above it.
 (2) High temperature: Monte Carlo == spherical model == first order -x_p x_q ... (alpha_1 = -x (1 - x) V / kT, binary).
 (3) fcc ternary on two shells: marginals of the Warren-Cowley tables, symmetry, Monte Carlo == spherical model at high T,
     and the mean-field instability against a direct evaluation.
 (4) B2-decorated bcc with an interaction between the two sublattices: marginals, and the sign of the response.
 (5) The file format of asr_driver::write_epi and the [[sro]] block.
"""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sro_alpha as sa   # noqa: E402

nfail = 0


def check(cond, name, info=''):
    global nfail
    print('  %s  %s %s' % ('ok  ' if cond else 'FAIL', name, info))
    if not cond:
        nfail += 1


def orbit(v):
    """All images of an integer vector under the 48 cubic operations."""
    out = set()
    for perm in ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)):
        for sg in np.ndindex(2, 2, 2):
            out.add(tuple(int((1 - 2 * sg[i]) * v[perm[i]]) for i in range(3)))
    return np.array(sorted(out))


def one_sub(lattice, species, conc, pairs):
    sites = [[0, 0, 0]] if lattice != 'fcc' else [[0, 0, 0]]
    if lattice == 'bcc':
        sites = [[0, 0, 0], [1, 1, 1]]
    if lattice == 'fcc':
        sites = [[0, 0, 0], [1, 1, 0], [1, 0, 1], [0, 1, 1]]
    return dict(lattice=lattice, alat=5.4, period=1,
                sublattices=[dict(name='a', species=species, conc=np.array(conc, float), sites=np.array(sites))],
                pairs=[dict(sub=(0, 0), shell=k + 1, vectors=orbit(v), u=np.array(u, float)) for k, (v, u) in enumerate(pairs)])


print('test_sro_alpha')
V = 2.0e-3                                             # Ry
tc = 1.5887 * V / sa.KB_RY                              # K

# ---- (1) bcc Ising transition
epi = one_sub('bcc', ['A', 'B'], [0.5, 0.5], [((1, 1, 1), [[0, 0], [0, V]])])
m = sa.Model(epi, 10)
lo = m.run(0.8 * tc, nsweep=1500, nequil=1500, seed=3)
hi = m.run(1.3 * tc, nsweep=1500, nequil=800, seed=3)
print('    bcc binary, T_c = %.0f K: alpha_1(0.8 T_c) = %.4f, alpha_1(1.3 T_c) = %.4f' % (tc, lo['alpha'][0][0, 1], hi['alpha'][0][0, 1]))
check(lo['alpha'][0][0, 1] < -0.55, 'ordered below the exact Ising T_c')
check(-0.25 < hi['alpha'][0][0, 1] < -0.08, 'short-range order only above T_c')
ts, qs = m.instability()
check(abs(ts / (2.0 * V / sa.KB_RY) - 1.0) < 1e-6 and np.allclose(sorted(qs), [0, 0, 1]),
      'mean-field instability: k T = 2 V at the H point', '(%.0f K, q = %s)' % (ts, qs))

# ---- (2) high temperature
t = 12.0 * V / sa.KB_RY
r = m.run(t, nsweep=3000, nequil=300, seed=5)
al, ok = m.kcm(t)
first = -0.25 * V / (sa.KB_RY * t)
print('    kT = 12 V: alpha_1 Monte Carlo %.5f, spherical model %.5f, first order %.5f' % (r['alpha'][0][0, 1], al[0][0, 1], first))
check(ok and abs(al[0][0, 1] - first) < 0.03 * abs(first), 'spherical model -> first order at high T')
check(abs(r['alpha'][0][0, 1] - al[0][0, 1]) < 1.5e-3, 'Monte Carlo == spherical model at high T')
check(abs(r['alpha'][0][0, 0] + r['alpha'][0][0, 1]) < 1e-12, 'binary marginal: alpha_AA = -alpha_AB at 50:50')

# ---- (3) fcc ternary, two shells
u1 = np.zeros((3, 3)); u1[1, 1] = 1.5e-3; u1[2, 2] = -1.0e-3; u1[1, 2] = u1[2, 1] = 0.4e-3
u2 = np.zeros((3, 3)); u2[1, 1] = -0.5e-3; u2[2, 2] = 0.3e-3; u2[1, 2] = u2[2, 1] = 0.1e-3
x = np.array([0.5, 0.3, 0.2])
epi3 = one_sub('fcc', ['A', 'B', 'C'], x, [((1, 1, 0), u1), ((2, 0, 0), u2)])
m3 = sa.Model(epi3, 8)
check(m3.nsite == 4 * 8**3 and len(m3.bonds[0][0]) == 12 * m3.nsite and len(m3.bonds[1][0]) == 6 * m3.nsite, 'fcc box: 12 and 6 bonds per site')
t = 3000.0
r = m3.run(t, nsweep=3000, nequil=300, seed=7)
al, ok = m3.kcm(t)
for k in range(2):
    a = r['alpha'][k]
    check(np.abs(a @ r['x']).max() < 1e-12 and np.abs(a - a.T).max() < 1e-12, 'shell %d: Monte Carlo marginals vanish, table symmetric' % (k + 1))
    check(np.abs(al[k] @ x).max() < 1e-12, 'shell %d: spherical-model marginals vanish' % (k + 1))
    d = np.abs(a - al[k]).max()
    check(ok and d < 4e-3, 'shell %d: Monte Carlo == spherical model at 3000 K' % (k + 1), '(max difference %.1e, largest alpha %.3f)' % (d, np.abs(al[k]).max()))
check(np.sign(r['alpha'][0][1, 1]) > 0 and np.sign(r['alpha'][0][2, 2]) < 0, 'U(B,B) > 0 makes B-B pairs rare, U(C,C) < 0 makes C-C pairs frequent')
e0 = sum(0.5 * len(p['vectors']) * (x @ (0.5 * (p['u'] + p['u'].T)) @ x) for p in epi3['pairs'])
rr = m3.run(1.0e7, nsweep=400, nequil=50, seed=9)
check(abs(rr['e'] - e0) < 2e-6, 'energy per site of the random alloy = sum_n (z_n / 2) x U x', '(%.3e vs %.3e Ry)' % (rr['e'], e0))
ts, qs = m3.instability(nk=8)
c = np.diag(x[1:]) - np.outer(x[1:], x[1:])
lam = []
for q in np.stack(np.meshgrid(*[np.arange(17)] * 3, indexing='ij'), -1).reshape(-1, 3) / 16.0:
    uq = sum(np.cos(2 * np.pi * p['vectors'] @ q).sum() * p['u'][1:, 1:] for p in epi3['pairs'])
    lam.append(np.linalg.eigvals(-c @ uq).real.max())
check(abs(ts - max(lam) / sa.KB_RY) < 1e-6 * abs(ts), 'instability temperature == largest eigenvalue of -C U(q)')

# ---- (4) two sublattices (B2 decoration of bcc), interaction between them
ub = np.array([[0.0, 0.0], [0.0, 1.0e-3]])
epi2 = dict(lattice='bcc', alat=5.4, period=1,
            sublattices=[dict(name='a', species=['A', 'B'], conc=np.array([0.6, 0.4]), sites=np.array([[0, 0, 0]])),
                         dict(name='b', species=['C', 'D'], conc=np.array([0.5, 0.5]), sites=np.array([[1, 1, 1]]))],
            pairs=[dict(sub=(0, 1), shell=1, vectors=orbit((1, 1, 1)), u=ub)])
m2 = sa.Model(epi2, 10)
r = m2.run(600.0, nsweep=2000, nequil=500, seed=11)
a = r['alpha'][0]
check(np.abs(a @ r['x'][2:]).max() < 1e-12 and np.abs(r['x'][:2] @ a).max() < 1e-12, 'two sublattices: both marginals vanish')
check(a[1, 1] > 0.02 and a[0, 0] > 0.0 and a[0, 1] < 0.0, 'U(B, D) > 0: B-D pairs rare, A-D pairs frequent', '(alpha_BD = %.4f)' % a[1, 1])
blk = sa.sro_block(m2, r['alpha'], 600.0)
check('sublattices = ["a", "b"]' in blk and blk.count('"]') == 1 and blk.count(', "') >= 4, '[[sro]] block for two sublattices lists all four pairs')

# ---- (5) file format
txt = '''name = "t"
lattice = "fcc"
alat = 6.7
period = 1
fermi_energy = -0.07
converged = true

[[sublattice]]
name = "fcc"
site_fraction = 1.0
species = ["Cu", "Ni"]
concentrations = [0.5, 0.5]
sites = [[0, 0, 0], [0, 1, 1], [1, 0, 1], [1, 1, 0]]

[[pair]]
sublattices = ["fcc", "fcc"]
shell = 1
vector = [1, 1, 0]
distance = 4.7376
multiplicity = 12
kernel = 0.2
vectors = [%s]
u_band = [
  [0.0, 0.0],
  [0.0, 1.0e-3],
]
u_es = [
  [0.0, 0.0],
  [0.0, 2.0e-4],
]
''' % ', '.join('[%d, %d, %d]' % tuple(v) for v in orbit((1, 1, 0)))
with tempfile.NamedTemporaryFile('w', suffix='.toml', delete=False) as f:
    f.write(txt)
e = sa.load_epi(f.name)
eb = sa.load_epi(f.name, band_only=True)
os.unlink(f.name)
check(abs(e['pairs'][0]['u'][1, 1] - 1.2e-3) < 1e-15 and abs(eb['pairs'][0]['u'][1, 1] - 1.0e-3) < 1e-15, 'epi.toml: band + screened Coulomb, or band only')
mf = sa.Model(e, 6)
r = mf.run(900.0, nsweep=600, nequil=200, seed=13)
blk = sa.sro_block(mf, r['alpha'], 900.0)
check(blk.startswith('[[sro]]') and 'pairs = [["Cu", "Ni", ' in blk and 'sublattices' not in blk, '[[sro]] block for one sublattice: the unlike pair only')
check(r['alpha'][0][0, 1] < 0.0, 'V(Cu, Ni) > 0: unlike neighbours favoured (alpha < 0)')

print('failed %d' % nfail)
sys.exit(1 if nfail else 0)
