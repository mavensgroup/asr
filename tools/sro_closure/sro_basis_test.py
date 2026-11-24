#!/usr/bin/env python3
"""Exact supercell solution of the basis-resolved charge-excess model against the branch-resolved SRO closure.

Tests multi-site sublattices:
  1. bcc folded into a 2-site B2 cell (must match 1-site primitive bcc on a commensurate mesh).
  2. hcp (2 sites per cell, 1 sublattice, z = 12 at ideal c/a = sqrt(8/3), z = 6 at c/a = 1.7).
  3. NiAs (4 sites per cell, 2 sublattices; ordered cations with db = 0, disordered anions on sublattice 2,
     z = 6 out-of-plane at c/a = 1.6 and z = 6 in-plane at c/a = 1.7).

On an L x L x L periodic supercell of the primitive cell, M_super is block-diagonalised by discrete Fourier
transform over the L^3 commensurate wavevectors k (with M(k = 0) projected onto the neutral cell subspace).
Canonical Metropolis Monte Carlo generates binary configurations on the active sublattice with pair coupling K kT
on the SRO shell, from ordering (K < 0) through random (K = 0) to phase separation (K > 0).
"""
from math import erfc, sqrt, pi
import numpy as np


def coulomb_k_basis(p, tau, kv, eta=None):
    """Ewald lattice Coulomb matrix M_bb'(k) in Ry/e (tau gauge). At k = 0 omits G = 0 and shifts uniform mode."""
    nb = tau.shape[1]
    om = abs(np.linalg.det(p))
    b = 2 * pi * np.linalg.inv(p).T
    et = (sqrt(pi) / (om / nb) ** (1 / 3)) if eta is None else eta
    dtmax = 0.0
    for b1 in range(nb):
        for b2 in range(b1 + 1, nb):
            dtmax = max(dtmax, np.linalg.norm(tau[:, b2] - tau[:, b1]))
    rmax = 7.0 / et
    gmax = 14.0 * et + np.linalg.norm(b, axis=0).max()
    nr = [int(np.ceil((rmax + dtmax) * np.linalg.norm(b[:, i]) / (2 * pi))) + 1 for i in range(3)]
    idx = np.array(np.meshgrid(*[np.arange(-k, k + 1) for k in nr], indexing='ij')).reshape(3, -1)
    tv = p @ idx
    tv = tv[:, np.linalg.norm(tv, axis=0) <= rmax + dtmax]
    ng = [int(np.ceil(gmax * np.linalg.norm(p[:, i]) / (2 * pi))) + 1 for i in range(3)]
    idx = np.array(np.meshgrid(*[np.arange(-k, k + 1) for k in ng], indexing='ij')).reshape(3, -1)
    gv = b @ idx
    gv = gv[:, np.linalg.norm(gv, axis=0) <= gmax]

    nk = kv.shape[1]
    mb = np.zeros((nb, nb, nk), dtype=complex)
    for b1 in range(nb):
        for b2 in range(b1, nb):
            d = tau[:, b2] - tau[:, b1]
            g_phase = np.exp(-1j * (gv.T @ d))
            vr = tv + d[:, None]
            rr = np.linalg.norm(vr, axis=0)
            keep = (rr > 1e-10) & (rr <= rmax)
            vr = vr[:, keep]
            erfc_r = np.array([erfc(et * x) / x for x in rr[keep]])
            for ik in range(nk):
                k = kv[:, ik]
                q = k[:, None] + gv
                q2 = (q * q).sum(0)
                if np.any(q2 < 1e-14):
                    sel = q2 >= 1e-14
                    s = 4 * pi / om * np.sum(np.exp(-q2[sel] / (4 * et ** 2)) / q2[sel] * g_phase[sel]) - pi / (om * et ** 2)
                else:
                    s = 4 * pi / om * np.sum(np.exp(-q2 / (4 * et ** 2)) / q2 * g_phase)
                sr = np.sum(erfc_r * np.exp(1j * (k @ vr)))
                val = s + sr
                if b1 == b2:
                    val = val.real - 2 * et / sqrt(pi)
                mb[b1, b2, ik] = 2.0 * val
                if b2 > b1:
                    mb[b2, b1, ik] = np.conj(mb[b1, b2, ik])
    for ik in range(nk):
        if np.linalg.norm(kv[:, ik]) < 1e-12:
            mb[:, :, ik] += 1e12 / nb
    return mb


def proto(name, coa=None):
    if name == 'bcc':
        avec = np.array([[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]).T
        frac = np.zeros((3, 1))
        label = np.array([1])
    elif name == 'b2_bcc':
        avec = np.eye(3)
        frac = np.array([[0.0, 0.5], [0.0, 0.5], [0.0, 0.5]])
        label = np.array([1, 1])
    elif name == 'hcp':
        ca = sqrt(8.0 / 3.0) if coa is None else coa
        avec = np.array([[1.0, -0.5, 0.0], [0.0, sqrt(3.0) / 2.0, 0.0], [0.0, 0.0, ca]])
        frac = np.array([[0.0, 1.0 / 3.0], [0.0, 2.0 / 3.0], [0.0, 0.5]])
        label = np.array([1, 1])
    elif name == 'nias':
        avec = np.array([[1.0, -0.5, 0.0], [0.0, sqrt(3.0) / 2.0, 0.0], [0.0, 0.0, coa]])
        frac = np.array([[0.0, 0.0, 1.0 / 3.0, 2.0 / 3.0],
                         [0.0, 0.0, 2.0 / 3.0, 1.0 / 3.0],
                         [0.0, 0.5, 0.25, 0.75]])
        label = np.array([1, 1, 2, 2])
    return avec, frac, avec @ frac, label


def sro_shell(p, tau, label, ksro):
    nb = tau.shape[1]
    sites = [b for b in range(nb) if (ksro == 0 or label[b] == ksro)]
    idx = np.array(np.meshgrid(*[np.arange(-2, 3)] * 3, indexing='ij')).reshape(3, -1)
    tv = p @ idx
    rmin = np.inf
    for b1 in sites:
        for b2 in sites:
            v = tv + (tau[:, b2] - tau[:, b1])[:, None]
            r = np.linalg.norm(v, axis=0)
            rmin = min(rmin, r[r > 1e-10].min())
    bonds = []
    for b1 in sites:
        for b2 in sites:
            v = tv + (tau[:, b2] - tau[:, b1])[:, None]
            r = np.linalg.norm(v, axis=0)
            for j in np.flatnonzero(np.abs(r - rmin) < 1e-7 * rmin):
                bonds.append((b1, b2, idx[:, j], v[:, j]))
    return sites, bonds, rmin


def zone_special_k(b):
    t1, h = 1.0 / 3.0, 0.5
    pts = [[0.0, 0.0, 0.0]]
    for i1 in (0, 1):
        for i2 in (0, 1):
            for i3 in (0, 1):
                if i1 or i2 or i3:
                    pts.append([i1 * h, i2 * h, i3 * h])
    pts.append([t1, t1, 0.0])
    pts.append([t1, t1, h])
    pts.append([0.5, 0.25, 0.0])
    return b @ np.array(pts).T


def branch_model(avec, tau0, label, alat, beta, nk=14, ksro=0):
    p, tau = avec * alat, tau0 * alat
    nb = tau.shape[1]
    b = 2 * pi * np.linalg.inv(p).T
    sites, bonds, r1_sh = sro_shell(p, tau, label, ksro)
    _, _, r1_nn = sro_shell(p, tau, label, 0)
    nb_sro = len(sites)
    z = len(bonds) // nb_sro
    smap = {s: i for i, s in enumerate(sites)}

    f = (np.array(np.meshgrid(*[np.arange(nk)] * 3, indexing='ij')).reshape(3, -1) + 0.5) / nk
    kv = b @ f
    mb = coulomb_k_basis(p, tau, kv)
    ksp = zone_special_k(b)
    msp = coulomb_k_basis(p, tau, ksp)

    meig = np.empty((nb, kv.shape[1]))
    u_sro = np.empty((nb_sro, nb, kv.shape[1]), dtype=complex)
    c1_br = np.empty((nb_sro, kv.shape[1]))
    v_br = np.empty((nb_sro, nb_sro, kv.shape[1]), dtype=complex)
    amin = 0.0
    for ik in range(kv.shape[1]):
        ev, u = np.linalg.eigh(mb[:, :, ik])
        meig[:, ik] = ev
        u_sro[:, :, ik] = u[sites, :]
        amin = max(amin, -ev.min())
        ck = np.zeros((nb_sro, nb_sro), dtype=complex)
        for (b1, b2, _, rv) in bonds:
            ck[smap[b1], smap[b2]] += np.exp(1j * (kv[:, ik] @ rv)) / z
        ec, vc = np.linalg.eigh(ck)
        c1_br[:, ik] = ec
        v_br[:, :, ik] = vc

    sp_ev = []
    sp_u = []
    sp_c1 = []
    sp_v = []
    for ik in range(ksp.shape[1]):
        ev, u = np.linalg.eigh(msp[:, :, ik])
        amin = max(amin, -ev.min())
        sp_ev.append(ev)
        sp_u.append(u[sites, :])
        ck = np.zeros((nb_sro, nb_sro), dtype=complex)
        for (b1, b2, _, rv) in bonds:
            ck[smap[b1], smap[b2]] += np.exp(1j * (ksp[:, ik] @ rv)) / z
        ec, vc = np.linalg.eigh(ck)
        sp_c1.append(ec)
        sp_v.append(vc)

    w0_eig = (np.abs(u_sro) ** 2).sum(0) / nb_sro  # [nb, nk]

    def beta_of(a):
        g = np.sum(w0_eig * (a / (a + meig))) / kv.shape[1]
        return a * (g - 1.0) * r1_nn / (2.0 * g)

    lo, hi = np.log(amin * 1e-9), np.log(1e9)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if beta_of(amin + np.exp(mid)) > beta:
            lo = mid
        else:
            hi = mid
    a = amin + np.exp(0.5 * (lo + hi))

    # Branch Green functions g_m(k) = v_m^H G_sro(k) v_m = sum_n |v_m^H u_sro(:,n)|^2 a / (a + ev_n)
    gk_br = np.empty((nb_sro, kv.shape[1]))
    for ik in range(kv.shape[1]):
        ov = np.abs(v_br[:, :, ik].conj().T @ u_sro[:, :, ik]) ** 2  # [nb_sro, nb]
        gk_br[:, ik] = ov @ (a / (a + meig[:, ik]))

    c1_flat = c1_br.ravel()
    gk_flat = gk_br.ravel()
    g00 = gk_flat.mean()
    g1 = (gk_flat * c1_flat).mean()

    # Special points: c1_ord and g_ord (largest mode G among ordering modes with minimal c1)
    c1_ord = min(c1_flat.min(), min(ec.min() for ec in sp_c1))
    g_ord_cands = []
    for ik in range(ksp.shape[1]):
        ov = np.abs(sp_v[ik].conj().T @ sp_u[ik]) ** 2
        g_sp = ov @ (a / (a + sp_ev[ik]))
        for m in range(nb_sro):
            if sp_c1[ik][m] < c1_ord + 1e-10:
                g_ord_cands.append(g_sp[m])
    for ik in range(kv.shape[1]):
        for m in range(nb_sro):
            if c1_br[m, ik] < c1_ord + 1e-10:
                g_ord_cands.append(gk_br[m, ik])
    # Ordering point: highest G (lowest effective M) among modes attaining c1_ord, matching screen_lattice
    g_ord = max(g_ord_cands)

    return dict(a=a, g00=g00, g1=g1, z=z, r1_nn=r1_nn, r1_sh=r1_sh,
                c1k=c1_flat, gk=gk_flat, c1_ord=c1_ord, g_ord=g_ord)


def clip_gamma(mdl, lam):
    c1, g, z, c1_ord, g_ord = mdl['c1k'], mdl['gk'], mdl['z'], mdl['c1_ord'], mdl['g_ord']
    al = lam / z
    if -1.0 <= lam <= -1.0 / c1_ord:
        return mdl['g00'] + mdl['g1'] * lam
    if lam > 0:
        c0_cl = c1[c1 < c1.max() - 1e-12].max()
        def pt(c0):
            w = np.maximum(0.0, c1 - c0)
            return np.sum(c1 * w) / np.sum(w), np.sum(g * w) / np.sum(w)
        acl, gcl = pt(c0_cl)
        if al >= acl:
            m = min(1.0, (al - acl) / (1.0 - acl))
            return (1.0 - m) * gcl
        lo, hi = c1_ord, c0_cl
    else:
        c0_ord = c1[c1 > c1.min() + 1e-12].min()
        def pt(c0):
            w = np.maximum(0.0, c0 - c1)
            return np.sum(c1 * w) / np.sum(w), np.sum(g * w) / np.sum(w)
        aord, gord_c = pt(c0_ord)
        if al <= aord:
            m = (aord - al) / max(aord - c1_ord, 1e-12)
            return (1.0 - m) * gord_c + m * g_ord
        lo, hi = c0_ord, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if pt(mid)[0] < al:
            lo = mid
        else:
            hi = mid
    return pt(0.5 * (lo + hi))[1]


class BasisSupercell:
    """Exact charge-excess solver and Metropolis MC on an L x L x L periodic supercell of any structure."""

    def __init__(self, name, coa, alat, L, ksro, a):
        avec, _, tau0, label = proto(name, coa)
        self.p, self.tau = avec * alat, tau0 * alat
        self.nb = self.tau.shape[1]
        self.L = L
        self.sites, bonds, _ = sro_shell(self.p, self.tau, label, ksro)
        self.nb_sro = len(self.sites)
        self.z = len(bonds) // self.nb_sro
        self.nsite_sro = self.nb_sro * (L ** 3)
        self.a = a
        smap = {s: i for i, s in enumerate(self.sites)}

        # Neighbour table on the SRO sublattice [nsite_sro, z]
        self.nbr = np.empty((self.nsite_sro, self.z), dtype=np.int32)
        per_b = {i: [] for i in range(self.nb_sro)}
        for (b1, b2, dn, _) in bonds:
            per_b[smap[b1]].append((smap[b2], dn))
        for i1 in range(L):
            for i2 in range(L):
                for i3 in range(L):
                    cell_id = (i1 * L + i2) * L + i3
                    for r1 in range(self.nb_sro):
                        u = cell_id * self.nb_sro + r1
                        for k, (r2, dn) in enumerate(per_b[r1]):
                            j1, j2, j3 = (i1 + dn[0]) % L, (i2 + dn[1]) % L, (i3 + dn[2]) % L
                            v = ((j1 * L + j2) * L + j3) * self.nb_sro + r2
                            self.nbr[u, k] = v

        # Commensurate k-mesh of the L x L x L supercell
        f = np.array(np.meshgrid(*[np.arange(L)] * 3, indexing='ij')).reshape(3, -1) / L
        b = 2 * pi * np.linalg.inv(self.p).T
        kv = b @ f
        mb = coulomb_k_basis(self.p, self.tau, kv)
        # Note: coulomb_k_basis is in the tau gauge exp(i k . (T + tau_b' - tau_b)).
        # Convert to periodic cell gauge M_tilde_bb'(k) = exp(-i k . (tau_b' - tau_b)) M_bb'(k) for FFT over cell indices!
        self.gsro_k = np.empty((L, L, L, self.nb_sro, self.nb_sro), dtype=complex)
        for ik in range(kv.shape[1]):
            k = kv[:, ik]
            ph = np.exp(-1j * (self.tau.T @ k))  # exp(-i k . tau_b)
            # M_tilde = diag(exp(i k.tau)) M diag(exp(-i k.tau))
            mtilde = (np.conj(ph)[:, None] * mb[:, :, ik]) * ph[None, :]
            gtilde = a * np.linalg.inv(a * np.eye(self.nb) + mtilde)
            i1, i2, i3 = int(round(f[0, ik] * L)), int(round(f[1, ik] * L)), int(round(f[2, ik] * L))
            self.gsro_k[i1, i2, i3] = gtilde[np.ix_(self.sites, self.sites)]

    def exact_gamma_binary(self, occ):
        """Exact binary Gamma eigenvalue on the supercell for occupation occ in {0, 1}^nsite_sro."""
        x1 = occ.mean()
        x0 = 1.0 - x1
        # Drive db = (+x0, -x1) / sqrt(x0 x1) has zero mean and unit x-variance: <db^2> = x1 x0^2/(x0 x1) + x0 x1^2/(x0 x1) = 1
        db = np.where(occ == 1, sqrt(x0 / x1), -sqrt(x1 / x0))
        db_grid = db.reshape(self.L, self.L, self.L, self.nb_sro)
        db_k = np.fft.fftn(db_grid, axes=(0, 1, 2))
        dq_k = np.einsum('ijkrs,ijks->ijkr', self.gsro_k, db_k)
        dq = np.fft.ifftn(dq_k, axes=(0, 1, 2)).real.ravel()
        # Project onto species subspace: Gamma_mode = <db * dq>
        gam_ex = np.mean(db * dq)
        # First-shell Warren-Cowley alpha_1: P(0 | 1) = fraction of 0 neighbours around 1 sites
        n1_nbr0 = np.sum((occ[:, None] == 1) & (occ[self.nbr] == 0)) / (self.nsite_sro * x1 * self.z)
        alpha1 = 1.0 - n1_nbr0 / x0
        return alpha1, gam_ex

    def mc_binary(self, x1, K, nsweep, seed):
        rng = np.random.default_rng(seed)
        n1 = int(round(x1 * self.nsite_sro))
        occ = np.zeros(self.nsite_sro, dtype=np.int8)
        occ[:n1] = 1
        rng.shuffle(occ)
        # E = K * (number of unlike first-shell pairs); K > 0 clustering, K < 0 ordering
        for _ in range(nsweep):
            iu = rng.integers(0, self.nsite_sro, size=self.nsite_sro)
            iv = rng.integers(0, self.nsite_sro, size=self.nsite_sro)
            ru = rng.random(self.nsite_sro)
            for t in range(self.nsite_sro):
                i, j = iu[t], iv[t]
                si, sj = occ[i], occ[j]
                if si == sj:
                    continue
                # Count unlike change: flipping i from si to sj changes unlike count at neighbours a != j
                nb_i = self.nbr[i]
                nb_j = self.nbr[j]
                # For binary {0, 1}: unlike(s, a) = (s != a). Change when s -> 1-s is (1 - 2*(s != a))
                de_i = np.sum(np.where(nb_i != j, np.where(occ[nb_i] != sj, 1, -1), 0))
                de_j = np.sum(np.where(nb_j != i, np.where(occ[nb_j] != si, 1, -1), 0))
                de = K * (de_i + de_j)
                if de <= 0.0 or ru[t] < np.exp(-de):
                    occ[i], occ[j] = sj, si
        return occ


def check_monotonic(mdl, label):
    z, c1_ord = mdl['z'], mdl['c1_ord']
    lams = np.linspace(z * c1_ord, z, 401)
    gams = [clip_gamma(mdl, l) for l in lams]
    lams_cl = np.linspace(-1.0, z, 301)
    gams_cl = [clip_gamma(mdl, l) for l in lams_cl]
    pos = all(g > 0 for g in gams[:-1]) and abs(gams[-1]) < 1e-14
    dec_all = all(gams[i + 1] <= gams[i] + 1e-12 for i in range(len(gams) - 1))
    dec_cl = all(gams_cl[i + 1] <= gams_cl[i] + 1e-12 for i in range(len(gams_cl) - 1))
    print(f"  {label}: c1_ord={c1_ord:+.4f}, g_ord={mdl['g_ord']:.4f}, G00={mdl['g00']:.4f}, z*G1={z*mdl['g1']:+.4f}, "
          f"pos={pos}, dec_all={dec_all}, dec_[-1,z]={dec_cl}, Gamma(z*c1_ord)={gams[0]:.4f}")


def main():
    alat = 2.85 / 0.529177210903
    print("=== 1. Reduction check: 2-site folded B2 bcc vs 1-site primitive bcc ===")
    avec1, _, tau1, lab1 = proto('bcc')
    avec2, _, tau2, lab2 = proto('b2_bcc')
    m1 = branch_model(avec1, tau1, lab1, alat, beta=1.0, nk=16, ksro=0)
    m2 = branch_model(avec2, tau2, lab2, alat, beta=1.0, nk=8, ksro=0)
    print(f"  1-site bcc (nk=16): a={m1['a']:.6f}, G00={m1['g00']:.6f}, G1={m1['g1']:+.6f}, c1_ord={m1['c1_ord']:+.4f}, g_ord={m1['g_ord']:.6f}")
    print(f"  2-site B2  (nk=8) : a={m2['a']:.6f}, G00={m2['g00']:.6f}, G1={m2['g1']:+.6f}, c1_ord={m2['c1_ord']:+.4f}, g_ord={m2['g_ord']:.6f}")
    check_monotonic(m2, "2-site B2 bcc beta=1.0")

    print("\n=== 2. Positivity and monotonicity on hcp and NiAs ===")
    for name, coa, ksro, beta in [('hcp', None, 0, 1.0), ('hcp', None, 0, 1.2),
                                  ('hcp', 1.7, 0, 1.0),
                                  ('nias', 1.6, 2, 1.0), ('nias', 1.6, 2, 1.2),
                                  ('nias', 1.7, 2, 1.0)]:
        avec, _, tau0, lab = proto(name, coa)
        mdl = branch_model(avec, tau0, lab, alat, beta=beta, nk=12, ksro=ksro)
        check_monotonic(mdl, f"{name} coa={coa} ksro={ksro} beta={beta}")

    print("\n=== 3. Exact supercell MC vs 1-shell linear vs branch-clip closure ===")
    for name, coa, ksro, beta, L, ks in [
        ('hcp', None, 0, 1.15, 8, (-0.50, -0.25, -0.10, 0.0, 0.10, 0.20, 0.30, 0.50)),
        ('nias', 1.6, 2, 1.15, 8, (-0.50, -0.25, -0.10, 0.0, 0.15, 0.30, 0.45, 0.70)),
    ]:
        avec, _, tau0, lab = proto(name, coa)
        mdl = branch_model(avec, tau0, lab, alat, beta=beta, nk=12, ksro=ksro)
        sc = BasisSupercell(name, coa, alat, L=L, ksro=ksro, a=mdl['a'])
        print(f"\n  {name} (coa={coa}, ksro={ksro}, beta={beta}, L={L}, {sc.nsite_sro} active sites): "
              f"a={mdl['a']:.4f}, G00={mdl['g00']:.4f}, z*G1={mdl['z']*mdl['g1']:+.4f}")
        print("      K   alpha_1   lambda |  Gamma_ex  Gamma_lin  Gamma_clip |   F_ex     F_lin    F_clip")
        for K in ks:
            occ = sc.mc_binary(0.5, K, nsweep=80, seed=100 + int((K + 2) * 100))
            al1, gex = sc.exact_gamma_binary(occ)
            lam = mdl['z'] * al1
            glin = mdl['g00'] + mdl['g1'] * lam
            gclp = clip_gamma(mdl, lam)
            fex = mdl['a'] * (1.0 / gex - 1.0)
            flin = mdl['a'] * (1.0 / glin - 1.0) if glin > 1e-6 else float('nan')
            fclp = mdl['a'] * (1.0 / gclp - 1.0)
            print(f"  {K:+5.2f}   {al1:+7.4f}  {lam:+7.3f} |   {gex:7.4f}    {glin:+7.4f}     {gclp:7.4f} | {fex:+7.3f}  {flin:+8.3f}   {fclp:+7.3f}")


if __name__ == '__main__':
    main()
