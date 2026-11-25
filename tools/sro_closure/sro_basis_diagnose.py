#!/usr/bin/env python3
"""Diagnosis of the site-averaged first-shell closure on sublattices with several sites per cell.

Before the branch-resolved closure, screen_structure had two defects on such sublattices.  This script shows both.
  1. The bond phases were taken relative to the basis offset tau_b' - tau_b, but coulomb_k_basis already carries the full
     bond phase exp(i k . R).  The basis offset was therefore removed twice, and bonds between different basis sites
     averaged towards zero: G1_wip against G1_true.
  2. The structure factor was averaged over the sites into one scalar c1(k) = (1 / nz) sum_j cos(k . R_j).  That is the
     in-phase vector's Rayleigh quotient of the shell matrix C(k), so <c1^2> = 1/z only when the in-phase vector is a branch
     of C(k) at almost every k.  Otherwise the branches mix and <c1^2> < 1/z, which breaks continuity at the edge of the
     linear range.  bcc in a two-site B2 cell is the case with equality (C(k) is real, so (1, 1) is always a branch); there
     only defect 1 shows.  hcp and the NiAs anions at c/a 1.6 show both.  At c/a 1.7 the NiAs anion shell has no bonds between
     different sites, so neither defect shows.
The branch values (eigenvalues of C(k), and g_m = v_m^H G v_m) satisfy <c1_m^2> = 1/z and <g_m c1_m> = G1 in every case.
The hardness a = 1.5 Ry is fixed (not derived from beta); the script is a diagnosis, not a supercell check
(see sro_basis_test.py for that).
"""
import numpy as np
from math import erfc, sqrt, pi


def coulomb_k_basis(p, tau, kv, eta=None):
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
            g_phase = np.exp(-1j * (gv.T @ d))  # [ng]
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
                    # Gamma point neutral-projected kernel: omit G=0, subtract pi/(om*et^2)
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
    return mb


def proto(name, coa=None):
    if name == 'hcp':
        ca = sqrt(8.0 / 3.0) if coa is None else coa
        avec = np.array([[1.0, -0.5, 0.0], [0.0, sqrt(3.0) / 2.0, 0.0], [0.0, 0.0, ca]])
        frac = np.array([[0.0, 1.0 / 3.0], [0.0, 2.0 / 3.0], [0.0, 0.5]])
        label = np.array([1, 1])
    elif name == 'nias':
        ca = coa
        avec = np.array([[1.0, -0.5, 0.0], [0.0, sqrt(3.0) / 2.0, 0.0], [0.0, 0.0, ca]])
        frac = np.array([[0.0, 0.0, 1.0 / 3.0, 2.0 / 3.0],
                         [0.0, 0.0, 2.0 / 3.0, 1.0 / 3.0],
                         [0.0, 0.5, 0.25, 0.75]])
        label = np.array([1, 1, 2, 2])
    elif name == 'b2_bcc':
        avec = np.eye(3)
        frac = np.array([[0.0, 0.5], [0.0, 0.5], [0.0, 0.5]])
        label = np.array([1, 1])
    return avec, avec @ frac, label


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
            r = r[r > 1e-10]
            rmin = min(rmin, r.min())
    bonds = []
    for b1 in sites:
        for b2 in sites:
            v = tv + (tau[:, b2] - tau[:, b1])[:, None]
            r = np.linalg.norm(v, axis=0)
            for j in np.flatnonzero(np.abs(r - rmin) < 1e-7 * rmin):
                bonds.append((b1, b2, v[:, j]))
    return sites, bonds, rmin


def diagnose_wip(name, coa, alat, nk, ksro=0):
    avec, tau0, label = proto(name, coa)
    p, tau = avec * alat, tau0 * alat
    b = 2 * pi * np.linalg.inv(p).T
    sites, bonds, r1_sh = sro_shell(p, tau, label, ksro)
    _, all_bonds, r1_nn = sro_shell(p, tau, label, 0)
    nb_sro = len(sites)
    z = len(bonds) // nb_sro
    f = (np.array(np.meshgrid(*[np.arange(nk)] * 3, indexing='ij')).reshape(3, -1) + 0.5) / nk
    kv = b @ f
    mb = coulomb_k_basis(p, tau, kv)
    # 1. Check WIP w1_eig bug vs true w1
    nb = tau.shape[1]
    g00_test = 0.0
    g1_wip = 0.0
    g1_true = 0.0
    a = 1.5
    for ik in range(kv.shape[1]):
        gk = a * np.linalg.inv(a * np.eye(nb) + mb[:, :, ik])
        for s in sites:
            g00_test += gk[s, s].real
        for (b1, b2, rv) in bonds:
            ph_wip = kv[:, ik] @ (rv - (tau[:, b2] - tau[:, b1]))
            ph_true = kv[:, ik] @ rv
            g1_wip += (gk[b1, b2] * np.exp(-1j * ph_wip)).real
            g1_true += (gk[b1, b2] * np.exp(-1j * ph_true)).real
    g00_test /= (nb_sro * kv.shape[1])
    g1_wip /= (len(bonds) * kv.shape[1])
    g1_true /= (len(bonds) * kv.shape[1])
    print(f"{name} (coa={coa}, ksro={ksro}): z={z}, G00={g00_test:.6f}, G1_wip={g1_wip:+.6f}, G1_true={g1_true:+.6f}")

    # 2. Compare scalar-averaged c1(k) vs branch-resolved C(k) eigenvalues
    c1_scalar = np.array([np.mean([np.cos(kv[:, ik] @ rv) for (_, _, rv) in bonds]) for ik in range(kv.shape[1])])
    c1_branch = []
    g_branch = []
    smap = {s: i for i, s in enumerate(sites)}
    for ik in range(kv.shape[1]):
        ck = np.zeros((nb_sro, nb_sro), dtype=complex)
        for (b1, b2, rv) in bonds:
            ck[smap[b1], smap[b2]] += np.exp(1j * (kv[:, ik] @ rv)) / z
        ev_c, vc = np.linalg.eigh(ck)
        gk = a * np.linalg.inv(a * np.eye(nb) + mb[:, :, ik])
        gsro = gk[np.ix_(sites, sites)]
        gb = np.real(np.diag(vc.conj().T @ gsro @ vc))
        c1_branch.extend(ev_c)
        g_branch.extend(gb)
    c1_branch = np.array(c1_branch)
    g_branch = np.array(g_branch)
    print(f"  scalar <c1^2> = {np.mean(c1_scalar**2):.6f}, 1/z = {1.0/z:.6f}, branch <c1^2> = {np.mean(c1_branch**2):.6f}")
    print(f"  scalar <g*c1> = {np.mean([np.mean(np.diag(a*np.linalg.inv(a*np.eye(nb)+mb[:,:,ik]))[sites])*c1_scalar[ik] for ik in range(kv.shape[1])]).real:+.6f}, "
          f"branch <g*c1> = {np.mean(g_branch * c1_branch):+.6f} (== G1_true {g1_true:+.6f})")
    print(f"  scalar c1 range [{c1_scalar.min():+.4f}, {c1_scalar.max():+.4f}], branch c1 range [{c1_branch.min():+.4f}, {c1_branch.max():+.4f}]")


if __name__ == '__main__':
    alat = 2.85 / 0.529177210903
    diagnose_wip('b2_bcc', None, alat, 10, ksro=0)
    diagnose_wip('hcp', None, alat, 10, ksro=0)
    diagnose_wip('nias', 1.6, alat, 10, ksro=2)
    diagnose_wip('nias', 1.7, alat, 10, ksro=2)
