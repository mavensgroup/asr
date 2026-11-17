#!/usr/bin/env python3
"""
Complete numerical verification script for Q1 - Q8 of docs/HANDOFF_spectral_argon.md:
1. Exact off-site scaling (Q2): Form 2A (frozen pair) and Form 2B (single-centre dressed target)
   for both H^(2) (route edc) and full ASA (route ed3), plus first-shell SRO (Q7).
2. Off-diagonal Lanczos extraction (Q3): Polarisation (a), Krylov projection (b), and Block Lanczos (c),
   testing exact Galerkin residual orders.
3. 3D Simple Cubic 3-species alloy (Q4, Q5, Q6, Q8):
   - A(k,E) from energy-dependent ASR (levels=3, rank=1, m=30 steps, hull tail) vs a sampled average (120 configurations)
     on the same 4 x 4 x 4 periodic cell (64 sites, nearest-neighbour hops).  Levels 3 covers the whole torus, so this tests
     rank 1 and sampling, not a lattice alloy.  The BZ sum rule on the commensurate k grid is an identity of the Fourier sum.
   - Five exact limits of Q8: Ordered limit, Equal-species limit, Dilute limit, DOS & Energy sum rules,
     Fejer window & Frobenius positive projection (Q6).
   - Test 5: exact enumeration with d and o replaced by their mean values (the error of dropping the off-diagonal
     disorder).  It is not an energy-independent recursion.
"""
import numpy as np
from itertools import product

np.set_printoptions(precision=6, suppress=True)


def augn_basis(conc):
    c = np.asarray(conc, dtype=float)
    n = len(c)
    sup = np.where(c > 0.0)[0]
    ns = len(sup)
    phi_sub = np.zeros((ns, ns), dtype=float)
    c_sub = c[sup]
    phi_sub[0, :] = 1.0
    for k in range(1, ns):
        v = np.zeros(ns, dtype=float)
        v[k] = 1.0
        for j in range(k):
            proj = np.sum(c_sub * v * phi_sub[j, :])
            v = v - proj * phi_sub[j, :]
        nrm = np.sqrt(np.sum(c_sub * v * v))
        phi_sub[k, :] = v / nrm
    phi = np.zeros((n, n), dtype=float)
    for l in range(ns):
        for j in range(ns):
            phi[l, sup[j]] = phi_sub[l, j]
    return phi


def local_weight_matrix(conc, par):
    phi = augn_basis(conc)
    c = np.asarray(conc, dtype=float)
    p = np.asarray(par, dtype=float)
    return (phi * (c * p)[None, :]) @ phi.T


def local_weight_matrix_complex(conc, par):
    phi = augn_basis(conc)
    c = np.asarray(conc, dtype=float)
    p = np.asarray(par, dtype=complex)
    return (phi * (c * p)[None, :]) @ phi.T


def build_ring_s(n_sites, t1=-0.25, t2=-0.06):
    S = np.zeros((n_sites, n_sites), dtype=float)
    for i in range(n_sites):
        S[i, (i + 1) % n_sites] = t1
        S[i, (i - 1) % n_sites] = t1
        S[i, (i + 2) % n_sites] = t2
        S[i, (i - 2) % n_sites] = t2
    return S


def h2_realization(S, cfg, enu, c_par, d_par, o_par):
    e = enu[cfg]
    c = c_par[cfg]
    d = d_par[cfg]
    o = o_par[cfg]
    h = np.diag(c) + (d[:, None] * S * d[None, :])
    return np.diag(e) + h - h @ np.diag(o) @ h


def ed3_realization(S, cfg, z, enu, c_par, d_par, o_par, pgam_par):
    e = enu[cfg]
    c = c_par[cfg]
    d = d_par[cfg]
    o = o_par[cfg]
    pg = pgam_par[cfg]
    p = pg + o * o
    eps = z - e
    h = np.diag(c) + (d[:, None] * S * d[None, :])
    I = np.eye(len(cfg))
    K = (
        np.diag(eps)
        - h
        + np.diag(eps * o) @ h
        + h @ np.diag(o * eps)
        - h @ np.diag(o) @ h
        + h @ np.diag(p * eps) @ h
    )
    G = np.linalg.inv(K)
    A0 = I + np.diag(o) @ h
    A2 = np.diag(np.sqrt(pg)) @ h
    G_asa = A0 @ G @ A0.T + A2 @ G @ A2.T
    return G, G_asa, h


def build_aug_configs(n_sites, n_sp, maxrank, flippable):
    cfgs = []
    for st in product(range(n_sp), repeat=n_sites):
        ok = True
        nflip = 0
        for r in range(n_sites):
            if st[r] > 0:
                if not flippable[r]:
                    ok = False
                    break
                nflip += 1
        if ok and nflip <= maxrank:
            cfgs.append(st)
    cfg_map = {st: idx for idx, st in enumerate(cfgs)}
    return cfgs, cfg_map


def build_scaled_aug_op_edc(S, conc_site, frozen_species, z, enu, c_par, d_par, o_par, central_sp, maxrank):
    n_sites = S.shape[0]
    n_sp = len(enu)
    flippable = [r not in frozen_species for r in range(n_sites)]
    cfgs, cfg_map = build_aug_configs(n_sites, n_sp, maxrank, flippable)
    n_cfg = len(cfgs)
    dim = n_sites * n_cfg

    dc = d_par[central_sp] ** 2
    ds = d_par ** 2
    g_sp = (dc / ds) * c_par
    q_sp = (ds / dc) * o_par
    ep_sp = z - (dc / ds) * (z - enu)
    Sp = dc * S

    Wg, Wq, Wep = [], [], []
    for r in range(n_sites):
        if r in frozen_species:
            sp = frozen_species[r]
            Wg.append(np.array([[g_sp[sp]]], dtype=float))
            Wq.append(np.array([[q_sp[sp]]], dtype=float))
            Wep.append(np.array([[ep_sp[sp]]], dtype=complex))
        else:
            Wg.append(local_weight_matrix(conc_site[r], g_sp))
            Wq.append(local_weight_matrix(conc_site[r], q_sp))
            Wep.append(local_weight_matrix_complex(conc_site[r], ep_sp))

    def build_rand_diag(W_list, dtype=complex):
        M = np.zeros((dim, dim), dtype=dtype)
        for ic, st in enumerate(cfgs):
            for r in range(n_sites):
                row = r * n_cfg + ic
                if r in frozen_species:
                    M[row, row] = W_list[r][0, 0]
                else:
                    l0 = st[r]
                    for k in range(n_sp):
                        st2 = list(st)
                        st2[r] = k
                        st2_t = tuple(st2)
                        if st2_t in cfg_map:
                            ic2 = cfg_map[st2_t]
                            col = r * n_cfg + ic2
                            M[row, col] = W_list[r][l0, k]
        return M

    G_op = build_rand_diag(Wg, dtype=float)
    Q_op = build_rand_diag(Wq, dtype=float)
    E_op = build_rand_diag(Wep, dtype=complex)
    S_op = np.kron(Sp, np.eye(n_cfg))
    hp = G_op + S_op
    Hz = E_op + hp - hp @ Q_op @ hp
    return Hz, hp, cfgs, cfg_map


def build_scaled_aug_op_ed3(S, conc_site, frozen_species, z, enu, c_par, d_par, o_par, pgam_par, central_sp, maxrank):
    n_sites = S.shape[0]
    n_sp = len(enu)
    flippable = [r not in frozen_species for r in range(n_sites)]
    cfgs, cfg_map = build_aug_configs(n_sites, n_sp, maxrank, flippable)
    n_cfg = len(cfgs)
    dim = n_sites * n_cfg

    dc = d_par[central_sp] ** 2
    ds = d_par ** 2
    p_par = pgam_par + o_par ** 2
    g_sp = (dc / ds) * c_par
    eps_sp = z - enu
    et_sp = (dc / ds) * eps_sp
    uu_sp = o_par * eps_sp
    qt_sp = (ds / dc) * (o_par - p_par * eps_sp)
    Sp = dc * S

    Wg, Wet, Wuu, Wqt = [], [], [], []
    for r in range(n_sites):
        if r in frozen_species:
            sp = frozen_species[r]
            Wg.append(np.array([[g_sp[sp]]], dtype=float))
            Wet.append(np.array([[et_sp[sp]]], dtype=complex))
            Wuu.append(np.array([[uu_sp[sp]]], dtype=complex))
            Wqt.append(np.array([[qt_sp[sp]]], dtype=complex))
        else:
            Wg.append(local_weight_matrix(conc_site[r], g_sp))
            Wet.append(local_weight_matrix_complex(conc_site[r], et_sp))
            Wuu.append(local_weight_matrix_complex(conc_site[r], uu_sp))
            Wqt.append(local_weight_matrix_complex(conc_site[r], qt_sp))

    def build_rand_diag(W_list, dtype=complex):
        M = np.zeros((dim, dim), dtype=dtype)
        for ic, st in enumerate(cfgs):
            for r in range(n_sites):
                row = r * n_cfg + ic
                if r in frozen_species:
                    M[row, row] = W_list[r][0, 0]
                else:
                    l0 = st[r]
                    for k in range(n_sp):
                        st2 = list(st)
                        st2[r] = k
                        st2_t = tuple(st2)
                        if st2_t in cfg_map:
                            ic2 = cfg_map[st2_t]
                            col = r * n_cfg + ic2
                            M[row, col] = W_list[r][l0, k]
        return M

    G_op = build_rand_diag(Wg, dtype=float)
    Et_op = build_rand_diag(Wet, dtype=complex)
    Uu_op = build_rand_diag(Wuu, dtype=complex)
    Qt_op = build_rand_diag(Wqt, dtype=complex)
    S_op = np.kron(Sp, np.eye(n_cfg))
    hp = G_op + S_op
    Kp = Et_op - hp + Uu_op @ hp + hp @ Uu_op - hp @ Qt_op @ hp
    return Kp, hp, cfgs, cfg_map


def csym_lanczos(Hz, v_start, m_steps, v_targets=None):
    n = len(v_start)
    u = v_start.astype(complex) / np.sqrt(np.sum(v_start * v_start))
    up = np.zeros(n, dtype=complex)
    a = np.zeros(m_steps, dtype=complex)
    b = np.zeros(m_steps, dtype=complex)
    U = np.zeros((n, m_steps), dtype=complex)
    m_used = 0
    for k in range(m_steps):
        U[:, k] = u
        w = Hz @ u
        a[k] = np.sum(u * w)
        m_used = k + 1
        if k == m_steps - 1:
            break
        w = w - a[k] * u
        if k > 0:
            w = w - b[k - 1] * up
        bk = np.sqrt(np.sum(w * w))
        b[k] = bk
        if abs(bk) < 1e-13:
            break
        up = u
        u = w / bk
    a = a[:m_used]
    b = b[:m_used]
    U = U[:, :m_used]
    C_proj = None
    if v_targets is not None:
        C_proj = v_targets.T @ U
    return a, b, C_proj, U


def csym_block_lanczos(Hz, V_start, m_steps, z):
    n, p = V_start.shape
    Q_list = []
    U = V_start.astype(complex).copy()
    for j in range(p):
        for k in range(j):
            U[:, j] -= (U[:, k] @ U[:, j]) * U[:, k]
        U[:, j] /= np.sqrt(U[:, j] @ U[:, j])
    for step in range(m_steps):
        Q_list.append(U)
        W = Hz @ U
        for Q_prev in Q_list:
            W = W - Q_prev @ (Q_prev.T @ W)
        U_next = W.copy()
        closed = False
        for j in range(p):
            for k in range(j):
                U_next[:, j] -= (U_next[:, k] @ U_next[:, j]) * U_next[:, k]
            nrm = np.sqrt(U_next[:, j] @ U_next[:, j])
            if abs(nrm) < 1e-12:
                closed = True
                break
            U_next[:, j] /= nrm
        if closed:
            break
        U = U_next
    Q_all = np.hstack(Q_list)
    T_blk = Q_all.T @ (z * np.eye(n) - Hz) @ Q_all
    rhs = Q_all.T @ V_start
    sol = np.linalg.solve(T_blk, rhs)
    G_2x2 = V_start.T @ Q_all @ sol
    R_blk = V_start - (z * np.eye(n) - Hz) @ (Q_all @ sol)
    return G_2x2, R_blk


def cf_hull(a, b, z, C_proj=None, use_tail=True):
    m = len(a)
    tail = 0.0 + 0.0j
    if use_tail and m >= 2 and abs(b[m - 1]) > 1e-12:
        T = np.diag(a) + np.diag(b[: m - 1], 1) + np.diag(b[: m - 1], -1)
        ev = np.linalg.eigvals(T)
        ilo = np.argmin(ev.real)
        ihi = np.argmax(ev.real)
        aa = 0.5 * (ev[ilo] + ev[ihi])
        bb = 0.25 * (ev[ihi] - ev[ilo])
        if abs(bb) > 1e-14:
            d = z - aa
            s = np.sqrt(d * d - 4.0 * bb * bb)
            g1 = (d - s) / (2.0 * bb * bb)
            g2 = (d + s) / (2.0 * bb * bb)
            tail = bb * bb * (g1 if abs(bb * g1) <= abs(bb * g2) else g2)
    gt = np.zeros(m, dtype=complex)
    gt[m - 1] = 1.0 / (z - a[m - 1] - tail)
    for j in range(m - 2, -1, -1):
        gt[j] = 1.0 / (z - a[j] - (b[j] ** 2) * gt[j + 1])
    g11 = gt[0]
    if C_proj is None:
        return g11
    prod_vec = np.zeros(m, dtype=complex)
    p = gt[0]
    for k in range(m):
        prod_vec[k] = p
        if k < m - 1:
            p = p * b[k] * gt[k + 1]
    return g11, C_proj @ prod_vec


def build_sc_cluster_aug(levels, n_sp, conc, L_per=6):
    sites = []
    half = L_per // 2
    for ix in range(-half, half):
        for iy in range(-half, half):
            for iz in range(-half, half):
                sites.append((ix, iy, iz))
    site_map = {pos: i for i, pos in enumerate(sites)}
    n_sites = len(sites)

    def wrap(c):
        return ((c + half) % L_per) - half

    nn_shifts = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    nbrs = [[] for _ in range(n_sites)]
    for i, (ix, iy, iz) in enumerate(sites):
        for dx, dy, dz in nn_shifts:
            p2 = (wrap(ix + dx), wrap(iy + dy), wrap(iz + dz))
            nbrs[i].append(site_map[p2])

    origin_idx = site_map[(0, 0, 0)]
    cfgs = [(0, 0)]
    for r in range(n_sites):
        if r == origin_idx:
            continue
        for k in range(1, n_sp):
            cfgs.append((r, k))
    cfg_map = {c: i for i, c in enumerate(cfgs)}

    start_node = (origin_idx, 0)
    node_level = {start_node: 0}
    queue = [start_node]
    idx = 0
    while idx < len(queue):
        r, ic = queue[idx]
        idx += 1
        lev = node_level[(r, ic)]
        if lev >= levels:
            continue
        f_site, f_st = cfgs[ic]
        r1_list = [r] + nbrs[r]
        for r1 in r1_list:
            r2_list = [r1] + nbrs[r1]
            for r2 in r2_list:
                visited = {r, r1, r2} - {origin_idx}
                cand_cfgs = [ic]
                if f_site == 0:
                    for v in visited:
                        for k in range(1, n_sp):
                            cand_cfgs.append(cfg_map[(v, k)])
                elif f_site in visited:
                    cand_cfgs.append(0)
                    for k in range(1, n_sp):
                        cand_cfgs.append(cfg_map[(f_site, k)])
                for ic2 in cand_cfgs:
                    nd2 = (r2, ic2)
                    if nd2 not in node_level:
                        node_level[nd2] = lev + 1
                        queue.append(nd2)

    for r_idx in range(n_sites):
        if (r_idx, 0) not in node_level:
            node_level[(r_idx, 0)] = levels
            queue.append((r_idx, 0))

    nodes = queue
    node_map = {nd: i for i, nd in enumerate(nodes)}
    return sites, site_map, nbrs, cfgs, cfg_map, nodes, node_map, origin_idx


def make_sc_aug_matvec(sites, nbrs, cfgs, cfg_map, nodes, node_map, origin_idx, conc, t_nn=-0.16):
    n_nodes = len(nodes)
    n_sp = len(conc)
    phi = augn_basis(conc)

    hop_rows, hop_cols = [], []
    for i, (r, ic) in enumerate(nodes):
        for r2 in nbrs[r]:
            nd2 = (r2, ic)
            if nd2 in node_map:
                hop_rows.append(i)
                hop_cols.append(node_map[nd2])
    hop_rows = np.array(hop_rows, dtype=int)
    hop_cols = np.array(hop_cols, dtype=int)

    node_site = np.array([nd[0] for nd in nodes], dtype=int)
    node_l0 = np.zeros(n_nodes, dtype=int)
    flip_rows, flip_cols, flip_l0, flip_k = [], [], [], []
    for i, (r, ic) in enumerate(nodes):
        if r == origin_idx:
            continue
        f_site, f_st = cfgs[ic]
        l0 = f_st if f_site == r else 0
        node_l0[i] = l0
        if f_site == 0 or f_site == r:
            for k in range(n_sp):
                ic2 = 0 if k == 0 else cfg_map[(r, k)]
                nd2 = (r, ic2)
                if nd2 in node_map:
                    flip_rows.append(i)
                    flip_cols.append(node_map[nd2])
                    flip_l0.append(l0)
                    flip_k.append(k)
        else:
            flip_rows.append(i)
            flip_cols.append(i)
            flip_l0.append(0)
            flip_k.append(0)

    flip_rows = np.array(flip_rows, dtype=int)
    flip_cols = np.array(flip_cols, dtype=int)
    flip_l0 = np.array(flip_l0, dtype=int)
    flip_k = np.array(flip_k, dtype=int)
    orig_mask = (node_site == origin_idx)

    return (hop_rows, hop_cols, flip_rows, flip_cols, flip_l0, flip_k, orig_mask, n_nodes, t_nn), phi


class FastAugOp3D:
    def __init__(self, op_data, phi, conc, g_sp, q_sp, ep_sp, dc, c_central):
        hop_rows, hop_cols, flip_rows, flip_cols, flip_l0, flip_k, orig_mask, n_nodes, t_nn = op_data
        self.hop_rows = hop_rows
        self.hop_cols = hop_cols
        self.flip_rows = flip_rows
        self.flip_cols = flip_cols
        self.orig_mask = orig_mask
        self.n_nodes = n_nodes
        self.t_scaled = dc * t_nn

        Wg = (phi * (conc * g_sp)[None, :]) @ phi.T
        Wq = (phi * (conc * q_sp)[None, :]) @ phi.T
        Wep = (phi * (conc * ep_sp)[None, :]) @ phi.T

        self.wg_links = Wg[flip_l0, flip_k]
        self.wq_links = Wq[flip_l0, flip_k]
        self.wep_links = Wep[flip_l0, flip_k]
        self.g_c = g_sp[c_central]
        self.q_c = q_sp[c_central]
        self.ep_c = ep_sp[c_central]

    def apply_rand(self, w_links, val_c, x):
        y = np.zeros(self.n_nodes, dtype=complex)
        y[self.orig_mask] = val_c * x[self.orig_mask]
        np.add.at(y, self.flip_rows, w_links * x[self.flip_cols])
        return y

    def apply_hp(self, x):
        y = self.apply_rand(self.wg_links, self.g_c, x)
        np.add.at(y, self.hop_rows, self.t_scaled * x[self.hop_cols])
        return y

    def __matmul__(self, x):
        h1 = self.apply_hp(x)
        qh1 = self.apply_rand(self.wq_links, self.q_c, h1)
        h2 = self.apply_hp(qh1)
        ex = self.apply_rand(self.wep_links, self.ep_c, x)
        return ex + h1 - h2


def run_all_tests():
    print("=" * 78)
    print("TEST 1: Exact Off-Site Scaling (Q2) on 6-Site 3-Species Ring (729 configs)")
    print("=" * 78)
    n_sites = 6
    n_sp = 3
    conc = np.array([0.50, 0.30, 0.20])
    enu = np.array([-0.12, -0.25, -0.05])
    c_par = np.array([0.08, -0.14, 0.21])
    d_par = np.array([0.32, 0.25, 0.35])   # D = (0.1024, 0.0625, 0.1225), ratio 1.96
    o_par = np.array([0.45, 0.80, 0.30])
    pgam_par = np.array([1.20, 2.10, 0.85])
    S = build_ring_s(n_sites, t1=-0.28, t2=-0.07)
    z = 0.05 + 0.04j

    G_edc_enum = np.zeros((n_sites, n_sites), dtype=complex)
    G_ed3_enum = np.zeros((n_sites, n_sites), dtype=complex)
    for cfg in product(range(n_sp), repeat=n_sites):
        cfg_arr = np.array(cfg, dtype=int)
        w = np.prod(conc[cfg_arr])
        H2 = h2_realization(S, cfg_arr, enu, c_par, d_par, o_par)
        G_edc_enum += w * np.linalg.inv(z * np.eye(n_sites) - H2)
        _, G_asa, _ = ed3_realization(S, cfg_arr, z, enu, c_par, d_par, o_par, pgam_par)
        G_ed3_enum += w * G_asa

    G_edc_pair_full = np.zeros(n_sites, dtype=complex)
    G_edc_pair_invT = np.zeros(n_sites, dtype=complex)
    G_edc_pair_rk1 = np.zeros(n_sites, dtype=complex)
    G_edc_pair_rk2 = np.zeros(n_sites, dtype=complex)
    G_ed3_pair_full = np.zeros(n_sites, dtype=complex)

    conc_all = [conc for _ in range(n_sites)]
    for c in range(n_sp):
        for rk, arr in [(n_sites, G_edc_pair_full), (1, G_edc_pair_rk1), (2, G_edc_pair_rk2)]:
            Hz, _, cfgs, _ = build_scaled_aug_op_edc(
                S, conc_all, {0: c}, z, enu, c_par, d_par, o_par, c, maxrank=rk
            )
            G_scaled = np.linalg.inv(z * np.eye(Hz.shape[0]) - Hz)
            arr[0] += conc[c] * G_scaled[0, 0]
            if rk == n_sites:
                G_edc_pair_invT[0] += conc[c] * G_scaled[0, 0]

        Kp, hp, cfgs, _ = build_scaled_aug_op_ed3(
            S, conc_all, {0: c}, z, enu, c_par, d_par, o_par, pgam_par, c, maxrank=n_sites
        )
        Kinv = np.linalg.inv(Kp)
        u0 = np.zeros(Kp.shape[0], dtype=complex); u0[0] = 1.0
        v0 = hp @ u0
        a0 = u0 + o_par[c] * v0
        G_ed3_pair_full[0] += conc[c] * (a0 @ Kinv @ a0 + pgam_par[c] * (v0 @ Kinv @ v0))

    for R in range(1, n_sites):
        for c in range(n_sp):
            for s in range(n_sp):
                P_cs = conc[c] * conc[s]
                T_Rs = d_par[c] / d_par[s]
                for rk, arr in [(n_sites, G_edc_pair_full), (1, G_edc_pair_rk1), (2, G_edc_pair_rk2)]:
                    Hz, _, cfgs, _ = build_scaled_aug_op_edc(
                        S, conc_all, {0: c, R: s}, z, enu, c_par, d_par, o_par, c, maxrank=rk
                    )
                    n_cfg = len(cfgs)
                    G_scaled = np.linalg.inv(z * np.eye(Hz.shape[0]) - Hz)
                    g_0R = G_scaled[0, R * n_cfg]
                    arr[R] += P_cs * T_Rs * g_0R
                    if rk == n_sites:
                        G_edc_pair_invT[R] += P_cs * (1.0 / T_Rs) * g_0R

                Kp, hp, cfgs, _ = build_scaled_aug_op_ed3(
                    S, conc_all, {0: c, R: s}, z, enu, c_par, d_par, o_par, pgam_par, c, maxrank=n_sites
                )
                n_cfg = len(cfgs)
                Kinv = np.linalg.inv(Kp)
                u0 = np.zeros(Kp.shape[0], dtype=complex); u0[0] = 1.0
                uR = np.zeros(Kp.shape[0], dtype=complex); uR[R * n_cfg] = 1.0
                v0 = hp @ u0
                vR = hp @ uR
                q_s = (d_par[s] / d_par[c]) ** 2 * o_par[s]
                a0 = u0 + o_par[c] * v0
                aR = uR + q_s * vR
                term1 = T_Rs * (a0 @ Kinv @ aR)
                term2 = (1.0 / T_Rs) * np.sqrt(pgam_par[c] * pgam_par[s]) * (v0 @ Kinv @ vR)
                G_ed3_pair_full[R] += P_cs * (term1 + term2)

    G_edc_single_full = np.zeros(n_sites, dtype=complex)
    G_edc_single_rk1 = np.zeros(n_sites, dtype=complex)
    G_edc_single_rk2 = np.zeros(n_sites, dtype=complex)
    phi = augn_basis(conc)
    for c in range(n_sp):
        T_vec = d_par[c] / d_par
        tau_c = phi @ (conc * T_vec)
        for rk, arr in [(n_sites, G_edc_single_full), (1, G_edc_single_rk1), (2, G_edc_single_rk2)]:
            Hz, _, cfgs, cfg_map = build_scaled_aug_op_edc(
                S, conc_all, {0: c}, z, enu, c_par, d_par, o_par, c, maxrank=rk
            )
            n_cfg = len(cfgs)
            G_scaled = np.linalg.inv(z * np.eye(Hz.shape[0]) - Hz)
            arr[0] += conc[c] * G_scaled[0, 0]
            for R in range(1, n_sites):
                b_vec = np.zeros(Hz.shape[0], dtype=complex)
                for k in range(n_sp):
                    st = [0] * n_sites
                    st[R] = k
                    st_t = tuple(st)
                    if st_t in cfg_map:
                        b_vec[R * n_cfg + cfg_map[st_t]] = tau_c[k]
                arr[R] += conc[c] * (G_scaled[0, :] @ b_vec)

    print(f"  Max |G_0R(Form 2A, full rank, +1 power T_R) - G_0R(enum)| = {np.max(np.abs(G_edc_pair_full - G_edc_enum[0, :])):.3e}")
    print(f"  Max |G_0R(Form 2A, full rank, -1 power T_R) - G_0R(enum)| = {np.max(np.abs(G_edc_pair_invT - G_edc_enum[0, :])):.3e}  <-- T_R^-1 FAILS")
    print(f"  Max |G_0R(Form 2B, full rank, single-centre) - G_0R(enum)| = {np.max(np.abs(G_edc_single_full - G_edc_enum[0, :])):.3e}")
    print(f"  Max |G_0R(ed3 full ASA, full rank)           - G_0R(enum)| = {np.max(np.abs(G_ed3_pair_full - G_ed3_enum[0, :])):.3e}")
    print(f"  Rank-1 max error |G_0R - G_enum|: Form 2A (pair) = {np.max(np.abs(G_edc_pair_rk1 - G_edc_enum[0, :])):.3e}, Form 2B (single) = {np.max(np.abs(G_edc_single_rk1 - G_edc_enum[0, :])):.3e}")
    print(f"  Rank-2 max error |G_0R - G_enum|: Form 2A (pair) = {np.max(np.abs(G_edc_pair_rk2 - G_edc_enum[0, :])):.3e}, Form 2B (single) = {np.max(np.abs(G_edc_single_rk2 - G_edc_enum[0, :])):.3e}")

    alpha_ab, alpha_ac, alpha_bc = -0.20, 0.10, 0.15
    alpha_mat = np.array([
        [0.0, alpha_ab, alpha_ac],
        [alpha_ab, 0.0, alpha_bc],
        [alpha_ac, alpha_bc, 0.0],
    ])
    for i in range(n_sp):
        alpha_mat[i, i] = -np.sum(conc * alpha_mat[i, :]) / conc[i]
    cond = conc[None, :] * (1.0 - alpha_mat)

    G_sro_enum = np.zeros(n_sites, dtype=complex)
    for cfg in product(range(n_sp), repeat=n_sites):
        cfg_arr = np.array(cfg, dtype=int)
        c0 = cfg_arr[0]
        w = conc[c0] * cond[c0, cfg_arr[1]] * cond[c0, cfg_arr[n_sites - 1]]
        for r in range(2, n_sites - 1):
            w *= conc[cfg_arr[r]]
        H2 = h2_realization(S, cfg_arr, enu, c_par, d_par, o_par)
        G_sro_enum += w * np.linalg.inv(z * np.eye(n_sites) - H2)[0, :]

    G_sro_form2b = np.zeros(n_sites, dtype=complex)
    for c in range(n_sp):
        conc_sro = [conc.copy() for _ in range(n_sites)]
        conc_sro[1] = cond[c, :]
        conc_sro[n_sites - 1] = cond[c, :]
        Hz, _, cfgs, cfg_map = build_scaled_aug_op_edc(
            S, conc_sro, {0: c}, z, enu, c_par, d_par, o_par, c, maxrank=n_sites
        )
        n_cfg = len(cfgs)
        G_scaled = np.linalg.inv(z * np.eye(Hz.shape[0]) - Hz)
        G_sro_form2b[0] += conc[c] * G_scaled[0, 0]
        T_vec = d_par[c] / d_par
        for R in range(1, n_sites):
            phi_R = augn_basis(conc_sro[R])
            tau_R = phi_R @ (conc_sro[R] * T_vec)
            b_vec = np.zeros(Hz.shape[0], dtype=complex)
            for k in range(n_sp):
                st = [0] * n_sites
                st[R] = k
                st_t = tuple(st)
                if st_t in cfg_map:
                    b_vec[R * n_cfg + cfg_map[st_t]] = tau_R[k]
            G_sro_form2b[R] += conc[c] * (G_scaled[0, :] @ b_vec)

    print(f"  Max |G_0R(SRO Form 2B, full rank)           - G_0R(MP enum)| = {np.max(np.abs(G_sro_form2b - G_sro_enum)):.3e}")

    print("\n" + "=" * 78)
    print("TEST 2: Off-Diagonal Extraction Accuracy Orders (Q3: Polarisation vs Projection)")
    print("=" * 78)
    Hz, _, cfgs, cfg_map = build_scaled_aug_op_edc(
        S, conc_all, {0: 0}, z, enu, c_par, d_par, o_par, 0, maxrank=1
    )
    n_cfg = len(cfgs)
    dim = Hz.shape[0]
    K_mat = z * np.eye(dim) - Hz
    K_inv = np.linalg.inv(K_mat)
    va = np.zeros(dim, dtype=complex); va[0] = 1.0
    vb = np.zeros(dim, dtype=complex); vb[1 * n_cfg] = 1.0
    g_ab_exact = vb @ K_inv @ va

    print("  steps m |  ||r_m||_2   | err Projection (Q3b) | err Polarisation (Q3a) | err Petrov-Galerkin")
    for m in [4, 6, 8, 10, 12, 16]:
        a_cf, b_cf, C_proj, U_a = csym_lanczos(Hz, va, m, v_targets=vb[:, None])
        _, f_proj = cf_hull(a_cf, b_cf, z, C_proj=C_proj, use_tail=False)
        g_proj = f_proj[0]
        T_a = z * np.eye(len(a_cf)) - (np.diag(a_cf) + np.diag(b_cf[:-1], 1) + np.diag(b_cf[:-1], -1))
        e1 = np.zeros(len(a_cf), dtype=complex); e1[0] = 1.0
        xm_a = U_a @ np.linalg.solve(T_a, e1)
        rm_a = va - K_mat @ xm_a
        rnorm = np.linalg.norm(rm_a)

        a_b, b_b, _, U_b = csym_lanczos(Hz, vb, m)
        T_b = z * np.eye(len(a_b)) - (np.diag(a_b) + np.diag(b_b[:-1], 1) + np.diag(b_b[:-1], -1))
        e1_b = np.zeros(len(a_b), dtype=complex); e1_b[0] = 1.0
        xm_b = U_b @ np.linalg.solve(T_b, e1_b)
        g_pg = g_proj + (xm_b @ rm_a)

        vp = va + vb
        vm = va - vb
        ap, bp, _, _ = csym_lanczos(Hz, vp, m)
        am, bm, _, _ = csym_lanczos(Hz, vm, m)
        gp = np.sum(vp * vp) * cf_hull(ap, bp, z, use_tail=False)
        gm = np.sum(vm * vm) * cf_hull(am, bm, z, use_tail=False)
        g_pol = 0.25 * (gp - gm)

        print(f"    {m:3d}   |  {rnorm:.3e}  |      {abs(g_proj - g_ab_exact):.3e}       |       {abs(g_pol - g_ab_exact):.3e}        |      {abs(g_pg - g_ab_exact):.3e}")

    print("\n" + "=" * 78)
    print("TEST 3: Block Recursion (Q3c) Second-Order Galerkin Error Identity")
    print("=" * 78)
    V_start = np.column_stack([va, vb])
    G_2x2_exact = V_start.T @ K_inv @ V_start
    for m_blk in [2, 3, 4, 6]:
        G_2x2, R_blk = csym_block_lanczos(Hz, V_start, m_blk, z)
        err_true = abs(G_2x2[0, 1] - G_2x2_exact[0, 1])
        err_id = abs((R_blk.T @ K_inv @ R_blk)[0, 1] - (G_2x2_exact[0, 1] - G_2x2[0, 1]))
        rnorm_blk = np.linalg.norm(R_blk)
        print(f"  block steps m={m_blk:2d} (2m={2*m_blk:2d} matvecs): ||R_m||={rnorm_blk:.3e}, |G_0R - G_blk|={err_true:.3e}, identity defect={err_id:.3e}")

    print("\n" + "=" * 78)
    print("TEST 4: 3D Simple Cubic 3-Species Alloy (Fe0.66 Ni0.14 Cr0.20 model, Q4-Q8)")
    print("=" * 78)
    conc_3d = np.array([0.66, 0.14, 0.20])
    enu_3d = np.array([-0.10, -0.22, -0.02])
    c_3d = np.array([0.04, -0.12, 0.16])
    d_3d = np.array([0.38, 0.30, 0.42])   # D = (0.1444, 0.0900, 0.1764), ratio D_Cr/D_Ni = 1.96
    o_3d = np.array([0.35, 0.60, 0.25])
    t_nn = -0.22

    levels = 3
    L_sc = 4
    sites, site_map, nbrs, cfgs_3d, cfg_map_3d, nodes, node_map, origin_idx = build_sc_cluster_aug(levels, 3, conc_3d, L_per=L_sc)
    n_nodes = len(nodes)
    print(f"  3D SC augmented space (L={L_sc}, levels={levels}, rank=1): {len(sites)} cluster sites, {n_nodes} reachable aug nodes")

    op_data, _ = make_sc_aug_matvec(sites, nbrs, cfgs_3d, cfg_map_3d, nodes, node_map, origin_idx, conc_3d, t_nn=t_nn)

    shell_sites = {}
    for idx_s, pos in enumerate(sites):
        r2 = pos[0]**2 + pos[1]**2 + pos[2]**2
        shell_sites.setdefault(r2, []).append((idx_s, pos))

    def compute_g0R_asr(z_val, enu_v, c_v, d_v, o_v, conc_v, m_steps=24, use_tail=True):
        phi_loc = augn_basis(conc_v)

        target_r_list = []
        for r2 in sorted(shell_sites.keys()):
            for idx_s, pos in shell_sites[r2]:
                target_r_list.append((idx_s, pos, r2))

        g_0R = {pos: 0.0 + 0.0j for _, pos, _ in target_r_list}
        v_start = np.zeros(n_nodes, dtype=complex)
        v_start[node_map[(origin_idx, 0)]] = 1.0

        for c in range(3):
            if conc_v[c] == 0.0:
                continue
            dc = d_v[c] ** 2
            ds = d_v ** 2
            g_sp = (dc / ds) * c_v
            q_sp = (ds / dc) * o_v
            ep_sp = z_val - (dc / ds) * (z_val - enu_v)
            Hz = FastAugOp3D(op_data, phi_loc, conc_v, g_sp, q_sp, ep_sp, dc, c)

            T_vec = d_v[c] / d_v
            tau_c = phi_loc @ (conc_v * T_vec)

            B_targets = np.zeros((n_nodes, len(target_r_list)), dtype=complex)
            for it, (idx_s, pos, r2) in enumerate(target_r_list):
                if idx_s == origin_idx:
                    B_targets[node_map[(origin_idx, 0)], it] = 1.0
                else:
                    for k in range(3):
                        ic = 0 if k == 0 else cfg_map_3d[(idx_s, k)]
                        if (idx_s, ic) in node_map:
                            B_targets[node_map[(idx_s, ic)], it] = tau_c[k]

            a_cf, b_cf, C_proj, _ = csym_lanczos(Hz, v_start, m_steps, v_targets=B_targets)
            _, f_proj = cf_hull(a_cf, b_cf, z_val, C_proj=C_proj, use_tail=use_tail)
            for it, (_, pos, _) in enumerate(target_r_list):
                g_0R[pos] += conc_v[c] * f_proj[it]
        return g_0R

    z_test = -0.05 + 0.025j
    g0R_ord = compute_g0R_asr(z_test, enu_3d, c_3d, d_3d, o_3d, np.array([1.0, 0.0, 0.0]), m_steps=30, use_tail=True)
    nk = L_sc  # exact k-grid of the L_sc^3 periodic cell
    kgrid = np.linspace(0, 2 * np.pi, nk, endpoint=False)
    kx, ky, kz = np.meshgrid(kgrid, kgrid, kgrid, indexing="ij")
    Sk = 2.0 * t_nn * (np.cos(kx) + np.cos(ky) + np.cos(kz))
    hk = c_3d[0] + (d_3d[0] ** 2) * Sk
    H2k_ord = enu_3d[0] + hk - o_3d[0] * (hk ** 2)
    Gk_ord = 1.0 / (z_test - H2k_ord)
    g00_kspace = np.mean(Gk_ord)
    g100_kspace = np.mean(Gk_ord * np.exp(-1j * kx))

    enu_eq = np.full(3, enu_3d[0])
    c_eq = np.full(3, c_3d[0])
    d_eq = np.full(3, d_3d[0])
    o_eq = np.full(3, o_3d[0])
    g0R_eq = compute_g0R_asr(z_test, enu_eq, c_eq, d_eq, o_eq, conc_3d, m_steps=30, use_tail=True)
    err_eq_ord = max(abs(g0R_eq[pos] - g0R_ord[pos]) for pos in g0R_ord)

    print(f"  Limit 1 (Ordered vs k-space): |G_00(ASR) - G_00(k)| = {abs(g0R_ord[(0,0,0)] - g00_kspace):.3e}, |G_100(ASR) - G_100(k)| = {abs(g0R_ord[(1,0,0)] - g100_kspace):.3e}")
    print(f"  Limit 2 (Equal species at x=(0.66,0.14,0.20) vs Ordered): max_R |G_0R(eq) - G_0R(ord)| = {err_eq_ord:.3e}")

    eps_dil = 1e-4
    g0R_dil = compute_g0R_asr(z_test, enu_3d, c_3d, d_3d, o_3d, np.array([1.0 - 2 * eps_dil, eps_dil, eps_dil]), m_steps=30, use_tail=True)
    diff_dil = abs(g0R_dil[(0, 0, 0)] - g0R_ord[(0, 0, 0)])
    print(f"  Limit 3 (Dilute eps={eps_dil:.1e}): |G_00(eps) - G_00(0)| = {diff_dil:.3e} = {diff_dil/eps_dil:.3f} * eps (linear O(eps) confirmed)")

    n_sc = L_sc ** 3
    sc_coords = np.array(sites, dtype=float)
    S_sc = np.zeros((n_sc, n_sc), dtype=float)
    for idx_a in range(n_sc):
        for idx_b in nbrs[idx_a]:
            S_sc[idx_a, idx_b] = t_nn

    k_pts = [
        ("Gamma", np.array([0.0, 0.0, 0.0])),
        ("Delta", np.array([np.pi / 2.0, 0.0, 0.0])),
        ("X",     np.array([np.pi, 0.0, 0.0])),
        ("M",     np.array([np.pi, np.pi, 0.0])),
        ("R",     np.array([np.pi, np.pi, np.pi])),
    ]
    e_grid = np.linspace(-0.45, 0.35, 9)
    eta = 0.03

    rng = np.random.default_rng(42)
    n_samples = 120
    sc_configs = [rng.choice(3, size=n_sc, p=conc_3d) for _ in range(n_samples)]
    H2_sc_list = [h2_realization(S_sc, cfg, enu_3d, c_3d, d_3d, o_3d) for cfg in sc_configs]

    bloch_vecs = {}
    for name, kv in k_pts:
        bloch_vecs[name] = np.exp(1j * (sc_coords @ kv)) / np.sqrt(n_sc)

    A_asr_sharp = np.zeros((len(k_pts), len(e_grid)))
    A_asr_fejer = np.zeros((len(k_pts), len(e_grid)))
    A_sc_exact = np.zeros((len(k_pts), len(e_grid)))
    dos_asr = np.zeros(len(e_grid))
    dos_bz_sum = np.zeros(len(e_grid))

    nk_bz = L_sc
    kbz = np.linspace(0, 2 * np.pi, nk_bz, endpoint=False)
    kx_b, ky_b, kz_b = np.meshgrid(kbz, kbz, kbz, indexing="ij")
    k_bz_vecs = np.column_stack([kx_b.ravel(), ky_b.ravel(), kz_b.ravel()])

    R_max_cut = np.sqrt(12.0)
    for ie, E_val in enumerate(e_grid):
        z_e = E_val + 1j * eta
        g0R = compute_g0R_asr(z_e, enu_3d, c_3d, d_3d, o_3d, conc_3d, m_steps=30, use_tail=True)
        dos_asr[ie] = -g0R[(0, 0, 0)].imag / np.pi

        acc_bz = 0.0
        for pos, gval in g0R.items():
            phase_bz = np.mean(np.exp(1j * (k_bz_vecs @ np.array(pos))))
            acc_bz += (-gval.imag / np.pi) * phase_bz.real
        dos_bz_sum[ie] = acc_bz

        for ik, (name, kv) in enumerate(k_pts):
            s_sharp = 0.0
            s_fejer = 0.0
            for pos, gval in g0R.items():
                r_len = np.linalg.norm(pos)
                ph = np.exp(1j * np.dot(kv, pos)).real
                im_g = -gval.imag / np.pi
                s_sharp += ph * im_g
                w_f = max(0.0, 1.0 - r_len / (R_max_cut + 0.5)) ** 2 * (1.0 + r_len / (2.0 * (R_max_cut + 0.5)))
                s_fejer += w_f * ph * im_g
            A_asr_sharp[ik, ie] = s_sharp
            A_asr_fejer[ik, ie] = s_fejer

        for H2_sc in H2_sc_list:
            G_sc = np.linalg.inv(z_e * np.eye(n_sc) - H2_sc)
            for ik, (name, _) in enumerate(k_pts):
                vk = bloch_vecs[name]
                gk_val = np.vdot(vk, G_sc @ vk)
                A_sc_exact[ik, ie] += (-gk_val.imag / np.pi) / n_samples

    print(f"  Limit 4 (BZ k-integral sum rule vs on-site DOS, an identity on the commensurate grid): max_E |<A(k,E)>_BZ - DOS(E)| = {np.max(np.abs(dos_bz_sum - dos_asr)):.3e}")
    print(f"  Limit 5 (Non-negativity on 5 k x 9 E, eta = {eta} Ry): min A_sharp = {np.min(A_asr_sharp):.4f}, min A_fejer = {np.min(A_asr_fejer):.4f}")
    rel_diff_sharp = np.linalg.norm(A_asr_sharp - A_sc_exact) / np.linalg.norm(A_sc_exact)
    rel_diff_fejer = np.linalg.norm(A_asr_fejer - A_sc_exact) / np.linalg.norm(A_sc_exact)
    print(f"  3D Alloy A(k,E) ASR vs {n_sc}-site cell, {n_samples}-sample Bloch average: L2 rel diff (sharp) = {rel_diff_sharp:.3f}, (Fejer, cut at {R_max_cut + 0.5:.2f} > largest distance {np.sqrt(12.0):.2f}) = {rel_diff_fejer:.3f}")
    print("\n  Sample A(k, E) values at E = -0.15 Ry (inside d-band) and E = +0.05 Ry:")
    for ik, (name, _) in enumerate(k_pts):
        print(f"    k = {name:5s}: E=-0.15 Ry -> ASR={A_asr_sharp[ik, 3]:.3f} (Fejer={A_asr_fejer[ik, 3]:.3f}, SC={A_sc_exact[ik, 3]:.3f});  E=+0.05 Ry -> ASR={A_asr_sharp[ik, 5]:.3f} (Fejer={A_asr_fejer[ik, 5]:.3f}, SC={A_sc_exact[ik, 5]:.3f})")

    print("\n" + "=" * 78)
    print("TEST 5: Mean off-diagonal parameters (exact enumeration) vs energy-dependent ASR (Q5); not an EI recursion")
    print("=" * 78)
    d_avg = np.full(3, np.sum(conc * d_par))
    o_avg = np.full(3, np.sum(conc * o_par))
    G_ei_diag = np.zeros((n_sites, n_sites), dtype=complex)
    for cfg in product(range(3), repeat=n_sites):
        cfg_arr = np.array(cfg, dtype=int)
        w = np.prod(conc[cfg_arr])
        H2_ei = h2_realization(S, cfg_arr, enu, c_par, d_avg, o_avg)
        G_ei_diag += w * np.linalg.inv(z * np.eye(n_sites) - H2_ei)

    for k_val, k_name in [(0.0, "k=0 (Gamma)"), (np.pi, "k=pi (X)")]:
        phases = np.exp(1j * k_val * np.arange(n_sites))
        A_exact_k = -np.sum(phases * G_edc_enum[0, :]).imag / np.pi
        A_ed_rk1_k = -np.sum(phases * G_edc_pair_rk1).imag / np.pi
        A_ed_s_rk1_k = -np.sum(phases * G_edc_single_rk1).imag / np.pi
        A_ei_diag_k = -np.sum(phases * G_ei_diag[0, :]).imag / np.pi
        print(f"  {k_name:12s}: A_exact={A_exact_k:.4f} | ED rank-1 (pair)={A_ed_rk1_k:.4f} (err {abs(A_ed_rk1_k-A_exact_k)/A_exact_k*100:.1f}%)"
              f" | ED rank-1 (single)={A_ed_s_rk1_k:.4f} (err {abs(A_ed_s_rk1_k-A_exact_k)/A_exact_k*100:.1f}%)"
              f" | EI (<d>,<o>)={A_ei_diag_k:.4f} (err {abs(A_ei_diag_k-A_exact_k)/A_exact_k*100:.1f}%)")


if __name__ == "__main__":
    run_all_tests()
