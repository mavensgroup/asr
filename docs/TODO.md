# TODO (ordered; updated each session)

## HANDOFF 2026-10-09 (read first; details in STATUS.md, plan in ROADMAP.md "Real-disorder roadmap")
**Decision (user, explicit and repeated):** the production method is the ENERGY-DEPENDENT ASR in the AMscf sense, implemented in AMscf's literal
form.  Not KKR-ASA, and NOT the energy-independent route (do not substitute it on equivalence grounds; ask first if a shortcut seems justified).
AMscf form (ASR/main.f90 ap1..ap13, ASR/hop_endep.f90): for central species c, H(E) = E - T(E - E_nu - h)T - T hoh T, T = Delta_c^1/2 Delta^-1/2:
linear P(E) = (E - C)/Delta, structure constants Delta_c^1/2 S Delta_c^1/2 (non-random), disorder only in the site-diagonal E-dependent term
(and in the scaled o), central site fixed to species c, one recursion per species, G_c(E) = <1|(E - H(E))^-1|1>.
Cross-checks only: EI H^(2)/H^(3) augmented-space route (E - H(E) = T (E - H^(2)) T, so identical G at infinite depth) and the exact KKR-ASA code
(`asr_pfun`, `potential_function = "exact"`, `kkr_check`).
**Workflow:** code sync via git hub `desktop:~/asr.git` (`asr_sync` / `asr_pull` in env.sh); desktop needs the SRM VPN; heavy runs on the
desktop in `~/asr_remote4/trial/*`; after module type changes `rm -rf build/ifx_*` (fpm keeps stale objects); never `pkill -f` a pattern that
matches the ssh command; one process per run directory.
**State:**
- Ordered systems validated: EI order 3 fixes CuZn E_F and FeCo Co moment (depth 64); exact KKR-ASA ED == same-model k-space KKR (D0 passed).
- Cross-check data for D2 (EI order 2, CuZn bcc, depth 4 rank 1, trial/d2): 50:50 matches Banerjee 2010 (Zn d 4 mRy, Cu d 21 mRy); the
  concentration trend is opposite to that paper, whose Fig. 3 is internally inconsistent (Cu0.75 panel holds +1.04 e below its "E_F").
  EI order-3 runs trial/d2/x{24,50,75}_o3 were running at handoff (cross-check only).
**Next (in order):**
1. [ ] Implement the AMscf energy-dependent ASR in the augmented-space operator: central site fixed to species c (non-flippable start site),
       Delta_c-scaled non-random S, site-diagonal E-dependent random term E(1 - Delta_c/Delta) + Delta_c C/Delta, scaled o term (-T hoh T);
       complex-symmetric Lanczos per contour energy (as in asr_endep) or AMscf's real seeds + spline; one recursion per species and orbital.
2. [ ] Gates: ordered limit == ordered ED/EI; disordered CuZn vs the EI cross-check at the same depth (equal at infinite depth; compare depth
       convergence and node counts); positive DOS; sum rules.
3. [ ] Production D2 on the AMscf ED route: CuZn x = 0.24 / 0.5 / 0.75; second independent concentration reference (KKR-CPA / supercells).
4. [ ] D3 electrostatics (SIM beta), D4 magnetic disorder (FeCo, FeV: Chakrabarti & Mookerjee 2005 Fe 1.12, V -0.11 muB) on the AMscf ED route.
5. [ ] SCF stall criterion (stop when dV, E_F and moments are stationary above tol).

## Now: finish M4 (own LMTO) in this order
1. [x] Spin polarisation (collinear) — DONE 2026-10-02: Fe 2.18, Ni 0.64, Cu 0.00 muB (`test_asa_mag`),
       SCF generalised to multiple classes + spin; `test_asa_ag` unchanged. Moments converge to Fe ~2.23, Ni ~0.64
       with smaller smearing (Gaussian bias ~2%).
2. [x] Multi-atom cells — DONE 2026-10-02: classes with own spheres, own Ewald Madelung (CsCl 1.762675, NaCl
       1.747565 to 4e-7), charge transfer, Anderson mixing of (smooth potential, E_nu). B2 FeCu from scratch:
       Fe 0.657/0.789/6.744, Cu 0.707/0.687/9.417 vs oracle 0.651/0.777/6.763, 0.706/0.677/9.426; E_F -0.1935 vs -0.1871.
       Still open: choosing sphere radii by a rule (neutrality / equal-volume); currently an input (oracle ratio in the test).
3. [x] (DONE) Fermi-level/BZ integration — linear tetrahedra (Bloechl weights) default, Gaussian kept. Fe moment
       2.204 (nk=12) -> 2.214 (16) -> ~2.23; Ni 0.603/0.619/0.628/0.634 at nk=12/16/20/28 -> Gaussian small-sigma 0.640.
       Original item: tetrahedron or Methfessel-Paxton (Gaussian smearing biases N(E_F)
       when E_F sits on a sharp d peak, seen in FeCu).
4. [x] (DONE: `asr example/ag_fcc.toml` runs SCF + params + DOS in 5 s; `test_lmto_input`) TOML `kind = "lmto"`: structure + species -> SCF -> parameters -> DOS (k-space and recursion).
5. [ ] Explain the 10-30 mRy absolute offsets vs Stuttgart Ag (E_nu, C, V(S)); see STATUS.

## Now (autonomous, decided 2026-10-02): M3 symmetry reduction of the augmented space
- [x] cubic point-group orbit reduction (48 ops): `asr_symmetry` + `build_augmented_operator_sym`; coefficients identical to
      the unreduced ones (3e-13 sc, 6e-13 bcc, 3e-11 fcc) with 32x / 39x / 46x fewer states (`test_augsym`).
- [ ] (running) depth 7-9 disordered terminator study (sc/bcc/fcc, off-diagonal disorder, concentration scan).

## M5 (alloys with LMTO Hamiltonian) — in progress
- [x] core operator `asr_augorb`, exact vs brute force.
- [ ] multi-shell neighbour graph + bond-block table from our screened S^alpha (`asr_strconst`), sc/bcc/fcc/B2.
- [x] pure-limit gate (`test_alloy_ordered`: coefficients == k-space to 1e-14).
- [ ] virtual-crystal / no-disorder limit (A=B parameters) same check; then a real binary alloy vs KPM/supercell truth.
- [x] allocation-free apply.
- [x] (DONE, `asr_augsym`, default in alloy_scf, TOML `kind="alloy"` driver) orbital-aware point-group reduction for `augorb` (little group of the start orbital; ~16x) so exact  [DONE: asr_augsym verified; next: use it in alloy SCF for depth>=4, TOML kind="alloy"]
      rank-limited generation reaches useful depth; drop flip-radius as an accuracy control (keep as a cheap exploratory option).
- [x] M6 core DONE for the non-magnetic, single-site level (`asr_alloy_scf`, ordered-limit gate passes). Remaining: real disorder runs, spin, SRO.
- [ ] (design record) M6 alloy SCF from recursion moments:
      * density from orthogonal-representation moments: rho_l = M0 phi^2 + 2 M1 phi phidot + M2 (phidot^2 + phi phiddot),
        M_k = int (E-E_nu)^k N_l(E) dE (phiddot already available from the variational equations);
      * species-resolved LDOS: <n_R G_RR> = transition element between u1=|R,{}> and M_R u1 = x u1 + s |R,{R}>;
        first column of (z-T)^-1 from the CF tail functions, [G]_k1 = G_11 prod_j (b_j g_{j+1});
      * energies via a complex contour with the terminated Green function; E_F from the electron count;
      * Madelung with species charges weighted by concentration (single-site level first; SRO/charge-correlation later);
      * gate: ordered limit (x=1) reproduces `asa_scf`; CuNi / CuZn / FeCr vs literature DOS.

## Next (set 2026-10-02)
- [ ] CuZn depth 4/5 (desktop): check Cu d / Zn d separation vs Tarafder 2006 (~0.37 Ry) and vs ordered B2 ASA (0.32 Ry).
- [ ] Validate the SIM (sign now screening) vs supercell / ordered-limit electrostatics; beta choice.
- [ ] Depth/terminator: ordered Ni needs 25-30 levels for its moment -> terminator quality (Luchini-Nex, multiband tails) on real alloy
      coefficients is the limit for magnetic Ni alloys; consider block/orbital-resolved (instead of per-orbital) recursion and larger depth on desktop.
- [ ] Heavy runs on the DESKTOP only (user instruction).

## Then
- [ ] s,p weights: Stuttgart moments include o,p terms - check FeCu s,p +3-5% excess with the new moment code.
- [ ] Combined correction / f downfolding (high-energy sp band; only if needed).
- [ ] M2 augmented space (binary), exact-average tests, then re-run the terminator comparison on disordered coefficients.
- [ ] Terminator workstream: read Luchini-Nex 1987 original if obtainable; multiband/gap-aware tails on real alloy coefficients.
- [ ] Symmetry reduction + performance (M3), HEA/SRO (M7), spectral/optics/exchange (M8), non-collinear/altermagnets (M9).

## Housekeeping
- [x] Explain the ~5% Fe moment excess of the alloy route: H^(2) is representation dependent (sigma 0.7 -> alpha_d 0.0168 vs standard 0.0107) and the
      kept S^alpha range rcl 1.2 is short.  sigma 0.64 + rcl 2.0 -> 2.240 vs ASA 2.226 (STATUS.md, 2026-10-06).
- [x] (2026-10-06, user approved) recursion route uses the standard TB-LMTO screening constants by default; with them the kept range rcl = 1.2 is
      converged (no cost increase).  Re-baselined on the desktop, 35/35 tests pass (STATUS.md).
- [x] B2 FeCo Co moment: third-order H (order = 3, h p h overlap term) + depth 64 -> Co 1.779 vs k-space ASA 1.774 (STATUS.md).
- [ ] B2 FeCo Fe moment 2.810 at order 3 / depth 64 vs 2.770 (k-space ASA), 2.75 (Elk); E_F 12 mRy below the ASA.
- [ ] Decide the default order (3 fixes CuZn E_F and FeCo Co; costs 2.5x per H application and ~6x augmented-space nodes at depth 3).
- [ ] Better sphere charges: second-order O^-1/2 (1 - X/2 + 3X^2/8) halves the charge loss; or Lowdin totals + sphere partition.
- [x] H^(2) sp error in B2 CuZn E_F: fixed by order = 3 (0.02 mRy at depth 16).
- [ ] Ni-type moments: depth >= 30 even with the better representation; ED route stability for magnetic metals (Fe ED SCF wanders at depth 8-16).
- [ ] Spin seed that polarises only the valence region (zeta0=0.4 over the core caused runaway).
- [ ] Desktop jobs: VPN needed (`nmcli con up id "SRM"`); results land in ~/asr_remote3/trial.
- [ ] fpm incremental builds can leave stale objects after a module signature change (seen: test_ordered_dos) -> do a clean
      rebuild (`rm -rf build/ifx_*`) when a signature changes; consider a `make check` wrapper that does it.
- [ ] array arguments aliasing reallocated host arrays (seen in the augmented-space builder): copy to locals before calling growers.
- [ ] `trial/` outputs are scratch; promote anything permanent to `example/` or `data/reference/`.
- [ ] Keep Stuttgart numbers as test oracles only (see memory: lmto-no-copying-rule).
