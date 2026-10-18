# Disordered magnetism: what is left (updated 2026-10-18)

Production route: asr_ed (AMscf form), levels 4, rank 1, recursion_steps 24, bp tail, SIM beta 1, VBH or PBE, fermi_kt 0.001 Ry.
Results: docs/STATUS.md, sections "Clean reruns of all disordered binaries" (q24) and "Fermi smearing, depth 48, SQS references, first fcc
alloys" (q25-q28).  Desktop: ~/asr_remote4/trial/q2x/SUMMARY.txt.

## CLOSED (2026-10-18): disordered binaries meet the 3-4% bar vs SQS
ASR (energy-dependent) vs our SQS with one potential per species (q32, STATUS "One potential per species SQS"): FeCo, FeNi (36 steps), NiCu
(36 steps), FeV Fe agree within ~0.2-5%; FeV V 0.07 muB off (small induced moment); FeRh Fe -4.8% (flat in steps) is the only open item.
Use 24 steps for strong moments, 36 for Ni-type.  Next project: partial (sublattice) disorder -> multicomponent -> SRO (HANDOFF below).

## HANDOFF (next thread): partial disorder, then multicomponent, then SRO (toward HEA)
Order (each with its paper gate, energy-dependent AMscf operator throughout, ed-is-target-route):
1. **Sublattice (partial) disorder**: concentrations per sublattice, some sublattices fully ordered.  Gate: Chakraborty, Chakrabarti, Mookerjee,
   cond-mat/0606737 (references/0606737v2.pdf): L2_1 Co2Fe0.4Cr0.6Al, a = 5.727 A, Co 8c, Al 4a, Fe/Cr random on 4b; Table 3 moments Co 1.075,
   Cr 0.972, Fe 2.617, Al -0.117 (their H is energy-independent second order + BP, 11-12 steps: expect ~3-4% agreement, not exact).  Also
   Co2CrAl (Co 0.760, Cr 1.548, Al -0.071) and Co2FeAl (Co 1.119, Fe 2.709, Al -0.101) as ordered end points; our own k-space LMTO for the
   ordered compounds and an SQS for the 4b disorder as the same-model references.  Existing pieces: ordered B2 mode (two sublattices,
   moments_ord), lattice = "cell" in the LMTO driver.
   PROGRESS 2026-10-19 (branch pd): general machinery done and gated (asr_augcls / asr_edc / asr_class_scf, `kind = "sublattice"`, STATUS
   "Sublattices"); k-space Co2FeAl / Co2CrAl done (Slater-Pauling 5 / 3); ASR runs + 80-atom SQS running in desktop ~/asr_pd/trial/q33.
2. **Multicomponent (N >= 3)**: occupation operator per site with N states (N-1 flips).  Gate: Alam & Mookerjee JPCM 21, 195503 (2009),
   Fe-Ni-Cr.  Node count grows fast: orbital-aware symmetry reduction + rank 1 are essential.
3. **SRO (Warren-Cowley alpha)**: gates Tarafder 2006 CuZn (alpha != 0; note the Banerjee 2010 disagreement in STATUS), then Alam 2009 ternary.
Earlier handoff (weak magnets + tail scan) is finished; its notes follow for the record.

## (finished) HANDOFF: weak magnets + tail scan
State: ASR agrees with our own SQS supercells for strong moments (FeCo; Fe in FeNi) but is wrong for weak / environment-sensitive ones:
| alloy | ASR r24 | SQS | error |
|---|---|---|---|
| fcc Ni0.75Cu0.25 | Ni 0.319 | 0.51 | -38% |
| fcc Fe0.5Ni0.5 | Ni 0.593 | 0.655 | -9% |
| bcc FeRh | Fe 2.866 (r48 2.877) | 2.99 | -0.11, not a depth effect |
| bcc FeV | Fe 1.848 (r48 1.890) / V -0.469 | 1.74 / -0.36 | ~0.1 too magnetic |
| pure fcc Ni (A = B) | 0.554 (r24), 0.609 (r36) | k-space 0.632 | chain resolution |
Three error sources (STATUS): (a) augmented-space truncation (levels, rank 1); (b) energy resolution of the chain (asymptotic width = whole spd
band ~2 Ry -> 24 steps resolve ~50-100 mRy; estimate, verify from the chain); (c) termination artefacts past the exactly generated region (hard
wall -> finite-cluster oscillation of a_n, b_n; the sqrt tail sits wherever step m lands).  A tail can only fix (c).

**Physics learned (2026-10-18, q29 / q30, STATUS "Tail scan" and "Contour path"):** the weak-moment error is the chain's energy resolution.
A d start orbital hybridises with s, p, so the chain's asymptotic width is the whole spd band (Ni: b_1 = 0.22, b_n ~ 0.62, edges -0.75 .. 1.77 Ry);
m steps resolve ~2.5 Ry / m.  A weak moment sits on a d-band edge at E_F (Ni minority d top); a coarse chain smears it across E_F -> moment
systematically LOW; strong moments (both d edges away from E_F) are immune.  The tail type and the contour path are irrelevant; more steps
help only ~1/m (A = B Ni at fixed potential: 12 steps 57%, 24 steps 74%, 64 steps 95% of the exact moment).  AMscf (12 steps, real-axis seeds,
BP) has no remedy and the published ASR moments (Ganguly 2011: Mn 0.6 muB below CPA) carry this error.  => the far environment must supply
the energy resolution: average-medium correction (item 7).  FeRh Fe (flat in m) and FeV (too magnetic) are NOT this effect.

7. [x] FAILED **Average-medium correction** (75e5732, `[asr] medium = "vca"`): G_c += G_med - CF_m(medium chain), medium = VCA of the scaled
       asr_ed operator, G_med from k space (exact).  Gate passed: A = B Ni identical at every m = 8..32, 0.5126 (kT 1 mRy, nk 60) vs tetrahedra
       0.5106.  q31 alloys (SCF from r24, kT 1 mRy): NiCu Ni unchanged at fixed potential (0.3185 vs 0.319), Cu a spurious -0.05 -> -0.08
       (killed at it 6, Ni 0.311); FeNi 2.531 / 0.616 (r36 plain: 2.555 / 0.641, SQS 2.535 / 0.662); FeRh 2.879 / 0.737 (unchanged).  The VCA
       medium has its d edge / spin asymmetry in the wrong place for an alloy site, so its truncation error is not the alloy's.  Not pursued
       further (a CPA medium would be the next step, but it is a lot of machinery for what 36 plain steps already do: NiCu -4%, FeNi -3%).
       Old codes checked (ASR_CORRECT / AMscf main, cgreen, pardos; thesis; Chakrabarti-Mookerjee EPJB): plain recursion, 12 steps, NN map,
       real-axis seeds + fit, chain-edge sqrt tail; no resolution remedy.  Their FeV 50-50: Fe 1.12 / V -0.11 (SQS 1.81 / -0.38).

Plan, in order (probe first, no SCF changes until the probe says so; production stays asr_ed):
1. [x] **Tail scan at fixed potential** (2026-10-18, q29, STATUS "Tail scan at fixed potential"): tail irrelevant (c negligible); Ni limited by
       energy resolution (b: ~1/m, exact needs m >= 64 and levels >= 5); FeRh Fe flat in m (not b, c) -> rank / levels.  Original spec:  Probe program: read a saved potential (pot.dat), build the asr_ed chain once to m_max = 48 per contour
       point, and evaluate the moments / occupations for every termination depth m = 8..48 and each tail (hull, bp, avg2) from the same chain
       (one chain contains every shorter truncation; cost = one chain per energy).  Also print b_n(E) at a few real energies near E_F and the tail
       edges per m.  Cases: pure fcc Ni (A = B; exact = k-space with the same parameters), NiCu r24 (q27/nicu_r24/pot.dat), FeRh r24 kT3
       (q25/ferh_vbh_r24_kt3/pot.dat), MnCr r24 (q23).  Output: moment vs m per tail -> how much is (b) vs (c).
2. [x] (q30: path independent to 1e-3 for Gamma >= 0.09 Ry; edges not a cause) **Path independence of the contour integral** at fixed potential: contour_points 32 vs 64, closed contour vs Fermi line (kT -> 0).  The
       tail edges are chosen per z by maxloc(real(eigenvalue)) of the complex chain (asr_endep cf_species) -> not guaranteed holomorphic in z.
       If the moments differ by > 1e-3: freeze the tail edges once per SCF iteration per species / orbital / spin.
3. [-] (not needed: q29 shows no oscillation about a stable mean) **Termination averaging** (if 1 shows oscillation about a stable mean): G(z) averaged over m in a window (e.g. 24..40), KPM-kernel analogue;
       analytic and sum-rule preserving, no extra Lanczos cost.  Gate: pure fcc Ni within 2% of k-space; NiCu / FeNi Ni towards SQS.
4. [-] (not needed: q29) **Two-band (period-2) terminator** only if 1 shows persistent period-2 oscillation of b_n (bcc d pseudogap).
5. [ ] **Rank 2** (FeRh, the cleanest case: Fe does not move with steps): first a node-count / memory probe on the desktop, then FeRh r24 rank 2
       vs SQS Fe 2.99.  Tests source (a).
6. [ ] **(b) dominates for Ni (q29)**: cost / feasibility of levels 5 + 64 steps for NiCu / FeNi (fcc alloy node count at levels 5);
       also a cheaper route to resolution (e.g. KPM-type moments or a narrower spectral window).  Original: more exact steps (levels 5 with the symmetry reduction, or a soft wall embedding the boundary in the average medium).

Pending on the desktop when this was written (2026-10-17 23:30): q27 nicu_r36, feni_r36 (restarts from r24), feni_sqs32_1; q28 (after q27)
FeCo / FeV 32-atom SQS reruns with mixing 0.1.

## Done (do not redo)
- [x] SIM sign and magnitude (326759b): V = -beta 2 qnet / R1.
- [x] SCF runaway: Anderson metric weighted by (r/(1+r))^2 + restart safeguard (ca48346); E_nu kept inside its band window (7bc52d2 ff.).
- [x] E_F search: bracketed false position (17416b2); species charges balance.
- [x] Antiparallel starts: [[species]] spin_sign (da1a24b); restarts: save_potential / start_potential / flip_spin (6692e95).
- [x] depth_tol (adaptive chain length) shown harmless at 12 / 24 / 48 steps (q16).
- [x] Our LMTO-ASA vs Questaal lm: band structure and charge transfer agree; Fe moment 1% low.
- [x] E_F noise floor: Fermi-smeared contour (cc2d513, b2367c2); kT 1 mRy is within ~0.003 muB of sharp, kT 3 mRy too coarse (q25, q27).
- [x] Independent disordered reference: SQS + lattice = "cell" / site_classes (faf2c78), sqs.py fcc triplet fix (cbd4250); results in STATUS.
- [x] External audit suggestion.md checked (STATUS "External audit"): no ASR bug in it; only gfac is real (tiny).

## Open, other
- [ ] Make fermi_kt = 0.001 the default for spin-polarised alloy runs (code default is still 0).
- [ ] k-space SCF with many classes: 3 of 4 32-atom SQS cells burst after reaching dV 1e-2..1e-3 (mixing 0.2).  If q28 (mixing 0.1) does not fix it,
      look at the Anderson history / restart with 32 classes and the E_F search in large cells.
- [ ] Generated depth (levels 5) at fixed steps for FeV and MnCr (covered by plan item 6).
- [ ] The 1% weaker magnetism vs Questaal (bcc Fe 2.252 vs 2.275, VWN): E_nu choice (~3 mRy lower in Questaal), core polarisation, Bloechl
      correction, scalar-relativistic density factor gfac (Questaal rhocor.f:298 / atomsr.f:2012; expected ~1e-4 muB).
- [ ] SCF start: FeCo overshoots E_F to +0.31 Ry and NiCu to +0.45 Ry (Cu moment 2.5 in it 1) from the atomic start; trace V_in / V_out in the
      alloy SCF, or start from pure-element potentials (start_potential).
- [ ] Magnetic state selection: parallel and antiparallel starts for Mn, Cr, V alloys; needs the alloy total energy (not implemented).
- [ ] Model limits (document, do not fix here): one collinear moment per species cannot describe sublattice AF, local-environment moment flips,
      or non-collinear order (Kaphle 2012, Ganguly 2011).
- [ ] References to re-check: FeV vs Chakrabarti-Mookerjee 2005 (theirs 1.12 / -0.11; our own SQS gives 1.74 / -0.36 in the same model, so their
      number does not fit this model); MnCr vs Banerjee-Mookerjee 2009 Fig. 1.

## After this (separate projects)
- Partial (sublattice) disorder: general lattice + basis graph, per-site configuration spaces (ordered sites never enter the cardinality
  sequence), automatic site symmetry, multi-class SCF with Ewald + SIM.  Target: Co2Fe0.4Cr0.6Al (Chakraborty 2006; experiment).
- HEA: n-component augmented space exists (asr_augn / asr_edn); needs the multi-class SCF and symmetry reduction from the item above.
- Kubo-Greenwood optical conductivity (Sadhukhan 2017): LMTO current operator + two-Green-function average (KPM route).
