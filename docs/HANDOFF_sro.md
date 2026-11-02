# HANDOFF: short-range order (new thread)

## Handover 2026-11-03 (start here)
State in one paragraph.  SRO electrostatics are solved.  The charge-excess model (`64bda88`, module `asr_screening`) gives
V = F dq with F = a [ (G00 + z G1 Ahat)^-1 - 1 ].  At the default beta = 1 it needs no fitting.  The per-atom supercells
(`trial/b45/ks_site`) confirm it: beta 1.00 to 1.02, lambda 0.20 to 0.25, hardness 1.27 (FINDINGS_LOG "B3 tested").  The ASR SRO
energy is within 5 to 12 % of those cells and the charge change within 10 to 15 %.  One error is left.  The ASR SRO hopping
(band) charge is 2.7 times too large: -0.0057 against -0.0021 e per unlike neighbour on Cu.  This makes the Zn - Cu band
separation shift 1.45 times too much (FINDINGS_LOG "Band-separation excess traced to the SRO hopping charge").

Running on the desktop (`~/asr_remote4/trial/b45/hop`, queues H1 and H2, started 2026-11-03 07:25, binary `bin/asr_64bda88`).
Every run is one iteration.  `SUMMARY.txt` logs start and end.  The charge is the output qnet in `run.log` (line `it  1 ... qnet=`).
| run | what it measures |
|---|---|
| `ks_a0_s1`, `ks_a0_s2`, `ks_a-0.15_s1`, `ks_a-0.15_s2`, `ks_a0.15_s1`, `ks_a-0.30_s1` | per-species SQS, free-atom potentials: pure hopping charge of the exact SRO structure |
| `at_L4_{rand,m015,p015,m027}` | ASR, same free-atom potentials, levels 4, 24 steps: same quantity in the ASR SRO model |
| `at_L4s36_{rand,m015}`, `at_L5_{rand,m015}` | the same at 36 steps and at levels 5 (truncation) |
| `fz_L4_{rand,m015,p015,m027}`, `fz_L5_{rand,m015}` | ASR from the converged random potential: hopping charge at the physical potential (cell value -0.0012 e at -0.148) |
How to read them.  Hopping charge = q_Cu(alpha) - q_Cu(0) at the same settings.  Per unlike neighbour divide by dn = 8 x 0.5 x
(-alpha) (0.594 at -0.148, 1.062 at -0.266).  If `at_L4` differs from the SQS by the same 2.7, the SRO model is the cause.
If `at_L5` or `at_L4s36` closes the gap, truncation is the cause.
Decision tree for the new thread:
1. Truncation is the cause: raise the defaults for SRO runs, after a cost check.
2. The SRO model is the cause: the shell atoms need their own conditional environment.  A neighbour of a Cu centre is itself a
   centre with an SRO shell.  Derive the next order of the Mookerjee-Prasad model (correlated second step) in general, prove its
   limits (alpha = 0, B2 at alpha = -1, the sum rules), then implement.  Do not fit a scaling to CuZn.
3. After that: the second alloy (fcc Fe-Ni or bcc Fe-Cr, the PRB target) for the full test; then the random-alloy ASR vs
   k-space gap (Zn d 0.057 Ry deeper in the ASR).
Analysis scripts (desktop `trial/b45/analysis`): `regress.py` (V and q fits), `qhop.py` (hopping fit), `cshift.py` (C_d shifts),
`ckappa.py`, `mom3.py` (DOS band centres), `cef.py` (model constants).  Run them from `trial/` or `trial/b45` as in FINDINGS_LOG.
Uncommitted on the laptop, on purpose: `sro_argon.md`, `sro_kimi.md`, `suggestion.md` (other models' derivations and notes).

Written 2026-10-30.  This thread works on ROADMAP phase 6 (short-range order) only.
Read first: this file, then `docs/HANDOFF.md` (working rules, desktop), `docs/STATUS.md`, `docs/ROADMAP.md` phase 6.
The full record of the SRO work so far is `docs/FINDINGS_LOG.md`, section "Short-range order: first-shell Warren-Cowley ...".
The original plan with gates B1 to B5 is `docs/archive/HANDOFF_hea_sro.md`, Part B.

## Model (implemented)
- First-shell Warren-Cowley SRO between a site and its first shell only (Mookerjee and Prasad PRB 48, 17724 (1993); Alam and
  Mookerjee JPCM 21, 195503 (2009), section 2.1).  P(c, s) = x_c y_s (1 - alpha_cs).  alpha > 0 means the pair avoids contact.
  alpha < 0 means the pair is favoured.
- Given the centre species c, the shell sites are independent with X_s^(c) = y_s (1 - alpha_cs).  Sites beyond the first shell are
  uncorrelated.  Shell sites are not correlated among themselves.  This is the model of the papers.
- Our route gets the centre projector for free.  `asr_ed` freezes the start site to the species of component c.  So the recursion of
  component c is the uncorrelated one, except that the first shell of the origin carries X^(c) instead of x.
- Component weights stay f_kappa x_c, because sum_c x_c X^(c) = x.
- Electrostatics: the charge-excess model with the conditional first shell (`asr_screening`, `64bda88`).  Without SRO it is the SIM.

## Code (committed on main 2026-10-30: `176a727` SRO, `f2f9ca7` spin-resolved DOS, `455a069` the user's docs)
| file | change |
|---|---|
| `src/asr_sro.f90` (new) | `sro_conditional`: X^(c) from x, y and alpha.  One sublattice: symmetric unlike pairs, diagonal derived.  Two sublattices: full table, both marginals checked.  X >= 0 checked. |
| `src/asr_augcls.f90` | `augcls_set_class_conc` resets one class's mt tables (node set unchanged).  `class_mt` builds the local basis on the species of non-zero weight only (a zero weight makes `augn_basis` singular). |
| `src/asr_class_scf.f90` | With SRO the operator has 2 ncls classes.  Class ncls + kappa' is the first shell of the origin on sublattice kappa' (minimum image).  `chain_worker` sets these classes to `sro_cond` of the frozen component.  Error stop for SRO between sublattices that are not first-shell neighbours. |
| `src/asr_driver.f90` | `[[sro]]` reader (`read_sro`), report block in `asr.out`.  Also the spin-resolved DOS columns of `dos.dat` (not SRO; committed separately as `f2f9ca7`). |
| `test/test_sro.f90` (new), `test/test_augcls.f90`, `test/test_run_input.f90` | unit gates below |
| `example/fe66ni14cr20_fcc_sro.toml` (new) | Alam triad 1-1-0 at levels 4, `hull` |

Input: `[[sro]] sublattices = ["A", "B"]` (optional with one sublattice), `shell = 1`, `pairs = [["Fe", "Ni", alpha], ...]`.

Tests on the merged main plus this work (laptop, 2026-10-30): `test_sro` 16/16, `test_augcls` 16/16, `test_run_input` 33/33,
`test_terminators` 11/11.  The full suite on committed main passed (41 programs, 306 checks; `docs/STATUS.md`).

Review on 2026-10-30 found no correctness bug.  One input fix went in before the commit: `read_sro` now rejects a repeated pair also
when its first alpha is 0 (`test_run_input` 34/34).  `suggestion.md` stays uncommitted.  `doc/build/` is ignored.

## Gates
| gate | status |
|---|---|
| B4.1 alpha = 0 equals the uncorrelated run (unit test; SCF gate G1: every printed number identical) | 🟢 |
| B4.2 brute force: ring of 5 with frozen centre and conditional shell (1.6e-15); fcc ternary D4h-reduced vs unreduced (8.3e-13) | 🟢 |
| B3 electrostatics with SRO: charge-excess model `64bda88` (FINDINGS_LOG "B3 re-derived", "B3 tested").  Per-atom cells confirm beta, lambda and the hardness; SRO energy within 5 to 12 %, charge change within 10 to 15 %.  The earlier options (SIM unchanged, `sro_shell`, (1 - beta) weight) are superseded. | 🟢 electrostatics; 🟠 hopping charge 2.7 times too large (test running) |
| B4.3 Tarafder et al. 2006 bcc CuZn: PRB version in `../references` (Fig. 5 SRO DOS, Fig. 6 band energy; digitised from the arXiv source, same figures) against our eta = 0.005 DOS (FINDINGS_LOG "Gate B4.3") | 🟢 B2 limit, Zn doublet, main Zn peak; 🟠 Cu band at moderate alpha; energies not comparable |
| B4.4 Alam 2009 Fig. 4 and 5 at levels 3 and 4: all claims hold with the paper's A and B swapped (Ni, Fe, Cr), all fail with its printed labels; our pair assignment checked; dE changes by at most 0.14 mRy from levels 3 to 4; our SRO effects about 10 times smaller at both levels.  Figure 2: the paper's Cr holds 1.6 e fewer electrons below E_F than ours (about 2 e charge transfer off Cr), so the gap comes from its random alloy (FINDINGS_LOG "Alam 2009 figure 2"). | 🟢 orderings and families; size gap explained (not ours) |
| B4.5 same-model references for bcc CuZn: per-species SQS (`ks/`, 6 cells) and per-atom cells (`ks_site/`, 5 cells), k-space, nk 5, gauss smear 0.005 (`tools/sqs.py --alpha ... --classes species`; site classes via `site_classes = true`).  DOS shape agrees; energies and charges agree with the per-atom cells after B3 | 🟢 done |
| SRO between two sublattices in an SCF run (🟢 2026-11-02: T1 and T2 identical to all digits, T3 converged in 35 iterations): T1 bcc one sublattice = B2 cell with equal sublattices (exact), T2 G1 on partial B2, T3 partial B2 Cu0.7Zn0.3 / Cu0.3Zn0.7 with alpha_CuZn = -0.05 (full SCF) | 🟠 queued (`trial/b45/twosub`) |

## Runs that exist (desktop `~/asr_remote4/trial/sro_bp`, 2026-10-24)
They were run from a plain copy `~/asr_sro`, retired on 2026-10-30.  Its binary and `SOURCE_MD5` are in `~/asr_remote4/trial/archive_clones/asr_sro`.  All runs used the
`bp` tail (the old default).  `SUMMARY.txt` has start, end and wall times.

Fe66Ni14Cr20 fcc NM, levels 3, 24 steps, `bp`.  Triad = (alpha_FeNi, alpha_FeCr, alpha_NiCr).  N(E_F) is the component total in
states/Ry/atom.
| run | iterations | E_F (Ry) | qnet Fe | Ni | Cr | N(E_F) Fe | Ni | Cr | E_total (Ry/site) |
|---|---|---|---|---|---|---|---|---|---|
| random (g1_nosro, 2 iterations from converged) | | -0.070292 | 0.0265 | 0.1460 | -0.1896 | 12.09 | 7.84 | 8.90 | -2522.70766986 |
| t110 | 15 | -0.071734 | 0.0160 | 0.1246 | -0.1401 | 12.38 | 9.67 | 9.59 | -2522.70724592 |
| t101 | 19 | -0.073827 | 0.0293 | 0.1117 | -0.1749 | 11.45 | 12.34 | 8.65 | -2522.70975395 |
| t011 | 16 | -0.071717 | 0.0097 | 0.1285 | -0.1219 | 12.90 | 9.48 | 9.89 | -2522.70509509 |
| t100 | 16 | -0.071886 | 0.0316 | 0.1242 | -0.1913 | 11.44 | 9.88 | 8.81 | -2522.70973062 |
| t010 | 13 | -0.070298 | 0.0113 | 0.1461 | -0.1395 | 12.94 | 7.81 | 9.62 | -2522.70511344 |
| t001 | 15 | -0.072053 | 0.0247 | 0.1284 | -0.1715 | 12.04 | 9.37 | 8.63 | -2522.70837832 |
| t111 | 18 | -0.073287 | 0.0138 | 0.1125 | -0.1244 | 12.37 | 12.30 | 9.90 | -2522.70669202 |

bcc CuZn 50:50, a = 2.85 A, levels 4, 24 steps, `bp`:
| run | alpha_CuZn | iterations | E_F (Ry) | qnet Cu | N(E_F) Cu | Zn | E_total (Ry/site) |
|---|---|---|---|---|---|---|---|
| cuzn_rand | 0 | 32 | 0.001302 | 0.1094 | 1.779 | 0.906 | -3446.83325127 |
| cuzn_m015 | -0.15 | 32 | 0.001115 | 0.1043 | 1.747 | 0.908 | -3446.83370998 |
| cuzn_m100 | -1 (first shell all unlike) | 33 | 0.002440 | 0.0790 | 1.571 | 0.920 | -3446.83479384 |

First observations (not checked further): ordering lowers the total energy (-0.46 mRy at -0.15, -1.54 mRy at -1).  The charge
transfer drops when the shell is ordered.  This energy trend comes with an electrostatic model that ignores the shell composition
(B3), so it is not yet a physical result.

## Next actions, in order
2026-11-03: B3 tested (FINDINGS_LOG "B3 tested").  The site-class cells give beta 1.00 to 1.02 and lambda 0.20 to 0.25, as the
charge-excess model predicts.  The new model (`64bda88`, beta 1) gives the SRO energy within 5 to 12 % and the charge change within
10 to 15 % of those cells.  Next: (a) find why the ASR band separation shifts 1.45 times too much; (b) the same test on a second
alloy (fcc Fe-Ni or bcc Fe-Cr); (c) the random-alloy ASR vs k-space gap (item 4 below).
2026-11-02: B3 derived and implemented (`22366b9`, FINDINGS_LOG "B3: screened shell"): V_p(shell) = -beta 2 (q_p - qbar)/R1 +
(1 - beta) 2 Delta_p/R1.  The species-class cells turned out to be a beta' = 0.5 model (finite size).  Running in `trial/b45`:
queue P (ASR beta = 0.495, `screen/`) and queue S (site-class cells, `ks_site/`).
1. Queue P vs the species-class cells: energies, Cu charge, Zn - Cu band separation at the four alpha.  Agreement tests the derivation.
2. Queue S: beta from the random site-class cell (fit V_i against q_i, as for the SIM), then the SRO energies and charges; compare
   with the ASR at that beta.  This sets the default sim_beta (now 1, not calibrated).
3. A second alloy (fcc, Fe-Ni or Fe-Cr) for the same two tests.
4. Open, independent of SRO: ASR vs k-space for the same random alloy (Zn d band 0.057 Ry deeper in the ASR, 62 mRy per atom).

## Rules that matter here
- All calculations on the desktop; unit tests on the laptop.  Code changes only on the laptop: commit, push to GitHub, pull on the desktop, then run.
- Production route `asr_ed` only.  Default tail `hull`; do not use `bp` for new runs (FINDINGS_LOG "Levels-4 SCF stall").
- A fix or a model choice must hold in general, not only for the test system (memory `general-fixes-not-patches`).
- Docs in simple sentences, no dashes.
