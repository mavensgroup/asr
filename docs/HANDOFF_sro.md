# HANDOFF: short-range order (new thread)

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
- Electrostatics are unchanged: SIM with a random first shell (gate B3 is open).

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
| B3 electrostatics with SRO.  Decided with the user 2026-10-31: keep SIM unchanged (option a) as the default.  The explicit shell charge sum_s X_s^(p) q_s - qbar (option b) becomes a sensitivity test once its double counting of the SIM screening charge is worked out.  Report observables that change little between the two.  SRO energies stay unphysical until then. | 🟢 decided (a); 🔴 test (b) |
| B4.3 Tarafder et al. 2006 bcc CuZn, alpha about -0.15: compare the SRO-induced change of the DOS, not absolute peaks (their bands are shifted by about 0.19 Ry) | 🟠 runs exist, not analysed |
| B4.4 Alam 2009 Fig. 4 and 5 at levels 3 and 4: all claims hold with the paper's A and B swapped (Ni, Fe, Cr), all fail with its printed labels; our pair assignment checked; dE changes by at most 0.14 mRy from levels 3 to 4; our SRO effects about 10 times smaller at both levels (FINDINGS_LOG "Gate B4.4"). | 🟢 orderings and families; 🟠 size gap open |
| B4.5 same-model reference: reverse Monte Carlo supercells at target alpha (new tool next to `tools/sqs.py`), k-space with one potential per species | 🔴 |
| SRO between two sublattices (for example L2_1) in an SCF run | 🔴 only the reader and `sro_conditional` are tested |

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
1. Done 2026-10-30: reviewed and committed; the four unit tests pass.  Main is not pushed yet (the user pushes).
2. Done 2026-10-31: the `hull` reruns are in `~/asr_remote4/trial/sro_hull`.  All converged, G1 passes, and `hull` agrees with `bp`
   (FINDINGS_LOG "SRO reruns with the `hull` tail").  Use the `hull` runs from now on.
3. B4.4 done 2026-10-31 at levels 3 and 4 (`trial/sro_hull`, `trial/sro_l4`; FINDINGS_LOG "Gate B4.4").  Open: the factor 10 in
   the size of the SRO effects.  It does not shrink from levels 3 to 4.  With the paper's A and B swapped, the grouping of the triads agrees, but the
   paper's Fe-Cr SRO peak is at +0.14 Ry and ours at -0.06 Ry, and the paper's is about ten times larger (FINDINGS_LOG "Alam 2009
   figures 4 and 5 against ours").  The paper's DOS is per spin.  Its Cr band sits at least 0.1 Ry above ours, so the figure 2
   comparison of the component DOS (HANDOFF item 2) comes first.
4. B4.3: the SRO-induced change of the CuZn DOS against Tarafder 2006 Fig. 3 (`../references/SRO-cuzn_tarafder2006.pdf`).
5. B3: put the two options (SIM unchanged vs explicit shell charge) to the user with numbers from the triads and CuZn.  An explicit
   shell term partly double counts the SIM screening charge.
6. B4.5: reverse Monte Carlo tool for supercells at target alpha, then k-space references with one potential per species.
7. Record every result in `docs/FINDINGS_LOG.md` and update the SRO rows of `docs/STATUS.md` and `docs/ROADMAP.md`.

## Rules that matter here
- All calculations on the desktop; unit tests on the laptop.  Code moves to the desktop by git (bundle or the user's push).
- Production route `asr_ed` only.  Default tail `hull`; do not use `bp` for new runs (FINDINGS_LOG "Levels-4 SCF stall").
- A fix or a model choice must hold in general, not only for the test system (memory `general-fixes-not-patches`).
- Docs in simple sentences, no dashes.
