# STATUS

Updated 2026-11-19.  `main` = `9383b29` (streamline, output-data, spectral window, parallel teams and hot loops merged).  The review fixes
(`docs/REVIEW_fable.md`) and the F4 analysis are on `main`.  The only other branch is `epi-sro` (unmerged).
🟢 done and verified · 🟠 in progress or partly verified · 🔴 not done · ⛔ blocker.
Only verified claims are marked 🟢.  Numbers, tables and history: `docs/FINDINGS_LOG.md` (section names in brackets).
Plan and order of work: `docs/ROADMAP.md`.  Where things are and how to run them: `docs/HANDOFF.md`.

## Production settings
Use these for new runs (energy-dependent ASR, `method = "asr"`):
| key | value | why |
|---|---|---|
| `[asr] levels` | 4 (3 for scans and HEAs) | HEAs meet the bar at levels 3 with 24 steps [HEA per-site reference] |
| `rank` | 1 | 🟠 the old evidence was 4-step chains [Queue 2 results].  At 24 steps and fixed potential (review 6.4): CuZn levels 3 no effect; FeNiCr levels 2 rank 1 is 0.6 / 1.6 / 2.4 % high in the Fe / Ni / Cr moments, rank 2 = rank 3; FeNiCr levels 3: rank 2 changes the moments by -0.5 / +0.2 / -2.6 % and the Cr charge by 12 %.  Self-consistent FeNiCr levels 3 (review 6.4, 2026-11-13): rank 2 lowers the moments by 6.3 / 4.0 / 11.5 % (Fe / Ni / Cr); charges equal to rank 1 within 0.0014 e.  Rank 1 is not converged for magnetic multicomponent alloys (outside the 3 to 4 % target).  CoCrFeNi levels 3 rank 2 SCF (`desk4/cfcn_L3_rank2`, iteration 17, dV 1.9e-5): Co 0.883, Cr -0.519, Fe 1.676, Ni 0.259 (rank 1: 0.914 / -0.619 / 1.751 / 0.266).  One rank-2 iteration plus the linear-response predictor gave 0.882 / -0.504 / 1.680 / 0.257 before the run ended (`docs/DERIVATION_scf_response.md` part C).  Rank 2 at levels 3: 500,607 nodes, 1.16 h per iteration on 20 threads, 15.6 GB.  Levels-2 SCF pairs (review 6.4, 2026-11-14) reproduce the levels-3 effect for FeNiCr (-6.3 / -4.6 / -13.7 %) at 1/8 of the cost; CoCrFeNi levels 2: Co -5.2, Cr -17.9, Fe -3.8, Ni -4.3 %.  bcc Fe0.875Cr0.125 levels 3: Fe -0.03 %, Cr -1.1 %, so rank 1 suffices there.  Charges converged at rank 1 in all cases (0.0026 e) |
| `recursion_steps` | 24 for strong moments, 36 for Ni-type weak moments | chain resolution [Tail scan at fixed potential].  🟠 not converged.  At a fixed potential the moments converge with the steps towards the Green function of the node set: flat from 64 (FeNiCr) to 96 (CoCrFeNi) steps; FeNiCr levels 2 and 3 have the same limit within 1 %.  At 24 steps the termination error is a few per cent: FeNiCr scatters by 3 to 5 % (Ni 14 %) without a trend, CoCrFeNi is 4 to 9 % low.  The SCF multiplies it by the softness of the moments: FeNiCr levels 2 SCF Fe 1.679 (24 steps), 1.471 (64), 1.515 (96) [`docs/DERIVATION_levels_steps.md`].  A step-converged default is not chosen yet |
| `ed_terminator` | `"hull"` (default since 3a47711) | `bp` puts a pole next to the contour [Levels-4 SCF stall] |
| `fermi_kt` | 0.001 | removes the N(E) noise floor, within 0.003 muB of sharp |
| `contour_points` | 32 | 32 = 64 to 1e-3 |
| `[scf] sim_beta` | 1 for bcc CuZn; 1.2 for fcc CoCrFeNi | corrected SIM [SIM sign and magnitude error].  The value depends on the alloy: per-atom cells give 1.0 (CuZn) and 1.2 (CoCrFeNi, review 6.6); the charge-excess model predicts both from the measured hardness.  Earlier HEA runs used 1.  With strong clustering SRO above about 1, the screened potential depends on the further-shell order, which the input does not fix [B3 closure] |
| `mixing`, `mixing_restart` | 0.2, 0.03 | stiff charge-transfer modes.  With `[[sro]]` the SRO potential is iterated in its own variable (the drive u = db, saved in `pot.dat`).  Since `1dfe1bc` u starts at zero when the start potential does not hold it; the SRO example converges from free atoms with 0.2 and no `mixing_restart` [B3 closure, the start of the drive] |
| `xc` | `"vbh"` (LDA) | the references use LDA |
| `route` | `"edc"` (default); `"ed3"` (full ASA) for charges and SRO energies, with `cluster_radius` 2.0 | ed3 at rcl 2.0 matches the per-atom SRO cells best [CuZn SRO, ed3 rcl 2.0] |
| `energy_block` | 0 (default since e8d193c) | blocks: 1.35x for CoCrFeNi but slower for a random binary; B2 ed3 at 96 steps differs by 1e-4 e between blocks 0 and 8, reproducibly with the review binary; not removed by F2 (review doubt 7) |
| `ef_linear_max` | 0 (default); 0.005 for long magnetic SCFs | FeNiCr spin 1.64x faster, moments within 1.5e-3 muB |

## Capabilities
| area | status | evidence |
|---|---|---|
| recursion core, ordered TB DOS vs k-space | 🟢 | `test_recursion_dense`, `test_ordered_dos` |
| own LMTO-ASA (radial, libxc, atom SCF, potential parameters, spin, Madelung, tetrahedra) | 🟢 | `test_asa_ag`, `test_asa_mag`, `test_asa_b2`; Questaal comparison [Our LMTO-ASA vs Questaal] |
| augmented space: binary, multi-orbital, N species, site classes; D4h class reduction | 🟢 | `test_augspace`, `test_augorb`, `test_augn`, `test_augcls` (exact vs brute force to 1e-12 or better) |
| energy-dependent ASR in AMscf form (`asr_ed`, `asr_edn`, `asr_edc`) | 🟢 | `test_asr_ed`, `test_asr_ed_n`, `test_asr_ed_sym` |
| class SCF: any sublattices and species, Madelung + SIM, exact electron count, Anderson with restart fix, total energy | 🟢 | binaries, Heuslers, HEAs below |
| OpenMP over independent chains | 🟢 | bit identical, l3 iteration 867 s to 398 s on 14 threads |
| one input format and driver, `asr.out` report | 🟢 | `test_run_input` |
| terminators `hull` (default), `bp`, `avg2` | 🟢 `hull`, `avg2` · 🟠 `bp` | `bp` is not safe on complex chains; `test_terminators` checks the Nevanlinna bound \|dG/dz\| <= -Im G / Im z [Levels-4 SCF stall]; a terminated energy-dependent chain is not causal in general (review F4).  Analysed 2026-11-13 (`docs/DERIVATION_f4.md`, on main since `bed1d80`): no tail is causal by construction; a non-causal pole disturbs the SCF only within about 10 kT of E_F; the code now counts the non-causal node values for `edc` too |
| short-range order, first shell (`asr_sro`, `[[sro]]`) | 🟢 | `test_sro` 16/16, `test_augcls` 16/16, SCF gate G1 [Short-range order]; all seven FeNiCr triads converge, including (1,1,1).  From free atoms since `1dfe1bc` (before: diverged; the first drive kept the free-atom charges of 5 e) [B3 closure, the start of the drive] |
| SRO electrostatics without the pole: minimal-tail closure of the first-shell structure factor, SCF in the drive u = db (`asr_screening`, merged `0af30ad`; drive starts at zero, `1dfe1bc`) | 🟢 | `test_screening` 70/70; exact charge-excess model on Monte Carlo supercells of 131072 sites (`tools/sro_closure`): within 1 % for random further shells [B3 closure] |
| k-space LMTO-ASA for ordered cells and SQS supercells | 🟢 | SQS references below |
| cross-check routes (EI H^(2)/H^(3), P(z) - S, exact KKR-ASA) | 🟢 kept for tests | not for production |
| N(E) probe `ef_scan`, chain dump `ef_scan_dump`, `bp` fallback counter | 🟢 | merged in 099ae3a |
| step profile from one chain dump (`tools/step_profile.py`): charges and moments for every chain length up to `recursion_steps`, `hull` tail, E_F set per length | 🟢 | the 24-step value equals the 24-step run to 4e-4 muB |
| linear response of the SCF from one chain dump (`tools/scf_response.py`): thermal DOS, on-site response, species response with single-site vertex weights, Stoner I from `asr.out`; predicts the self-consistent change from one fixed-potential iteration | 🟢 | FeNiCr levels 2, 3 (rank 2; levels 3 to 4) within 1.5 points of the moment in all nine entries, including the Ni sign reversal; CoCrFeNi levels 2, 3 within 1 point for Co, Fe, Ni, Cr 2 to 3 points [`docs/DERIVATION_scf_response.md` part C] |
| full generalised ASA route `ed3` (`asr_ed3`, two chains per energy) | 🟢 | gates 1 to 4b [DERIVATION_ed3]; merge check on 9 examples (dos identical, potentials to 1.6e-14) |
| block kernels over contour energies (`asr_augblk`, `energy_block`) | 🟠 | `test_augcls_blk`; random edc bitwise equal for blocks 0/2/4/8; B2 ed3 differs by 1e-4 (E_F -0.16219 vs -0.16214): with the review binary (F2, F3) 0.12916 vs 0.12928 e, bit-reproducible; a rounding sensitivity of long ordered chains, not a kernel bug; cause of the size open (review doubt 7) |
| one contour per iteration for small E_F steps (`ef_linear_max`) | 🟢 opt-in | CuZn 90 to 84 passes, same qnet; FeNiCr spin 84 to 66 passes, 32000 s to 19539 s |
| MPI over the chain tasks (`asr_mpi`, build with `-DUSE_MPI`) | 🟢 | `test_mpi`; serial vs 1/2/4 ranks bitwise for CuZn B2, random edc and ed3; spin FeNiCr SRO within 3.4e-14 (serial run-to-run spread 7.7e-15, caused by `orbital_rep`, fixed in review F3); node scaling below |

## Validation against same-model references (3 to 4 % bar)
| system | result | status |
|---|---|---|
| bcc Fe, fcc Ni, B2 FeCo, B2 CuZn ordered | ED route = same-model KKR-ASA to < 1 mRy, 0.01 muB | 🟢 |
| bcc FeCo, FeV; fcc FeNi, NiCu (36 steps) vs SQS (one potential per species) | within 0.2 to 5 % | 🟢 |
| bcc FeRh | Fe 4.8 % low, flat in steps | 🟠 |
| L2_1 Co2FeAl, Co2CrAl, Co2Fe0.4Cr0.6Al | Fe, Cr, Al within 4 %; Co 5 to 8 % low (36 steps) | 🟢 Co explained (chain resolution) |
| fcc Fe66Ni14Cr20 levels 3 NM, spin; levels 4 NM | converged; l4 needs `hull` | 🟢 |
| fcc Fe66Ni14Cr20 SRO triads (Alam and Mookerjee), levels 3 NM, closure (`desk9`) | all seven converge (the old form diverged for (1,1,1)); a species with only its own kind as neighbours is exactly neutral; the other charges agree with the old form to 0.003 e.  From free atoms (`1dfe1bc`): levels 2, all seven (`desk18`), charges equal to `desk9` within 0.003 e; levels 3, the example (1,1,0) in 40 iterations and (1,1,1) in 64, both equal to `desk9`; the other five from free atoms in `desk24` (29 to 43 iterations), energies equal to `desk9` within 5e-7 Ry, charges within 1e-5 e | 🟢 [B3 closure]; 🟢 levels 3 from free atoms, all seven |
| CoCrFeNi levels 3 spin (hull), 24 steps, vs 40-atom SQS with one potential per site, H^(2) | sim_beta 1: moments Co -4.2 %, Fe -2.8 %, Ni -3.7 %, Cr +22 %, but charges 7 to 12 % low.  sim_beta 1.2 (the cell's value, review 6.6): charges within 1 to 6 %, moments Co -5.5 %, Fe -2.3 %, Ni -5.6 %, Cr +26 %.  One cell gives about 10 % sampling error per species | 🟠 the beta-1 agreement was partly compensation; the moment error is in the recursion (rank, levels, steps), not the electrostatics |
| CoCrFeNi per-site cell, full ASA against H^(2) | moments change by less than 3 %, Ni by 7 %; `ed3` is not needed for random HEA moments | 🟢 |
| CoCrFeNi levels 3 at 36 and 48 steps | moments rise past the reference (36 steps: Co +7 %, Fe +8 %, Cr +40 %).  Cause: at the fixed potential the 24-step moments are 4 to 9 % below the step limit (levels 2, 192 steps: Co 1.005, Cr -0.646, Fe 1.817, Ni 0.300), and the SCF gain is about 4.  The step-converged rank-1 moments at that potential lie above the per-site cell for Co (+4 %), Cr (+32 %) and Ni (+6 %); rank 2 lowers them by 3 to 18 %.  The 24-step agreement rests on two errors of opposite sign [`docs/DERIVATION_levels_steps.md`] | 🟢 cause found; 🔴 step-converged rank-2 run |
| bcc CuZn SRO, `ed3` rcl 2.0, converged (smriti, 4 nodes) vs per-atom cells.  Check 1 is done: the cell reruns `a0_s1`, `a-0.30_s1` with the fixed code give the table values to 3e-7 Ry, so the corrected cell energies stand.  Check 2 waits in the smriti queue: `ed3` levels 4 at +-0.148.  The `ed3` entries at -0.266 predate the SRO closure and need a rerun | q_Cu (random) 0.1023 vs 0.102.  dE (mRy) at -0.148 / +0.148 / -0.266.  **Corrected form** (on main since F1 `bca3a15`): ed3 levels 3 -1.14 / +1.27 / -2.29, levels 4 at -0.266 -1.73; cells -1.07..-1.10 / +1.10 / -1.77.  **Old form** (E_band = q1 + E_nu q0, before F1): ed3 levels 3 -1.33 / +1.53 / -2.64, levels 4 at -0.266 -2.11; cells -1.44..-1.49 / +1.48 / -2.30.  edc levels 4 (not affected by F1) with the closure (`desk9`): -1.30 / +1.53 / -2.10, dq_Cu +0.010 / -0.011 / +0.015 (before the closure -1.28 / +1.55 / -2.16).  dq_Cu (both forms): ed3 +0.0126 / -0.0130 / +0.0225 (levels 4 at -0.266 the same), cells +0.012 / -0.012 / +0.019.  Band-centre shift at -0.148: ed3 +9.5 mRy vs cells +6.5 (edc +9.1) | 🟠 energies.  Old form: ed3 levels 3 -8 to -11 % / +3 % / +15 %, levels 4 -8 %.  Corrected form: ed3 levels 3 +4 to +7 % / +15 % / +30 %, levels 4 at -0.266 -2 %.  🟢 dq_Cu at +-0.148; 🟠 dq_Cu at -0.266 overshoots 18 % at levels 3 and 4 (that cell has shell-2/3 correlations 0.07/0.05); band shift 1.45x as edc |
| CoCrFeMnNi levels 3 spin (hull) vs the per-species 40-atom SQS | Fe -10 %, Co -28 %, Ni -13 %, Cr +7.5 %, Mn -2.3 %; that reference is superseded, the per-site cell is not run | 🟠 |
| Alam and Mookerjee 2009 Fig. 2 (Fe66Ni14Cr20 DOS) | levels 3 shapes agree, Ni-d peak 0.08 Ry higher | 🟠 levels 4 comparison open |

## Blockers and open problems
| item | status |
|---|---|
| cost: CoCrFeNi levels 4 takes 2.7 h per iteration, CoCrFeMnNi more | ⛔ for the HEA SRO series; HEA levels 4 and 5 are low priority now.  Not for binaries: bcc Fe0.875Cr0.125 levels 4 with 48 steps converges in about 80 min (`fecr/asr/L4s48_*`) |
| HEA reference: the per-species SQS has charges about 1.4 times smaller than the per-site cell, so it is not a reference for charges or moments | 🟢 found; 🔴 more per-site cells |
| SRO electrostatics (B3): charge-excess model (`64bda88`, beta 1, no fit) gives SRO energy within 5 to 12 % and charge change within 10 to 15 % of per-atom supercells; fit there gives beta 1.0, lambda 0.21 as predicted.  The first-shell form had a pole for strong clustering (FeNiCr triads blew up); the minimal-tail closure removes it and keeps the first-shell form wherever that is realizable [B3 closure] | 🟢 charges; 🟢 pole removed; 🟠 energies: `edc` 19 to 39 % above the corrected cells (row "bcc CuZn SRO").  The per-atom H^(2) cells give -2.37 mRy with `sigma` 0.7 (`desk4/h2cells`) and -2.34 mRy with `std_alpha = true` (`desk25`); the representation moves the SRO energy by 1 %, so `edc` (-2.10) is 10 % below its own-Hamiltonian cell and the excess over the full-ASA cells is the H^(2) Hamiltonian; 🟠 band separation shifts 1.45 times too much |
| SRO gates B4.3 (Tarafder CuZn), B4.4 (Alam Fig. 4 and 5), B4.5 (RMC supercells) | 🟠 B4.4 done (the size gap comes from Alam's random alloy); B4.3 done (B2 limit and main Zn peak agree, Cu band differs); B4.5: DOS shape agrees, energies and band positions do not |
| unit tests on `sro-drive-start` = main `1dfe1bc` (laptop, ifx, 2026-11-16/17): `fpm test` runs 38 programs, 374 checks, 0 failed.  13 slow tests (SCF, enumeration; 121 checks) are fpm examples in `test_slow/`, run by `asr_test_slow`; all passed on `95fe365`; after the drive change only `test_scf_guards` (10/10, the one slow test that runs the class SCF) was rerun.  The 17 diagnostic `probe_*` programs are fpm examples in `probe/` and take arguments | 🟢 |
| gfortran 14 + OpenBLAS 0.3.28 build (smriti, CentOS 7): 6 examples equal to ifx/MKL to all printed digits; unit tests pass except the `test_augcls` three-class ring check (err 18.7).  Cause: the test's own configuration count under gfortran 14, not OpenBLAS; fixed with d54f018 (unallocated species parameters) in review 3.6; the ring check passes with gfortran 14.  All LAPACK calls have explicit interfaces (`asr_lapack`, `651fb66`): gfortran 16 with -Werror=implicit-interface compiles all 155 sources | 🟢 (on main; CPC requirement) |
| MPI shutdown crash: one of two production runs (m015) ended with a segfault 3 s after writing all output; a rerun did not repeat it | 🟠 watch |
| FeRh Fe -4.8 % (rank 2 untested; rank is now a live suspect, review 6.4) | 🟠 |
| SCF guards (review 6.5, `36c10b1`): a run with E_nu held at its band window in the last iteration is reported as not converged; valence weight below the contour bottom is measured and flagged above 1e-3 e; `test_scf_guards` 10/10 | 🟢 |
| bcc Fe0.875Cr0.125 SRO: ASR runs at alpha = 0, +-0.10 (levels 3, 4; 48 steps) and per-atom k-space cells at alpha = 0, +-0.1, two seeds each (desktop `~/asr_remote4/trial/fecr`, 5 to 9 Nov) | 🔴 comparison not done |
| `example/fe66ni14cr20_fcc.toml` says that the non-magnetic SCF does not converge at levels 4 (it alternates).  This file says it converged with `hull`, as do ROADMAP and FINDINGS_LOG [Levels-4 SCF stall].  The run `desk11/rand_fa` was stopped by its cap in the first iteration, so it does not decide this | 🟠 source of the comment unknown |

## Known discrepancies (documented, not blocking)
- H^(2) depends on the screening constants at third order.  The recursion uses `standard_alpha()`; `[kspace] second_order` uses the
  hard-core alpha of `sigma` unless `std_alpha = true`.  B2 CuZn (a = 2.85 A): k-space H^(2) moves by 3.0 mRy per atom between the two
  sets.  The recursion was 2.2 mRy above k-space with the default set and is 0.3 mRy below it (nk 32) with `std_alpha`.  Steps,
  levels, supercell and cluster radius move the recursion energy by less than 0.5 mRy [B2 CuZn: the 2.2 mRy between recursion and k-space].  Every
  internal check of H^(2) recursion against k-space needs `std_alpha = true`.
- Ag parameters: 10 to 30 mRy absolute offsets vs Stuttgart (E_nu, C, V(S)) [Self-consistent ASA fcc Ag].
- bcc Fe moment 1 % below Questaal (2.252 vs 2.275, VWN) [Our LMTO-ASA vs Questaal].
- Weak moments at an edge near E_F (Ni, Co in half-metals) are low at 24 steps: chain energy resolution, about 1/m [Tail scan].  This does not hold for the HEAs at levels 3: against the per-site cell they meet the bar at 24 steps, and more steps overshoot [HEA per-site reference].
- Disordered MnCr and FeV references from the old ASR papers are not reproduced; our own SQS are the gates [MnCr vs Banerjee and
  Mookerjee; One potential per species SQS].
- `asr_alloy_scf` P(z) - S route: unstable antisymmetric charge mode in its A = B test (cross-check route only).
- Total energy: Kohn-Sham form, Fermi smearing entropy not included; checked on fcc Cu (a0 3.546 A, B 199 GPa, route not recorded).
  Before F1 (`bca3a15`, now on main) the k-space full ASA and `ed3` band energies lacked E_nu p q2 per channel (about -77 mRy per
  CuZn site).  fcc Cu k-space (nk 20 tetrahedra): a0 3.523 A, B 172 GPa without F1; 3.563 A, 185 GPa with it.  The 3.546 A /
  199 GPa above is not reproduced by either binary, so it came from another route or other settings.  `edc` is not affected.

## Cost (one iteration, desktop, 28 cores shared)
| system | levels | nodes | wall per iteration | memory |
|---|---|---|---|---|
| Fe66Ni14Cr20 NM | 4 | 64,567 | 11 min (28 threads) | 1.5 to 3 GB |
| CoCrFeNi | 3 / 4 | 23,681 / 96,763 | 7 to 12 min (14 threads) / 2.7 h | 0.6 / 2.3 GB |
| CoCrFeMnNi | 3 / 4 | 31,535 / 128,959 | 11 to 12 min (14 threads) / not measured | 0.8 / about 3 GB |
The operator build is serial: 16 min for the quinary at levels 4.

Cluster (smriti, 48 cores per node, MPI one rank of 48 threads per node), random bcc CuZn `ed3` levels 3, rcl 2.0:
| nodes | one iteration from the converged potential (s) | speed-up | efficiency |
|---|---|---|---|
| 1 | 4811 | 1.00 | 1.00 |
| 2 | 2588 | 1.86 | 0.93 |
| 4 | 1551 | 3.10 | 0.78 |
| 8 | 1053 | 4.57 | 0.57 |
Full SCF on 4 nodes: random 14 iterations in 5.9 h, m015 15 iterations in 6.1 h (desktop: 19.4 h for one iteration on 6 threads).
