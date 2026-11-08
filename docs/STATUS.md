# STATUS

Updated 2026-11-08.  🟢 done and verified · 🟠 in progress or partly verified · 🔴 not done · ⛔ blocker.
Only verified claims are marked 🟢.  Numbers, tables and history: `docs/FINDINGS_LOG.md` (section names in brackets).
Plan and order of work: `docs/ROADMAP.md`.  Where things are and how to run them: `docs/HANDOFF.md`.

## Production settings
Use these for new runs (energy-dependent ASR, `method = "asr"`):
| key | value | why |
|---|---|---|
| `[asr] levels` | 4 (3 for scans and HEAs) | HEAs meet the bar at levels 3 with 24 steps [HEA per-site reference] |
| `rank` | 1 | rank 2 and 3 change nothing at levels 4 [Queue 2 results] |
| `recursion_steps` | 24 for strong moments, 36 for Ni-type weak moments | chain resolution [Tail scan at fixed potential] |
| `ed_terminator` | `"hull"` (default since 3a47711) | `bp` puts a pole next to the contour [Levels-4 SCF stall] |
| `fermi_kt` | 0.001 | removes the N(E) noise floor, within 0.003 muB of sharp |
| `contour_points` | 32 | 32 = 64 to 1e-3 |
| `[scf] sim_beta` | 1 | corrected SIM [SIM sign and magnitude error] |
| `mixing`, `mixing_restart` | 0.2, 0.03 | stiff charge-transfer modes |
| `xc` | `"vbh"` (LDA) | the references use LDA |

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
| terminators `hull` (default), `bp`, `avg2` | 🟢 `hull`, `avg2` · 🟠 `bp` | `bp` is not safe on complex chains; `test_terminators` checks the Nevanlinna bound \|dG/dz\| <= -Im G / Im z [Levels-4 SCF stall] |
| short-range order, first shell (`asr_sro`, `[[sro]]`) | 🟠 | committed `176a727`; `test_sro` 16/16, `test_augcls` 16/16, SCF gate G1 [Short-range order]; plan in `docs/HANDOFF_sro.md` |
| k-space LMTO-ASA for ordered cells and SQS supercells | 🟢 | SQS references below |
| cross-check routes (EI H^(2)/H^(3), P(z) - S, exact KKR-ASA) | 🟢 kept for tests | not for production |
| N(E) probe `ef_scan`, chain dump `ef_scan_dump`, `bp` fallback counter | 🟢 | merged in 099ae3a |

## Validation against same-model references (3 to 4 % bar)
| system | result | status |
|---|---|---|
| bcc Fe, fcc Ni, B2 FeCo, B2 CuZn ordered | ED route = same-model KKR-ASA to < 1 mRy, 0.01 muB | 🟢 |
| bcc FeCo, FeV; fcc FeNi, NiCu (36 steps) vs SQS (one potential per species) | within 0.2 to 5 % | 🟢 |
| bcc FeRh | Fe 4.8 % low, flat in steps | 🟠 |
| L2_1 Co2FeAl, Co2CrAl, Co2Fe0.4Cr0.6Al | Fe, Cr, Al within 4 %; Co 5 to 8 % low (36 steps) | 🟢 Co explained (chain resolution) |
| fcc Fe66Ni14Cr20 levels 3 NM, spin; levels 4 NM | converged; l4 needs `hull` | 🟢 |
| CoCrFeNi levels 3 spin (hull), 24 steps, vs 40-atom SQS with one potential per site, H^(2) | Co -4.2 %, Fe -2.8 %, Ni -3.7 %, Cr +22 %; one cell gives about 10 % sampling error per species | 🟢 Co, Fe, Ni; 🟠 Cr needs more cells |
| CoCrFeNi per-site cell, full ASA against H^(2) | moments change by less than 3 %, Ni by 7 %; `ed3` is not needed for random HEA moments | 🟢 |
| CoCrFeNi levels 3 at 36 and 48 steps | moments rise past the reference (36 steps: Co +7 %, Fe +8 %, Cr +40 %); no convergence in steps at fixed levels 3 | 🟠 cause open |
| CoCrFeMnNi levels 3 spin (hull) vs the per-species 40-atom SQS | Fe -10 %, Co -28 %, Ni -13 %, Cr +7.5 %, Mn -2.3 %; that reference is superseded, the per-site cell is not run | 🟠 |
| Alam and Mookerjee 2009 Fig. 2 (Fe66Ni14Cr20 DOS) | levels 3 shapes agree, Ni-d peak 0.08 Ry higher | 🟠 levels 4 comparison open |

## Blockers and open problems
| item | status |
|---|---|
| cost: CoCrFeNi levels 4 takes 2.7 h per iteration, CoCrFeMnNi more | ⛔ for the SRO series and the Fe-Cr series; HEA levels 4 and 5 are low priority now |
| HEA reference: the per-species SQS has charges about 1.4 times smaller than the per-site cell, so it is not a reference for charges or moments | 🟢 found; 🔴 more per-site cells |
| SRO electrostatics (B3): charge-excess model (`64bda88`, beta 1, no fit) gives SRO energy within 5 to 12 % and charge change within 10 to 15 % of per-atom supercells; fit there gives beta 1.0, lambda 0.21 as predicted | 🟢 energies and charges; 🟠 band separation shifts 1.45 times too much |
| SRO gates B4.3 (Tarafder CuZn), B4.4 (Alam Fig. 4 and 5), B4.5 (RMC supercells) | 🟠 B4.4 done (the size gap comes from Alam's random alloy); B4.3 done (B2 limit and main Zn peak agree, Cu band differs); B4.5: DOS shape agrees, energies and band positions do not |
| full unit-test suite on merged main `3986fc2` (desktop, 2026-10-30, about 75 min): 41 test programs, 306 checks, 0 failed; 6 `probe_*` programs stop because they need arguments or local reference files (not tests) | 🟢 |
| gfortran + OpenBLAS build never tried | 🔴 (CPC requirement) |
| FeRh Fe -4.8 % (rank 2 untested) | 🟠 low priority |

## Known discrepancies (documented, not blocking)
- Ag parameters: 10 to 30 mRy absolute offsets vs Stuttgart (E_nu, C, V(S)) [Self-consistent ASA fcc Ag].
- bcc Fe moment 1 % below Questaal (2.252 vs 2.275, VWN) [Our LMTO-ASA vs Questaal].
- Weak moments at an edge near E_F (Ni, Co in half-metals) are low at 24 steps: chain energy resolution, about 1/m [Tail scan].  This does not hold for the HEAs at levels 3: against the per-site cell they meet the bar at 24 steps, and more steps overshoot [HEA per-site reference].
- Disordered MnCr and FeV references from the old ASR papers are not reproduced; our own SQS are the gates [MnCr vs Banerjee and
  Mookerjee; One potential per species SQS].
- `asr_alloy_scf` P(z) - S route: unstable antisymmetric charge mode in its A = B test (cross-check route only).
- Total energy: Kohn-Sham form, Fermi smearing entropy not included; checked on fcc Cu (a0 3.546 A, B 199 GPa).

## Cost (one iteration, desktop, 28 cores shared)
| system | levels | nodes | wall per iteration | memory |
|---|---|---|---|---|
| Fe66Ni14Cr20 NM | 4 | 64,567 | 11 min (28 threads) | 1.5 to 3 GB |
| CoCrFeNi | 3 / 4 | 23,681 / 96,763 | 7 to 12 min (14 threads) / 2.7 h | 0.6 / 2.3 GB |
| CoCrFeMnNi | 3 / 4 | 31,535 / 128,959 | 11 to 12 min (14 threads) / not measured | 0.8 / about 3 GB |
The operator build is serial: 16 min for the quinary at levels 4.
