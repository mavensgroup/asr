# STATUS

Updated 2026-10-30.  🟢 done and verified · 🟠 in progress or partly verified · 🔴 not done · ⛔ blocker.
Only verified claims are marked 🟢.  Numbers, tables and history: `docs/FINDINGS_LOG.md` (section names in brackets).
Plan and order of work: `docs/ROADMAP.md`.  Where things are and how to run them: `docs/HANDOFF.md`.

## Production settings
Use these for new runs (energy-dependent ASR, `method = "asr"`):
| key | value | why |
|---|---|---|
| `[asr] levels` | 4 (3 for scans) | Co, Ni moments in HEAs need 4 or 5 |
| `rank` | 1 | rank 2 and 3 change nothing at levels 4 [Queue 2 results] |
| `recursion_steps` | 24 for strong moments, 36 for Ni-type weak moments | chain resolution [Tail scan at fixed potential] |
| `ed_terminator` | `"hull"` (the code default is still `"bp"`) | `bp` puts a pole next to the contour [Levels-4 SCF stall] |
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
| terminators `hull`, `bp`, `avg2` | 🟢 `hull`, `avg2` · 🟠 `bp` | `bp` is not safe on complex chains [Levels-4 SCF stall] |
| short-range order, first shell (`asr_sro`, `[[sro]]`) | 🟠 | uncommitted; `test_sro` 16/16, `test_augcls` 16/16, SCF gate G1 [Short-range order] |
| k-space LMTO-ASA for ordered cells and SQS supercells | 🟢 | SQS references below |
| cross-check routes (EI H^(2)/H^(3), P(z) - S, exact KKR-ASA) | 🟢 kept for tests | not for production |
| N(E) probe `ef_scan`, chain dump `ef_scan_dump` | 🟢 | branch `l4stall`, not merged |

## Validation against same-model references (3 to 4 % bar)
| system | result | status |
|---|---|---|
| bcc Fe, fcc Ni, B2 FeCo, B2 CuZn ordered | ED route = same-model KKR-ASA to < 1 mRy, 0.01 muB | 🟢 |
| bcc FeCo, FeV; fcc FeNi, NiCu (36 steps) vs SQS (one potential per species) | within 0.2 to 5 % | 🟢 |
| bcc FeRh | Fe 4.8 % low, flat in steps | 🟠 |
| L2_1 Co2FeAl, Co2CrAl, Co2Fe0.4Cr0.6Al | Fe, Cr, Al within 4 %; Co 5 to 8 % low (36 steps) | 🟢 Co explained (chain resolution) |
| fcc Fe66Ni14Cr20 levels 3 NM, spin; levels 4 NM | converged; l4 needs `hull` | 🟢 |
| CoCrFeNi levels 3 spin (hull) vs 40-atom SQS | Fe -3.6 %, Co -9.6 %, Ni -10 %, Cr 27 % too large | 🟠 Co, Ni need levels 4 or 5 |
| CoCrFeMnNi levels 3 spin (hull) vs 40-atom SQS | Fe -10 %, Co -28 %, Ni -13 %, Cr +7.5 %, Mn -2.3 % | 🟠 same |
| Alam and Mookerjee 2009 Fig. 2 (Fe66Ni14Cr20 DOS) | levels 3 shapes agree, Ni-d peak 0.08 Ry higher | 🟠 levels 4 comparison open |

## Blockers and open problems
| item | status |
|---|---|
| default tail is still `bp` in the code; `l4stall` not merged | ⛔ for every new run until fixed (set `ed_terminator = "hull"`) |
| cost: CoCrFeNi levels 4 takes 2.7 h per iteration, CoCrFeMnNi more | ⛔ for HEA levels 4 and 5, SRO series, Fe-Cr series |
| SRO electrostatics (B3): SIM with a correlated shell, no literature prescription | 🔴 decision needed |
| SRO gates B4.3 (Tarafder CuZn), B4.4 (Alam Fig. 4 and 5), B4.5 (RMC supercells) | 🔴 |
| full unit-test suite: last complete run 2026-10-09 (39 programs, desktop); later threads ran subsets; full run takes more than 20 min | 🟠 |
| gfortran + OpenBLAS build never tried | 🔴 (CPC requirement) |
| FeRh Fe -4.8 % (rank 2 untested) | 🟠 low priority |

## Known discrepancies (documented, not blocking)
- Ag parameters: 10 to 30 mRy absolute offsets vs Stuttgart (E_nu, C, V(S)) [Self-consistent ASA fcc Ag].
- bcc Fe moment 1 % below Questaal (2.252 vs 2.275, VWN) [Our LMTO-ASA vs Questaal].
- Weak moments at an edge near E_F (Ni, Co in half-metals and HEAs) are low at 24 steps: chain energy resolution, about 1/m [Tail scan].
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
