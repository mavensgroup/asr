# HANDOFF (next thread): l4 SCF stall, first HEA results, speed work (2026-10-28)

Continues `docs/HANDOFF_hea_a1.md`. Read that and `docs/HANDOFF_hea_sro.md` first (rules, Part A and B plan).
Production route is unchanged: energy-dependent asr_ed class SCF (memory `ed-is-target-route`).

## State of the code
- Local and GitHub `main` are at `574c467` (pushed on 2026-10-21). Commits of the 2026-10-21 thread:
  - `b727ff4` OpenMP over the independent chains (`chain_greens` / `chain_worker` in `src/asr_class_scf.f90`), Kohn-Sham total
    energy per site (printed every iteration, in `res%etot`, in `asr.out`), `save_potential` written every iteration, Anderson
    safeguard changes (`src/asr_mixing.f90`), `[scf] mixing_restart`.
  - `5c8075e` the first Anderson step (no history) also uses `mixing_restart`; test in `test/test_mixing.f90`.
  - `ef89271` `tools/sqs.py` handles N species (pair probabilities P_ss'(k) = x_s x_s'); the binary output is unchanged.
  - `574c467` the saved potential carries E_F (a restart starts its E_F search there), one E_F pass when the Newton step is
    below 1 mRy, per-iteration diagnostics (weighted residual norm next to dV, where the maximum dV sits, E_nu residual).
- Unit tests run: `test_mixing` (6/6) and `test_run_input` (26/26). The full `fpm test` was not run (more than 20 min).
- The user's own uncommitted files are untouched: `.gitignore`, `README.md`, `src/asr_kinds.f90`, `doc/`, `ford.md`,
  `suggestion.md`, `docs/JIJ_PLAN.md`, `docs/PUBLICATION_PLAN.md`.
- Desktop binaries: `~/asr_remote4` is at `5c8075e` (binary `build/ifx_A08D5D447AF67060/app/asr`). `~/asr_q` has `574c467`
  (binary `~/asr_q/bin/asr_574c467`). Pull with `git pull --ff-only github main` (https remote) or a git bundle over scp.

## Total energy
- The energy is the Kohn-Sham form with the output density in the input potential. The Fermi smearing entropy (kT 1 mRy) is not included.
- Check: fcc Cu, 5 lattice constants, levels 4: a0 = 3.546 A, B = 199 GPa (usual LDA overbinding, so the energy is sane).
- Near self-consistency the energy fluctuates by about 1e-4 Ry per site in the stalled l4 run, so it is only useful once the SCF converges.

## Open problem 1: Fe66Ni14Cr20 levels 4 does not converge
Runs (all in `~/asr_remote4/trial/a1/`):
- `l4` (binary 5c8075e): 20 iterations from the levels-3 potential, then 12 more from its own potential. dV stays near 4e-3 to 6e-3.
- `l4b` (binary 574c467, from the stall potential): 40 iterations, last dV 4.0e-3, not converged. Queue summary in
  `trial/hea/SUMMARY.txt`.
- The SCF does not crawl. It jumps between two states with no Anderson restarts and a full history of 6:
  state A: qnet Cr -0.185, E_total -2522.70758; state B: qnet Cr -0.200, E_total -2522.70745.
- Old failure (before `b727ff4`): 15 restarts in 35 iterations, each linear step multiplied dV by about 6. The restart fix removed that.

Diagnostics `a1/diag` (D1 to D4, each one iteration from the stall potential, first step = start + 0.03 r, so the residual is recovered
from the saved potential; script `diag_resid.py`):
- D1 and D2 (same input): residuals identical to 1e-13. The map is deterministic.
- D3 and D4: Cr potential shifted by 1e-5 and 1e-4 Ry. The response ratio 4 over 3 should be about 10 for a linear map. It is 5.5 (Fe),
  6.7 (Ni), 5.7 (Cr). So the map has a break at the 1e-5 to 1e-4 Ry scale. The largest residual sits at r = 0.7 to 0.8 Bohr (d region).
- `a1/diag0`: the same three runs with `depth_tol = 0` (all chains to 24 levels). **Identical numbers** (dV 6.007e-3, same ratios).
  **The adaptive recursion depth is ruled out.** One full-depth iteration takes about 38 min on 28 threads.

Remaining suspects, in order:
1. **E_F search path.** The old binary stops at |N - N_elec| < 2e-4 and then interpolates the moments linearly in the bracket. A small
   potential shift can change the number of passes and the bracket, so the interpolation error jumps. Test: rerun D1, D3, D4 with a
   tight stopping tolerance (for example 1e-7, or a fixed number of passes). The 574c467 binary has a different rule (one pass when the
   Newton step is below 1 mRy), so also run D1, D3, D4 with it and compare.
2. **Terminator fallback.** `cf_species` in `src/asr_endep.f90` (`bp_fixed_point`): if 200 fixed-point iterations or `zgeev` fail, it
   falls back to the hull tail without any message. Add a counter (like `ed_nchain`) and print it per iteration. Then run one iteration
   with `ed_terminator = "hull"` and compare the ratio.
3. Anderson itself: check whether a plain small linear step (mixing 0.02 to 0.05, no Anderson) converges slowly but surely. The stiff
   mode has gain about 8 (response/shift), so beta must be below 2/9 for plain mixing.
4. The Fermi contour (`fermi_kt` 0.001, 32 points): try 48 points in one diagnostic iteration.

## HEA first results (levels 3, 24 steps, ASR vs 40-atom SQS in k-space, one potential per species)
Inputs: `trial/hea/{cfn,cfmn}_{nm,sp}_l3/in.toml` (a = 3.57 A CoCrFeNi, 3.60 A CoCrFeMnNi), SQS in `trial/hea/sqs_*`.
Spin runs start from atomic potentials, Cr and Mn with `spin_sign = -1`, zeta0 0.1. All four ASR runs and both references converged.

| moment (muB) | ASR CoCrFeNi | SQS | ASR CoCrFeMnNi | SQS |
|---|---|---|---|---|
| Fe | 1.751 | 1.809 | 1.553 | 1.723 |
| Co | 0.934 | 1.025 | 0.622 | 0.854 |
| Ni | 0.274 | 0.302 | 0.139 | 0.157 |
| Cr | -0.599 | -0.471 | -0.213 | -0.214 |
| Mn | | | -0.587 | -0.599 |
| average per site | 0.590 | 0.666 | 0.303 | 0.384 |

- Meets the 3 to 4 % bar: Fe (CoCrFeNi), Mn and Cr (CoCrFeMnNi). Weak moments (Co, Ni) are 9 to 27 % low: chain resolution at
  levels 3 (memory `asr-resolution-physics`). Cr in CoCrFeNi is 27 % too large. Each SQS is one 40-atom cell, so it has its own scatter.
- Total energies (Ry/site): CoCrFeNi NM -2615.34917, spin -2615.35057 (spin lower by 1.40 mRy); CoCrFeMnNi NM -2555.05804,
  spin -2555.05872 (lower by 0.68 mRy).
- ASR cost: 10 to 13 h for 58 to 66 iterations (14 threads). k-space SQS: 40 to 65 min.
- Not extracted yet: the Cr d sum rule (open item A5). Iteration 1 of every run prints it, for example in `trial/hea/cost/*.log`.
  Values at E_F = 6 Ry (d column): CoCrFeNi Co 0.983, Cr 1.008, Fe 1.065, Ni 0.934; CoCrFeMnNi Co 0.970, Cr 1.019, Fe 1.047, Mn 1.026, Ni 0.926. Cr is not
  the outlier here, so A5 looks harmless for these alloys.
- Alam and Mookerjee comparison (Ni-Fe peak splitting, Cr spikes at levels 4) is still open because l4 does not converge.

## Cost table (one iteration from atomic start, partly with another job on the other 14 cores)
| system | levels | nodes | wall | memory |
|---|---|---|---|---|
| Fe66Ni14Cr20 | 4 | 64,567 | about 11 min later iterations (13 min with the first E_F passes) | 1.5 GB |
| CoCrFeNi | 3 / 4 | 23,681 / 96,763 | 35 min / 2.7 h | 0.6 / 2.3 GB |
| CoCrFeMnNi | 3 / 4 | 31,535 / 128,959 | 63 min / not finished | 0.8 / about 3 GB |

The operator build is serial. For the quinary at levels 4 it takes about 16 min once per run.

## Next steps (suggested)
1. Run the two cheap tests of suspects 1 and 2 on the idle desktop (D1, D3, D4 style, about 40 min each). Fix what they show.
2. Finish l4, then close the Alam and Mookerjee comparison (Ni-Fe splitting levels 3 to 4, Cr spikes).
3. Speed: fewer E_F-search contours (one full contour, then a short arc between the old and new E_F; Cauchy gives the same
   moments), a better first E_F guess (574c467 only helps restarts), O_h instead of D4h (up to 3x fewer nodes), parallel or cached operator build.
4. HEA at levels 4 or 5 for the weak moments once the speed work is in. Then more SQS cells for error bars.
5. Part B (short-range order) as in `HANDOFF_hea_sro.md`.
6. Total energy ranking of competing magnetic states (needs a converged SCF and a finer energy check).

## Practical notes
- The desktop rebooted once between 2026-10-21 and 2026-10-23 and the remote-control session `asr` (tmux) did not come back.
  There is no tmux server there now. Restart it if the user wants phone access (memory `desktop-rc-tmux`).
- The SRM VPN drops every few hours: `timeout 60 nmcli con up id "SRM"` (the passwd-file warning is harmless). Wrap ssh in `timeout`.
- Launch pattern: `setsid nohup ... < /dev/null &` inside the ssh command; write the PID from the same shell. A first-step
  mistake this time: `echo $! > pid` after an `&&` chain recorded the bash PID, not the program PID. Check with `pgrep -af app/asr`.
- Kill by PID only. Never `pkill -f` a pattern that also appears in your own command line.
- A long ssh command that runs the build in the foreground can hit the 600 s tool limit. Build first, launch second.
- After a module signature change: `rm -rf build/ifx_*`. All calculations run on the desktop (memory `run-on-desktop`).
- Writing style for text the user reads or for docs: simple sentences, no dashes (memory `writing-style-simple-sentences`).
