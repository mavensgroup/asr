# HANDOFF (next thread): Part A1 run, SCF mixing cycle, parallelism plan (2026-10-21)

Continues `docs/HANDOFF_hea_sro.md` (Parts A and B, rules).  Read that first; this file records what the 2026-10-20/21 thread did
and what is open.  Production route unchanged: energy-dependent asr_ed class SCF (memory `ed-is-target-route`).

## State at handoff
- Code: local main `e7d6c49` (one commit on top of `736bc75`), **not pushed**; the user pushes, then the desktop
  `~/asr_remote4` does `git pull --ff-only github main` (`asr_sync` in env.sh).  The desktop binary is still `736bc75`.
- `e7d6c49`: DOS file fix + new example + test entry.
  - `src/asr_driver.f90` (writer, ~l.478): `res%edos` is per spin channel; for `spin = 1` the file now carries 2/nspin times it,
    as the header ("summed over spin") always claimed.  The SCF was never affected (it uses `fs = 2/nsp`).  DOS files written by
    the old binary at spin = 1 (all `trial/a1/l3`, `l4` output) are **per spin: multiply by 2**.
  - `example/fe66ni14cr20_fcc.toml`: fcc, one sublattice Fe/Ni/Cr 0.66/0.14/0.20, a = 3.59 A (experimental austenite; the paper
    gives no number), non-magnetic, levels 4, 24 steps, 32 contour points.  Read-checked in `test/test_run_input.f90` (26/26).
- Unit tests: all suites that ran passed; the full `fpm test` exceeds 20 min (`test_asr_ed_sym` alone is long), so the suites after
  it were not run in this thread.  `probe_*` programs abort with "end-of-file during read" without arguments (expected, not tests).
- User's uncommitted edits (README, .gitignore, `asr_kinds.f90` comment, `doc/`, `ford.md`, `suggestion.md`) were left alone.

## Desktop runs (`~/asr_remote4/trial/a1/`, binary `build/ifx_A08D5D447AF67060/app/asr` at 736bc75)
| dir | setup | status |
|---|---|---|
| `l3` | levels 3, spin 1 | **converged**, 27 it, ~10 min/it on 14 threads, nodes 15827 |
| `l3_spin` | levels 3, spin 2, zeta0 0.1, Cr `spin_sign = -1` | **converged**, 42 it, ~30 min/it (14 threads, other job running) |
| `l4` | levels 4, spin 1, `start_potential = ../l3/...` | **running, PID 2617298, not converging** (see SCF cycle below); nodes 64567 |

Results (qnet = electrons - Z, E in Ry):
- `l3` non-magnetic: E_F -0.0703; qnet Fe +0.026, Ni +0.146, Cr -0.189.  d sum rules (M0/5 at 6 Ry) Fe 1.04, Ni 0.93, Cr 1.02.
  d occupations (per-spin file x 2) Fe 6.2, Ni 8.1, Cr 4.4; N(E_F) 21.7 states/Ry/atom.
- `l3_spin`: ferrimagnet, m Fe +1.370, Ni +0.256, Cr -0.931 muB (cell average 0.754 muB/atom); E_F -0.0661; charges as non-magnetic
  within 0.01 e.  **No total energy in the class SCF**, so NM vs ferri cannot be ranked yet.
- `l4` (cycling): qnet Cr between -0.185 and -0.227, E_F -0.0708 to -0.0711; depth changes little vs levels 3.

Comparison with Alam & Mookerjee 2009 Fig. 2 (`../references/SRO-1_alam.pdf` p. 8; both E_F ~ -0.07 Ry absolute):
- agrees: order Ni < Fe < Cr, Cr split off above E_F, charge Cr -> Ni (their "Cr pushed up, Fe/Ni down"), total DOS main peak with
  a low-energy shoulder (paper -0.21 / -0.33, ours -0.19 / -0.31).
- differs: Ni-d peak -0.29 (paper) vs -0.21 (ours), so Ni overlaps Fe more in ours; Fe-d flat up to E_F, so N(E_F)/peak ~0.8 vs ~0.5.
  Their H is energy-independent 2nd order + BP, and asr2 found published ASR results of that group disagreeing by up to 0.19 Ry, so
  50-80 mRy is not a verdict.  Check whether the Ni-Fe splitting moves with depth (l3 -> l4) before judging.
- their DOS axis is probably per spin (peak heights match our per-spin values); the paper does not say.
- small spikes in Cr-d at E_F + 0.18-0.20 Ry (terminator?); check at levels 4.

## Open problem 1: SCF limit cycle at levels 4 (blocks l4)
Log (`l4/run.log`): dV falls to ~3e-4 .. 2e-3, the safeguard in `asr_mixing::anderson_step` (active only for dV > 1e-3, set in
`asr_class_scf.f90:325`) restarts the history (residual > 2x previous and/or step against the residual), the next step is the plain
linear one `xin + beta r` with beta = 0.2, and dV jumps ~8x (1.5e-3 -> 1.3e-2).  15 restarts in 35 iterations.  Reading: linear mixing
at 0.2 is unstable for a charge-sloshing mode (Cr charge +-0.02 e, through the SIM term); Anderson with intact history controls it.
l3 escaped because it got below 1e-3 early (safeguard off) and never restarted again.
- `save_potential` is written **only at the end** (`asr_class_scf.f90:355`), so l4 has no intermediate potential: a restart can only
  start again from `../l3/fe66ni14cr20_pot.dat`.
- Options (user had not decided at handoff): (a) kill PID 2617298, rerun l4 from the l3 potential with `mixing = 0.05`; (b) code:
  after a restart keep the newest (dx, dr) pair or take a damped step instead of a full linear one, and write the potential every
  iteration (cheap, also makes restarts possible); (c) both.  The same safeguard is used by `asr_asa` and `asr_alloy_scf`
  (see the open item in HANDOFF_hea_sro.md).

## Open problem 2: cost (user wants parallelism; "dont implement now" on 2026-10-21, plan agreed in discussion)
Measured: fcc ternary nodes 15827 (levels 3) / 64567 (levels 4) vs 8648 for the largest Heusler component.  The group used is the
class-preserving subgroup of **D4h (order 16)**, not O_h (48): up to ~3x fewer nodes for cubic sublattices if O_h were used.
Where the time goes (`asr_class_scf.f90` `ed_moments`, l.392-415): OpenMP only inside the matvec (`asr_augcls.f90:158-271`, fork/join
per Lanczos step over 15k-64k nodes, memory bound); the loop spin x component x orbital(4) x contour point (~33) is serial, and each
point is an independent chain (H depends on E): ~800 chains for spin 2, ternary.  The whole contour is redone at every E_F search step
(`l.212-214`, up to 16; 6 in iteration 1).  Plan, in order:
1. task-level OpenMP over the independent chains (per-thread copy of the mutable op state: params, frozen site, start orbital, `wk`;
   graph tables shared; reduce mom/dosef);
2. fewer E_F-search contours (one contour + short correction near E_F, or reuse the DOS(E_F) Newton step late in the SCF);
3. MPI (not coarrays) only for multi-node: same task farm, one MPI_Allreduce of the moment arrays; hybrid MPI x OpenMP.  Ask the
   user whether a cluster is available (asr2 has `hpc_convergence/`).
SRO in our frozen-centre route does not change the node or chain count (only the shell-class concentrations per component); the HEA
cost driver is N components x (N-1) local states per flipped site.

## asr2 clues for Part B (all user's own code, reusable; `../asr2`)
- `src/asr_sro.f90:46-216` `build_sro_shell` (derived diagonal, symmetry, alpha <= 1, collective ordering budget with
  `admissible_pair_range`), `sro_from_pair_counts` (supercell histogram -> alpha, for the RMC reference); tests `test/test_sro.f90`.
- Reference data: `data/reference/sro/fe66ni14cr20_fig5_moments.dat` (Alam Fig. 5 M2/M4 at E = 0), `data/reference/cuzn/tarafder2006_*`,
  digitisers `tools/digitize_{alam2009_fig5,tarafder2006_fig3}.py`.  Docs: `docs/MULTICOMPONENT_SRO.md` (sec. 3 constraints, 11-12
  gates and their failures), `docs/TERNARY_HEA_ASR_BRIEF.md`.
- Conventions: alpha > 0 clustering; triad (1,1,1) -> conditional = identity; triad pair order Fe-Ni, Fe-Cr, Ni-Cr.
- Our route needs only the branch c (centre frozen): no correlated move model, no gauge test.  **Triads give X_s^(c) = 0**: the
  shell class must drop that species per component (augn_basis on a measure with a zero weight is singular).
- asr2 lessons: depth 1-3 band energies flip sign of the SRO preference (not converged); compare SRO-induced changes and the Cu-Zn
  d splitting (0.358 Tarafder / 0.365 Banerjee) rather than absolute peaks; the "alloy of Cu with itself" calibration isolates the
  disorder-route error.  asr2's FeNiCr depth-1 M2 ordering failed and its M2 scale (0.31-0.35) differs from the paper's (~0.20):
  pin the paper's moment definition before B4.4.

## Next steps (suggested)
1. Decide/fix the l4 cycle (problem 1), finish l4, compare Ni-Fe splitting and Cr spikes with l3.
2. Total energy in the class SCF (needed to rank NM / ferri here and the competing HEA states in A3).
3. A2: `tools/sqs.py` to N species; k-space SQS reference for the magnetic ternary (compare l3_spin moments at the 3-4 % bar).
4. A4 / problem 2 when the user says go: task-level OpenMP first, O_h reduction.
5. Then Part B (port `build_sro_shell`, shell class per component).

## Practical notes from this thread
- **All calculations on the desktop** (memory `run-on-desktop`); I wrongly ran timing probes on the laptop once.
- The SRM VPN drops every few hours: `timeout 60 nmcli con up id "SRM"` reconnects (a passwd-file warning is harmless).  Wrap ssh in
  `timeout` so a dropped link does not hang the shell.
- Launch pattern that works: `ssh desktop 'cd DIR; source ~/asr_remote4/env.sh; FORT_BUFFERED=false OMP_NUM_THREADS=14 setsid nohup BIN
  in.toml > run.log 2>&1 < /dev/null & echo $! > pid'` with `;` separators.  Chaining with `&&` before the `&` backgrounds the whole
  chain and the ssh session hangs (pid file then holds the bash PID).  The binary needs `env.sh` sourced (MKL libraries).
- Never `pkill -f` a pattern that also appears in your own command line: it kills the calling shell (happened twice).
- Fortran stdout to a pipe is buffered; set `FORT_BUFFERED=false` to see progress.
