# HANDOFF: start here

Updated 2026-10-30.  This is the only current handoff.  Old handoffs are in `docs/archive/`.
Read order for a new thread: this file, `docs/STATUS.md` (dashboard), `docs/ROADMAP.md` (order of work).  Details and numbers are in
`docs/FINDINGS_LOG.md`.

## Code state
| where | what |
|---|---|
| `main` | local: this docs cleanup on top of `3ca0c9a`, not pushed; GitHub: check with `git log origin/main -1` (last known push `574c467`) |
| main working tree, uncommitted | SRO work: `src/asr_sro.f90`, `test/test_sro.f90`, `example/fe66ni14cr20_fcc_sro.toml`, changes in `asr_augcls`, `asr_class_scf`, `asr_driver`, `test_augcls`, `test_run_input`; the user's own edits in `README.md`, `.gitignore`, `src/asr_kinds.f90`, `doc/`, `ford.md`, `suggestion.md`, `docs/JIJ_PLAN.md`, `docs/PUBLICATION_PLAN.md` |
| branch `l4stall` (worktree `../asr_l4stall`), not pushed | `93aa84e` `[asr] ef_scan` probe (N(E) at a fixed potential, split by contour part) and the `bp` fallback counter `ed_nfall`; `03ac4e5` `ef_scan_dump` (chain coefficients for offline tail tests) |
| stash in `../asr_l4stall` | rejected `bp` repair "continued from the real chain" (fails the fine continuity test; FINDINGS_LOG 2026-10-30) |
| default tail in the code | still `bp` (`asr_driver`, `asr_class_scf`).  Production runs must set `ed_terminator = "hull"` until the default changes |

## Desktop state (`ssh desktop`, idle on 2026-10-30)
| clone | commit | use |
|---|---|---|
| `~/asr_remote4` | `5c8075e` | main run tree, `env.sh`, all `trial/` directories |
| `~/asr_q` | `574c467` | binary `bin/asr_574c467` used by the HEA l3 `bp` queue |
| `~/asr_l4` | `03ac4e5` (`l4stall`) | binary of all 2026-10-29/30 runs (`build/ifx_A08D5D447AF67060/app/asr`) |
| `~/asr_spindos` | `03ac4e5` (`l4stall`) | not created by the l4-stall thread; origin unknown, leave it |
| `~/asr_pd`, `~/asr_base`, `~/asr_v1`, `~/asr_v2` | older | disposable diagnostic clones |

Recent runs (`~/asr_remote4/trial/`):
- `a1/l4_hull`: Fe66Ni14Cr20 levels 4 NM, `hull`, converged in 6 iterations (DOS in `fe66ni14cr20_dos.dat`).
- `hea/*_l3_hull`: CoCrFeNi and CoCrFeMnNi levels 3, NM and spin, `hull`, all converged; `hea/*_l3` are the `bp` runs.
- `a1/scan/{bp,hull,avg2,bpdump}`: N(E) scans; `bpdump/chains.dat` (25 MB) with `tailtest.py` (offline replica of `cf_species`).

## Next actions
1. Make `hull` the default `ed_terminator`; keep `bp` as an option with a warning in the manual; add a unit test that N(E) is smooth on
   a fixed chain family (data: Cr e_g chain at node 36 in `bpdump/chains.dat`); run the unit tests; merge `l4stall`.
2. Review and commit the uncommitted SRO work (with the user).
3. Alam and Mookerjee Fig. 2 comparison from `a1/l4_hull` (Ni-Fe splitting from levels 3 to 4, Cr spikes).
4. Speed work (ROADMAP phase 5).

## Working rules
- All calculations on the desktop, unit tests may run on the laptop.  The SRM VPN drops every few hours:
  `timeout 60 nmcli con up id "SRM"`.  Wrap every ssh in `timeout`.
- Code moves by git.  The desktop pulls the https remote `github` (`git pull --ff-only github main`).  Unpushed branches go over as a
  bundle: `git bundle create x.bundle 574c467..<branch>`, `scp`, `git fetch /tmp/x.bundle <branch>:<branch>`.  Only the user pushes.
- Build on the desktop with `source ~/asr_remote4/env.sh; asr_build` (plain `fpm build` misses the MKL link flags).  After a module
  signature change: `rm -rf build/ifx_*`.
- Launch: `setsid nohup bash run.sh < /dev/null > run.log 2>&1 &` then `echo $! > pid` in the same shell.  `run.sh` sources `env.sh`
  and sets `FORT_BUFFERED=false OMP_NUM_THREADS=<n> MKL_NUM_THREADS=1`.  Check with `pgrep -af app/asr`.
- Kill by PID only.  Never `pkill -f` a pattern that also appears in your own command line.
- A long build in the foreground of one ssh command can hit the 600 s tool limit: build first, launch second.
- Production route only `asr_ed` (memory `ed-is-target-route`); LMTO parameters only from our code.
- Docs and text for the user: simple sentences, no dashes (memory `writing-style-simple-sentences`).
