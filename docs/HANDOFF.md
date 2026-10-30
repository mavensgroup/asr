# HANDOFF: start here

Updated 2026-10-30.  This is the only current handoff.  Old handoffs are in `docs/archive/`.
Read order for a new thread: this file, `docs/STATUS.md` (dashboard), `docs/ROADMAP.md` (order of work).  Details and numbers are in
`docs/FINDINGS_LOG.md`.  A thread on short-range order also reads `docs/HANDOFF_sro.md`.

## Code state
| where | what |
|---|---|
| `main` | local `dd435bb` (merge of `l4stall`, SRO, spin-resolved DOS, docs), not pushed; GitHub `origin/main` is `574c467`.  The laptop has one checkout (`asr_new`) and one remote (`origin`) |
| main, committed 2026-10-30 | `176a727` SRO (`docs/HANDOFF_sro.md`), `f2f9ca7` spin-resolved `dos.dat` columns, `455a069` the user's manual and docs; only `suggestion.md` is left uncommitted |
| `l4stall` (merged into main; branch and worktree deleted 2026-10-30) | `93aa84e` `[asr] ef_scan` probe and `bp` fallback counter `ed_nfall`; `03ac4e5` `ef_scan_dump`; `3a47711` default tail `hull`, Nevanlinna-bound test |
| `stash@{0}` in `asr_new` | rejected `bp` repair "continued from the real chain" (fails the fine continuity test; FINDINGS_LOG 2026-10-30) |
| default tail in the code | `hull` (since `3a47711`); the uncommitted SRO example and `doc/manual` now say `hull` too |

## Desktop state (`ssh desktop`, idle on 2026-10-30)
| clone | commit | use |
|---|---|---|
| `~/asr_remote4` | `3986fc2` (merged main, by bundle, 2026-10-30); pulls `dd435bb` and later once main is on GitHub | the only `asr_new` checkout on the desktop: run tree, `env.sh`, all `trial/` directories; full test run in `trial/fulltest` |
| `~/asr_remote4/trial/archive_clones` | | archive of the clones retired on 2026-10-30 (`asr_q`, `asr_l4`, `asr_spindos`, `asr_pd`, `asr_base`, `asr_v1`, `asr_v2`, `asr_sro`, `asr_remote`, `asr_remote2`, `asr_remote3`): per clone `HEAD.txt`, `uncommitted.patch`, untracked run files; binaries `asr_q/bin/asr_574c467` and `asr_l4/asr_03ac4e5` (the binary of the 2026-10-29/30 runs) |
| `~/asr_remote4/trial/sro_bp` | | the SRO runs of 2026-10-24 (`bp` tail; `docs/HANDOFF_sro.md`) |
| `~/asr2`, `~/asr2_omp` | | the older project `mavens-group/asr2`, not this code |

The hub `~/asr.git` is deleted (its branches are in main).  `~/asr_remote4` has one remote, `github`.

Recent runs (`~/asr_remote4/trial/`):
- `a1/l4_hull`: Fe66Ni14Cr20 levels 4 NM, `hull`, converged in 6 iterations (DOS in `fe66ni14cr20_dos.dat`).
- `hea/*_l3_hull`: CoCrFeNi and CoCrFeMnNi levels 3, NM and spin, `hull`, all converged; `hea/*_l3` are the `bp` runs.
- `a1/scan/{bp,hull,avg2,bpdump}`: N(E) scans; `bpdump/chains.dat` (25 MB) with `tailtest.py` (offline replica of `cf_species`).

## Next actions
1. SRO gates (`docs/HANDOFF_sro.md`).  The SRO code is committed.
2. Alam and Mookerjee Fig. 2 comparison from `a1/l4_hull` (Ni-Fe splitting from levels 3 to 4, Cr spikes).
3. Speed work (ROADMAP phase 5).
4. Push main (the user: `git push origin main`; GitHub still had `574c467` on 2026-10-30).  Then `~/asr_remote4` pulls with `git pull --ff-only github main`.

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
