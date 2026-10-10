# Handoff: D4h point-group reduction of the n-component augmented space (written 2026-10-10)

Start a NEW Claude Code session for this (Opus): `cd ~/Devel/ASR/asr_new && claude --model opus`, then tell it:
"Read docs/HANDOFF_augn_sym.md and do the task described there."  Run `asr_pull` first (code moves only via the git hub, see memory
desktop-git-sync).  The other (Sonnet) session keeps working on the desktop queues, the SCF driver and docs: do NOT edit
`src/asr_alloy_scf.f90`, `src/asr_alloy_driver.f90`, `asr_ed.f90`, `asr_augorb.f90`, `asr_augsym.f90` (binary, validated, in production).

## Task
Add the cubic point-group (D4h, 16 operations, as in `asr_augsym`) reduction to the n-component operator `src/asr_augn.f90`, as a NEW module
`src/asr_augnsym.f90` (do not modify the validated unreduced `asr_augn`).  Target: same Lanczos coefficients and Green function as the unreduced
operator, with ~15-40x fewer nodes on fcc/bcc, for any number of species ns, with a non-flippable and/or frozen start site (central species).

## State (commit 7cff90c)
- `asr_augn`: node = (site R, flipped sites with states 1..ns-1); config entry code = site*ns + state (sorted).  Local basis phi_0 = 1.., rank-one
  M^(s)_lk = x_s phi_l(s) phi_k(s); `apply_rand` uses W_lk = sum_s X_s M^(s)_lk (`mt(l,k,s)`), `flipidx(k,i)` = node with R's own state set to k
  (self for k = own state, 0 if absent/not allowed), `lst(i)` = own state, `frozen(i)` = species of a frozen site, `hop(slot,i)`.  H = E_nu + h - h o h.
  `augn_rand_c` applies complex per-species parameters (used by `asr_edn` for the energy-dependent diagonal).
- `asr_edn` (n-species energy-dependent ASR on `asr_augn`): `edn_setup`, `edn_green`; it needs `op%apply` (real, e = 0), `augn_rand_c`, `op%nn`.
- Tests that pass: `test_augn` (ring brute force ns = 3, frozen site, ns = 2 == asr_augorb), `test_asr_ed_n` (3 species ED vs dense).
- Template to copy the SYMMETRY machinery from: `src/asr_augsym.f90` (+ its test `test/test_augsym_orb.f90`, `test/test_edam_sym.f90`/`test_asr_ed_sym.f90` for the
  frozen-site ED check).  Read its header comment: canonical representatives under the group, amplitudes psi(g node) = chi(g) D(g)^T psi(rep),
  normalised coordinates w = sqrt(n_node) psi, links store (representative index, operation g), `btab(l,l2,slot,g)` = S block times D(g)^T,
  `canon_node(perm, site, cfg, ...)` with `perm(site, iop)` from `d4h_ops(coord, 2*lsuper, perm, mats)` and `orbital_rep(mats, dmat)` (asr_symmetry).

## What changes relative to asr_augsym
1. Species labels are group-invariant scalars: the group permutes SITES only.  A config entry code = site*ns + state maps to perm(site)*ns + state; re-sort
   after mapping (`canon_node` must map and sort codes, ties => orbit size as before).  The start site must be invariant (osz = 1) as in asr_augsym.
2. The flip table has ns entries per node (`flipidx(k,i)`, k = own-state target), each a (representative, g) link; the own-state k = lst gives the node itself
   (g = identity).  `md`/`s` of the binary code are replaced by W_lk lookups exactly as in `asr_augn%apply_rand` (use `mt`).
3. Frozen/non-flippable site and complex-parameter application (`augn_rand_c`) must exist in the reduced module too (same signatures as asr_augn:
   `augn_freeze_site`, `augn_rand`, `augn_rand_c` -> name them augnsym_*), because `asr_edn` will be switched to the reduced operator by the other session.
4. Keep the start-orbital character chi(g) = D(g)[L0,L0] (`augsym_set_start` analogue).  Orbital block structure (9 orbitals) is the same as asr_augsym.

## Gates (write `test/test_augn_sym.f90`; all must pass before commit)
- fcc supercell like `test_augsym_orb` (depth 2-3, rank 1), ns = 3 (conc e.g. 0.3/0.5/0.2) with random but physical-sized parameters (use `build_lmto_alloy`
  blocks, species parameters of three made-up potpar_t like in test_asr_ed_n): reduced vs unreduced chain coefficients to 1e-9 for the s, p_z, t2g, e_g
  start orbitals; also with the start site frozen to each species (reduced `edn_green` vs unreduced `edn_green`, a_k, b_k, G to 1e-9).
- ns = 2: reduced `asr_augnsym` == `asr_augsym` coefficients (binary regression), same node counts.
- Report node counts and the reduction ratio.  Keep the unit tests small (the first symmetry test run took >10 min at depth 4, rank 2: use depth <= 3, rank 1).

## Rules / gotchas
- Heavy runs only on the desktop (VPN: `nmcli con up id "SRM"`), but queues q2/q3 occupy it until ~tomorrow: run only small tests locally, no big jobs on the desktop
  without asking.  Never `pkill -f` a pattern (kills your own shell); `timeout` cannot run shell functions (asr_test etc.): call `fpm test --compiler ifx --flag "$ASR_FFLAGS" --link-flag "$ASR_LDFLAGS" --target X`
  directly or without timeout; `rm -rf build/ifx_*` after module type changes; fpm only builds library modules used by a test/app (a new module compiles when a test uses it).
- Build/test: `source env.sh; asr_build; asr_test --target test_augn_sym`.
- Do not change numerical results of the existing binary route.  Commit with the attribution lines from the session reminder; `asr_sync` pushes to the hub and fast-forwards
  the desktop copy; tell the user when done (what passed, node reduction ratios) and update docs/STATUS.md (section "n-component augmented space") and docs/TODO.md item 3b(i).
- After this: the other session will (ii) generalise `alloy_scf` to ns species using the reduced operator, (iii) a 3-component bcc test, (iv) short-range order.
