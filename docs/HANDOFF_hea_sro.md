# HANDOFF (next thread): multicomponent sublattice disorder and short-range order, toward HEA (2026-10-20)

Two parts, in order: **A** multicomponent disorder on sublattices (uncorrelated), then **B** short-range order (Warren-Cowley).
Production route throughout: the energy-dependent AMscf operator (asr_ed form, central site frozen to the species of the component;
memory `ed-is-target-route`).  Accuracy bar: 3-4 % vs the same-model reference (memory `asr-accuracy-target`).

## State at handoff (all verified; details in docs/STATUS.md, sections "Sublattices" and "One input format")
- Code: main at `0e484b4` (local commits `525ffad`, `0e484b4`; the user pushes to GitHub mavensgroup/asr, then the desktop
  `~/asr_remote4` does `git pull --ff-only github main`).  Desktop scratch clone with these sources: `~/asr_pd`, runs in `trial/q33-q37`.
- **One program, one input format** (`src/asr_driver.f90`): decorated sc/bcc/fcc + `[[sublattice]]` (positions, species, concentrations)
  or an explicit cell + `[[site]]` (k-space only).  `[calculation] method = "asr"` (default, any structure) | `"kspace"` (ordered only;
  primitive-cell reduction, one ASA class per sublattice).  No `kind` key.  New features go into this format (memory `single-input-format`).
- **Site-class operator** `asr_augcls`: every class (sublattice) has its own species list (N >= 1) and concentrations; local basis
  phi_l orthonormal in <f,g> = sum_s x_s f g (`asr_augn::augn_basis`); symmetry reduction by the class-preserving subgroup of D4h.
  Gates (test_augcls 11/11): ring of 4 with 3 classes vs exact average 2.7e-15; decorated bcc reduced vs unreduced 4e-13; binary class vs
  augsym 3.4e-15; **ternary class vs augn 1.8e-12 (nodes 188625 -> 13600)**.  The N-species machinery exists; it has not been run in an SCF.
- `asr_edc` (energy-dependent Green function of one component), `asr_class_scf` (components p = (sublattice, species), weights
  w_p = f_kappa x_s, Madelung of the sublattice averages + SIM of the deviation, exact electron count, Anderson with the restart fix).
- SCF stops at tol: r36 Co2Fe0.4Cr0.6Al converged in 26 iterations (dV 1.9e-6).  Moments Co 0.942 / Fe 2.773 / Cr 1.222 / Al -0.082 vs
  SQS 1.017 / 2.849 / 1.186 / -0.086 (Co low = chain resolution at the minority-gap edge, memory `asr-resolution-physics`).
- Cost: that run is ~17 min per SCF iteration on 14 threads (4 components, 36 steps, nodes 8648 / 4639 / 4516).  HEA will have more
  components; measure before launching long runs.

## Part A: multicomponent sublattice disorder (no SRO)
A1. **Gate: Alam & Mookerjee, JPCM 21, 195503 (2009)** (`../references/SRO-1_alam.pdf`): fcc Fe66Ni14Cr20 (stainless steel) and
    Fe80-xNixCr20 (x = 14, 17); also Fe0.05Ni0.05Cr0.90 and Fe0.90Ni0.05Cr0.05 (their Fig. 3).  **Non-magnetic** (they say magnetism is not
    analysed), DOS only: Fig. 2 (partial and total DOS, E_F marked); lattice constant "at the total-energy minimum" (no number in the
    text: use experimental austenite a ~ 3.59 A and say so).  Their H is energy-independent second order + BP: compare shapes and
    peak positions, not digits.  Input: one fcc sublattice, `centering = "F"`, species ["Fe", "Ni", "Cr"], concentrations.
A2. **Same-model reference**: SQS with one potential per species, k-space (`method = "kspace"`, `lattice = "cell"`, `site_classes = false`).
    `tools/sqs.py` is **binary only**: extend it to N species (match the pair probabilities P_ss'(shell) = x_s x_s' for all species pairs on
    the first shells, plus the decoration option it already has).  Then compare moments/DOS for a magnetic ternary (e.g. fcc FeNiCr with
    spin, or bcc FeCoCr) at the 3-4 % bar.
A3. **HEA targets** after the gate: equiatomic fcc CoCrFeNi and CoCrFeMnNi (spin 2; Cr / Mn antiparallel starts via `spin_sign`; expect
    competing magnetic states), and B2-type AlCoCrFeNi (bcc graph, two sublattices, several species on each: sublattice + multicomponent
    together).
A4. **Performance**: node count and time per iteration vs N (3, 4, 5) at levels 4, rank 1, 24 / 36 steps.  Check where OpenMP acts
    (inside `apply`?) and whether components x orbitals x contour nodes can run in parallel instead.
A5. **Watch**: the Cr d sum rule in the Heusler runs is 0.907 at E_F = 6 Ry (Al d 1.067, others ~1.0).  Not understood; check it on a
    Cr-rich case before trusting Cr-heavy HEA results.

## Part B: short-range order
B1. **Formalism** (Mookerjee & Prasad PRB 48, 17724 (1993); Alam & Mookerjee 2009 section 2.1 and appendix): pairwise correlations between
    a central site R0 and its first shell only.  Given the central occupation j, the shell sites are independent with conditional
    probabilities
        P(s at R_k | c at R0) = X_s^(c):  X_s^(c) = x_s (1 - alpha_cs) for s != c,   X_c^(c) = x_c + sum_{s != c} x_s alpha_cs,
    and the augmented operator of a shell site is  N~_Rk = sum_j P_R0^(j) (x) N_Rk^(j)  (projector on the centre's state times the
    conditional operator).  Homogeneity, sum_c x_c X_s^(c) = x_s, needs symmetric alpha_cs = alpha_sc, and X >= 0 bounds alpha.
    Sites beyond the first shell: uncorrelated (x_s).  Ternary: three parameters alpha_AB, alpha_BC, alpha_AC (their Table 1: alpha = 1
    means that pair avoids NN contact).
B2. **The simplification our route gets for free**: asr_ed freezes the start site to the component's species c, so P_R0^(j) -> delta_jc.
    In the recursion for component c, the first-shell sites of the start site are just **uncorrelated with concentrations X^(c)** instead
    of x.  In code that is a shell class around the start site whose concentrations depend on the component:
    - concentrations enter `asr_augcls` only through the `mt` tables (`build_augcls_operator`, `augn_basis(mdl%conc(...))`, around line
      428); the node set does not depend on them -> add a routine that recomputes `mt` for given class concentrations and call it per
      component before the recursion (like `augcls_set_params` / `augcls_freeze_site` in `ed_moments`);
    - give the first-shell sites of the origin their own class id in `cls(:)` (same species list as their sublattice).  The shell is
      invariant under the site group, so the symmetry reduction stays valid; `class_ops` must see the shell class (class-preserving
      subgroup unchanged for a full shell).
    - Partial disorder: the shell of a site may lie on another sublattice (L2_1: the NN of 4b are 8c).  Define alpha per (sublattice
      pair, species pair), first shell only.
    State the model limit in the docs: only centre-shell correlations (shell sites uncorrelated among themselves), as in the papers.
B3. **Electrostatics with SRO**: the SIM term -beta 2 (q_p - qbar) / R1 assumes a random first shell.  With SRO the shell composition of
    component p is X^(p): decide (with a reference) whether to keep SIM or use the explicit shell charge sum_s X_s^(p) q_s.  Component
    weights w_p stay f_kappa x_s.
B4. **Gates**, in order:
    1. alpha = 0 gives exactly the uncorrelated result (unit test, same numbers);
    2. brute force: the ring / small-cluster test of test_augcls with the conditional product measure for a frozen centre;
    3. Tarafder et al. 2006, CuZn bcc 50:50, NN alpha ~ -0.15 to -0.18 (`../references/SRO-cuzn_tarafder2006.pdf`, only 4 pages, DOS
       Fig. 3).  Their band positions are shifted by ~0.19 Ry (non-self-consistent; docs/ENERGY_DEPENDENT_ASR_AND_CUZN_FINDINGS.md):
       compare the SRO-induced *change*, not absolute peaks;
    4. Alam 2009 Fig. 4 (Fe66Ni14Cr20 with alpha triads 1-1-0, 1-0-1, 0-1-1, 1-0-0, 0-1-0, 0-0-1, 1-1-1): DOS shapes and the ordering of the
       second moments M2 (Fig. 5: M2(110) < M2(101) < M2(011));
    5. same-model reference: supercells with target alpha by reverse Monte Carlo (new tool next to sqs.py), k-space with one potential
       per species.
B5. **Input**: e.g. `[[sro]] shell = 1, pairs = [["Fe", "Ni", -0.1], ...]` (plus sublattice names when there are several); reject alpha
    outside the X >= 0 bounds and asymmetric input.  `docs/ROADMAP.md` M7 already plans an `[sro]` section.

## Rules (from memory; keep)
- Calculations run on the desktop (`ssh desktop`, VPN `nmcli con up id "SRM"` when it drops), not on the laptop; unit tests may run
  locally.  Launch with `setsid nohup ... < /dev/null &` (plain `&` inside ssh hangs).  Kill by PID, never `pkill -f`.
- Never swap in the energy-independent route for production; LMTO parameters are generated by our code (no copying from AMscf /
  Questaal).
- After a module signature change: `rm -rf build/ifx_*` (fpm kept a stale test binary on 2026-10-19 and it segfaulted).
- Pushing to the public GitHub needs the user (or the user pushes).  Never put secrets in tracked files.

## Open items carried over (not part of this handoff)
- k-space (`asr_asa`) and binary `asr_alloy_scf` still apply the Anderson restart safeguards at all dV (only class_scf switches them off
  below 1e-3; the ring fix applies everywhere).
- `asr_alloy_scf` P(z) - S route: unstable antisymmetric charge mode in its A = B test SCF (STATUS).
- FeRh Fe -4.8 % vs SQS; README and `doc/manual` still describe the removed `kind` values (the user's uncommitted edits).
