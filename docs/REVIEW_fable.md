# Review of code and physics (Fable 5.2 session, 2026-11-10)

Branch `review-fable`, off `main` e8d193c.  Not merged.  Pushed to `origin` so that the desktop and smriti could build
it.  Task: `docs/HANDOFF_review_fable.md`.
The code work and the unit tests were done on the laptop.  Runs: desktop `~/asr_review` (F4 on 2026-11-11, `trial/f4`;
section 4 jobs in `trial/desk1`), smriti gfortran 14 jobs 398595, 398596 and 398614 (3.6).  The smriti SRO results
(`~/asr_trial/mpi/prod`) and the desktop per-atom cells (`~/asr_remote4/trial/b45/ks_site`) were read for F1.
References checked on 2026-11-11: the SRO measure of 3.1 is the one of Tarafder et al. PRB 74, 144204 (2006), Sec. IV,
and of Alam and Mookerjee JPCM 21, 195503 (2009) (`../references`).

Marks: **proven** means a derivation plus a test or a number from a run here.  **suspected** means an argument without a
decisive test.

## 1. Findings, most severe first

### F1. Full-ASA band energy lost $E_\nu p^\gamma q_2$ per channel (proven, fixed, commit bca3a15)
- Where: `src/asr_asa.f90` (band term of `esum`), `src/asr_class_scf.f90` (band term of `esum`, route `ed3`).
- What: both used $E_\mathrm{band} = \sum_l (q_1 + E_\nu q_0)$.  In the generalised problem $H z = \varepsilon O z$ the $q_k$
  are matrix elements.  The sphere charge of a channel is $q_0 + p^\gamma q_2$.  For every eigenvector
  $\varepsilon z^T O z = z^T H z = \sum_{RL} [q_1 + E_\nu (q_0 + p^\gamma q_2)]$, because all factors next to $h$ in $H$ are
  diagonal.  The same follows in one sphere from $H\phi = E_\nu\phi$, $H\dot\phi = E_\nu\dot\phi + \phi$.
- Evidence: `test/test_band_energy.f90` compares the moment form with the occupied eigenvalue sum of the same bands and
  weights (B2 \ce{FeCu}; tetrahedra, Gaussian, spin).  Before: off by 0.52 to 0.55 Ry per cell.  After: $4 \times 10^{-14}$ Ry.
  The second-order branch was exact before and after ($2 \times 10^{-14}$).
- Effect on the \ce{CuZn} SRO energies (2026-11-11).  The change of an SRO energy is $\Delta \sum_p w_p \sum_l E_\nu p^\gamma q_2$
  between two SRO states.  It was evaluated from the printed $E_\nu$, $p^\gamma$ and $q_2$ of the converged runs.  In a converged
  run the printed table matches the moments (`DERIVATION_ed3.md` section 9 item 1).  Check: for class Zn1 of `a0_s1` the table
  gives qnet $-0.01318$ against the printed $-0.01316$.  `ed3` (smriti `~/asr_trial/mpi/prod`): the correction is $-76.2$ to
  $-76.8$ mRy per site.  Per-atom cells (desktop `~/asr_remote4/trial/b45/ks_site`, 128 sites, 384 table rows each): $-76.72$
  (`a0_s1`), $-76.35$ and $-76.33$ (`a-0.15_s1`, `_s2`), $-77.10$ (`a0.15_s1`), $-76.18$ (`a-0.30_s1`) mRy per atom.  Cell dE is
  taken against `a0_s1`.
  | levels | $\alpha$ | `ed3` dE old (mRy) | `ed3` dE new | cells old | cells new | new ratio `ed3`/cells | dq_Cu ASR | dq_Cu cells |
  |---|---|---|---|---|---|---|---|---|
  | 3 | $-0.148$ | $-1.331$ | $-1.142$ | $-1.44$, $-1.49$ | $-1.066$, $-1.098$ | 1.04 to 1.07 | $+0.0126$ | $+0.012$ |
  | 3 | $+0.148$ | $+1.528$ | $+1.269$ | $+1.48$ | $+1.099$ | 1.15 | $-0.0129$ | $-0.012$ |
  | 3 | $-0.266$ | $-2.643$ | $-2.287$ | $-2.30$ | $-1.766$ | 1.30 | $+0.0225$ | $+0.019$ |
  | 4 | $-0.266$ | $-2.106$ | $-1.726$ | $-2.30$ | $-1.766$ | 0.98 | $+0.0225$ | $+0.019$ |
  The correction moves the cell energies by 23 to 26 %, and `ed3` by 13 to 17 %.  So the old agreement at $\pm 0.148$ was
  partly the same error on both sides.  With the correct form `ed3` overshoots by 4 to 15 % at levels 3.  At $-0.266$ levels 3
  overshoots by 30 % and levels 4 agrees to 2 %.  So the energy at strong order is not converged in levels at levels 3.  The
  charge is converged there (dq_Cu is the same at levels 3 and 4).  `edc` (energy moments, not changed by F1) gives $-1.28$,
  $+1.55$ and $-2.16$ mRy, 17 to 41 % above the corrected cells.
- Check of the table method (done 2026-11-11): the cells `a0_s1` and `a-0.30_s1` were rerun with the fixed binary c12eb4f
  (desktop `~/asr_review/trial/desk1/ks_site`).  They give $-3446.84486201$ and $-3446.84662776$ Ry per atom.  The table values
  are $-3446.84486167$ and $-3446.84662747$.  The difference is $3 \times 10^{-7}$ Ry in both, so dE is $-1.7658$ mRy both
  ways.  The corrected cell column stands.
- Effect: total energies of the k-space full ASA and of `ed3` change.  Charges, $E_F$ and potentials do not.  The default
  route `edc` does not change.  fcc \ce{Cu}, k-space full ASA, VBH, nk 20 tetrahedra, nine lattice constants 3.44 to 3.64 A,
  third-order fit in $V^{-2/3}$ (desktop `~/asr_review/trial/desk1/cu_eos`, 2026-11-11): `main` e8d193c gives $a_0$ 3.523 A and
  $B$ 172 GPa; c12eb4f gives $a_0$ 3.563 A and $B$ 185 GPa.  So F1 moves $a_0$ by +0.040 A (1.1 %) and $B$ by +13 GPa.  The
  energy at fixed $a$ moves by $-44$ to $-49$ mRy.  The fcc \ce{Cu} check in `STATUS.md` ($a_0$ 3.546 A, $B$ 199 GPa) is not
  reproduced by the old binary either, so it used another route or other settings.
- Docs: `docs/DERIVATION_ed3.md` section 9 item 2, `doc/manual/theory/scf_electrostatics.md` ("Total energy").

### F2. Reduced operator: rounding leaves the symmetric subspace and is amplified (proven, fixed, commit f19cb05)
- Where: `src/asr_augcls.f90` (new `augcls_project_c`, stabiliser tables), `src/asr_augblk.f90` (`blk_project`),
  chains in `src/asr_edc.f90` and `src/asr_ed3.f90`.
- What: a reduced vector is valid only if $\psi(r) = \chi(g) D(g)^T \psi(r)$ for every $g$ in the stabiliser of node $r$.
  The reduced operator is symmetric on this subspace $V$ only.  Outside $V$ it is another matrix with a relative asymmetry
  of 0.1 to 0.4.  The Lanczos recurrence amplifies the rounding part outside $V$.
- Evidence (ordered B2 on a bcc supercell of 128 sites, 24 reduced nodes, dense solve on the unreduced operator as
  reference): the part of the `ed3` Lanczos vector outside $V$ is $5 \times 10^{-16}$ at 8 steps, $1 \times 10^{-11}$ at 24,
  $2 \times 10^{-5}$ at 36, 3 % at 48 and up to 17 % later.  The `ed3` blocks then stall $4 \times 10^{-5}$ from the dense
  result at any length.  Full reorthogonalisation does not help ($1 \times 10^{-4}$).  With the projection the chain reaches
  $1 \times 10^{-14}$.  `test/test_augcls_project.f90` (14 checks) holds the projector properties, the dense comparison and
  the block comparison.
- Effect on production: on a random fcc ternary operator (levels 2, 1704 reduced nodes) the projection changes $G$ by less
  than $10^{-10}$ up to 48 steps.  So random alloys at 24 to 48 steps are not affected.  Ordered and partly ordered cells
  with long chains are (gate 2 used 96 and 160 steps on B2).  Production B2 `ed3` at 96 steps (doubt 7): the block
  difference stays at $10^{-4}$ e with the projection, so F2 is not its cause on that case.
- Not fixed: the older reduced operator `asr_augsym` (routes kept as cross-checks) has the same structure and no projector.

### F3. Run-to-run differences at rounding level came from `orbital_rep` (proven, fixed, commit 5545a36)
- Where: `src/asr_symmetry.f90`, `orbital_rep`.
- What: the orbital matrices were fitted on points from an unseeded `random_seed()`.  Entries of a monomial matrix were
  $1 \pm 10^{-15}$ and changed in every run.  They multiply every reduced operator application.
- Evidence: `test_orbital_rep` has two new checks (entries exactly $0, \pm 1$; same bits in a second call).  Both failed
  before and pass now.  A 96-step `ed3` chain differed between two runs on one thread before the fix.  After it all printed
  digits agree between runs and between 1 and 4 threads.
- This resolves doubt 8.  It is not a thread-order effect.  No OpenMP sum is involved.

### F4. The terminated energy-dependent chain is not causal in general (proven by example, not fixed)
- What: the $m$-step chain is the bilinear Galerkin approximant of $(zM - A)^{-1}$ on a Krylov space that depends on $z$
  ($M = T^2$, $A = T H^{(2)} T$).  It is a real rational function of $z$.  Its poles need not be real.  The `hull` tail adds
  complex band edges, so its branch cut can lie above $z$.
- Evidence: scan of an fcc ternary test operator (levels 2, three central species, four start orbitals, 801 energies,
  $\mathrm{Im}\,z$ = 0.02, 0.005, 0.001; 24, 36, 48 steps).  The scan was stopped after 28 of the 36 cases.  With `hull`: 6 of these 28 have negative DOS points; the worst
  is $-2.1$ states/Ry at $E = 0.330$ ($\mathrm{Im}\,z$ = 0.005, 48 steps, p orbital), 1.1 % of the weight.  Without a tail:
  similar isolated points, and $\mathrm{Im}\,G = +6.9$ at 72 steps in one case.
- Consequence: the dips and negative values reported for the \ce{Fe66Ni14Cr20} SRO DOS have a known mechanism.
- Production case (desktop, 2026-11-11, binary `0ee51d6`, so with the projection of F2): fcc \ce{Fe66Ni14Cr20}, triad (1, 1, 0),
  spin, levels 3, 24 steps, `hull`.  DOS at $\eta$ = 0.02 (321 energies) and $\eta$ = 0.005 (801 energies, $-0.6$ to 0.4 Ry):
  no negative value in any of the 30 columns (components, $l$, spin).  So the report is not reproduced in this setting.  Caveat:
  the SCF from free atoms did not converge in 99 iterations (dV 0.05 to 0.17 at the end, Cr charge $-4.8$ e after iteration 1).
  A second SCF from its potential reached dV $1$ to $5 \times 10^{-4}$ and stalled.  The DOS is at that potential.  Control
  (desktop `trial/f4/ctl_main`, `ctl_review`): `main` e8d193c gives the same iteration 1 as this branch to all printed digits
  (E_F search, E_F $-0.467557$, qnet 1.3098 / 0.6566 / $-4.7820$, dV 78.01, E $-2450.32088120$).  So the unstable SCF comes
  from this input with a free-atom start and spin, not from the branch.  A non-magnetic start potential does not help.
  `read_potentials` copies a one-spin file into both spin channels, and the seed `zeta0` acts only on the free atoms.  So
  such a start is spin symmetric, and the SCF keeps it non-magnetic.  The next try is a smaller mixing (0.1, restart 0.03)
  from `trial/f4/scf/pot.dat`.  It needs about 150 iterations (about 8 h), so it waits for smriti (input
  `trial/desk1/fnc_spin.toml.for_smriti`).
- A causal alternative exists in exact form: the same $G_c$ is the resolvent of the pencil $(A, M)$ with real, fixed
  operators.  A Lanczos chain in the $M$ inner product has real coefficients and is Herglotz at every length.  It is the
  energy-independent route on $H^{(2)}$ with a frozen centre.  The ground rules exclude a swap, so nothing was changed.
- Docs: `doc/manual/theory/energy_dependent.md` ("Causality of a truncated chain").

### F5. Comment and doc errors corrected
- `energy_dependent.md` said the identity $E - H(E) = T (E - H^{(2)}) T$ holds "at infinite depth".  It is exact for every
  configuration.  Only the agreement of the two chains needs convergence.
- `asr_class_scf.f90` says the first-order $E_F$ step has an error $O(\mathrm{step}^2)$ in the moments.  It is first order:
  the step uses the DOS at $E_F + i\,$`ed_eta` (Lorentzian, 0.01 Ry), while the contour derivative is the thermally
  broadened DOS ($kT$ = 0.001 Ry).  The electron count stays exact.  The comment is corrected (see doubt 6).

## 2. The ranked doubts

| # | doubt | status | reason |
|---|---|---|---|
| 1 | SRO model and band positions | narrowed | see 3.1 |
| 2 | `ed3` at strong order | narrowed | levels 3 to 4 at $-0.266$: dq_Cu unchanged (+0.0225 against +0.019 in the cell), so the 18 % charge overshoot is not truncation; with F1 on both sides dE is $-2.29$ (levels 3) and $-1.73$ (levels 4) against $-1.77$ in the cell, so the energy overshoot at levels 3 is truncation in levels (table in F1) |
| 3 | energy-dependent form | resolved | identity proven by algebra (section 3.2); code checked line by line; existing ring tests agree to $10^{-15}$ |
| 4 | `ed3` derivation | resolved with one correction | Results 1 to 5 re-derived; band energy corrected (F1); gates judged in 3.3 |
| 5 | total energy | narrowed | F1 fixed; functional form and double counting checked (3.4); 0.1 mRy claim not yet tested by a run |
| 6 | `ef_linear_max` | narrowed | same fixed point proven; the step term is excluded, so the 1.5e-3 $\mu_B$ is the distance to the fixed point at tolerance 1e-5 (3.5) |
| 7 | block kernels | narrowed, not resolved | not a kernel bug: the single-energy chain differs from itself by the same amount under a $10^{-15}$ shift of $z$; on the 128-site test cell block equals single to $10^{-14}$ at 144 steps after F2.  Production B2 `ed3` (supercell 16, 96 steps, 64 points, one iteration; desktop `trial/desk1/f2_b2`, 2026-11-11): with c12eb4f eb0 and eb8 still differ (qnet 0.12916 / 0.12928, E_F $-0.162172$ / $-0.162199$).  A repeat of eb0 is bit identical, so this is a fixed rounding sensitivity of the long ordered chains, about $10^{-4}$ e.  F2 does not remove it there.  The old binaries were not even reproducible: `main` and `5473149` give qnet 0.12935 and 0.12904 for the same eb0 input (F3) |
| 8 | thread-order dependence | resolved | F3 |
| 9 | OpenBLAS ring test | resolved | not OpenBLAS and not the library: the test's own configuration count was wrong under gfortran 14; fixed in the test, passes there (3.6) |
| 10 | MPI shutdown crash | not reached | `asr_mpi.F90` read; no defect found by reading; no test possible here |
| 11 | terminators and causality | narrowed | F4: violations proven on a test operator; the FeNiCr SRO spin DOS at levels 3 has no negative value |
| 12 | electrostatics (B3) | resolved by derivation | 3.7 |
| 13 to 15 | lower priority | not reached | F2 does not explain 13 on the evidence here (no effect in a random operator up to 48 steps) |

## 3. Details

### 3.1 SRO (doubts 1 and 2)
- The conditional-shell construction in `asr_augcls` and `asr_class_scf` is the Mookerjee-Prasad measure exactly: centre
  frozen, first shell with $X^{(c)}$, all other sites with $x$.  The brute-force ring test proves the operator.
- What the model gets exactly.  Any quantity that is linear in the species of single sites, given the centre, needs only
  the pair probabilities.  Two cases: the first energy moment of the local DOS of $H^{(2)}$ (it involves sites within the
  $S$ range only), and the mean potential of the charge-excess model ($V_0 = \sum_j W_{0j} b_j$).  For these the model is
  exact if the pair correlations of all shells are right.  It has $\alpha_n = 0$ for $n \ge 2$.
- So two error sources are separate: missing pair correlations of further shells, and missing correlations among the
  shell sites.  The first is cheap to add in the same operator: one more class per shell with its own conditional
  concentrations, and $F = (G_{00} + \sum_n G_n \hat A_n)^{-1} (W_{00} + \sum_n W_n \hat A_n)$ in `asr_screening`.
  The $-0.266$ reference cell has $\alpha_2 = 0.07$, $\alpha_3 = 0.05$.  A bare estimate of the shell-2 term on the Cu
  potential is $z_2 (2/R_2) \alpha_2 q_\mathrm{Cu} = 0.234\,\alpha_2$ Ry before screening.  With the fitted weight 0.11 this
  is 1.8 mRy on Cu and the opposite on Zn at $\alpha_2 = 0.07$.  It reduces the band separation, as observed.  This is
  **suspected** as a main part of the overshoot at $-0.266$.  Test: run the ASR with shell 2 and 3 classes at the cell's
  own $\alpha_2$, $\alpha_3$.
- A correlated next order is derivable.  On a tree the pair-correlated measure is a Markov chain,
  $P = \prod_i X(s_i | s_{\pi(i)})$.  Use local bases that are orthonormal for each parent species.  Then a site-diagonal
  parameter at a site of generation $n$ has matrix elements that also change the local states of its ancestors:
  generation 1 gives the present single-site weights, generation 2 gives two-site weights
  $\sum_{s_j} X(s_j|c)\,\phi_{l_j}\phi_{l_j'}(s_j)\,W^{(s_j)}_{l_k l_k'}$.  The node set is unchanged.  The operator needs
  pair-flip links.  On a real lattice a second-shell site has several parents (four in bcc), and in fcc the shell sites
  are neighbours of each other.  A closure is then needed (for example the Kirkwood product), and the weights depend on
  all parents.  This is a new operator feature.  It was not implemented.
- The factor 1.45 at $\pm 0.148$ is not explained here.  The log's own analysis says $dC = dV_\mathrm{ext} - U\,dq$ is a
  small difference of large terms.  With `ed3` the charge change agrees to 5 %, so $dV_\mathrm{ext}$ must be about 15 %
  larger in the ASR than in the cells.  This needs the cell data on the desktop.

### 3.2 The energy-dependent form (doubt 3)
$T h T = (D_c/D_R) c + d_c S d_c = g + S'$.  $T (E - E_\nu) T = (D_c/D_R)(E - E_\nu)$.
$T h o h T = h' (T^{-1} o T^{-1}) h' = h' q h'$ with $q = o D_R/D_c$.  So $E - H(E) = T (E - H^{(2)}) T$ with
$H(E) = [E - D_c (E - E_\nu)/D_R] + h' - h' q h'$.  `edc_setup` sets exactly $g$, $|d_c|$, $q$ (`src/asr_edc.f90:48-50`),
and `edc_green` adds $e_s(z)$ (`:72`).  $|d_c|$ instead of $d_c$ is a diagonal sign transformation.  $T = d_c/d_R$ with
signed $d$ makes $T d$ constant for any sign convention of $\sqrt\Delta$.  $T_{R_0} = 1$ gives
$G_c = \langle R_0 | (E - H^{(2)})^{-1} | R_0 \rangle$ on every configuration with $c$ on $R_0$.

### 3.3 The gate redefinitions of `ed3` (doubt 4)
- Gate 2 (section 10a): justified.  The kept $S$ range is a convergence parameter.  `rcl` 2.0 moved `edc` onto its own
  k-space reference as well (0.1221 against 0.1219).  The limits were not changed.
- Gate 3 (section 10b): justified with a caveat.  The reweighted model at fixed $E_F$ is not neutral, and the ASR is.  The
  correction was tested on a prediction for `edc` made before the run (1.6 %).  The caveat is the open $E_F$ shift
  (+0.00064 Ry against about +0.0015 Ry), and the pass holds only at supercell 16 with 64 contour points.

### 3.4 Total energy (doubt 5)
- Kinetic term: $\sum\varepsilon - \int n \tilde v_\mathrm{in}$ with the potential that the bands used, Madelung shift
  included on both sides.  Correct after F1.
- Hartree: `poisson_hartree` with $\tfrac12 \int n V_H$, nuclear term cancelled as documented.  Madelung:
  $\tfrac12 q_p V^\mathrm{M}_p$ with the output charges.  With SRO the cross term vanishes because
  $\sum_{p \in \kappa} x_p (F dq)_p = 0$ (marginal condition of $\hat A$).
- Smearing entropy: about $\tfrac{\pi^2}{6} N(E_F) (kT)^2$, 0.03 mRy at $kT$ = 1 mRy and 20 states/Ry.  Differences are
  smaller.
- Open: no run here tests 0.1 mRy.  The functional is variational, so the error of a difference is second order in the
  density error.

### 3.5 `ef_linear_max` (doubt 6)
The $E_F$ carried to the next iteration is $E_F + (N_\mathrm{elec} - N)/D_\eta$.  At a joint fixed point of potential and
$E_F$ the step is zero, so $N(E_F) = N_\mathrm{elec}$ and the moments are those of a full contour.  The default path has
the same fixed point.  The iteration contracts if $|1 - N'/D_\eta| < 1$.  The moments at a non-zero last step carry an
error of first order, $|\mathrm{step}| \times$ (species and spin DOS mismatch between $\eta$ and $kT$ broadening).  The
SCF stops on $dV$ only.  The logs (desktop `~/asr_ed3/trial/speed/s3_e33fe0b`, 2026-11-11) settle which term it is.  The run
with `ef_linear_max` 0.005 stopped at iteration 42 with dV $9.8 \times 10^{-6}$ and a last step of $-3.2 \times 10^{-7}$ Ry.
That step times a DOS of order 10 states/Ry per spin gives less than $10^{-5}$ $\mu_B$, so the step term is far too small.
Its previous iteration had dV $5.3 \times 10^{-5}$.  The default run stopped at iteration 57 with dV $4.5 \times 10^{-6}$, after
$1.5 \times 10^{-5}$ at iteration 56.  So both runs stop at different points of a slow approach to one fixed point.  The 1.5e-3 $\mu_B$ (Fe 1.34077 against 1.33929,
Cr $-0.91661$ against $-0.91534$) is that distance at tolerance 1e-5, a soft magnetic mode.  Charges agree to 3e-5 e.  Not
tested: a rerun at tolerance 1e-6 (long, a smriti job).

### 3.6 OpenBLAS ring test (doubt 9)
- The reference is well conditioned: $\| K^{-1} \| \le \max(D_s/D_c) / \mathrm{Im}\,z \approx 6$, condition number below 100.
- The recursion side calls no BLAS or LAPACK routine whose result it uses.
- Here: gfortran 16.2 with OpenBLAS 0.3.34 and with reference LAPACK (FlexiBLAS) pass, 40 of 40 random seeds, largest error
  $8 \times 10^{-15}$.
- smriti `ringtest/RESULT.txt`: the gfortran 14 build fails with err 18.7 both with OpenBLAS and with MKL.  So the library is
  not the cause.  The node counts of all six ring cases equal those of the passing builds.
- Cause (resolved 2026-11-11, smriti jobs 398595, 398596, 398614).  Per-case output showed that the exact reference was the
  recursion value times a constant: 0.048, 1.31 and 6.6 for start classes 1, 2, 3, the same for s and d orbitals and at every
  energy.  So the configuration weights of the reference did not sum to 1, and the library was right.  The configuration count
  came from an implied-do inside `product()`; gfortran 14 evaluated it wrongly.  Commit 4c70b9f computes it in a plain loop and
  adds a check that the weights sum to 1.  gfortran 14, optimised and with `-fcheck=all -finit-real=snan`: all ring checks pass
  (weights $2 \times 10^{-16}$, ring $3 \times 10^{-15}$).  The `associate` change of 919f73f was not the cause; it is harmless.
- The same gfortran 14 suite found a real defect: `build_augcls_operator` assigned the species parameters also when they were
  not allocated (undefined behaviour; ifx tolerates it, gfortran 14 crashed in `test_augcls_blk`).  Fixed in d54f018;
  `test_augcls_blk` passes 18/18 with both gfortran 14 builds.  Full gfortran 14 suite (398595, before these two fixes): the
  other 48 programs pass, including `test_band_energy`, `test_augcls_project` and `test_orbital_rep`.
- The ring test now uses a fixed seed (919f73f).
- Branch `lapack-interfaces`: it builds with gfortran 16 here, and that build was used for this test.  No defect seen.

### 3.7 Electrostatics (doubt 12)
$q = G b$, $V = W b$ with $G = a/(a + M)$, $W = a(1 - G)$.  Single-site elimination gives
$dq = (G_{00} + G_1 \hat A) db$ and $V = (W_{00} + W_1 \hat A) db$.  Both matrices are functions of $\hat A$ and commute, so
the code's $F$ equals $V dq^{-1}$.  $W_1 = -a G_1$ because $\langle \cos k R_1 \rangle = 0$.  Without SRO
$F = -2\beta/R_1$, with the sign of `scf_electrostatics.md` (electron excess lowers its own potential).  The single-site
energy is the exact average of the model: the fluctuation terms of $\tfrac{a}{2}(q - b)^2$ and of $\tfrac12 q V$ cancel
($+\mathrm{var}(V)/2a$ and $-\mathrm{var}(V)/2a$).

## 4. What needs the desktop or smriti
Rule (user, 2026-11-11): the desktop takes jobs of 2 to 3 h.  Long jobs and jobs that need MPI wait for smriti.
1. Effect of F1.  Done: the correction of both sides of the SRO comparison from the tables (F1 table), and fcc \ce{Cu} $E(a)$
   (F1), and the reruns of `a0_s1` and `a-0.30_s1` with the fixed code (F1, 6.1).
2. Done: F2 on B2 `ed3` at 96 steps, `energy_block` 0 and 8, new and `main` binary, and a repeat (doubt 7).  The gap of
   $10^{-4}$ e remains.  It is reproducible with c12eb4f.
3. Done: the F4 control shows that `main` has the same unstable start (F4).  Done on the desktop (2026-11-12,
   `~/asr_review/trial/desk2/fnc_sro_spin`, binary 36c10b1): the spin SRO SCF with mixing 0.1 from `trial/f4/scf/pot.dat`
   converged in 98 iterations (dV $6.6 \times 10^{-6}$; no clamp event; weight below the contour $1.4 \times 10^{-9}$ e).  E_F
   $-0.072941$ Ry, qnet Fe / Ni / Cr 0.0002 / 0.0366 / $-0.0264$, moments 0.351 / 0.093 / 0.014 $\mu_B$.  It is the same
   solution that the stalled second SCF of F4 approached (same charges).  The DOS at $\eta$ = 0.005 (801 energies) has no
   negative entry in any of its 27 columns; the smallest value is 0.032 states/Ry.  So the reported negative DOS is not
   reproduced at the converged potential either.  With Fe surrounded by Fe only (triad 1, 1, 0) the Fe moment falls from
   1.34 (random alloy) to 0.35 $\mu_B$, close to fcc Fe at this lattice constant; Cr is almost non-magnetic.
4. Read on 2026-11-11: `ringtest/RESULT.txt` (3.6) and the levels-4 results (F1 table, doubt 2).  The gfortran 14 suite of this
   branch ran (3.6).
5. Done: the two FeNiCr logs for doubt 6 (3.5).
6. Submitted on smriti (`~/asr_trial/review/m015_L4_r20`, `p015_L4_r20`, two nodes each, binary 36c10b1): `ed3` levels 4 at
   $\pm 0.148$, because levels 3 is not converged for the SRO energy (F1 table).
7. Done on the desktop (`~/asr_review/trial/desk2/crfeconi_b12`): CoCrFeNi levels 3, 24 steps, `sim_beta` 1.2 (6.6).
8. smriti state on 2026-11-12: all 110 `cpu` nodes are down or drained.  The four jobs (398637 to 398640) wait.  Items 3 and 7
   therefore also run on the desktop (`~/asr_review/trial/desk2/fnc_sro_spin`, `crfeconi_b12`, binary 36c10b1, 14 threads
   each, started 2026-11-12).  The levels-4 `ed3` runs of item 6 do not fit the desktop: they took 9 to 11 h on 4 nodes of
   48 cores, which is about 3 days each on 28 cores.

## 5. Tests
New: `test_band_energy` (4 checks), `test_augcls_project` (14 checks), two checks in `test_orbital_rep`.
Full suite on this branch at f19cb05 (ifx, laptop, 2026-11-10): 50 test programs, 435 checks, 0 failed.
Full suite at 36c10b1, with `test_scf_guards` (ifx, desktop `~/asr_review/trial/desk1/suite2/suite.log`, 2026-11-12): 51 test
programs, 445 checks, 0 failed (same fpm exit code from the `probe_*` programs).
Full suite at c12eb4f, after d54f018 and the ring-test changes (ifx, desktop `~/asr_review/trial/desk1/suite/suite.log`,
2026-11-11, 47 min on 14 threads): 50 test programs, 436 checks, 0 failed.  fpm exits with code 24 because the six `probe_*`
programs stop without their arguments or local reference files, as before (`STATUS.md`).  They are not tests.

## 6. Scrutiny: does the implementation match the theory, and what is ad hoc (2026-11-11)

Read in full, line by line against the theory documents: `asr_class_scf`, `asr_edc`, `asr_ed3`, `asr_augcls` (operator, builder,
projector), `asr_sro`, `asr_screening`, `asr_endep::cf_species` (terminator), `asr_mixing`, `clamp_enu`.  Not read again here:
the radial solver, the atom SCF, the xc wrapper, the structure constants and the contour nodes.  Those are covered by their own
tests and by the k-space comparisons.

### 6.1 Parts where the code is the theory (no free choice)
| part | theory | code | check |
|---|---|---|---|
| `edc` operator | $E - H(E) = T (E - H^{(2)}) T$ | `edc_setup` sets $g$, $\lvert d_c \rvert$, $q$; `edc_green` adds $e_s(z)$; `augcls_apply` is $e + h - h o h$ | identity to $10^{-15}$ on random matrices (10d of `DERIVATION_ed3.md`); ring tests |
| `ed3` operator | $K' = \tilde\varepsilon - h' + u h' + h' u - h' \tilde Q h'$ | `ed3_kmul`, `ed3_coefs` | same; `test_asr_ed3`, `test_asr_ed3_cls` |
| `ed3` moments | $q_0, q_1, q_2$ from chains $\lvert a\rangle, \lvert v\rangle$ | `ed3_blocks`, `ed_moments` | same |
| $o^\alpha$, $p^\alpha$ | $1/o^\alpha = C - E_\nu - \Delta/(\gamma - \alpha)$, $p^\alpha = p^\gamma + (o^\alpha)^2$ | `to_alpha`, `o_alpha` | re-derived here; equals the code |
| density | `edc`: $M_0 \phi^2 + 2 M_1 \phi\dot\phi + M_2 (\dot\phi^2 + \phi\ddot\phi)$, charge $M_0$; `ed3`: $q_0 \phi^2 + 2 q_1 \phi\dot\phi + q_2 \dot\phi^2$, charge $q_0 + p^\gamma q_2$ | `class_scf` lines 381 to 385 with `phi2` | sphere integrals re-derived |
| band energy | $\sum_l [M_1 + E_\nu \times \text{channel charge}]$ | `class_scf` line 380 | `test_band_energy`; cell rerun below |
| total energy | $\sum\varepsilon - \int n \tilde v_\mathrm{in} + \tfrac12 \int n V_H + \int n \varepsilon_{xc} + \tfrac12 q V^\mathrm{M}$ | `class_scf` lines 370 to 434 | 3.4 |
| SRO probabilities | $X^{(c)}_s = y_s (1 - \alpha_{cs})$, derived diagonal, marginals | `asr_sro` | equals Tarafder 2006 Sec. IV and Alam 2009 |
| SRO measure | centre frozen, first shell conditional, rest random | shell classes in `class_scf`, `augcls_set_class_conc` | brute-force ring test |
| SRO electrostatics | $F = (G_{00} + G_1 \hat A)^{-1} (W_{00} + W_1 \hat A)$ | `asr_screening` | 3.7 |

Direct check of F1 on a production cell (desktop `trial/desk1/ks_site/a0_s1`, binary c12eb4f): the rerun gives
$-3446.84486201$ Ry per atom.  The table correction predicted $-3446.84486167$.  They agree to $3 \times 10^{-7}$ Ry.  E_F is
unchanged ($-0.03542808$ Ry).

### 6.2 Safeguards that act on the path only (not on the fixed point)
| item | where | why it does not change the result |
|---|---|---|
| step limits on the potential (0.15 Ry) and on $E_\nu$ (0.08 Ry) | `class_scf` 470, 475 | zero step at the fixed point |
| $E_\nu$ update $E_\nu + 0.7\, M_1/M_0$ | 444 | fixed point is $M_1 = 0$ for any factor: $E_\nu$ at the centre of gravity of the occupied channel |
| $E_\nu$ band-window clamp (`clamp_enu`, 5 % margin) | 476 | checked in the logs: it fired in iterations 15 to 20 of 31 in the `edc` \ce{CuZn} run and never at convergence; never in the `ed3` \ce{CuZn} or the spin \ce{FeNiCr} run.  At convergence Zn d has $q_1 = -0.0002$, the unclamped fixed point |
| E_F search limits, DOS floor, 16 passes | 310 to 345 | the count is made exact afterwards (line 349); at convergence the step is $10^{-5}$ Ry or less (3.5) |
| Anderson safeguards and restarts | `asr_mixing` | $x_\mathrm{next} = x_\mathrm{in}$ when the residual is zero |
| start $E_F = \max E_{\nu,d} + 0.35$ | 274 | start value only |

If the clamp is ever active in the last iteration, the result depends on the 5 % margin.  The code prints a line each time.  A
converged run with that line in its last iteration must not be used.  This rule is not enforced by the code.

### 6.3 Approximations with a control parameter (these decide the accuracy)
| approximation | parameter, production value | converged? | evidence |
|---|---|---|---|
| node set: support of $H^{L}$ from the start orbital | `levels` 3 or 4 | not shown for energies; charges yes for \ce{CuZn} | SRO dE at $-0.266$: $-2.29$ (levels 3) against $-1.73$ (levels 4) |
| at most `rank` flipped sites in a configuration | `rank` 1 | see 6.4 | AMscf keeps up to 3 |
| Lanczos steps on that finite node set (hard wall) | `recursion_steps` 24 or 36 | **no**: it does not converge at fixed levels | CoCrFeNi: 36 and 48 steps overshoot (`STATUS.md`) |
| square-root tail from the extreme eigenvalues of the finite complex chain | `ed_terminator` hull | heuristic; not causal in general | F4 |
| early exit of a chain | `depth_tol` $10^{-5}$ | numerical tolerance | gate runs use 0 |
| kept range of $S^\alpha$ | `cluster_radius` 1.2 (2.0 for charges) | 1.2 is not converged: 0.01 e | gate 2 |
| contour | 32 points, $kT$ 1 mRy | yes to $10^{-3}$ for random alloys; not for ordered cells | gate 2 item 3 |
| one potential per component, equal spheres | none | model | SQS with one potential per species |
| screened Madelung of charge fluctuations | `sim_beta` 1 | model parameter | fitted value 1.0 for \ce{CuZn} only |
| first-shell SRO model | none | model | 3.1 |

Three of these need a plain statement.

1. **`recursion_steps` is not a convergence parameter.**  The node set is exact for the first `levels` steps only.  Later steps
   run on a finite, truncated space with a hard wall.  More steps resolve the spectrum of that truncated space.  They do not
   approach the alloy.  So "24 steps for strong moments, 36 for weak moments" is a tuned resolution, chosen where the result
   agrees with the reference.  The log already separates this ("termination artefacts beyond the exact region").  A converged
   result needs levels and steps to grow together.  That is not affordable now (cost blocker in `STATUS.md`).
2. **The sum-rule print is a weak guard.**  It is $1.000 \pm 0.005$ for \ce{CuZn} with 64 points, but 1.02 to 1.04 for
   \ce{FeNiCr} with 32 points, because its contour runs to 6 Ry.  So it cannot detect a loss of 1 % of the weight below the
   contour bottom ($\min E_\nu - 1.2$ Ry).
3. **`sim_beta` = 1 is a model value.**  It is the first-shell screening model.  The per-atom \ce{CuZn} cells give 1.0.  No
   HEA cell was used to check it.  The per-site CoCrFeNi cell has the data (charge against potential per site).

### 6.4 Rank test (new, desktop `~/asr_review/trial/desk1/rank`, binary c12eb4f)
Production uses `rank` 1.  The only evidence was \ce{CuZn} with 4-step chains (FINDINGS_LOG "Queue 2 results": rank 1 = 2 = 3
to three digits).  Production runs 24 to 36 steps.  FeRh was once suspected of a rank error and dismissed on that evidence.
Test: one iteration from one converged potential, 24 steps, `hull`, supercell 8, 32 contour points.  Only `rank` changes.

| case | quantity | rank 1 | rank 2 | rank 3 | nodes (rank 1 / 2 / 3) |
|---|---|---|---|---|---|
| \ce{CuZn} `edc`, levels 2 | qnet Cu | 0.1145 | 0.1119 | 0.1126 | 652 / 3089 / 6606 |
| \ce{CuZn} `edc`, levels 3 | qnet Cu | 0.1096 | 0.1096 | not run | 4555 / 58620 |
| \ce{CuZn} `edc`, levels 3 | E_F (Ry) | 0.002028 | 0.002419 | | |
| \ce{Fe66Ni14Cr20} spin, levels 2 | moments Fe / Ni / Cr ($\mu_B$) | 1.3772 / 0.2702 / $-0.9254$ | 1.3683 / 0.2659 / $-0.9030$ | 1.3689 / 0.2661 / $-0.9032$ | 2143 / 22958 / 95650 |
| \ce{Fe66Ni14Cr20} spin, levels 2 | qnet Fe / Ni / Cr | 0.0139 / 0.1807 / $-0.1722$ | 0.0079 / 0.1817 / $-0.1532$ | 0.0078 / 0.1815 / $-0.1529$ | |
| \ce{Fe66Ni14Cr20} spin, levels 3 | moments Fe / Ni / Cr ($\mu_B$) | 1.3406 / 0.2482 / $-0.9164$ | 1.3339 / 0.2487 / $-0.8930$ | not run | 15827 / 500607 |
| \ce{Fe66Ni14Cr20} spin, levels 3 | qnet Fe / Ni / Cr | 0.0245 / 0.1552 / $-0.1896$ | 0.0175 / 0.1547 / $-0.1661$ | | |

Reading.  These are one-iteration responses at a fixed potential, so self-consistency can enlarge them.
- The rank series converges at rank 2: rank 3 equals rank 2 to 0.05 % in the moments (\ce{FeNiCr}, levels 2).
- Rank 1 is not converged at levels 2: 2 % in the \ce{CuZn} charge, 0.6 / 1.6 / 2.4 % in the Fe / Ni / Cr moments, and 0.019 e
  (11 %) in the Cr charge.  Rank 2 lowers the moments.  This is the expected direction: rank 1 propagates in an averaged medium
  around one flipped site, so it broadens the bands too little.
- \ce{CuZn} at levels 3 shows no rank effect in the charge (below $10^{-4}$ e) and 0.4 mRy in E_F.
- \ce{FeNiCr} at levels 3 (production levels): rank 2 changes the Fe / Ni / Cr moments by $-0.5$ / $+0.2$ / $-2.6$ % and the Cr
  charge by 0.024 e (12 %).  dV after this one iteration is 0.014 Ry, against $4 \times 10^{-6}$ Ry in the rank-1 control.  So
  rank 1 is not converged for this alloy at production settings.  The self-consistent size is not measured.  A rank-2
  iteration costs about 1.9 h on 14 desktop threads (500,607 nodes, 11.5 GB).
- Consequence: for magnetic multicomponent alloys the published accuracy of rank 1 is not established.  The Cr moment and the
  Cr charge are the sensitive quantities.  FeRh (Fe 4.8 % low, flat in steps) should be rerun at rank 2.

### 6.5 Guards added (commit 36c10b1)
1. `E_\nu` clamp.  `class_scf` and `asa_scf` count the clamp events.  A run that reaches the dV tolerance with $E_\nu$ held at
   its window is reported as not converged, with a warning in the log and in `asr.out`.  The margin is the named variable
   `enu_window_margin` (0.05, unchanged).
2. Weight below the contour bottom.  `class_scf` integrates $G$ on a closed path around the 3 Ry below the contour bottom (20
   energies, about a third of one E_F pass) at the first and the final potential.  $G$ is analytic there when no state lies
   below, so the integral vanishes.  A value above $10^{-3}$ e per site is flagged.  This replaces the sum rule as the guard.
   The sum-rule print stays, labelled as a rough check.
3. `test_scf_guards` (9 checks; the first unit test that runs `class_scf`; fcc \ce{Cu}): the clamp is silent and the weight
   below is $10^{-8}$ e at production settings; a looser margin (0.01) gives the same bits; a margin of 0.999 is flagged and
   not converged; a contour bottom inside the band is flagged (4.8 e missing).
No default and no result changes.

### 6.6 `sim_beta` for the HEA, from the per-site cell (analysis, no run)
Data: the 40-atom CoCrFeNi cells with one potential per site (desktop `~/asr_ed3/trial/hea_sqs/site_full`, `site_h2`).  The
Madelung potential of every site was computed here from the positions and the converged charges (Ewald,
`scratchpad/hea/beta_fit.py`), then regressed on the charge.
| quantity | full ASA cell | $H^{(2)}$ cell |
|---|---|---|
| $\beta$ from the species means, $\langle V \rangle_s = -2 \beta \langle q \rangle_s / R_1$: Co / Cr / Ni | 1.17 / 1.21 / 1.21 | 1.17 / 1.21 / 1.20 |
| hardness $a = -dV/dq$ within one species (Ry/e): Co / Cr / Fe / Ni | 0.93 / 0.92 / 0.90 / 1.02 | 0.94 / 0.92 / 0.90 / 1.03 |
| $\beta$ that the model of `asr_screening` gives for $a$ = 0.90 to 0.96 Ry/e (fcc, this lattice constant) | 1.18 to 1.22 | |
| $\beta$ from one line through all 40 sites | 1.38 | 1.37 |
Fe has $\langle q \rangle = -0.005$ e, so its species-mean value is undefined.  The species means are the quantity that the
single-potential-per-species SCF uses.  So the cell gives $\beta = 1.2$ for fcc CoCrFeNi, and the charge-excess model predicts
the same value from the measured hardness.  The model is confirmed.  The input value 1.0 is 17 % too small for this alloy.
For bcc \ce{CuZn} the same analysis gave 1.0 (FINDINGS_LOG "B3 tested").  So `sim_beta` is alloy dependent, and 1.0 is not a
universal default.  At $\beta$ = 1 the ASR charges of CoCrFeNi are 7 to 12 % smaller than the cell's species means (Co 0.083
against 0.089, Cr $-0.211$ against $-0.239$, Ni 0.137 against 0.155).

Run at `sim_beta` 1.2 (desktop `~/asr_review/trial/desk2/crfeconi_b12`, 19 iterations from the $\beta$ = 1 potential, converged,
no clamp event).  Reference: the per-site cell with $H^{(2)}$ (`site_h2`, the same Hamiltonian as `edc`), species means.
| species | qnet $\beta$ 1 | qnet $\beta$ 1.2 | qnet cell | moment $\beta$ 1 | moment $\beta$ 1.2 | moment cell |
|---|---|---|---|---|---|---|
| Co | 0.0830 ($-7$ %) | 0.0947 (+6 %) | 0.0892 | 0.927 ($-4.2$ %) | 0.914 ($-5.5$ %) | 0.967 |
| Cr | $-0.2114$ ($-11$ %) | $-0.2402$ (+1 %) | $-0.2382$ | $-0.599$ (+22 %) | $-0.619$ (+26 %) | $-0.490$ |
| Fe | $-0.0082$ | $-0.0090$ | $-0.0064$ | 1.743 ($-2.8$ %) | 1.751 ($-2.3$ %) | 1.793 |
| Ni | 0.1367 ($-12$ %) | 0.1544 ($-1$ %) | 0.1554 | 0.272 ($-3.7$ %) | 0.266 ($-5.6$ %) | 0.282 |
So $\beta$ = 1.2 puts the charges on the cell (Cr and Ni within 1 %, Co within 6 %).  The moments do not improve: Co and Ni move
from 4 % to 5.5 % low, and Cr stays 22 to 26 % too large in magnitude.  The moment error is therefore not electrostatic.  It
comes from the recursion (rank, levels and steps; 6.3, 6.4).  The rank test points the right way: rank 2 reduces the Cr moment
of \ce{FeNiCr} by 2.6 to 3 %.  The `STATUS.md` row "CoCrFeNi levels 3 spin" (Co, Fe, Ni within 3 to 4 %) holds only at $\beta$ =
1, where the charges are 7 to 12 % off.  So that agreement was partly a compensation between electrostatics and the recursion.
A parameter-free form exists: the hardness is the on-site response $-dV_p/dq_p$ of the recursion itself, and $\beta$ follows
from it through `screen_lattice`.  This is not implemented.
