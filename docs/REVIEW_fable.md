# Review of code and physics (Fable 5.2 session, 2026-11-10)

Branch `review-fable`, off `main` e8d193c.  Not merged, not pushed.  Task: `docs/HANDOFF_review_fable.md`.
All work was done on the laptop (unit tests and small probes).  The desktop and smriti were not used.  smriti was not
reachable (connection timed out), so its result files were not read.

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
- Effect: total energies of the k-space full ASA and of `ed3` change.  Charges, $E_F$ and potentials do not.  The default
  route `edc` does not change.  **Not yet measured**: fcc \ce{Cu} $a_0$ and $B$ (k-space), and the \ce{CuZn} SRO energies of
  `ed3` and of the per-atom cells (`STATUS.md` row "bcc CuZn SRO").  Both sides of that comparison used the old form, so
  both must be recomputed.  The change of an SRO energy is $\Delta \sum_p w_p \sum_l E_\nu p^\gamma q_2$ between two SRO
  states.  It can be evaluated from the printed $q_2$, $p^\gamma$ and $E_\nu$ of the finished runs without new runs.
- Docs: `docs/DERIVATION_ed3.md` section 9 item 2, `doc/manual/theory/scf_electrostatics.md` ("Total energy").

### F2. Reduced operator: rounding leaves the symmetric subspace and is amplified (proven, fixed)
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
  with long chains are (gate 2 used 96 and 160 steps on B2).  Not measured on a production run.
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
- Consequence: the dips and negative values reported for the \ce{Fe66Ni14Cr20} SRO DOS have a known mechanism.  The
  production case itself was not reproduced (it needs the desktop).
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
| 2 | `ed3` at strong order | not resolved | F1 changes both sides of the energy comparison; levels-4 runs not read; see 3.1 |
| 3 | energy-dependent form | resolved | identity proven by algebra (section 3.2); code checked line by line; existing ring tests agree to $10^{-15}$ |
| 4 | `ed3` derivation | resolved with one correction | Results 1 to 5 re-derived; band energy corrected (F1); gates judged in 3.3 |
| 5 | total energy | narrowed | F1 fixed; functional form and double counting checked (3.4); 0.1 mRy claim not yet tested by a run |
| 6 | `ef_linear_max` | narrowed | same fixed point proven; size of the 1.5e-3 $\mu_B$ not attributed (3.5) |
| 7 | block kernels | resolved | not a kernel bug: the single-energy chain differs from itself by the same amount under a $10^{-15}$ shift of $z$; the persistent part was F2; block equals single to $10^{-14}$ at 144 steps after F2 |
| 8 | thread-order dependence | resolved | F3 |
| 9 | OpenBLAS ring test | narrowed | 3.6 |
| 10 | MPI shutdown crash | not reached | `asr_mpi.F90` read; no defect found by reading; no test possible here |
| 11 | terminators and causality | narrowed | F4: violations proven on a test operator; FeNiCr case not reproduced |
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
SCF stops on $dV$ only.  Whether 1.5e-3 $\mu_B$ is this term or the distance to the fixed point along a soft magnetic
mode needs the two logs (the printed "shift" of the last iteration) or a rerun at a tighter tolerance.

### 3.6 OpenBLAS ring test (doubt 9)
- The reference is well conditioned: $\| K^{-1} \| \le \max(D_s/D_c) / \mathrm{Im}\,z \approx 6$, condition number below 100.
- The recursion side calls no BLAS or LAPACK routine whose result it uses.
- Here: gfortran 16.2 with OpenBLAS 0.3.34 and with reference LAPACK (FlexiBLAS) pass, 40 of 40 random seeds, largest error
  $8 \times 10^{-15}$.  So the failure belongs to the gfortran 14 / OpenBLAS 0.3.28 build.  `RESULT.txt` on smriti decides
  between compiler and library; it was not read.
- The test uses an unseeded `random_seed()`.  It should get a fixed seed.  Not changed.
- Branch `lapack-interfaces`: it builds with gfortran 16 here, and that build was used for this test.  No defect seen.

### 3.7 Electrostatics (doubt 12)
$q = G b$, $V = W b$ with $G = a/(a + M)$, $W = a(1 - G)$.  Single-site elimination gives
$dq = (G_{00} + G_1 \hat A) db$ and $V = (W_{00} + W_1 \hat A) db$.  Both matrices are functions of $\hat A$ and commute, so
the code's $F$ equals $V dq^{-1}$.  $W_1 = -a G_1$ because $\langle \cos k R_1 \rangle = 0$.  Without SRO
$F = -2\beta/R_1$, with the sign of `scf_electrostatics.md` (electron excess lowers its own potential).  The single-site
energy is the exact average of the model: the fluctuation terms of $\tfrac{a}{2}(q - b)^2$ and of $\tfrac12 q V$ cancel
($+\mathrm{var}(V)/2a$ and $-\mathrm{var}(V)/2a$).

## 4. What needs the desktop or smriti
1. Effect of F1: fcc \ce{Cu} $a_0$ and $B$; \ce{CuZn} SRO energies of `ed3` and of the per-atom cells (from existing output).
2. Effect of F2 on a B2 `ed3` run at 96 steps, with `energy_block` 0 and 8 (doubt 7 on the production case).
3. \ce{Fe66Ni14Cr20} SRO DOS with spin: count of negative DOS points (doubt 11).
4. `ringtest/RESULT.txt` and the levels-4 results on smriti (doubts 9 and 2).
5. The two FeNiCr logs for doubt 6.

## 5. Tests
New: `test_band_energy` (4 checks), `test_augcls_project` (14 checks), two checks in `test_orbital_rep`.
Full suite on this branch at f19cb05 (ifx, laptop, 2026-11-10): 50 test programs, 435 checks, 0 failed.  The only later
source change is a comment in `asr_class_scf.f90`; the build passes.
