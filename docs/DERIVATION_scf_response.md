# Linear response of the SCF: rank amplification and the hardness (derivation)

Written 2026-11-13 (Fable 5.2).  Two analytic items of the handover: why self-consistency enlarges the rank effect, with a
cheaper way to reach rank 2 (part A), and a parameter-free hardness for `sim_beta` (part B).  Part A is proven and partly
supported by data.  Part B is derived and not tested.  No code is changed.

## A. Rank and self-consistency

### A.1 The amplification
Let $x$ be the SCF variable (potentials and $E_\nu$).  Let $F_r(x)$ be the output of one iteration at rank $r$.  The fixed
points are $x_r = F_r(x_r)$.  Define the defect $d(x) = F_2(x) - F_1(x)$ and the Jacobian $J_r = F_r'$.

**Result A1.**  To first order in the defect, $x_2 - x_1 = (1 - J_2)^{-1}\, d(x_1)$.

Proof.  $0 = F_2(x_2) - x_2 = F_1(x_1) + d(x_1) + J_2 (x_2 - x_1) - x_2 = d(x_1) - (1 - J_2)(x_2 - x_1)$.  $\square$

$d(x_1)$ is the one-iteration change at the fixed rank-1 potential.  This is what REVIEW 6.4 measured first.  The
self-consistent change is this defect times the response $(1 - J)^{-1}$ of the SCF.  The response is large in a direction
where $J$ has an eigenvalue close to 1.  In a magnetic alloy this is the longitudinal moment: a larger exchange splitting
gives a larger output moment.  The gain $1/(1 - \lambda)$ is the Stoner enhancement of the alloy at its working point.

Data, \ce{Fe66Ni14Cr20}, levels 3 (REVIEW 6.4).
| quantity | defect at fixed potential | self-consistent change | gain |
|---|---|---|---|
| moment Fe ($\mu_B$) | $-0.0067$ | $-0.0843$ | 12.6 |
| moment Ni ($\mu_B$) | $+0.0005$ | $-0.0099$ | sign change |
| moment Cr ($\mu_B$) | $+0.0234$ | $+0.1053$ | 4.5 |
| charge Cr (e) | $+0.0235$ | $+0.0009$ | 0.04 |

The charges are stiff: the electrostatic feedback has a large negative eigenvalue of $J$, so the gain is far below 1.  The
moments are soft.  A second data set shows the same softness without any change of rank.  The levels-2 rank-1 SCF from the
levels-3 potential (desktop `~/asr_review/trial/desk4/cal/fnc_L2_r1`, 2026-11-13) ends at moments 1.679 / 0.333 / $-1.133$.
The first iteration at that potential gave 1.377 / 0.270 / $-0.925$ (REVIEW 6.4), against 1.341 / 0.248 / $-0.916$ at
levels 3.  The gains are 9, 4 and 24.  This change is too large for a linear estimate, but it shows the same soft mode.

Consequences.
1. In this alloy a change of the recursion by 0.5 % at fixed potential becomes 5 % in the self-consistent moments.  This
   holds for every approximation of REVIEW 6.3: rank, levels, steps and the tail.  The 3 to 4 % target for the moments
   therefore asks for a fixed-potential accuracy of a few tenths of a per cent.
2. A comparison of two methods is cleaner at a common fixed potential.  It then measures the method.  The self-consistent
   comparison measures the method times the softness of the alloy.  Both numbers should be reported.
3. The same holds for the reference.  A supercell of the same alloy has the same soft mode, so its moments are equally
   sensitive to its own parameters (cell size, sampling, sphere radii).

### A.2 Rank 2 at the cost of a few rank-2 iterations
**Result A2 (two-level SCF).**  Define the outer iteration: $x^{(k+1)}$ is the fixed point of the map
$x \mapsto F_1(x) + d(x^{(k)})$.  Then
(a) every fixed point of the outer iteration is a rank-2 fixed point, with no approximation;
(b) the error obeys $e^{(k+1)} = (1 - J_1)^{-1} (J_2 - J_1)\, e^{(k)}$ to first order.

Proof.  (a) At the outer fixed point $x = F_1(x) + d(x) = F_2(x)$.  (b) Subtract $x_2 = F_1(x_2) + d(x_2)$ from
$x^{(k+1)} = F_1(x^{(k+1)}) + d(x^{(k)})$ and linearise.  $\square$

Each outer step costs one rank-2 iteration (for $d$) and one rank-1 SCF.  The soft mode is relaxed inside the rank-1 SCF,
where an iteration is 20 to 40 times cheaper.  The outer iteration converges when the rank correction depends weakly on
the potential.  The first outer step is Result A1 with the full nonlinear rank-1 response.

In the code the defect is a small object: the difference of the output moments $q_0, q_1, q_2$ per component, $l$ and spin
between the two ranks at the same input potential, each at its own $E_F$.  Both sets are neutral, so the defect carries no
net charge and the $E_F$ condition of the inner SCF is unchanged.

Cost estimate for CoCrFeNi at levels 3 (desktop, 10 threads; scaled from the \ce{FeNiCr} timings, not yet measured): one
rank-2 iteration about 7 h, one rank-1 SCF about 3 h.  Three outer steps are about 30 h.  The direct rank-2 SCF needs about
24 iterations, about 7 days.

Status.  Not implemented.  It needs an input option for a moment correction file and an output of the moments after one
iteration.  The test is \ce{Fe66Ni14Cr20} at levels 2, where the direct rank-2 SCF costs 2 to 3 h and is being run
(`desk4/cal/fnc_L2_r2`).

## B. A parameter-free hardness for `sim_beta`

### B.1 The model
`src/asr_screening.f90` uses the charge-excess model $q_i = b_i - (V_i - \mu)/a$ with one hardness $a$.  `sim_beta` fixes $a$
through $\beta = -W_{00} R_1 / (2 G_{00})$.  REVIEW 6.6 found $a$ = 0.90 to 1.03 Ry/e in the per-site CoCrFeNi cell and
$\beta$ = 1.2 from it.  So $a$ is a property of the alloy.  It can be computed.

### B.2 The bare on-site response from the chains that already exist
**Result B1.**  Shift the potential of one site by a constant $\delta V$.  In the `edc` form with this site as the frozen
centre the Green function of the centre becomes $G_L(z; \delta V) = 1 / (1/G_L(z) - \delta V)$ for every start orbital $L$.
The bare response of the site charge at fixed $\mu$ is
$$\chi_p = -\frac{\partial q_p}{\partial\, \delta V} = \frac{1}{\pi}\, \mathrm{Im} \sum_j w_j \sum_{L, \sigma} G_{p L \sigma}(z_j)^2 .$$
It needs no new chain.

Proof.  A constant shift moves $C$ and $E_\nu$ of the site together and leaves $\Delta$ and $\gamma$.  So $h$ is unchanged and
$H^{(2)} \to H^{(2)} + \delta V P_c$.  $T = 1$ on the centre, so $K \to K - \delta V P_c$.  The perturbation is a multiple of the
unit matrix on the centre.  The on-site block of $G$ is diagonal in the four orbital classes on a cubic site.  The Dyson
equation then closes on the centre and gives the first formula.  Its derivative at $\delta V = 0$ is $G_L^2$.  For one level
$G = 1/(z - \varepsilon)$ the formula gives $\chi = -f'(\varepsilon) > 0$.  $\square$

### B.3 The hardness
The site also responds inside its own sphere.  Let $U_p$ be the shift of the potential felt by the valence states per added
electron at a fixed external potential (Hartree and exchange-correlation inside the sphere).  An external shift $dV$ then
gives $dq = -\chi_p (dV + U_p\, dq)$.

**Result B2.**  $a_p = -dV/dq = 1/\chi_p + U_p$.

Limits.  $\chi_p \to 0$ (no states near $E_F$): $a \to \infty$, the charges are frozen, and $\beta \to 0$.  This is the
unscreened limit of the model.  A large $\chi_p$ gives $a \to U_p$.  The model needs $a > a_{\min}$ of the lattice
(`screen_t%a_min`).  $U_p$ is of the order of 1 Ry/e, above $a_{\min}$ (0.76 Ry/e for bcc \ce{CuZn}).

### B.4 What is needed for a test
1. $\chi_p$: one sum over the contour nodes of an iteration.  `tools/f4_dump.py` prints it from an `ef_scan_dump` file, so
   the F4 runs of 2026-11-13 give $\chi$ of Fe, Ni and Cr in \ce{Fe66Ni14Cr20} without a further run.
2. $U_p$: a finite difference in the sphere solver.  Add a small valence charge with the shape of the output density,
   rebuild the potential, and take the shift of the band centres $C_l$ weighted with the channel responses.
3. Compare $a_p$ with 0.90 to 1.03 Ry/e (CoCrFeNi cell) and the resulting $\beta$ with 1.2 (fcc CoCrFeNi) and 1.0 (bcc
   \ce{CuZn}).  The model has one $a$ for all sites.  The spread of $a_p$ over the species is its error bar.

Status.  Derived.  Not tested.  It becomes a default only after the test in item 3.
