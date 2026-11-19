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

## C. A predictor of the self-consistent change from one fixed-potential iteration

Written 2026-11-16 (Fable 5.2).  Part A gives $x_2 - x_1 = (1 - J)^{-1} d(x_1)$.  This part gives a model of $(1 - J)^{-1}$ in
the moments and charges of the components.  It uses only chains that the code already runs.  It has no fitted number.
Tool: `tools/scf_response.py` (input: an `ef_scan_dump` file with three scan energies and the `asr.out` of the same potential).

### C.1 Three response functions per component and spin
Let $q_{p\sigma}$ be the charge of component $p$ in spin $\sigma$.  The potential enters the chains of one spin only.

1. Uniform shift of the potential of spin $\sigma$ on all sites.  It equals a shift of $E_F$.  The response is the thermal
   density of states, $n_{p\sigma} = \partial q_{p\sigma} / \partial E_F = -\frac{1}{\pi}\, \mathrm{Im} \sum_j w_j\, G'_{p\sigma}(z_j)$.
2. Shift on the site alone.  Result B1: $\chi_{p\sigma} = \frac{1}{\pi}\, \mathrm{Im} \sum_j w_j \sum_L G_{pL\sigma}(z_j)^2$.
3. Shift on all other sites.  By 1 and 2 its response is $n_{p\sigma} - \chi_{p\sigma} = \frac{1}{\pi}\, \mathrm{Im} \sum_j w_j \sum_L
   (-G' - G^2)$.  This is exact.  It is negative when the site takes charge from its neighbours ($\chi > n$).

A shift on the other sites of one species $p'$ needs a closure.  In the single-site picture a species shift $\delta v_{p'}$
changes the effective medium by $\delta\sigma(z) = \sum_{p'} w_{p'} \kappa_{p'}(z)\, \delta v_{p'}$ with the vertex weights
$$\kappa_{p'}(z) = \frac{G_{p'}(z)^2}{\sum_{p''} w_{p''} G_{p''}(z)^2}, \qquad \sum_{p'} w_{p'} \kappa_{p'} = 1 .$$
This is the single-site (CPA) derivative of the coherent potential, $\partial\sigma / \partial\varepsilon_{p'} = w_{p'} D_{p'}^2 /
\sum w D^2$ with $D_p = G_p / \bar G$.  The average $\bar G$ cancels.  So
$$B_{pp'\sigma} = w_{p'}\, \frac{1}{\pi}\, \mathrm{Im} \sum_j w_j \sum_L \big(-G'_{pL\sigma} - G_{pL\sigma}^2\big)\, \kappa_{p'L\sigma}(z_j),
\qquad \sum_{p'} B_{pp'\sigma} = n_{p\sigma} - \chi_{p\sigma} .$$
The sum rule holds exactly (checked to $10^{-14}$).  Equal weights $\kappa = 1$ give the simpler closure
$B_{pp'} = (n_p - \chi_p)\, w_{p'}$.

$G'$ at a node comes from the scan: the nodes move with $E_F$, so $G' = [G(z_j^+) - G(z_j^-)] / (z_j^+ - z_j^-)$.  The nodes
that do not move lie far below the bands, where $G' = -G^2$.  The node sum of $-G'/\pi$ equals the scan difference of the
charges to $5 \times 10^{-3}$ states/Ry.

### C.2 The model
A change $(\delta m_p, \delta q_p)$ of the output moments and charges changes the input potential by
$$\delta v_{p\sigma} = -\sigma\, \tfrac12 I_p\, \delta m_p + \kappa_q\, \delta q_p, \qquad I_p = \frac{C_{d\downarrow} - C_{d\uparrow}}{m_p} .$$
$I_p$ is the Stoner parameter of the component.  The code prints $C_d$ and $m$ in `asr.out`.  The values are 0.068 (Fe),
0.082 (Ni), 0.061 (Cr) and 0.070 (Co) Ry/$\mu_B$ in every state looked at.  $\kappa_q$ is the hardness of the component
charge, about 1 Ry/e ($U$ of the sphere minus $2\beta/R_1$).  The output then changes by
$$\delta n_{p\sigma} = -\chi_{p\sigma}\, \delta v_{p\sigma} - \sum_{p'} B_{pp'\sigma}\, \delta v_{p'\sigma} + n_{p\sigma}\, \delta\mu .$$

**Result C1.**  Let $(\delta m^0, \delta q^0)$ be the change after one iteration at the fixed potential.  To first order the
self-consistent change solves the linear system
$$\delta m_p = \delta m^0_p + \sum_\sigma \sigma\, \delta n_{p\sigma}, \qquad \delta q_p = \delta q^0_p + \sum_\sigma \delta n_{p\sigma},
\qquad \sum_p w_p\, \delta q_p = 0 .$$
It has $2 n_c + 1$ unknowns ($\delta m$, $\delta q$, $\delta\mu$).

This is Result A1 with $J$ restricted to rigid shifts.  Limits.  (a) One component, no charge change: $\delta m = \delta m^0
/ (1 - I\, \tilde n)$ with $\tilde n = 2 n_\uparrow n_\downarrow / (n_\uparrow + n_\downarrow)$.  This is the Stoner
enhancement at a fixed charge.  (b) $\kappa_q \to \infty$: $\delta q = 0$, the charges are frozen.  For $\kappa_q \chi \gg 1$ the
charge change is $\delta q^0 / (\kappa_q \chi)$, about 0.04 $\delta q^0$, as measured (part A).  The moments do not depend on
$\kappa_q$ between 0.5 and 5 Ry/e (fourth digit).  (c) $B = 0$: every component has its own Stoner factor with $\chi$ in
place of $n$.

### C.3 Test
Data.  $n$, $\chi$, $B$ from chain dumps at the converged rank-1 potentials (desktop `desk4/f4scan/fnc_rand_ef`,
`desk15/fnc_L2`, `desk15/cfcn_L2`, `desk15/cfcn_L3`; levels as stated, 24 steps, `hull`).  \ce{Fe66Ni14Cr20}, levels 3:

| | Fe up / down | Ni up / down | Cr up / down |
|---|---|---|---|
| $n$ (1/Ry) | 12.18 / 10.99 | 7.46 / 7.20 | 8.84 / 9.25 |
| $\chi$ (1/Ry) | 10.82 / 15.56 | 6.00 / 8.36 | 11.50 / 12.82 |

Gain matrix $\partial\, \delta m / \partial\, \delta m^0$ (rows and columns Fe, Ni, Cr): $(6.7,\ 0.1,\ -1.0;\ 0.6,\ 2.4,\ -0.2;\ -3.8,\ -0.2,\ 4.0)$.
Fe is soft.  Cr follows Fe with the opposite sign.  Ni is stiff.

Changes in per cent of $|m_p|$.  "Fixed" is the measured change after one iteration.  "SCF" is the measured
self-consistent change.  The three closures are the predictions.

| case | species | fixed | on-site ($B = 0$) | equal weights | vertex weights | SCF |
|---|---|---|---|---|---|---|
| \ce{Fe66Ni14Cr20}, levels 3, rank 1 to 2 | Fe | $-0.5$ | $-5.4$ | $-4.1$ | $-5.1$ | $-6.3$ |
| | Ni | $+0.2$ | $-0.4$ | $-0.2$ | $-2.8$ | $-4.0$ |
| | Cr | $+2.6$ | $+10.5$ | $+12.1$ | $+13.6$ | $+11.5$ |
| \ce{Fe66Ni14Cr20}, rank 1, levels 3 to 4 | Fe | $+0.6$ | $+5.3$ | $+3.1$ | $+3.5$ | $+5.0$ |
| | Ni | $+1.8$ | $+7.1$ | $+6.9$ | $+8.0$ | $+8.6$ |
| | Cr | $+0.2$ | $+0.4$ | $-2.1$ | $-2.7$ | $-3.6$ |
| \ce{Fe66Ni14Cr20}, levels 2, rank 1 to 2 | Fe | $-0.8$ | $-4.8$ | $-4.2$ | $-5.5$ | $-6.3$ |
| | Ni | $-1.5$ | $-3.7$ | $-3.4$ | $-3.6$ | $-4.6$ |
| | Cr | $+2.7$ | $+10.4$ | $+12.2$ | $+14.9$ | $+13.7$ |

A positive change of Cr is a smaller magnitude (the Cr moment is negative).  The levels-2 line uses the change at the
levels-2 rank-1 fixed point (`desk15/fnc_L2_r2_it1`: moments 1.6653 / 0.3281 / $-1.1026$ against 1.6785 / 0.3332 / $-1.1334$).
The pair in `desk4/cal` started from the levels-3 potential, which is not that fixed point.

Reading.
- The vertex closure gives every sign.  This includes the reversal of Ni under rank 2 and the Cr change under levels 4.
  Both come from the coupling to Fe.  The on-site model has no such coupling and misses both.
- The magnitudes are within 1.5 points of $|m|$ in all nine entries (0.7 to 0.8 of the SCF change for Fe).
- The predicted charge changes are below 0.001 e.  The measured ones are below 0.0022 e.
- The error bar of the predictor follows from its own gain.  A gain $g$ turns a relative error $\epsilon$ of $I \chi$ into
  $g\, \epsilon$.  For Fe ($g$ about 7) an error of 3 % in $I$ or in the response gives 20 % in the predicted change.  The
  rigid-shift form of $I$ (the $d$ splitting for all channels) is of that accuracy.  So 20 to 30 % of the change is the
  honest error bar.  It is 1 to 2 points of the moment here.

### C.4 Prediction for CoCrFeNi, levels 3, rank 2
Input: the converged rank-1 potential (`desk2/crfeconi_b12`, $\beta$ 1.2; moments Co 0.9142, Cr $-0.6185$, Fe 1.7511,
Ni 0.2664) and the first iteration of the rank-2 run (`desk4/cfcn_L3_rank2`: 0.9137 / $-0.5934$ / 1.7372 / 0.2652; charges
0.0813 / $-0.2123$ / $-0.0210$ / 0.1519 against 0.0947 / $-0.2402$ / $-0.0090$ / 0.1544).  Dump: `desk15/cfcn_L3`.

| | Co | Cr | Fe | Ni |
|---|---|---|---|---|
| fixed-potential change (% of $|m|$) | $-0.1$ | $+4.1$ | $-0.8$ | $-0.5$ |
| predicted self-consistent change | $-3.5$ | $+18.5$ | $-4.0$ | $-3.7$ |
| predicted moments ($\mu_B$) | 0.882 | $-0.504$ | 1.680 | 0.257 |
| run at iteration 15 (dV $1.0 \times 10^{-4}$, not converged) | 0.8845 | $-0.5184$ | 1.6769 | 0.2591 |

The predicted charge changes are below 0.001 e.  The run was at iteration 15 when this was written, and its moments were
known.  The prediction does not use them.  Error bar: 30 % of each change, so Cr $-0.50 \pm 0.035$.

### C.5 Use
1. One rank-2 iteration at the rank-1 potential plus one three-energy dump at rank 1 give the self-consistent rank-2 moments
   to 1 to 2 points.  For CoCrFeNi at levels 3 this is one rank-2 iteration (5.4 h) and a dump (47 min) instead of about 20
   iterations (4.5 days).
2. The same gain matrix serves every small change of the method at that state: rank, levels, steps, the tail.  It turns a
   fixed-potential error into the self-consistent error (see `DERIVATION_levels_steps.md`).
3. The gain matrix is a property of the alloy state.  It says which moments are soft before a costly run is started.
