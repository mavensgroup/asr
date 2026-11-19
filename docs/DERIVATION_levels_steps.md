# `recursion_steps` at fixed `levels`: what converges, and what the short chain costs (derivation and tests)

Written 2026-11-16 (Fable 5.2).  Question of the handover: CoCrFeNi at levels 3 meets its reference at 24 steps, and 36 and
48 steps give larger moments.  Why, and how can levels and steps be converged together?  Statements marked **Result** are
proven here and checked numerically (section 5).  The production data are in section 3.  No production code is changed.
Tools: `tools/step_profile.py` (new), `tools/scf_response.py` (new), `tools/f4_dump.py`.

Short answers.
1. More steps do not run away.  At a fixed potential the moments converge with the chain length.  They are flat from
   about 64 steps (\ce{Fe66Ni14Cr20}) and 96 steps (CoCrFeNi) on.  For \ce{Fe66Ni14Cr20} the limit is the same at levels 2
   and levels 3 within 1 %.  For CoCrFeNi the levels-3 chain of 64 steps lies 2 to 5 % from the levels-2 limit and is not
   yet flat.
2. The chain of 24 to 36 steps is not converged.  In \ce{Fe66Ni14Cr20} its moments scatter by 3 to 5 % with the chain
   length (Ni up to 14 %), without a trend.  In CoCrFeNi the 24-step moments are low by 4 to 9 % against the limit, at
   both levels.  The SCF multiplies the fixed-potential error by the gain of the alloy (up to 7 for Fe, about 4 in
   CoCrFeNi).  This is the "overshoot" of 36 and 48 steps: the longer chain removes the low bias, and the SCF enlarges
   the change.
3. The limit of many steps at fixed levels is the Green function of the finite node set (Result S2).  It is exact for the
   first `levels` + 1 steps (Result S1).  The error of both truncations is one identity (Result S3).  In this alloy the
   node set is not the weak point.  The chain length is.
4. Error bar the code can print: the spread of the moments over the last chain lengths (free, from the same chains), times
   the gain matrix of `DERIVATION_scf_response.md` part C.  Recommendation in section 6.

## 1. Setting
Route `edc`.  $K(z) = zM - A$ on the full augmented space, with $M > 0$ diagonal and $M e_1 = e_1$ (`DERIVATION_f4.md`
section 1).  `levels` = $L$ keeps the node set $V_L$: every node that $L$ applications of $A$ reach from the start node.
One application of $H^{(2)}$ moves by up to two hops, so $V_L$ holds the sites within $2L$ hops and their flipped
configurations (at most `rank` flipped sites).  $P$ is the projector on $V_L$, $K_L = P K P$ the truncated operator, and
$G_L(z) = e_1^T K_L(z)^{-1} e_1$.  The code runs $m$ = `recursion_steps` steps of the complex-symmetric Lanczos chain of
$K_L$ and closes it with a tail.  $G_{m,L}$ is the tail-free value.  $G_\infty$ is the exact function of the full space
(for the given `rank`).

## 2. Results

**Result S1 (exact head).**  The coefficients $a_1 \dots a_{L+1}$ and $b_1 \dots b_L$ of the truncated chain are those of the
full operator at every $z$.

Proof.  $u_k \in \mathrm{span}\{K^j e_1, j < k\}$ has its support within level $k - 1$.  For $k \le L + 1$ this lies in $V_L$,
and $K u_k$ agrees with $K_L u_k$ on $V_L$.  $a_k = u_k^T K u_k$ needs $K u_k$ on the support of $u_k$ only.  $b_k$ needs the
norm of $K u_k$ on level $k$, which lies in $V_L$ for $k \le L$.  $\square$

So at levels 3 a chain of 24 steps has four exact $a$ and three exact $b$.  The other twenty levels are those of the node
set.

**Result S2 (limit of many steps).**  For $m \to \dim V_L$ the chain value tends to $G_L(z)$, with any tail.  $G_L$ is a
Herglotz function: a finite sum of real poles with positive weights that add up to 1.

Proof.  The tail-free value is the Galerkin approximant on the Krylov space of $K_L$ (F4 Result 1a), which grows to the
cyclic subspace of $e_1$.  The tail enters through $g_{1m}^2$ (F4 section 4), and $g_{1m} \propto \rho_m \to 0$.  $K_L$ is the
pencil $z M_L - A_L$ with $M_L > 0$ and $M_L e_1 = e_1$.  $\square$

So `recursion_steps` is a convergence parameter at fixed levels.  Its limit is the centre of the finite node set with a hard
wall, not the alloy.  Whether that limit is close to the alloy is the question of `levels` alone.

**Result S3 (one error identity for both truncations).**  Let $x_m \in V_L$ be the Galerkin solution of the $m$-step chain
and $r = e_1 - K x_m$ its residual in the full space.  Then $r = r_c + r_\ell$ with
$$r_c = -\rho_m\, q_{m+1} \in V_L \ \ (\text{chain residual}), \qquad r_\ell = (1 - P)\, A\, x_m \perp V_L \ \ (\text{leak out of the node set}),$$
$$G_\infty(z) - G_{m,L}(z) = r^T K(z)^{-1} r, \qquad |G_\infty - G_{m,L}| \le \frac{\|r_c\|_{M^{-1}}^2 + \|r_\ell\|_{M^{-1}}^2}{\mathrm{Im}\, z} .$$
(a) $r_\ell = 0$ for $m \le L$.  (b) For $m \to \infty$ at fixed $L$, $r_c \to 0$ and $r_\ell \to (1 - P) A K_L^{-1} e_1$, which
does not vanish.

Proof.  $M$ is diagonal, so $(1 - P) K P = -(1 - P) A P$.  Inside $V_L$ the Lanczos relation gives $P r = -\rho_m q_{m+1}$
(F4 Result 2).  $x_m^T r = 0$: $r_c$ is orthogonal to the Krylov space in the bilinear form, and $r_\ell$ to all of $V_L$.
Then $G_\infty - G_{m,L} = e_1^T K^{-1} r = (x_m + K^{-1} r)^T r = r^T K^{-1} r$.  The bound follows as in F4 Result 2.  The two
parts have disjoint supports, so the $M^{-1}$ norms add.  (a) $x_m$ lies within level $m - 1 < L$.  $\square$

Reading.  Steps reduce $r_c$.  Only levels reduce $r_\ell$.  Steps are converged for a given node set when $\|r_c\|$ is
below $\|r_\ell\|$.  As a bound the identity is as loose as F4 Result 2 near the real axis.  Its use is the structure:
there are two independent errors, and the data of section 3 say which one is large.

## 3. Production data: the step profile
One dump holds the chain coefficients of every contour node.  The first $m'$ coefficients are the $m'$-step chain.  So one
run at $m$ steps gives the charges and moments for every $m' \le m$ with the code's own `hull` tail
(`tools/step_profile.py`; $E_F$ is set again for every $m'$).  The check: $m'$ = 24 reproduces the 24-step run
(1.3410 / 0.2484 / $-0.9161$ against 1.3406 / 0.2482 / $-0.9164$).

Runs: desktop `desk14`, binary `asr_1dd6058`.  Random fcc \ce{Fe66Ni14Cr20}, spin, rank 1, `hull`, 32 contour points,
supercell 8.  One fixed potential for all: the converged levels-3, 24-step potential.  Moments Fe / Ni / Cr in $\mu_B$.

| chain length $m'$ | levels 3, $kT$ 1 mRy | levels 2, $kT$ 1 mRy | levels 2, $kT$ 5 mRy |
|---|---|---|---|
| 14 | 1.346 / 0.224 / $-0.903$ | 1.367 / 0.261 / $-0.902$ | 1.364 / 0.258 / $-0.900$ |
| 18 | 1.295 / 0.254 / $-0.921$ | 1.316 / 0.246 / $-0.920$ | 1.315 / 0.244 / $-0.914$ |
| 22 | 1.349 / 0.266 / $-0.913$ | 1.383 / 0.253 / $-0.913$ | 1.377 / 0.251 / $-0.909$ |
| 24 | 1.341 / 0.248 / $-0.916$ | 1.377 / 0.270 / $-0.925$ | 1.372 / 0.266 / $-0.925$ |
| 28 | 1.391 / 0.263 / $-0.899$ | 1.394 / 0.262 / $-0.904$ | 1.386 / 0.258 / $-0.902$ |
| 32 | 1.414 / 0.226 / $-0.912$ | 1.403 / 0.232 / $-0.913$ | 1.395 / 0.230 / $-0.909$ |
| 36 | 1.409 / 0.238 / $-0.916$ | 1.431 / 0.208 / $-0.916$ | 1.418 / 0.211 / $-0.915$ |
| 48 | 1.364 / 0.284 / $-0.896$ | 1.367 / 0.269 / $-0.900$ | 1.366 / 0.264 / $-0.900$ |
| 64 | 1.372 / 0.264 / $-0.903$ | 1.360 / 0.245 / $-0.899$ | 1.365 / 0.245 / $-0.901$ |
| 80 | 1.375 / 0.260 / $-0.900$ | 1.364 / 0.242 / $-0.903$ | 1.367 / 0.246 / $-0.909$ |
| 96 | 1.377 / 0.252 / $-0.904$ | 1.371 / 0.255 / $-0.900$ | 1.3670 / 0.2525 / $-0.9028$ |
| 128 | | 1.373 / 0.250 / $-0.898$ | 1.3682 / 0.2527 / $-0.9026$ |
| 160 | | 1.374 / 0.254 / $-0.893$ | 1.3679 / 0.2536 / $-0.9020$ |
| 192 | | 1.378 / 0.253 / $-0.895$ | 1.3681 / 0.2538 / $-0.9021$ |
| limit | 1.377(7) / 0.257(4) / $-0.904(3)$ | 1.372(6) / 0.254(6) / $-0.898(4)$ | 1.368 / 0.254 / $-0.902$ |

The limits are the means over $m' \ge 56$ (levels 3) and $m' \ge 96$ (levels 2), with the spread in the last digit.

Reading.
1. The moments converge with the chain length (Result S2).  At $kT$ = 5 mRy they are converged to four digits from 112
   steps.  At $kT$ = 1 mRy a spread of 0.5 % stays at levels 2 up to 192 steps.  That part is the discrete spectrum of the
   node set, seen at the lowest Matsubara node ($\pi kT$ = 3 mRy).  It is the only sign of the hard wall in these data.
2. The limits of levels 2 and levels 3 agree within 0.4 % (Fe), 1 % (Ni) and 0.6 % (Cr).  So $r_\ell$ of Result S3 is not
   the large term here.  The alloy scatters strongly, and the centre does not feel a wall at four hops.
3. Below about 40 steps the moments scatter with the chain length: Fe between 1.30 and 1.41 ($\pm 4$ %), Ni between 0.21 and
   0.28 ($\pm 14$ %), Cr within $\pm 1.5$ %.  The scatter has the same size at levels 2 and 3 and at $kT$ = 1 and 5 mRy.  It
   has no trend, and neighbouring lengths differ in sign.  It is the termination error of a short chain: the `hull` tail
   stands for all structure below the resolution of the chain, about $W/m$.
4. The 24-step value of levels 3 lies 2.6 % (Fe) and 3 % (Ni) below the limit and 1.3 % above it in $|m_{\mathrm{Cr}}|$.  The
   24-step value of levels 2 lies 0.4 % (Fe), 6 % (Ni) and 3 % (Cr) above.  Both are draws from the scatter.
5. Raising $kT$ does not remove the scatter of the short chain (hypothesis of F4 section 5a: not confirmed in that form).
   It removes the spread of the long chain.

CoCrFeNi, spin, rank 1, $\beta$ 1.0, the converged levels-3 24-step potential of `~/asr_ed3/trial/hea_steps` (desktop
`desk19`).  Moments Co / Cr / Fe / Ni in $\mu_B$.  The 24-step value of levels 3 reproduces the run (`hea_steps/ctl_s24`).

| chain length $m'$ | levels 3 | levels 2 |
|---|---|---|
| 24 | 0.927 / $-0.599$ / 1.743 / 0.272 | 0.963 / $-0.606$ / 1.793 / 0.269 |
| 36 | 0.944 / $-0.612$ / 1.804 / 0.256 | 0.933 / $-0.642$ / 1.824 / 0.231 |
| 48 | 1.005 / $-0.613$ / 1.802 / 0.284 | 1.009 / $-0.625$ / 1.818 / 0.277 |
| 64 | 0.983 / $-0.625$ / 1.814 / 0.286 | 0.987 / $-0.641$ / 1.820 / 0.294 |
| 96 | | 1.006 / $-0.646$ / 1.809 / 0.297 |
| 128 | | 1.002 / $-0.652$ / 1.818 / 0.302 |
| 192 | | 1.002 / $-0.643$ / 1.817 / 0.302 |
| limit | not reached | 1.005(4) / $-0.646(5)$ / 1.817(4) / 0.300(3) |

Here the short chain is biased, not only noisy.  At levels 3 the 24-step moments lie 8 % (Co), 7 % (Cr), 4 % (Fe) and 9 %
(Ni) below the levels-2 limit.  The per-site H^(2) cell has 0.967 / $-0.490$ / 1.793 / 0.282.  So the step-converged
moments at this fixed potential already lie above the cell for Co (+4 %), Cr (+32 %) and Ni (+6 %).

## 4. What the SCF makes of it
The SCF multiplies a fixed-potential change by the gain matrix of the alloy (`DERIVATION_scf_response.md` part C;
\ce{Fe66Ni14Cr20}: Fe 6.7, Ni 2.4, Cr 4.0 on the diagonal, Cr follows Fe with $-3.8$).  So a termination error of 3 % in
the Fe moment becomes about 15 to 20 % at self-consistency.  This one mechanism accounts for the following.

| observation (self-consistent, 24 steps unless stated) | size | reading |
|---|---|---|
| \ce{Fe66Ni14Cr20}: levels 2 against levels 3 | Fe 1.679 against 1.339 (+25 %) | the two 24-step values differ by 2.7 % at the fixed potential (1.377, 1.341); their limits agree |
| \ce{Fe66Ni14Cr20}: levels 4 against levels 3 | +5 / +9 / +3 % | fixed potential: +0.6 / +1.8 / $-0.2$ %, inside the scatter |
| CoCrFeNi levels 3: 36 against 24 steps | Co +12, Cr +14, Fe +11, Ni $-2$ % | fixed potential: +1.7 / +2.1 / +3.5 / $-5.9$ % (`hea_steps`), the size of the scatter |
| weak moments "need 36 steps" (Ni in \ce{Ni_{0.75}Cu_{0.25}}) | 38 % low at 24 steps | the weakest moment has the largest relative scatter |

Prediction (made before the run ended).  Corrected after the run: see section 7.  Take the limit of section 3 as the fixed-potential change of
\ce{Fe66Ni14Cr20} against the 24-step chain of levels 3: +0.034 / +0.008 / +0.014 $\mu_B$.  The gain matrix gives
+0.22 / +0.04 / $-0.07$ $\mu_B$.  So a step-converged SCF should end near 1.56 / 0.29 / $-0.99$ $\mu_B$, at levels 2 as at
levels 3 (linear response, error bar 30 % of the change).  Test: `desk20/FNC_L2_s96` (section 7).

Consequence for the comparison with the per-site cell.  CoCrFeNi at levels 3, rank 1 and 24 steps is within 4 % of the
cell for Co, Fe and Ni.  Rank 2 lowers the moments by 3 to 4 % (Co, Fe, Ni) and 16 to 18 % (Cr).  36 steps raise them by 11
to 14 %.  So the agreement at (rank 1, 24 steps) rests on two errors of opposite sign: rank 1 is too high, and this draw of
the short chain is too low.  The converged statement needs both: rank 2 and a step-converged chain.

## 5. Numerical check of S1 to S3 (laptop, `numpy`)
Simple cubic cluster of 1331 sites, random site energies and hoppings, random diagonal $M$ with $M_{11} = 1$, node set of
levels 3 (63 sites), three energies ($\mathrm{Im}\, z$ = 0.1, 0.02, 0.005).  S1: $a_1 \dots a_4$ and $b_1 \dots b_3$ agree with
the full chain to $10^{-10}$, and $a_5$, $b_4$ differ.  S3: the identity holds to $10^{-14}$ for $m$ = 3 to 40, $r_\ell = 0$
for $m$ = 3, and the bound holds.  S2: $|G_m - G_L|$ falls to $5 \times 10^{-4}$ at $m$ = 40 for $\mathrm{Im}\, z$ = 0.1, while
$\|r_\ell\|$ stays at 1.3.

## 6. Recommendation
1. Say in the paper that `recursion_steps` converges at fixed `levels` and that the short chain has a termination error of
   a few per cent at a fixed potential.  REVIEW 6.3 item 1 ("more steps resolve the truncated space, they do not approach
   the alloy") is half right: the limit is the node set, but for \ce{Fe66Ni14Cr20} the node set of levels 2 is already
   within 1 % of levels 3 (CoCrFeNi: 2 to 5 % at 64 steps, where the levels-3 chain is not yet flat).  How the paper presents this is the user's call.
2. Converge the chain first, then test levels.  At levels 2 a chain of 96 steps costs about half of levels 3 with 24 steps
   (2143 against 15827 nodes).  Rank 2 at levels 2 has 22958 nodes, against 500607 at levels 3.  So "levels 2, rank 2, about
   96 steps" may be the converged and affordable setting for the soft alloys.  This needs the runs of section 7 first.
3. Error bar from data the code has: the spread of the moments over the last chain lengths $m' = m - 8 \dots m$ at the
   converged potential (`tools/step_profile.py` on a three-energy dump), times the gain matrix (`tools/scf_response.py`).
   For \ce{Fe66Ni14Cr20} at 24 steps: spread 0.03 / 0.01 / 0.01 $\mu_B$, so $\pm 0.2$ / $\pm 0.03$ / $\pm 0.1$ $\mu_B$ at
   self-consistency.  The code could print the spread in every iteration at no cost in chains (a change of `cf_species`;
   not done, it needs the user's OK).
4. A cheap partial cure, also not implemented: average $G$ over the last chain lengths.  The mean of $G_{m'}$ over
   $m' = m - K + 1 \dots m$ keeps weight 1, is causal where every term is, and has the same limit as the chain.  On the
   profiles above the mean over $m'$ = 20 to 28 lies within 1 % (Fe), 3 % (Ni) and 1.2 % (Cr) of the limit, against 2.6 %,
   6 % and 3 % for single lengths.  It reduces the scatter by a factor 2 to 3.  It does not replace the long chain.
5. $kT$ = 1 mRy is below the level spacing of the node set of levels 2 (spread 0.5 % at 192 steps).  $kT$ = 5 mRy removes
   that spread.  A step-converged production run should check both.

## 7. Self-consistent test, and what is still open
\ce{Fe66Ni14Cr20}, spin, levels 2, rank 1, from the levels-3 potential (desktop `desk20`):

| steps | Fe | Ni | Cr | state |
|---|---|---|---|---|
| 24 (`desk4/cal/fnc_L2_r1`) | 1.679 | 0.333 | $-1.134$ | converged |
| 64 | 1.471 | 0.304 | $-0.975$ | stopped by the cap at iteration 40, dV $7 \times 10^{-5}$ |
| 96 | 1.515 | 0.332 | $-1.012$ | converged |
| levels 3, 24 steps | 1.339 | 0.248 | $-0.915$ | converged |
| levels 4, 24 steps | 1.406 | 0.270 | $-0.948$ | converged |

Reading.  The self-consistent moments move by 10 % (Fe) from 24 to 96 steps.  64 and 96 steps still differ by 3 % (Fe)
and 8 % (Ni).  That is the fixed-potential difference of these two lengths (0.8 % and 4 %) times the gain.  So the
self-consistent moments need 96 steps or more before they are stable to about 3 %.  The linear prediction of section 4
(1.56 / 0.29 / $-0.99$) is right for Fe and Cr within its error bar and wrong for Ni, which does not fall.  The levels-2
and levels-3 SCF values at 24 steps differ by 25 % only because their two 24-step chains differ at the fixed potential.

Open.
- A step-converged SCF at levels 3, to confirm that levels 2 and 3 agree self-consistently as they do at the fixed potential.
- CoCrFeNi: a step-converged rank-2 SCF at levels 2 (about 96 steps), and its comparison with the per-site cell.  The data
  above say that rank 1 with converged steps overshoots the cell for Co, Cr and Ni, and that rank 2 lowers the moments by
  3 to 18 %.  The 24-step agreement rests on these two errors of opposite sign.
