# Causality of the truncated energy-dependent chain (finding F4): derivation

Written 2026-11-13 (Fable 5.2).  It answers the three F4 questions of `docs/REVIEW_fable.md`.  All statements marked
**Result** are proven here.  Each one was also checked on small real pencils (section 8).  The production check (chain dumps
of \ce{Fe66Ni14Cr20}) is set up in section 6.3.  Nothing in the production code is changed by this document.

Short answers.
1. Does the SCF see F4?  The contour sum does not need an analytic approximant.  It needs the approximant to be accurate at
   the nodes.  A non-causal pole changes the charge only through its ordinary residue.  It harms the SCF in one case: when
   it lies within about 10 kT of $E_F$ and below the contour line.  Then $N(\mu)$ is not monotone there (Results 9 and 10).
2. Is there a Herglotz termination inside the $H(E)$ form?  Not as a tail.  A tail acts on the last level only.  The exact
   condition under which every causal tail gives a causal $G$ is a $2 \times 2$ test on the head of the chain (Result 7).  It
   fails often.  Real band edges make the tail causal, but they do not make $G$ causal.  A causal completion of the chain
   data always exists (Result 6), but it needs a subspace that grows exponentially with the length.  The causal Krylov method
   for the same pencil is the Lanczos chain in the $M$ inner product.  That is the energy-independent route.
3. Can the violation be bounded and reported?  Yes, in three ways: an exact error identity with a rigorous bound (Result 2),
   the $2 \times 2$ test (Result 7), and the count of nodes with $-\mathrm{Im}\,G \le 0$.  Section 7 says which one the code
   should print.

## 1. Setting
Route `edc` (`src/asr_edc.f90`).  For a central species $c$ the secular operator is
$$K(z) = z - H(z) = zM - A, \qquad M = T^2,\quad A = T H^{(2)} T,\quad T_R = (\Delta_c/\Delta_R)^{1/2}.$$
$M$ and $A$ are real symmetric operators on the augmented space.  $M$ is positive definite.  The start vector $e_1$ is the
centre with the average configuration.  The centre is frozen to $c$, so $T = 1$ there and $M e_1 = e_1$.  The exact function
is $G(z) = e_1^T K(z)^{-1} e_1$.  Route `ed3` has the same structure with $M \to O'$ (positive definite, DERIVATION_ed3 Result 3)
and $A \to H'$.  Results 2, 4, 6, 7, 9 and 10 use only that $A$ and $M$ are real symmetric and $M > 0$.  So they hold for
`ed3` as well.  Results 1c and 3 also use $M e_1 = e_1$.

The exact $G$ is a Herglotz function.  With $A x_n = \varepsilon_n M x_n$ and $x_n^T M x_m = \delta_{nm}$,
$G(z) = \sum_n (e_1^T x_n)^2 / (z - \varepsilon_n)$.  The weights are positive and sum to $e_1^T M^{-1} e_1 = 1$.

The code runs the complex-symmetric Lanczos chain of $H(z)$ from $e_1$ for every contour node $z$.  The vectors are
normalised with the bilinear form $u^T u = 1$ (no complex conjugate).  The coefficients are $a_k(z) = u_k^T H(z) u_k$ and
$b_k(z)$.  The tridiagonal matrix is $T_m(z)$.  The code returns $G_m(z) = [(z - T_m(z) - t\, e_m e_m^T)^{-1}]_{11}$ with the
tail $t$ of `cf_species`.  I write $\tilde T_m(z) = z - T_m(z) = Q_m^T K(z) Q_m$ for the chain of $K$.  Its diagonal is
$z - a_k$ and its off-diagonal is $-b_k$.

## 2. What the m-step chain is

**Result 1 (structure).**
(a) $G_m(z)$ without a tail is the bilinear Galerkin approximant of $e_1^T K^{-1} e_1$ on the Krylov space
$V_m(z) = \mathrm{span}\{e_1, K e_1, \dots, K^{m-1} e_1\}$.
(b) With $s_j(z) = e_1^T K(z)^j e_1$, a real polynomial of degree at most $j$,
$G_m = \mathbf{s}^T \mathcal{H}^{-1} \mathbf{s}$, $\mathbf{s} = (s_0, \dots, s_{m-1})$, $\mathcal{H}_{ij} = s_{i+j+1}$.  So $G_m$ is a
rational function of $z$ with real coefficients.  Its poles are real or come in conjugate pairs.
(c) The first level does not depend on $z$: $a_1 = A_{11}$ and $b_1^2 = \sum_{j \ne 1} A_{j1}^2$.  The second vector $u_2$ is real
and fixed, and $z - a_2(z) = m_2 z - \tilde a_2$ with $m_2 = u_2^T M u_2 > 0$.  So $G_1$ and $G_2$ are Herglotz.
(d) The third vector is $u_3 \propto z v + w$ with the real vectors $v = (M - m_2) u_2$ and
$w = -(A - \tilde a_2) u_2 - \beta_1 e_1$, $\beta_1 = -u_2^T A e_1$.
$b_2^2 = z^2 |v|^2 + 2 z\, v \cdot w + |w|^2$ has a conjugate pair of zeros unless $v \parallel w$.  From $m = 3$ on the chain is
complex in an essential way.

Proof.  (a) $x_m = Q_m \tilde T_m^{-1} e_1$ solves $Q_m^T (K x_m - e_1) = 0$, and $G_m = e_1^T x_m$.  (b) Take the basis
$K^j e_1$.  The Galerkin matrix is $(K^i e_1)^T K (K^j e_1) = s_{i+j+1}$ because $K$ is symmetric.  (c) $M e_1 = e_1$ gives
$K e_1 = z e_1 - A e_1$, so $K_{11} = z - A_{11}$ and the part of $K e_1$ orthogonal to $e_1$ is $-(A e_1 - A_{11} e_1)$.  (d) follows
by one more step and the Cauchy-Schwarz inequality.  $\square$

**Result 2 (exact error and bound).**  Let $q_{m+1}$ be the next Lanczos vector and $\rho_m = b_m\, [\tilde T_m^{-1}]_{1m}$.  Then
$$G(z) - G_m(z) = \rho_m^2\; q_{m+1}^T K(z)^{-1} q_{m+1}, \qquad
|G(z) - G_m(z)| \le \frac{|\rho_m|^2\; q_{m+1}^H M^{-1} q_{m+1}}{\mathrm{Im}\, z}.$$
So the exact $G(z)$ lies in a disc around the tail-free value.  The radius is computable from the chain and one extra norm.

Proof.  The residual is $r_m = e_1 - K x_m$.  The Lanczos relation $K Q_m = Q_m \tilde T_m + b_m q_{m+1} e_m^T$ (up to the sign of
$b_m$) gives $r_m = -\rho_m q_{m+1}$.  Galerkin means $v^T r_m = 0$ for all $v \in V_m$.  Then
$G - G_m = e_1^T K^{-1} r_m = (x_m + K^{-1} r_m)^T r_m = r_m^T K^{-1} r_m$.  For the bound put $w = K^{-1} q$.  $A$ and $M$ are real
symmetric, so $\mathrm{Im}(w^H K w) = \mathrm{Im} z\; w^H M w$.  Also $w^H K w = w^H q$.  With the Cauchy-Schwarz inequality in the
$M$ norms, $\mathrm{Im} z\, \|w\|_M^2 \le \|w\|_M \|q\|_{M^{-1}}$.  So $\|w\|_M \le \|q\|_{M^{-1}} / \mathrm{Im} z$ and
$|q^T w| \le \|q\|_{M^{-1}} \|w\|_M$.  $\square$

Remarks.  $|\rho_m|$ falls like $(1 + \mathrm{Im} z / W)^{-m}$ for a band of half width $W$.  At the main line of the Fermi contour
($\mathrm{Im}\, z = 0.1$ Ry) 24 steps give a few per cent.  At the lowest Matsubara node ($\mathrm{Im}\, z = \pi kT = 0.003$ Ry) the
bound is of no use.  There the value of $G$ rests on the tail model.  This is the precise meaning of "`recursion_steps` is a
resolution" (REVIEW 6.3).  In the tests the bound is 15 to 400 times the true error.

**Result 3 (weight).**  If the remainder of the fraction below level 1 vanishes for $z \to \infty$ (the generic case), then
$z\, G_m(z) \to 1$.  The residues of $G_m$ sum to 1.  A non-causal pole moves weight.  It does not create or destroy weight.

Proof.  Result 1c: $G_m = 1 / (z - A_{11} - b_1^2 g_2(z))$ with constant $A_{11}$ and $b_1$.  $\square$

## 3. Why causality fails

**Result 4 (real axis).**  For real $E$ the chain is real.  Let $\kappa(E)$ be an eigenvalue of $\tilde T_m(E)$ with the
normalised eigenvector $x$.  It is a Ritz value of $K(E)$ on $V_m(E)$.  A real pole $E_p$ of $G_m$ is a zero of some $\kappa$.  Its
residue is $x_1^2 / \kappa'(E_p)$.  The slope is
$$\kappa'(E) = x^T (Q_m^T M Q_m)\, x + 2\, b_m x_m\, (c^T x), \qquad c = (dQ_m/dE)^T q_{m+1}.$$
The first term is positive and at least the smallest eigenvalue $\mu_{\min}$ of $M$.  The second term is the rotation of the
Krylov space with the energy.  Its factor $b_m |x_m|$ is the residual norm of the Ritz pair.  A pole of negative weight
needs $2\, b_m |x_m|\, |c^T x| > \mu_{\min}$.  So it needs a Ritz value that is not converged.  A Ritz pair with a small
residual gives a pole of positive weight.  When a pole of negative weight meets a pole of positive weight they leave the
real axis as a conjugate pair.

Proof.  $d\tilde T_m/dE = Q_m^T M Q_m + [\tilde T_m, \Omega] + b_m (c\, e_m^T + e_m c^T)$ with the antisymmetric
$\Omega = Q_m^T dQ_m/dE$.  The commutator has zero expectation in an eigenvector.  Near $E_p$,
$G_m \approx x_1^2 / (\kappa'(E_p)(E - E_p))$.  $\square$

For the full operator the second term is absent and every eigenvalue of $K(E)$ rises with $E$.  This is the Herglotz
property.  In an alloy all Ritz values inside the band are unconverged, because the spectrum is continuous.  So the
violation cannot be excluded there.  It is also the region that the tail is meant to model.

**Result 5 (minimal example).**  $N = 4$, $m = 3$, $M = \mathrm{diag}(1, 1, 1, 4)$,
$$A = \begin{pmatrix} 0 & 0 & 1 & 1 \\ 0 & 0 & -2 & 0 \\ 1 & -2 & -2 & 0 \\ 1 & 0 & 0 & 2 \end{pmatrix}.$$
The exact poles are $-3.4533$, $-0.2665$, $0.7290$, $1.4907$.  $G_3$ has the real poles $-2.9029$, $-0.2101$, $1.1694$ and the pair
$1.5551 \pm 0.2715\, i$ with the residue $-0.0827 \mp 0.0474\, i$.  The residues sum to 1.000000.  The three real poles carry
1.1654 and the pair carries $-0.1654$.  So F4 is a property of the exact arithmetic.  It is not a rounding effect.  $m = 3$ is
the smallest length at which it can occur (Result 1c).

**Result 6 (a causal completion always exists).**  Let $W_m$ be the real space spanned by all words of length at most $m - 1$
in $A$ and $M$ applied to $e_1$.  Let $(A_W, M_W)$ be the Galerkin projection of the pencil on $W_m$.  Then the $m$-step chain
of $(A_W, M_W)$ equals the $m$-step chain of $(A, M)$ at every $z$.  And $G_W(z) = e_1^T (z M_W - A_W)^{-1} e_1$ is Herglotz with
weight 1.

Proof.  $K^j e_1 \in W_m$ for $j \le m - 1$.  By induction the compressed operator gives the same vectors.  The coefficients
$a_1 \dots a_m$ and $b_1 \dots b_{m-1}$ are bilinear forms of these vectors with $K$.  $M_W$ is positive definite and
$M_W e_1 = e_1$.  $\square$

Consequence.  The chain data at all energies never contradict causality.  The violation comes from closing the chain
after $m$ steps at each $z$ separately.  The dimension of $W_m$ grows up to $2^m - 1$, so this completion is not a
method.

## 4. Tails (question 2)
A tail replaces the missing remainder by a number $t$ on the last level.  By the Sherman-Morrison formula
$$G(t) = g_{11} + \frac{g_{1m}^2\, t}{1 - g_{mm}\, t}, \qquad g_{ij} = [\tilde T_m(z)^{-1}]_{ij}.$$
Call a tail value causal when $\mathrm{Im}\, t \le 0$.  A square-root tail with real band edges is causal.

**Result 7 (Weyl disc test).**  Fix $z$ with $\mathrm{Im} z > 0$.  Every causal tail value gives $-\mathrm{Im}\, G > 0$ if and only if
the real symmetric matrix
$$\mathcal{W}_m(z) = -\mathrm{Im} \begin{pmatrix} g_{11} & g_{1m} \\ g_{1m} & g_{mm} \end{pmatrix}$$
is positive definite.

Proof.  Sufficiency.  Let $\mathcal{G}$ be the $2 \times 2$ corner of $\tilde T_m^{-1}$.  Then
$G(t) = [(\mathcal{G}^{-1} - t\, e_2 e_2^T)^{-1}]_{11}$.  For a complex symmetric matrix
$\mathrm{Im}(\mathcal{G}^{-1}) = -\mathcal{G}^{-1} (\mathrm{Im}\, \mathcal{G})\, \mathcal{G}^{-H}$.  So $\mathcal{W} > 0$ gives
$\mathrm{Im}(\mathcal{G}^{-1} - t e_2 e_2^T) > 0$ for $\mathrm{Im}\, t \le 0$, and the inverse of such a matrix has a negative
definite imaginary part.  Necessity.  If $-\mathrm{Im}\, g_{11} \le 0$, $t = 0$ fails.  If $\mathrm{Im}\, g_{mm} \ge 0$, the pole
$t = 1/g_{mm}$ of the map is a causal tail value.  Otherwise the real $t$ axis maps to a circle with the centre
$g_{11} + i g_{1m}^2 / (2\, \mathrm{Im}\, g_{mm})$ and the radius $|g_{1m}|^2 / (2 |\mathrm{Im}\, g_{mm}|)$.  Its top has the imaginary part
$\mathrm{Im}\, g_{11} + (\mathrm{Im}\, g_{1m})^2 / |\mathrm{Im}\, g_{mm}|$.  This is non-negative exactly when
$\det \mathcal{W} \le 0$.  $\square$

For a chain with real coefficients $\mathcal{W}_m > 0$ always holds.  These are the nested Weyl discs.  For the
energy-dependent chain the disc can leave the lower half plane.

**Result 8 (the candidates).**
1. *Tail band edges forced real.*  The tail is then causal.  By Result 7 this gives a causal $G$ only where
   $\mathcal{W}_m > 0$.  In the tests $\mathcal{W}_m$ fails at 26 % of the points.  The `hull` tail gives 19 negative values
   there, and the same tail with real edges still gives 11.  No proof is possible, because Result 7 is an equivalence.
2. *A tail from the real pencil, used only as terminator.*  It is one more causal tail value.  The same statement holds.
3. *A positive fit of the computed node values.*  By Pick's theorem the values $G(z_j)$ have a Herglotz interpolant if and
   only if the matrix $P_{jk} = -(G(z_j) - \overline{G(z_k)}) / (z_j - \bar z_k)$ is positive semidefinite.  The set of such data
   is convex and contains the exact values.  So the Euclidean projection onto it is causal by construction and never
   increases the distance to the exact values.  This is the only candidate with a proof.  But it only removes the symptom.
   In the tests the Pick violation is $4 \times 10^{-5}$ to $4 \times 10^{-3}$ of the largest eigenvalue, and the projection
   changes the electron count by less than the truncation error itself.
4. *What is causal by construction.*  Galerkin projection of the pencil on any fixed real subspace $V$ is Herglotz:
   $e_1^T (z M_V - A_V)^{-1} e_1$ with $M_V > 0$.  The choice $V = K_m(M^{-1} A, e_1)$ is optimal for all $z$ at once.  It is the
   Lanczos chain in the $M$ inner product, which is the energy-independent route with a frozen centre.  The ground rules
   exclude it for production.  For `ed3` it is not practical in any case, because it needs $O'^{-1}$ on a vector in every
   step.
5. *Pointwise positivity at full order.*  Hermitian Galerkin on the real span of $V_m(z)$ and its complex conjugate gives
   $-\mathrm{Im}\, G(z) > 0$ at every $z$ and keeps the second-order accuracy of the bilinear form.  It is not analytic in $z$.
   It needs the $2m$ real vectors in memory and has no tail theory.  In the tests it is not more accurate than the chain.
   It is not recommended.

So the answer to question 2 is no for tails, and yes only in forms that leave the chain (items 3 and 4).

## 5. What the SCF sees (question 1)

**Result 9 (the contour needs accuracy, not analyticity).**  The exact $G$ is analytic in the upper half plane.  So the
contour sum of the exact $G$ is the exact thermal electron count up to the quadrature error.  For any approximant $\tilde G$
$$\delta N = -\frac{1}{\pi}\, \mathrm{Im} \sum_j w_j \big(\tilde G(z_j) - G(z_j)\big).$$
The poles of $\tilde G$ do not enter.  Only its error at the nodes enters.  The sentence in the header of
`src/asr_contour.f90` ("analytic there, as the terminated continued fraction is") is wrong for the energy-dependent chain,
but the contour formula does not depend on it.

**Result 10 (the electron count of the tail-free chain).**  On the Fermi contour of `contour_nodes_fermi` (line at
$\mathrm{Im}\, z = \Gamma$, Matsubara nodes below it), with $f(z) = 1 / (e^{(z - \mu)/kT} + 1)$ and simple poles $z_p$ with
residues $r_p$,
$$N(\mu) = \sum_{p\ \mathrm{real}} r_p\, f(E_p) \;+\; 2\, \mathrm{Re} \sum_{0 < \mathrm{Im} z_p < \Gamma} r_p\, f(z_p).$$

Proof.  Close the path with the real axis.  The enclosed poles of $f G_m$ are the Matsubara poles (residue $-kT\, G_m$) and
the poles of $G_m$ (residue $r_p f(z_p)$).  The Matsubara terms cancel against the Matsubara nodes of the contour.  On the
real axis a conjugate pair contributes a real function.  $\square$

Consequences.
- A pole pair far below $\mu$ ($\mathrm{Re}\, z_p < \mu - 10\, kT$) has $f = 1$.  It enters with $2\, \mathrm{Re}\, r_p$, which is the
  weight it carries in the sum rule.  The charge is then exact for the approximant.  The error is the ordinary truncation
  error of Result 9.
- A pole pair far above $\mu$ does not enter.
- A pole pair within about $10\, kT$ of $\mu$ enters with the complex Fermi function.  Its phase is $\mathrm{Im}\, z_p / kT$, and
  $|f|$ is large next to a Matsubara pole.  Then $N(\mu)$ changes by the weight of the pair within a few $kT$ and can decrease.
  A real pole of negative weight does the same.  In the minimal example $N$ falls from 1.165 to 1.000 across the pair.
- A pole above the line ($\mathrm{Im}\, z_p > \Gamma$) is not counted.  Its weight is missing from the charge.  The CPC path test
  (moments change by at most $10^{-3}$ between $\Gamma$ = 0.09 and 0.4 Ry) bounds the weight between these two lines.

So the SCF sees F4 as a non-monotone $N(\mu)$ in a window of about $\pm 10\, kT$ = $\pm 0.01$ Ry around a non-causal feature
that lies below the contour line.  The stall of the `bp` tail at levels 4 had this form: a pole next to the contour and a
steep $N(E_F)$ (comment at `ed_tail` in `src/asr_endep.f90`).  Away from $E_F$ a non-causal feature is part of the truncation error and has no separate effect.  The
DOS notch of random \ce{Fe66Ni14Cr20} (Cr up d, $E - E_F = +0.11$ Ry) is 110 $kT$ above $E_F$.  It does not enter the charges.

## 5a. Consequence for `recursion_steps`
Result 2 gives the first rigorous measure of the chain length.  On the main line of the contour the tail-free error is
below $|\rho_m|^2\, q^H M^{-1} q / \mathrm{Im}\, z$, and it falls geometrically with the number of steps.  At the Matsubara
nodes next to the real axis no number of steps on a finite node set converges, because the chain then resolves the discrete
spectrum of that set (REVIEW 6.3 item 1).  The lowest node is at $\pi kT$ = 0.003 Ry.  A chain of $m$ steps on a band of half
width $W$ resolves about $W/m$, which is 0.02 Ry at 24 steps.  So at $kT$ = 1 mRy the electron count near $E_F$ rests on the tail
for every affordable length.  This suggests one way to make the length a convergence parameter: raise $kT$ with the
resolution, so that $\pi kT\, m / W$ stays above a fixed number, and extrapolate in $kT$.  This is a hypothesis.  The test is
the moments of CoCrFeNi against 24, 36 and 48 steps at $kT$ = 5 mRy.  If the overshoot at 36 and 48 steps disappears, the
hypothesis holds.

## 6. Tools and the production check
6.1 `tools/f4_dump.py` reads a file written by `[asr] ef_scan_dump` and reports, per scan energy: the nodes with
$-\mathrm{Im}\, G \le 0$ with the `hull` tail and without a tail, the nodes where $\mathcal{W}_m$ is not positive definite,
the nodes where the `hull` tail has $\mathrm{Im}\, t > 0$, the electron count and its increments, and the phase vortices of $G$
on the grid of scan energies and Matsubara frequencies.  A vortex of charge $-1$ is a pole and $+1$ is a zero.  A Herglotz
function has neither in the upper half plane.

6.2 Test of the tool.  A synthetic dump of the minimal example (Result 5) gives one pole vortex in the cell
$E \in [1.5551, 1.5611]$, $\mathrm{Im}\, z \in [0.264, 0.289]$.  The pole is at $1.5551 + 0.2715\, i$.  On a coarse grid (cells
larger than the distance between the pole and the nearest zero) the vortex is found one cell off.

6.3 Production runs (desktop `~/asr_review/trial/desk4/f4scan`, queued 2026-11-13 behind the calibration set).  Route `edc`,
levels 3, 24 steps, `hull`, converged potentials, 33 scan energies with step 0.002 Ry around $E_F$:
`fnc_rand_ef` (random \ce{Fe66Ni14Cr20}, spin), `fnc_sro_ef` (SRO triad (1, 1, 0)), and `fnc_rand_notch` (17 energies around
$E_F + 0.11$ Ry, the DOS notch).  Results: to be added.

## 7. What the code can report (question 3)
| quantity | cost | meaning | use |
|---|---|---|---|
| nodes with $-\mathrm{Im}\, G \le 0$ | none | a causal $G$ has none | print for `edc` as for `ed3`; warn when not zero |
| $\mathcal{W}_m(z)$ not positive definite | $O(m)$ per chain | the head does not guarantee causality for causal tails | verbose only; it is too strict as a warning (26 % in the tests against 0.3 % real violations) |
| disc radius of Result 2 | one norm and one $M^{-1}$ per chain | rigorous bound on the tail-free error | a convergence measure for `recursion_steps` on the main line |
| $N(\mu)$ increments (`ef_scan`) | one contour per energy | negative thermal DOS at $E_F$ | diagnostic when an SCF stalls |
| smallest eigenvalue of the Pick matrix | one small eigenproblem per chain | lower bound on the node error: $\max_j |\delta G_j| \ge -\lambda_{\min} / (2S)$, $S^2 = \sum_{jk} |z_j - \bar z_k|^{-2}$ | weak in the tests; not proposed |

## 8. Numerical checks (laptop, `numpy`; simple cubic clusters of 343 sites with three species)
| statement | result |
|---|---|
| Result 2, identity | holds to $10^{-14}$ at 16 points ($m$ = 8 to 36, $\mathrm{Im}\, z$ = 0.003 to 0.1) |
| Result 2, bound | holds everywhere; 15 to 400 times the true error |
| F4 in exact arithmetic | isolated $-\mathrm{Im}\, G \le 0$ points with and without full reorthogonalisation, identical counts |
| Result 4, slope formula | equals the finite difference to 5 digits; 1.4 % ($m$ = 12) and 2.7 % ($m$ = 24) of the Ritz values fall with $E$; their residuals are 5 to 7 times the median |
| Result 5 | rational function from the Hankel form equals the chain to $10^{-9}$; residues sum to 1.000000 |
| Result 6 | chain of the word-space pencil equals the full chain to $10^{-14}$ ($m$ = 3, 4, 5; dimensions 4, 8, 15) |
| Result 7, sufficiency | 6416 points, 200 random causal tails each: no violation where $\mathcal{W}_m > 0$; all 19 (`hull`) and 11 (real edges) negative values lie where $\mathcal{W}_m$ fails |
| Result 7, necessity | at all 119 failing points of one scan the constructed real tail gives $\mathrm{Im}\, G \ge 0$ |
| Result 10 | contour sum equals the pole formula to 6 digits at 7 values of $\mu$ across the pair |
| `hull` tail | $\mathrm{Im}\, t > 0$ at 39 % of the points |

## 9. Proposed changes (each needs the user's OK)
1. Correct the header of `src/asr_contour.f90` (comment only) as in Result 9.
2. Count the nodes with $-\mathrm{Im}\, G \le 0$ for route `edc` too and print a warning when the count is not zero.
3. Manual, "Causality of a truncated chain": replace "There is no proof of causality for any tail" by Results 7 and 10.
4. No change of the tail.  Real edges are not a fix (Result 8 item 1).
