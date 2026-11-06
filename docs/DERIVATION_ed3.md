# Full-ASA energy-dependent recursion (route `ed3`): derivation

Branch `ed3`, written 2026-11-04.  Status: the operator algebra is proven.  The moment extraction is specified.  Gate 1 is
passed (2026-11-04, `test/test_asr_ed3.f90` 14/14, `test/test_asr_ed3_cls.f90` 9/9).  The production routine is not written
yet.  Open items are listed in section 9.  Main stays on the production route (`asr_edc`) until gates 2 to 4 pass.
Update 2026-11-06: the production routine is `src/asr_ed3.f90` with `[asr] route = "ed3"`.  The default route stays `edc`.

## 1. Why
The production route `asr_edc` is exactly the orthogonal second-order Hamiltonian $H^{(2)} = E_\nu + h - h o h$.  It has no
$h p h$ term.  The k-space reference solves the full generalised ASA problem.  The diagnostic `[kspace] second_order = true`
(`59b00fc`) gives k-space with the recursion's own Hamiltonian.  At identical free-atom potentials (bcc \ce{CuZn}, a = 2.85 A,
one iteration; FINDINGS_LOG and `trial/b45/h2`):

| quantity | k-space full ASA | k-space $H^{(2)}$ | ASR ($H^{(2)}$) |
|---|---|---|---|
| ordered B2: $q_\mathrm{Cu}$ (e), $E_F$ (Ry) | 0.1193, -0.1623 | 0.1267, -0.1124 | 0.1163, -0.1275 |
| 128-site SQS: $q_\mathrm{Cu}(0)$ | 0.1459 | 0.1776 | 0.1687 (24 steps), 0.1772 (36 steps) |
| SRO hopping charge at $\alpha_\mathrm{CuZn} = -0.148$ | -0.0050 | -0.0103 | -0.0111 (24), -0.0102 (36) |

At $\alpha = 0$ the augmented space reproduces $H^{(2)}$ to about 1 %.  $H^{(2)}$ doubles the SRO hopping response of the
correlated cells.  So the recursion needs the full generalised problem.

The last row compares two structures.  The ASR column uses the first-shell Mookerjee-Prasad model.  The k-space columns use the
correlated cells.  So its agreement between -0.0102 and -0.0103 does not show that the model is exact.  The exact model value in
the full ASA is -0.0072 +- 0.0006 e (FINDINGS_LOG "MP reweighting"), 1.45 times the correlated value.  The route here keeps that
model, so it cannot remove this factor.

The exact model value with $H^{(2)}$ is now measured (2026-11-06, `trial/b45/mp2`; the same 10 cells as `trial/b45/mp`, with
`second_order = true`).  It is the binomial average of a quadratic fit, with a statistical error of about 10 %.

| $\alpha_\mathrm{CuZn}$ | exact model, $H^{(2)}$ | ASR (levels 4, 24 steps) | exact model, full ASA | correlated cells, full ASA |
|---|---|---|---|---|
| -0.148 | -0.0112 | -0.0111 | -0.0072 | -0.0050 |
| +0.148 | +0.0111 | +0.0111 | +0.0071 | +0.0050 |
| -0.266 | -0.0201 | -0.0199 | -0.0131 | -0.0083 |

So the ASR reproduces its own model at $H^{(2)}$ within 1 % at all three $\alpha$.  The old factor of 1.55 between the ASR and
the full-ASA model is therefore all $H^{(2)}$.  The route has no hidden error here.  The target of this route is -0.0072.  The
correlated cells with $H^{(2)}$ give -0.0103, so the ratio of model to correlated value is 1.08 at $H^{(2)}$ and 1.45 in the full
ASA.  So $H^{(2)}$ hides most of the model error, and this route will expose it.

## 2. Notation (as in the code)
All potential parameters are diagonal in site $R$ and orbital $L$.  The screening $\alpha$ is the standard TB-LMTO set.  It is
the same for every species.  So the screened structure constants $S = S^\alpha$ do not depend on the configuration.
- $E_\nu$, $c = C^\alpha - E_\nu$, $D = \Delta^\alpha$, $d = D^{1/2}$, $o = o^\alpha$ (`o_alpha`), $p = p^\alpha$ (`to_alpha`).
- `to_alpha` (src/asr_lmto_ham.f90:126) gives $p^\alpha = p^\gamma - (o^\gamma)^2 + (o^\alpha)^2$.  From the gamma
  representation $o^\gamma = 0$.  So $p = p^\gamma + o^2$.  Here $p^\gamma = \langle \dot\phi | \dot\phi \rangle > 0$.
- $h = c + d S d$.  It is real symmetric.
- $\varepsilon(E) = E - E_\nu$, diagonal.
- Generalised problem of `bands_core` (src/asr_asa.f90, `zhegv` branch):
  $O = 1 + o h + h o + h p h$ and $H = E_\nu + h + E_\nu o h + h o E_\nu + h (o + p E_\nu) h$.  In the code the linear term is
  $(\tfrac12 + E_{\nu,i} o_i + \tfrac12 + E_{\nu,j} o_j)\, h_{ij}$: one $h$, not two.

## 3. The secular matrix and its factorisation
**Result 1.**  $K(E) = E O - H = \varepsilon - h + \varepsilon o h + h o \varepsilon - h o h + h p \varepsilon h$.

Proof.  The zeroth order is $E - E_\nu = \varepsilon$.  The terms linear in $h$ are $E (o_i + o_j) h_{ij} - (1 + E_{\nu,i} o_i +
E_{\nu,j} o_j) h_{ij} = -h_{ij} + \varepsilon_i o_i h_{ij} + h_{ij} o_j \varepsilon_j$.  The quadratic terms are
$E\, h p h - h (o + p E_\nu) h = -h o h + h p \varepsilon h$, because the middle factors are diagonal.  Every ordering is
kept.  Per-orbital $E_\nu$ is allowed.

**Result 2 (factorisation).**  With $p = p^\gamma + o^2$,
$$K(E) = (1 + h o)\,[\varepsilon (1 + o h) - h] + h\, p^\gamma \varepsilon\, h .$$
Proof.  Expand the product: $\varepsilon + \varepsilon o h - h + h o \varepsilon + h o \varepsilon o h - h o h$.  The term
$h o \varepsilon o h$ equals $h\, o^2 \varepsilon\, h$.  Adding $h p^\gamma \varepsilon h$ gives $h (o^2 + p^\gamma) \varepsilon h
= h p \varepsilon h$.  This is Result 1.

**Result 3.**  $O = (1 + h o)(1 + o h) + h p^\gamma h$.  Both terms are positive semidefinite.  So $O$ is positive definite
when $1 + o h$ is invertible.  Then the eigenvalues of the pencil are real.  Contour integration of $G = K^{-1}$ is valid.

## 4. Central-species scaling
Fix the central site 1 to species $c$.  Let $T = \mathrm{diag}((D_c / D_R)^{1/2})$.  $D_c$ is one number per orbital.  So
$T d = D_c^{1/2}$ on every site.  $T_1 = 1$.

**Result 4.**  $h' = T h T = g + S'$ with $g_R = (D_c / D_R)\, c_R$ (random, diagonal) and $S' = D_c^{1/2} S D_c^{1/2}$
(not random).  Then
$$K'(E) = T K T = \tilde\varepsilon - h' + u h' + h' u - h' \tilde Q h' ,$$
$$\tilde\varepsilon_R = \frac{D_c}{D_R}\,\varepsilon_R, \qquad u_R = o_R\, \varepsilon_R, \qquad
\tilde Q_R = \frac{D_R}{D_c}\,(o_R - p_R\, \varepsilon_R) = q_R - \tilde p_R\, \varepsilon_R .$$
$q_R = o_R D_R / D_c$ is the `asr_edc` value.  $\tilde p_R = p_R D_R / D_c$.

Proof.  $T$ is diagonal.  So $T (\varepsilon o h) T = \varepsilon o\, h'$ and $T (h X h) T = h' (T^{-1} X T^{-1}) h'$ for any
diagonal $X$.  Use $T^{-2} = D_R / D_c$.

Properties.
- Every random factor is site diagonal and sits around the non-random $S'$.  This is the structure the augmented space
  needs.  It holds for any number of species and sublattices, for spin (per channel), and for per-orbital $E_\nu$.
- Every coefficient is linear in $E$.  The operator has no poles.  The pole problem of the old $P(z) - S$ route
  (`asr_endep`) cannot occur.
- $K'$ is symmetric.  At complex $E$ it is complex symmetric.  The site-diagonal factors depend on the species and on $l$
  only.  So the point-group reduction of `asr_augcls` and the SRO shell classes apply unchanged.
- $\langle 1 | K'^{-1} | 1 \rangle = \langle 1 | K^{-1} | 1 \rangle = G_{11}$, because $T_1 = 1$.

**Canonical form** (for the implementation).  Insert $h' = g + S'$:
$$K' = D^\mathrm{diag} - (V S' + S' V) - S' \tilde Q S' ,$$
$$D^\mathrm{diag}_s = \frac{D_c}{D_s}\,[\varepsilon_s - c_s + 2 o_s c_s \varepsilon_s - c_s^2 (o_s - p_s \varepsilon_s)], \qquad
V_s = (\tfrac12 + c_s o_s) - (o_s + c_s p_s)\,\varepsilon_s .$$
The single $-S'$ term splits into two halves.  So $V_s$ starts with $\tfrac12$, not 1.  The $D_c / D_s$ factors cancel in
$g \tilde Q = c (o - p \varepsilon)$.  So $V$ has no such factor.

**Matrix-vector product.**  With $H_\mathrm{eff} = E - K'$ and $e = E - \tilde\varepsilon$ (the `asr_edc` diagonal):
$$H_\mathrm{eff}\,\Psi = e\,\Psi + (1 - u)\, h'\Psi - h'\,[\,u\,\Psi - \tilde Q\, h'\Psi\,] .$$
It needs two applications of $S'$ per Lanczos step.  The current route also needs two.

## 5. The moments a sphere charge needs
In real space $G = K^{-1} = \sum_\alpha z_\alpha z_\alpha^T / (E - E_\alpha)$ with $z_\alpha^T O z_\alpha = 1$.
`bands_moments` uses, per site and orbital, the occupied sums of $\mathrm{xm0} = z_1^2$, $\mathrm{xm1} = z_1 (h z)_1$ and
$\mathrm{xm2} = (h z)_1^2$.  So the recursion must give three central blocks:
$$G_{11},\qquad (h G)_{11},\qquad (h G h)_{11} .$$
Then, as in `bands_moments` and `asa_scf`:
- $q_0 = \mathrm{xm0} + 2 o\, \mathrm{xm1} + o^2 \mathrm{xm2}$, $q_1 = \mathrm{xm1} + o\, \mathrm{xm2}$, $q_2 = \mathrm{xm2}$;
- density $q_0 \phi^2 + 2 q_1 \phi \dot\phi + q_2 \dot\phi^2$;
- sphere charge $q_0 + p^\gamma q_2 = \mathrm{xm0} + 2 o\, \mathrm{xm1} + p\, \mathrm{xm2}$ (with the code's $p$).

The sphere integral uses $\int \phi^2 = 1$, $\int \phi \dot\phi = 0$ and $\int \dot\phi^2 = p^\gamma$.  The current route
uses the density $M_0 \phi^2 + 2 M_1 \phi\dot\phi + M_2 (\dot\phi^2 + \phi \ddot\phi)$.  That form integrates to $M_0$.  Both
forms are consistent once the moments mean the same thing.  Here they do not: in `ed3` the moments are matrix elements, not
energy moments of one Green function.

## 6. Exact extraction: two chains per energy
Since $G = T K'^{-1} T$, $h T = T^{-1} h'$ and $T_1 = 1$:
$$(h G)_{11} = \langle v | K'^{-1} | u \rangle, \qquad (h G h)_{11} = \langle v | K'^{-1} | v \rangle, \qquad
|u\rangle = |1 L\rangle, \quad |v\rangle = h' |u\rangle = g_c |u\rangle + |w\rangle, \quad |w\rangle = S' |u\rangle .$$

**Result 5.**  $|w\rangle$ does not depend on the configuration.  $g_c$ belongs to the frozen centre and $S'$ is not random.
So $|w\rangle$ lies in the reference sector of the augmented space (no flipped site).  $S'$ commutes with the site group.
So $|w\rangle$ transforms like $|1 L\rangle$ and lives in the same reduced basis.

Recipe per (central species, start orbital, contour energy):
1. Chain from $|u\rangle$: $G_{uu}$ from the continued fraction with the `hull` tail, as now.
2. In the same chain, $G_{wu} = \sum_m (w^T v_m)\, [J^{-1}]_{m1}$.  Here $v_m$ are the Lanczos vectors and $J$ is the
   terminated tridiagonal matrix.  The bracket is the bilinear form of the complex-symmetric Lanczos.  The dot products
   are cheap, because $w$ has support on the first shell only.
3. Chain from $|w\rangle$: $G_{ww}$ from its continued fraction.  The same projection onto $u$ gives $G_{uw}$.  By symmetry
   $G_{uw} = G_{wu}$.  This gives a free consistency check of step 2.
4. $(h G)_{11} = g_c G_{uu} + G_{wu}$ and $(h G h)_{11} = g_c^2 G_{uu} + 2 g_c G_{wu} + G_{ww}$.

Cost: two chains instead of one, each with two $S'$ applications per step.  So about twice the current route.

**Production choice (2026-11-06): chains from $|a\rangle$ and $|v\rangle$.**  The recipe above is exact, but its sphere charge
depends on the cross term.  A diagonal Lanczos element has an error of second order in the residual.  A projected off-diagonal
element has an error of first order.  With $|a\rangle = |u\rangle + o_c |v\rangle$ the gamma-representation moments are
$$q_0 = \langle a | K'^{-1} | a \rangle, \qquad q_1 = \langle v | K'^{-1} | a \rangle, \qquad q_2 = \langle v | K'^{-1} | v \rangle .$$
Both vectors lie in the reference sector, as in Result 5.  So $q_0$ and $q_2$ are diagonal elements of two chains.  The sphere
charge $q_0 + p^\gamma q_2$ needs no cross term, and it is a positive sum of two Herglotz functions.  Only $q_1$ needs the
projection.  The two chains give it twice, and the code uses the mean.  The cost is the same two chains.

Measured before any gate run, on the d channels of bcc CuZn parameters with a conditional shell (levels 4, 1414 nodes, 20 energies).
Mean relative error of the charge combination against 240-step chains:

| steps | $\mathrm{Im}\,z$ | chains $u, w$ | chains $a, v$ |
|---|---|---|---|
| 24 | 0.03 | 6.8e-3 | 5.4e-3 |
| 24 | 0.10 | 8.1e-4 | 1.5e-4 |
| 24 | 0.30 | 8.2e-5 | 1.4e-7 |
| 36 | 0.03 | 1.9e-3 | 1.1e-3 |
| 36 | 0.10 | 1.7e-4 | 2.0e-5 |
| 36 | 0.30 | 4.4e-7 | 2.0e-10 |

The two choices agree at 240 steps to 6e-7.  With chains from $u, w$ the block matrix of $-\mathrm{Im}\,G$ failed to be positive
semidefinite on terminated chains at 8 of 72 energies at $\mathrm{Im}\,z = 0.01$ and 24 steps.  With chains from $a, v$ it held at
all 20 test energies from $\mathrm{Im}\,z = 0.01$ to 0.3 (`test_asr_ed3_prod`).  The two sets of energies differ, so this is
not a strict comparison.

Why one chain is not enough.  Row 1 of $K' G' = 1$ reads
$\tilde\varepsilon_1 G_{11} - (1 - u_1)(h' G')_{11} - (h' u G')_{11} + (h' \tilde Q h' G')_{11} = 1$.
The last two terms carry $u_R$ and $\tilde Q_R$ of the neighbours.  They are new unknowns.  So this identity does not fix
$(h G)_{11}$ from $G_{11}$.  We know no other closure for $p^\gamma \neq 0$.

## 7. Single-site dressing (variant to measure, not for production)
**Result 6.**  If $p^\gamma = 0$, Result 2 gives $K = (1 + h o)[\varepsilon (1 + o h) - h]$.  On every eigenvector
$h z = \varepsilon (1 + o h) z$, since $1 + h o$ is invertible.  On the centre this gives
$(h z)_1 = \varepsilon_1 z_1 / (1 - \varepsilon_1 o_1)$.  Then $\mathrm{xm1}$ and $\mathrm{xm2}$ follow from $G_{11}$ with the
single-site factors $\varepsilon_1 / (1 - \varepsilon_1 o_1)$ and its square.  One chain would be enough.

With $p^\gamma \neq 0$ the defect is $r = h z - \varepsilon (1 + o h) z = (1 + h o)^{-1} h\, p^\gamma \varepsilon\, h z$.  Its size
relative to $(h z)_1$ is of order $p^\gamma \varepsilon^2$.  This is not small for d bands.  The B2 run prints
$p^\gamma = 2.27$ Ry$^{-2}$ for Cu d.  There $p^\gamma q_2 = 0.142$ e.  For $\varepsilon = 0.1$ to 0.2 Ry the defect is 2 to 9 %.
Through $p\, \mathrm{xm2}$ and $2 o\, \mathrm{xm1}$ this is about 0.01 to 0.03 e.  That is the size of the effects to be fixed.
So the variant is only a candidate.  Gate 1 must measure its error.

The factor $1/(1 - \varepsilon_1 o_1)$ has a pole at $E_\nu + 1/o_1$.  With the standard $\alpha$ and the B2 parameters,
$o_d \approx 1.0$ Ry$^{-1}$ for Cu and Zn.  So the pole sits about 1 Ry above $E_\nu$ (Zn d about +0.2 Ry, Cu d about +0.46 Ry).
For s, $o \approx -0.43$ Ry$^{-1}$, so the pole sits about 2.3 Ry below $E_\nu$.  Check both against the actual contour
before using the variant.

**Second-order limit in closed form.**  At $p^\gamma = 0$ Result 2 gives $K^{-1} = [\varepsilon (1 + o h) - h]^{-1} (1 + h o)^{-1}$.
So $(1 + o h)\, G\, (1 + h o) = [\varepsilon - h (1 + o h)^{-1}]^{-1}$ exactly.  Its site-1 block is
$G_{11} + o_1 (h G)_{11} + (G h)_{11}\, o_1 + o_1 (h G h)_{11}\, o_1$, the $q_0$ combination.  Expanding
$h (1 + o h)^{-1} = h - h o h + \dots$ gives $[E - H^{(2)}]^{-1}$.  So $H^{(2)}$ is the second-order truncation of the
physical-amplitude resolvent.  Checked here.

## 8. Limits
- **Ordered limit.**  With one species per sublattice the augmented space is the reference sector only.  Then $K$ is the
  `bands_core` pencil term by term.  The ASR differs from k-space only by the $S$ range (`rcl` 1.2 alat against $4w$ in
  k-space) and the chain truncation.
- **Second-order limit.**  It is not term by term.  At $p = 0$ the new route is the generalised problem with first-order
  overlap in the non-orthogonal basis.  $H^{(2)}$ is its orthonormal form to second order.  With $u = 0$ and $p = 0$,
  $K' = \tilde\varepsilon - h' - h' q h'$, while `asr_edc` has $\tilde\varepsilon - h' + h' q h'$.  The routes agree only in the
  physical moments, to second order.  Compare them numerically.
- **No SRO.**  The operator does not depend on $\alpha_\mathrm{SRO}$.  Only the shell class concentrations do.  At
  $\alpha_\mathrm{SRO} = 0$ they are $x$.  Homogeneity $\sum_c x_c X^{(c)} = x$ is unchanged.
- **Electron count.**  Over all sites, $\sum_{E_\alpha \le E_F} z_\alpha^T O z_\alpha$ is the number of occupied states.  So the
  total count from the sphere charges matches k-space.
- **Per-orbital weight.**  It is not 1.  $\sum_\alpha Q_1 = [O^{-1} + o_1 (h O^{-1} + O^{-1} h) + p_1\, h O^{-1} h]_{11}$, with
  $o_1$, $p_1$ of the centre.  This is not $(O O^{-1})_{11}$, because $O$ carries $o$ and $p$ of the intermediate sites.  Only
  the trace over all sites gives the number of states.  The present sum-rule print ($M_0 = 1$ at 6 Ry) needs a new target.
- **Representation.**  The generalised problem should not depend on $\alpha$ apart from the $S$ cut-off.  This is to be
  tested (gate 4), not assumed.

## 9. Open items
1. **The k-space oracle's own bookkeeping: closed (2026-11-04).**  For B2 after one iteration, $\sum_l (q_0 + p^\gamma q_2)$
   from the printed table gave qnet 0.1086 against the printed 0.1193.  Cause: after the SCF loop `asa_scf` recomputes the
   printed parameters, $p$ included, at the mixed new potential and the updated $E_\nu$ (src/asr_asa.f90:327-337).  The
   printed $q_0, q_1, q_2$ are the moments of the last iteration.  After one iteration they do not match.  In a converged
   run they do: the 128-site cell `trial/b45/ks/a0_s1` gives Cu $\sum_l (q_0 + p^\gamma q_2)$ = 11.08176, so qnet 0.08176
   (printed 0.08175), and Zn -0.08176 (printed -0.08175).  So the oracle is exact.  Compare qnet, not table sums, for
   one-iteration runs.
2. **Band energy.**  `asa_scf` uses $q_1 + E_\nu q_0$ per channel.  At convergence $q_1 = 0$ ($E_\nu$ at the band centre).  That treats $q_1$ as an energy moment.  In `ed3` the
   band energy should be the occupied sum of $E_\alpha z_\alpha^T O z_\alpha$, a contour integral of
   $E\,[G_{11} + 2 o (h G)_{11} + p (h G h)_{11}]$.  Decided (user, 2026-11-05): use the k-space form, $q_1 + E_\nu q_0$ per
   channel, so that the total energy matches the oracle.
3. **Tail of the cross term: closed for the charge.**  `cf_species` already gives the projection with the tail.  The cross term is
   weak near the axis, so the production chains use $|a\rangle, |v\rangle$ (section 6) and the charge does not need it.  Only
   $q_1$ does.  The two chains give it twice, and each SCF iteration prints the largest difference.
4. **E_F search.**  It needs $N(E)$ from the sphere charges at each trial $E_F$.  All three blocks are needed there too.
5. **Cost.**  Two chains per energy.  Measure against `asr_edc`.
6. **Storage of $w$: closed.**  $|w\rangle = S' |1 L\rangle$ covers the whole $S'$ range (bcc with `rcl` 1.2 alat: 8 + 6
   sites), not only the first shell.  Build it as $h'|u\rangle - g_c |u\rangle$ with the operator itself.  Then it is in the
   normalised orbit coordinates of the reduced operator.  Gate 1 B2 confirms it: reduced and unreduced agree for all four
   blocks to 6e-11 relative, and $|w|$ agrees to 2e-16.

## 10. Gates
Gate 1 results (2026-11-04).  Part A (`test_asr_ed3`, binary ring, `asr_augorb`): complex symmetry 8e-17; the three blocks
against the exact average of $K = E O - H$ 8e-15, 2e-15, 3e-16; two chains against the dense solve 2e-15; the
$q \to -q$ regression against $\langle [E - H^{(2)}]^{-1}_{11} \rangle$ 2e-15; canonical form 9e-16 (with $V$ starting at 1
instead of $\tfrac12$ it fails at 8e-2); single-site dressing exact at $p^\gamma = 0$ (3e-16), and 2.5 to 5 % off at
$p^\gamma \varepsilon^2 \approx 0.055$.  Part B (`test_asr_ed3_cls`, `asr_augcls`): ring of 5 with four species, ordered,
binary and ternary classes and two conditional shell classes (zero weights included): all checks at 1e-13 or better;
fcc ternary with the conditional shell class, D4h-reduced against unreduced: 6e-11 relative.  Zero-weight local states
have empty rows and columns (`class_mt`); a dense check must put 1 on their diagonal.

1. **Unit test (laptop).**  A small ring or cluster with random species, unequal $E_\nu$ and SRO
   ($\alpha$ = -0.2, 0, +0.2).  Brute-force configuration average of $G_{11}$, $(h G)_{11}$, $(h G h)_{11}$ from dense
   $K = E O - H$ against the augmented-space chains run to full length.  Target $10^{-12}$.  The same test measures the
   error of the single-site variant (section 7).  A second check: the canonical form equals the factored matrix-vector
   product to machine precision.  A third check: with $u = 0$, $p = 0$ and $q \to -q$, $K'$ is the `asr_edc` operator,
   so the chain must reproduce `asr_edc` to machine precision.
2. **Ordered B2 \ce{CuZn} at free atoms** (desktop).  Against k-space full ASA: $q_\mathrm{Cu}$ 0.1193, $E_F$ -0.1623 Ry.
   The agreement should be as good as $H^{(2)}$ against k-space $H^{(2)}$ (0.010 e in qnet, 0.007 e per channel).  Match
   the $S$ range.  Compare qnet, not table sums (open item 1).
3. **128-site SQS at free atoms** (desktop, `trial/b45/hop/ks_*`).  $q_\mathrm{Cu}(0)$ near 0.146.  The hopping charge at
   $\alpha = -0.148$ near the exact first-shell model value, -0.0072 +- 0.0006 e.  The correlated value -0.0050 e is not the
   target, because this route keeps the first-shell model.  Also the Zn d band position against the k-space DOS.
4. **Representation and SCF.**  Shift $\alpha$ by 20 %: the charges should not move beyond the $S$ cut-off effect.  A full
   \ce{CuZn} SCF should converge as well as the current route.  Report the cost per iteration.

## 10a. Gate 2 result (2026-11-06): not passed as set
Binary `ab50499`, desktop `~/asr_ed3/trial/gate2`, `gate2d`, `gate2e`, `gate2f`.  B2 \ce{CuZn} at free atoms, one iteration.
The cleanest oracle is k-space full ASA at nk 32, Gaussian 1 mRy: qnet 0.1193, E_F -0.1617 Ry, E -3446.74519 Ry.

| run | qnet Cu | E_F (Ry) | E + 3446 (Ry) |
|---|---|---|---|
| `ed3`, supercell 8, 36 steps, 32 contour points (gate input) | 0.1217 | -0.1833 | -0.74669 |
| `ed3`, supercell 16, 96 steps, 32 points | 0.1310 | -0.1608 | -0.74367 |
| `ed3`, supercell 16, 96 steps, 64 points | 0.1292 | -0.1623 | -0.74457 |
| `ed3`, supercell 16, 96 steps, 128 points | 0.1295 | -0.1622 | -0.74455 |
| k-space, nk 32, Gaussian 1 mRy | 0.1193 | -0.1617 | -0.74519 |
| k-space, nk 32, tetrahedra | 0.1183 | -0.1610 | -0.74458 |

Against the limits: E_F passes at the converged settings (0.5 mRy), and the total energy agrees to 0.6 mRy.  The charge does not pass: qnet is 0.0102 high (limit 0.010), and the channels at 128 points are Cu s -0.0042, p -0.0063, d +0.0206, Zn s -0.0042, p -0.0102, d +0.0042 (limit 0.007).

What the diagnostics fixed:
1. Representation is not the cause.  k-space full ASA at standard alpha against sigma 0.7 ($\alpha_d$ differs by 57 %): E_F 0.02 mRy, qnet 0.0001 e, E 0.04 mRy.  The same change in k-space $H^{(2)}$: E_F 1.3 mRy, qnet 0.005 e, E 3.8 mRy.
2. The E_F gap of the gate input is finite size.  An ordered run on a periodic graph has a discrete spectrum, like a Gamma-centred k-mesh of the supercell.  E_F moved from -0.1833 to -0.1608 with supercell 8 to 16.  k-space at nk 8 is off by 11 mRy in the same way.  `edc` closes from about 36 mRy to about 3 mRy.
3. 32 contour points are too few for ordered cells.  With 128 points the weighted sum rule at 6 Ry is 1.0000 (0.9991 at 64, 0.953 at 32).  This confirms the electron count of section 8 in the production code.
4. Steps, tail, depth tolerance and Fermi smearing do not explain the charge: 96 and 160 steps agree to 1e-4 e at supercell 16.

Open: the residual transfer of about 0.01 e from Zn to Cu.  `edc` has the same sign, +0.006 e against k-space $H^{(2)}$ at standard alpha.  So the origin may be shared by both routes, but this is not shown.  The unit tests exclude the operator algebra.  Untested suspect: $S$ range and inversion cluster at the converged settings (`rcl` 2.0 was only run at supercell 8 with 32 points).

**Update 2026-11-07: the residual is the kept $S$ range, and gate 2 passes at converged settings** (`trial/gate2g`; supercell 16,
96 steps, 64 contour points).  `ed3` at `rcl` 2.0: qnet 0.1186, E_F -0.1621 Ry, E -3446.74575 Ry, channels Cu -0.0007, +0.0002,
-0.0002, Zn -0.0010, +0.0025, -0.0009.  All meet the limits.  A larger inversion radius at `rcl` 1.2 changes nothing (qnet
0.1293), so the kept range matters and the inversion cluster does not.  `edc` at `rcl` 2.0: qnet 0.1221 against 0.1219 in k-space
$H^{(2)}$ at standard alpha.  So the residual was shared by both routes.  `rcl` 1.2 is not converged for this charge: it
moves about 0.01 e from Zn to Cu.  The cost at `rcl` 2.0 is 7.4 times that of `rcl` 1.2 for `ed3` and 4.4 times for `edc`.

## 11. Provenance
The first draft was written by another model from a prompt of this thread.  A second draft by the same model took over the
corrections below and added the closed form of the second-order limit and the $q \to -q$ regression check; both are checked
here.  The $H$ of section 2 in the prompt had a doubled linear $h$; it is fixed above.  Checked here: Results 1 to 4, the canonical
form, the matrix-vector product, Result 5, Result 6 and the limits.  Corrected: the factor $\tfrac12$ in $V$ (the draft had
1); the second-order limit (the draft forced a sign to match `asr_edc`); the one-chain claim (not shown for $p^\gamma \ne 0$);
the size of $p^\gamma$ for d bands (the draft assumed 0.05 or less); the pole of the dressing factor (the draft placed it
5 Ry away); the per-orbital sum rule (the draft claimed 1).
