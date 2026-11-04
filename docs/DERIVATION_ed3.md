# Full-ASA energy-dependent recursion (route `ed3`): derivation

Branch `ed3`, written 2026-11-04.  Status: the operator algebra is proven.  The moment extraction is specified.  Nothing is
implemented.  Open items are listed in section 9.  Main stays on the production route (`asr_edc`) until the gates of
section 10 pass.

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

The augmented space reproduces $H^{(2)}$ to about 1 %.  The error is $H^{(2)}$ itself.  It doubles the SRO hopping response.
So the recursion needs the full generalised problem.

## 2. Notation (as in the code)
All potential parameters are diagonal in site $R$ and orbital $L$.  The screening $\alpha$ is the standard TB-LMTO set.  It is
the same for every species.  So the screened structure constants $S = S^\alpha$ do not depend on the configuration.
- $E_\nu$, $c = C^\alpha - E_\nu$, $D = \Delta^\alpha$, $d = D^{1/2}$, $o = o^\alpha$ (`o_alpha`), $p = p^\alpha$ (`to_alpha`).
- `to_alpha` (src/asr_lmto_ham.f90:126) gives $p^\alpha = p^\gamma - (o^\gamma)^2 + (o^\alpha)^2$.  From the gamma
  representation $o^\gamma = 0$.  So $p = p^\gamma + o^2$.  Here $p^\gamma = \langle \dot\phi | \dot\phi \rangle > 0$.
- $h = c + d S d$.  It is real symmetric.
- $\varepsilon(E) = E - E_\nu$, diagonal.
- Generalised problem of `bands_core` (src/asr_asa.f90, `zhegv` branch):
  $O = 1 + o h + h o + h p h$ and $H = E_\nu + (1 + E_\nu o) h + h (1 + o E_\nu) + h (o + p E_\nu) h$.

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
1. **The k-space oracle's own bookkeeping.**  For B2 at free atoms, $\sum_l (q_0 + p^\gamma q_2)$ from the printed table gives
   qnet 0.1086.  The run prints 0.1193.  The gap of 0.011 e must be explained before k-space serves as a 0.01 e oracle.
   In `second_order` mode the printed occupations add up to qnet exactly.  Suspects: the small-component terms in the
   $\dot\phi^2$ norm, or a different $p$ in the print.
2. **Band energy.**  `asa_scf` uses $q_1 + E_\nu q_0$ per channel.  That treats $q_1$ as an energy moment.  In `ed3` the
   band energy should be the occupied sum of $E_\alpha z_\alpha^T O z_\alpha$, a contour integral of
   $E\,[G_{11} + 2 o (h G)_{11} + p (h G h)_{11}]$.  Decide which form to use.  The total energy must match the k-space
   form for the oracle.
3. **Tail of the cross term.**  $G_{wu}$ uses the full column $[J^{-1}]_{m1}$.  Check that the `hull` tail treats it as well
   as $G_{uu}$.  The two chains give $G_{wu}$ twice; their difference measures this.
4. **E_F search.**  It needs $N(E)$ from the sphere charges at each trial $E_F$.  All three blocks are needed there too.
5. **Cost.**  Two chains per energy.  Measure against `asr_edc`.

## 10. Gates
1. **Unit test (laptop).**  A small ring or cluster with random species, unequal $E_\nu$ and SRO
   ($\alpha$ = -0.2, 0, +0.2).  Brute-force configuration average of $G_{11}$, $(h G)_{11}$, $(h G h)_{11}$ from dense
   $K = E O - H$ against the augmented-space chains run to full length.  Target $10^{-12}$.  The same test measures the
   error of the single-site variant (section 7).  A second check: the canonical form equals the factored matrix-vector
   product to machine precision.
2. **Ordered B2 \ce{CuZn} at free atoms** (desktop).  Against k-space full ASA: $q_\mathrm{Cu}$ 0.1193, $E_F$ -0.1623 Ry.
   The agreement should be as good as $H^{(2)}$ against k-space $H^{(2)}$ (0.007 e per channel).  Match the $S$ range.
3. **128-site SQS at free atoms** (desktop, `trial/b45/hop/ks_*`).  $q_\mathrm{Cu}(0)$ near 0.146 and the hopping charge near
   -0.0050 e at $\alpha = -0.148$.  Also the Zn d band position against the k-space DOS.
4. **Representation and SCF.**  Shift $\alpha$ by 20 %: the charges should not move beyond the $S$ cut-off effect.  A full
   \ce{CuZn} SCF should converge as well as the current route.  Report the cost per iteration.

## 11. Provenance
The first draft was written by another model from a prompt of this thread.  Checked here: Results 1 to 4, the canonical
form, the matrix-vector product, Result 5, Result 6 and the limits.  Corrected: the factor $\tfrac12$ in $V$ (the draft had
1); the second-order limit (the draft forced a sign to match `asr_edc`); the one-chain claim (not shown for $p^\gamma \ne 0$);
the size of $p^\gamma$ for d bands (the draft assumed 0.05 or less); the pole of the dressing factor (the draft placed it
5 Ry away); the per-orbital sum rule (the draft claimed 1).
