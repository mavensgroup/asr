---
title: Recursion and terminators
---

# Recursion and terminators

## Lanczos / Haydock recursion

For a Hermitian operator \(H\) and a start vector \(u_1\), the Lanczos process
builds an orthonormal basis \(\{u_n\}\) in which \(H\) is tridiagonal:

\[
H|u_n\rangle=b_{n-1}|u_{n-1}\rangle+a_n|u_n\rangle+b_n|u_{n+1}\rangle .
\]

Vectors are normalised, so Haydock's \(\beta_n^2\) (unnormalised form) is
\(b_n^2\) here. Only a routine `apply(H, v)` is required; this is the abstract
interface of [[asr_operator]] and the reason **the same recursion code serves
tight-binding, LMTO and augmented-space operators** — the operator is matrix-free.
The implementation is [[asr_recursion]].

!!! note "Real symmetric today"
    The recursion is real-symmetric and scalar. Block and complex-Hermitian
    versions are needed for degenerate orbitals and non-collinear magnetism
    and are on the [roadmap](../appendix/roadmap.html).

## Continued fraction

The diagonal element of the resolvent is

\[
G_{11}(z)=\cfrac{1}{z-a_1-\cfrac{b_1^2}{z-a_2-\cfrac{b_2^2}{z-a_3-\dots}}} ,
\qquad \mathrm{LDOS}(E)=-\tfrac1\pi\,\mathrm{Im}\,G_{11}(E+i\eta) .
\]

[[asr_greens]] evaluates it from the bottom up with a *tail* \(t(z)\)
replacing the unknown part after level \(n\). If \(t\) is wrong, the
unresolved structure of the spectrum is missing: **depth, not the tail, is the
lever** — see below.

## Terminators

A terminator is a pair *(explicit extended coefficients, asymptotic tail)*
([[asr_terminator]]). Asymptotes \(a_\infty,b_\infty\):

| Asymptote | Rule |
|---|---|
| `last` | the last computed \((a_n,b_n)\) |
| `avg2` | average of the last two levels |
| Beer–Pettifor (`bp`) | \((A,B)\) such that the denominator of the terminated fraction vanishes exactly at both band edges \(A\pm2B\) — no loss of spectral weight; solved by fixed-point iteration from the `avg2` start |
| period 2 | levels repeat with period 2; the link \(b_n\) is exact |

The square-root tail of \((a_\infty,b_\infty)\) is

\[
t(z)=\frac{(z-a_\infty)-\sqrt{(z-a_\infty)^2-4b_\infty^2}}{2}.
\]

Joiners: *butt* (none) or *Luchini–Nex* (linear interpolation of the coefficients from
\((a_n,b_n)\) to \((a_\infty,b_\infty)\) over `ninterp` levels, then held constant).

!!! warning "Luchini–Nex is unverified"
    The implementation follows the published *description* via secondary
    sources. It has **not** been compared with the original paper. Treat results
    that depend on this joiner with care.

The energy-dependent route uses `hull`, `bp` and `avg2` tails ([[asr_endep]]), the
`ed_terminator` key. The default is `hull`.

| Tail | Rule on the chain at one contour energy |
|---|---|
| `hull` | The band edges are the extreme Ritz values of the chain itself. There is no fixed point and no end-site shift. |
| `bp` | The Beer–Pettifor fixed point. The last diagonal element is shifted to \(a_m - B\). |
| `avg2` | The average of the last two levels. |

!!! warning "Why `bp` is not the default"
    On an energy-dependent chain the coefficients are complex. The `bp` end-site
    shift can then split off a bound state. Its pole can lie very close to a
    contour node. The Green function of that chain stays continuous, but it
    becomes very steep in energy. In the Fe66Ni14Cr20 fcc run at levels 4 this
    made \(N(E_F)\) and the Cr charge jump. The SCF stalled for 72 iterations.
    With `hull` the same run converged in 6 iterations.

    The exact \(G_c(z)\) is the resolvent of a real symmetric operator, so it
    obeys the Nevanlinna bound \(|dG/dz| \le -\mathrm{Im}\,G/\mathrm{Im}\,z\).
    `test_terminators` checks this bound on two real Cr e_g chains. `hull`
    reaches 0.41 of the bound. `bp` exceeds it by a factor of about 1000. No
    general repair of `bp` was found (`docs/FINDINGS_LOG.md`, "Levels-4 SCF
    stall"). The `bp` tail is kept for comparison only.

    `bp` was first chosen because `hull` is too narrow on very short chains
    (4 to 6 steps). Production chains have 24 or 36 steps. There all tails agree
    within 0.005 muB.

## What the measurements say

Measured against exact k-space truth (max \(|\Delta\mathrm{LDOS}|\),
`app/term_bench`, `test_terminators`):

* sc, fcc: all single-band terminators agree within a few percent once depth ≥ 16;
  depth limits accuracy.
* bcc: nothing helps; ≈ 5·10⁻³ at depth 60 for every one-band tail (van Hove
  singularity at band centre, oscillating \(b_n^2\)). *Open problem.*
* gapped B2 compound: period-2 tail is 1–3 orders of magnitude better than any
  one-band tail (e.g. 1.9·10⁻² vs 2.7·10⁻⁵ at depth 60).
* Real d band (ordered Ag, LMTO): no terminator rescues shallow depth; ≈ 1 %
  LDOS error needs depth ≈ 40. The augmented space reaches only ~10–20 levels.

The consequence for the disordered code: use the highest `levels` the machine
allows, check LDOS convergence in `levels`, and treat features narrower than the
depth-limited resolution as unreliable. See the
[Convergence guide](../guides/convergence.html).

## Contour integration

Occupied moments are evaluated on the complex plane, where the continued fraction is
analytic and smooth:

\[
M_k=-\frac1\pi\,\mathrm{Im}\int_C (z-E_\nu)^k\,G(z)\,dz
=\int_{E_\text{bot}}^{E_F}(E-E_\nu)^k N(E)\,dE ,
\]

where \(C\) is a semicircle in the upper half-plane from \(E_\text{bot}\) to
\(E_F\), integrated by Gauss–Legendre in the angle ([[asr_contour]],
`contour_points`). For a finite temperature `fermi_kt` the contour is
extended by Matsubara poles (`contour_nodes_fermi`, `fermi_poles`).

## Reference methods

Independent of the continued-fraction machinery the code provides:

* [[asr_kspace]]: exact k-space Green function of ordered cubic NN bands,
* [[asr_kpm]]: Kernel Polynomial Method (Chebyshev, stochastic trace, Jackson
  kernel) used as ground truth for large disordered supercells,
* [[asr_lmto_realization]]: one random realisation of the LMTO alloy on a
  periodic supercell, for comparison with the augmented-space recursion.
