---
title: TB-LMTO-ASA
---

# TB-LMTO-ASA

This is the one Hamiltonian that `asr` implements at present. The method itself
([Theory](../theory/index.html)) does not depend on it; see
[Hamiltonians](index.html) for the interface.

The Hamiltonian is the second-order **tight-binding linear muffin-tin orbital**
Hamiltonian in the **atomic-sphere approximation** for an \(spd\) basis, in the
screened (\(\alpha\)) representation. *Every* parameter is generated inside
the code; nothing is read from tables, and nothing from other codes is used as
an input (Stuttgart's AMscf and Questaal are reference *code* only).

## From atom to potential parameters

For each species and each \(l\):

1. **Radial equation** ([[asr_radial]]). A scalar-relativistic
   Koelling–Harmon solver on a shifted logarithmic mesh
   \(r_i=b(e^{a(i-1)}-1)\), in Rydberg units (\(\hbar^2/2m=1\), \(e^2=2\),
   \(c=274.07\ldots\)):

   \[
   P'=MQ+\tfrac Pr,\qquad
   Q'=-\tfrac Qr+\Big[\tfrac{l(l+1)}{Mr^2}+V-E\Big]P,\qquad M=1+\tfrac{E-V}{c^2}.
   \]

   As \(c\to\infty\) this is the Schrödinger equation. The energy derivatives
   \(\dot\phi\) and \(\ddot\phi\) are obtained **exactly** from the variational
   equations, not by finite differences.
   `relativistic = false` selects the non-relativistic limit.

2. **Parameters** ([[asr_potpar]]) at the sphere radius \(S\), through Wronskians of
   the bare \(K\) and \(J\) functions with the normalised \(\phi\) and
   \(\dot\phi\) (\(\langle\phi|\dot\phi\rangle=0\)):

   \[
   C=E_\nu-\frac{W\{K,\phi\}}{W\{K,\dot\phi\}},\quad
   \sqrt\Delta=-\sqrt{\frac{(S/w)^{2l+1}}{2S}}\;\frac1{W\{K,\dot\phi\}},\quad
   \gamma=(S/w)^{2l+1}\frac{W\{J,\dot\phi\}}{W\{K,\dot\phi\}},\quad
   p=\langle\dot\phi|\dot\phi\rangle ,
   \]

   with \(w\) the average Wigner–Seitz radius.

## Structure constants

[[asr_strconst]] and [[asr_lmto_cell]] build the bare canonical (\(\kappa=0\))
structure constants \(S^0\) on a finite cluster around each atom and screen them
with the Andersen–Jepsen transformation

\[
S^\alpha=S^0\,(1-\alpha S^0)^{-1},\qquad
\alpha_l=\frac{\sigma^{2l+1}}{2(2l+1)}\ \text{(hard-core)},
\]

with `sigma` = 0.7 by default, or the *standard* TB-LMTO set
\(\alpha=(0.3485,\,0.05303,\,0.010714)\) for s, p, d. The real spherical
harmonics and the Gaunt coefficients ([[asr_harmonics]]) are computed by exact quadrature, so
nothing depends on tabulated Clebsch–Gordan data.

!!! note "Orbital order"
    Within an l-shell the order is \(m=-l..l\): s | p(y,z,x) | d(xy, yz, 3z²−r²,
    xz, x²−y²) — the Stuttgart/INFO order.

## The Hamiltonian

In the \(\alpha\) representation

\[
h=(C-E_\nu)+\sqrt\Delta\,S^\alpha\sqrt\Delta,\qquad
H^{(2)}=E_\nu+h-h\,o\,h,
\]

with \(o^\alpha=\dfrac{\gamma-\alpha}{(\gamma-\alpha)(E_\nu-C)-\Delta}\)
(as in the Stuttgart `makpph`). For an ordered crystal the \(9\times9\) block
\(H(\mathbf k)\) is built in [[asr_lmto_ham]] (one atom) or [[asr_lmto_cell]] (general
cell, per-species parameters and screening) and integrated over the Brillouin
zone with the linear tetrahedron method ([[asr_tetra]]).

For disordered alloys the structure constants must not depend on the local
composition, so **one common screening** is used for all species
([[asr_lmto_alloy]]); species enter only through the site-diagonal random
parameters \(c\), \(\sqrt\Delta\), \(o\), \(E_\nu\), which become operators
\(X_B+(X_A-X_B)\,\mathcal M_R\) in the augmented space. The matrix-free
operator application is [[asr_augorb]] and the point-group reduced version is
[[asr_augsym]]. [[asr_lmto_realization]] applies the same Hamiltonian to a *single
numerical realisation* and optionally implements the third-order ASA
Hamiltonian:

\[
H^{(3)}=H^{(2)}+h\,o\,h\,o\,h-\tfrac12\,(h\,p\,h\,h+h\,h\,p\,h),\quad p=p^\gamma .
\]

!!! warning "Conventions that were each a bug once"
    1. Stuttgart `ctrl` files store \(c\), \(\sqrt\Delta\), \(p\) in the **gamma
       representation**; they must be transformed to the structure-constant
       \(\alpha\) with `pptrns` before use. Skipping the transformation placed
       the p band at +3 Ry.
    2. \(S\) enters with sign **+**; the sign of \(\sqrt\Delta(p)\) is a pure gauge.
    3. A finite-cluster \(S^\alpha\) row is Hermitian only to truncation error, so the
       Bloch sum is Hermitianised.

## What is *not* included

* No combined correction and no downfolding; the basis is s, p, d.
* Collinear spin only; non-collinear magnetism needs the complex-Hermitian
  recursion first.
* Spherical (ASA) potentials; no full-potential terms.
