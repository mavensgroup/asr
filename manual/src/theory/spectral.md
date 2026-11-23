---
title: Bloch spectral function
---

# Bloch spectral function

The configuration-averaged Bloch spectral function \(A(\mathbf k,E)\) is the band structure of
a disordered alloy. This page explains how the code obtains it from the recursion
([[asr_spectral]], [[asr_edc]] `edc_green_proj`). The [spectral guide](../guides/spectral.html)
explains the input and the output. The method works on route `edc`.

## Definition

Let \(\langle G_{0L,RL'}(z)\rangle\) be the configuration average of the physical Green function
between orbital \(L\) at the origin and orbital \(L'\) at site \(R\). For a translation-invariant
average of a real symmetric operator, the spectral matrix is a Fourier sum of the imaginary parts:

\[
A_{LL'}(\mathbf k,E)=\sum_{R}e^{i\mathbf k\cdot\mathbf R}
\left(-\frac1\pi\,\mathrm{Im}\,\langle G_{0L,RL'}(E+i\eta)\rangle\right) .
\]

The proof uses two facts. \(G\) is symmetric in its indices for every configuration. The average
is invariant under lattice translations. The code reports the trace over \(L\), per graph site and
per spin, summed over sublattices with their weights \(f_\kappa\).

The exact \(A\) has three properties. The first two hold at every \(\eta>0\).

* **Positivity.** The matrix \(A(\mathbf k,E)\) is Hermitian positive semidefinite for every
  configuration. An average keeps this.
* **DOS sum rule.** The Brillouin-zone average of \(A_{LL}\) is the on-site DOS of orbital
  \(L\). The phases leave only \(R=0\).
* **Ordered limit.** In an ordered crystal \(A\) becomes \(\sum_n w_{n\mathbf k}\,
  \delta(E-E_{n\mathbf k})\) as \(\eta\to0\). For a decoration with several sites in the cell this
  is the unfolded band structure.

The unfolding uses the Bloch vectors of the graph lattice. It needs every sublattice to be one
lattice of translations. Random alloys, B2 and L2\(_1\) qualify. L1\(_2\) does not, because its
face sites are related by rotations only.

## Off-site elements from the SCF chains

The SCF chains of `edc` run on the scaled operator of a frozen centre ([Energy-dependent
ASR](energy_dependent.html)). For every configuration
\(E-H^{(c)}(E)=T\,(E-H^{(2)})\,T\). Inverting gives

\[
G^{(2)}=T\,(E-H^{(c)}(E))^{-1}\,T ,\qquad T_R(s)=\left(\frac{D_{c,L'}}{D_{s,L'}}\right)^{1/2} .
\]

\(T=1\) on the centre. On the target site the factor is \(T_R(s)\) to the power \(+1\). So the
off-site average is

\[
\langle G^{(2)}_{0L,RL'}\rangle=\sum_c x_c\,
\langle 0L|\,(E-H^{(c)}(E))^{-1}\,|\widetilde{RL'}^{(c)}\rangle .
\]

The target vector is dressed with the scaling of site \(R\) in augmented space:

\[
|\widetilde{RL'}^{(c)}\rangle=\sum_k\tau_k\,|RL',\{R:k\}\rangle ,\qquad
\tau_k=\sum_s p_s(R)\,T_R(s)\,\phi_k(s) .
\]

Here \(\phi_k\) is the local basis of the site and \(p_s(R)\) its concentration. The target has at
most one flipped site. The centre alone is frozen, so the operator keeps the full site group. The
chain from \(|0L\rangle\) gives the vector \((E-H^{(c)})^{-1}|0L\rangle\) everywhere. So one chain
per centre species, start orbital and energy gives every target site at once. The spectral chains
use the operator reduced with D\(_{2h}\), the diagonal operations of the site group. There all
nine real orbitals are one dimensional. So every start orbital is allowed and no rotation between
orbitals is needed.

## Accuracy and checks

A diagonal element of a Lanczos chain has an error of second order in the residual. An off-site
element by projection has an error of first order. So off-site elements need longer chains than
the on-site one. `[spectral] steps` sets that length without repeating the SCF.

The truncated augmented space does not give \(G_{0L,RL'}=G_{0L',-RL}\) exactly. The two come
from different chains, one from each origin. The code uses their mean and prints the largest
difference as the reciprocity defect. It also prints the change when the chains are cut four
steps earlier.

## Short-range order

With first-shell SRO the centre fixes the conditional concentrations of its first shell,
\(X_s^{(c)}=x_s(1-\alpha_{cs})\). Every other site keeps \(x_s\). So a target in the first shell
is dressed with \(p_s(R)=X_s^{(c)}\), and a target further out with \(x_s\). The marginal
condition \(\sum_cx_c\alpha_{cs}=0\) makes every site carry \(x_s\) on average. The symmetrised
mean of the two origins makes \(A\) Hermitian at every \(\mathbf k\). At \(\alpha=0\) the
definition reduces to the one above.

The model has limits. It changes the first-shell bonding and the disorder broadening. It moves
band centres through the self-consistent charges. It cannot produce superlattice peaks or zone
folding at an ordering vector, because SRO beyond the first shell is zero. Long-range order needs
sublattices in the decoration.

## The window

The Fourier sum stops at `rmax`. A sharp stop rings. For an ordered metal at
\(\eta=0.01\) Ry the sharp sum dips to \(-30\,\%\) of its peak. The default window is the overlap
volume of two balls of diameter \(r_\text{max}\) at distance \(R\), divided by the volume of one
ball:

\[
w(R)=\left(1-\frac{R}{r_\text{max}}\right)^2\left(1+\frac{R}{2r_\text{max}}\right),\qquad
R<r_\text{max} .
\]

Its Fourier transform is the squared transform of a ball, so it is non-negative. Samples of such
a function on a lattice keep a non-negative Fourier series, by Poisson summation. So the windowed
\(A\) is the exact \(A\) convolved in \(\mathbf k\) with a non-negative kernel of width about
\(1/r_\text{max}\). It is non-negative wherever the real-space elements are exact. \(w(0)=1\), so
the DOS sum rule and the energy integral at each \(\mathbf k\) are unchanged. As \(r_\text{max}\)
grows, \(w\to1\) and \(A\) tends to the exact \(A\). The \(\mathbf k\) resolution is therefore
about \(1/r_\text{max}\).
