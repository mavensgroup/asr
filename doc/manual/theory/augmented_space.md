---
title: Augmented space
---

# Augmented space

## The problem

For a substitutional alloy each site \(R\) carries a random variable
\(n_R\in\{0,1\}\) (A or B) with \(\Pr(n_R=1)=x\). A Hamiltonian
\(H(\{n_R\})\) depends on the whole configuration, and a physical quantity is the
*configuration average* of the Green function

\[
\langle\langle G(z)\rangle\rangle=\Big\langle\!\Big\langle (z-H(\{n\}))^{-1}\Big\rangle\!\Big\rangle .
\]

The CPA replaces this by a single-site self-energy. The augmented-space
formalism (Mookerjee 1973) instead gives the **exact** average, which in
addition contains the effect of the environment (clustering, short-range
order, off-diagonal disorder).

## The theorem

Each random variable with a probability density admits an
*augmented space* \(\Phi_R\): a local Hilbert space in which \(n_R\) becomes a
**Hermitian operator** \(\mathcal M_R\) whose spectral density equals the
probability density of \(n_R\). The averaged Green function is then the
*diagonal element of an ordinary resolvent* in the enlarged space
\(\Psi=\mathcal H\otimes\prod_R\Phi_R\):

\[
\langle\langle G(z)\rangle\rangle=\langle R,\varnothing|\,(z-\widetilde H)^{-1}\,|R,\varnothing\rangle ,
\qquad
\widetilde H = H\big(\{n_R\to\mathcal M_R\}\big).
\]

Here \(|\varnothing\rangle\) is the *reference* configuration state
(all sites in the "average" state).

## Binary alloys

For a binary site, \(\Phi_R\) is two-dimensional with basis
\(|\varnothing\rangle\) (A/B average state) and \(|R\rangle\) (state with the site
*flipped*). The operator for \(n_R\) is

\[
\mathcal M_R = x\,\mathcal P_R^{\varnothing} + y\,\mathcal P_R^{\downarrow}
+ \sqrt{xy}\,\big(\mathcal T^{+}_R+\mathcal T^{-}_R\big),\qquad y=1-x ,
\]

where \(\mathcal P\) are projectors on the two local states and \(\mathcal T^\pm\)
flip the site; the flip amplitude is \(\sqrt{xy}\). Any site-diagonal random
parameter \(X_R=X_B+(X_A-X_B)\,n_R\) becomes \(X_B+(X_A-X_B)\,\mathcal M_R\)
(implemented as `X(n) = X_B + (X_A - X_B) n` in [[asr_lmto_alloy]]).

A **state** of the augmented space is a *node*: the site index plus the set
of sites currently flipped. This is the data structure in [[asr_augspace]],
which uses hashed states.

!!! note "Verified exactly"
    For finite rings and clusters the recursion equals the exact average over
    all \(2^N\) configurations to \(10^{-14}\) (`test_augspace`, `test_augorb`).

## More than two components (HEA)

For \(n_s\) species with probabilities \(x_s\) the local space is the
function space on the species set with inner product
\(\langle f,g\rangle=\sum_s x_s f(s)g(s)\). An orthonormal basis
\(\phi_0=1,\phi_1,\dots,\phi_{n_s-1}\) (Gram–Schmidt on the species
indicators) gives, for the indicator of species \(s\), the rank-one matrix

\[
M^{(s)}_{lk}=x_s\,\phi_l(s)\,\phi_k(s).
\]

A site-diagonal parameter \(X_R=\sum_s X_s[s_R=s]\) becomes
\(W_{lk}=\sum_s X_s M^{(s)}_{lk}\), and \(n_s=2\) reproduces the binary case
exactly. Implemented in [[asr_augn]]; for different species sets per site class
(sublattice disorder, ordered classes that never flip) see [[asr_augcls]].

## Why the node set is finite in practice

\(\widetilde H\) acts locally: a hop moves the current site, and a random
parameter can flip the site it acts on. Starting from \(|R_0,\varnothing\rangle\),
level \(n\) of the Lanczos recursion reaches only nodes within \(n\) hops,
with flips only on the visited sites. The node count grows quickly with depth and
is controlled by:

* `levels` — the depth,
* `rank` — the maximum number of simultaneously flipped sites; states beyond
  rank are folded into an average medium,
* `flip_radius` — only sites near the origin may flip,
* `supercell` — periodic graph size (so that no hopping leaves the graph),
* `symmetry` — the point group (next section).

## Symmetry reduction

The recursion is started on an orbital of the origin site, and
\(\widetilde H\) commutes with the point-group operations \(U_g\) that fix
the start site. If the start vector transforms as
\(U_g u_1=\chi(g)u_1\) with \(\chi(g)=D(g)[L_0,L_0]=\pm1\), it is enough to store
**canonical representatives** of each orbit of nodes. A non-representative node
\(f=g^{-1}r\) satisfies

\[
\psi(f)=\chi(g)\,D(g)^{T}\psi(r),
\]

and the vectors are kept in normalised coordinates
\(w=\sqrt{n_\text{node}}\,\psi\) (\(n_\text{node}\): orbit size) so that
the ordinary Lanczos recursion applies unchanged. The reduction is **exact**
(not an approximation) and reduces the node count by ≈ 40×.

| Module | Group | Scope |
|---|---|---|
| [[asr_symmetry]] | O\(_h\) (48) site permutations on periodic sc/bcc/fcc supercells | geometry |
| [[asr_augsym]] | D\(_{4h}\) (16) | binary, start orbital s, p\(_z\), t\(_{2g}\), e\(_g\) |
| [[asr_augcls]] | any supplied group fixing the start site and class | multi-class, multi-species |

The four start orbitals (s, p\(_z\), xy, 3z²−r²) are those with a
well-defined character under D\(_{4h}\), and the l-resolved DOS is assembled
from them with weights 1, 3, 3, 2 (the degeneracies of s, p, t\(_{2g}\), e\(_g\)).

## Species-resolved Green functions

The total average \(\langle\langle G_{RR}\rangle\rangle\) is a mix of A and B on
site \(R\). The *partial* functions follow from the occupation operator
\(n_R\to\mathcal M_R\):

\[
\langle\langle n_R\,G_{RL,RL}\rangle\rangle=\langle u_1|\mathcal M_R\,G|u_1\rangle ,
\qquad G_A=\frac{\langle\langle nG\rangle\rangle}{x},\qquad
G_B=\frac{\langle\langle G\rangle\rangle-xG_A}{1-x}.
\]

With the Lanczos basis \(\{u_k\}\) and \(v=\mathcal M_Ru_1\) this equals
\(\sum_kc_k[G]_{k1}\), \(c_k=\langle u_k|v\rangle\), and the first column
of the tridiagonal resolvent follows from the continued fraction tails
(\(g_j=1/(z-a_j-b_j^2g_{j+1})\), \([G]_{k1}=g_1\prod_{j<k}b_jg_{j+1}\)). This is
[[asr_partial]]. Because \(n_R\) is a projector the energy-dependent route
needs only the transition element \(\langle n_Rg_{RR}\rangle\) in addition
(see [Energy-dependent ASR](energy_dependent.html)).
