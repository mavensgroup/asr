---
title: Energy-dependent ASR
---

# Energy-dependent ASR

The recursion needs an operator. For the TB-LMTO-ASA Hamiltonian
([Hamiltonians](../hamiltonian/index.html)) the operator depends on energy, and
the code offers two production forms and keeps older forms as modules for
tests and cross-checks. Only the production forms can be selected in the
input file, with `[asr] route`.

| `route` | Operator | Recursions | Module |
|---|---|---|---|
| `edc` (default) | AMscf form \(H(E)=E-T(E-E_\nu-h)T-ThohT\), the orthogonal \(H^{(2)}\) | One chain per (central species, start orbital, contour energy); start site frozen | [[asr_edc]], [[asr_ed]], [[asr_edn]] |
| `ed3` | Full generalised ASA problem \(K(E)=EO-H\), including the \(hph\) term | As `edc`, with two start vectors per chain | [[asr_ed3]] |

| Reference form (modules only) | Operator | Module |
|---|---|---|
| Energy-independent | \(H^{(2)}=E_\nu+h-hoh\), a fixed operator; stored Lanczos basis, \(G\) at any energy | [[asr_augsym]], [[asr_alloy_scf]] |
| \(P(z)-S\) | \(A(z)=P^\alpha(z)-S^\alpha\), complex-symmetric chain per energy | [[asr_endep]] |

!!! tip "Which route should I use?"
    Use `edc`, the production route. `ed3` adds the \(hph\) term that `edc` omits
    and agrees better with the k-space full-ASA reference where that term
    matters (derivation [below](#ed3-the-full-generalised-asa-problem)). It costs more per
    iteration. The energy-independent and \(P(z)-S\) forms are exact
    cross-checks of `edc` at convergence but have no input switch.

## Energy-independent (reference form)

\(H^{(2)}\) is a fixed operator. The augmented-space operator is built once per
SCF iteration, four recursions per spin run (orbitals s, p\(_z\), xy, 3z²−r²),
the Lanczos bases are stored, and the species Green functions are formed from
\(\langle n_RG\rangle\) via the stored basis ([[asr_partial]]). Moments \(M_k\)
on a contour give \(E_F\), the charge and the density

\[
\rho_l=M_0\,\phi^2+2M_1\,\phi\dot\phi+M_2\,(\dot\phi^2+\phi\ddot\phi) .
\]

## \(P(z)-S\) (reference form)

At each complex energy \(z\) the recursion runs on

\[
A(z)=P^\alpha(z)-S^\alpha,\qquad
P^\alpha(z)=\frac{z-C}{\Delta+(\gamma-\alpha)(z-C)} ,
\]

with \((C,\Delta,\gamma)\) of the species on each site. The physical
(orthonormal, gamma-representation) Green function is

\[
G=\lambda+\mu\,g^\alpha\mu,\quad
\mu=\frac{\sqrt\Delta}{\Delta+\delta(z-C)},\quad
\lambda=\frac{\delta}{\Delta+\delta(z-C)},\quad
\delta=\gamma-\alpha .
\]

\(A(z)\) is complex-symmetric, so the chain uses the complex-symmetric Lanczos
process ([[asr_edn]] for \(n\) species). The exact potential function
from the radial equation at each \(z\) ([[asr_pfun]]) can replace the two-parameter \(P\) in this form; then \(\lambda=C'/C\) is used so that the poles at the real
zeros of \(C\) cancel and the function is analytic elsewhere (the textbook
\(-P''/(2P')\) has spurious poles in the upper half-plane at the complex zeros
of \(P'\) of an exact, meromorphic \(P\)).

Because \(n_R\) is a projector, \(n_R\mu_R^2=\mu_A^2n_R\), the species-resolved
functions need only the average \(g_{RR}\) and the transition element
\(\langle n_Rg_{RR}\rangle\) from one chain per energy.

## `edc`: the AMscf form

For a **central species** \(c\) (the start site frozen to \(c\)), one recursion
per species and orbital:

\[
H(E)=E-T(E-E_\nu-h)T-ThohT,\qquad
T_R=\big(\Delta_c/\Delta_R\big)^{1/2},\qquad
G_c(E)=\langle R_0|(E-H(E))^{-1}|R_0\rangle .
\]

Since \(Td=\sqrt{\Delta_c}\) is constant, the structure constants appear only as the
**non-random** \(S'=\Delta_c^{1/2}S\Delta_c^{1/2}\), and all disorder sits in
site-diagonal terms:

\[
H(E)=\big[E-\Delta_c(E-E_{\nu,R})/\Delta_R\big]+h'-h'qh' ,\qquad
h'=g+S',\ \ g_R=\tfrac{\Delta_c}{\Delta_R}(C_R-E_{\nu,R}),\ \ q_R=\tfrac{o_R\Delta_R}{\Delta_c}.
\]

The identity \(E-H(E)=T\,(E-H^{(2)})\,T\) is exact for every configuration. Proof:
\(T h T=h'\), \(T(E-E_\nu)T=\Delta_c(E-E_\nu)/\Delta_R\) and
\(T\,h\,o\,h\,T=h'\,(T^{-1}oT^{-1})\,h'=h'qh'\). Here \(T_R=d_c/d_R\) with the signed
\(d=\sqrt\Delta\), so \(T d\) is the same on every site for any sign convention. The code
uses \(|d_c|\). This differs by a diagonal sign matrix per orbital and does not change a
diagonal element. \(T=1\) on the frozen centre, so
\(G_c(E)=\langle R_0|(E-H^{(2)})^{-1}|R_0\rangle\) averaged over all configurations with
species \(c\) on \(R_0\). This is the species-resolved Green function. The
energy-independent route with a frozen centre gives the same function. The two routes are
exact cross-checks of one another at convergence. At a finite number of steps they are
different approximants.

## `ed3`: the full generalised ASA problem

`edc` works with the orthogonal \(H^{(2)}\). It has no \(hph\) term. The k-space reference
solves the full generalised ASA problem. `ed3` brings the recursion to that problem
([[asr_ed3]]). It keeps the frozen centre and the scaling of `edc`.

### The secular matrix

All potential parameters are diagonal in site and orbital. The screening \(\alpha\) is the
same for every species, so \(S=S^\alpha\) does not depend on the configuration. Write
\(c=C^\alpha-E_\nu\), \(D=\Delta^\alpha\), \(d=D^{1/2}\), \(o=o^\alpha\),
\(p=p^\alpha=p^\gamma+o^2\), \(h=c+dSd\) and \(\varepsilon=E-E_\nu\). The overlap and the
Hamiltonian of the ASA problem are

\[
O=1+oh+ho+hph,\qquad H=E_\nu+h+E_\nu oh+hoE_\nu+h(o+pE_\nu)h .
\]

Their secular matrix is

\[
K(E)=EO-H=\varepsilon-h+\varepsilon oh+ho\varepsilon-hoh+hp\varepsilon h .
\]

It factorises as \(K=(1+ho)\,[\varepsilon(1+oh)-h]+h\,p^\gamma\varepsilon\,h\). The overlap is
\(O=(1+ho)(1+oh)+hp^\gamma h\). Both terms are positive semidefinite, so \(O\) is positive
definite when \(1+oh\) is invertible. The eigenvalues of the pencil are then real, and the
contour integral of \(G=K^{-1}\) is valid.

### Central-species scaling

Fix the centre to species \(c\) and scale with \(T_R=(D_c/D_R)^{1/2}\), as in `edc`. Then
\(h'=ThT=g+S'\) with the random diagonal \(g_R=(D_c/D_R)\,c_R\) and the non-random
\(S'=D_c^{1/2}SD_c^{1/2}\). The scaled secular matrix is

\[
K'(E)=TKT=\tilde\varepsilon-h'+uh'+h'u-h'\tilde Qh' ,
\]

\[
\tilde\varepsilon_R=\frac{D_c}{D_R}\,\varepsilon_R,\qquad u_R=o_R\varepsilon_R,\qquad
\tilde Q_R=\frac{D_R}{D_c}\,(o_R-p_R\varepsilon_R) .
\]

The proof uses only that \(T\) is diagonal: \(T(hXh)T=h'(T^{-1}XT^{-1})h'\) for any diagonal
\(X\), and \(T^{-2}=D_R/D_c\). The scaled matrix has four properties.

* Every random factor is site diagonal. It sits around the non-random \(S'\). This is the
  structure the augmented space needs. It holds for any number of species and sublattices,
  for each spin, and for an \(E_\nu\) per orbital.
* Every coefficient is linear in \(E\). So \(K'\) has no poles. The pole problem of the
  \(P(z)-S\) form cannot occur.
* \(K'\) is symmetric. At complex \(E\) it is complex symmetric. Its diagonal factors depend
  on species and \(l\) only, so the point-group reduction and the SRO shell classes apply
  unchanged.
* \(T=1\) on the centre. So \(\langle 1|K'^{-1}|1\rangle=G_{11}\).

The code applies \(K'\) in its canonical form
\(K'=D^\text{diag}-(VS'+S'V)-S'\tilde QS'\). Here
\(V_s=(\tfrac12+c_so_s)-(o_s+c_sp_s)\,\varepsilon_s\). The single \(-S'\) term splits into two
halves, so \(V\) starts with \(\tfrac12\). One Lanczos step needs two applications of \(S'\),
as in `edc`.

### Moments from two chains

A sphere charge needs three central blocks, \(G_{11}\), \((hG)_{11}\) and \((hGh)_{11}\). The
energy moments of one Green function are not enough. Take \(|u\rangle=|1L\rangle\),
\(|v\rangle=h'|u\rangle\) and \(|a\rangle=|u\rangle+o_c|v\rangle\). Since \(T=1\) on the
centre, the moments in the gamma representation are

\[
q_0=\langle a|K'^{-1}|a\rangle,\qquad q_1=\langle v|K'^{-1}|a\rangle,\qquad
q_2=\langle v|K'^{-1}|v\rangle .
\]

\(|v\rangle=g_c|u\rangle+S'|u\rangle\) does not depend on the configuration. So both vectors lie
in the reference sector of the augmented space. They also lie in the same reduced basis. One
complex-symmetric Lanczos chain starts from each. \(q_0\) and \(q_2\) are diagonal elements.
Their error is of second order in the Lanczos residual. \(q_1\) comes from projecting each chain
onto the other start vector. That error is of first order. The two chains give \(q_1\) twice.
The code uses the mean and prints the largest difference.

The sphere charge is \(q_0+p^\gamma q_2\). It needs no cross term. It is a positive sum of two
Herglotz functions. The density is \(q_0\phi^2+2q_1\phi\dot\phi+q_2\dot\phi^2\). The band energy
of a channel is \(q_1+E_\nu\,(q_0+p^\gamma q_2)\). This follows from \(E_\alpha z^TOz=z^THz\) for
every eigenvector. The form \(q_1+E_\nu q_0\) holds only for the energy moments of `edc`. In
the generalised problem it would lose \(E_\nu p^\gamma q_2\) per channel.

One chain is not enough. With chains from \(|u\rangle\) and \(S'|u\rangle\) instead, the charge
needs the cross term \((hG)_{11}\), which carries the first-order error. On terminated chains
that choice also lost the positivity of the \(2\times 2\) block of \(-\mathrm{Im}\,G\) at some
energies. The chains from \(|a\rangle\) and \(|v\rangle\) kept it at every test energy. When
\(p^\gamma=0\), one chain would suffice, because \((hz)_1\) then follows from \(z_1\) on every
eigenvector. For d bands \(p^\gamma\) is not small, so the code does not use this shortcut.

### Limits

* **Ordered limit.** With one species per sublattice the augmented space is the reference
  sector alone. Then \(K\) is the k-space pencil term by term. The two differ only by the
  range of \(S\) and by the chain length.
* **Second-order limit.** At \(p^\gamma=0\),
  \((1+oh)\,G\,(1+ho)=[\varepsilon-h(1+oh)^{-1}]^{-1}\) exactly. Expanding
  \(h(1+oh)^{-1}=h-hoh+\dots\) gives \([E-H^{(2)}]^{-1}\). So \(H^{(2)}\) is the second-order
  truncation, and `ed3` and `edc` agree in the physical moments to second order.
* **Electron count.** The sum of the sphere charges over all sites is the number of occupied
  states. The weight of a single orbital is not 1, because \(O\) carries \(o\) and \(p\) of the
  neighbours. Only the trace over all sites counts states.
* **Representation.** The generalised problem does not depend on \(\alpha\) apart from the
  range of \(S\). `ed3` has about half the representation dependence of `edc`.

## The reduced operator and its symmetric subspace

The point-group reduction stores one representative per orbit of nodes. A reduced vector is
valid only in the *symmetric subspace* \(V\): for every operation \(g\) of the stabiliser
of a representative node \(r\),

\[
\psi(r)=\chi(g)\,D(g)^{T}\psi(r).
\]

The reduced operator maps \(V\) into \(V\) and is symmetric on \(V\). On the rest of the
coordinate space it is a different matrix, and that matrix is not symmetric. A Lanczos chain
starts in \(V\). Rounding moves it out of \(V\) by \(10^{-16}\) per step. The recurrence
then amplifies this part. In an ordered B2 cell the part outside \(V\) is \(10^{-11}\) of
the vector after 24 steps, \(10^{-5}\) after 36 steps and several percent after 48 steps.
The `ed3` chains then miss the exact result by \(10^{-5}\) to \(10^{-4}\) however long they
run. So every chain applies the projector

\[
(P\psi)(r)=\frac1{|S_r|}\sum_{g\in S_r}\chi(g)\,D(g)^{T}\psi(r)
\]

to each new Lanczos vector ([[asr_augcls]] `augcls_project_c`, [[asr_augblk]] `blk_project`).
\(P\) is a real symmetric projector, and the operator maps \(V\) into \(V\). So the chain with
\(P\) is the Lanczos chain of the operator restricted to \(V\). With it the reduced
chains equal a dense solve on the unreduced operator to \(10^{-14}\) (`test_augcls_project`).
In random alloys at 24 to 48 steps the projection changes nothing visible.

## Causality of a truncated chain

The exact \(G_c(z)\) is a Herglotz function: \(-\mathrm{Im}\,G_c>0\) for
\(\mathrm{Im}\,z>0\). An energy-independent Lanczos chain keeps this property at every
length, because its fraction has real coefficients. The energy-dependent chain does not. Write
\(z-H(z)=zM-A\) with \(M=T^2\) and \(A=T H^{(2)} T\), both real. The chain of \(m\) steps
gives the bilinear Galerkin approximant on the Krylov space of \(zM-A\). That space depends on
\(z\). The approximant is a rational function of \(z\) with real coefficients. Its poles
are the zeros of a real polynomial, and they need not be real. A pole in the upper half plane
makes \(-\mathrm{Im}\,G\) negative near it. Such points occur with and without a tail. A
scan of an fcc ternary test operator (three central species, four start orbitals, 801 energies,
24 to 48 steps) found them at isolated energies for \(\mathrm{Im}\,z\) = 0.001 to 0.02. With
the `hull` tail the worst case was a DOS of \(-2.1\) states/Ry at one energy
(\(\mathrm{Im}\,z\) = 0.005, 48 steps), with 1 % of the weight negative. The tail adds a second
mechanism: its band edges are complex, and its branch cut can lie above \(z\).

Four statements are proven.

1. The chain is causal for one and two steps. A pole off the real axis can occur from three
   steps on. On the real axis a pole of negative weight belongs to a Ritz value that is not
   converged. Converged Ritz values give poles of positive weight.
2. The residues of the chain without a tail sum to one. A non-causal pole moves weight. It
   does not create or destroy it.
3. No tail can repair the chain in general. Let \(g_{ij}\) be the elements of
   \((z-T_m(z))^{-1}\). Every tail value with \(\mathrm{Im}\,t\le 0\) gives
   \(-\mathrm{Im}\,G>0\) if and only if the \(2\times 2\) matrix
   \(-\mathrm{Im}\,[g_{11}, g_{1m}; g_{1m}, g_{mm}]\) is positive definite. A tail with real
   band edges has \(\mathrm{Im}\,t\le 0\), but the matrix is often not positive definite.
4. The contour sum does not need an analytic approximant. It needs the chain to be accurate at
   the nodes. A pole pair between the contour line and the real axis enters the electron
   count with \(2\,\mathrm{Re}\,[r_p f(z_p)]\), where \(f\) is the Fermi function at the complex
   pole. Far below \(E_F\) this is the ordinary weight of the pair. Within about \(10\,kT\) of
   \(E_F\) the count is then not monotone in \(E_F\), and the SCF can stall.

Both routes print the number of contour node values that are not causal, and `asr.out`
carries a warning when the number is not zero in the last iteration. Treat a DOS with
negative values as not converged in the number of steps.

## Tails, depth and cost

| Setting | Effect |
|---|---|
| `ed_terminator = "hull" \| "bp" \| "avg2"` | Tail of each energy's chain; default `hull`, which is the safe choice ([Recursion and terminators](recursion.html)) |
| `depth_tol` | A chain stops early (from level 8, every 4 levels) when \(G\) changes by less than this (relative); `levels` is the maximum |
| `contour_points` | Number of energies (and therefore chains) |
| `fermi_kt` | Smeared occupations; makes \(N(E_F)\) smooth in the SCF |

Statistics on how many chains and levels were used are kept in
[[asr_chain]] (`ed_nchain`, `ed_nlev`).

## An exact reference: k-space

`method = "kspace"` ([Input: k-space](../input/kspace.html)) solves the same
Hamiltonian for an *ordered* structure by the Brillouin-zone integral. It is the
exact reference of the recursion at identical potential parameters.
