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
| `ed3` | Full generalised ASA problem \(K(E)=EO-H\), including the \(hph\) term | As `edc`, with two start vectors per chain | [[asr_ed3]] (derivation in `docs/DERIVATION_ed3.md`) |

| Reference form (modules only) | Operator | Module |
|---|---|---|
| Energy-independent | \(H^{(2)}=E_\nu+h-hoh\), a fixed operator; stored Lanczos basis, \(G\) at any energy | [[asr_augsym]], [[asr_alloy_scf]] |
| \(P(z)-S\) | \(A(z)=P^\alpha(z)-S^\alpha\), complex-symmetric chain per energy | [[asr_endep]] |

!!! tip "Which route should I use?"
    Use `edc`, the production route. `ed3` adds the \(hph\) term that `edc` omits
    and agrees better with the k-space full-ASA reference where that term
    matters (see `docs/DERIVATION_ed3.md`, section 1). It costs more per
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

Four statements are proven (`docs/DERIVATION_f4.md`).

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
[[asr_endep]] (`ed_nchain`, `ed_nlev`).

## An exact reference: k-space

`method = "kspace"` ([Input: k-space](../input/kspace.html)) solves the same
Hamiltonian for an *ordered* structure by the Brillouin-zone integral. It is the
exact reference of the recursion at identical potential parameters.
