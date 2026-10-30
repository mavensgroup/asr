---
title: Energy-dependent ASR
---

# Three recursion routes

The `kind = "alloy"` and `"sublattice"` calculations select how the TB-LMTO
Hamiltonian is turned into a recursion with `[asr] hamiltonian`.

| `hamiltonian` | Operator | Recursions | Module |
|---|---|---|---|
| `energy-independent` | \(H^{(2)}=E_\nu+h-hoh\) (or \(H^{(3)}\)); a *fixed* operator | One per start orbital (s, p\(_z\), t\(_{2g}\), e\(_g\)) and spin; **stored Lanczos basis**, G at any energy | [[asr_augsym]], [[asr_alloy_scf]] |
| `energy-dependent` | \(A(z)=P^\alpha(z)-S^\alpha\), random diagonal \(P^\alpha(z)\) | One complex-symmetric chain **per complex energy** \(z\) | [[asr_endep]] |
| `asr_ed` | AMscf form \(H(E)=E-T(E-E_\nu-h)T-ThohT\) | One recursion per central species, start site frozen | [[asr_ed]], [[asr_edn]], [[asr_edc]] |

`kind = "sublattice"` always uses the `asr_ed` form (via [[asr_edc]] and
[[asr_class_scf]]).

!!! tip "Which one should I use?"
    The energy-**independent** route is cheapest (one set of coefficients gives the
    DOS at every energy) and is the natural first choice, but its
    \(H^{(2)}\) is only an approximation of the energy dependence of the
    potential function. The **energy-dependent** forms keep the exact \(z\)-dependence
    (linear \(P(E)=(E-C)/\Delta\) in `asr_ed`) and are the production method for
    real alloys. Use the two against each other as a cross-check.

## Energy-independent

\(H^{(2)}\) is expanded to a chosen `order` in \(h\). The augmented-space
operator is built once per SCF iteration, four recursions per spin run
(orbitals s, p\(_z\), xy, 3z²−r²), the Lanczos bases are stored, and the species
Green functions are formed from \(\langle n_RG\rangle\) via the stored basis
([[asr_partial]]). Moments \(M_k\) on a contour give \(E_F\), the charge and the density

\[
\rho_l=M_0\,\phi^2+2M_1\,\phi\dot\phi+M_2\,(\dot\phi^2+\phi\ddot\phi)
\]

(`density = "expand"`; with `"stuttgart"` the last term is
\(M_2(\dot\phi^2-p\phi^2)\)).

## Energy-dependent: \(P(z)-S\)

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
process ([[asr_edn]] for \(n\) species). `potential_function = "exact"` replaces
the two-parameter \(P\) by the exact potential function from the radial equation at
each \(z\) ([[asr_pfun]]); \(\lambda=C'/C\) is used so that the poles at the real
zeros of \(C\) cancel and the function is analytic elsewhere (the textbook
\(-P''/(2P')\) has spurious poles in the upper half-plane at the complex zeros
of \(P'\) of an exact, meromorphic \(P\)).

Because \(n_R\) is a projector, \(n_R\mu_R^2=\mu_A^2n_R\), the species-resolved
functions need only the average \(g_{RR}\) and the transition element
\(\langle n_Rg_{RR}\rangle\) from one chain per energy.

## `asr_ed`: the AMscf form

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

At infinite depth \(E-H(E)=T\,(E-H^{(2)})\,T\), so \(G_c\) equals the
species-resolved energy-independent function. The two routes are therefore
exact cross-checks of one another at convergence.

## Tails, depth and cost

| Setting | Effect |
|---|---|
| `ed_terminator = "hull" \| "bp" \| "avg2"` | Tail of each energy's chain; default `hull`, which is the safe choice ([Recursion and terminators](recursion.html)) |
| `depth_tol` | A chain stops early (from level 8, every 4 levels) when \(G\) changes by less than this (relative); `levels` is the maximum |
| `contour_points` | Number of energies (and therefore chains) |
| `fermi_kt` | Smeared occupations; makes \(N(E_F)\) smooth in the SCF |
| `tail_scan` | Diagnostic: moments against termination depth, each tail at its own \(E_F\) |

Statistics on how many chains and levels were used are kept in
[[asr_endep]] (`ed_nchain`, `ed_nlev`).

## An exact reference: KKR-ASA

`potential_function = "exact"` together with `kkr_check = nk` runs, after the SCF,
a same-model k-space KKR-ASA calculation of the *ordered* system on an \(n_k^3\)
mesh, an independent test of the recursion route with exactly the same potential
function.
