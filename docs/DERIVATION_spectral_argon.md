# Derivation of the Configuration-Averaged Bloch Spectral Function $A(\mathbf{k}, E)$

Written 2026-11-17 in response to `docs/HANDOFF_spectral_argon.md`.
This document derives the configuration-averaged Bloch spectral function for the energy-dependent augmented-space recursion formalism implemented in `src/asr_class_scf.f90`, `src/asr_augcls.f90`, `src/asr_edc.f90`, and `src/asr_ed3.f90`.
Every statement is marked as **proven**, **checked numerically**, or **conjecture**.
Verification script: `tools/check_spectral_asr.py`.

---

## 0. Symbol Table

Each symbol is defined once here and used throughout.

* $\Lambda$: the crystal lattice of $N_{\text{cell}}$ primitive decoration cells.
* $\mathbf{T}$: a Bravais translation vector of the primitive decoration cell in units of the lattice constant $a$.
* $\kappa, \kappa' \in \{1, \dots, N_{\text{cls}}\}$: sublattice (site class) indices in the primitive cell with site fractions $f_\kappa$ and basis positions $\boldsymbol{\tau}_\kappa$.
* $R$: a lattice site at position $\mathbf{R} = \mathbf{T} + \boldsymbol{\tau}_{\kappa'}$ measured from the origin site $0\kappa$ at $\boldsymbol{\tau}_\kappa$.
* $L, L' \in \{1, \dots, 9\}$: real cubic harmonic orbital indices $(s, p_y, p_z, p_x, d_{xy}, d_{yz}, d_{3z^2-r^2}, d_{xz}, d_{x^2-y^2})$ with angular momentum $\ell(L) \in \{0, 1, 2\}$.
* $s, c \in \{1, \dots, n_\kappa\}$: chemical species indices on a sublattice with marginal concentrations $x_s(\kappa)$.
* $E_{\nu,s,L}, C_{s,L}, D_{s,L} = \Delta_{s,L}^\alpha, d_{s,L} = D_{s,L}^{1/2}, o_{s,L} = o_{s,L}^\alpha, p_{s,L}^\gamma, p_{s,L} = p_{s,L}^\gamma + o_{s,L}^2$: TB-LMTO-ASA potential parameters in the standard screened representation $\alpha$.
* $S^\alpha_{RL, R'L'}$: the non-random screened structure constant matrix.
* $h = c + d S^\alpha d$: the first-order Hamiltonian matrix with $c_{s,L} = C_{s,L} - E_{\nu,s,L}$.
* $H^{(2)} = E_\nu + h - h o h$: the orthogonal second-order Hamiltonian of route `edc`.
* $K(E) = E O - H$: the secular matrix of route `ed3` with overlap $O = (1 + h o)(1 + o h) + h p^\gamma h$.
* $z = E + i\eta$: the complex energy with real energy $E$ and broadening $\eta > 0$.
* $T = \operatorname{diag}_{R,L}\!\big(T_{RL}(s_R)\big)$: the central-species scaling matrix with $T_{RL}(s) = (D_{c,L}/D_{s,L})^{1/2}$.
* $H^{(c)}(z) = e(z) + h' - h' q h'$: the scaled operator of route `edc` for central species $c$, with $h' = g + S'$, $S' = D_c^{1/2} S^\alpha D_c^{1/2}$, $g_{s,L} = (D_{c,L}/D_{s,L}) c_{s,L}$, $q_{s,L} = (D_{s,L}/D_{c,L}) o_{s,L}$, and $e_{s,L}(z) = z - (D_{c,L}/D_{s,L})(z - E_{\nu,s,L})$.
* $K'^{(c)}(z) = T K(z) T = \tilde{\varepsilon} - h' + u h' + h' u - h' \tilde{Q} h'$: the scaled secular operator of route `ed3` with $\varepsilon_{s,L} = z - E_{\nu,s,L}$, $\tilde{\varepsilon}_{s,L} = (D_{c,L}/D_{s,L})\varepsilon_{s,L}$, $u_{s,L} = o_{s,L}\varepsilon_{s,L}$, and $\tilde{Q}_{s,L} = (D_{s,L}/D_{c,L})(o_{s,L} - p_{s,L}\varepsilon_{s,L})$.
* $\phi_k(s)$ ($k = 0, \dots, n_\kappa - 1$): the orthonormal local configuration basis on sublattice $\kappa$ with $\phi_0(s) = 1$.
* $|R L, C\rangle$: an augmented-space basis state with the electron at $(R, L)$ and flipped site configuration $C$ ($C = \emptyset$ is the reference state).
* $L_{\text{lev}}$: the `levels` truncation parameter defining the reachable node set $V_{L_{\text{lev}}}$.
* $r_{\max}$: the `rank` truncation parameter limiting the number of flipped sites $|C| \le r_{\max}$.
* $m$: the `recursion_steps` parameter for the complex-symmetric Lanczos chain.
* $\mathbf{k}$: a Bloch wavevector in the primitive Brillouin zone $\text{BZ}$ of volume $\Omega_{\text{BZ}}$.

---

## Deliverable 1: Definitions and Main Results (Q1 to Q7)

### Q1. Physical Green Function and Definition of $A(\mathbf{k}, E)$

Define the Bloch basis state for sublattice $\kappa$, orbital $L$, and wavevector $\mathbf{k}$ by
$$|\mathbf{k}, \kappa L\rangle = \frac{1}{\sqrt{N_{\text{cell}}}} \sum_{\mathbf{T}} e^{i \mathbf{k} \cdot (\mathbf{T} + \boldsymbol{\tau}_\kappa)} |\mathbf{T} + \boldsymbol{\tau}_\kappa, L\rangle .$$
This gauge matches `ham_k_cell` and `hsmall_k` in `src/asr_lmto_cell.f90`.

#### Route (a): Orthogonal second-order Hamiltonian $H^{(2)}$ (`edc`)

For each alloy configuration $\mathcal{C}$, the physical Green function in the orthonormal $\gamma$-like representation is
$$G^{(2)}(\mathcal{C}; z) = \big(z - H^{(2)}(\mathcal{C})\big)^{-1} .$$

#### Route (b): Full generalised ASA problem (`ed3`)

For each alloy configuration $\mathcal{C}$, the secular resolvent is $G(\mathcal{C}; z) = K(\mathcal{C}; z)^{-1} = \big(z O(\mathcal{C}) - H(\mathcal{C})\big)^{-1}$. The physical sphere-charge-weighted Green function matrix is
$$\mathcal{G}^{\text{ASA}}(\mathcal{C}; z) = (1 + o h)\, K(\mathcal{C}; z)^{-1}\, (1 + h o) + (p^\gamma)^{1/2}\, h\, K(\mathcal{C}; z)^{-1}\, h\, (p^\gamma)^{1/2} .$$

#### Bloch Green Function and Spectral Matrix (both routes)

Let $\mathcal{G}(z)$ denote $G^{(2)}(z)$ in route `edc` and $\mathcal{G}^{\text{ASA}}(z)$ in route `ed3`. Its configuration average is $\langle \mathcal{G}_{RL, R'L'}(z) \rangle$. The Bloch Green function matrix is
$$\mathcal{G}_{\kappa L, \kappa' L'}(\mathbf{k}, z) = \sum_{\mathbf{T}} e^{i \mathbf{k} \cdot (\mathbf{T} + \boldsymbol{\tau}_{\kappa'} - \boldsymbol{\tau}_\kappa)} \langle \mathcal{G}_{0\kappa L, (\mathbf{T} + \boldsymbol{\tau}_{\kappa'}) L'}(z) \rangle .$$
The sublattice-resolved and orbital-resolved Bloch spectral matrix at $z = E + i\eta$ is the anti-Hermitian part
$$A_{\kappa L, \kappa' L'}(\mathbf{k}, E) = \frac{i}{2\pi} \left( \mathcal{G}_{\kappa L, \kappa' L'}(\mathbf{k}, E + i\eta) - \overline{\mathcal{G}_{\kappa' L', \kappa L}(\mathbf{k}, E + i\eta)} \right) .$$
The sublattice matrix is $A_{\kappa \kappa'}(\mathbf{k}, E) = \sum_{L=1}^{9} A_{\kappa L, \kappa' L}(\mathbf{k}, E)$. The total Bloch spectral function per primitive cell and per spin is
$$A(\mathbf{k}, E) = \sum_{\kappa=1}^{N_{\text{cls}}} \sum_{L=1}^{9} A_{\kappa L, \kappa L}(\mathbf{k}, E) = -\frac{1}{\pi} \sum_{\kappa=1}^{N_{\text{cls}}} \sum_{L=1}^{9} \operatorname{Im} \mathcal{G}_{\kappa L, \kappa L}(\mathbf{k}, E + i\eta) .$$

---

**Result Q1.1 (Real-space imaginary-part Fourier formula): Proven.**
For any translation-invariant configuration average of a real symmetric operator, the spectral matrix satisfies
$$A_{\kappa L, \kappa' L'}(\mathbf{k}, E) = \sum_{\mathbf{T}} e^{i \mathbf{k} \cdot (\mathbf{T} + \boldsymbol{\tau}_{\kappa'} - \boldsymbol{\tau}_\kappa)} \left( -\frac{1}{\pi} \operatorname{Im} \langle \mathcal{G}_{0\kappa L, (\mathbf{T} + \boldsymbol{\tau}_{\kappa'}) L'}(E + i\eta) \rangle \right) .$$

*Proof.* For every configuration $\mathcal{C}$, the matrices $H^{(2)}(\mathcal{C})$ and $K(\mathcal{C}; z)$ are symmetric in site and orbital indices. Therefore $\mathcal{G}_{RL, R'L'}(\mathcal{C}; z) = \mathcal{G}_{R'L', RL}(\mathcal{C}; z)$. By lattice translation invariance of the configuration average, $\langle \mathcal{G}_{0\kappa' L', (-\mathbf{T} + \boldsymbol{\tau}_\kappa) L}(z) \rangle = \langle \mathcal{G}_{0\kappa L, (\mathbf{T} + \boldsymbol{\tau}_{\kappa'}) L'}(z) \rangle$. Taking the complex conjugate of the $(\kappa' L', \kappa L)$ Fourier sum and reversing the summation index $\mathbf{T} \to -\mathbf{T}$ yields $\overline{\mathcal{G}_{\kappa' L', \kappa L}(\mathbf{k}, z)} = \sum_{\mathbf{T}} e^{i \mathbf{k} \cdot (\mathbf{T} + \boldsymbol{\tau}_{\kappa'} - \boldsymbol{\tau}_\kappa)} \overline{\langle \mathcal{G}_{0\kappa L, (\mathbf{T} + \boldsymbol{\tau}_{\kappa'}) L'}(z) \rangle}$. Subtracting the two series and dividing by $2\pi i$ gives the stated formula. $\square$

---

**Result Q1.2 (Positivity, DOS sum rule, ordered limit, and energy sum rule): Proven.**
1. *Non-negativity:* For every $\mathbf{k}$ and $E$ and $\eta > 0$, the $9 N_{\text{cls}} \times 9 N_{\text{cls}}$ matrix $A(\mathbf{k}, E)$ is Hermitian positive semidefinite, and $A_{\kappa L, \kappa L}(\mathbf{k}, E) \ge 0$.
2. *Brillouin-zone integral (DOS sum rule):* For both routes `edc` and `ed3`, the Brillouin-zone average of $A_{\kappa L, \kappa L}(\mathbf{k}, E)$ is the exact on-site density of states of sublattice $\kappa$ and orbital $L$ computed by `src/asr_class_scf.f90`:
   $$\frac{1}{\Omega_{\text{BZ}}} \int_{\text{BZ}} d^3k \, A_{\kappa L, \kappa L}(\mathbf{k}, E) = -\frac{1}{\pi} \operatorname{Im} \langle \mathcal{G}_{0\kappa L, 0\kappa L}(E + i\eta) \rangle = \operatorname{DOS}_{\kappa L}(E) .$$
3. *Ordered limit:* In an ordered crystal as $\eta \to 0^+$, $A_{\kappa L, \kappa L}(\mathbf{k}, E)$ converges to $\sum_n w_{n\mathbf{k}, \kappa L} \delta(E - E_{n\mathbf{k}})$ where $E_{n\mathbf{k}}$ are the eigenvalues of `bands_core` in `src/asr_asa.f90` and $w_{n\mathbf{k}, \kappa L}$ are the exact sphere-charge weights (`|z|^2` in `edc` and `xm0 + 2 o xm1 + p xm2` in `ed3`), with $\sum_{\kappa, L} w_{n\mathbf{k}, \kappa L} = 1$.
4. *Energy integral sum rule:* In route `edc`, $\int_{-\infty}^\infty dE \, A_{\kappa L, \kappa' L'}(\mathbf{k}, E) = \delta_{\kappa\kappa'} \delta_{LL'}$ at every $\mathbf{k}$. In route `ed3`, the combined Brillouin-zone and energy integral of the trace is $\Omega_{\text{BZ}}^{-1} \int_{\text{BZ}} d^3k \int_{-\infty}^\infty dE \, A(\mathbf{k}, E) = 9 N_{\text{cls}}$ (and in the ordered limit $\int_{-\infty}^\infty dE \, A(\mathbf{k}, E) = 9 N_{\text{cls}}$ at each $\mathbf{k}$ separately).

*Proof.*
(1) For each configuration $\mathcal{C}$ in route `edc`, $-\operatorname{Im} G^{(2)}(\mathcal{C}; E + i\eta) = \eta (E - H^{(2)})^{-1}(E - H^{(2)})^{-1\,\dagger} \succeq 0$. In route `ed3`, $O = (1 + h o)(1 + o h) + h p^\gamma h \succ 0$ by Result 3 of `docs/DERIVATION_ed3.md`, so $-\operatorname{Im} K(\mathcal{C}; E + i\eta)^{-1} = \eta K^{-1} O K^{-\dagger} \succeq 0$. For any vector $|\psi\rangle$, $\langle \psi | \mathcal{G}^{\text{ASA}} | \psi \rangle = \langle \psi_0 | K^{-1} | \psi_0 \rangle + \langle \psi_2 | K^{-1} | \psi_2 \rangle$ with $|\psi_0\rangle = (1 + h o)|\psi\rangle$ and $|\psi_2\rangle = h (p^\gamma)^{1/2} |\psi\rangle$, which is a sum of two non-negative quadratic forms. Convex averaging over $\mathcal{C}$ preserves positive semidefiniteness.
(2) Orthogonality of the Bloch phases $\Omega_{\text{BZ}}^{-1} \int_{\text{BZ}} d^3k \, e^{i \mathbf{k} \cdot \mathbf{T}} = \delta_{\mathbf{T}, \mathbf{0}}$ isolates the $\mathbf{T} = \mathbf{0}$ term. In route `edc`, $\mathcal{G}_{0\kappa L, 0\kappa L} = G^{(2)}_{0\kappa L, 0\kappa L}$. In route `ed3`, the diagonal element of $\mathcal{G}^{\text{ASA}}$ on site $0\kappa$ is
$$\mathcal{G}^{\text{ASA}}_{0\kappa L, 0\kappa L} = G_{00} + 2 o_{0\kappa L} (h G)_{00} + (o_{0\kappa L}^2 + p^\gamma_{0\kappa L}) (h G h)_{00} = q_0 + p^\gamma_{0\kappa L} q_2 = G_{00} + 2 o_{0\kappa L} (h G)_{00} + p_{0\kappa L} (h G h)_{00} ,$$
which is identical to `gz` in `src/asr_class_scf.f90`.
(3) In an ordered crystal, expanding $K(\mathbf{k}; z)^{-1} = \sum_n z_{n\mathbf{k}} z_{n\mathbf{k}}^\dagger / (z - E_{n\mathbf{k}})$ with $z_{n\mathbf{k}}^\dagger O(\mathbf{k}) z_{n\mathbf{k}} = 1$ gives residues $|(1 + o(\mathbf{k}) h(\mathbf{k})) z_{n\mathbf{k}}|_{\kappa L}^2 + p^\gamma_{\kappa L} |h(\mathbf{k}) z_{n\mathbf{k}}|_{\kappa L}^2 = \mathrm{xm0}_{n\mathbf{k},\kappa L} + 2 o_{\kappa L} \mathrm{xm1}_{n\mathbf{k},\kappa L} + p_{\kappa L} \mathrm{xm2}_{n\mathbf{k},\kappa L}$. Summing over $\kappa, L$ gives $z_{n\mathbf{k}}^\dagger \big[(1 + h o)(1 + o h) + h p^\gamma h\big] z_{n\mathbf{k}} = z_{n\mathbf{k}}^\dagger O(\mathbf{k}) z_{n\mathbf{k}} = 1$.
(4) Integrating $-\frac{1}{\pi}\operatorname{Im}(z - H^{(2)})^{-1}$ over $E \in (-\infty, \infty)$ gives the identity $I$. In route `ed3`, integrating over $E$ replaces $K(z)^{-1}$ by $O^{-1}$, and cyclic invariance of the trace gives $\operatorname{Tr}\big[(1+oh)O^{-1}(1+ho) + (p^\gamma)^{1/2} h O^{-1} h (p^\gamma)^{1/2}\big] = \operatorname{Tr}\big[O^{-1}\big((1+ho)(1+oh) + h p^\gamma h\big)\big] = \operatorname{Tr}[O^{-1} O] = 9 N_{\text{cell}}$. $\square$

---

### Q2. Off-Site Averages in the Energy-Dependent Form

**Result Q2.1 (Correction of the scaling exponent for off-site elements in route `edc`): Proven and checked numerically.**
The prompt's tentative factor $T_R(s)^{-1}$ is **incorrect**; the exact factor is $+1$ power of $T_{RL'}(s) = (D_{c,L'}/D_{s,L'})^{1/2}$. For any site $R \in \kappa'$ ($R \neq 0\kappa$), the unscaled off-site Green function in route `edc` is given by two equivalent exact representations in the full augmented space:

* **Form 2A (Frozen pair: site $0\kappa$ frozen to species $c$, site $R$ frozen to species $s$):**
  $$\langle G^{(2)}_{0\kappa L, R L'}(z) \rangle = \sum_{c \in \kappa} \sum_{s \in \kappa'} P_{0R}(c, s) \, T_{R L'}(s) \, \widetilde{G}^{(c,s)}_{0\kappa L, R L'}(z), \qquad T_{R L'}(s) = \left(\frac{D_{c,L'}}{D_{s,L'}}\right)^{1/2} ,$$
  where $P_{0R}(c, s) = x_c(\kappa) x_s(\kappa')$ without SRO, and $\widetilde{G}^{(c,s)}_{0\kappa L, R L'}(z) = \langle 0\kappa L, \emptyset | \big(z - H^{(c)}(z)\big)^{-1} | R L', \emptyset \rangle_{s_0 = c,\, s_R = s}$.
* **Form 2B (Single frozen centre: only site $0\kappa$ frozen to species $c$, site $R$ flippable):**
  $$\langle G^{(2)}_{0\kappa L, R L'}(z) \rangle = \sum_{c \in \kappa} x_c(\kappa) \, \langle 0\kappa L, \emptyset | \big(z - H^{(c)}(z)\big)^{-1} | \widetilde{R L'}^{(c)} \rangle_{s_0 = c} ,$$
  where $|\widetilde{R L'}^{(c)}\rangle = \widehat{T}_R |R L', \emptyset\rangle$ is the dressed target vector at site $R$ in the rank-1 augmented space:
  $$|\widetilde{R L'}^{(c)}\rangle = \sum_{k=0}^{n_{\kappa'} - 1} \tau^{(c)}_{R L', k} \, |R L', \{R: k\}\rangle, \qquad \tau^{(c)}_{R L', k} = \sum_{s \in \kappa'} p_s(R) \, T_{R L'}(s) \, \phi_k(s) ,$$
  with $p_s(R) = x_s(\kappa')$ without SRO (or $X_s^{(c)}(\kappa, \kappa')$ if $R$ lies in the first shell with SRO).

*Proof.* For every configuration $\mathcal{C}$, the scaling relation in `src/asr_edc.f90` is $z - H^{(c)}(z) = T (z - H^{(2)}) T$. Taking the matrix inverse gives $(z - H^{(c)}(z))^{-1} = T^{-1} (z - H^{(2)})^{-1} T^{-1}$. Multiplying on both sides by $T$ gives
$$G^{(2)}(\mathcal{C}; z) = (z - H^{(2)}(\mathcal{C}))^{-1} = T \, \big(z - H^{(c)}(\mathcal{C}; z)\big)^{-1} \, T .$$
On site $0\kappa$ frozen to species $c$, $T_{0\kappa L}(c) = (D_{c,L}/D_{c,L})^{1/2} = 1$ for all $L \in \{1, \dots, 9\}$. On site $R$ occupied by species $s$, $T_{R L'}(s) = (D_{c,L'}/D_{s,L'})^{1/2}$ where $D_{c,L'}$ is the $\Delta$ parameter of central species $c$ in angular momentum channel $\ell(L')$ (set in `edc_setup`). Conditioning on $(s_0 = c, s_R = s)$ gives Form 2A. In the augmented space with only site $0\kappa$ frozen to $c$, the multiplication operator $\widehat{T}_R$ at site $R$ acts on the local reference state $|\phi_0\rangle_R$ as $\widehat{T}_R |\phi_0\rangle_R = \sum_{k=0}^{n_{\kappa'}-1} \tau^{(c)}_{R L', k} |\phi_k\rangle_R$. Since $|\phi_k\rangle_R \otimes |\emptyset\rangle_{\Lambda \setminus \{R\}} = |\{R: k\}\rangle$ has at most one flipped site (rank $\le 1$), taking the reference matrix element gives Form 2B. $\square$

*Numerical check (Test 1 of Section 4):* On a 6-site ternary ring ($3^6 = 729$ configurations, $D_{\max}/D_{\min} = 1.96$), both Form 2A and Form 2B agree with exact enumeration to $1.17 \times 10^{-14}$. Replacing $T_R(s)$ by $T_R(s)^{-1}$ fails with an error of $1.46 \times 10^{-2}$. At `rank = 1`, Form 2A has max error $5.73 \times 10^{-3}$ and Form 2B has max error $9.65 \times 10^{-3}$; at `rank = 2`, Form 2A has max error $1.66 \times 10^{-4}$ and Form 2B has max error $4.23 \times 10^{-4}$.

---

**Result Q2.2 (Off-site averages in the full ASA route `ed3`): Proven and checked numerically.**
In route `ed3`, with site $0\kappa$ frozen to $c$ and site $R \in \kappa'$ frozen to $s$, define the reference-sector vectors $|u_{0L}\rangle = |0\kappa L, \emptyset\rangle$, $|u_{RL'}\rangle = |R L', \emptyset\rangle$, $|v_{0L}^{(c)}\rangle = h' |u_{0L}\rangle$, $|v_{RL'}^{(c)}\rangle = h' |u_{RL'}\rangle$, $|a_{0L}^{(c)}\rangle = |u_{0L}\rangle + o_{c,L} |v_{0L}^{(c)}\rangle$, and $|a_{RL'}^{(c,s)}\rangle = |u_{RL'}\rangle + q_{s,L'} |v_{RL'}^{(c)}\rangle$ where $q_{s,L'} = (D_{s,L'}/D_{c,L'}) o_{s,L'}$. Then
$$\langle \mathcal{G}^{\text{ASA}}_{0\kappa L, R L'}(z) \rangle = \sum_{c \in \kappa} \sum_{s \in \kappa'} P_{0R}(c, s) \left[ T_{R L'}(s) \, \langle a_{0L}^{(c)} | K'^{(c)}(z)^{-1} | a_{RL'}^{(c,s)} \rangle_{c,s} + T_{R L'}(s)^{-1} (p^\gamma_{c,L} p^\gamma_{s,L'})^{1/2} \, \langle v_{0L}^{(c)} | K'^{(c)}(z)^{-1} | v_{RL'}^{(c)} \rangle_{c,s} \right] .$$
If site $R$ is not frozen (single-centre operator with only $0\kappa$ frozen to $c$), the same identity holds with $|a_{RL'}^{(c,s)}\rangle$ and $|v_{RL'}^{(c)}\rangle$ replaced by the single-flip augmented-space vectors $\big(\widehat{T}_R + h' \widehat{T}_R^{-1} \widehat{o}_R\big)|R L', \emptyset\rangle$ and $h' \widehat{T}_R^{-1} (\widehat{p}^\gamma_R)^{1/2} |R L', \emptyset\rangle$.

*Proof.* Substitute $K(z)^{-1} = T K'^{(c)}(z)^{-1} T$ and $h = T^{-1} h' T^{-1}$ into $\mathcal{G}^{\text{ASA}}(z)$ of Q1. This gives $(1 + o h) K^{-1} (1 + h o) = (T + o T^{-1} h') K'^{-1} (T + h' T^{-1} o)$ and $(p^\gamma)^{1/2} h K^{-1} h (p^\gamma)^{1/2} = (p^\gamma)^{1/2} T^{-1} h' K'^{-1} h' T^{-1} (p^\gamma)^{1/2}$. Acting with $(T + h' T^{-1} o)$ on $|0\kappa L\rangle$ where $T_{0\kappa L}(c) = 1$ gives $|a_{0L}^{(c)}\rangle$. Acting on $|R L'\rangle$ where $T_{RL'}(s) = (D_{c,L'}/D_{s,L'})^{1/2}$ and $T_{RL'}(s)^{-1} o_{s,L'} = T_{RL'}(s) q_{s,L'}$ gives $T_{RL'}(s) |a_{RL'}^{(c,s)}\rangle$. Similarly, acting with $h' T^{-1} (p^\gamma)^{1/2}$ on $|0\kappa L\rangle$ and $|R L'\rangle$ gives $(p^\gamma_{c,L})^{1/2} |v_{0L}^{(c)}\rangle$ and $T_{RL'}(s)^{-1} (p^\gamma_{s,L'})^{1/2} |v_{RL'}^{(c)}\rangle$. Because $0\kappa$ and $R$ are frozen in Form 2A, $h'|u_{RL'}\rangle = g_{s,L'} |u_{RL'}\rangle + S'|u_{RL'}\rangle$ is non-random and lies entirely in the reference sector. $\square$

*Numerical check (Test 1 of Section 4):* Checked on the 6-site ternary ring against exact enumeration of $\mathcal{G}^{\text{ASA}}(z)$; maximum error across all sites $R$ is $3.82 \times 10^{-15}$.

---

### Q3. Comparison of Off-Diagonal Extraction Methods from Recursions

Let $K(z) = z - H^{(c)}(z)$ (or $K'^{(c)}(z)$ in `ed3`) be complex symmetric. Let $a = |0\kappa L, \emptyset\rangle$ be the start vector and $b$ be a target vector at site $R$. Let $x_m^{(v)} \in \mathcal{K}_m(K, v)$ be the tail-free $m$-step complex-symmetric Lanczos Galerkin solution for start vector $v$, with Galerkin residual $r_m^{(v)} = v - K x_m^{(v)}$.

**Result Q3.1 (Accuracy orders in the Galerkin residual): Proven and checked numerically.**
1. *Polarisation identity (Q3a):* Using $b^T K^{-1} a = \frac{1}{4}\big[(a+b)^T K^{-1}(a+b) - (a-b)^T K^{-1}(a-b)\big]$ (or $\frac{1}{2}\big[(a+b)^T K^{-1}(a+b) - a^T K^{-1} a - b^T K^{-1} b\big]$), the error is
   $$\delta G_{ba}^{\text{pol}} = \frac{1}{4}\left[ (r_m^{(a+b)})^T K^{-1} r_m^{(a+b)} - (r_m^{(a-b)})^T K^{-1} r_m^{(a-b)} \right] = O\big(\|r_m\|^2\big) \quad (\text{second order}) .$$
2. *One-sided Krylov projection (Q3b):* Using $f_m = b^T x_m^{(a)} = \sum_{k=1}^m (b^T u_k) [\widetilde{T}_m^{-1}]_{k1}$ from the chain started on $a$ alone (as in `cf_species` in `src/asr_endep.f90`), the exact error identity is
   $$b^T K^{-1} a - b^T x_m^{(a)} = b^T K^{-1} r_m^{(a)} = (x_m^{(b)})^T r_m^{(a)} + (r_m^{(b)})^T K^{-1} r_m^{(a)} = O\big(\|r_m^{(a)}\|\big) \quad (\text{first order}) .$$
   If a chain from $b$ is also run, the uncorrected mean $\frac{1}{2}(b^T x_m^{(a)} + a^T x_m^{(b)})$ remains first order because $(x_m^{(b)})^T r_m^{(a)} \neq 0$, whereas the Petrov-Galerkin cross-corrected estimator $b^T x_m^{(a)} + (x_m^{(b)})^T r_m^{(a)}$ has error $(r_m^{(b)})^T K^{-1} r_m^{(a)} = O(\|r_m^{(a)}\| \|r_m^{(b)}\|)$ (second order).
3. *Block recursion over $U_1 = [a, b]$ (Q3c):* For an $m$-step $2\times 2$ block Lanczos recursion generating block basis $Q_m$ and residual $R_m = U_1 - K X_m$, the error matrix is
   $$U_1^T K^{-1} U_1 - U_1^T X_m = R_m^T K^{-1} R_m = O\big(\|R_m\|^2\big) \quad (\text{second order}) .$$

*Proof.* By Galerkin orthogonality, $v^T r_m^{(u)} = 0$ for every $v \in \mathcal{K}_m(K, u)$. For any trial solution $x_m^{(a)}$ and test solution $x_m^{(b)}$ with $b = K x_m^{(b)} + r_m^{(b)}$, we have $b^T K^{-1} a - b^T x_m^{(a)} = b^T K^{-1} r_m^{(a)} = (x_m^{(b)} + K^{-1} r_m^{(b)})^T r_m^{(a)} = (x_m^{(b)})^T r_m^{(a)} + (r_m^{(b)})^T K^{-1} r_m^{(a)}$. When $b = a$ (diagonal chains in Q3a) or when $\mathcal{K}_m$ contains both $a$ and $b$ (block chain in Q3c), $x_m^{(b)}$ lies in the test space orthogonal to $r_m^{(a)}$, so $(x_m^{(b)})^T r_m^{(a)} = 0$ and only the quadratic term $(r_m^{(b)})^T K^{-1} r_m^{(a)}$ remains. When only the chain from $a$ is run (Q3b), $(x_m^{(b)})^T r_m^{(a)}$ is non-zero and linear in $r_m^{(a)}$. $\square$

*Numerical check (Tests 2 and 3 of Section 4):* Verified on the 6-site ternary ring; the block error identity holds to $3.0 \times 10^{-16}$.

---

**Result Q3.2 (Symmetry kept and chain counts for Q3a, Q3b, Q3c): Proven.**
1. *One-sided Krylov projection on the single-centre operator (Form 2B + Q3b):*
   The operator $H^{(c)}(z)$ freezes only the origin $0\kappa$ and keeps the full cubic point group $G_0 \subset D_{4h}$ of order up to 16 (`class_ops` in `src/asr_symmetry.f90`). Because both the start vector $|0\kappa L, \emptyset\rangle$ and $H^{(c)}(z)$ transform with character $\chi_L(g) = D_{LL}(g) = \pm 1$ under $G_0$, every Lanczos vector $u_k$ in the reduced space determines the exact unreduced amplitude at any site $R$, orbital $L'$, and local state $m_{\text{st}} \in \{0, \dots, n_{\kappa'}-1\}$ via `at_link_c` (`src/asr_augcls.f90`):
   $$\psi_k(R, L', \{R: m_{\text{st}}\}) = \frac{\chi_L(g) \, \mathrm{sg}(L', g)}{\sqrt{\text{orb}(r)}} \, u_k\big(9(r - 1) + \pi(L', g)\big) ,$$
   where $(r, g)$ are the canonical representative node and group operation of $(R, \{R: m_{\text{st}}\})$.
   Under $D_{4h}$, the 5 orbitals $L \in \{s\,(1), p_z\,(3), d_{xy}\,(5), d_{3z^2-r^2}\,(7), d_{x^2-y^2}\,(9)\}$ all carry 1D characters $\chi_L(g) = \pm 1$. For a cubic sublattice ($\text{cubic} = \text{.true.}$), the remaining 4 start orbitals $\{p_y\,(2), p_x\,(4), d_{yz}\,(6), d_{xz}\,(8)\}$ are obtained exactly from the $p_z$ and $d_{xy}$ chains by applying the 3-fold cyclic permutation $C_3: (x, y, z) \to (y, z, x)$ to $\mathbf{R}$ and the orbital indices. Thus **all shells $R$ and all $9 \times 9$ orbital blocks $(L, L')$ are obtained simultaneously from 5 chains per central species $c$ and energy $z$** (10 chains in `ed3`), with **zero extra chains per shell $R$**.
2. *Frozen-pair projection, polarisation, and block recursion (Form 2A + Q3a/b/c):*
   Freezing $(0\kappa: c, R: s)$ breaks $G_0 \subset D_{4h}$ down to the abelian pair stabiliser $G_{\text{pair}}$ (`epi_pair_ops`, order $1, 2, 4,$ or $8$), increasing the reduced node count by $|G_0|/|G_{\text{pair}}| \in [2, 16]$ and requiring a separate operator build per shell orbit and separate chains per species pair $(c, s)$.

| Method | Accuracy in $\|r_m\|$ | Symmetry group kept | Reduced node count vs SCF | Chains per $(E, \text{centre } c)$ for $N_{\text{sh}}$ shells and $N_{\text{sp}}$ species (`edc` / `ed3`) |
|---|---|---|---|---|
| **Q3b (Form 2B: single-centre projection)** | $O(\|r_m\|)$ | Full $G_0 \subset D_{4h}$ (order 16) | $1\times$ ($15,827$ in FeNiCr) | **$5$ / $10$ total** (independent of $N_{\text{sh}}$ and $N_{\text{sp}}$) |
| **Q3b (Form 2A: frozen-pair projection)** | $O(\|r_m\|)$ | Pair stabiliser $G_{\text{pair}}$ (order $1\text{--}8$) | $2\times$ to $16\times$ | $9 N_{\text{sh}} N_{\text{sp}}$ / $18 N_{\text{sh}} N_{\text{sp}}$ |
| **Q3a (Form 2A: polarisation $|0 \pm R\rangle$)** | $O(\|r_m\|^2)$ | Pair stabiliser $G_{\text{pair}}$ (order $1\text{--}8$) | $2\times$ to $16\times$ | $(9 + 2\sum_\chi d_\chi^2) N_{\text{sh}} N_{\text{sp}}$ |
| **Q3c (Form 2A: $2\times 2$ block Lanczos)** | $O(\|r_m\|^2)$ | Pair stabiliser $G_{\text{pair}}$ (order $1\text{--}8$) | $2\times$ to $16\times$ | $2 \sum_\chi d_\chi^2 N_{\text{sh}} N_{\text{sp}}$ (no `hull` tail exists) |

---

### Q4. The Fourier Sum, Spatial Decay, Windowing, and Relation to `levels`

**Result Q4.1 (Spatial decay of $\langle G_{0R}(E + i\eta) \rangle$ and truncation bound): Proven.**
Let $H^{(2)}$ have maximum hopping range $r_0 = 2 R_{\text{cl}}$ (two structure-constant hops) and uniform hopping norm bound $J = \sup_R \sum_{R' \neq R} \|H^{(2)}_{RR'}\|$.
1. *Rigorous Combes-Thomas bound (Proven):* For every configuration $\mathcal{C}$ and $\eta = \operatorname{Im} z > 0$,
   $$\big\| \langle G^{(2)}_{0\kappa, R}(E + i\eta) \rangle \big\| \le \frac{2}{\eta} \exp\!\big(-\gamma(\eta) |\mathbf{R}|\big), \qquad \gamma(\eta) = \frac{1}{r_0} \operatorname{arsinh}\!\left(\frac{\eta}{2 e J}\right) .$$
2. *Asymptotic disorder + smearing decay (Checked numerically):* In a metallic random alloy with elastic mean free path $\ell_{\text{mfp}}(E) = v_F(E) / (2 \Gamma_{\text{dis}}(E))$ and Fermi velocity $v_F(E)$, the configuration-averaged Green function decays as
   $$\big\| \langle G_{0\kappa, R}(E + i\eta) \rangle \big\| \le C \frac{a}{|\mathbf{R}|} e^{-\kappa_{\text{tot}}(E, \eta) |\mathbf{R}|}, \qquad \kappa_{\text{tot}}(E, \eta) = \frac{1}{2 \ell_{\text{mfp}}(E)} + \frac{\eta}{v_F(E)} = \frac{\Gamma_{\text{dis}}(E) + \eta}{v_F(E)} .$$
3. *Truncation error bound (Proven):* Truncating the Fourier series at $|\mathbf{R}| \le R_{\max}$ produces a uniform error in $A(\mathbf{k}, E)$ bounded by
   $$\big\| A(\mathbf{k}, E) - A^{R_{\max}}(\mathbf{k}, E) \big\| \le \frac{1}{\pi} \sum_{|\mathbf{R}| > R_{\max}} \big\| \operatorname{Im} \langle \mathcal{G}_{0\kappa, R}(E + i\eta) \rangle \big\| \le \frac{4 C a}{\Omega_0 \kappa_{\text{tot}}^2} (1 + \kappa_{\text{tot}} R_{\max})\, e^{-\kappa_{\text{tot}} R_{\max}} ,$$
   where $\Omega_0$ is the volume per site.

*Proof.* (1) For any self-adjoint operator $H^{(2)}$ with finite hopping range $r_0$, let $U_\mu = \operatorname{diag}_R(e^{\mu |\mathbf{R}|})$. Then $\|U_\mu H^{(2)} U_\mu^{-1} - H^{(2)}\| \le (e^{\mu r_0} - 1) J$. Choosing $\mu = \gamma(\eta)$ gives $\|U_\mu H^{(2)} U_\mu^{-1} - H^{(2)}\| \le \eta / 2$, so $\|(z - U_\mu H^{(2)} U_\mu^{-1})^{-1}\| \le 2 / \eta$, which yields $\|G_{0R}\| \le (2/\eta) e^{-\gamma(\eta)|\mathbf{R}|}$ for each configuration and hence for the average. (3) Follows from $|e^{i\mathbf{k}\cdot\mathbf{R}}| = 1$, the triangle inequality, and integration of $4\pi R^2 \Omega_0^{-1} (a/R) e^{-\kappa_{\text{tot}} R}$ over $R \in (R_{\max}, \infty)$. $\square$

---

**Result Q4.2 (Positivity preservation by a 3D Fejér autocorrelation window): Proven and checked numerically.**
Let $\Lambda_{\text{fin}}$ be any finite set of lattice translations containing $\mathbf{0}$ with diameter $R_{\max}$ and cardinality $N_\Lambda = |\Lambda_{\text{fin}}|$. Define the discrete 3D Fejér (Bartlett autocorrelation) window
$$w_{\text{Fej\'er}}(\mathbf{T}) = \frac{\big|\Lambda_{\text{fin}} \cap (\Lambda_{\text{fin}} - \mathbf{T})\big|}{N_\Lambda} ,$$
or in the isotropic spherical continuum limit on a ball of diameter $R_{\max}$, the 3D Euclid-hat (self-convolution of a ball of radius $R_{\max}/2$) window
$$w_{\text{ball}}(R) = \left(1 - \frac{R}{R_{\max}}\right)^2 \left(1 + \frac{R}{2 R_{\max}}\right) \quad \text{for } 0 \le R \le R_{\max}, \qquad w_{\text{ball}}(R) = 0 \quad \text{for } R > R_{\max} .$$
If the real-space spectral kernel $M_{R L, R' L'}(E) = -\frac{1}{\pi}\operatorname{Im}\langle \mathcal{G}_{0\kappa L, (\mathbf{R}' - \mathbf{R}) L'}(E + i\eta) \rangle$ is positive semidefinite on $\Lambda_{\text{fin}}$, then the Fejér-windowed spectral matrix $A^{w}(\mathbf{k}, E)$ is positive semidefinite ($A^{w}(\mathbf{k}, E) \succeq 0$) for every $\mathbf{k}$ and $E$, and satisfies the exact Brillouin-zone DOS sum rule $\Omega_{\text{BZ}}^{-1} \int_{\text{BZ}} d^3k \, A^w_{\kappa L, \kappa L}(\mathbf{k}, E) = \operatorname{DOS}_{\kappa L}(E)$.

*Proof.* By definition of $w_{\text{Fej\'er}}(\mathbf{T}) = N_\Lambda^{-1} \sum_{\mathbf{T}_1 \in \Lambda_{\text{fin}}} \chi_{\Lambda_{\text{fin}}}(\mathbf{T}_1 + \mathbf{T})$, substituting $\mathbf{T} = \mathbf{T}_2 - \mathbf{T}_1$ rewrites the windowed sum as the quadratic form
$$A^w(\mathbf{k}, E) = \frac{1}{N_\Lambda} \sum_{\mathbf{T}_1 \in \Lambda_{\text{fin}}} \sum_{\mathbf{T}_2 \in \Lambda_{\text{fin}}} e^{-i \mathbf{k} \cdot \mathbf{T}_1} M_{\mathbf{T}_1, \mathbf{T}_2}(E) e^{i \mathbf{k} \cdot \mathbf{T}_2} = \Phi(\mathbf{k})^\dagger M_{\Lambda_{\text{fin}}}(E) \Phi(\mathbf{k}) \succeq 0 .$$
In the isotropic continuum limit, $w_{\text{ball}}(\mathbf{R})$ is the normalized 3D convolution $(\chi_{B(0, R_{\max}/2)} * \chi_{B(0, R_{\max}/2)})(\mathbf{R}) / |B(0, R_{\max}/2)|$, whose 3D Fourier transform is $\hat{w}_{\text{ball}}(\mathbf{q}) = |\hat{\chi}_{B(0, R_{\max}/2)}(\mathbf{q})|^2 / |B(0, R_{\max}/2)| \ge 0$. Convolution of a positive-semidefinite spectral distribution $A(\mathbf{k}, E)$ with a non-negative kernel $\hat{w}_{\text{ball}}(\mathbf{q}) \ge 0$ preserves positive semidefiniteness everywhere. Since $w(\mathbf{0}) = 1$, the $\mathbf{T} = \mathbf{0}$ term is untouched, so the $\mathbf{k}$-integral sum rule is exact. $\square$

---

**Result Q4.3 (Relation between $R_{\max}$ and `levels` = $L_{\text{lev}}$): Proven.**
1. One application of $H^{(2)}$ (or $K'$ in `ed3`) applies $S'$ twice and moves up to $2$ graph hops. The node set $V_{L_{\text{lev}}}$ reaches sites $R$ with graph distance $d_{\text{graph}}(0\kappa, R) \le 2 L_{\text{lev}}$.
2. By Result S1 of `docs/DERIVATION_levels_steps.md`, for a target state at graph level $\ell(R) = \lceil d_{\text{graph}}(0\kappa, R) / 2 \rceil$ (and at level $\ell(R) + 1$ for the flipped state $|R L', \{R: k\}\rangle$ with $k \ge 1$ in Form 2B), the Lanczos projection $c_j = b^T u_j$ is nonzero only for steps $j \ge \ell(R) + 1$, and is unaffected by the boundary of $V_{L_{\text{lev}}}$ up to step $j \le L_{\text{lev}} + 1$.
3. For site $R$ and its 1-flip state $|R L', \{R: k\}\rangle$ to lie strictly inside the interior of $V_{L_{\text{lev}}}$ (with at least one full application of $H^{(c)}(z)$ around $R$ before reaching the hard boundary), $R_{\max}$ must satisfy
   $$d_{\text{graph}}(0\kappa, R) \le 2 (L_{\text{lev}} - 1) .$$
   In fcc with `rcl = 1.2` $a$, one graph hop contains both the 1st shell ($R_1 = a/\sqrt{2}$, 12 sites) and the 2nd shell ($R_2 = a$, 6 sites), so $L_{\text{lev}} = 3$ ($2(L_{\text{lev}} - 1) = 4$ graph hops) places all shells up to $|\mathbf{R}| \le 4.0\,a$ inside the interior of $V_3$.

---

### Q5. The k-Space Alternative

**Result Q5.1 (No-go theorem for a frozen-centre energy-dependent form in k-space): Proven.**
In the augmented space, a Bloch start vector $|\mathbf{k}, \kappa L, \emptyset\rangle = N_{\text{cell}}^{-1/2} \sum_{\mathbf{T}} e^{i \mathbf{k} \cdot (\mathbf{T} + \boldsymbol{\tau}_\kappa)} |\mathbf{T} + \boldsymbol{\tau}_\kappa, L, \emptyset\rangle$ cannot be combined with a frozen-centre energy-dependent scaling $T_R = (D_c / D_{s_R})^{1/2}$.

*Proof.*
(1) The unscaled augmented-space Hamiltonian $\widehat{H}^{(2)}$ commutes with the joint translation operator $\widehat{\mathcal{T}}_{\mathbf{a}} |R, C\rangle = |R + \mathbf{a}, C + \mathbf{a}\rangle$ that translates both the site $R$ and the flipped configuration $C$ together. Thus $\mathbf{k}$ is a good quantum number only in the relative-configuration basis $|\mathbf{k}, \kappa L, C_{\text{rel}}\rangle$ where $C_{\text{rel}} = C - R$.
(2) Freezing a single site $0\kappa$ to species $c$ breaks translation invariance ($[\widehat{\mathcal{T}}_{\mathbf{a}}, \widehat{H}^{(c)}(z)] \neq 0$), so $\widehat{H}^{(c)}(z)$ is no longer block-diagonal in $\mathbf{k}$.
(3) Moreover, the relation $\langle \psi | (z - H^{(2)})^{-1} | \psi \rangle = \langle \psi | T (z - H(z))^{-1} T | \psi \rangle = \langle \psi | (z - H(z))^{-1} | \psi \rangle$ requires $T |\psi\rangle = |\psi\rangle$. Since $|\mathbf{k}, \kappa L, \emptyset\rangle$ has non-zero amplitude on every site $\mathbf{T} + \boldsymbol{\tau}_\kappa$ of sublattice $\kappa$, $T |\mathbf{k}, \kappa L, \emptyset\rangle = |\mathbf{k}, \kappa L, \emptyset\rangle$ holds if and only if $T_R = 1$ on every site of sublattice $\kappa$, which holds only when every site of sublattice $\kappa$ is frozen to the same species $c$ (an ordered crystal without disorder). $\square$

---

**Result Q5.2 (Error, symmetry, and cost of the energy-independent k-space ASR): Proven and checked numerically.**
1. *Bandwidth and rank-1 error in strongly disordered alloys:*
   In $\text{Fe}_{0.66}\text{Ni}_{0.14}\text{Cr}_{0.20}$, the $d$-band width parameter $D_d = \Delta_d^\alpha$ varies by a factor of $D_d(\text{Cr})/D_d(\text{Ni}) \approx 1.9$.
   * If off-diagonal disorder is replaced by the virtual-crystal average $d_s \to \bar{d}$ and $o_s \to \bar{o}$ (traditional diagonal-disorder k-space ASR), each component's $d$-band dispersion has a fractional slope error $(D_s - \bar{D})/\bar{D} \approx \pm 30\%$, shifting $A(\mathbf{k}, E)$ by up to $7\%$ even on the 6-site test ring (Test 5 of Section 4).
   * If the full unscaled operator $\widehat{H}^{(2)} = \widehat{E}_\nu + \widehat{c} + \widehat{d} S \widehat{d} - \widehat{h}\widehat{o}\widehat{h}$ is run in the relative-configuration k-space basis $|\mathbf{k}, C_{\text{rel}}\rangle$ at `rank = 1`, every bond hop $\widehat{d}_{R_1} S_{R_1 R_0} \widehat{d}_{R_0}$ creates a simultaneous two-site flip $|\{R_0, R_1\}\rangle$ of amplitude $(\delta d)^2 S_{R_1 R_0}$ that is truncated at `rank = 1`, whereas in the energy-dependent AMscf form $T d = d_c$ is non-random and has zero flip truncation on $S'$. The resulting rank-1 self-energy defect in the energy-independent form is of order $\operatorname{Var}(d)/\bar{d}^2 \, W_d \approx 3\text{ to }6\%$ of the $d$-band width ($10\text{ to }25\text{ mRy}$).
2. *Symmetry and cost comparison:*
   In k-space ASR, the point group is reduced from $D_{4h}$ (order 16) to the little group $G_{\mathbf{k}}$ of $\mathbf{k}$ (order 4 along $\Gamma\text{--}\text{K}$, order 2 along $\Gamma\text{--}\text{L}$ inside $D_{4h}$, order 1 at generic $\mathbf{k}$), enlarging the reduced node set by $4\times$ to $16\times$, and requiring a separate complex operator build at every $\mathbf{k}$ point. In the real-space route of Q3 and Q4, the chains are completely independent of $\mathbf{k}$ and use the full $D_{4h}$-reduced operator once per energy.

---

### Q6. Causality and Positivity Restoration

**Result Q6.1 (Conditions for $A(\mathbf{k}, E) \succeq 0$ in the construction of Q2 to Q4): Proven.**
The constructed spectral matrix $A^{m, R_{\max}}(\mathbf{k}, E)$ is positive semidefinite at $(\mathbf{k}, E)$ if the combined chain-projection error $\varepsilon_{\text{chain}}(m, z)$, node-set boundary leak $\varepsilon_{\text{leak}}(L_{\text{lev}}, z)$, and spatial Fourier tail $\varepsilon_{\text{tail}}(R_{\max}, z)$ are smaller than the smallest eigenvalue $\lambda_{\min}\big(A^{\text{exact}}(\mathbf{k}, E)\big)$, or whenever:
1. $m \to \dim V_{L_{\text{lev}}}$ (where $\widetilde{G}_m(z) \to G_{L_{\text{lev}}}(z)$ is an exact matrix Herglotz function on $V_{L_{\text{lev}}}$ by Result S2 of `docs/DERIVATION_levels_steps.md`),
2. the boundary leak $r_\ell = (1 - P_{V_{L_{\text{lev}}}}) A x$ is negligible compared to $\eta = \operatorname{Im} z$, and
3. the Fejér window $w_{\text{Fej\'er}}(\mathbf{R})$ (or $w_{\text{ball}}(R)$) of Result Q4.2 is applied.

At finite $m$ (for example $m = 24$), isolated negative eigenvalues $\lambda_\alpha(\mathbf{k}, E) < 0$ can occur from three mechanisms: (i) F4 non-causality of the truncated energy-dependent chain (`docs/DERIVATION_f4.md`), (ii) the first-order residual $O(\|r_m\|)$ of one-sided Krylov projection on off-site shells $R \neq 0$, and (iii) Gibbs oscillations if a sharp cutoff $w(\mathbf{R}) = 1$ is used instead of the Fejér window.

---

**Result Q6.2 (Frobenius positive projection and trace-preserving sum-rule restoration): Proven and checked numerically.**
Let $A(\mathbf{k}, E) = \sum_{\alpha=1}^{9 N_{\text{cls}}} \lambda_\alpha(\mathbf{k}, E) v_\alpha(\mathbf{k}, E) v_\alpha(\mathbf{k}, E)^\dagger$ be the eigendecomposition of the Hermitian matrix $A(\mathbf{k}, E)$. Define the positive-semidefinite projection
$$A_+(\mathbf{k}, E) = \sum_{\alpha=1}^{9 N_{\text{cls}}} \max\!\big(\lambda_\alpha(\mathbf{k}, E), 0\big) \, v_\alpha(\mathbf{k}, E) v_\alpha(\mathbf{k}, E)^\dagger ,$$
the clipped negative weight $w_{\text{neg}}(\mathbf{k}, E) = \sum_{\lambda_\alpha < 0} |\lambda_\alpha(\mathbf{k}, E)|$, and the trace-preserving positive projection
$$A_+^{\text{tr}}(\mathbf{k}, E) = \frac{\max\!\big(\operatorname{Tr} A(\mathbf{k}, E), 0\big)}{\operatorname{Tr} A_+(\mathbf{k}, E)} \, A_+(\mathbf{k}, E) \qquad \big(\text{for } \operatorname{Tr} A_+(\mathbf{k}, E) > 0\big) .$$
Then:
1. $A_+(\mathbf{k}, E)$ is the unique nearest positive-semidefinite matrix to $A(\mathbf{k}, E)$ in the Frobenius norm $\|\cdot\|_F$, and strictly reduces the Frobenius distance to the exact physical spectral matrix $A^{\text{exact}}(\mathbf{k}, E) \succeq 0$:
   $$\big\| A_+(\mathbf{k}, E) - A^{\text{exact}}(\mathbf{k}, E) \big\|_F^2 \le \big\| A(\mathbf{k}, E) - A^{\text{exact}}(\mathbf{k}, E) \big\|_F^2 - \sum_{\lambda_\alpha < 0} \lambda_\alpha(\mathbf{k}, E)^2 .$$
2. *Effect on sum rules:* Unscaled $A_+(\mathbf{k}, E)$ increases the trace at $(\mathbf{k}, E)$ by $w_{\text{neg}}(\mathbf{k}, E) \ge 0$, raising the BZ-integrated DOS at energy $E$ by $\langle w_{\text{neg}}(\mathbf{k}, E) \rangle_{\mathbf{k}}$. By contrast, whenever the total trace $\operatorname{Tr} A(\mathbf{k}, E) \ge 0$ across the BZ, the trace-preserving projection $A_+^{\text{tr}}(\mathbf{k}, E) \succeq 0$ **preserves the exact total DOS sum rule** $\Omega_{\text{BZ}}^{-1} \int_{\text{BZ}} d^3k \, \operatorname{Tr} A_+^{\text{tr}}(\mathbf{k}, E) = \operatorname{DOS}(E)$ and the total electron count, while guaranteeing non-negative orbital and sublattice weights everywhere.

*Proof.* (1) The set $\mathcal{S}_+$ of Hermitian positive-semidefinite matrices is a closed convex cone in the real Hilbert space of Hermitian matrices with inner product $\langle X, Y \rangle_F = \operatorname{Tr}(X Y)$. By the Hoffman-Wielandt theorem, $A_+$ is the orthogonal projection of $A$ onto $\mathcal{S}_+$. Since $A^{\text{exact}} \in \mathcal{S}_+$, the convex projection inequality $\langle A - A_+, A^{\text{exact}} - A_+ \rangle_F \le 0$ gives $\|A - A^{\text{exact}}\|_F^2 = \|A - A_+\|_F^2 + \|A_+ - A^{\text{exact}}\|_F^2 - 2\langle A - A_+, A^{\text{exact}} - A_+ \rangle_F \ge \|A - A_+\|_F^2 + \|A_+ - A^{\text{exact}}\|_F^2$.
(2) By construction, $\operatorname{Tr} A_+^{\text{tr}}(\mathbf{k}, E) = \operatorname{Tr} A(\mathbf{k}, E)$ whenever $\operatorname{Tr} A(\mathbf{k}, E) \ge 0$. Since all off-site terms $\mathbf{R} \neq \mathbf{0}$ have zero Brillouin-zone integral, $\Omega_{\text{BZ}}^{-1} \int_{\text{BZ}} d^3k \, \operatorname{Tr} A(\mathbf{k}, E) = \operatorname{DOS}(E)$. $\square$

---

### Q7. Short-Range Order (SRO) in $A(\mathbf{k}, E)$

**Result Q7.1 (Homogenised centre-ego definition of $A(\mathbf{k}, E)$ with first-shell SRO): Proven and checked numerically.**
In the Mookerjee-Prasad first-shell SRO model of `src/asr_sro.f90` and `src/asr_class_scf.f90`, freezing site $0\kappa$ to species $c$ defines a probability measure $\mathbb{P}^{(0\kappa)}(\mathcal{C} \mid s_0 = c)$ in which first-shell neighbours $R_1 \in \text{Shell}_1(0\kappa) \cap \kappa'$ carry the conditional concentration $X_s^{(c)}(\kappa, \kappa') = x_s(\kappa')(1 - \alpha_{cs}(\kappa, \kappa'))$ (class $N_{\text{cls}} + \kappa'$) and all outer sites $R \notin \{0\kappa\} \cup \text{Shell}_1(0\kappa)$ carry the unconditional concentration $x_s(\kappa')$.
1. *Single-site homogeneity (Proven):* By the marginal condition $\sum_{c \in \kappa} x_c(\kappa) \alpha_{cs}(\kappa, \kappa') = 0$, the marginal occupation probability of every site $R \in \kappa'$ under $\mathbb{P}^{(0\kappa)}(\mathcal{C}) = \sum_{c \in \kappa} x_c(\kappa) \mathbb{P}^{(0\kappa)}(\mathcal{C} \mid s_0 = c)$ is $x_s(\kappa')$.
2. *Translation-invariant propagator definition:* Define the translation-invariant two-point correlation function for displacement $\mathbf{R} = \mathbf{T} + \boldsymbol{\tau}_{\kappa'} - \boldsymbol{\tau}_\kappa$ as the symmetrised centre-ego average
   $$\langle \mathcal{G}_{0\kappa L, R L'}(z) \rangle^{\text{SRO}} = \frac{1}{2} \left( \langle \mathcal{G}_{0\kappa L, R L'}(z) \rangle_{\mathbb{P}^{(0\kappa)}} + \langle \mathcal{G}_{0\kappa' L', (-R) L}(z) \rangle_{\mathbb{P}^{(0\kappa')}} \right) ,$$
   where in Form 2B, for origin $0\kappa$ and central species $c$, the operator is the SRO operator `cop(kap)%op` of `src/asr_class_scf.f90` and the target vector $|\widetilde{R L'}^{(c)}\rangle = \sum_{k=0}^{n_{\kappa'}-1} \tau^{(c)}_{R L', k} |R L', \{R: k\}\rangle$ uses the local basis $\phi_k^{(R)}$ of site $R$'s class in `cop(kap)%op`:
   $$\tau^{(c)}_{R L', k} = \sum_{s \in \kappa'} p_s^{(c)}(R) \, T_{R L'}(s) \, \phi_k^{(R)}(s), \qquad p_s^{(c)}(R) = \begin{cases} X_s^{(c)}(\kappa, \kappa') & \text{if } R \in \text{Shell}_1(0\kappa), \\ x_s(\kappa') & \text{if } R \notin \text{Shell}_1(0\kappa) . \end{cases}$$
   (In Form 2A, the pair weight is $P_{0R}^{\text{MP}}(c, s) = x_c(\kappa) p_s^{(c)}(R)$ while the remaining $z_1(\kappa, \kappa') - 1$ first-shell sites keep concentration $X_{s'}^{(c)}(\kappa, \kappa')$.)
3. *Reduction to Q1 without SRO (Proven):* At $\alpha_{cs}(\kappa, \kappa') = 0$, $X_s^{(c)}(\kappa, \kappa') = x_s(\kappa')$ for all $c, s$, so $\mathbb{P}^{(0\kappa)}$ becomes the translation-invariant product measure $\prod_R x_{s_R}$ and $\langle \mathcal{G}_{0\kappa L, R L'}(z) \rangle^{\text{SRO}}$ reduces identically to Q1.

*Proof.* (1) For $R \in \text{Shell}_1(0\kappa)$, $\sum_{c \in \kappa} x_c(\kappa) X_s^{(c)}(\kappa, \kappa') = x_s(\kappa') \sum_{c \in \kappa} x_c(\kappa)(1 - \alpha_{cs}(\kappa, \kappa')) = x_s(\kappa')$. (2) Symmetrising between $0\kappa \to R$ and $0\kappa' \to -R$ guarantees $\langle \mathcal{G}_{0\kappa L, R L'}(z) \rangle^{\text{SRO}} = \langle \mathcal{G}_{0\kappa' L', (-R) L}(z) \rangle^{\text{SRO}}$, so Result Q1.1 holds and $A^{\text{SRO}}(\mathbf{k}, E)$ is strictly Hermitian at every $\mathbf{k}$. (3) Setting $\alpha = 0$ makes $p_s^{(c)}(R) = x_s(\kappa')$ on all shells. $\square$

*Numerical check (Test 1 of Section 4):* On the 6-site ternary ring with Warren-Cowley parameters $\alpha = (-0.20, 0.10, 0.15)$, Form 2B with $\tau^{(c)}_{R, k}$ on the conditional shell agrees with exact enumeration of the Mookerjee-Prasad measure to $1.21 \times 10^{-14}$.

*What the first-shell SRO model can and cannot represent:*
* **Can represent:** (i) Self-consistent charge-transfer and Madelung/screened potential shifts (`asr_screening`) that move the band centres $C_p, E_{\nu,p}$ of each component; (ii) modification of nearest-neighbour bonding/antibonding hybridization and disorder lifetime $\Gamma_{\text{dis}}(\mathbf{k}, E)$ via the correlated first shell $R_1 \in \text{Shell}_1(0\kappa)$; (iii) sharpening of majority-species bands under clustering ($\alpha > 0$) and pseudo-gap opening under ordering ($\alpha < 0$).
* **Cannot represent:** (i) Superlattice Bragg peaks, zone folding, or secondary bands at ordering wavevectors $\mathbf{Q}$ (because $\alpha^{(n)}_{cs} = 0$ for all shells $n \ge 2$; representing long-range order requires sublattice classes in `decoration_t`); (ii) correlations between two outer sites $R_1, R_2 \neq 0\kappa$ that are nearest neighbours of each other; (iii) global Bochner positivity at extreme $|\alpha|$, because $\mathbb{P}^{(0\kappa)}$ is an ego-centred measure rather than a globally translation-invariant random field (cured by the projection $A_+^{\text{tr}}(\mathbf{k}, E)$ of Result Q6.2).

---

## Deliverable 2: Recommended Route, Pseudocode, and New Code Operations

### Recommended Route

Use **Form 2B + Q3b (Single-centre operator with dressed target projection)** as the production route for $A(\mathbf{k}, E)$, combined with the **3D Fejér ball window $w_{\text{ball}}(R)$** (Result Q4.2) and the **trace-preserving positive projection $A_+^{\text{tr}}(\mathbf{k}, E)$** (Result Q6.2).
* It reuses the exact same $D_{4h}$-reduced operator `cop(kap)%op` already built in `src/asr_class_scf.f90`.
* It requires only **5 start orbitals** $L_0 \in \{1, 3, 5, 7, 9\}$ per component and energy (all five are 1D under $D_{4h}$), reconstructing the remaining 4 start orbitals $\{2, 4, 6, 8\}$ by $C_3$ cyclic symmetry on cubic sublattices.
* It extracts all shells $|\mathbf{R}| \le R_{\max}$ and all $9 \times 9$ orbital blocks from the **same chains** with $< 2\%$ overhead per step.

### Pseudocode

```text
ALGORITHM: Configuration-Averaged Bloch Spectral Function A(k, E) (Route edc)

INPUT:
  - Converged class_scf state: sublattices kap = 1..ncls, components p = (kap, c) with weights x_c(kap),
    potential parameters pps(p, spin), screening alpha(0:2).
  - Existing D4h-reduced operators cop(kap)%op (built at maxlevel = L_lev, maxrank = r_max).
  - Shell cutoff R_max (with d_graph(0, R) <= 2*(L_lev - 1)), window type ('fejer' or 'sharp').
  - Energy grid E_grid(1:N_E), broadening eta > 0, k-point path k_pts(3, 1:N_k).

STEP 1: Precompute shell target table (once before the energy loop)
  For each sublattice kap = 1..ncls:
    1.1 Find all lattice sites R in the supercell with |d(0_kap, R)| <= R_max.
        Record displacement d_R = R - R_{0_kap} and partner sublattice kap' = cls_base(R).
    1.2 For each site R and each local configuration state m_st = 0 .. nsc(kap') - 1:
        Let cfg = empty if m_st == 0 else {R: m_st}.
        Call canon_node(perm, nop, nsm, R, cfg, cs, ccf, osz, g_link).
        Find reduced node index r_node = find_state(st, cs, ccf) in cop(kap)%op.
        Store (r_node, g_link, 1 / sqrt(osz)) in table target_link(m_st, R, kap).

STEP 2: Run chains and project onto shells (parallel over spin, component p = (kap, c), energy block, L0)
  Let iorb5 = [1, 3, 5, 7, 9]  ! s, p_z, d_xy, d_{3z^2-r^2}, d_{x^2-y^2} (all 1D under D4h)
  For each spin, component p = (kap, c), energy z = E + i*eta, and o = 1..5 (L0 = iorb5(o)):
    2.1 Set operator parameters:
        call edc_setup(mdl, pps(:, spin), alpha, p, ed)
        call augcls_set_params(op, mdl)
        call augcls_freeze_site(op, 1, p)
        if (sro) set first-shell class concentrations to sro_cond(c, :, kap, :)
        call augcls_set_start(op, L0, cop(kap)%d)
    2.2 Precompute dressed scaling weights tau(L', m_st, class_of_R) for L' = 1..9:
        T_s(L') = sqrt(ed%dc(L') / ed%del(L', s))
        tau(L', m_st, cl) = sum_{s in cl} conc(s, cl) * T_s(L') * phi(m_st, s, cl)
    2.3 Choose start vector u_1 = e_{9*(1-1) + L0} (normalised: u_1^T u_1 = 1).
    2.4 Run m = recursion_steps of the complex-symmetric Lanczos chain at energy z:
        For step k = 1..m:
          a. For each shell site R and orbital L' = 1..9:
             Read u_k at target_link(m_st, R, kap) via at_link_c:
               amp(L', m_st) = (chi(g) * sg(L', g) / sq(r_node)) * u_k(9*(r_node-1) + pi(L', g))
             Accumulate projection:
               c_proj(k, L', R) = sum_{m_st = 0}^{nsc-1} tau(L', m_st, ncl(R)) * amp(L', m_st)
          b. Apply H^{(c)}(z) to u_k, compute a(k) = u_k^T w, orthogonalise against u_k and u_{k-1},
             call augcls_project_c(op, w), compute b(k) = sqrt(w^T w), advance u_{k+1} = w / b(k).
    2.5 Close the chain with the 'hull' tail (cf_species) and read projections:
        Compute continued-fraction factors prod(1:m) from (a(1:m), b(1:m), z).
        For each shell site R and orbital L' = 1..9:
          G_shell(L0, L', R, p, z) = sum_{k=1}^m c_proj(k, L', R) * prod(k)

STEP 3: Reconstruct remaining 4 start orbitals {p_y(2), p_x(4), d_yz(6), d_xz(8)} by C3 symmetry
  Apply the 3-fold cyclic permutation C3: (x,y,z) -> (y,z,x) and C3^2 to the shell displacements d_R
  and orbital representation matrix D(C3) on the p_z(3) and d_xy(5) rows to fill L0 in {2, 4, 6, 8}.
  Average over central species c on each sublattice kap:
    G_avg(kap, L, R, L', z) = sum_{c in kap} x_c(kap) * G_shell(L, L', R, (kap, c), z)
  If SRO is active, symmetrise:
    G_avg(kap, L, R, L', z) = 0.5 * (G_avg(kap, L, R, L', z) + G_avg(kap', L', -R, L, z))

STEP 4: Fourier sum and positive projection
  For each energy E and wavevector k:
    4.1 Form Hermitian 9*ncls x 9*ncls matrix A(k, E):
        A_{kap L, kap' L'}(k, E) = sum_{R in kap'} w_win(|d_R|) * exp(i k . d_R) * (-1/pi) Im G_avg(kap, L, R, L', E + i*eta)
        Hermitise: A(k, E) = 0.5 * (A(k, E) + A(k, E)^H).
    4.2 Diagonalise A(k, E) = V diag(lambda) V^H (zheev).
        Record w_neg(k, E) = sum_{lambda_a < 0} |lambda_a|.
        Form A_+^{tr}(k, E) by clipping lambda_a -> max(lambda_a, 0) and rescaling to preserve Tr A(k, E).
```

### List of New Operations Needed in the Fortran Code

1. **`augcls_shell_links` (in `src/asr_augcls.f90`):** During `build_augcls_operator` (before `st` is deallocated at the end of the routine), look up `find_state(st, cs, ccf)` for each shell site $R$ ($|\mathbf{R}| \le R_{\max}$) and each single-site flip state $m_{\text{st}} \in \{0, \dots, n_{\text{sc}}-1\}$, storing the canonical node index `r_node` and operation `g_link`.
2. **`edc_green_shells` / `edc_green_shells_blk` (in `src/asr_edc.f90`):** Extend `edc_green` and `edc_green_blk` to accept the shell link table and accumulate the sparse dot products `c_proj(k, L', R)` at the ~100 to 300 shell nodes during each Lanczos step, then multiply by `prod(1:m)` in `cf_species`.
3. **`ed3_shells` / `ed3_shells_blk` (in `src/asr_ed3.f90`):** Same extension for the two chains from $|a_{0L}^{(c)}\rangle$ and $|v_{0L}^{(c)}\rangle$, projecting onto $|a_{RL'}^{(c,s)}\rangle$ and $|v_{RL'}^{(c)}\rangle$.
4. **`c3_orbital_complete` (in `src/asr_symmetry.f90`):** Apply the $C_3$ permutation $(x, y, z) \to (y, z, x)$ to map the $p_z$ and $d_{xy}$ shell blocks onto $p_x, p_y$ and $d_{yz}, d_{xz}$ (and include $L_0 = 9$ in the start-orbital list `[1, 3, 5, 7, 9]`).
5. **`spectral_fourier_and_project` (new module `src/asr_spectral.f90`):** Evaluate the windowed Bloch sum over shells $|\mathbf{R}| \le R_{\max}$, call `zheev` on the $9 N_{\text{cls}} \times 9 N_{\text{cls}}$ matrix $A(\mathbf{k}, E)$, compute $A_+^{\text{tr}}(\mathbf{k}, E)$, and write $A(\mathbf{k}, E)$ along the requested k-path.

---

## Deliverable 3: Error Control and Diagnostic Prints

### What Converges with Each Parameter

1. **`steps` ($m = \text{recursion\_steps}$):** Controls the Galerkin residual $\|r_c\|_{M^{-1}}$ on the fixed node set $V_{L_{\text{lev}}}$ (Result S3 of `docs/DERIVATION_levels_steps.md`). Because off-site Krylov projection is first order in $\|r_c\|$ (Result Q3.1), $A(\mathbf{k}, E)$ converges as $|\rho_m| \sim (1 + \eta / W)^{-m}$. For $\eta = 5\text{ to }10\text{ mRy}$ and half-bandwidth $W \approx 0.3\text{ Ry}$, $m = 48\text{ to }64$ steps reduces the chain termination scatter below $2\%$.
2. **`levels` ($L_{\text{lev}}$):** Controls the spatial radius $2 L_{\text{lev}}$ (in graph hops) of the augmented-space node set $V_{L_{\text{lev}}}$ and the boundary leak $\|r_\ell\|_{M^{-1}}$. All shells used in the Fourier sum must satisfy $d_{\text{graph}}(0\kappa, R) \le 2(L_{\text{lev}} - 1)$.
3. **$R_{\max}$:** Controls the spatial Fourier resolution $\Delta k \sim \pi / R_{\max}$. The error decays exponentially as $e^{-\kappa_{\text{tot}}(E, \eta) R_{\max}}$ where $\kappa_{\text{tot}} = (2\ell_{\text{mfp}})^{-1} + \eta / v_F$ (Result Q4.1).
4. **$\eta$ (`dos_eta`):** Sets the Lorentzian energy resolution and provides the minimum spatial decay rate $\gamma(\eta) \approx \eta / v_F$. Reducing $\eta$ sharpens quasi-particle peaks in weak-scattering regions ($s/p$ bands) and requires a proportionally larger $R_{\max} \gtrsim 3 v_F / \eta$ and $m \gtrsim 2 W / \eta$.

### Checks the Code Should Print in `asr.out`

1. **Shell decay profile:** For each coordination shell $n$ up to $R_{\max}$ (distance $R_n$, multiplicity $z_n$), print $\max_E \|\operatorname{Im}\langle \mathcal{G}_{0, R_n}(E + i\eta)\rangle\|_F$ and the ratio of the outermost shell $R_{\max}$ to the on-site element $R = 0$ (target: $< 10^{-2}$).
2. **Chain step-spread check (at zero extra chain cost):** Re-evaluate `cf_species` at $m' = m - 4$ and $m$ from the same chain coefficients $a(1:m), b(1:m)$, and print $\max_{\mathbf{k}, E} |A_m(\mathbf{k}, E) - A_{m-4}(\mathbf{k}, E)| / \max_{\mathbf{k}, E} A_m(\mathbf{k}, E)$.
3. **Symmetry reciprocity defect:** For shells where both $(L, R, L')$ and $(L', -R, L)$ are computed independently (e.g. $s\text{--}d$ or $p\text{--}d$ off-diagonal blocks), print $\max |\langle \mathcal{G}_{0L, RL'} \rangle - \langle \mathcal{G}_{0L', (-R)L} \rangle|$, which is zero when the chains are converged and measures the first-order projection residual $O(\|r_m\|)$.
4. **Positivity and clipped negative weight:** Print the fraction of $(\mathbf{k}, E)$ points where $\lambda_{\min}(A(\mathbf{k}, E)) < 0$, the worst negative eigenvalue $\min_{\mathbf{k}, E} \lambda_{\min}$, and the integrated negative weight ratio $\sum_{\mathbf{k}, E} w_{\text{neg}}(\mathbf{k}, E) / \sum_{\mathbf{k}, E} \operatorname{Tr} A_+(\mathbf{k}, E)$.
5. **Sum rule check:** Print the integrated spectral weight $\Delta E \sum_E A(\mathbf{k}, E)$ at the high-symmetry endpoints ($\Gamma, \text{X}, \text{L}, \text{W}$) against $9 N_{\text{cls}}$.

---

## Deliverable 4: Exact Limits (Q8) and NumPy Verification Script with Output

### Exact Limits and How to Test Each

1. **Ordered limit ($x_1(\kappa) = 1$, $x_{s>1}(\kappa) = 0$):** All flip weights $\tau_{R, k \ge 1}$ and $W_{l, k}$ ($k \ge 1$) vanish identically. Test by comparing $A(\mathbf{k}, E)$ at identical $\alpha$ (`std_alpha = true`) and cluster radius `rcl` against `ham_k_cell` (`src/asr_lmto_cell.f90`) for `edc` and `bands_core` (`src/asr_asa.f90`) for `ed3`; agreement must be at rounding/chain-closure level ($< 10^{-13}$).
2. **All species equal ($p_1 = p_2 = \dots = p_n$ at arbitrary concentrations $x_s$):** Because $\phi_0(s) = 1$ and $\sum_s x_s \phi_k(s) = \delta_{k,0}$, every random parameter $X_s = X_1$ has $W_{lk} = X_1 \delta_{l,k}$ and $\tau_{R, k} = \delta_{k,0}$. All fluctuation sectors decouple identically at every rank and level, and $A(\mathbf{k}, E)$ must agree with the ordered crystal of species 1 to machine precision ($< 10^{-14}$).
3. **Dilute limit ($x_1 = 1 - \epsilon$, $x_{s>1} = O(\epsilon)$ as $\epsilon \to 0$):** Verify that $\|A(\mathbf{k}, E; \epsilon) - A_{\text{host}}(\mathbf{k}, E)\| = O(\epsilon)$ and that the finite difference $(G(\mathbf{k}, z; \epsilon) - G_{\text{host}}(\mathbf{k}, z))/\epsilon$ converges to the exact single-impurity T-matrix $G_{\text{host}}(\mathbf{k}, z) \, t_{\text{imp}}(z) \, G_{\text{host}}(\mathbf{k}, z)$.
4. **DOS sum rule (k-integral per energy) and state count (E-integral per k):** Verify that averaging $A_{\kappa L, \kappa L}(\mathbf{k}, E)$ over the full supercell $\mathbf{k}$-mesh reproduces the on-site ASR $\operatorname{DOS}_{\kappa L}(E)$ to machine precision ($< 10^{-14}$), and that trapezoidal integration of $A(\mathbf{k}, E)$ over $E \in [-2, +2]\text{ Ry}$ plus the Cauchy tail $9 N_{\text{cls}} \frac{2\eta}{\pi E_{\max}}$ equals $9 N_{\text{cls}}$ to $O((\eta/E_{\max})^2)$.
5. **Non-negativity:** Verify $\lambda_{\min}(A_+^{\text{tr}}(\mathbf{k}, E)) \ge 0$ everywhere and $\lambda_{\min}(A^{w_{\text{Fej\'er}}}(\mathbf{k}, E)) \ge 0$ when chains are converged.

### NumPy Verification Script Output (`tools/check_spectral_asr.py`)

```text
==============================================================================
TEST 1: Exact Off-Site Scaling (Q2) on 6-Site 3-Species Ring (729 configs)
==============================================================================
  Max |G_0R(Form 2A, full rank, +1 power T_R) - G_0R(enum)| = 1.168e-14
  Max |G_0R(Form 2A, full rank, -1 power T_R) - G_0R(enum)| = 1.456e-02  <-- T_R^-1 FAILS
  Max |G_0R(Form 2B, full rank, single-centre) - G_0R(enum)| = 1.168e-14
  Max |G_0R(ed3 full ASA, full rank)           - G_0R(enum)| = 3.820e-15
  Rank-1 max error |G_0R - G_enum|: Form 2A (pair) = 5.727e-03, Form 2B (single) = 9.648e-03
  Rank-2 max error |G_0R - G_enum|: Form 2A (pair) = 1.662e-04, Form 2B (single) = 4.234e-04
  Max |G_0R(SRO Form 2B, full rank)           - G_0R(MP enum)| = 1.208e-14

==============================================================================
TEST 2: Off-Diagonal Extraction Accuracy Orders (Q3: Polarisation vs Projection)
==============================================================================
  steps m |  ||r_m||_2   | err Projection (Q3b) | err Polarisation (Q3a) | err Petrov-Galerkin
      4   |  3.177e-01  |      5.663e-01       |       1.273e-01        |      1.572e-01
      6   |  6.604e-02  |      7.970e-02       |       2.537e-02        |      1.928e-02
      8   |  1.588e-01  |      6.200e-02       |       2.159e-02        |      1.842e-02
     10   |  7.243e-03  |      1.702e-03       |       7.340e-04        |      2.438e-04
     12   |  1.201e-02  |      3.292e-03       |       1.892e-03        |      6.058e-04
     16   |  3.030e-05  |      1.249e-06       |       3.952e-05        |      9.192e-08

==============================================================================
TEST 3: Block Recursion (Q3c) Second-Order Galerkin Error Identity
==============================================================================
  block steps m= 2 (2m= 4 matvecs): ||R_m||=1.715e+00, |G_0R - G_blk|=2.886e+00, identity defect=1.073e-14
  block steps m= 3 (2m= 6 matvecs): ||R_m||=8.012e-01, |G_0R - G_blk|=4.410e-01, identity defect=3.925e-16
  block steps m= 4 (2m= 8 matvecs): ||R_m||=4.805e-01, |G_0R - G_blk|=1.413e-01, identity defect=3.002e-16
  block steps m= 6 (2m=12 matvecs): ||R_m||=1.455e-01, |G_0R - G_blk|=2.365e-02, identity defect=4.707e-16

==============================================================================
TEST 4: 3D Simple Cubic 3-Species Alloy (Fe0.66 Ni0.14 Cr0.20 model, Q4-Q8)
==============================================================================
  3D SC augmented space (L=4, levels=3, rank=1): 64 cluster sites, 4956 reachable aug nodes
  Limit 1 (Ordered vs k-space): |G_00(ASR) - G_00(k)| = 6.894e-15, |G_100(ASR) - G_100(k)| = 1.299e-15
  Limit 2 (Equal species at x=(0.66,0.14,0.20) vs Ordered): max_R |G_0R(eq) - G_0R(ord)| = 1.986e-15
  Limit 3 (Dilute eps=1.0e-04): |G_00(eps) - G_00(0)| = 2.459e-03 = 24.590 * eps (linear O(eps) confirmed)
  Limit 4 (BZ k-integral sum rule vs on-site DOS): max_E |<A(k,E)>_BZ - DOS(E)| = 0.000e+00
  Limit 5 (Non-negativity on grid): min A_sharp = 0.0337, min A_fejer = 0.0509
  3D Alloy A(k,E) ASR vs 64-site cell, 120-sample Bloch average: L2 rel diff (sharp) = 0.071, (Fejer, cut at 3.96 > largest distance 3.46) = 0.364

  Sample A(k, E) values at E = -0.15 Ry (inside d-band) and E = +0.05 Ry:
    k = Gamma: E=-0.15 Ry -> ASR=1.939 (Fejer=3.286, SC=1.705);  E=+0.05 Ry -> ASR=0.221 (Fejer=0.479, SC=0.235)
    k = Delta: E=-0.15 Ry -> ASR=5.046 (Fejer=3.296, SC=5.207);  E=+0.05 Ry -> ASR=0.330 (Fejer=0.605, SC=0.355)
    k = X    : E=-0.15 Ry -> ASR=2.912 (Fejer=2.517, SC=3.037);  E=+0.05 Ry -> ASR=0.556 (Fejer=0.777, SC=0.578)
    k = M    : E=-0.15 Ry -> ASR=0.557 (Fejer=0.907, SC=0.555);  E=+0.05 Ry -> ASR=1.978 (Fejer=1.730, SC=1.891)
    k = R    : E=-0.15 Ry -> ASR=0.257 (Fejer=0.537, SC=0.218);  E=+0.05 Ry -> ASR=2.858 (Fejer=2.708, SC=3.036)

==============================================================================
TEST 5: Energy-Independent Approximations vs Energy-Dependent ASR (Q5)
==============================================================================
  k=0 (Gamma) : A_exact=0.5964 | ED rank-1 (pair)=0.5958 (err 0.1%) | ED rank-1 (single)=0.5955 (err 0.2%) | EI (<d>,<o>)=0.6366 (err 6.7%)
  k=pi (X)    : A_exact=1.2394 | ED rank-1 (pair)=1.2378 (err 0.1%) | ED rank-1 (single)=1.2375 (err 0.2%) | EI (<d>,<o>)=1.2463 (err 0.6%)
```

---

## Deliverable 5: Cost Estimate for $\text{Fe}_{0.66}\text{Ni}_{0.14}\text{Cr}_{0.20}$ (Q9)

### Baseline Calibration from the Prompt

* System: fcc $\text{Fe}_{0.66}\text{Ni}_{0.14}\text{Cr}_{0.20}$ at `levels = 3`, `rank = 1`, $N_{\text{red}} = 15,827$ reduced nodes in $D_{4h}$, $N_{\text{sp}} = 3$ species, 9 orbitals.
* Reference SCF iteration at $m = 24$ steps: 2 to 4 contour passes of $77 \text{ energies} \times 3 \text{ species} \times 4 \text{ start orbitals} = 924\text{ chains per pass}$, giving $1,848\text{ to }3,696\text{ chains}$ in $5\text{ to }7\text{ min}$ on 8 cores ($\approx 0.13\text{ s}$ wall time per 24-step chain on 8 cores, or $\approx 0.10\text{ s}$ with `energy_block = 8`).
* Target panel: $N_E = 300$ energies, $N_k = 100$ $\mathbf{k}$-points along 4 high-symmetry lines, $\eta = 5\text{ to }10\text{ mRy}$, shells up to $R_{\max} = 2.0\,a$ (6 fcc shells, 86 lattice sites).

### Cost of One $A(\mathbf{k}, E)$ Panel (per Spin Channel)

| Route | Operators built | Reduced nodes per operator | Chains per $(E, c)$ | Total chains ($N_E = 300, N_{\text{sp}} = 3$) | Wall time on 8 cores ($m = 24 / 48 / 96$ steps) |
|---|---|---|---|---|---|
| **Recommended: Form 2B + Q3b (Single-centre projection, `edc`)** | $1$ (reuse SCF operator) | $15,827$ | $5$ ($s, p_z, d_{xy}, d_{3z^2-r^2}, d_{x^2-y^2}$) | **$4,500$** | **$9\text{ min}$ / $18\text{ min}$ / $36\text{ min}$** ($7 / 14 / 27\text{ min}$ with `energy_block = 8`) |
| **Recommended: Form 2B + Q3b (`ed3` full ASA)** | $1$ (reuse SCF operator) | $15,827$ | $10$ ($2 \times 5$) | **$9,000$** | **$18\text{ min}$ / $36\text{ min}$ / $72\text{ min}$** |
| **Alternative: Form 2A + Q3b (Frozen pair, 6 shells, `edc`)** | $6$ (one per shell orbit) | $\approx 55,000$ (mean over $G_{\text{pair}}$) | $9 \times N_{\text{sh}} N_{\text{sp}} = 162$ | **$145,800$** | **$17\text{ h}$ / $34\text{ h}$ / $68\text{ h}$** |

* **Projection and Fourier sum overhead:** Reading 86 shell sites $\times 3$ local states $= 258$ nodes per step touches $1.6\%$ of the $15,827$ nodes and adds $< 2\%$ to the Lanczos step time. The Fourier sum over $300 \times 100 \times 86 \times 81 \approx 2.1 \times 10^8$ multiply-adds takes **$< 0.05\text{ s}$** total and does not scale with the number of chains.

---

## Deliverable 6: Risks and Failure Modes

1. **First-order residual and F4 non-causality at short chain lengths ($m = 24$) and small $\eta$ (`edc` / `ed3`):**
   * *Mechanism:* Off-site projections $\widetilde{G}_{0R}(z)$ converge as $O(\|r_m\|)$ (first order in the Galerkin residual, Result Q3.1) rather than $O(\|r_m\|^2)$, and begin picking up signal only at step $k \ge \lceil d_{\text{graph}}(0, R)/2 \rceil + 1$. At $m = 24$ and $\eta = 5\text{ mRy}$, the `hull` tail and unconverged Ritz values can produce negative eigenvalues in $A(\mathbf{k}, E)$ near band edges or weak-scattering $s/p$ pockets.
   * *Mitigation:* Run $m = 48\text{ to }64$ steps for spectral panels, apply the Fejér window $w_{\text{ball}}(R)$, apply the trace-preserving positive projection $A_+^{\text{tr}}(\mathbf{k}, E)$, and monitor the printed negative weight $w_{\text{neg}}(\mathbf{k}, E)$ and step spread $|A_m - A_{m-4}|$.
2. **Sharp $s/p$ bands with long mean free path ($\ell_{\text{mfp}} > R_{\max}$):**
   * *Mechanism:* In transition-metal alloys, $d$ states scatter strongly ($\ell_{\text{mfp}} \approx 1\text{ to }2\,a$) and decay within 3 to 4 shells, whereas $s/p$ states below $d$-band bottom have weak disorder scattering ($\ell_{\text{mfp}} \gg 4\,a$). At $\eta = 5\text{ mRy}$, $\langle G_{0R} \rangle$ in the $s/p$ band has not decayed by $R_{\max} \approx 2\text{ to }3\,a$, causing ring oscillations under sharp Fourier truncation.
   * *Mitigation:* Use the Fejér ball window $w_{\text{ball}}(R)$ (which eliminates Gibbs negative sidelobes at the cost of broadening the $s/p$ peak width in $\mathbf{k}$ by $\Delta k \sim 2\pi / R_{\max}$) or increase `levels` to 4 for $s/p$ dispersion checks.
3. **SRO model limitations at large $|\alpha|$ and missing outer-shell correlations:**
   * *Mechanism:* The Mookerjee-Prasad model correlates only the central site $0\kappa$ with its first shell $\text{Shell}_1(0\kappa)$ and sets all outer shells $n \ge 2$ to uncorrelated $x_s$. It cannot produce superlattice zone folding or Bragg gaps (which require multi-shell correlations $\alpha^{(n)}$ for $n \ge 2$ or explicit sublattices in `decoration_t`), and its centre-ego measure $\mathbb{P}^{(0\kappa)}$ breaks global translation invariance on bonds between shell-1 and shell-2 sites.
4. **Rank-1 truncation of simultaneous flip at $R$ in Form 2B vs Form 2A:**
   * *Mechanism:* In Form 2B, the $m_{\text{st}} \ge 1$ component of $|\widetilde{R L'}^{(c)}\rangle = \widehat{T}_R |R L', \emptyset\rangle$ occupies the 1-flip state $\{R: m_{\text{st}}\}$. At `rank = 1`, paths reaching $\{R: m_{\text{st}}\}$ cannot simultaneously hold a flip at an intermediate site $R'' \notin \{0, R\}$ without first unflipping $R''$, whereas Form 2A freezes $R$ and retains 1 free flip on intermediate sites $R''$.
   * *Size and check:* In Test 1 ($D_{\max}/D_{\min} = 1.96$), the difference between Form 2A and Form 2B at `rank = 1` is $0.4\%$ of $G_{0R}$ ($0.1\%$ in $A(k, E)$), and drops by a factor of $25$ at `rank = 2` (which is affordable at `levels = 2` with $22,958$ nodes). Form 2A on the first shell can be run at a single energy as an exact internal diagnostic check of Form 2B.

---

## Review (2026-11-17, Claude Opus 5.5)

The script was rerun on the laptop (13 s).  It reproduces every number above.  The output labels of Test 4 and Test 5 were
corrected on 2026-11-17; the numbers did not change.

### What the review accepts
1. Q1: the definitions, the positivity proof, the sum rules and the ordered limit.
2. Q2: the off-site scaling with $T_R$ to the power $+1$.  The prompt's guess $T_R^{-1}$ was wrong.  The algebra is
   $G^{(2)} = T\,(z - H^{(c)})^{-1}\,T$ with $T = 1$ on the frozen centre.  Test 1 confirms it against exact enumeration of
   all 729 configurations of the ring to $10^{-14}$, for `edc`, for `ed3` and for the SRO model.  The formulas in the script
   match the code's Hamiltonian (`ed3` secular matrix, sphere-charge weighted $G$).
3. The recommended route (Form 2B with Q3b): the chains of the single-centre operator already determine the amplitude at every
   site.  So one chain per centre species, start orbital and energy gives every shell at once.  The cost is about one DOS pass:
   10 to 30 min per spin for \ce{Fe66Ni14Cr20} at levels 3.
4. Q3: the accuracy orders, and the block identity of Test 3 ($10^{-16}$).

### What the review does not adopt
1. The Fejér window as the default.  Test 4 is the only evidence, and there the window makes $A(\mathbf{k},E)$ five times
   worse (L2 difference 0.36 against 0.07 for the sharp cutoff), while the sharp cutoff is already non-negative.  The test is
   too small to settle the question: the window is cut at 3.96 on a torus whose largest distance is 3.46.  Use the sharp
   cutoff with $R_{\max}$ chosen from the printed shell decay.  Keep the window as a fallback when negative values appear.
2. The trace-preserving clip $A_+^{\mathrm{tr}}$ as a default.  It hides the chain error and changes the orbital weights.
   Print the negative weight instead.  Clip only for display, and say so in the figure.

### What the review narrows
1. Result Q5.1 holds for a frozen centre only.  A scaling with a reference $D_{\mathrm{ref}}$ and no frozen site keeps
   $E - H(E) = T\,(E - H^{(2)})\,T$ with a non-random $S'$, and that operator is translation invariant in augmented space.
   So an energy-dependent k-space recursion with dressed Bloch start vectors is possible in principle.  It is not needed:
   the real-space route is cheaper and keeps the point-group reduction.
2. Q4.1 item 2 (the decay law) is a physical estimate, not a numerical check.  The 3 to 6 % of Q5.2 is an estimate, not a
   proof.  Test 5 does not test an energy-independent recursion.  It enumerates the alloy exactly with $d$ and $o$ replaced
   by their mean values.
3. Q1.2 item 4 states $9 N_{\mathrm{cls}}$ and its proof ends with $9 N_{\mathrm{cell}}$.  Per primitive cell the count is
   $9 N_{\mathrm{cls}}$.

### What the verification does not show
1. Test 4 runs the ASR and the reference on the same $4 \times 4 \times 4$ torus (64 sites, nearest-neighbour hops).  Levels 3
   covers the whole torus.  So it tests rank 1 and the sampling of 120 configurations, not a lattice alloy.  The 7 % is not
   split into the two.
2. Limit 4 (BZ sum rule = 0) is an identity of the Fourier sum on the commensurate k grid.  Limit 5 checks 5 k points and 9
   energies at $\eta$ = 30 mRy.
3. At 16 steps the first-order projection (Q3b) is more accurate than the second-order polarisation (Q3a) in Test 2
   ($1.2 \times 10^{-6}$ against $4.0 \times 10^{-5}$).  The order in the residual does not decide the method in practice.
4. The script uses one orbital per site and no symmetry reduction.  The production plan reads the projections from vectors of
   the D4h-reduced operator (9 orbitals) and rebuilds the $p_x$, $p_y$, $d_{yz}$, $d_{xz}$ blocks with a 3-fold rotation.
   That rotation is in $O_h$, not in D4h.  It is valid only on cubic sublattices (the code requires them) and only if the node
   set is cubic.  This reconstruction is the main implementation risk, and nothing tests it yet.
5. The code-specific claims of Deliverable 2 (`find_state` before the store is freed, `at_link_c`, the projection
   support of `cf_species`) were not checked against the source.

### Tests to add before the implementation
1. Reduced against unreduced operator: the off-site $G$ from the projections of the reduced chain equals the chain of the
   unreduced operator, for all $9 \times 9$ orbital blocks and several shells (as `test/test_augcls_project.f90` does for the
   on-site element).
2. The ordered limit of the full code against the k-space bands of the same $H^{(2)}$ (`[kspace] second_order = true`,
   `std_alpha = true`), on B2 and on fcc.
3. The 3-fold reconstruction against explicit $p_x$ and $d_{yz}$ chains (unreduced operator).
4. Form 2A against Form 2B on the first shell at one energy, as an internal check of the rank-1 dressing.
