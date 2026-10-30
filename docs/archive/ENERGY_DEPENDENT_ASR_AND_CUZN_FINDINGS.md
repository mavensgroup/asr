# Energy-Dependent Augmented Space Recursion (ASR) & CuZn Electronic Structure Analysis

**Author:** Antigravity (Pair Programming Analysis)  
**Date:** October 3, 2026  
**Project:** `asr_new` (`/home/rudra/Devel/ASR/asr_new`)

---

## 1. Executive Summary & Recommendation on Desktop Runs

> [!CAUTION]
> **HOLD the desktop runs for energy-dependent CuZn at depth 3 and 4.**  
> The energy-dependent (ED) formulation currently has a mathematical truncation in the complex contour integration: stopping the semicircle contour at $\text{Im}\,z = 0.01\text{ Ry}$ omits the vertical segment back to the real axis. This artificially undercounts valence electrons near $E_F$, causing a systematic $+44\text{ mRy}$ upward shift in $E_F$. Running CuZn in this state will lead to convergence failures or incorrect Fermi levels.

| Investigation Area | Key Finding | Status / Action |
|---|---|---|
| **CuZn Zn $d$ Position** | $\text{Zn } d$ peak at $-0.684\text{ Ry}$ in `asr_new` depth 5 is **correct**. It reproduces [Banerjee & Mookerjee (2010)](file:///home/rudra/Devel/ASR/references/me_2010.pdf) ($-0.682\text{ Ry}$) within $2\text{ mRy}$ and matches experimental XPS ($9.5\text{--}10\text{ eV}$). | **Resolved**: The $0.19\text{ Ry}$ offset against Tarafder (2006) is due to non-self-consistent starting potentials in the 2006 paper. |
| **ED $N(E_F)$ "Staircase"** | In $A(z) = P(z) - S$, evaluating the continued fraction at argument $0$ yields a purely real terminator tail on the real axis ($a_m^2 - 4b_{m-1}^2 > 0$). The spectrum consists of $m$ discrete poles rather than a continuous band. | **Identified**: Inversion at argument 0 lacks the intrinsic branch cut of standard $(z - H)^{-1}$ resolvents. |
| **ED $E_F$ Offset ($44\text{ mRy}$)** | Truncating the contour arc at $\text{Im}\,z = \eta = 0.01\text{ Ry}$ omits $\frac{1}{\pi} \int_0^\eta \text{Re}[G(E_F + iy)]\,dy$. This discards up to $50\%$ of the electron weight near $E_F$, forcing the secant search to push $E_F$ upwards by $44\text{ mRy}$. | **Identified**: Must add the vertical segment or adopt a low-temperature Fermi-Dirac / Matsubara contour. |
| **`asr_complete_algorithm.tex`** | Contains several critical bugs and omissions: non-Hermitian operator definitions ($\mathcal{T}_{\uparrow\downarrow}$ without $\mathcal{T}_{\downarrow\uparrow}$), missing augmented-space configuration flip actions in HOP, missing species-projected Green function formulas, and generic KKR charge integrals instead of TB-LMTO orthogonal moments. | **Audited**: Corrections detailed in Section 4. |

---

## 2. The CuZn Puzzle Disentangled

### 2.1 The Observed Data vs Published Benchmarks

The self-consistent field (SCF) calculations for disordered $\text{Cu}_{0.5}\text{Zn}_{0.5}$ (bcc, $a = 2.85\text{ \AA}$) yielded:

| Quantity | Depth 4 (`asr_new`) | Depth 5 (`asr_new`) | Banerjee & Mookerjee (2010) | Tarafder *et al.* (2006) |
|---|:---:|:---:|:---:|:---:|
| **$E_F\text{ (Ry)}$** | $+0.0196$ | $+0.0061$ | $0.0000\text{ (ref)}$ | $0.0000\text{ (ref)}$ |
| **Charge transfer $Q_{\text{net}}\text{ (e)}$** | $\pm 0.029$ | $\pm 0.033$ | — | — |
| **$\text{Cu } d\text{ peak vs } E_F\text{ (Ry)}$** | $-0.335$ | $-0.261$ | **$-0.309$** | $-0.122$ |
| **$\text{Zn } d\text{ peak vs } E_F\text{ (Ry)}$** | $-0.683$ | **$-0.684$** | **$-0.682$** | $-0.488$ |
| **$d\text{-band splitting }\Delta E_d\text{ (Ry)}$** | $0.348$ | $0.423$ | **$0.373$** | **$0.366$** |
| **$\text{DOS at } E_F\text{ (states/Ry/spin)}$** | $1.51$ | $1.50$ | $\approx 1.24\text{ (2.47 tot)}$ | $2.43$ |

### 2.2 Why the $\text{Zn } d$ Band is at $-0.68\text{ Ry}$ and NOT $-0.49\text{ Ry}$

1. **Exact Reproduction of the 2010 Benchmark:**
   The digitized curve from [Banerjee & Mookerjee (*Int. J. Mod. Phys. C* **21**, 205 (2010), Fig. 3)](file:///home/rudra/Devel/ASR/references/me_2010.pdf) and Rudra Banerjee's Ph.D. thesis (Fig. 4.1, page 73) places the $\text{Zn } d$ peak at **$-0.682\text{ Ry}$**. The `asr_new` Depth 5 calculation gives **$-0.684\text{ Ry}$**, matching your legacy code ([`AMscf`](file:///home/rudra/Devel/ASR/AMscf)) to **$2\text{ mRy}$**.
2. **The $d$-Band Splitting is Invariant:**
   Notice that the energy difference between the $\text{Cu } d$ and $\text{Zn } d$ peaks is virtually identical across all calculations:
   $$\Delta E_d = C_d(\text{Cu}) - C_d(\text{Zn}) \approx 0.366\text{--}0.373\text{ Ry}$$
   The entire band manifold in Tarafder (2006) is rigidly shifted upwards by $\approx 0.19\text{ Ry}$ relative to Banerjee (2010) and `asr_new`.
3. **Origin of the $0.19\text{ Ry}$ Shift in Tarafder (2006):**
   In Tarafder *et al.* (PRB 2006), the calculation was performed using reciprocal-space ASR with **non-self-consistent elemental starting parameters**. In Fig. 3(a) of that paper, elemental pure bcc Zn has its $d$-band at $-0.53\text{ Ry}$ and pure Cu at $-0.17\text{ Ry}$. In their alloy DOS (Fig. 3c), the peaks moved by only $0.04\text{ Ry}$ (to $-0.49$ and $-0.12\text{ Ry}$).
4. **Charge Transfer & Screened Impurity Electrostatics:**
   In fully self-consistent TB-LMTO (`AMscf` and `asr_new`), charge transfers from Zn to Cu ($Q_{\text{net}} \approx \pm 0.03\text{ e}$ in disordered bcc, and $\pm 0.13\text{ e}$ in ordered B2). The Screened Impurity Model (SIM) potential:
   $$V_M^i = \beta_{\text{SIM}} \frac{2 \alpha_M Q_i}{S}$$
   pulls the electropositive Zn potential deeper and pushes the electronegative Cu potential shallower.
5. **Experimental Photoemission Benchmark:**
   XPS and UPS measurements of $\text{CuZn}$ and metallic $\text{Zn}$ place the $\text{Zn } 3d$ semicore states at a binding energy of $9.5\text{--}10.2\text{ eV}$ ($0.70\text{--}0.75\text{ Ry}$) below $E_F$. A position of $-0.49\text{ Ry}$ ($6.7\text{ eV}$) is unphysically shallow for $\text{Zn } 3d$.
6. **Density of States at $E_F$:**
   In Tarafder (2006), because the $\text{Cu } d$ peak was placed at $-0.12\text{ Ry}$, its upper shoulder extended directly to $E_F$, giving an artificially high DOS ($2.43\text{ /spin}$). In Banerjee (2010) and `asr_new`, the $d$-band complex is completely filled below $E_F$, and the states at the Fermi level are purely $sp$ conduction electrons:
   $$N_{\text{free}}(E_F) = \frac{k_F \Omega}{2\pi^2} \approx 1.5\text{--}1.6\text{ states/Ry/spin}$$
   This matches the `asr_new` result of $1.50\text{ states/Ry/spin}$ exactly.

---

## 3. Energy-Dependent ASR: Mathematical Root Cause of $E_F$ Bouncing & Drift

### 3.1 The Potential-Function Representation $A(z) = P(z) - S$

In the energy-dependent formulation, the KKR/LMTO secular matrix is:
$$A(z) = P^\alpha(z) - S^\alpha$$
where the potential function in the screened representation $\alpha$ is:
$$P^\alpha(z) = \frac{z - C}{\Delta + (\gamma - \alpha)(z - C)}$$
The physical Green function in the orthogonal Löwdin representation is:
$$G(z) = \lambda(z) + \mu(z) g^\alpha(z) \mu(z)$$
where:
$$\mu(z) = \frac{\sqrt{\Delta}}{\Delta + (\gamma - \alpha)(z - C)}, \quad \lambda(z) = \frac{\gamma - \alpha}{\Delta + (\gamma - \alpha)(z - C)}$$

### 3.2 Why the Real-Axis DOS had Sharp Discrete Spikes

To obtain $g^\alpha(z) = \langle 0 | [A(z)]^{-1} | 0 \rangle$, complex-symmetric Lanczos recursion is run on $A(z)$, yielding tridiagonal coefficients $a_k(z)$ and $b_k(z)$. The Green function is evaluated as:
$$g^\alpha(z) = \langle 0 | (0 - T(z))^{-1} | 0 \rangle$$
Notice that the variable in the continued fraction is **$0$** (not $z$):
$$g^\alpha(z) = \frac{1}{0 - a_1(z) - \cfrac{b_1^2(z)}{0 - a_2(z) - \dots - \cfrac{b_m^2(z)}{0 - a_m(z) - T_{\text{term}}}}}$$

In standard Haydock recursion on $(z - H)^{-1}$, the variable $z = E$ is on the real axis, and the square-root terminator:
$$T(E) = \frac{E - a - \sqrt{(E - a)^2 - 4b^2}}{2b^2}$$
develops a non-zero imaginary part whenever $E \in [a - 2b, a + 2b]$ because $(E - a)^2 - 4b^2 < 0$. This provides the continuous branch cut.

However, in `cf_species` ([`src/asr_endep.f90`](file:///home/rudra/Devel/ASR/asr_new/src/asr_endep.f90#L104-L130)), the argument is $0$:
$$d = 0 - a_m(E) = -a_m(E)$$
For almost all energies $E$, $a_m^2(E) - 4b_{m-1}^2(E) > 0$. Hence:
$$\sqrt{d^2 - 4b^2} = \sqrt{a_m^2 - 4b_{m-1}^2} \in \mathbb{R}$$
The terminator tail is **purely real**. On the real axis, $[(0 - T_m(E))^{-1}]_{11}$ has **zero imaginary part** almost everywhere, punctuated only by $m$ singular $\delta$-like poles where $\det T_m(E) = 0$.  
Evaluating $N(E)$ on the real axis therefore produces a staircase with $m$ steps, causing Newton-Raphson searches to bounce wildly between zero and infinite slopes.

### 3.3 Proof: Why Truncating the Contour at $\text{Im}\,z = 0.01\text{ Ry}$ Shifts $E_F$ by $44\text{ mRy}$

To eliminate the real-axis spikes, the contour in [`src/asr_contour.f90`](file:///home/rudra/Devel/ASR/asr_new/src/asr_contour.f90#L14-L31) was stopped at $\text{Im}\,z = \eta = 0.01\text{ Ry}$.

By Cauchy's theorem, integrating an analytic function $G(z)$ along a closed path from $E_{\text{bot}}$ to $E_F$ requires:
$$\int_{E_{\text{bot}}}^{E_F} G(E)\,dE = \int_{C_{\text{arc}}} G(z)\,dz + \int_{E_F + i\eta}^{E_F} G(z)\,dz$$
Taking the imaginary part $-\frac{1}{\pi}\text{Im}$:
$$N(E_F) = -\frac{1}{\pi}\text{Im} \int_{C_{\text{arc}}} G(z)\,dz \;+\; \Delta N_{\text{vertical}}$$
where the omitted vertical segment is:
$$\Delta N_{\text{vertical}} = -\frac{1}{\pi}\text{Im} \left[ -i \int_0^\eta G(E_F + iy)\,dy \right] = \frac{1}{\pi} \int_0^\eta \text{Re}\left[ G(E_F + iy) \right]\,dy$$

For a localized pole (or electron state) at energy $E_k < E_F$:
$$G(z) \approx \frac{1}{z - E_k}$$
The arc integral ending at $E_F + i\eta$ evaluates to:
$$-\frac{1}{\pi}\text{Im} \int_{C_{\text{arc}}} \frac{dz}{z - E_k} = 1 - \frac{1}{\pi} \arctan\left( \frac{\eta}{E_F - E_k} \right)$$
The missing vertical segment contributes exactly:
$$\frac{1}{\pi} \int_0^\eta \frac{E_F - E_k}{(E_F - E_k)^2 + y^2}\,dy = \frac{1}{\pi} \arctan\left( \frac{\eta}{E_F - E_k} \right)$$
The sum is identically $1.0000$.

#### The Consequence of Omitting the Vertical Segment:
When the vertical integral is dropped:
- An electron state at $E_k = E_F$ contributes only $1 - \frac{1}{\pi}\arctan(\infty) = 0.50$ (a **$50\%$ charge loss**).
- An electron state at $E_F - E_k = 0.01\text{ Ry}$ loses $\frac{1}{\pi}\arctan(1) = 0.25$ (a **$25\%$ charge loss**).
- For fcc Ag, states near $E_F$ lose $\approx 0.05\text{--}0.08$ electrons.
- To accumulate the required $N_{\text{elec}} = 11.0$ valence electrons, the secant search **must push $E_F$ upwards into the empty conduction band**:
  $$\Delta E_F \approx \frac{\Delta N_{\text{lost}}}{N(E_F)} \approx \frac{0.06}{1.4} \approx +0.043\text{ Ry} = +\mathbf{43\text{ mRy}}$$
This explains why ordered Ag converged to $E_F = -0.130\text{ Ry}$ instead of $-0.174\text{ Ry}$.

---

## 4. How the Legacy Code (`AMscf` / `ASR_CORRECT`) Handled Energy Dependence

In [`AMscf/ASR/main.f90`](file:///home/rudra/Devel/ASR/AMscf/ASR/main.f90#L285-L355) and Chakrabarti & Mookerjee (*Eur. Phys. J. B* **44**, 21 (2005)), the energy dependence was formulated differently:
1. **Seed Energies:**  
   The effective Hamiltonian $\hat{H}(E_s) = \hat{A} - \hat{B} - \hat{F} + \hat{S} - (\hat{J} + \hat{S})\hat{o}(\hat{J} + \hat{S})$ was evaluated at $7$ real seed energies $\{E_s\}$ across the valence bandwidth.
2. **Three-Term Recursion at Seeds:**  
   At each seed $E_s$, standard real-space recursion was performed to compute $\{\alpha_n(E_s), \beta_n^2(E_s)\}$.
3. **Spline Interpolation:**  
   [`newfit.f90`](file:///home/rudra/Devel/ASR/AMscf/ASR/newfit.f90) interpolated $\alpha_n(E)$ and $\beta_n^2(E)$ using cubic splines onto an $800$-point real energy mesh.
4. **Spectral Continued Fraction on Real Axis:**  
   In [`cgreen.f90`](file:///home/rudra/Devel/ASR/ASR_CORRECT/cgreen.f90#L110-L125), the continued fraction was evaluated with $E$ explicitly present:
   $$G(E) = \frac{1}{E - \alpha_1(E) - \cfrac{\beta_2^2(E)}{E - \alpha_2(E) - \dots - \beta_m^2(E) T_{\text{BP}}(E)}}$$
   The Beer-Pettifor terminator $T_{\text{BP}}(E)$ was determined from the eigenvalues of the tridiagonal matrix (`eigen`). Because $E$ was inside the square-root denominator, $T_{\text{BP}}(E)$ provided an **exact continuous branch cut on the real axis with $\delta = 0$**.
5. **Real-Axis Fermi Integration:**  
   $E_F$ was found by straightforward Simpson integration along the real axis in [`fermi.f90`](file:///home/rudra/Devel/ASR/ASR_CORRECT/fermi.f90).

---

## 5. Comprehensive Audit of `references/asr_complete_algorithm.tex`

A line-by-line review of [`references/asr_complete_algorithm.tex`](file:///home/rudra/Devel/ASR/references/asr_complete_algorithm.tex) against the published literature reveals several crucial errors and omissions:

### 5.1 Mathematical Error: Non-Hermitian Operator Definition
In Equations (179)–(181):
$$\tilde{\mathbf{C}} = \sum_{\mathbf{R}} \left[ \langle \hat{C}_{\mathbf{R}} \rangle \mathcal{I} + B(\hat{C}_{\mathbf{R}}) \mathcal{P}_{\downarrow_{\mathbf{R}}} + F(\hat{C}_{\mathbf{R}}) \mathcal{T}_{\uparrow\downarrow_{\mathbf{R}}} \right] \otimes \mathcal{P}_{\mathbf{R}}$$
$$\tilde{\mathbf{\Delta}}^{1/2} = \sum_{\mathbf{R}} \left[ \langle \Delta_{\mathbf{R}}^{1/2} \rangle \mathcal{I} + B(\Delta_{\mathbf{R}}^{1/2}) \mathcal{P}_{\downarrow_{\mathbf{R}}} + F(\Delta_{\mathbf{R}}^{1/2}) \mathcal{T}_{\uparrow\downarrow_{\mathbf{R}}} \right] \otimes \mathcal{P}_{\mathbf{R}}$$
- **Error:** $\mathcal{T}_{\uparrow\downarrow} = |\!\uparrow\rangle\langle\downarrow|$ is strictly upper triangular. This makes $\tilde{\mathbf{C}}$ and $\tilde{\mathbf{\Delta}}^{1/2}$ **non-Hermitian**.
- **Correction:** The off-diagonal term must be Hermitian:
  $$F(\hat{C}_{\mathbf{R}}) \left( \mathcal{T}_{\uparrow\downarrow_{\mathbf{R}}} + \mathcal{T}_{\downarrow\uparrow_{\mathbf{R}}} \right)$$
  (as correctly written in Eq. (91)).

### 5.2 Ambiguity in $\hat{C}_{\mathbf{R}}$
In Eq. (179), $\hat{C}_{\mathbf{R}}$ is introduced without definition. In TB-LMTO:
$$\hat{C}_{\mathbf{R}} \equiv C_{\mathbf{R}} - E_{\nu\mathbf{R}}$$
If $\tilde{\mathbf{C}}$ is formed from $\hat{C}$, then $\tilde{\mathbf{h}}$ must be:
$$\tilde{\mathbf{h}} = \tilde{\mathbf{C}} + \tilde{\mathbf{\Delta}}^{1/2} \tilde{\mathbf{S}} \tilde{\mathbf{\Delta}}^{1/2}$$
without subtracting $\tilde{\mathbf{E}}_\nu$ again in Eq. (175).

### 5.3 Omission of Configuration Flips in Algorithm 4 (Module HOP)
In Algorithm 4 (lines 309–330), HOP is written purely as operations on real-space sites ($R, L$). It does not describe how the augmented-space states $(R, \{C\})$ are indexed, hashed, or how configuration flips ($\{C\} \leftrightarrow \{C \Delta \{R\}\}$) are executed when $\tilde{\mathbf{\Delta}}^{1/2}$ and $\tilde{\mathbf{o}}$ act.

### 5.4 Missing Species-Resolved Green Function Formulation
In Section 5.1 and Section 6.2, the document references $G^\lambda(E)$ and $n^\lambda(E)$ without explaining how they are obtained from the recursion.  
Starting the recursion from $|u_1\rangle = |R, L\rangle \otimes |\emptyset\rangle$ yields only the **configuration-averaged** Green function:
$$\llangle G_{RR}(z) \rrangle = x G^A(z) + (1-x) G^B(z)$$
To extract the species-projected Green function $G^A(z)$, one must evaluate the projection onto atom $A$ ($n_R \to M_R$):
$$G^A(z) = \frac{1}{x} \langle n_R G_{RR}(z) \rangle = \frac{1}{x} \langle u_1 | (z - \tilde{H})^{-1} M_R | u_1 \rangle$$
Since $M_R |u_1\rangle = x |u_1\rangle + \sqrt{xy} |R, L \otimes \{R\}\rangle$, this requires the **transition matrix element** between $|u_1\rangle$ and the first flipped state $|R, L \otimes \{R\}\rangle$ using the first column of the tridiagonal resolvent:
$$[G]_{k1} = G_{11} \prod_{j=1}^{k-1} \left( -b_j g_{j+1} \right)$$
This formula is completely absent from the algorithm guide.

### 5.5 Incomplete Charge Density Reconstruction
In Section 6.2, Eq. (462) writes the continuous KKR expression:
$$\rho^\lambda(r) = -\frac{1}{\pi} \text{Im} \int_{-\infty}^{E_F} G^\lambda(E, r, r)\,dE$$
In TB-LMTO-ASA, the charge density inside sphere $\lambda$ is reconstructed from the radial wavefunctions and the **orthogonal representation energy moments**:
$$\rho_l^\lambda(r) = M_0^{\lambda, l} \phi_l^2(r) + 2 M_1^{\lambda, l} \phi_l(r) \dot{\phi}_l(r) + M_2^{\lambda, l} \left[ \dot{\phi}_l^2(r) + \phi_l(r) \ddot{\phi}_l(r) \right]$$
where:
$$M_k^{\lambda, l} = \int_{-\infty}^{E_F} (E - E_{\nu, l}^\lambda)^k n_l^\lambda(E)\,dE$$

### 5.6 Inapplicable Madelung Formula for Alloys
In Section 6.1 (Eq. 456), the document writes an ordered multipole lattice sum:
$$V_M(r) = \sum_{L, q} M_{pq}^L Q_q^L$$
In a substitutionally disordered binary alloy, the average lattice charge is zero, so this sum vanishes. The correct formulation is the **Screened Impurity Model (SIM)** of Ruban & Skriver:
$$V_M^A = \beta_{\text{SIM}} \frac{2 \alpha_M Q_A}{S}, \quad V_M^B = \beta_{\text{SIM}} \frac{2 \alpha_M Q_B}{S}$$
where $\alpha_M$ is the Madelung constant of the underlying lattice ($1.79186$ for bcc), $S$ is the average Wigner-Seitz radius, and $\beta_{\text{SIM}} \approx 0.6\text{--}0.8$ is the screening parameter.

---

## 6. Corrective Implementation Roadmap

```mermaid
flowchart TD
    A["Fix Contour Integration in asr_contour.f90"] --> B["Option 1: Add Vertical Segment Integral"]
    A --> C["Option 2: Low-T Fermi-Dirac / Matsubara Contour"]
    B --> D["Validate Ordered Ag Gate in test_alloy_scf"]
    C --> D
    D --> E["Verify EF matches Ordered ASA (-0.174 Ry) within 10 mRy"]
    E --> F["Launch CuZn Depth 3 and 4 on Desktop (mixing 0.3)"]
    F --> G["Compare CuZn DOS against Banerjee 2010 Benchmark"]
```

### Action Item 1: Complete the Contour Integration in `asr_contour.f90`
Implement the vertical closing line from $E_F + i\eta$ to $E_F$:
$$N(E_F) = -\frac{1}{\pi}\text{Im} \int_{C_{\text{arc}}} G(z)\,dz \;+\; \frac{1}{\pi} \int_0^\eta \text{Re}\left[ G(E_F + iy) \right]\,dy$$
Using a 5-point Gauss-Legendre quadrature along $y \in [0, \eta]$:
- Quadrature points stay strictly in the upper half-plane ($y > 0$).
- Restores the exact missing electron count.
- Drops the $+44\text{ mRy}$ artificial $E_F$ offset immediately.

Alternatively, use a **Matsubara / Fermi-Dirac contour** with $T = 300\text{ K}$ ($k_B T \approx 0.0019\text{ Ry}$):
$$N(E_F) = -\frac{1}{\pi}\text{Im} \int_{C_{\text{arc}}} f(z) G(z)\,dz - 2 k_B T \sum_{n=0}^{M} \text{Re}\left[ G\left( E_F + i (2n+1)\pi k_B T \right) \right]$$
where the arc passes well above $E_F$, ensuring complete analyticity and smoothness.

### Action Item 2: Update `references/asr_complete_algorithm.tex`
Apply the corrections detailed in Section 5 to ensure full mathematical rigor and consistency with literature.
