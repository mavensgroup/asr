---
title: SCF and electrostatics
---

# SCF and electrostatics

## One iteration

```text
potentials V_p(r)  ──▶  radial solver  ──▶  φ, φ̇, φ̈  ──▶  C, √Δ, p, γ, E_ν
      ▲                                                          │
      │                                                          ▼
 Hartree + xc + Madelung                      Green function / band energies
      ▲                                                          │
      └──────  density ρ_l = M0φ² + 2M1φφ̇ + M2(φ̇²+φφ̈)  ◀── moments M_k, E_F
```

* **Core.** Relaxed core in the sphere (Hankel tail at \(S\)); the shells below the
  valence `konfig` are re-solved every iteration ([[asr_asa]],
  `core_shells_for`). The spherical free-atom stack
  ([[asr_atom]]) is used for the starting potential and as a validation
  of the radial/Poisson/XC stack.
* **Hartree.** Poisson solve with a neutral-sphere boundary (`poisson_hartree` in
  [[asr_atom]]); the net charge adds the Madelung shift.
* **xc.** libxc through its C ABI ([[asr_xc]]); LDA and PBE GGA, spin
  polarised. Densities in e/bohr³, results in Ry. GGA on a spherical
  density uses \(\sigma=(dn/dr)^2\) and
  \(v=\partial e/\partial n-r^{-2}\,d/dr[r^2\,2(\partial e/\partial\sigma)\,dn/dr]\).
* **Fermi level.** For `lmto`, from the tetrahedron/Gauss integration and the
  electron count; for alloys, a *common* \(E_F\) from \(\sum_pw_pN_p=\sum_pw_pZ_{\text{val},p}\).
* **Mixing.** Anderson (Pulay/DIIS-type) mixing of the smooth potentials and
  \(E_\nu\) ([[asr_mixing]]).

## Anderson mixing and its safeguard

With residual \(r_k=F(x_k)-x_k\),

\[
x_{k+1}=x_k+\beta r_k-\sum_j\gamma_j(\Delta x_j+\beta\,\Delta r_j),\qquad
\gamma=\arg\min\big\|r_k-\sum_j\gamma_j\Delta r_j\big\| ,
\]

with history length `m`. The history is **dropped** (plain step \(x_k+\beta r_k\))
when the Anderson step points against the current residual, \((x_{k+1}-x_k)\cdot r_k\le0\), or
when \(|r_k|\) grew by more than a factor 2. The safeguard exists because in
strongly nonlinear early iterations (magnetic starts) extrapolation otherwise
marched the potential away at the safety-clamp rate.

!!! tip "If the SCF oscillates"
    Lower `mixing` (0.3 → 0.1–0.2), and for magnetic systems lower `zeta0`. See
    [Troubleshooting](../guides/troubleshooting.html).

## Electrostatics in disordered alloys

In a random alloy a sphere's net charge \(q_\text{net}=N-Z\) is partly screened by
the surrounding spheres. The **Screened Impurity Model** (Ruban & Skriver, PRB 66,
024201 (2002)) screens the charge by the first coordination shell (radius \(R_1\)):

\[
V=-\beta_\text{SIM}\,\frac{2\,q_\text{net}}{R_1}\quad(\text{Ry}),
\]

added to the vacuum-boundary sphere Hartree potential, which already holds
\(+2q_\text{net}/S\). `sim_beta` is \(\beta_\text{SIM}\); `1` is full first-shell
screening, `0` switches the screening off (isolated sphere charges).

For **sublattices** ([[asr_class_scf]]) the potential of component \(p\) on
sublattice \(\kappa\) is

\[
V_p=\sum_{\kappa'}\Phi(\kappa,\kappa')\,\bar q_{\kappa'}-\beta_\text{SIM}\frac{2\,(q_p-\bar q_\kappa)}{R_1},
\qquad \bar q_\kappa=\sum_{p\in\kappa}x_pq_p ,
\]

the Madelung potential of the sublattice-average charges over the decoration
cell (Ewald, [[asr_ewald]]), plus the screened-impurity term for the deviation from
the sublattice average. For a binary random alloy \(\bar q=0\) and this reduces to
the SIM above; for an ordered sublattice there is no SIM term.

The screened term comes from a charge-excess model of the graph lattice ([[asr_screening]]): every site responds to its
own Madelung potential with one hardness \(a\), and `sim_beta` fixes \(a\). Without SRO it is the SIM term above. With
first-shell **short-range order** the charge deviations respond to their drive \(\delta b\) as
\(\delta q=\Gamma\,\delta b\) and \(V=a(1-\Gamma)\,\delta b\), so \(V=F\,\delta q\) with \(F=a(\Gamma^{-1}-1)\).
\(\Gamma\) is the lattice average of the screening function times the SRO structure factor. The first shell alone
gives \(\Gamma=G_{00}+G_1\hat A\). For strong clustering this has a zero, because the truncated structure factor turns
negative. The code keeps the first-shell structure factor wherever it is non-negative, which matches the ASR's own model
of uncorrelated sites beyond the first shell. Beyond that range it clips the structure factor at zero and renormalises it,
so it stays positive with as little order on further shells as possible. \(\Gamma\) then stays positive. On a Bravais lattice it vanishes only for complete segregation, where the
segregated species must stay neutral (\(F=\infty\)). The SCF iterates \(\delta b\) itself, so this limit is finite. \(\delta b\) is saved in `save_potential`. A start
without it (free atoms, or a potential of a run without SRO) starts it at zero, because such a potential holds no SRO part.
`asr.out` lists the SRO modes with \(\Gamma\), the first-shell value and \(F\).

**Sublattices with several sites per cell.** On a sublattice \(\kappa\) with \(n>1\) sites per primitive cell (hcp,
either sublattice of NiAs) the first shell is an \(n\times n\) Hermitian matrix

\[
C_{rr'}(\mathbf k)=\frac1z\sum_{j:\,r\to r'}e^{i\mathbf k\cdot\mathbf R_j},
\]

summed over the shell bonds \(\mathbf R_j\) from site \(r\) to site \(r'\). Its eigenvalues \(c_{1,m}(\mathbf k)\) are the
branches of a concentration wave. Branch \(m\) has the screening function
\(g_m(\mathbf k)=v_m^\dagger(\mathbf k)\,G_{\kappa\kappa}(\mathbf k)\,v_m(\mathbf k)\). Here \(v_m\) is the eigenvector, and
\(G_{\kappa\kappa}\) is the block of \(G=a\,(a+M)^{-1}\) on the sites of the sublattice. The closure clips the structure
factor on the branches. For \(n=1\) it is the closure above. The averages \(\langle c_{1,m}\rangle=0\),
\(\langle c_{1,m}^2\rangle=1/z\), \(\langle g_m\rangle=G_{00}\) and \(\langle g_m c_{1,m}\rangle=G_1\) hold for every \(n\).
So the linear range and its edges do not change. \(\Gamma\) stays positive. It need not decrease with \(\alpha_1\), because
\(\mathrm d\Gamma/\mathrm d\alpha_1\) is the covariance of \(g_m\) and \(c_{1,m}\) over the clipped tail divided by the variance
of \(c_{1,m}\). In NiAs with \(c/a=1.6\) the anion shell lies out of the basal plane. There \(\Gamma\) rises from \(G_{00}\) on
the ordering side, has a maximum, and falls below \(G_{00}\) at \(\alpha_1=-1\). The exact supercell model does the same. At complete clustering
(\(\alpha_1=1\)) the clipped structure factor ends evenly on the branches with \(c_{1,m}=1\). So \(\Gamma\) ends at the
average of \(g_m\) over them. When the shell joins every site of the cell in three dimensions, as on a Bravais lattice or in
ideal hcp, the only such branch is the uniform mode, and \(\Gamma=0\). On a sublattice smaller than the cell the other
sublattices screen. In NiAs a Cu-rich and a Zn-rich anion region each stay neutral through the charges of their cations, so
the anion charges need not vanish. With a shell inside the basal plane any stacking of pure planes has \(\alpha_1=1\). In both
cases the end value is positive, as in the exact supercell model.

For the **ordered B2** real-space route the two-site Madelung shifts are exact
(Ewald sums with a neutralising background: \(V_i=\sum_j\Phi(i,j)\,q_j\) in
Hartree, ×2 for Ry).

## Total energy

The total energy is the Kohn-Sham functional of the output density in the input potential. Per
component (Ry),

\[
E_p=E_\text{band}+\sum\varepsilon_\text{core}-\int n\,\tilde v_\text{in}
+\tfrac12\int n\,V_H[n]+\int n\,\varepsilon_{xc}[n]+\tfrac12\,q_p\,V^\text{M}_p ,
\qquad \tilde v=v+2Z/r .
\]

The first three terms are the kinetic energy plus \(\int n\,2Z/r\). That nuclear term
cancels against the sphere electrostatics. \(V^\text{M}_p\) is the Madelung and screening
potential of the output charges. The Fermi smearing entropy is not included. With
\(kT=1\) mRy this changes the energy by about \(\tfrac{\pi^2}{6}N(E_F)(kT)^2\), a few
hundredths of a mRy.

The band energy must be the occupied eigenvalue sum of the Hamiltonian that the route solves.

* Routes with an orthogonal Hamiltonian (`edc`, `second_order`): the \(M_k\) are energy
  moments, the sphere charge of a channel is \(M_0\), and
  \(E_\text{band}=\sum_l (M_1+E_\nu M_0)\).
* Full ASA (`[kspace]` default and `route = "ed3"`): the generalised problem
  \(Hz=\varepsilon Oz\) has \(\varepsilon\,z^TOz=z^THz\). All factors next to \(h\) in
  \(H\) are diagonal, so

  \[
  z^THz=\sum_{RL}\big[q_1+E_\nu\,(q_0+p^\gamma q_2)\big].
  \]

  The band energy of a channel is \(q_1\) plus \(E_\nu\) times its sphere charge
  \(q_0+p^\gamma q_2\). The same follows inside one sphere from \(H\phi=E_\nu\phi\) and
  \(H\dot\phi=E_\nu\dot\phi+\phi\). `test_band_energy` checks the identity against the
  eigenvalue sum of the bands to \(10^{-13}\) Ry.

## Charge definitions

In the production `asr` routes, `edc` uses the orthogonal \(M_0\) (`q0`) sphere charges, while `ed3` and `kspace` use the full ASA sphere charges \(q_0 + p\,q_2\). (The legacy reference module [[asr_alloy_scf]] also supported a `charges = "lowdin" | "sphere"` switch for order-3 fixed-operator tests.)

## Convergence measure

The change of \(rV/(1+r)\) between iterations (Ry) — *not* the unweighted
potential, which is dominated by \(r\sim10^{-5}\) Bohr at the nucleus.
