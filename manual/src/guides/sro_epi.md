---
title: Short-range order and EPIs
---

# Short-range order and effective pair interactions

`asr` supports two complementary directions for chemical short-range order (SRO):

1. **Prescribed first-shell SRO (`[[sro]]`):** fix Warren–Cowley parameters \(\alpha_{ij}\) on the first coordination shell and solve the electronic structure and electrostatics self-consistently.
2. **Effective pair interactions (`[epi]`):** compute the pair interactions \(U(c,s)\) from the converged random alloy so you can study ordering tendencies or feed them into a Monte Carlo / cluster-variation calculation.

---

## 1. Prescribed first-shell SRO (`[[sro]]`)

### How the model works

In the augmented-space recursion of a central species \(i\), the first coordination shell of the frozen centre carries the conditional probabilities

\[
P(j \mid i) = c_j\,(1 - \alpha_{ij}),
\]

where \(c_j\) is the sublattice concentration of species \(j\) ([[asr_sro]]). Sites in the first shell are independent of each other, and sites beyond the first shell are uncorrelated.

- \(\alpha_{ij} > 0\) (\(i \neq j\)): species \(i\) and \(j\) avoid contact (clustering / segregation); \(\alpha_{ij} = 1\) means they are never nearest neighbours.
- \(\alpha_{ij} < 0\) (\(i \neq j\)): unlike neighbours are preferred (chemical ordering).
- \(\alpha_{ij} = 0\): random (uncorrelated) first shell.

On one sublattice, you specify only the unlike pairs; the diagonal \(\alpha_{ii}\) is derived automatically from \(\sum_j P(j \mid i) = 1\). With two sublattices, `sublattices = ["A", "B"]` is required, `[a, b, alpha]` gives species `a` on the first sublattice and `b` on the second, unset pairs default to `0`, and both marginals must vanish.

### Worked example: fcc \(\text{Fe}_{66}\text{Ni}_{14}\text{Cr}_{20}\)

`example/fe66ni14cr20_fcc_sro.toml` sets complete Fe segregation from Ni and Cr (\(\alpha_{\text{FeNi}} = 1\), \(\alpha_{\text{FeCr}} = 1\)) while leaving Ni–Cr random (\(\alpha_{\text{NiCr}} = 0\)):

```toml
[[sublattice]]
name = "fcc"
positions = [[0.0, 0.0, 0.0]]
species = ["Fe", "Ni", "Cr"]
concentrations = [0.66, 0.14, 0.20]

[[sro]]
pairs = [["Fe", "Ni", 1.0], ["Fe", "Cr", 1.0], ["Ni", "Cr", 0.0]]
```

In the reference output (`example/reference/fe66ni14cr20_fcc_sro/asr.out`), the `=== short-range order ===` section reports:

```text
 conditional probability P(j|i) of species j on a first-shell site, centre species i
   centre    P(Fe|i)    P(Ni|i)    P(Cr|i)
   ---------------------------------------
   Fe         1.0000     0.0000     0.0000
   Ni         0.0000     0.8000     0.2000
   Cr         0.0000     0.1400     0.8600

 screening modes of the SRO closure
   mode  lambda(Ahat)  Gamma = dq/db  G00 + G1 lambda    F [Ry/e]
   --------------------------------------------------------------
      1        0.0000        1.35649          1.35649     -0.4169
      2        7.9200        0.45535          0.48673      1.8976
      3       12.0000        0.00000          0.03867    infinite
 F infinite: Gamma = 0, full clustering on that mode; its charge deviation vanishes
```

Notice mode 3 (\(\lambda = 12\), full Fe clustering): \(\Gamma = 0\) and \(F = \infty\), so the Fe charge deviation vanishes and Fe stays neutral without any electrostatic singularity. The SCF iterates the drive \(u = \delta b\) directly (saved in `save_potential`), keeping strong clustering stable.

!!! warning "Strong clustering and `sim_beta`"
    With `sim_beta` above about `1`, the result for strong clustering depends on the order beyond the first shell, which `[[sro]]` does not fix. Keep `sim_beta = 1.0` unless you have a specific reason to vary it.

---

## 2. Effective pair interactions (`[epi]`)

### What `[epi]` computes

After the SCF converges, `[epi]` computes the mixed second difference of the grand potential \(\Omega\) at frozen potentials ([[asr_epi]]):

\[
U(c, s) = \Omega(c, s) - \Omega(c_0, s) - \Omega(c, s_0) + \Omega(c_0, s_0),
\]

where \(c_0, s_0\) are the first species of the two sublattices. On one sublattice, the exchange interaction of an unlike pair \((p, q)\) is

\[
V(p, q) = U(p, p) + U(q, q) - U(p, q) - U(q, p),
\]

where \(V > 0\) favours unlike neighbours (ordering) and \(V < 0\) favours like neighbours (clustering). Both `route = "edc"` and `route = "ed3"` are supported.

### Worked example: bcc \(\text{Cu}_{50}\text{Zn}_{50}\)

`example/cuzn_alloy_epi.toml` computes pair interactions over the first 3 coordination shells (8, 6, and 12 sites on bcc):

```toml
[epi]
shells = 3               # coordination shells (8, 6 and 12 sites on bcc)
lambda_points = 4        # nodes of the species path Cu -> Zn
file = "cuzn_epi.toml"
```

### Checking accuracy and reducing cost

1. **Asymmetry check:** \(U(p, q)\) and \(U(q, p)\) are evaluated from different chains (path site and frozen site swapped). `asr.out` prints their difference as a direct error measure of the chain length and path quadrature.
2. **Path quadrature (`lambda_points`):**
   - `4` (default): ~0.1 % accuracy even for chemically dissimilar species.
   - `2`: ~0.01 % for neighbours in a transition series, at half the cost.
   - `1`: midpoint rule only; use for a quick survey, not final numbers.
3. **`symmetric = true`:** On one sublattice, computes only \(s \ge c\) and copies \(U(s, c) = U(c, s)\). Cuts chains by 25 % (4 species) to 30 % (5 species), but disables the asymmetry check—verify one shell without it first.
4. **Energy blocking (`[asr] energy_block`):** Shares the same contour-energy blocking kernels (`asr_augblk`) as the SCF.

### Reading `epi.toml`

`epi.toml` is a standalone TOML file containing run metadata (`alat` in Bohr, `fermi_energy` in Ry), `[[sublattice]]` definitions, and one `[[pair]]` table per pair type with `shell`, `vector` (units of `alat / 2`), `distance` (Bohr), `multiplicity`, `kernel` \(W(R)\), `u_band`, and `u_es` (both in Ry). The total interaction matrix is `u_band + u_es`.

```python
import tomllib

with open("cuzn_epi.toml", "rb") as f:
    epi = tomllib.load(f)

for pair in epi["pair"]:
    ub = pair["u_band"]
    ue = pair["u_es"]
    u = [[ub[i][j] + ue[i][j] for j in range(len(ub[0]))] for i in range(len(ub))]
    v_ab = u[0][0] + u[1][1] - u[0][1] - u[1][0]
    print(f"shell {pair['shell']} (d={pair['distance']:.3f} Bohr): V_AB = {v_ab * 1e3:.3f} mRy")
```
