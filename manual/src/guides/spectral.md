---
title: Bloch spectral functions
---

# Bloch spectral functions (`[spectral]`)

The configuration-averaged Bloch spectral function \(A(\mathbf k, E)\) is the disordered analogue of a band structure: at each \(\mathbf k\) along a high-symmetry path it shows broadened quasiparticle peaks whose widths reflect disorder scattering ([[asr_spectral]]).

In `asr`, `[spectral]` runs **after** the SCF converges and is available on `method = "asr"` with `route = "edc"`. Every sublattice must be a single lattice of translations (random alloys, B2, L2\(_1\); not L1\(_2\)).

---

## 1. Setting up `[spectral]`

Here is `example/cuzn_alloy_spectral.toml` for random bcc \(\text{Cu}_{50}\text{Zn}_{50}\):

```toml
[spectral]
kpath = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
labels = ["G", "H", "N", "G", "P"]
nk = 20
emin = -0.9
emax = 0.5
npoints = 141
eta = 0.01
rmax = 2.0
steps = 48               # off-site elements converge later than the on-site one
```

### Key choices

- **`kpath` coordinates:** Cartesian coordinates in units of \(2\pi / a\) of the **graph lattice** (`alat`), not fractional coordinates of the primitive cell.
- **`rmax`:** `rmax` (in units of `alat`) sets the real-space cutoff of the Fourier sum. The site cluster of the recursion always holds the ball of `rmax`. The targets must also lie inside the node set that `levels` generates; a run stops if one does not, and warns when the targets reach beyond its interior.
- **`window` (`"ball"` vs `"sharp"`):**
  - `"ball"` (default): weights displacements by \(w(R) = (1 - R/r_{\max})^2 (1 + R/(2r_{\max}))\). Its Fourier transform is non-negative, so \(A(\mathbf k, E) \ge 0\) wherever the real-space elements are converged, while preserving the BZ average and energy integral.
  - `"sharp"`: sets \(w(R) = 1\) up to `rmax`. Gives narrower \(\mathbf k\)-features at fixed `rmax`, but introduces Gibbs ringing (which can drive \(A < 0\), e.g. \(-30\,\%\) on ordered Cu at \(\eta = 0.01\text{ Ry}\)).
- **`steps`:** Off-site Green function elements converge more slowly with Lanczos depth than the on-site element. Setting `steps` larger than the SCF `recursion_steps` (e.g. `48`) improves off-site accuracy without repeating the SCF.

---

## 2. Built-in quality checks in `asr.out`

When `[spectral]` is active, `asr.out` reports a `=== Bloch spectral function ===` section with:

- **Target hop distance vs interior of the node set:** checks that the augmented-space node set covers `rmax`.
- **Reciprocity defect:** difference between the two origins of a pair \(\langle G_{0L,RL'}\rangle\) and \(\langle G_{0L',-RL}\rangle\), which come from separate chains.
- **Truncation check:** change in \(G\) when the chains are cut 4 steps earlier.
- **Negativity count:** number of \((\mathbf k, E)\) points with \(A < 0\) (indicates chains that are too short, or ringing when `window = "sharp"`).
- **Spatial decay:** decay of off-site \(\mathrm{Im}\,G\) with shell distance.

---

## 3. Plotting `spectral.dat` (and `.bands`)

`spectral.dat` (default filename) contains one block of energies per \(\mathbf k\)-point separated by blank lines. For an ordered single-species structure, `<file>.bands` is also written with the exact eigenvalues of \(H^{(2)}(\mathbf k)\) for direct overlay.

### Gnuplot

```gnuplot
set pm3d map
set xlabel "k-path (2*pi/a)"
set ylabel "E - E_F (Ry)"
# Replace EF with the Fermi energy from the header of spectral.dat
EF = -0.05
splot "spectral.dat" using 1:($5 - EF):6 with pm3d notitle
```

### Python (`numpy` + `matplotlib`)

```python
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("spectral.dat")
# Number of energies per k-point from the header (# mesh : nk k points x ne energies)
ne = 141
nk = len(data) // ne

dist = data[::ne, 0]
E = data[:ne, 4]
A = data[:, 5].reshape(nk, ne).T

plt.pcolormesh(dist, E, A, shading="auto", cmap="inferno")
plt.xlabel(r"$k$-path ($2\pi/a$)")
plt.ylabel(r"$E$ (Ry)")
plt.colorbar(label=r"$A(\mathbf{k}, E)$ (states/Ry)")
plt.tight_layout()
plt.savefig("spectral.png", dpi=200)
```
