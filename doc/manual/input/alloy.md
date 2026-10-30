---
title: kind = "alloy"
---

# `kind = "alloy"` — binary random alloy (and ordered B2)

Self-consistent binary substitutional alloy A\(_x\)B\(_{1-x}\) on an fcc or bcc
lattice, using the augmented-space recursion. Species- and l-resolved DOS are
written. The ordered CsCl (`b2`) structure is available as the *non-random*
counterpart: it uses the same machinery without augmented space (nothing is
random), a real-space recursion on a periodic bcc supercell with exact two-site
Madelung shifts. Driver: [[asr_alloy_driver]]; SCF: [[asr_alloy_scf]].

```toml
[calculation]
name = "CuZn bcc alloy"
kind = "alloy"

[structure]
lattice = "bcc"
alat = 2.85              # Å

[[species]]
symbol = "Cu"
z = 29
nval = 11
konfig = [4, 4, 3, 4]
occupations = [1.0, 10.0]
concentration = 0.5      # x_A; species B gets 1 - x_A

[[species]]
symbol = "Zn"
z = 30
nval = 12
konfig = [4, 4, 3, 4]
occupations = [2.0, 10.0]

[scf]
spin = 1
xc = "vbh"
sim_beta = 1.0
max_iterations = 60

[asr]
levels = 3               # recursion depth — convergence must be checked
rank = 2                 # max flipped sites per augmented-space state
supercell = 8            # graph size, conventional cells
symmetry = true          # exact D4h reduction
hamiltonian = "energy-independent"

[dos]
emin = -0.9
emax = 0.7
npoints = 801
eta = 0.02

[output]
dos_file = "cuzn_alloy_dos.dat"
```

## `[structure]`

| Key | Meaning |
|---|---|
| `lattice` | `fcc`, `bcc` (random alloy) or `b2` (ordered, 50:50, A on the corner and B on the body centre) |
| `alat` | Lattice constant in Å (> 0). For `b2` this is the cubic cell. |

## `[[species]]`

**Exactly two** entries. Keys as in the [overview](index.html), plus:

* `concentration` — required on the **first** species for `fcc`/`bcc`; strictly
  between 0 and 1. If given on both, they must sum to 1. For `b2` omit it (or give
  0.5).

## `[scf]`

Common keys plus the alloy-specific ones:

| Key | Default | Meaning |
|---|---|---|
| `spin`, `xc`, `sigma`, `mixing`, `relativistic` | see overview | |
| `tolerance` | `1e-5` | |
| `max_iterations` | `60` | |
| `zeta0` | `0.1` | |
| `sim_beta` | `1.0` | Screened Impurity Model strength. The net charge of a sphere is screened by its first shell, \(V = -\beta_\text{SIM}\,2\,q_\text{net}/R_1\) (Ry). `1` is full first-shell screening. See [Electrostatics](../theory/scf_electrostatics.html). |
| `save_potential` | `""` | Write the converged species potentials and \(E_\nu\) to this file (restart, fixed-potential diagnostics) |
| `start_potential` | `""` | Start from the potentials in this file instead of free atoms |
| `flip_spin` | `""` | With `start_potential` and `spin = 2`: `"A"`, `"B"` or `"AB"` swaps the spin channels of that species in the start potential (antiparallel starts). Requires `start_potential`. |

## `[asr]` — the recursion

This is the most important table. Convergence strategy is in the
[Convergence guide](../guides/convergence.html).

### Size and truncation of the augmented space

| Key | Default | Meaning |
|---|---|---|
| `levels` | `4` | Depth of the generated node set / Lanczos levels (≥ 2). **The main accuracy lever.** |
| `recursion_steps` | `0` | Lanczos steps run on the generated node set; `0` means `= levels`. If non-zero it must be ≥ `levels`; the node set then acts as a hard wall. |
| `rank` | `2` | Maximum number of simultaneously flipped sites in an augmented-space state (≥ 0). Sites beyond are treated as an average medium. |
| `supercell` | `8` | Size of the periodic graph in conventional cells (≥ 2 when `symmetry = true`). |
| `flip_radius` | `99` | Only sites within this graph distance of the origin may flip. |
| `symmetry` | `true` | Exact D\(_{4h}\) point-group reduction of the augmented space (≈ 40× fewer nodes). **Required** by the energy-dependent route for random alloys. |

### Hamiltonian and screening

| Key | Default | Values | Meaning |
|---|---|---|---|
| `hamiltonian` | `energy-independent` | `energy-independent`, `energy-dependent`, `asr_ed` | Which recursion is run — see [Energy-dependent ASR](../theory/energy_dependent.html). |
| `order` | `2` | 2, 3 | Order of the energy-independent orthonormal ASA Hamiltonian: 2 is \(E_\nu + h - hoh\); 3 adds \(hohoh - (hphh + hhph)/2\) |
| `order_switch` | `1e-3` | | Order 3: iterate with order 2 until \(\Delta V\) < this value, then switch (cubic terms are unreliable far from self-consistency) |
| `screening` | `auto` | `auto`, `standard`, `hardcore` | `standard` uses the TB-LMTO \(\alpha\) = (0.3485, 0.05303, 0.010714); `hardcore` uses `sigma`; `auto` = `standard` for the energy-independent route and `hardcore` for `energy-dependent` (but `standard` for `asr_ed`). |
| `cluster_radius` | 1.2 | alat | Kept range of \(S^\alpha\) |
| `inversion_radius` | 0 | alat | Cluster used to invert \(1 - S^0\alpha\); `0` (or ≤ `cluster_radius`) means invert on the kept cluster |
| `charges` | `lowdin` | `lowdin`, `sphere` | Order 3: orthonormal-orbital weights vs. ASA sphere charges (the latter loses ≈ 0.03 e per cell at this order) |
| `density` | `expand` | `expand`, `stuttgart` | \(M_2\) term of the density: `expand` = \(M_2(\dot\phi^2+\phi\ddot\phi)\); `stuttgart` = \(M_2(\dot\phi^2-p\phi^2)\) |
| `terminator` | `avg2` | `avg2`, `bp` | Tail of the energy-independent route: average of the last two coefficients, or Beer–Pettifor |

### Energy integration

| Key | Default | Meaning |
|---|---|---|
| `contour_points` | `60` | Gauss–Legendre points on the complex-energy contour |
| `fermi_kt` | `0` | Energy-dependent routes: electronic temperature (Ry). `> 0` gives a Fermi-smeared contour with smooth \(N(E_F)\); `0` is a sharp \(E_F\). |
| `fermi_poles` | `0` | Matsubara poles enclosed by the Fermi contour (path height \(2n\pi kT\)); `0` selects them automatically (≈ 0.1 Ry). Range 0–200. |

### Energy-dependent route only

| Key | Default | Meaning |
|---|---|---|
| `potential_function` | `parametrised` | `exact`: \(P_l(z)\) is solved from the radial equation at each contour energy ([[asr_pfun]]) |
| `exact_density` | `true` | With `exact`: density from \(-\tfrac1\pi\,\mathrm{Im}\int[G\phi(z)^2+T(z)]\,dz\); `false` uses energy moments of the exact \(G\) |
| `ed_terminator` | `hull` | Tail of the chain: `hull`, `bp`, `avg2` |
| `depth_tol` | `1e-5` | Stop a chain when \(G\) changes by less than this (`levels` is the maximum) |
| `kkr_check` | `0` | Ordered ED: after the SCF run a same-model k-space KKR-ASA check on an \(n_k^3\) mesh (`0` = off) |
| `tail_scan` | `0` | **Diagnostic** (`asr_ed`): tabulate moments against termination depth 8 … `tail_scan` and tail type, then stop |

!!! danger "Validity constraints enforced by the reader"
    * Energy-dependent recursion for a random alloy **needs** `symmetry = true`.
    * `recursion_steps` must be `0` or ≥ `levels`.
    * `fermi_poles` is limited to 0–200; `order` to 2 or 3.
    * `lattice = "b2"` forbids a `concentration` other than 0.5.

## `[dos]` and `[output]`

| Key | Default | Meaning |
|---|---|---|
| `emin`, `emax` | −1, 1 | Energy window (Ry), `emax > emin` |
| `npoints` | 801 | Grid points |
| `eta` | `0.02` | Imaginary part of the energy for the energy-independent route |
| `dos_file` | `dos.dat` | |
| `out_file` | `asr.out` | |

## Examples shipped with the code

| File | Demonstrates |
|---|---|
| `cuzn_alloy.toml` | Random bcc Cu\(_{50}\)Zn\(_{50}\), energy-independent \(H^{(2)}\) |
| `cuzn_alloy_ed.toml` | Same system, energy-dependent \(P(z)-S\) recursion |
| `cuzn_alloy_asr_ed.toml` | Same, AMscf-form energy-dependent recursion (`asr_ed`) |
| `cuzn_b2_realspace.toml` | Ordered B2 CuZn by real-space recursion |
| `cuzn_b2_ed.toml` | Ordered B2 CuZn, energy-dependent |
| `feco_b2_realspace.toml` | Ferromagnetic ordered B2 FeCo by real-space recursion |
