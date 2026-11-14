---
title: Recursion and disorder (method = "asr")
---

# Input: the recursion

`method = "asr"` is the default. The keys below are read only for it.

## `[asr]`

### Augmented space

| Key | Default | Meaning |
|---|---|---|
| `levels` | `4` | Recursion depth, ≥ 2 |
| `recursion_steps` | `0` | Lanczos steps; `0` = `levels`; otherwise ≥ `levels` |
| `rank` | `1` | Maximum number of simultaneously flipped sites in an augmented-space state (≥ 0). Sites beyond are an average medium. Rank 2 is a refinement; for soft-moment alloys it changes moments by several percent ([Convergence](../guides/convergence.html)) |
| `supercell` | `8` | Graph size in conventional cells; a multiple of `period` |
| `flip_radius` | `99` | Only sites within this distance may flip |
| `cluster_radius`, `inversion_radius` | code defaults | Real-space cluster for the structure constants and for their inversion |

### Energy integration

| Key | Default | Meaning |
|---|---|---|
| `contour_points` | `32` | Gauss–Legendre nodes on the complex-energy contour |
| `fermi_kt` | `0.001` | Electronic temperature (Ry); gives a smooth \(N(E_F)\) |
| `fermi_poles` | `0` | Matsubara poles enclosed (0 = automatic, range 0 to 200) |
| `ed_terminator` | `hull` | Tail of each chain: `hull`, `bp`, `avg2`. Keep `hull` ([Recursion and terminators](../theory/recursion.html)) |
| `depth_tol` | code default | A chain stops early when \(G\) changes by less than this; `levels` is the maximum |

### Route (Hamiltonian-dependent)

| Key | Default | Meaning |
|---|---|---|
| `route` | `edc` | `edc`: AMscf energy-dependent form of \(H^{(2)}\). `ed3`: the full generalised ASA problem including the \(hph\) term. See [Energy-dependent ASR](../theory/energy_dependent.html) |

### Speed

| Key | Default | Meaning |
|---|---|---|
| `energy_block` | `0` | Contour energies per chain block. `0` runs one chain per energy. `8` runs eight in one block, 1.35× faster for CoCrFeNi at about 6× the memory. Same result within the chain's length error |
| `ef_linear_max` | `0` | Ry. \(E_F\) steps up to this size are taken to first order without a new contour (exact at self-consistency); `0` = off |

An MPI build (`-DUSE_MPI`) splits the chain tasks over ranks with no input change.

### Diagnostics

| Key | Default | Meaning |
|---|---|---|
| `ef_scan`, `ef_scan_step`, `ef_scan_center` | `0` | Scan \(N(E)\) around \(E_F\) at the start potential |
| `ef_scan_dump` | — | File: dump \(N(E)\) at the start potential, then stop |
| `alpha_scale` | `1` | Multiply the screening constants by this factor |

### Constraints enforced by the reader

* `levels` ≥ 2; `recursion_steps` is 0 or ≥ `levels`.
* `supercell` is a multiple of `period`.
* `fermi_poles` in 0 to 200; `route` is `edc` or `ed3`.
* `energy_block` ≥ 0; `ef_linear_max` ≥ 0; `alpha_scale` > 0.

## `[[sro]]`

First-shell Warren–Cowley short-range order (decorated lattice, `method = "asr"`).
The first shell of the frozen centre carries the conditional concentrations
\(P(s\mid c)\) ([[asr_sro]]); shell sites are independent of each other, and sites
beyond the first shell are uncorrelated.

| Key | Meaning |
|---|---|
| `sublattices` | `["A", "B"]`; optional with one sublattice |
| `shell` | `1` |
| `pairs` | `[["Fe", "Ni", alpha], …]`. \(\alpha>0\): the pair avoids contact (clustering); \(\alpha=1\): never neighbours; \(\alpha<0\): ordering |

With one sublattice the unlike pairs are symmetric and the diagonal is derived.
With two sublattices, `[a, b, alpha]` is species `a` on the first and `b` on the
second; unset pairs are 0 and both marginals must vanish. Example:
`example/fe66ni14cr20_fcc_sro.toml`.
