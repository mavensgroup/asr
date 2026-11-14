---
title: k-space reference (method = "kspace")
---

# Input: the k-space reference

`method = "kspace"` solves the same Hamiltonian by Brillouin-zone integration
for an **ordered** structure (one species per sublattice). It is the exact
reference of the recursion. Disorder needs `method = "asr"`.

## `[kspace]`

| Key | Default | Meaning |
|---|---|---|
| `nk` | `12` | SCF k mesh \(n_k^3\); ≥ 4 |
| `bz` | `tetra` | `tetra`: linear tetrahedron method ([[asr_tetra]]); `gauss`: Gaussian smearing |
| `smear` | `0.01` | Gaussian width, Ry (`bz = "gauss"`) |
| `sigma` | `0.7` | Hard-core screening radius of the \(\alpha\) representation |
| `trace` | `false` | Verbose iteration trace |
| `dos_nk` | `24` | k mesh of the DOS, independent of the SCF mesh |
| `sphere_ratio` | `1` | B2 only: sphere-radius ratio \(r_A/r_B\) (1 = equal spheres) |

Diagnostics:

| Key | Default | Meaning |
|---|---|---|
| `second_order` | `false` | Use the orthogonal \(H^{(2)}\) that the `edc` route uses, instead of the full generalised ASA problem |
| `alpha_scale` | `1` | Multiply the screening constants by this factor |
| `std_alpha` | `false` | Use the standard \(\alpha\) set of the recursion instead of the hard-core \(\alpha\) of `sigma` |

!!! warning "Sharp Fermi-level features"
    When \(E_F\) sits on a sharp d peak, Gaussian smearing biases \(N(E_F)\) upward
    (measured: 19.35 electrons at 6 mRy against 19.11 at 1.5 mRy for B2 FeCu).
    Prefer `bz = "tetra"`.

`[scf]` keys `sim_beta`, `start_potential`, `save_potential`, `mixing_restart` and
`[asr]` are rejected with `method = "kspace"`.
