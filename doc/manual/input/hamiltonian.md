---
title: Hamiltonian (TB-LMTO-ASA)
---

# Input: Hamiltonian

These keys belong to the Hamiltonian, not to the method. Today the only
Hamiltonian is [TB-LMTO-ASA](../hamiltonian/tb_lmto_asa.html), so they are
documented for it. A second Hamiltonian would bring its own species block
and a selector key ([Adding a Hamiltonian](../hamiltonian/adding.html)).

## `[[species]]`

The atomic data from which the code generates all potential parameters itself
(radial solver, free-atom start). Any number of entries, unique symbols.

| Key | Required | Meaning |
|---|---|---|
| `symbol` | yes | Label; matches `species` in `[[sublattice]]` / `[[site]]` |
| `z` | yes | Atomic number |
| `nval` | yes | Number of valence electrons |
| `konfig` | no | `[ns, np, nd, nf]`: valence principal quantum numbers, e.g. `[4, 4, 3, 4]` for Fe |
| `occupations` | no | `[n_s, n_d]` of the neutral-atom start, e.g. `[2.0, 6.0]` for Fe |
| `spin_sign` | no (`1`) | `1` or `-1`: sign of the start polarisation of this species |

Concentrations are not here; they live in `[[sublattice]]`.

## `[scf]`

| Key | Default | Meaning |
|---|---|---|
| `spin` | `1` | `1` non-magnetic, `2` collinear spin-polarised |
| `xc` | `vbh` | libxc functional name (`vbh`: von Barth–Hedin LDA) |
| `relativistic` | `true` | Scalar-relativistic radial solver; `false` is the non-relativistic limit |
| `zeta0` | `0.1` (`0.4` for `kspace`) | Initial spin polarisation (seed) |
| `mixing` | `0.3` | Anderson mixing parameter |
| `tolerance` | `1e-5` (`1e-6` for `kspace`) | Convergence of the potential, Ry |
| `max_iterations` | `60` (`100` for `kspace`) | Iteration limit |

Keys read with `method = "asr"` only:

| Key | Default | Meaning |
|---|---|---|
| `sim_beta` | `1.0` | Hardness of the screened-impurity (SIM) charge-excess model, \(0\le\beta<\) the lattice maximum (about 1.39 on bcc). The net charge of a sphere is screened by its first shell, \(V=-\beta\,2q_\text{net}/R_1\) (Ry). With `[[sro]]` it also sets the screened potential of the correlated first shell. See [Electrostatics](../theory/scf_electrostatics.html) |
| `mixing_restart` | = `mixing` | Linear step after an Anderson restart; smaller for stiff charge-transfer maps |
| `save_potential` | — | File written every iteration |
| `start_potential` | — | Start from a saved potential |

An exit with code 3 means the SCF did not converge within `max_iterations`.
