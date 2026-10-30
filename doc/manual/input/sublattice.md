---
title: kind = "sublattice"
---

# `kind = "sublattice"` — decorated lattices

Ordered compounds, **partial (sublattice) disorder** and **multicomponent (high-entropy)**
alloys on one input format. Sites of an sc/bcc/fcc *graph* are assigned to
**sublattices** by a *decoration cell*; each sublattice carries its own list of
species and concentrations. The SCF uses the energy-dependent augmented-space
recursion (AMscf form). Driver: [[asr_class_driver]]; SCF: [[asr_class_scf]];
geometry: [[asr_sublattice]]; operator: [[asr_augcls]], [[asr_edc]].

## Idea

The graph lattice is the lattice on which hopping is defined (e.g. bcc with
\(a_\text{bcc}=a/2\) for an L2\(_1\) Heusler cell). A decoration cell of
`period`³ conventional graph cubes contains the sublattice sites, given in
**fractional coordinates of the decoration cell**. Every graph site of the
decoration cell must belong to **exactly one** sublattice.

A *component* \(p=(\kappa,s)\) — sublattice \(\kappa\), species \(s\) — has its
own potential, \(E_\nu\) and potential parameters, and weight \(w_p=f_\kappa
x_{s,\kappa}\) per site. All spheres have one common radius.

## Example: L2\(_1\) Co\(_2\)Fe\(_{0.4}\)Cr\(_{0.6}\)Al

Co on 8c, Al on 4a, Fe/Cr random on 4b:

```toml
[calculation]
name = "co2fe04cr06al L2_1, kind sublattice"
kind = "sublattice"

[structure]
lattice = "bcc"          # the graph
alat = 2.8635            # Å, of the graph lattice (= a/2 for a = 5.727 Å)
period = 2               # decoration cell = 2×2×2 graph cubes
centering = "F"          # translations of the cell

[[sublattice]]
name = "8c"
positions = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]
species = ["Co"]

[[sublattice]]
name = "4a"
positions = [[0.0, 0.0, 0.0]]
species = ["Al"]

[[sublattice]]
name = "4b"
positions = [[0.5, 0.5, 0.5]]
species = ["Fe", "Cr"]
concentrations = [0.4, 0.6]

[[species]]
symbol = "Co"
z = 27
nval = 9
konfig = [4, 4, 3, 4]
occupations = [2.0, 7.0]
# … Fe, Cr, Al likewise
```

## `[structure]`

| Key | Default | Meaning |
|---|---|---|
| `lattice` | `bcc` | Graph lattice: `sc`, `bcc`, `fcc` |
| `alat` | — (> 0) | Lattice constant of the **graph** lattice, Å |
| `period` | `1` | Decoration cell = `period`³ conventional cubes (≥ 1) |
| `centering` | `P` | `P`, `I` or `F`: translations of the decoration cell applied to every position |

## `[[sublattice]]`

| Key | Required | Meaning |
|---|---|---|
| `name` | no (default `S`) | Label |
| `positions` | yes | `[[x,y,z], …]`, fractional coordinates of the decoration cell |
| `species` | yes | List of species symbols (each a `[[species]]`); no duplicates within a sublattice |
| `concentrations` | with > 1 species | Same length as `species`; each > 0; sum to 1. May be omitted for a single species. |

## `[[species]]`

As for [`alloy`](alloy.html) but any number, with unique symbols. No `concentration`
key — concentrations live in `[[sublattice]]`.

## `[scf]`

`spin`, `xc`, `sim_beta`, `mixing`, `tolerance` (1e-5), `max_iterations` (60),
`zeta0` (0.1), `relativistic`, `start_potential`, `save_potential`. There is no
`sigma` and no `flip_spin` for this kind.

## `[asr]`

| Key | Default | Meaning |
|---|---|---|
| `levels` | `4` | Recursion depth |
| `recursion_steps` | `0` | Lanczos steps (`0` = `levels`) |
| `rank` | `1` | Max flipped sites per state |
| `supercell` | `8` | Graph size in conventional cells; **a multiple of `period`** |
| `flip_radius` | `99` | Only sites within this distance may flip |
| `contour_points` | `32` | Contour nodes |
| `cluster_radius`, `inversion_radius` | code defaults | As in [alloy](alloy.html) |
| `ed_terminator` | `hull` | `hull`, `bp`, `avg2`; keep `hull` (see [Recursion and terminators](../theory/recursion.html)) |
| `fermi_kt` | `0.001` | Ry; Fermi-smeared contour |
| `fermi_poles` | `0` | Automatic |
| `depth_tol` | code default | Chain stopping tolerance |

## `[dos]`, `[output]`

`emin` (−1), `emax` (1), `npoints` (**401**), `eta` (0.02); `dos_file`, `out_file`.

!!! tip "Cubic site symmetry"
    The start orbitals are s, p\(_z\), t\(_{2g}\) and e\(_g\); the code checks that
    each site has the cubic symmetry this assumes.

## Special quasirandom structures

For ordered *supercell* treatments of disorder (a cross-check against the
recursion), see [SQS generation](../guides/sqs.html) and `kind = "lmto"` with
`lattice = "cell"`.
