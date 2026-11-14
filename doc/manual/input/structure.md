---
title: Structure and disorder
---

# Structure and disorder

Nothing on this page depends on the Hamiltonian. It describes a periodic
structure, the sublattices of its decoration cell, and the species and
concentrations on each sublattice. The same blocks run under both methods.

## `[structure]`

| Key | Default | Meaning |
|---|---|---|
| `lattice` | `bcc` | The **graph lattice** `sc`, `bcc` or `fcc`: every site of it belongs to one sublattice. `cell` (general cell, `method = "kspace"` only) takes `[[site]]` entries |
| `alat` | — (required, > 0) | Lattice constant of the graph lattice, Å |
| `period` | `1` | Decoration cell = `period`³ conventional cubes (≥ 1) |
| `centering` | `P` | `P`, `I` or `F`: translations of the decoration cell applied to every position |
| `vectors` | — | `cell` only: three lattice vectors `[[a1], [a2], [a3]]` in `alat` units |
| `site_classes` | `false` | `cell` only: one potential per site instead of per species (local moments and charges, as an SQS supercell needs) |

`period` / `centering` and `lattice = "cell"` exclude each other.

## `[[sublattice]]` (decorated lattice)

| Key | Required | Meaning |
|---|---|---|
| `name` | no (`S`) | Label |
| `positions` | yes | `[[x, y, z], …]`, fractional coordinates of the decoration cell |
| `species` | yes | Species symbols (each a `[[species]]`); no duplicates within a sublattice |
| `concentrations` | with > 1 species | Same length as `species`; each > 0; sum to 1. May be omitted for one species |

One species on every sublattice is an ordered compound. Several species on a
sublattice is disorder on that sublattice. Any mix of the two is allowed
(partial disorder), with any number of species (binary to high-entropy).

### Examples

Pure bcc Fe: one sublattice, `centering = "I"`.

```toml
[structure]
lattice = "bcc"
alat = 2.8665
centering = "I"

[[sublattice]]
name = "bcc"
positions = [[0.0, 0.0, 0.0]]
species = ["Fe"]
```

Ordered B2 FeCo: `centering = "P"`, two sublattices (cube corner and centre).

Random fcc \ce{Fe66Ni14Cr20}: `lattice = "fcc"`, `centering = "F"`, one sublattice
with `species = ["Fe", "Ni", "Cr"]` and `concentrations = [0.66, 0.14, 0.20]`
(`example/fe66ni14cr20_fcc.toml`).

Partially disordered L2\(_1\) \ce{Co2Fe_{0.4}Cr_{0.6}Al}: the cell is a decoration of the
bcc graph with \(a_\text{bcc}=a/2\), so `period = 2`, `centering = "F"`; Co on
8c, Al on 4a, Fe/Cr random on 4b (`example/co2fe04cr06al_l21.toml`).

!!! tip "Cubic site symmetry"
    The start orbitals are s, p\(_z\), t\(_{2g}\) and e\(_g\); the code checks that
    each site has the cubic symmetry this assumes.

## `[[site]]` (general cell, k-space only)

`lattice = "cell"` takes one `[[site]]` per atom: `position = [x, y, z]`
(Cartesian, `alat` units) and `species` (a `[[species]]` symbol). Use it for
special quasirandom structures ([SQS generation](../guides/sqs.html)) as an
ordered-supercell cross-check of the recursion.

```toml
[structure]
lattice = "cell"
alat = 2.86
vectors = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
site_classes = true

[[site]]
position = [0.0, 0.0, 0.0]
species = "Fe"
```
