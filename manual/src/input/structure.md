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
| `lattice` | `bcc` | A cubic **graph lattice** `sc`, `bcc` or `fcc` (every site of it belongs to one sublattice), a named prototype `b2`, `hex`, `hcp` or `nias` ([Hexagonal structures](#hexagonal-structures)), or `cell`: explicit `vectors` with `[[sublattice]]` positions, or (`method = "kspace"` only) with `[[site]]` entries |
| `alat` | — (required, > 0) | Lattice constant of the graph lattice, Å |
| `c_over_a` | — | Hexagonal prototypes: \(c/a\). Required for `hex` and `nias`; `hcp` defaults to the ideal \(\sqrt{8/3}\). `coa` is the same key; `c` gives \(c\) in Å instead |
| `period` | `1` | Decoration cell = `period`³ conventional cubes (≥ 1) |
| `centering` | `P` | `P`, `I` or `F`: translations of the decoration cell applied to every position |
| `vectors` | — | `cell` only: three lattice vectors `[[a1], [a2], [a3]]` in `alat` units; `[[sublattice]]` positions are then fractional in these vectors |
| `site_classes` | `false` | `cell` only: one potential per site instead of per species (local moments and charges, as an SQS supercell needs) |

`period` / `centering` and `lattice = "cell"` exclude each other.

## `[[sublattice]]` (decorated lattice)

| Key | Required | Meaning |
|---|---|---|
| `name` | no (`S`) | Label |
| `positions` | cubic lattices and `cell` | `[[x, y, z], …]`, fractional coordinates of the decoration cell (or of the prototype cell). Optional for a named prototype: see below |
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

Random fcc \(\text{Fe}_{66}\text{Ni}_{14}\text{Cr}_{20}\): `lattice = "fcc"`, `centering = "F"`, one sublattice
with `species = ["Fe", "Ni", "Cr"]` and `concentrations = [0.66, 0.14, 0.20]`
(`example/fe66ni14cr20_fcc.toml`).

Partially disordered L2\(_1\) \(\text{Co}_2\text{Fe}_{0.4}\text{Cr}_{0.6}\text{Al}\): the cell is a decoration of the
bcc graph with \(a_\text{bcc}=a/2\), so `period = 2`, `centering = "F"`; Co on
8c, Al on 4a, Fe/Cr random on 4b (`example/co2fe04cr06al_l21.toml`).

!!! tip "Site symmetry"
    The space group, the site groups and the start orbitals come from the structure. A cubic site starts from s, p\(_z\),
    t\(_{2g}\) and e\(_g\) as before; a site of lower symmetry starts from the orbital groups its site group allows.

## Hexagonal structures

The prototypes use the cell \(\mathbf a_1=(1,0,0)\), \(\mathbf a_2=(-\tfrac12,\tfrac{\sqrt3}{2},0)\),
\(\mathbf a_3=(0,0,c/a)\) in units of `alat` (= \(a\)). Basis sites in fractional coordinates:

| `lattice` | Sites | Sublattices without `positions` |
|---|---|---|
| `hex` | \((0,0,0)\) | one |
| `hcp` | \((0,0,0)\), \((\tfrac13,\tfrac23,\tfrac12)\) | one sublattice holds both sites; or two, one site each |
| `nias` | cation 2a: \((0,0,0)\), \((0,0,\tfrac12)\); anion 2c: \((\tfrac13,\tfrac23,\tfrac14)\), \((\tfrac23,\tfrac13,\tfrac34)\) | first `[[sublattice]]` = cation, second = anion; or four, one site each |

Without `positions` the sublattices take the prototype's sites in this order. With `positions` (given on every sublattice) each
position must be a site of the prototype, and every site must belong to one sublattice; this splits the sites in any other way.
A structure outside the list takes `lattice = "cell"` with `vectors` and `[[sublattice]]` positions.

```toml
[structure]
lattice = "nias"
alat = 3.619            # Angstrom
c_over_a = 1.391

[[sublattice]]
name = "cation"
species = ["Ni"]

[[sublattice]]
name = "anion"
species = ["As"]
```

!!! warning "Validated range of the hexagonal structures"
    The ordered limit agrees with `method = "kspace"` (same Hamiltonian) to the accuracy of the recursion itself, as for the
    cubic structures: hcp Co and Cu, and NiAs-type NiAs and ferromagnetic MnSb, at levels 4 to 6 and 24 to 48 recursion steps
    (charges to 0.007 e, moments to about 1 %). Antiferromagnetic order with the two cation sites on separate sublattices,
    vacancies and empty spheres, and short-range order on more than one sublattice are not validated yet.

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
