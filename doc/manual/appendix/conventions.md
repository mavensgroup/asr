---
title: Units and conventions
---

# Units and conventions

## Units

| Quantity | In the input | Inside the code | In output files |
|---|---|---|---|
| Lattice constant `alat` | Å | Bohr (÷ 0.529177210903) | Bohr and Å in `asr.out` |
| Cartesian site positions | units of `alat` | same | same |
| Energies (`emin`, `emax`, `eta`, `smear`, `fermi_kt`, `mixing` parameters) | Ry | Ry (\(\hbar^2/2m=1\), \(e^2=2\)) | Ry |
| Charges | — | electrons | electrons |
| Moments | — | | μ\(_B\) |
| Density of states | — | | states / Ry / atom (per spin channel for `lmto`, summed over spin for `alloy`/`sublattice`) |

Libxc works in Hartree and e/bohr³; the wrapper ([[asr_xc]]) converts to Rydberg.
The speed of light is \(c=274.071979\) Ry atomic units ([[asr_radial]]).

## Orbitals

Order within a shell is \(m=-l,\dots,l\): s | p(y, z, x) | d(xy, yz, 3z²−r², xz, x²−y²).
Nine orbitals per site (s, p, d). The \(\mathrm{D}_{4h}\) start orbitals of the
augmented-space recursion are s, p\(_z\), t\(_{2g}\) (xy), e\(_g\) (3z²−r²) with
degeneracy weights 1, 3, 3, 2.

## Gamma representation

Potential parameters \((C,\sqrt\Delta,p,\gamma)\) printed in `asr.out` and
`params_file` are in the **gamma representation** (\(\alpha=\gamma\)). To use them
with structure constants \(S^\alpha\) of another screening, transform with the
standard (`pptrns`) conversion; the code does this internally and the sign of
\(\sqrt\Delta\) is a gauge.

## Lattices

| Name | Meaning | Where |
|---|---|---|
| `sc`, `bcc`, `fcc` | One-atom cubic Bravais lattices | `model`, `lmto`, `alloy` (random: bcc, fcc), `sublattice` (graph) |
| `b2` | CsCl, two atoms: A at (0,0,0), B at (½,½,½) | `lmto`, `alloy` (ordered) |
| `cell` | General cell, `vectors` and `[[site]]` | `lmto` |

Integer half-lattice coordinates on the graphs: sc all even; bcc all of one
parity; fcc coordinate sum even ([[asr_lattice_toy]], [[asr_lattices]]).

## Spin

`spin = 1`: restricted. `spin = 2`: collinear (up and down). In the alloy
quantities, `m` is the moment per atom of that species in μ\(_B\).

## Naming

* *Class* — an inequivalent ASA sphere (own radius, \(Z\), core, valence).
* *Species* — a chemical element entry.
* *Component* — a (sublattice, species) pair.
* *Node* — a basis state of the augmented space: a site plus its flipped sites.
