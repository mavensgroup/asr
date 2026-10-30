---
title: kind = "lmto"
---

# `kind = "lmto"` — ordered TB-LMTO-ASA

Self-consistent tight-binding LMTO in the atomic-sphere approximation for an
**ordered** crystal: the code generates the potential parameters
(\(E_\nu, C, \sqrt\Delta, p, \gamma\)) from its own radial solver, iterates the
charge density to self-consistency with a k-space integration, then writes the
parameters and an l-resolved DOS. Driver: [[asr_lmto_driver]]; the SCF itself is
[[asr_asa]].

## Minimal examples


```toml
[calculation]
name = "fcc Ag"
kind = "lmto"

[structure]
lattice = "fcc"
alat = 4.0900            # Å

[[species]]
symbol = "Ag"
z = 47
nval = 11
konfig = [5, 5, 4, 4]
occupations = [1.0, 10.0]

[scf]
spin = 1
xc = "vbh"
nk = 12
bz = "tetra"
mixing = 0.3
tolerance = 1.0e-6
max_iterations = 100

[dos]
emin = -0.9
emax = 1.0
npoints = 801
smear = 0.01
nk = 24
```

Other ready-to-run files: `fe_bcc.toml` (ferromagnetic Fe, `spin = 2`),
`fecu_b2.toml` (two classes, Madelung energy, charge transfer), `feco_b2.toml`
(ferromagnetic B2 FeCo).

## `[structure]`

| Key | Default | Meaning |
|---|---|---|
| `lattice` | `fcc` | `fcc`, `bcc`, `sc` (one `[[species]]`); `b2` (CsCl, **exactly two** `[[species]]`); `cell` (general) |
| `alat` | — (required, > 0) | Lattice constant in Å |
| `ratio` | `1.0` | `b2` only: sphere-radius ratio \(r_A/r_B\). Spheres always fill the cell volume. |
| `vectors` | — | `cell` only: three lattice vectors `[[a1],[a2],[a3]]` in `alat` units |
| `site_classes` | `false` | `cell` only: `true` gives **one ASA class per site** (local moments and charges — what an SQS supercell needs); `false` gives one class per species |

### General cells and SQS supercells

```toml
[structure]
lattice = "cell"
alat = 2.86
vectors = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
site_classes = true

[[site]]
position = [0.0, 0.0, 0.0]   # Cartesian, alat units
species = "Fe"

[[site]]
position = [0.5, 0.5, 0.5]
species = "Co"
```

Every `[[site]]` needs `position` (3 numbers) and a `species` equal to the
`symbol` of one `[[species]]`. Special quasirandom structures can be generated
with [`tools/sqs.py`](../guides/sqs.html), which emits exactly these tables.

## `[scf]`

All common keys from the [overview](index.html) plus:

| Key | Default | Meaning |
|---|---|---|
| `nk` | `12` | k-mesh \(n_k^3\) of the SCF; ≥ 4 |
| `bz` | `tetra` | `tetra`: linear tetrahedron method ([[asr_tetra]]); `gauss`: Gaussian smearing |
| `smear` | `0.01` | Gaussian width in Ry (when `bz = "gauss"`) |
| `trace` | `false` | Verbose iteration trace |

!!! warning "Sharp Fermi-level features"
    When \(E_F\) sits on a sharp d peak, Gaussian smearing biases \(N(E_F)\) upward
    (measured: 19.35 electrons at 6 mRy against 19.11 at 1.5 mRy for B2 FeCu).
    Prefer `bz = "tetra"`.

## `[dos]`

| Key | Default | Meaning |
|---|---|---|
| `emin`, `emax` | −1, 1 | Window in Ry (`emax > emin`) |
| `npoints` | `801` | Grid points |
| `smear` | `0.01` | Gaussian broadening of the DOS |
| `nk` | `24` | k-mesh of the DOS (independent of the SCF mesh) |

## `[output]`

`dos_file` (`dos.dat`), `params_file` (`params.dat`), `out_file` (`asr.out`).

## Console output

```text
asr lmto  fcc Ag
  lattice fcc  a(Bohr)=  7.72898  atoms/cell=1  xc=vbh  spin=1
  SCF converged=T  iterations=12  dV= 9.82E-07  E_F= -0.173718 Ry
  class Ag        S(Bohr)=  3.0205  qnet= -0.00000  moment=  0.00000
  wrote ag_fcc_params.dat
  wrote ag_fcc_dos.dat
  wrote asr.out
```

## Validation

Ag reproduces the Stuttgart parameters \((\sqrt\Delta, \gamma, p)\) to ≈ 1 %,
Fe gives 2.2 μ<sub>B</sub>, Ni 0.64 μ<sub>B</sub>, Cu 0, and B2 FeCu occupations agree to
0.01–0.02 e. Details: [Validation](../validation.html).
