---
title: Quick start
---

# Quick start

This tutorial walks through **four calculations** with the one input format of `asr`, all by the augmented-space recursion: an ordered compound, a ferromagnetic ordered compound, a random binary alloy, and a site-disordered Heusler compound. The Hamiltonian is the same in all four (TB-LMTO-ASA, the only one implemented; see [Hamiltonians](hamiltonian/index.html)). Only the structure changes.

---

## 0. Setup and Scratch Directory

Make sure you have installed the prerequisites and built `asr` (see the [Installation](installation.html) guide):

```bash
git clone https://github.com/mavensgroup/asr.git
cd asr
source env.sh
asr_build
```

All example inputs live in `example/`. Because `asr` always writes output files (`asr.out`, `*_dos.dat`, `*_params.dat`, saved potentials) **into the same directory as the input file**, copy the examples into a scratch directory first to keep your working tree clean:

```bash
mkdir -p /tmp/run
cp example/*.toml /tmp/run/
export OMP_NUM_THREADS=16
export OMP_STACKSIZE=256M
```

In the commands below, `asr` refers to the installed binary (`asr_run -- <file>` works identically from the repository root).

---

## 1. Ordered Compound: B2 CuZn (≈ 2 min)

An ordered crystal needs no special method. Every sublattice holds one species, so there is no disorder to average, and the recursion runs on the ordered real-space Hamiltonian. `asr` runs a complete self-consistent TB-LMTO-ASA calculation from scratch: it solves the free-atom scalar-relativistic equation for the starting potential, computes the screened structure constants \(S^\alpha\), integrates the recursion Green function on a complex-energy contour, updates the sphere potentials via libxc + Poisson, and mixes to self-consistency with Anderson mixing.

### Input (`cuzn_b2.toml`)

```toml
# Ordered B2 CuZn by the energy-dependent recursion.
[calculation]
name = "B2 CuZn"

[structure]
lattice = "bcc"
alat = 2.92055
period = 1
centering = "P"          # two sublattices: cube corner and cube centre

[[sublattice]]
name = "a"
positions = [[0.0, 0.0, 0.0]]
species = ["Cu"]

[[sublattice]]
name = "b"
positions = [[0.5, 0.5, 0.5]]
species = ["Zn"]

[[species]]
symbol = "Cu"
z = 29
nval = 11
konfig = [4, 4, 3, 4]
occupations = [1.0, 10.0]

[[species]]
symbol = "Zn"
z = 30
nval = 12
konfig = [4, 4, 3, 4]
occupations = [2.0, 10.0]

[scf]
spin = 1
xc = "vbh"
mixing = 0.3
tolerance = 1.0e-5
max_iterations = 80

[asr]
recursion_steps = 24
supercell = 12
contour_points = 32

[dos]
emin = -0.9
emax = 0.6
npoints = 401
eta = 0.01

[output]
dos_file = "cuzn_b2_dos.dat"
```

### Run & Output

```bash
asr /tmp/run/cuzn_b2.toml
```

The terminal prints a short summary (SCF convergence, \(E_F\), per-class charge and moment) and the files written. The run converges in 37 iterations, in about 2 min on 8 threads.

Inspect `/tmp/run/asr.out`. The converged potential parameters and \(l\)-resolved charges are summarized in:

```text
=== converged results (per component atom; cell = decoration cell) ========
 Fermi energy                   -0.04550604 Ry
 class 1  Cu@a   Z=  29.0  valence=11  S=   2.71742 Bohr   atoms/cell=1
   qnet (electrons in sphere - Z)      0.13477 e
   spin  l      E_nu         C        sqrt(Delta)       p         gamma        occ(q0)
      1  s   -0.527072   -0.441189    -0.398659    0.063845    0.421305     0.72072
      1  p   -0.374622    0.520164     0.388781    0.031766    0.111313     0.83541
      1  d   -0.311871   -0.305992     0.091900    2.931075   -0.003933     9.57863
 class 2  Zn@b   Z=  30.0  valence=12  S=   2.71742 Bohr   atoms/cell=1
   qnet (electrons in sphere - Z)     -0.13477 e
```

Charge flows from Zn to Cu (0.135 e per atom), and the cell stays neutral.

!!! tip "Plot the DOS"
    `/tmp/run/cuzn_b2_dos.dat` holds the total and class/\(l\)-projected densities of states. The columns are described in [Output files](output.html).

!!! note "Harmless IEEE floating-point notices on exit"
    Intel `ifx` binaries may print `forrtl: warning (525): IEEE_DIVIDE_BY_ZERO is signaling` (or `527: IEEE_OVERFLOW`) upon program termination. These are standard Fortran runtime summaries indicating that an IEEE floating-point status flag was raised during execution; the calculation still exits with status `0` and all numerical results are unaffected.

---

## 2. Ferromagnetic Ordered Compound: B2 FeCo (`spin = 2`, ≈ 2.5 min)

`spin = 2` switches on spin polarisation. Everything else works as in the previous section. The input is `example/feco_b2_asr.toml`:

```toml
# Ordered B2 FeCo, ferromagnetic, by the energy-dependent recursion (default method).
[calculation]
name = "B2 FeCo FM, recursion"

[structure]
lattice = "bcc"
alat = 2.857             # Angstrom
period = 1
centering = "P"          # two sublattices: cube corner and cube centre

[[sublattice]]
name = "a"
positions = [[0.0, 0.0, 0.0]]
species = ["Fe"]

[[sublattice]]
name = "b"
positions = [[0.5, 0.5, 0.5]]
species = ["Co"]

[[species]]
symbol = "Fe"
z = 26
nval = 8
konfig = [4, 4, 3, 4]
occupations = [2.0, 6.0]

[[species]]
symbol = "Co"
z = 27
nval = 9
konfig = [4, 4, 3, 4]
occupations = [2.0, 7.0]

[scf]
spin = 2
xc = "vbh"
mixing = 0.2
zeta0 = 0.1
tolerance = 1.0e-5
max_iterations = 200

[asr]
recursion_steps = 24     # Lanczos steps (ordered: one component per sublattice)
supercell = 12
contour_points = 32

[dos]
emin = -0.9
emax = 0.6
npoints = 601
eta = 0.01

[output]
dos_file = "feco_b2_asr_dos.dat"
```

---

## 3. Random Binary Alloy: bcc Cu\(_{50}\)Zn\(_{50}\)

The only change against an ordered crystal is a sublattice with several species and their concentrations. `asr` then averages the Green function exactly over the disorder in the Mookerjee augmented space, with point-group (\(O_h / D_{4h}\)) symmetry reduction of the configuration space, a complex-energy contour at each SCF step, and Screened Impurity Model (SIM) electrostatics.

### Input (`cuzn_alloy.toml`)

```toml
# Random bcc Cu0.5Zn0.5 (A2) by the energy-dependent augmented-space recursion: one sublattice with two species.
[calculation]
name = "CuZn bcc alloy"

[structure]
lattice = "bcc"
alat = 2.85              # Angstrom
centering = "I"

[[sublattice]]
name = "bcc"
positions = [[0.0, 0.0, 0.0]]
species = ["Cu", "Zn"]
concentrations = [0.5, 0.5]

[[species]]
symbol = "Cu"
z = 29
nval = 11
konfig = [4, 4, 3, 4]
occupations = [1.0, 10.0]

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
mixing = 0.3
max_iterations = 60

[asr]
levels = 4               # augmented-space generations
recursion_steps = 24     # Lanczos steps; convergence must be checked
rank = 1                 # max flipped sites per augmented-space state
supercell = 8            # graph size in conventional cells
contour_points = 32

[dos]
emin = -0.9
emax = 0.7
npoints = 401
eta = 0.02

[output]
dos_file = "cuzn_alloy_dos.dat"
```

### Run & Output

```bash
asr /tmp/run/cuzn_alloy.toml
```

The terminal prints the SCF history and, per species, the net charge and moment. `asr.out` records the recursion settings, the number of augmented-space nodes and the SIM parameter. `cuzn_alloy_dos.dat` holds the total and species/\(l\)-projected densities of states ([Output files](output.html)).

!!! warning "Check the recursion depth"
    `levels`, `recursion_steps` and `rank` control the accuracy. Too few steps truncate the continued fraction early and can produce oscillations or small negative values near band edges. Follow the [Convergence guide](guides/convergence.html) before trusting a result.

---

## 4. Sublattice Disorder & Multicomponent Alloys: L2\(_1\) Co\(_2\)Fe\(_{0.4}\)Cr\(_{0.6}\)Al

A crystal with **several inequivalent Wyckoff sublattices**, some ordered and some disordered, or with three or more alloy components, uses the same blocks. Each `[[sublattice]]` gives its fractional coordinates on the decorated Bravais graph, the occupying `species` and, if disordered, their `concentrations`. See `example/co2fe04cr06al_l21.toml`, which models the partially disordered L2\(_1\) Heusler alloy Co\(_2\)Fe\(_{0.4}\)Cr\(_{0.6}\)Al:

```toml
[structure]
lattice = "bcc"
alat = 2.8635
period = 2
centering = "F"

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
```

Co (`8c`) and Al (`4a`) occupy ordered sublattices while Fe and Cr share `4b` randomly in a 40:60 ratio. Running:

```bash
asr /tmp/run/co2fe04cr06al_l21.toml
```

converges the spin-polarised multi-sublattice SCF with species-resolved magnetic moments. The input uses production recursion settings (`levels = 4`, `recursion_steps = 24`, `spin = 2`) and takes several minutes on a 16-core workstation.

For the species-dependent first-shell correlations of a high-entropy alloy see `example/fe66ni14cr20_fcc_sro.toml` and the `[[sro]]` section of the [input reference](input/asr.html).

---

## 5. Where to Go Next

| Goal | Manual Chapter |
|---|---|
| Look up every TOML section and key | [Input reference](input/index.html) |
| See which part of the input and code belongs to the Hamiltonian | [Hamiltonians](hamiltonian/index.html) |
| Understand every table and column in `asr.out` and `*_dos.dat` | [Output files](output.html) |
| Choose `levels`, `recursion_steps`, `rank`, and terminators for publication-grade runs | [Convergence guide](guides/convergence.html) |
| Converge ferromagnetic, ferrimagnetic, and half-metallic alloys | [Magnetism guide](guides/magnetism.html) |
| Cross-check ASR results against Special Quasirandom Structures (SQS) | [SQS cross-checks](guides/sqs.html) |
| Understand the augmented-space formalism and the TB-LMTO-ASA equations | [Theory](theory/index.html), [TB-LMTO-ASA](hamiltonian/tb_lmto_asa.html) |
| Diagnose input errors or SCF convergence issues | [Troubleshooting](guides/troubleshooting.html) |
