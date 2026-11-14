---
title: Quick start
---

# Quick start

This tutorial walks through **four calculations** with the one input format of `asr`: an ordered crystal by the exact k-space reference, the same kind of crystal by the augmented-space recursion, a random binary alloy, and a site-disordered Heusler compound. The Hamiltonian is the same in all four (TB-LMTO-ASA, the only one implemented; see [Hamiltonians](hamiltonian/index.html)). Only the structure and the `method` change.

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

## 1. Ordered Crystal by k-space: fcc Ag (`method = "kspace"`, ≈ 9 s)

With `method = "kspace"`, `asr` runs a complete self-consistent TB-LMTO-ASA calculation from scratch: it solves the free-atom scalar-relativistic equation for the starting potential, computes screened canonical structure constants \(S^\alpha\), integrates the Brillouin zone by the Blöchl tetrahedron method, updates the sphere potential via libxc + Poisson, and mixes to self-consistency with Anderson mixing. This is the exact reference for ordered structures; the recursion in the next section solves the same model.

### Input (`ag_fcc.toml`)

```toml
# fcc Ag, k-space TB-LMTO-ASA (the exact reference of the ordered crystal).
[calculation]
name = "fcc Ag"
method = "kspace"

[structure]
lattice = "fcc"
alat = 4.0900            # Angstrom (experimental)
centering = "F"          # one sublattice: the fcc lattice itself

[[sublattice]]
name = "fcc"
positions = [[0.0, 0.0, 0.0]]
species = ["Ag"]

[[species]]
symbol = "Ag"
z = 47
nval = 11
konfig = [5, 5, 4, 4]
occupations = [1.0, 10.0]

[scf]
spin = 1
xc = "vbh"
mixing = 0.3
tolerance = 1.0e-6
max_iterations = 100

[kspace]
nk = 12
bz = "tetra"
dos_nk = 24

[dos]
emin = -0.9
emax = 1.0
npoints = 801
eta = 0.01

[output]
dos_file = "ag_fcc_dos.dat"
params_file = "ag_fcc_params.dat"
```

### Run & Output

```bash
asr /tmp/run/ag_fcc.toml
```

The terminal prints a short summary (SCF convergence, \(E_F\), per-class charge and moment) and the files written.

Inspect `/tmp/run/asr.out`. The converged potential parameters and \(l\)-resolved charges are summarized in:

```text
=== converged results =====================================================
 Fermi energy                   -0.17371795 Ry
 class 1  Ag   Z=  47.0  valence=11  S=   3.02046 Bohr   atoms/cell=1
   spin  l      E_nu         C        sqrt(Delta)       p         gamma        occ(q0)
      1  s   -0.533991   -0.449414     0.364910    0.081764    0.424921     0.68472
      1  p   -0.437422    0.425053    -0.365166    0.036272    0.115341     0.66416
      1  d   -0.485924   -0.488864    -0.104418    2.061327    0.005199     9.50980
```

These `s/p/d` occupations (`0.685 / 0.664 / 9.510`) and potential parameters (\(\sqrt{\Delta}\), \(\gamma\), \(p\)) agree with the Stuttgart TB-LMTO-ASA reference to \(\approx 1\%\) without importing any external parameters (see [Validation](validation.html)).

!!! tip "Plot the l-resolved DOS"
    Plot the Ag \(4d\) band from `/tmp/run/ag_fcc_dos.dat` using `gnuplot`:

        gnuplot -p -e "plot '/tmp/run/ag_fcc_dos.dat' u 1:4 w l lw 2 t 'Ag 4d'"

    Columns are `E  s  p  d` (states/Ry/atom per spin channel), and `E_F` is printed in the second header comment line.

!!! tip "Spin-polarized and multi-class ordered crystals"
    Try running ferromagnetic bcc Fe (`asr /tmp/run/fe_bcc.toml`, which converges to \(m \approx 2.21\,\mu_B\)) and ordered B2 FeCu (`asr /tmp/run/fecu_b2.toml`, which demonstrates a two-class crystal with Ewald Madelung electrostatics and inter-sublattice charge transfer).

!!! note "Harmless IEEE floating-point notices on exit"
    Intel `ifx` binaries may print `forrtl: warning (525): IEEE_DIVIDE_BY_ZERO is signaling` (or `527: IEEE_OVERFLOW`) upon program termination. These are standard Fortran runtime summaries indicating that an IEEE floating-point status flag was raised during execution; the calculation still exits with status `0` and all numerical results are unaffected.

---

## 2. The Same Structure by Recursion: B2 FeCo (default `method = "asr"`)

Leave out `method` and the same structure block runs by the energy-dependent augmented-space recursion, the default solver. For an ordered crystal there is no disorder to average, and the recursion reproduces the k-space result to the accuracy of the chain length. Compare `example/feco_b2_asr.toml` with `example/feco_b2_kspace.toml`: the `[structure]`, `[[sublattice]]` and `[[species]]` blocks are identical, and only `method`, `[kspace]` (k mesh) and `[asr]` (recursion) differ. This pairing is the basic validation of the recursion ([Validation](validation.html)).

```toml
# Ordered B2 FeCo, ferromagnetic, by the energy-dependent recursion (default method).  k-space: example/feco_b2_kspace.toml.
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
