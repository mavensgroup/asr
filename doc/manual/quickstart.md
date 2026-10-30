---
title: Quick start
---

# Quick start

This tutorial walks through **four short calculations**—one for each workflow supported by `asr` (`model`, `lmto`, `alloy`, and `sublattice`)—starting from a simple nearest-neighbour tight-binding band and progressing to a self-consistent disordered alloy and a site-disordered Heusler compound.

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

## 1. Haydock Recursion on a Model Band (`kind = "model"`, 1 s)

The simplest calculation constructs a nearest-neighbour tight-binding \(s\)-band on a real-space cluster, tridiagonalizes it with the Lanczos (Haydock) recursion, terminates the continued fraction with a square-root analytic tail, and evaluates the local density of states (LDOS).

### Input (`fcc_ordered.toml`)

```toml
[calculation]
name = "fcc NN s-band"
kind = "model"

[structure]
lattice = "fcc"

[model]
hopping = 1.0

[recursion]
nlevels = 40
terminator = "sqrt-avg2"

[dos]
emin = -14.0
emax = 14.0
eta = 0.05
npoints = 801

[output]
dos_file = "fcc_ordered_dos.dat"
```

### Run & Output

```bash
asr /tmp/run/fcc_ordered.toml
```

```text
asr 0.1.0  fcc NN s-band
input: /tmp/run/fcc_ordered.toml
wrote /tmp/run/fcc_ordered_dos.dat
```

`/tmp/run/fcc_ordered_dos.dat` is a two-column table (`E  LDOS`) for the 12-nearest-neighbour fcc band (`[-12t, +4t]`), computed from 40 Lanczos steps and a square-root terminator. This is the core pipeline—*linear operator → Lanczos coefficients \((a_n, b_n)\) → terminated continued fraction → Green function \(\to\) LDOS*—that every first-principles calculation in `asr` builds upon. Try lowering `nlevels = 10` to observe finite-depth oscillations vs terminator smoothing (see [Theory: recursion](theory/recursion.html)).

---

## 2. Self-Consistent Ordered Crystal: fcc Ag (`kind = "lmto"`, ≈ 9 s)

For `kind = "lmto"`, `asr` runs a complete self-consistent TB-LMTO-ASA calculation from scratch: it solves the free-atom Dirac-like scalar-relativistic equation for the starting potential, computes screened canonical structure constants \(S^\alpha\), integrates the Brillouin zone by the Blöchl tetrahedron method, updates the sphere potential via libxc + Poisson, and mixes to self-consistency with Anderson mixing.

### Input (`ag_fcc.toml`)

```toml
[calculation]
name = "fcc Ag"
kind = "lmto"

[structure]
lattice = "fcc"
alat = 4.0900            # Angstrom (experimental)

[[species]]
symbol = "Ag"
z = 47
nval = 11
konfig = [5, 5, 4, 4]    # valence principal quantum numbers s,p,d,f; lower states are core
occupations = [1.0, 10.0] # neutral-atom valence s and d occupations (starting potential only)

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

[output]
dos_file = "ag_fcc_dos.dat"
params_file = "ag_fcc_params.dat"
```

### Run & Output

```bash
asr /tmp/run/ag_fcc.toml
```

```text
asr lmto  fcc Ag
  lattice fcc  a(Bohr)=  7.72898  atoms/cell=1  xc=vbh  spin=1
  SCF converged=T  iterations=12  dV= 9.82E-07  E_F= -0.173718 Ry
  class Ag        S(Bohr)=  3.0205  qnet= -0.00000  moment=  0.00000
  wrote ag_fcc_params.dat
  wrote ag_fcc_dos.dat
  wrote asr.out
```

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

## 3. Random Binary Alloy: bcc Cu\(_{50}\)Zn\(_{50}\) (`kind = "alloy"`, ≈ 30 s)

For `kind = "alloy"`, `asr` constructs the joint real-space \(\otimes\) configuration-space Hamiltonian in the Mookerjee augmented space, applies exact point-group (\(O_h / D_{4h}\)) symmetry reduction on the disorder configurations, integrates the configuration-averaged Green function along a semicircular complex energy contour at each SCF step, and includes Screened Impurity Model (SIM) Madelung electrostatics.

### Input (`cuzn_alloy.toml`)

```toml
[calculation]
name = "CuZn bcc alloy"
kind = "alloy"

[structure]
lattice = "bcc"
alat = 2.85              # Angstrom

[[species]]
symbol = "Cu"
z = 29
nval = 11
konfig = [4, 4, 3, 4]
occupations = [1.0, 10.0]
concentration = 0.5

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
max_iterations = 60

[asr]
levels = 3               # Lanczos levels (depth); raise for production runs
rank = 2                 # max flipped sites per augmented-space configuration
supercell = 8            # real-space graph size in conventional cells
symmetry = true          # exact D4h symmetry reduction of configuration space
hamiltonian = "energy-independent"
```

### Run & Output

```bash
asr /tmp/run/cuzn_alloy.toml
```

```text
  SCF converged=T  iterations=30  E_F=  0.017089 Ry
  species Cu  qnet=  0.06709  moment=  0.00000
  species Zn  qnet= -0.06709  moment=  0.00000
```

In `/tmp/run/asr.out`:

```text
  recursion      energy-independent          levels=3  max rank=2  supercell=8
  augmented-space nodes 63175   D4h symmetry=T   SIM beta= 1.000   mixing= 0.300
```

Here `63175` symmetry-reduced augmented-space states were generated, and Cu transfers \(0.067\,e\) to Zn under the Screened Impurity Model (`sim_beta = 1.0`). The output file `cuzn_alloy_dos.dat` contains the total and species/\(l\)-projected densities of states (`E  total  A_s A_p A_d  B_s B_p B_d  total_spin1  total_spin2`).

!!! warning "Coarse tutorial settings vs production convergence"
    `cuzn_alloy.toml` sets `levels = 3` so it finishes in ~30 seconds. At `levels = 3`, the continued fraction is truncated very early and partial \(l\)-resolved DOS curves can exhibit unphysical oscillations or small negative values near band edges. For production alloy calculations, use the energy-dependent route (`example/cuzn_alloy_asr_ed.toml` or `example/cuzn_alloy_ed.toml`) with deeper recursion (`levels = 4`, `recursion_steps = 24`–`36`, `rank = 1`, `ed_terminator = "hull"`), and consult the [Convergence guide](guides/convergence.html).

---

## 4. Sublattice Disorder & Multicomponent Alloys: L2\(_1\) Co\(_2\)Fe\(_{0.4}\)Cr\(_{0.6}\)Al (`kind = "sublattice"`)

When a crystal has **multiple inequivalent Wyckoff sublattices**—some ordered and some chemically disordered, or with three or more alloy components—use `kind = "sublattice"`. Each `[[sublattice]]` block specifies its fractional coordinates on the decorated Bravais graph, the list of occupying `species`, and (if disordered) their `concentrations`.

Look at `example/co2fe04cr06al_l21.toml`, which models the partially disordered L2\(_1\) Heusler alloy Co\(_2\)Fe\(_{0.4}\)Cr\(_{0.6}\)Al:

```toml
[calculation]
name = "co2fe04cr06al L2_1, kind sublattice"
kind = "sublattice"

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

Here Co (`8c`) and Al (`4a`) occupy ordered sublattices while Fe and Cr share the `4b` sublattice randomly in a 40:60 ratio. Running:

```bash
asr /tmp/run/co2fe04cr06al_l21.toml
```

converges the spin-polarized multi-sublattice SCF with sublattice-averaged Ewald Madelung potentials plus intra-sublattice SIM screening, yielding species-resolved magnetic moments (\(\mu_{\text{Co}} \approx 0.94\,\mu_B\), \(\mu_{\text{Fe}} \approx 2.75\,\mu_B\), \(\mu_{\text{Cr}} \approx 1.24\,\mu_B\), \(\mu_{\text{Al}} \approx -0.08\,\mu_B\)). Note that this input uses production recursion settings (`levels = 4`, `recursion_steps = 24`, `spin = 2`) and takes several minutes on a 16-core workstation.

---

## 5. Where to Go Next

| Goal | Manual Chapter |
|---|---|
| Look up every TOML section and key | [Input reference](input/index.html) (`model`, `lmto`, `alloy`, `sublattice`) |
| Understand every table and column in `asr.out` and `*_dos.dat` | [Output files](output.html) |
| Choose `levels`, `recursion_steps`, `rank`, and terminators for publication-grade runs | [Convergence guide](guides/convergence.html) |
| Converge ferromagnetic, ferrimagnetic, and half-metallic alloys | [Magnetism guide](guides/magnetism.html) |
| Cross-check ASR results against Special Quasirandom Structures (SQS) | [SQS cross-checks](guides/sqs.html) |
| Understand the augmented-space formalism and TB-LMTO-ASA equations | [Theory](theory/index.html) |
| Diagnose input errors or SCF convergence issues | [Troubleshooting](guides/troubleshooting.html) |
