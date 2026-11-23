---
title: Input reference
ordered_subpage: structure.md
                 hamiltonian.md
                 asr.md
                 kspace.md
---

# Input reference

An `asr` run is described by **one TOML file** in **one format**. The file is
parsed by [toml-f](https://github.com/toml-f/toml-f) and validated by the reader
in [[asr_driver]]. The old `[calculation] kind = ...` input (`model`, `lmto`,
`alloy`, `sublattice`) is gone; the reader rejects `kind` with a message.

## Rules that apply to every input

1. **Strict keys.** Unknown top-level sections and unknown keys inside a section
   are errors (`asr: input error: unknown key "…" in [scf]`). A typo can never
   silently fall back to a default.
2. **Missing sections are empty sections.** Every key has a default, so an absent
   `[dos]` or `[output]` uses the defaults listed in the tables.
3. **Bad values are errors** and are caught *before* any computation starts.
4. **Paths** in `[output]`, `save_potential` and `start_potential` are relative to
   the directory of the *input file*, not to the current directory.
5. **Units.** `alat` is in Å. Energies (`emin`, `emax`, `eta`, `smear`,
   `fermi_kt`) are in Ry. Positions in `[[site]]` are Cartesian in units of the
   lattice constant. See [Units and conventions](../appendix/conventions.html).

## Which section is for what

| Section | Content | Depends on the Hamiltonian? | Page |
|---|---|---|---|
| `[calculation]` | `name`, `method` | no | below |
| `[structure]`, `[[sublattice]]`, `[[site]]` | geometry, sublattices, species and concentrations | no | [Structure](structure.html) |
| `[[sro]]` | first-shell short-range order | no (the method) | [Recursion](asr.html) |
| `[epi]` | effective pair interactions after the SCF | no (the method) | [Recursion](asr.html) |
| `[asr]` | recursion, contour, route, diagnostics | the method, except `route` | [Recursion](asr.html) |
| `[[species]]` | atomic data from which the parameters are generated | **yes** | [Hamiltonian](hamiltonian.html) |
| `[scf]` | xc, spin, mixing, convergence, screening model | partly | [Hamiltonian](hamiltonian.html) |
| `[kspace]` | Brillouin-zone reference | **yes** | [k-space](kspace.html) |
| `[dos]`, `[output]` | energy window, files | no | below |
| `[spectral]` | Bloch spectral function A(k, E) after the SCF | no | below |

Only [TB-LMTO-ASA](../hamiltonian/tb_lmto_asa.html) is implemented, so there is
no key to select a Hamiltonian ([Hamiltonians](../hamiltonian/index.html)).

## `[calculation]`

| Key | Default | Meaning |
|---|---|---|
| `name` | `asr` | Free text, echoed in output headers |
| `method` | `asr` | `asr`: augmented-space recursion, any decorated structure, ordered or disordered (the default, also for an ordered crystal). `kspace`: Brillouin-zone solution of the same Hamiltonian, ordered structures only; an internal test tool, not a production route |

The same structure block runs under both methods.  The test inputs in `data/kspace/` show the `[kspace]` section.

## `[dos]`

| Key | Default | Meaning |
|---|---|---|
| `emin`, `emax` | −1, 1 | Energy window in Ry; `emax > emin` |
| `npoints` | `401` (`801` for `kspace`) | Grid points |
| `eta` | `0.02` (`0.01` for `kspace`) | Recursion: imaginary part of the energy. k-space: Gaussian smearing |

## `[output]`

| Key | Default | Meaning |
|---|---|---|
| `dos_file` | `dos.dat` | DOS table |
| `out_file` | `asr.out` | Report |
| `params_file` | `params.dat` | Potential parameters, `method = "kspace"` only |

See [Output files](../output.html).

## `[spectral]`

The configuration-averaged Bloch spectral function \(A(\mathbf k,E)\), the band structure of the disordered alloy, at the
final potential.  Route `edc` only.  It is unfolded on the Brillouin zone of the graph lattice (bcc, fcc, sc) and given per
graph site and per spin.  Every sublattice must be one lattice of translations (random alloys, B2, L2\(_1\); not L1\(_2\)).

| Key | Default | Meaning |
|---|---|---|
| `kpath` | required | Corners of the path, Cartesian, units \(2\pi/a\) of the graph lattice, at least two |
| `labels` | empty | One label per corner |
| `nk` | `30` | Points per segment |
| `emin`, `emax` | those of `[dos]` | Energy window (Ry) |
| `npoints` | `201` | Energies |
| `eta` | `0.01` | Imaginary part of the energy (Ry) |
| `rmax` | `2.0` | Largest displacement in the Fourier sum, units of \(a\); below `supercell / 2` |
| `window` | `"ball"` | Weight of each displacement in the Fourier sum: `"ball"` or `"sharp"` (below) |
| `steps` | `recursion_steps` | Length of the spectral chains |
| `file` | `spectral.dat` | Output table, one block per k point |

Cost: one chain per spin, component, orbital (9) and energy on the operator reduced with D2h, about twice the nodes of the
SCF operator.  `asr.out` lists the checks: the largest graph distance of a target against the interior of the node set, the
difference between the two origins of a pair, the change when the chains are cut four steps earlier, the number of points with
\(A<0\) (chains that are too short, or the ringing of the sharp window), and the decay of \(\mathrm{Im}\,G\) with the shell.

The window.  The Fourier sum stops at `rmax`.  With `window = "sharp"` every displacement up to `rmax` has weight 1.  The sum
then rings: for ordered fcc Cu at `eta` 0.01 Ry, \(A\) falls to \(-30\,\%\) of its peak.  The default `"ball"` weights a
displacement \(R\) by \(w(R) = (1-R/r_{\max})^2\,(1+R/(2r_{\max}))\), the overlap of two balls of diameter \(r_{\max}\).  Its Fourier
transform is non-negative.  So the result is the exact \(A\) smoothed in \(\mathbf k\) with a positive kernel of width about
\(1/r_{\max}\): it cannot ring below zero, and its BZ average and its energy integral at each \(\mathbf k\) are unchanged.  The
price is resolution in \(\mathbf k\).  Raise `rmax` (and `supercell`) for sharper bands.

The energy integral of \(A\) at a k point counts the bands inside the window, at most 9 per spin.  The file sums both spins for
`spin = 1`.

## Skeleton

```toml
[calculation]
name = "my run"

[structure]            # geometry: Hamiltonian-independent
lattice = "bcc"
alat = 2.86

[[sublattice]]
name = "a"
positions = [[0.0, 0.0, 0.0]]
species = ["Fe"]

[[species]]            # atomic data: Hamiltonian-specific (TB-LMTO-ASA)
symbol = "Fe"
z = 26
nval = 8
konfig = [4, 4, 3, 4]
occupations = [2.0, 6.0]

[scf]
spin = 2

[asr]                  # the method
levels = 4

[dos]
emin = -0.9
emax = 0.8
```
