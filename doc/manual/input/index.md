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
| `[asr]` | recursion, contour, route, diagnostics | the method, except `route` | [Recursion](asr.html) |
| `[[species]]` | atomic data from which the parameters are generated | **yes** | [Hamiltonian](hamiltonian.html) |
| `[scf]` | xc, spin, mixing, convergence, screening model | partly | [Hamiltonian](hamiltonian.html) |
| `[kspace]` | Brillouin-zone reference | **yes** | [k-space](kspace.html) |
| `[dos]`, `[output]` | energy window, files | no | below |

Only [TB-LMTO-ASA](../hamiltonian/tb_lmto_asa.html) is implemented, so there is
no key to select a Hamiltonian ([Hamiltonians](../hamiltonian/index.html)).

## `[calculation]`

| Key | Default | Meaning |
|---|---|---|
| `name` | `asr` | Free text, echoed in output headers |
| `method` | `asr` | `asr`: augmented-space recursion, any decorated structure, ordered or disordered (the default, also for an ordered crystal). `kspace`: Brillouin-zone solution of the same Hamiltonian, ordered structures only; the exact reference |

The same structure block runs under both methods: delete `method = "kspace"` and
`[kspace]` from `example/fe_bcc.toml` and the recursion runs it.

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
