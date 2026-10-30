---
title: Input reference
ordered_subpage: model.md
                 lmto.md
                 alloy.md
                 sublattice.md
---

# Input reference

An `asr` run is described by **one TOML file**. The file is parsed by
[toml-f](https://github.com/toml-f/toml-f) and validated by the readers in
[[asr_input]], [[asr_input_lmto]], [[asr_alloy_driver]] and [[asr_class_driver]].

## Rules that apply to every kind

1. **Strict keys.** Unknown top-level sections and unknown keys inside a section
   are errors (`asr: input error: unknown key "…" in [scf]`). A typo can never
   silently fall back to a default.
2. **Missing sections are empty sections.** Every key has a default, so an absent
   `[dos]` or `[output]` simply uses the defaults listed in the tables.
3. **Bad values are errors** and are caught *before* any computation starts
   (for example `[scf] spin must be 1 or 2`, `[dos] emax must exceed emin`).
4. **Paths** in `[output]` (and `save_potential` / `start_potential`) are
   relative to the directory of the *input file*, not to the current directory.
5. **Units.** `alat` is in Å. Energies (`emin`, `emax`, `eta`, `smear`, `fermi_kt`)
   are in Ry. Positions in `[[site]]` are Cartesian in units of the lattice
   constant. See [Units and conventions](../appendix/conventions.html).

## Choosing the kind

```toml
[calculation]
name = "my run"          # free text, echoed in output headers
kind = "lmto"            # model | lmto | alloy | sublattice
```

`kind` defaults to `model` when absent.


| You have … | Use |
|---|---|
| a toy s-band to test the recursion/terminator | [`model`](model.html) |
| an ordered crystal (any lattice and basis), or an SQS supercell | [`lmto`](lmto.html) |
| a random binary alloy on fcc/bcc, or the ordered B2 counterpart | [`alloy`](alloy.html) |
| an ordered compound, partial sublattice disorder, or ≥ 3 components | [`sublattice`](sublattice.html) |

## Blocks shared by the LMTO-based kinds

The blocks `[[species]]`, `[scf]`, `[dos]` and `[output]` appear in `lmto`,
`alloy` and `sublattice` files with minor differences, listed on each page.

### `[[species]]`

One table per chemical species. Everything the code needs to start a self-consistent
calculation is given here — nothing is tabulated.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `symbol` | string | `"X"` | Label used in reports and to refer to the species from `[[site]]` / `[[sublattice]]`. Must be unique. |
| `z` | number | — (required, > 0) | Atomic number (nuclear charge) |
| `nval` | integer | — (required, > 0) | Number of valence electrons |
| `konfig` | `[ns, np, nd, nf]` | `[4,4,3,4]` | Principal quantum number of the **lowest valence** state of each \(l\). States below are treated as core. |
| `occupations` | `[n_s, n_d]` | `[2.0, 0.0]` | Neutral-atom valence s and d occupations. Used **only** to build the starting potential. The rest of `nval` goes to p. |
| `spin_sign` | ±1 | `1` | Sign of the starting spin polarisation (`zeta0`) of this species. `-1` gives an antiparallel start. Any value is reduced to its sign. |
| `concentration` | number | — | `alloy` only, first species only. See [alloy](alloy.html). |

!!! tip "Choosing `konfig` and `occupations`"
    For a 3d transition metal use `konfig = [4, 4, 3, 4]` and `occupations =
    [2.0, 6.0]` (Fe: 4s²3d⁶). For a 4d metal such as Ag use `[5, 5, 4, 4]` and
    `[1.0, 10.0]`. For sp metals with no d valence, give the lowest *d* shell the
    code may treat as valence, e.g. Al: `konfig = [3, 3, 3, 4]`, `occupations =
    [2.0, 0.0]`.

### `[scf]` (common keys)

| Key | Type | Default | Meaning |
|---|---|---|---|
| `spin` | 1 or 2 | `1` | 1: non-magnetic, 2: collinear spin-polarised |
| `xc` | string | `"vbh"` | `vbh` (von Barth–Hedin LDA), `pw`, `pwmod`, `vwn`, `pbe` (GGA). Any other name aborts in [[asr_xc]]. |
| `mixing` | number | 0.3 | Linear/Anderson mixing parameter β ([[asr_mixing]]) |
| `tolerance` | number | `1e-6` (`lmto`), `1e-5` (others) | Convergence of the weighted potential change, in Ry |
| `max_iterations` | integer | 100 (`lmto`), 60 (others) | SCF iteration limit |
| `zeta0` | number | 0.4 (`lmto`), 0.1 (others) | Initial spin polarisation of the atomic start (`spin = 2`) |
| `relativistic` | bool | `true` | Scalar-relativistic radial equation (Koelling–Harmon). `false` is Schrödinger. |
| `sigma` | number | 0.7 | Hard-core screening radius of the screened structure constants (`lmto`, `alloy`) |

!!! warning "Magnetic starts"
    A large `zeta0` is applied to the whole atomic density *including the core* and
    can make the SCF run away (observed for Fe with 0.4 in the alloy route). The
    alloy/sublattice default of 0.1 is deliberate. See the
    [magnetism guide](../guides/magnetism.html).

### `[output]`

| Key | Default | Written by |
|---|---|---|
| `dos_file` | `dos.dat` | all kinds |
| `params_file` | `params.dat` | `lmto` only |
| `out_file` | `asr.out` | `lmto`, `alloy`, `sublattice` |

See [Output files](../output.html) for the formats.
