---
title: kind = "model"
---

# `kind = "model"` — ordered tight-binding band

The simplest workflow: a nearest-neighbour single-orbital band on an sc, bcc or
fcc lattice. The local density of states at the origin site is computed by
**recursion** on a finite cluster ([[asr_lattice_toy]]) closed with a
square-root terminator ([[asr_greens]]). It exists to test and demonstrate the
recursion core and is the fastest way to see a terminated continued fraction at work.

```toml
[calculation]
name = "fcc NN s-band"
kind = "model"

[structure]
lattice = "fcc"          # sc | bcc | fcc

[model]
hopping = 1.0            # t, in the energy unit of the output

[recursion]
nlevels = 40             # number of (a_n, b_n) levels
terminator = "sqrt-avg2" # sqrt-last | sqrt-avg2

[dos]
emin = -14.0
emax = 14.0
eta = 0.05
npoints = 801

[output]
dos_file = "fcc_ordered_dos.dat"
```

## Keys

| Section | Key | Default | Meaning |
|---|---|---|---|
| `[calculation]` | `name` | `asr` | Label |
| | `kind` | `model` | Must be `model` |
| `[structure]` | `lattice` | `sc` | `sc`, `bcc` or `fcc` |
| `[model]` | `hopping` | `1.0` | Nearest-neighbour hopping *t*. The band is \(\varepsilon(\mathbf k)\) of the lattice, so the half-width is \(zt\) (z = 6, 8, 12) |
| `[recursion]` | `nlevels` | `40` | Recursion depth, ≥ 2 |
| | `terminator` | `sqrt-avg2` | Tail from the last coefficients: `sqrt-last` uses \((a_n, b_n)\) of the last level; `sqrt-avg2` averages the last two |
| `[dos]` | `emin`, `emax` | −8, 8 | Energy window; `emax > emin` |
| | `eta` | `0.05` | Imaginary part of \(z = E + i\eta\); must be > 0 |
| | `npoints` | `801` | Grid points, ≥ 2 |
| `[output]` | `dos_file` | `dos.dat` | Output file |

For an ordered NN band the first recursion coefficients are known exactly:
\(a_1 = 0\) and \(b_1^2 = z\) (the coordination number). The test-suite checks this.

## Output

Two comment lines followed by two columns `E  LDOS(E)`:

```text
# fcc NN s-band
# nlevels=40 terminator=sqrt-avg2 eta= 0.0500
 -1.40000000E+01  1.28233610E-04
 …
```

!!! note "Accuracy"
    Accuracy is limited by *depth*, not by the tail. Measurements and the open
    bcc van Hove problem are on the [Validation](../validation.html) page.
    The `term_bench` program ([[term_bench]]) compares every terminator against
    exact k-space results.
