---
title: Troubleshooting
---

# Troubleshooting

## Input errors (exit code 1)

Messages have the form `asr: input error: <message>`. Frequent ones:

| Message | Cause and fix |
|---|---|
| `cannot read '<file>': …` | File missing or invalid TOML |
| `unknown key "k" in [scf]` / `in top level` | Typo, or a key that belongs to the other `method` (`[asr]` needs `method = "asr"`, `[kspace]` needs `"kspace"`). Check the [reference](../input/index.html). |
| `[structure] alat (Angstrom) must be > 0` | `alat` missing |
| `[[species]] needs z and nval` | Both are required |
| `[calculation] kind is no longer read` | Old input. Use the one format: delete `kind`, set `method` if you want k-space ([Input](../input/index.html)) |
| `lattice "cell" needs [structure] vectors` / `[[site]] entries` | Missing blocks |
| `[[site]] species "X" is not a [[species]] symbol` | Spelling |
| `method = "kspace" needs an ordered structure` | Disorder runs with `method = "asr"` |
| `lattice = "cell" (explicit sites) runs with … method = "kspace" only` | Explicit-site cells are k-space only |
| `concentrations must sum to 1` | |
| `[asr] recursion_steps must be 0 or >= levels` | |
| `[[sublattice]] … concentrations must be > 0 and sum to 1` | |
| `[[sublattice]] with several species needs concentrations` | |

!!! note "Unknown `xc`"
    An unrecognised `xc` name is *not* caught by the reader; it aborts with
    `xc_init: unknown functional` from [[asr_xc]]. Use `vbh`, `pw`, `pwmod`,
    `vwn` or `pbe`.

## Exit code 3: SCF not converged

The report ends `status … SCF NOT CONVERGED` (see the `SCF history`). Try, in
order:

1. Lower `mixing` (0.3 → 0.2 → 0.1).
2. Magnetic runs: lower `zeta0` (0.4 → 0.1). See [Magnetism](magnetism.html).
3. Raise `max_iterations`.
4. Energy-dependent routes: set `fermi_kt` (e.g. 0.001).
5. Restart from a good potential with `save_potential` / `start_potential`.
6. If \(E_F\) sits on a very sharp d peak, use `bz = "tetra"` in `lmto`.

## Symptoms

| Symptom | Likely cause |
|---|---|
| Negative l-resolved DOS | Recursion too shallow → [Convergence](convergence.html) |
| Charge transfer changes sign or oscillates over many iterations | Normal in the first 15–20 iterations of the alloy SCF (see the shipped CuZn example, which reaches 0.067 e only after ~20 iterations); persistent oscillation → `sim_beta`, `mixing` |
| Fermi level drifts between iterations | Sharp \(N(E)\); use `fermi_kt`, tetrahedra |
| \(E_F\) differs from another code by ≈ 10 mRy | `charges = "sphere"` at order 3; compare using `lowdin` |
| s,p occupations 3–5 % above a Stuttgart reference | Known open item: Stuttgart moments include \(o\) and \(p\) terms, ours are plain orthonormal \(|u|^2\) (untested hypothesis) |
| Unexplained link/runtime errors after editing a type or procedure signature | fpm kept stale objects: `rm -rf build/ifx_*` and rebuild |
| `IEEE_DIVIDE_BY_ZERO is signaling` on exit | Runtime flag notice; see [Quick start](../quickstart.html) |

## Getting help from the program

* `asr.out` records the full input and effective settings (screening \(\alpha\),
  node counts, depth). Attach it to bug reports.
* Diagnostics for the energy-dependent route: `tail_scan`, `kkr_check`.
* Keep `trial/` for scratch runs; it is ignored by git.
