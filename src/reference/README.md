# Reference layer

The program (`app/main.f90`) does not use any module in this directory. They are the earlier generations of the
method and independent numerics. The tests and the probes use them to cross-check the production code. A test may
check that the class operator with the `edc` route reproduces the binary route. Or it may check that the recursion
equals an exact configuration average or a dense solve.

| Module | What it is | Production counterpart |
|---|---|---|
| `asr_alloy_scf` | binary self-consistent alloy | `asr_class_scf` |
| `asr_augorb` | binary multi-orbital augmented-space operator | `asr_augcls` |
| `asr_augsym` | binary operator reduced by the point group | `asr_augcls` (with `asr_symmetry`) |
| `asr_augn` | n-species operator on one sublattice | `asr_augcls` |
| `asr_ed` | binary energy-dependent route (AMscf form) | `asr_edc` |
| `asr_edn` | n-species energy-dependent route | `asr_edc` |
| `asr_endep` | binary potential-function routes; re-exports the chain machinery of `asr_chain` | `asr_chain` |
| `asr_lmto_realization` | ordered real-space operator | `asr_augcls` |
| `asr_recursion`, `asr_greens`, `asr_terminator`, `asr_kpm`, `asr_partial` | plain Lanczos, Green functions, terminators, kernel polynomials, partial disorder | inside `asr_edc`, `asr_ed3`, `asr_chain` |
| `asr_pfun`, `asr_kspace`, `asr_lattice_toy` | potential functions, toy k-space and lattice models | none (tests only) |

A reference module may use production modules. A production module never uses a reference module. The rule keeps
the program small. Run `tools/check_layers.py` to check it.
