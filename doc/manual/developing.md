---
title: Developing
---

# Developing

## Workflow

Clone the repository from GitHub:

```bash
git clone https://github.com/mavensgroup/asr.git
cd asr
```

Then load the environment and build:

```bash
source env.sh
asr_build --tests         # build library, app, tests
asr_test test_augspace    # one test (fpm test NAME)
asr_test                  # everything, ≈ 3 min on 16 cores
```

* **After a module *signature* change**: `rm -rf build/ifx_*` (fpm's incremental
  build can keep stale objects).
* **Scratch work** goes in `trial/` (git-ignored; keep inputs, never commit
  outputs). Anything worth keeping moves to `example/` (inputs) or
  `data/reference/` (truth data).
* **Code repository and sync**: The public upstream repository is
  [https://github.com/mavensgroup/asr](https://github.com/mavensgroup/asr).
  `env.sh` provides `asr_sync` (push `main` to GitHub and fast-forward the desktop
  working copy) and `asr_pull`.
* **Heavy runs** (deep recursions) go to the desktop; one process per run
  directory.

## Writing a test

Tests are programs `test/test_<topic>.f90` picked up automatically
(`auto-tests = true`). Use the in-repo harness [[asr_testing]]:

```fortran
!> One line stating what is verified and against what independent truth.
program test_example
  use asr_kinds, only: dp
  use asr_testing
  implicit none
  print '(a)', 'test_example'
  call check(1 + 1 == 2, 'arithmetic still works')
  call check_close(0.1_dp + 0.2_dp, 0.3_dp, 1.0e-12_dp, 'sum')
  call report()           ! prints "passed n, failed m." and stops with 1 on failure
end program test_example
```

The convention of the project: **test against an independent truth** — dense
diagonalisation, exact enumeration of all configurations, k-space sums, NIST
atoms, a different recursion route. A test that only checks that a function
reproduces itself is not accepted as validation. Exploratory scripts are named
`probe_*` and do not assert.

## Adding a new Hamiltonian

The recursion needs only `apply`. Extend the abstract type of [[asr_operator]]:

```fortran
type, extends(linop) :: my_op
contains
  procedure :: apply => my_apply   ! y = H x
  procedure :: dim   => my_dim
end type my_op
```

then call `lanczos(op, v, nlevels, a, b, m)` ([[asr_recursion]]) and
`cf_green(a, b, n, z, tail)` ([[asr_greens]]). For an augmented-space operator
also override `mapply` (the configuration projector \(\mathcal M\) used by the
species-resolved Green function).

## Adding an input key

1. Add the key to the `allowed` list passed to `check_keys` for its section — the
   reader rejects everything not listed.
2. `get_value(child, 'key', variable, default)`.
3. Validate the value at the end of the reader with a clear message.
4. Mirror it in the module header comment, add a row to the
   [input reference](input/index.html), and add a case to the matching
   `test_*_input`.

## Documentation conventions (FORD)

The API reference is built by [FORD](https://github.com/Fortran-FOSS-Programmers/ford)
from the source comments; this manual is the set of pages in `doc/manual/`.

| Marker | Where | Use |
|---|---|---|
| `!>` | line(s) **before** a module, type, procedure | Description (Markdown, MathJax) |
| `!!` | **after** an argument or component declaration | Description of that item |

```fortran
!> Linear interpolation of the table at x.
!> Extrapolation is an error.
function interp(x) result(y)
  real(dp), intent(in) :: x   !! abscissa, Rydberg
  real(dp) :: y               !! interpolated value
```

* Math in comments and pages: inline `\( … \)`, display `\[ … \]` or `$$ … $$`.
* Link to API entities from pages and comments with `[[name]]`, e.g. `[[asr_augspace]]`,
  `[[asr_augspace(module):build_augspace]]` for a specific procedure.
* The first `!>` block of a module is its summary on the module list: make
  the first sentence stand alone.
* Admonitions: `!!! note`, `!!! warning`, `!!! tip`, `!!! danger` with an
  indented body.

### Manual pages

Each page in `doc/manual/` begins with front matter:

```markdown
---
title: Page title
ordered_subpage: child1.md
                 child2.md
---
```

`ordered_subpage` fixes the order of children; a directory is a page via its
`index.md`. Cross-page links are ordinary relative Markdown links with the
`.html` extension (`[text](input/alloy.html)`).

### Building

```bash
ford ford.md              # doc/build/index.html
```

Check the output for `Warning:` lines about unresolved `[[…]]` links and fix them.

## Code style

* Free-form, `implicit none`, `private` by default with explicit `public` lists.
* `real(dp)` from [[asr_kinds]]; Rydberg atomic units internally, Bohr lengths.
* Input errors are returned as a message string (`msg`), never `stop` in a reader.
* Keep the stale-object rule in mind when changing `type` definitions.
* Do not paste numbers from other codes into the source: parameters are generated.
