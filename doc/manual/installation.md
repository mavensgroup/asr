---
title: Installation
---

# Installation

This guide covers everything needed to install dependencies, obtain the `asr` source code, compile the executable and test suite, and configure your runtime environment on a Linux workstation or high-performance computing (HPC) cluster.

---

## 1. System Requirements & Dependencies

`asr` is written in standard-conforming modern Fortran (Fortran 2008/2018) and uses the [Fortran Package Manager (`fpm`)](https://fpm.fortran-lang.org) for dependency resolution, compilation, and testing.

### Required Dependencies

| Component | Minimum / Target Version | Role in `asr` | How it is provided |
|---|---|---|---|
| **Fortran compiler** | Intel `ifx` ≥ 2024.0 *(primary)* | Compiles Fortran 2008/2018 source and OpenMP regions | Intel oneAPI HPC Toolkit |
| **BLAS & LAPACK** | Intel oneMKL *(or OpenBLAS / Netlib)* | Dense linear algebra (`dgesv`, `zheev`, `dsterf`) and threaded kernels | Bundled with Intel oneAPI (`-qmkl=parallel`) |
| **OpenMP** | OpenMP 4.5+ | Parallel complex-energy contour and Brillouin-zone loops | Bundled with compiler (`-qopenmp`) |
| **[fpm](https://fpm.fortran-lang.org)** | ≥ 0.10 | Build system, dependency fetcher, and test runner | Standalone binary, Conda, Pip, or Spack |
| **[libxc](https://libxc.gitlab.io)** | **≥ 7.0** (C library) | LDA (`vbh`, `vwn`, `pw`, `pwmod`) and GGA (`pbe`) exchange–correlation | System package (`dnf`/`apt`), Conda, Spack, or built from source |
| **[toml-f](https://github.com/toml-f/toml-f)** | `v0.5.2` | Validated TOML input parser | **Fetched automatically** by `fpm` on first build (`build/dependencies/toml-f`) |
| **Git** | ≥ 2.0 | Cloning `asr` and fetching `toml-f` | System package manager |

### Optional Dependencies

| Component | Version | Purpose |
|---|---|---|
| **Python 3 + NumPy** | Python ≥ 3.9 | Running the Special Quasirandom Structure generator (`tools/sqs.py`) |
| **[FORD](https://github.com/Fortran-FOSS-Programmers/ford)** | ≥ 7.0 | Building this HTML user manual and API reference (`ford ford.md`) |
| **Graphviz (`dot`)** | any recent | Generating module dependency and call graphs in FORD documentation |

!!! note "Why no Fortran `.mod` files are needed for `libxc`"
    `asr` calls `libxc` exclusively through its stable **C ABI** using Fortran's intrinsic `iso_c_binding` (see [[asr_xc]]). You do **not** need `libxc` Fortran interface modules (`xc_f03_lib_m.mod`), and you do **not** need to recompile `libxc` with `ifx`. Any system or CMake build of `libxc` providing `libxc.so` (or `libxc.a`) is sufficient.

---

## 2. Installing the Dependencies

Choose the installation method below that matches your system (Linux distribution package manager, Conda/Mamba, Spack, or building from source).

### 2.1 Intel oneAPI Compiler (`ifx`) and oneMKL

All development, validation benchmarks, and regression tests target Intel's LLVM-based Fortran compiler (`ifx`) paired with Intel oneMKL.

* **Workstations (Linux installer / repositories)**:
  Install the free **Intel oneAPI Base Toolkit** (for oneMKL) and **Intel oneAPI HPC Toolkit** (for `ifx`) from the [Intel oneAPI downloads page](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html) or via your distribution's package manager:

```bash
# Fedora / RHEL / Rocky Linux (after adding the Intel oneAPI yum/dnf repo)
sudo dnf install intel-oneapi-compiler-fortran intel-oneapi-mkl-devel

# Ubuntu / Debian (after adding the Intel oneAPI apt repo)
sudo apt install intel-oneapi-compiler-fortran intel-oneapi-mkl-devel
```

* **HPC Clusters (Environment Modules / Lmod)**:
  Load the compiler and MKL modules provided by your cluster administrator:

```bash
module load intel-oneapi-compilers intel-oneapi-mkl
```

Verify that `ifx` is available after loading the environment:

```bash
source /opt/intel/oneapi/setvars.sh
ifx --version
```

### 2.2 Exchange–Correlation Library (`libxc` ≥ 7)

`fpm.toml` links against `-lxc` (`link = ["xc"]`). Ensure the `libxc` shared or static library (`libxc.so` / `libxc.a`) is installed and visible to the linker.

#### Option A: System Package Manager (Fedora 40+ / Rolling Releases)

Recent Fedora and rolling-release distributions ship `libxc` ≥ 7 directly:

```bash
# Fedora
sudo dnf install libxc libxc-devel

# Verify installed version (should be >= 7.0)
pkg-config --modversion libxc
```

#### Option B: Conda / Mamba (Cross-platform, No Root Required)

```bash
conda install -c conda-forge "libxc>=7.0"
export LIBRARY_PATH="$CONDA_PREFIX/lib:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH}"
```

#### Option C: Building `libxc` ≥ 7 from Source (Ubuntu/Debian LTS, RHEL, or HPC Clusters)

If your Linux distribution ships an older `libxc` (5.x or 6.x) or you do not have root privileges, build `libxc` 7.x into `~/.local` using CMake (takes under a minute):

```bash
curl -L https://gitlab.com/libxc/libxc/-/archive/7.0.0/libxc-7.0.0.tar.bz2 | tar xj
cmake -B libxc-7.0.0/build -S libxc-7.0.0 \
      -DCMAKE_INSTALL_PREFIX="$HOME/.local" \
      -DBUILD_SHARED_LIBS=ON \
      -DENABLE_FORTRAN=OFF \
      -DBUILD_TESTING=OFF
cmake --build libxc-7.0.0/build -j
cmake --install libxc-7.0.0/build
```

Then expose `~/.local/lib` (or `~/.local/lib64`) in your shell profile (`~/.bashrc`):

```bash
export LIBRARY_PATH="$HOME/.local/lib:$HOME/.local/lib64:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="$HOME/.local/lib:$HOME/.local/lib64:${LD_LIBRARY_PATH}"
```

### 2.3 Fortran Package Manager (`fpm`)

Install `fpm` using any of the following methods:

```bash
# Method 1: Single static Linux x86_64 binary (no dependencies required)
mkdir -p "$HOME/.local/bin"
curl -Lo "$HOME/.local/bin/fpm" \
  https://github.com/fortran-lang/fpm/releases/download/v0.12.0/fpm-0.12.0-linux-x86_64-gcc-12
chmod +x "$HOME/.local/bin/fpm"

# Method 2: Conda / Mamba
conda install -c conda-forge fpm

# Method 3: Python pip / pipx
pip install --user fpm

# Method 4: Spack (HPC clusters)
spack install fpm
```

Verify with `fpm --version`. Ensure `$HOME/.local/bin` is on your `PATH`.

---

## 3. Obtaining the Source Code

The official `asr` repository is hosted on GitHub at **<https://github.com/mavensgroup/asr>**.

Clone the repository with Git:

```bash
# Clone via HTTPS
git clone https://github.com/mavensgroup/asr.git
cd asr

# Or clone via SSH (for contributors with SSH keys configured)
git clone git@github.com:mavensgroup/asr.git
cd asr
```

### Repository Layout

| Directory / File | Contents |
|---|---|
| `fpm.toml` | Package manifest (`asr` metadata, `libxc` link target, `toml-f` dependency) |
| `env.sh` | Environment setup script (loads Intel oneAPI, sets compiler/linker flags, defines helper functions) |
| `src/` | Library modules (`asr_radial`, `asr_xc`, `asr_potpar`, `asr_strconst`, `asr_asa`, `asr_augspace`, `asr_ed`, `asr_endep`, …) |
| `app/main.f90` | Command-line executable entry point (`asr <input.toml>`) |
| `test/` | Automated verification and regression test programs (`test_*.f90`) |
| `example/` | Ready-to-run TOML input files covering all four calculation kinds (`model`, `lmto`, `alloy`, `sublattice`) |
| `data/reference/` | Independent benchmark reference data (Stuttgart LMTO, NIST, exact configurations) |
| `tools/sqs.py` | Special Quasirandom Structure (SQS) generator script |
| `doc/` & `ford.md` | FORD user manual (`doc/manual/`), custom stylesheet (`doc/css/asr.css`), and generated HTML (`doc/build/`) |

---

## 4. Building `asr`

### 4.1 Recommended Build Using `env.sh`

From the root of the cloned `asr` repository, source `env.sh` and run `asr_build`:

```bash
source env.sh          # loads Intel oneAPI (ifx, MKL) and defines fpm wrappers
asr_build --tests      # compiles libasr, the asr executable, and all test programs
```

On the first run, `fpm` automatically clones `toml-f` (`v0.5.2`) into `build/dependencies/toml-f` and compiles it alongside `asr`.

#### What `env.sh` configures

`env.sh` sources `/opt/intel/oneapi/setvars.sh` (edit line 2 of `env.sh` if your oneAPI installation lives elsewhere or is loaded via `module load`) and defines three shell functions (`asr_build`, `asr_test`, `asr_run`) with the following compiler and linker flags:

```text
FPM_FC      = ifx
ASR_FFLAGS  = -O2 -g -traceback -qopenmp -qmkl=parallel -free -standard-semantics
              -warn all -fp-model=precise -heap-arrays 1024
ASR_LDFLAGS = -qopenmp -qmkl=parallel
```

| Flag | Why it is used |
|---|---|
| `-O2 -g -traceback` | Full optimization while preserving symbolic line-number tracebacks on runtime errors |
| `-qopenmp` | Enables OpenMP multi-threading across complex contour points and k-points |
| `-qmkl=parallel` | Links threaded Intel oneMKL for BLAS/LAPACK (`zheev`, `dgesv`, `dsterf`) |
| `-free -standard-semantics` | Free-form Fortran with strict F2018 standard semantics (including value-initialized module variables and reallocation on assignment) |
| `-fp-model=precise` | Prevents unsafe floating-point associative rearrangements that would break reproducibility of Lanczos recursion and SCF mixing |
| `-heap-arrays 1024` | Allocates temporary arrays larger than 1 MB on the heap rather than the thread stack, preventing stack overflows in deep augmented-space recursions |

### 4.2 Installing the `asr` Binary onto Your `PATH`

To run `asr <input.toml>` directly from any working directory without invoking `fpm`, install the compiled executable into `~/.local/bin` (or any prefix on your `PATH`):

```bash
source env.sh
fpm install --compiler ifx \
            --flag "$ASR_FFLAGS" \
            --link-flag "$ASR_LDFLAGS" \
            --prefix "$HOME/.local"
```

This places the binary at `~/.local/bin/asr`. Verify with:

```bash
which asr
```

### 4.3 Building on Air-Gapped HPC Nodes (Offline Build)

Compute or compile nodes without external internet access cannot fetch `toml-f` from GitHub during `fpm build`. On a login node with internet access (or on your workstation before syncing), pre-fetch the dependency once:

```bash
# On a machine/node with internet access:
git clone --depth 1 --branch v0.5.2 https://github.com/toml-f/toml-f.git build/dependencies/toml-f
```

Once `build/dependencies/toml-f` exists, `asr_build` runs completely offline.

!!! warning "Stale build objects after a module signature change"
    `fpm`'s incremental dependency scanner tracks module `.mod` files, but can occasionally retain stale `.o` objects if a procedure interface or derived-type layout changes during development. If you encounter unexplained link or runtime errors after modifying code in `src/`, wipe the build cache and recompile:

        rm -rf build/ifx_*
        asr_build --tests

---

## 5. Verifying the Installation

Always run the automated test suite after building on a new machine to verify the radial solver, `libxc` linkage, structure constants, k-space tetrahedron integration, and augmented-space recursion:

```bash
export OMP_NUM_THREADS=16
asr_test
```

The full test suite runs in **≈ 3 minutes on 16 cores**. You can also run individual test executables by name for a fast check (a few seconds):

```bash
asr_test test_radial       # hydrogenic levels to 1e-10 Ry; He/Ne vs NIST LDA to 3e-7 Ha
asr_test test_augspace     # augmented-space recursion vs exact 2^N configuration average (1e-14)
asr_test test_asa          # self-consistent fcc Ag, bcc Fe, fcc Cu/Ni, and B2 FeCu
```

Every test program prints `passed N, failed 0.` and exits with status `0`. Next, run a 1-second end-to-end smoke test on an example input:

```bash
mkdir -p /tmp/asr_smoke && cp example/fcc_ordered.toml /tmp/asr_smoke/
asr_run -- /tmp/asr_smoke/fcc_ordered.toml
```

Expected output:

```text
asr 0.1.0  fcc NN s-band
input: /tmp/asr_smoke/fcc_ordered.toml
wrote /tmp/asr_smoke/fcc_ordered_dos.dat
```

---

## 6. Runtime Environment & Command-Line Usage

### Command-Line Syntax and Exit Codes

```bash
asr <path/to/input.toml>
# or via fpm:
asr_run -- <path/to/input.toml>
```

`asr` accepts **one** positional argument: the path to a TOML input file. All output files (`asr.out`, `*_dos.dat`, `*_params.dat`, saved potentials) are written **in the same directory as the input file**.

| Exit Code | Meaning |
|---|---|
| `0` | Calculation finished normally and (for self-consistent kinds) SCF converged within `tolerance` |
| `1` | Input validation error — diagnostic printed to stdout as `asr: input error: <reason>` |
| `2` | Invalid command-line invocation (not called with exactly one argument) |
| `3` | SCF loop reached `max_iterations` without meeting `tolerance` |

### OpenMP Threading and Stack Size

Energy contour integrations (`[asr] contour_points`), k-space meshes, and energy-dependent Lanczos recursions are parallelized with OpenMP. Before running calculations, set the thread count and give each OpenMP worker thread adequate stack space:

```bash
export OMP_NUM_THREADS=16      # set to physical core count
export OMP_STACKSIZE=256M      # recommended for deep augmented-space recursions
ulimit -s unlimited            # lift primary thread stack limit
```

The active thread count is recorded at the top of `asr.out`.

---

## 7. Building This Documentation

To build the HTML user manual and API reference locally from the source tree:

```bash
pip install --user ford
# Optional: install Graphviz ('dot') if you want call/module dependency graphs
# sudo dnf install graphviz   (Fedora)
# sudo apt install graphviz   (Ubuntu/Debian)

ford ford.md
xdg-open doc/build/index.html
```

Set `graph: true` in `ford.md` if Graphviz is installed and you wish to generate interactive dependency diagrams across modules and procedures.

---

## 8. Troubleshooting Installation Issues

| Symptom | Cause | Fix |
|---|---|---|
| `ifx: command not found` | Intel oneAPI environment not loaded | Run `source /opt/intel/oneapi/setvars.sh` (or `source env.sh`), or load your cluster's oneAPI module |
| `ld: cannot find -lxc` | Linker cannot find `libxc.so` or `libxc.a` | Install `libxc-devel` (`dnf`) / `libxc-dev` (`apt`), or append the `libxc` library directory to `LIBRARY_PATH` |
| `libxc.so.15: cannot open shared object file` | Runtime loader cannot find `libxc.so` | Append the `libxc` library directory (e.g. `~/.local/lib64`) to `LD_LIBRARY_PATH` |
| `forrtl: severe (174): SIGSEGV` on large runs | Thread stack exhausted in Lanczos recursion | Run `ulimit -s unlimited`, set `export OMP_STACKSIZE=256M`, and keep `-heap-arrays 1024` |
| Link/runtime errors after pulling commits | Stale `.mod` or `.o` files in `build/` | Clean and rebuild: `rm -rf build/ifx_* && asr_build --tests` |

Next step: head to the **[Quick start](quickstart.html)** to run four progressively richer calculations (`model`, `lmto`, `alloy`, and `sublattice`).
