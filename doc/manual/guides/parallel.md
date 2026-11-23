---
title: Parallelism and performance
---

# Parallelism and performance

`asr` uses a hybrid **MPI + OpenMP** model for `method = "asr"` ([[asr_mpi]], [[asr_chain]]). The same source builds either the serial/OpenMP executable or the MPI+OpenMP executable when compiled with `-DUSE_MPI`.

---

## 1. How work is distributed

A chain task is one `(start orbital, species, spin, contour energy)` tuple, giving

\[
N_{\text{tasks}} = 4\,N_{\text{species}}\,N_{\text{spin}}\,N_E
\]

independent Lanczos chains per SCF pass (for example, 256 chains for a non-magnetic binary alloy with 32 contour points).

1. **Across MPI ranks:** When built with `-DUSE_MPI`, chain tasks are split across MPI ranks with no change to the input file. Only rank `0` prints to stdout and writes output files.
2. **Within a rank (OpenMP):** While each OpenMP thread has multiple chains, each thread executes whole chains independently.
3. **Nested OpenMP teams:** When a rank has fewer chain tasks than OpenMP threads (e.g. across many nodes, or during `ef_linear_max` / single-energy passes), teams of threads share the matrix-vector products of individual chains (`chains: ... teams of ... threads` in the log). Results are independent of the split.

!!! warning "Do not launch a serial build under `mpirun`"
    If you compile without `-DUSE_MPI` and run with `mpirun`, every rank would run the full calculation and overwrite the same output files. `asr` detects this and stops immediately.

---

## 2. Recommended runtime environment

Set the math library to **single-threaded** mode so BLAS/LAPACK calls inside OpenMP chain loops do not oversubscribe cores:

```bash
export OMP_NUM_THREADS=16          # physical cores per rank (e.g. one rank per node)
export OMP_STACKSIZE=256M          # prevent stack overflow in deep recursions
export MKL_NUM_THREADS=1           # avoid oversubscription from oneMKL
export OPENBLAS_NUM_THREADS=1      # same if using OpenBLAS
ulimit -s unlimited
```

For Intel OpenMP when nested teams are active:

```bash
export KMP_HOT_TEAMS_MAX_LEVEL=2
```

If you pin threads to cores, use `spread` at the outer level and `close` at the inner level so each outer thread keeps its own core partition for its team:

```bash
export OMP_PLACES=cores
export OMP_PROC_BIND=spread,close
```

---

## 3. Building with MPI (`-DUSE_MPI`)

`asr_mpi` includes `mpif.h` directly and declares explicit interfaces, so it does not depend on compiler-specific `mpi.mod` files. Pass your MPI include and library flags alongside `ASR_FFLAGS` / `ASR_LDFLAGS`:

```bash
source env.sh
# OpenMPI example:
fpm build --compiler ifx \
  --flag "$ASR_FFLAGS -DUSE_MPI -I/path/to/openmpi/include" \
  --link-flag "$ASR_LDFLAGS -L/path/to/openmpi/lib -lmpi_mpifh -lmpi"

# MPICH / Intel MPI: replace -lmpi_mpifh -lmpi with -lmpifort -lmpi
```

Check load balance in `asr.out` (and `timing_s` in `asr.json`): the timing section reports `min`, `mean`, `max`, and `max/mean` of the chain wall time across ranks. A `max/mean` noticeably above `1.0` indicates rank imbalance.

---

## 4. Algorithmic speedups in `[asr]`

Two `[asr]` keys reduce runtime without changing the physics:

| Key | Default | Effect | Trade-off |
|---|---|---|---|
| `energy_block` | `0` | Packs `N` contour energies per chain block using `asr_augblk` (also used by `[epi]`). `energy_block = 8` measured **1.35× faster** on CoCrFeNi | ~6× more memory for `energy_block = 8`; rounding-level differences within the chain-length error |
| `ef_linear_max` | `0.0` | Takes \(E_F\) steps up to this size (Ry) to first order without re-running the contour (exact at self-consistency) | Avoids extra contour passes during the \(E_F\) root search |
| `depth_tol` | `1e-5` | Stops individual chains early (checked from level 8 every 4 levels) once \(G\) converges within `depth_tol` | `levels` becomes an upper bound; check `ed_nlev` in `asr.out` |

Note that `method = "kspace"` has no MPI split: if run under MPI, rank `0` executes the k-space calculation while other ranks wait.
