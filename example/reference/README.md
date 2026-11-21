# Reference outputs of the examples

Each directory holds the output of one input in `example/`: the text report `asr.out`, the summary `asr.json` and the
density-of-states file.  A rerun with release v1.0.1 should reproduce them apart from time stamps, timings and the version line.
The runs predate the version number 1.0.1: release v1.0.1 is commit db797e3 plus this directory and the version number.
Host names and working directories have been removed.

Program: main at commit db797e3 (`git describe`: v1.0.0-7-gdb797e3), built with GCC 14.2.0, OpenBLAS and Open MPI 5.0.8.
Computer: two nodes of the PARAM Smriti cluster, one MPI rank of 48 OpenMP threads per node.

| Input | Route | levels | steps | spin | SCF iterations | Wall time (s) |
|---|---|---|---|---|---|---|
| `cuzn_b2.toml` | edc | 4 | 24 | 1 | 37 | 20 |
| `cuzn_alloy.toml` | edc | 4 | 24 | 1 | 34 | 942 |
| `feco_b2_asr.toml` | edc | 4 | 24 | 2 | 28 | 25 |
| `co2fe04cr06al_l21.toml` | edc | 4 | 24 | 2 | 47 | 927 |
| `fe66ni14cr20_fcc.toml` | edc | 3 | 24 | 1 | 27 | 780 |
| `fe66ni14cr20_fcc_sro.toml` | edc | 3 | 24 | 1 | 40 | 1155 |

On a workstation (Intel Core i7-14700, 28 OpenMP threads, Intel ifx, same commit), `cuzn_b2` takes 26 s, `feco_b2_asr`
41 s and `cuzn_alloy` 60 min, with the same iteration counts and Fermi energies as above.
