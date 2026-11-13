# Numerical checks of docs/DERIVATION_f4.md

Small real pencils $(A, M)$ on simple cubic clusters (`f4lab.py`), `numpy` and `scipy` only.  Run each script in this directory.

| script | checks | time |
|---|---|---|
| `f4_check1.py` | Result 2 (error identity and bound); first scan for $-\mathrm{Im}\,G \le 0$ | 20 s |
| `f4_check2.py` | F4 with and without full reorthogonalisation (exact arithmetic) | 2 min |
| `f4_check3.py` | Result 7, sufficiency; counts for the `hull` tail and for real edges | 45 s |
| `f4_check4.py` | Result 7, necessity; Result 5 (minimal example, writes `minimal_example.npy`); Result 10 | 2 s |
| `f4_check5.py` | Result 6 (word-space pencil); Pick matrix and positive-measure projection | 5 s |
| `f4_check6.py` | Result 4 (slope of the Ritz values); pointwise positive alternatives | 4 s |
| `mk_synth_dump.py` | synthetic `ef_scan_dump` file of the minimal example for `tools/f4_dump.py` (run `f4_check4.py` first) | 1 s |
