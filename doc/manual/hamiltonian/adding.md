---
title: Adding a Hamiltonian
---

# Adding a Hamiltonian

This page lists what has to be written to run the augmented-space recursion on
a Hamiltonian other than TB-LMTO-ASA. It is a guide for developers, not a
promise of a feature.

1. **Parameters.** A routine that gives, for each species and spin, the
   site-diagonal parameters of the new model (the analogue of
   \(E_\nu, C, \sqrt\Delta, o\) of [TB-LMTO-ASA](tb_lmto_asa.html)).
2. **Hopping blocks.** One block per neighbour slot of the Bravais graph,
   common to all species. For a model with species-dependent hopping, the
   method does not apply as it stands.
3. **Operator.** Either extend `linop` ([[asr_operator]]) for a fixed operator,
   or provide the \(z\)-dependent form as [[asr_edc]] does for the ASA.
4. **Class description.** Fill `class_alloy_t` ([[asr_augcls]]) with the
   parameters and the species sets of the site classes.
5. **Density rule.** Give the map from contour moments to charge and moment
   per species (the analogue of the \(\phi^2\), \(\phi\dot\phi\) expansion
   in [Energy-dependent ASR](../theory/energy_dependent.html)).
6. **Input.** Add the species keys to the strict reader in [[asr_driver]] and a
   selector key under `[calculation]`.
7. **Tests.** The two checks every Hamiltonian must pass are the exact
   \(2^N\) configuration average and agreement of the A = B alloy with the
   ordered result (see [Validation](../validation.html)).
