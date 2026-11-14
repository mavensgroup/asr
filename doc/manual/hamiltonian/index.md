---
title: Hamiltonians
ordered_subpage: tb_lmto_asa.md
                 adding.md
---

# Hamiltonians

**asr is a method, not a Hamiltonian.** The augmented-space recursion needs an
operator on a graph whose site-diagonal parts depend on the local species and
whose off-diagonal parts do not. It does not care where that operator comes
from. This chapter separates the two things.

| Layer | Question it answers | Where |
|---|---|---|
| Method | How is the disorder average taken, and how is the Green function obtained? | [Theory](../theory/index.html) |
| Hamiltonian | What are the site parameters, the hopping blocks and the energy dependence? | this chapter |

## What the method asks of a Hamiltonian

A Hamiltonian plugs in if it supplies the following.

1. **A graph.** Sites on a periodic lattice, and for every site the same list of
   neighbour slots. The hopping blocks belong to the slots and are the same for
   every site and every species (see below).
2. **A local basis.** A fixed number of orbitals per site, with a point-group
   representation if the symmetry reduction of the augmented space is wanted
   ([[asr_symmetry]], [[asr_augsym]]).
3. **Species-dependent site parameters.** Per species and spin, the site-diagonal
   quantities that become random operators in the augmented space.
4. **Species-independent hopping.** The off-diagonal part must not depend on the
   local composition. This is what makes the augmented-space operator exact and
   finite.
5. **An energy dependence.** Either none, so that a fixed operator can be
   applied (`apply` of [[asr_operator]]), or a rule that gives the operator at a
   complex energy \(z\) (one chain per contour node).
6. **A density rule.** How the Green function moments on the contour turn into
   charge and moment per species, and from there into the next potential
   (self-consistency).

Items 1 and 5 enter through the abstract operator `linop` of [[asr_operator]].
Items 2 to 4 enter through the site-class description of [[asr_augcls]]. Item 6
is the SCF layer ([[asr_class_scf]]).

## What is implemented

| Hamiltonian | Status | Basis | Page |
|---|---|---|---|
| **TB-LMTO-ASA** (second-order, screened \(\alpha\) representation, collinear spin) | implemented, the only one | spd, 9 orbitals per site | [TB-LMTO-ASA](tb_lmto_asa.html) |
| Wannier / tight-binding tables | not implemented | any | [Adding a Hamiltonian](adding.html) |
| Full-potential or non-collinear LMTO | not implemented | | |

!!! note "There is no `hamiltonian` key in the input file"
    Because only one Hamiltonian exists, the input file does not select it. The
    species block (`z`, `nval`, `konfig`, `occupations`), the `[scf]` keys `xc`,
    `relativistic` and `sim_beta` and the `[kspace]` keys `sigma`, `alpha_scale`,
    `std_alpha` and `second_order` are the **Hamiltonian-specific** part of the
    input. They are collected on [Input: Hamiltonian](../input/hamiltonian.html),
    so that the rest of the input (structure, disorder, recursion, output) can
    be read as Hamiltonian-independent. When a second Hamiltonian is added, a
    selector key will be introduced and the Hamiltonian page of the input
    reference will split by Hamiltonian.

## Hamiltonian-independent and Hamiltonian-specific code

| Hamiltonian-independent | Hamiltonian-specific (TB-LMTO-ASA) |
|---|---|
| [[asr_operator]], [[asr_recursion]], [[asr_greens]], [[asr_terminator]] | [[asr_radial]], [[asr_potpar]], [[asr_atom]], [[asr_xc]] |
| [[asr_contour]], [[asr_partial]], [[asr_mixing]] | [[asr_strconst]], [[asr_harmonics]], [[asr_ewald]] |
| [[asr_augspace]], [[asr_augn]], [[asr_symmetry]], [[asr_lattices]], [[asr_sublattice]] | [[asr_lmto_ham]], [[asr_lmto_cell]], [[asr_lmto_alloy]], [[asr_pfun]] |
| [[asr_augcls]] (operator build from class parameters) | [[asr_edc]], [[asr_ed]], [[asr_edn]], [[asr_ed3]] (the H(E) forms of the ASA) |
| [[asr_class_scf]] (SCF loop) | the potential, density and electrostatics steps called from it |

The boundary is not yet perfectly clean. For example [[asr_class_scf]] still calls
the ASA density and Madelung steps directly, and [[asr_augcls]] assumes the
\(h = c + dSd\) structure of the LMTO operator. The table shows the intended
split and the current state of the modules.
