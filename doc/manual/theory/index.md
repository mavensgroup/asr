---
title: Theory
ordered_subpage: augmented_space.md
                 recursion.md
                 lmto.md
                 energy_dependent.md
                 scf_electrostatics.md
---

# Theory

This chapter explains what the code computes and why the numerical methods
work, at the depth needed to interpret results and read the source. It is not a
substitute for the original papers (see [References](../appendix/references.html)).

| Page | Question it answers |
|---|---|
| [Augmented space](augmented_space.html) | How can disorder be averaged *exactly* with a single Hamiltonian? |
| [Recursion and terminators](recursion.html) | How is the Green function extracted, and why does depth matter? |
| [TB-LMTO-ASA](lmto.html) | What is the Hamiltonian, and which parameters does the code generate? |
| [Energy-dependent ASR](energy_dependent.html) | What are the three recursion routes (`energy-independent`, `energy-dependent`, `asr_ed`)? |
| [SCF and electrostatics](scf_electrostatics.html) | How are charges, Fermi level, Madelung energy and mixing handled? |

## The pipeline

```text
 input.toml ─▶ species ─▶ radial solver ─▶ potential parameters (E_ν, C, √Δ, p, γ)
                                              │
        structure constants S^α  ◀────────────┘
                 │
                 ▼
   operator H (ordered: H(k) | disordered: augmented-space operator on a graph)
                 │
        Lanczos / Haydock recursion  →  (a_n, b_n)
                 │
     terminated continued fraction  →  G(z)  →  DOS, moments M_k over a contour
                 │
   charge density → Hartree + xc + Madelung → mixing ──┐
                 ▲                                       │
                 └───────────── self-consistency ◀───────┘
```
