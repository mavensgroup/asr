---
project: asr
author: Rudra Banerjee
author_description: Assistant Professor, Computational Condensed Matter, SRM Institute of Science and Technology (SRMIST). Leads the <a href="https://mavens-group.github.io/">MAVENs group</a>: disordered and high-entropy alloys, 2D and quantum materials, magnetocaloric and spintronic compounds.<br><a href="https://orcid.org/0000-0002-5028-4844">ORCID</a> · <a href="https://scholar.google.com/citations?user=FU6cIKAAAAAJ">Scholar</a>
website: https://mavens-group.github.io/author/rudra-banerjee/
github: https://github.com/rudrab
linkedin: https://www.linkedin.com/in/rudra-banerjee-a6bb99b7
license: mit
version: 0.1.0
src_dir: ./src
         ./app
output_dir: ./doc/ford_api
exclude_dir: ./build
             ./trial
             ./test
media_dir: ./doc/media
favicon: ./doc/media/favicon.png
extensions: f90
display: public
         protected
source: true
graph: false
search: true
sort: alpha
proc_internals: false
md_extensions: markdown.extensions.toc
               markdown.extensions.sane_lists
               markdown.extensions.admonition
css: ./doc/css/asr.css
creation_date: %Y-%m-%d
print_creation_date: true
parallel: 4
quiet: false
warn: false
---

# Fortran 2008/2018 API Reference

Welcome to the internal source code and API reference for **asr** (Augmented-Space Recursion Suite). This documentation is extracted directly from the Fortran docstrings (`!>`, `!!`) in `src/` and `app/` using FORD.

## API Navigation

| Section | Description | Direct Link |
| :--- | :--- | :--- |
| **Modules** | All 52 Fortran modules implementing recursion, TB-LMTO-ASA, and SCF | [Browse Modules](lists/modules.html) |
| **Derived Types** | Public & internal data types (`potpar_t`, `augorb_op`, `cell_t`, etc.) | [Browse Types](lists/types.html) |
| **Procedures** | Mathematical solvers, Hamiltonians, Lanczos recursion, and utilities | [Browse Procedures](lists/procedures.html) |
| **Source Files** | Syntax-highlighted Fortran source code | [Browse Source](lists/files.html) |
| **Programs** | Executable drivers (`asr` CLI entrypoint) | [Browse Programs](lists/programs.html) |

---
&larr; [Return to asr User Manual](../manual.html)

