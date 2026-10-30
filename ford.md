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
page_dir: ./doc/manual
output_dir: ./doc/build
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

<div id="cover-page" class="scientific-cover">
  <!-- Masthead -->
  <header class="cover-masthead">
    <div class="masthead-brand">
      <span class="masthead-dot"></span>
      <span class="masthead-lab">The MAVENs Group &bull; SRMIST</span>
    </div>
    <div class="masthead-meta">
      <span>Computational Condensed Matter</span>
      <span class="sep">&bull;</span>
      <span>v0.1.0 (2026)</span>
    </div>
  </header>

  <!-- Main Stage: Asymmetric Split Grid -->
  <main class="cover-stage">
    <!-- Left Column: Physics & Geometry -->
    <section class="stage-visual">
      <div class="visual-card">
        <div class="visual-frame">
          <img src="media/asr.svg" class="cover-logo" alt="ASR Logo">
        </div>
        <p class="theory-caption">Exact disorder averaging in Mookerjee augmented space via Haydock-Lanczos recursion</p>
        <div class="theory-tags">
          <span class="tag">Beyond CPA</span>
          <span class="tag">TB-LMTO-ASA</span>
          <span class="tag">High-Entropy Alloys</span>
          <span class="tag">Koelling–Harmon</span>
        </div>
      </div>
    </section>

    <!-- Right Column: Identity & Access Portals -->
    <section class="stage-content">
      <div class="identity-header">
        <div class="super-title">First-Principles Disordered Electronic Structure</div>
        <h1 class="main-title">asr</h1>
        <p class="main-subtitle">Augmented-Space Recursion Suite</p>
        <p class="main-abstract">
          A self-contained electronic structure code for random alloys and disordered materials.
          Integrates exact configuration averaging, scalar-relativistic radial solvers, screened LMTO structure constants,
          and self-consistent ASA electrostatics without empirical potential parameters.
        </p>
      </div>

      <!-- Curated Scientific Portals -->
      <nav class="portal-grid" aria-label="Documentation Portals">
        <a href="page/index.html" class="portal-card portal-card-primary">
          <div class="portal-icon"><i class="fa fa-book"></i></div>
          <div class="portal-text">
            <span class="portal-title">User Manual &amp; Theory</span>
            <span class="portal-desc">Installation, input syntax, workflows, and physics guides</span>
          </div>
          <span class="portal-arrow">&rarr;</span>
        </a>

        <a href="lists/modules.html" class="portal-card">
          <div class="portal-icon"><i class="fa fa-sitemap"></i></div>
          <div class="portal-text">
            <span class="portal-title">API &amp; Architecture Reference</span>
            <span class="portal-desc">Fortran modules, procedures, and derived type definitions</span>
          </div>
          <span class="portal-arrow">&rarr;</span>
        </a>

        <a href="https://github.com/mavensgroup/asr" class="portal-card">
          <div class="portal-icon"><i class="fa fa-github"></i></div>
          <div class="portal-text">
            <span class="portal-title">Source Code Repository</span>
            <span class="portal-desc">GitHub repository, validation suite, and issue tracker</span>
          </div>
          <span class="portal-arrow">&rarr;</span>
        </a>
      </nav>
    </section>
  </main>

  <!-- Colophon / Academic Provenance -->
  <footer class="cover-colophon">
    <div class="colophon-author">
      Developed by <a href="https://mavens-group.github.io/author/rudra-banerjee/"><strong>Rudra Banerjee</strong></a> &bull; <a href="https://mavens-group.github.io/">The MAVENs Group</a> &bull; Department of Physics &amp; Nanotechnology, SRMIST
    </div>
    <div class="colophon-specs">
      <span>Fortran 2008/2018</span> &bull; <span>OpenMP</span> &bull; <span>Libxc &ge; 7.0</span> &bull; <a href="https://opensource.org/licenses/MIT">MIT License</a>
    </div>
  </footer>
</div>
