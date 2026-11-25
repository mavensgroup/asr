!> Self-consistent alloy with sublattices (site classes) from the energy-dependent augmented-space recursion (asr_ed form,
!> asr_augcls + asr_edc): ordered compounds, partial (sublattice) disorder and multicomponent alloys.  Collinear spin.
!>
!> A component p = (sublattice kappa, chemical species s) has its own potential, E_nu and potential parameters; its weight per
!> site is w_p = f_kappa x_{s,kappa} (f_kappa: fraction of the sites on sublattice kappa).  Equal ASA spheres (one common S).
!> Per iteration:  component potentials -> potential parameters -> for every sublattice one class operator with the origin on
!> that sublattice; for every component p of it the start site is frozen to p and the operator re-parametrised with p's Delta
!> (asr_edc_setup) -> G_p(z) on the Fermi contour, one start orbital per orbital group of the site (cubic: s, p_z, t2g, e_g)
!> -> moments M_k (common E_F from sum_p w_p N_p = sum_p w_p Z_val,p) -> density M0 phi^2 + 2 M1 phi phidot + M2 (phidot^2 +
!> phi phiddot), relaxed core, Hartree, xc ->
!>     V_p = sum_kappa' Phi(kappa, kappa') qbar_kappa'  +  sum_p' F(p, p') (q_p' - qbar_kappa'),   qbar_kappa = sum_{p in kappa} x_p q_p:
!> Madelung of the sublattice-average charges over the decoration cell (Ewald) plus the screened potential of the deviations
!> from the sublattice averages, F from asr_screening (charge-excess model of the graph lattice, hardness solved from sim_beta).
!> Without SRO F = -2 sim_beta / R1 (the screened-impurity model; ordered sublattice: no deviation).  With first-shell SRO,
!> Ahat(p, p') = z(kappa, kappa') (X_s'^(p) - x_s') (the conditional shell; derivation in asr_screening), and
!> Gamma = dq/db is the closure of asr_screening on every eigenmode of Ahat.  The SRO potential is
!> then iterated in the drive u = db of the charge-excess model (asr_screening::screen_update, u is part of the mixing vector):
!>     u_out = dq + (I - Gamma) u_in,   V = a (I - Gamma) u_out.
!> The fixed point is V = F dq, F = a (Gamma^-1 - I), and dq = 0 on a mode with Gamma = 0 (full clustering, F infinite).
!> u starts at zero unless start_potential holds it: a potential without an SRO part has no drive.
!> Without SRO the potential is F dq as before (bit for bit).
!> -> Anderson mixing of the smooth potentials and E_nu.
!> Total energy per site (Ry, Kohn-Sham form: output density, input potential), E = sum_p w_p E_p with
!>     E_p = E_band + sum eps_core - int n v~_in + 1/2 int n V_H[n] + int n eps_xc[n] + 1/2 q_p vm_p,
!> E_band = sum_l (M1_l + E_nu,l M0_l) (M0 = sphere charge of the channel), v~ = v + 2Z/r (the nuclear terms -2Z int n/r of T and of the sphere electrostatics cancel),
!> vm_p = Madelung + screened fluctuation shift from the output charges (1/2 q vm = Madelung energy of the sublattice charges +
!> 1/2 sum w dq F dq; diag(w) F is symmetric, so this is variational).  The Fermi-smearing entropy (k T ~ 1 mRy) is not included.
!> The chains (spin x component x start orbital x energy) are independent (H depends on E) and run in parallel (OpenMP), one
!> per task; every thread keeps its own copy of the class operator.
!>
!> route = 'ed3' (asr_ed3; manual/src/theory/energy_dependent.md, section ed3): the full generalised ASA problem instead of H^(2).  The chains give the central
!> blocks xm0 = G_11, xm1 = (hG)_11, xm2 = (hGh)_11 and the moments are those of k-space (bands_moments): q0 = xm0 + 2 o xm1 + o^2 xm2,
!> q1 = xm1 + o xm2, q2 = xm2, density q0 phi^2 + 2 q1 phi phidot + q2 phidot^2, channel charge q0 + p^gamma q2 = xm0 + 2 o xm1 + p xm2.
!> E_band = sum_l (q1 + E_nu,l (q0 + p^gamma q2)) as in asa_scf: the occupied eigenvalue sum.  E_F comes from the channel charges.  gz(...) then holds G + 2 o hG + p hGh, whose
!> -Im / pi is the DOS of the sphere charge, so dos_final and ef_scan need no other change.
module asr_class_scf
  use asr_kinds, only: dp, pi
  use asr_lapack, only: zheev
  use asr_radial
  use asr_xc
  use asr_atom, only: cumint, poisson_hartree
  use asr_potpar, only: potpar_l_t, potpar_at_enu, clamp_enu
  use asr_lmto_ham, only: potpar_t, to_alpha, o_alpha
  use asr_asa, only: asa_class_t, asa_input_t, asa_class_result_t, atomic_start
  use asr_strconst, only: standard_alpha
  use asr_kinds, only: i8
  use asr_runstats
  use asr_augspace, only: graph_t, graph_distances
  use asr_augcls
  use asr_edc
  use asr_ed3
  use asr_augblk, only: blk_work_t, blk_work_init, blk_nbmax
  use asr_chain, only: ed_depth_tol, ed_nchain, ed_nlev, ed_nfall, ed_tail, chain_teams, chain_vec_par
  !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads, omp_get_max_active_levels, omp_set_max_active_levels
  use asr_crystal, only: structure_t, ws_radius
  use asr_crystal_symmetry, only: symop_t, space_group, is_monomial, orbital_matrix, orbital_groups, reduction_ops
  use asr_crystal_graph, only: crystal_graph_t, crystal_graph, graph_site_ops
  use asr_lmto_alloy, only: structure_blocks
  use asr_sublattice
  use asr_contour, only: contour_nodes_closed, contour_nodes_fermi, fermi_npole, fermi_nodes_count, gauss_legendre
  use asr_mixing, only: anderson_t, anderson_init, anderson_step
  use asr_screening, only: screen_t, screen_structure, screen_matrix, screen_update, screen_project, mode_f, screen_structure_shells
  use asr_epi
  use asr_mpi, only: par_rank, par_size, par_root, par_sum_z, par_sum_i8, par_sum_i
  use asr_spectral, only: spectral_input_t, spectral_kpath, spectral_structure_ok, spectral_origin_classes, spectral_unfold, &
                          spectral_dressing, spectral_window
  implicit none
  private
  public :: class_scf_input_t, class_scf_result_t, class_scf

  type :: class_scf_input_t
    type(asa_class_t), allocatable :: spec(:)          !! chemical species
    type(structure_t) :: st                            !! crystal structure (primitive vectors, basis sites, sublattice labels)
    type(decoration_t) :: dec                          !! optional cubic decoration (used if st%nb == 0)
    integer, allocatable :: nsc(:), spc(:, :)           !! per sublattice: number of species, chemical species spc(k, sublattice)
    real(dp), allocatable :: conc(:, :)                 !! conc(k, sublattice)
    logical :: sro = .false.                            !! first-shell short-range order (asr_sro)
    real(dp), allocatable :: sro_cond(:, :, :, :)       !! sro_cond(k, k', kappa, kappa') = P(species k' on a first-shell site of
                                                        !! sublattice kappa' | species k at a site of kappa)
    real(dp) :: alat = 0.0_dp                          !! graph lattice constant (Bohr)
    character(len=8) :: lattice = 'bcc', xc = 'vbh'
    real(dp) :: rcl = 1.2_dp, rinv = 0.0_dp            !! kept S^alpha range, inversion cluster radius (alat)
    integer :: nlevels = 4, nrec = 0, maxrank = 1, flip_radius = 99, ncont = 32
    real(dp) :: beta = 0.3_dp, tol = 1.0e-5_dp, sim_beta = 1.0_dp, zeta0 = 0.1_dp
    real(dp) :: beta_restart = 0.0_dp                  !! linear step after an Anderson restart (0: beta)
    integer :: maxit = 60, mix_history = 6, nspin = 1
    logical :: relativistic = .true.
    real(dp) :: ed_eta = 0.01_dp, fermi_kt = 0.001_dp, depth_tol = 1.0e-5_dp
    integer :: fermi_poles = 0
    character(len=8) :: ed_tail = 'hull'                !! asr_chain%ed_tail; 'bp' is not safe on complex chains
    character(len=8) :: route = 'edc'                   !! 'edc' (H^(2), asr_edc) or 'ed3' (full generalised ASA problem, asr_ed3)
    real(dp) :: alpha_scale = 1.0_dp                    !! diagnostic: screening constants = alpha_scale * standard_alpha()
    integer :: eblock = 0                               !! contour energies per chain block (asr_augblk); 0 (default): one energy per chain, real kernels
    real(dp) :: ef_linear_max = 0.0_dp                  !! E_F steps up to this size (Ry) are first order, without a new contour; 0: off
    integer :: ef_scan = 0                             !! diagnostic: N(E) on 2 ef_scan + 1 energies at the start potential, then stop
    real(dp) :: ef_scan_step = 2.0e-5_dp, ef_scan_center = huge(1.0_dp)   !! scan step (Ry), centre (default: start E_F)
    character(len=256) :: ef_scan_dump = ''            !! diagnostic: file for the chain coefficients of every scan energy
    integer :: ndos = 0
    real(dp) :: dos_emin = -1.0_dp, dos_emax = 1.0_dp, dos_eta = 0.02_dp
    character(len=256) :: pot_in = '', pot_out = ''
    real(dp) :: ebot_offset = 1.2_dp                   !! the contour starts ebot_offset (Ry) below the lowest E_nu; res%wbelow checks it
    integer :: epi_shells = 0                          !! effective pair interactions (asr_epi) on this many shells after the SCF; 0: off
    integer :: epi_nlam = 4                            !! Gauss-Legendre nodes of the species path
    logical :: epi_sym = .false.                       !! one sublattice: compute U(c, s) for s >= c only, U(s, c) = U(c, s)
    type(spectral_input_t) :: spectral                 !! Bloch spectral function after the SCF (asr_spectral; route edc)
  end type class_scf_input_t

  type :: class_scf_result_t
    integer :: ncomp = 0, nspin = 1, iter = 0
    integer, allocatable :: comp_cls(:), comp_spec(:), nodes(:)
    real(dp), allocatable :: comp_w(:), comp_x(:)       !! weight per site f_kappa x; concentration on the sublattice x
    type(asa_class_result_t), allocatable :: comp(:)
    real(dp) :: ef = 0.0_dp, dv = 0.0_dp, moment_site = 0.0_dp, etot = 0.0_dp   !! etot: total energy per site (Ry)
    logical :: converged = .false.
    real(dp), allocatable :: hist_dv(:), hist_ef(:), hist_moment(:, :), hist_qnet(:, :), hist_etot(:)
    real(dp), allocatable :: edos_e(:), edos(:, :, :, :)   !! DOS(e, l, component, spin) per atom
    !> report data, from the last iteration.  eterm(n_et): the terms of etot per site (asr_runstats et_*; their sum is etot to
    !> rounding); ecomp(n_et, component): the same per atom of the component (eterm = sum_p w_p ecomp(:, p)).
    real(dp) :: eterm(n_et) = 0.0_dp
    real(dp), allocatable :: ecomp(:, :)
    !> dosef(component, spin, 0:2): DOS at E_F + i ed_eta per atom and per spin channel, s p d (states/Ry; nspin = 1 holds one
    !> channel, so the total is twice the sum); ed3: the DOS of the channel charge q0 + p q2.  qcell: the cell charge
    !> sum_p w_p qnet_p (e per site) that the Madelung sum removes (the residual of the E_F search; 0 = neutral).
    real(dp), allocatable :: dosef(:, :, :)
    real(dp) :: qcell = 0.0_dp
    !> per iteration, rows 1:iter: E_F passes (contour evaluations), wall seconds, largest change of E_nu before mixing (Ry)
    integer, allocatable :: hist_npass(:)
    real(dp), allocatable :: hist_wall(:), hist_denu(:)
    !> last iteration: mean recursion depth used per chain (edc) and bp tail fallbacks; Anderson restarts in the run, history used
    real(dp) :: depth_mean = 0.0_dp
    integer :: nfall = 0, mix_nrestart = 0, mix_nh = 0
    type(run_timer_t) :: timer                            !! wall / CPU time per part of the run (asr_runstats)
    real(dp), allocatable :: chain_rank_wall(:)           !! wall time of each rank's own chains, summed over the run (s)
    integer :: chain_split(2) = [0, 1]                    !! thread split of the last chain pass: teams (nout), threads per team
    type(screen_t) :: screen                               !! screened fluctuation electrostatics (asr_screening)
    real(dp), allocatable :: fscr(:, :)                    !! F(p, p'): potential of component p per unit deviation of p' (Ry/e);
                                                           !! huge if an SRO mode has Gamma = 0 (F infinite)
    real(dp), allocatable :: gscr(:, :), wscr(:, :)        !! Gamma = dq/db and W = V/db of the screening model
    real(dp), allocatable :: mode_lam(:), mode_gam(:)      !! SRO eigenmodes: eigenvalue of Ahat, Gamma of the mode
    logical :: uform = .false.                             !! SRO potential iterated in the drive u = db (screen_update)
    real(dp), allocatable :: udrive(:)                     !! u of the last iteration (e)
    real(dp), allocatable :: qch(:, :, :)                  !! qch(component, spin, l): sphere charge per channel, one spin channel
    !> E_nu channels held at their band window (clamp_enu) in the last iteration / in the whole run.  nclamp > 0: the E_nu fixed
    !> point lies outside the window, the result depends on the window margin, and converged is .false.
    integer :: nclamp = 0, nclamp_run = 0
    !> valence electrons per site in the 3 Ry below the contour bottom, at the first (wbelow1) and the final potential.  The
    !> occupied moments miss this weight.  It must vanish (|wbelow| < wbelow_tol); below_ok is that test.
    real(dp) :: wbelow = 0.0_dp, wbelow1 = 0.0_dp
    logical :: below_ok = .true.
    !> contour node values (node, start orbital, component, spin) of the last moment evaluation, and how many of them are not
    !> causal: -Im G <= 0 (edc) or the block -Im G not positive semidefinite (ed3).  A Herglotz function has none.  The
    !> truncated energy-dependent chain can have some (manual/src/theory/energy_dependent.md, causality); the charges then carry the chain error there.
    integer :: nnode = 0, nnoncausal = 0
    !> radius of the site clusters of the recursion (asr_crystal_graph; alat units) and their number of sites per sublattice
    real(dp) :: cluster_radius = 0.0_dp
    integer, allocatable :: cluster_sites(:)
    !> effective pair interactions (asr_epi), one entry per pair type; epi_wbelow: largest valence weight of the path site
    !> (electrons) in the 3 Ry below the contour bottom over all path points (must vanish, as wbelow)
    type(epi_table_t), allocatable :: epi(:)
    real(dp) :: epi_wbelow = 0.0_dp
    !> Bloch spectral function (asr_spectral), per graph site and per spin channel: spec_a(energy, k, l, spin) summed over the
    !> origin classes with their site fractions, spec_ac(energy, k, class, spin) the part of each origin class (sum = total).
    real(dp), allocatable :: spec_k(:, :), spec_dist(:), spec_e(:), spec_a(:, :, :, :), spec_ac(:, :, :, :)
    integer, allocatable :: spec_icorner(:)
    !> diagnostics: shell radii (alat) and the largest |Im G| of each shell relative to the on-site one; the largest difference
    !> of G_{0L,RL} and G_{0L,-R L} from the two origins (relative); the change of G when the chains are cut 4 steps earlier
    !> (relative; off-site and on-site); the number of (E, k, spin) points with A < 0 and the most negative A relative to the largest;
    !> the energy integral of A at the first and the last k point (per spin; 9 for a window that holds every band)
    real(dp), allocatable :: spec_rsh(:), spec_decay(:), spec_sum(:, :)
    !> ordered structure (one species per sublattice): eigenvalues (9 * nb, k, spin) and unfolding weights W_n(k) of H^(2)(k) at
    !> the path points, from the same structure constant blocks and parameters as the recursion (exact unfolded A of the operator)
    real(dp), allocatable :: spec_bands(:, :, :), spec_bwgt(:, :, :)
    real(dp) :: spec_rdef = 0.0_dp, spec_spread = 0.0_dp, spec_spread0 = 0.0_dp, spec_amin = 0.0_dp
    integer :: spec_nneg = 0, spec_nt = 0, spec_maxhop = 0, spec_hoplim = 0, spec_nsteps = 0
    integer, allocatable :: spec_nodes(:)
  end type class_scf_result_t

  integer, parameter :: nvert = 8
  integer, parameter :: lorb9(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]
  real(dp), parameter, public :: wbelow_tol = 1.0e-3_dp     !! electrons per site

  type :: class_op_t
    type(augcls_op) :: op
    real(dp), allocatable :: d(:, :, :)
  end type class_op_t

  !> The site cluster of one sublattice (asr_crystal_graph: start site on the sublattice) with what the operators need
  type :: sub_graph_t
    type(crystal_graph_t) :: cg
    logical, allocatable :: flip(:)                       !! flippable sites (graph distance <= flip_radius, not the start)
    integer, allocatable :: dist(:)                       !! graph distance from the start site
    integer, allocatable :: cls(:)                        !! sublattice of every site
    integer, allocatable :: stype(:)                      !! structure-constant table of every site (translation type)
    integer, allocatable :: perm(:, :)                    !! site operations (graph_site_ops): images of every site
    real(dp), allocatable :: rot(:, :, :)                 !! their Cartesian rotations
  end type sub_graph_t

contains

  subroutine class_scf(inp, res, verbose)
    type(class_scf_input_t), intent(in) :: inp
    type(class_scf_result_t), intent(out) :: res
    logical, intent(in), optional :: verbose
    type(xc_lda_t) :: xf
    type(asa_input_t) :: ai
    type(potpar_t), allocatable :: pps(:, :)                ! (component, spin)
    type(potpar_l_t) :: ppl
    type(structure_t) :: st
    type(symop_t), allocatable :: sops(:)
    type(sub_graph_t), allocatable :: sgr(:)
    type(class_alloy_t) :: mdl
    type(class_op_t), allocatable :: cop(:)
    type(edn_t) :: ed
    type(anderson_t) :: mixer
    type(asa_class_result_t), allocatable :: cr(:)
    integer, allocatable :: offx(:), cls(:), pcl(:), psp(:), styp(:), gidx(:), red(:)
    integer, allocatable :: norb(:), iorbk(:, :), lorbk(:, :)   !! start orbitals per sublattice: count, orbital, l
    real(dp), allocatable :: worbk(:, :), dgo(:, :, :)          !! their weights (orbital group sizes); orbital matrices
    logical, allocatable :: mono(:)
    real(dp), allocatable :: y(:, :, :, :, :, :), enu(:, :, :), mom(:, :, :, :), dosef(:, :, :), w(:), xp(:), phimad(:, :), qbar(:)
    real(dp), allocatable :: sfrac(:)
    integer :: ksro
    complex(dp), allocatable :: dmp_a(:, :, :, :, :), dmp_b(:, :, :, :, :)   ! ef_scan_dump: chain coefficients (level, node, orb, comp, spin)
    integer, allocatable :: dmp_m(:, :, :, :)
    real(dp), allocatable :: xin(:), xout(:), xnext(:), wmix(:), nval(:, :), ncore(:, :), ntot(:), vh(:), vxc(:, :), rho(:, :), exc(:)
    real(dp), allocatable :: tmp(:), cum(:), yc(:, :, :), vtmp(:, :), etmp(:, :), esum(:)
    real(dp) :: alpha(0:2), c, avw, rmin, fs, nelec, ef, nn, dvc, dv, e, elo, ehi, flo, fhi, vm, qcell, rlo, rhi, etot, ep, br
    real(dp) :: ef_in, r2n, denu, rdv, efshift
    integer :: npass, efmode, pdv, sgdv, kdv
    real(dp), allocatable :: molo(:, :, :, :), mohi(:, :, :, :), dslo(:, :, :), dshi(:, :, :)
    integer, allocatable :: pk(:), zsh(:, :)
    real(dp), allocatable :: ahat(:, :)
    real(dp), allocatable :: dqv(:), ucur(:), uout(:), vfl(:)
    logical :: uload
    real(dp) :: tw0, tc0, tw_it
    real(dp), allocatable :: ecp(:, :)
    integer :: offu
    character(len=:), allocatable :: smsg
    integer :: p2
    integer :: nsp, ncomp, ncls, p, kap, k, n, nmax, it, i, sg, l, nx, jf, nfe, side, ngen, nlev, is, nod, nmc, nomax, b0
    integer :: grp(9), ngrp, rep(9), wgt(9)
    logical :: vb, hlo, hhi, core_ok, ed3
    integer :: teams_said(2)                    !! the last thread split of chain_greens that was printed (nout, nin)
    real(dp), allocatable :: qch(:, :, :), dsm(:, :, :, :), qclo(:, :, :), qchi(:, :, :), dmlo(:, :, :, :), dmhi(:, :, :, :)
    real(dp) :: phi2, asym_max, psd_max
    integer :: npsd, nnode

    ed3 = trim(inp%route) == 'ed3'
    if (ed3 .and. len_trim(inp%ef_scan_dump) > 0) error stop 'class_scf: ef_scan_dump is not available with route = ed3'
    phi2 = merge(0.0_dp, 1.0_dp, ed3)                  ! weight of phi phiddot in the density
    vb = .false.; if (present(verbose)) vb = verbose
    teams_said = 0
    npsd = 0; nnode = 0
    nsp = inp%nspin; fs = 2.0_dp / real(nsp, dp)
    c = c_nonrel; if (inp%relativistic) c = c_ry
    if (inp%st%nb > 0) then
      st = inp%st
    else
      call decoration_structure(inp%dec, st, smsg)
      if (len(smsg) > 0) error stop smsg
    end if
    ncls = st%nlabel
    allocate (sfrac(ncls))
    do kap = 1, ncls
      sfrac(kap) = real(count(st%label == kap), dp) / real(st%nb, dp)
    end do
    call space_group(st, sops)
    avw = ws_radius(st) * inp%alat
    alpha = standard_alpha() * inp%alpha_scale
    ksro = 0
    if (inp%sro) then
      do kap = 1, ncls
        if (any(abs(inp%sro_cond(1:inp%nsc(kap), 1:inp%nsc(kap), kap, kap) - &
                    spread(inp%conc(1:inp%nsc(kap), kap), 1, inp%nsc(kap))) > 1.0e-12_dp)) ksro = kap
      end do
    end if
    call screen_structure(st, inp%alat, inp%sim_beta, res%screen, smsg, sub_sro=ksro)
    if (len(smsg) > 0) error stop smsg
    ngen = inp%nlevels
    nlev = inp%nlevels; if (inp%nrec > 0) nlev = inp%nrec
    tw0 = wall_seconds(); tc0 = cpu_seconds()
    call timer_start(res%timer, tb_setup)
    ! ---- components
    ncomp = sum(inp%nsc)
    allocate (pcl(ncomp), psp(ncomp), pk(ncomp), w(ncomp), xp(ncomp))
    p = 0
    do kap = 1, ncls
      do k = 1, inp%nsc(kap)
        p = p + 1; pcl(p) = kap; pk(p) = k; psp(p) = inp%spc(k, kap); xp(p) = inp%conc(k, kap)
        w(p) = sfrac(kap) * inp%conc(k, kap)
      end do
    end do
    nelec = sum([(w(p) * real(inp%spec(psp(p))%nval, dp), p=1, ncomp)])
    allocate (ucur(ncomp), uout(ncomp), vfl(ncomp), dqv(ncomp)); ucur = 0.0_dp; uout = 0.0_dp; vfl = 0.0_dp
    uload = .false.                                    ! u read from start_potential; otherwise u starts at zero
    ! operator classes: the sublattices, and with SRO also ncls + kappa' = first-shell sites of the origin on sublattice kappa'
    ! (same species; their concentrations are set per central component in chain_worker)
    nmc = ncls; if (inp%sro) nmc = 2 * ncls
    mdl%nl = 9; mdl%nsp = ncomp; mdl%ncls = nmc
    allocate (mdl%nsc(nmc), mdl%spc(maxval(inp%nsc), nmc), mdl%conc(maxval(inp%nsc), nmc))
    mdl%spc = 0; mdl%conc = 0.0_dp
    p = 0
    do kap = 1, ncls
      mdl%nsc(kap) = inp%nsc(kap)
      do k = 1, inp%nsc(kap)
        p = p + 1; mdl%spc(k, kap) = p; mdl%conc(k, kap) = inp%conc(k, kap)
      end do
    end do
    if (inp%sro) then
      mdl%nsc(ncls + 1:) = mdl%nsc(1:ncls); mdl%spc(:, ncls + 1:) = mdl%spc(:, 1:ncls); mdl%conc(:, ncls + 1:) = mdl%conc(:, 1:ncls)
    end if
    call structure_madelung(st, inp%alat, phimad)
    allocate (qbar(ncls))
    call xc_init(xf, inp%xc)
    ! ---- meshes and atomic starts (equal spheres)
    allocate (cr(ncomp), enu(0:2, 2, ncomp), pps(ncomp, 2)); enu = 0.0_dp
    nmax = 0
    do p = 1, ncomp
      call make_radmesh(inp%spec(psp(p))%z, avw, 0.012_dp, cr(p)%mesh, c)
      n = cr(p)%mesh%n; nmax = max(nmax, n)
      allocate (cr(p)%v(n, nsp), cr(p)%nrad(n, nsp))
      ai%alat = inp%alat; ai%xc = inp%xc; ai%relativistic = inp%relativistic; ai%nspin = nsp; ai%zeta0 = inp%zeta0
      if (allocated(ai%cls)) deallocate (ai%cls)
      allocate (ai%cls(1)); ai%cls(1) = inp%spec(psp(p))
      allocate (vtmp(n, nsp), etmp(0:2, nsp))
      call timer_stop(res%timer, tb_setup); call timer_start(res%timer, tb_atom)
      call atomic_start(ai, 1, cr(p)%mesh, c, xf, nsp, vtmp, etmp)
      call timer_stop(res%timer, tb_atom); call timer_start(res%timer, tb_setup)
      cr(p)%v = vtmp
      do sg = 1, nsp; enu(:, sg, p) = etmp(:, sg); end do
      deallocate (vtmp, etmp)
    end do
    ef_in = huge(1.0_dp)
    if (len_trim(inp%pot_in) > 0) call read_potentials()
    allocate (y(2, 0:2, nmax, 0:2, ncomp, 2)); y = 0.0_dp
    call timer_stop(res%timer, tb_setup)
    call all_params()
    call timer_start(res%timer, tb_setup)
    ! ---- site clusters (one per sublattice, start site on it; asr_crystal_graph), S blocks per translation type, one class
    ! operator per sublattice.  The cluster holds every site the closure reaches in ngen levels (and the spectral targets).
    rmin = 0.0_dp; if (inp%spectral%on) rmin = inp%spectral%rmax
    allocate (sgr(ncls))
    do kap = 1, ncls
      b0 = findloc(st%label, kap, 1)
      call crystal_graph(st, sops, inp%rcl, ngen, sgr(kap)%cg, smsg, start=b0, rmin=rmin)
      if (len(smsg) > 0) error stop smsg
    end do
    call structure_blocks(st, sgr(1)%cg, max(inp%rcl, inp%rinv), alpha, mdl%block, mdl%block0, styp)
    allocate (mdl%c(9, ncomp), mdl%d(9, ncomp), mdl%o(9, ncomp), mdl%e(9, ncomp))
    call edc_setup(mdl, pps(:, 1), alpha, 1, ed)
    allocate (cop(ncls), res%nodes(ncls), zsh(ncls, ncls)); zsh = 0
    allocate (norb(ncls), iorbk(9, ncls), lorbk(9, ncls), worbk(9, ncls)); iorbk = 0; lorbk = 0; worbk = 0.0_dp
    do kap = 1, ncls
      associate (sk => sgr(kap))
        call graph_distances(sk%cg%g, 1, sk%dist)
        sk%flip = sk%dist <= inp%flip_radius; sk%flip(1) = .false.
        sk%cls = st%label(sk%cg%basis)
        sk%stype = styp(sk%cg%basis)
        call graph_site_ops(sk%cg, st, sops, sk%perm, sk%rot, gidx)
      end associate
      cls = sgr(kap)%cls
      if (inp%sro) then
        call shell_classes(kap, cls)
        do k = 1, ncls; zsh(kap, k) = count(cls == ncls + k); end do
        do k = 1, ncls                                 ! SRO with a sublattice that is not in the first shell: input error
          if (.not. any(cls == ncls + k)) then
            do i = 1, inp%nsc(kap)
              if (any(abs(inp%sro_cond(i, 1:inp%nsc(k), kap, k) - inp%conc(1:inp%nsc(k), k)) > 1.0e-12_dp)) &
                error stop 'class_scf: [[sro]] for a sublattice pair that are not first-shell neighbours'
            end do
          else if (vb) then
            do i = 1, inp%nsc(kap)
              print '(a,i0,a,i0,a,i0,a,i0,a,*(f8.4))', '  SRO: sublattice ', kap, ' component ', i, ', first-shell sites ', &
                count(cls == ncls + k), ' on sublattice ', k, ': P(s | centre) =', inp%sro_cond(i, 1:inp%nsc(k), kap, k)
            end do
          end if
        end do
      end if
      ! start orbitals: one per orbital group of the class-preserving site group, weight = group size; the operator is reduced
      ! by its monomial operations under which every start orbital is one dimensional (asr_crystal_symmetry)
      call class_group(kap, cls, gidx, dgo, mono)
      call orbital_groups(dgo, mono, grp, ngrp, rep, wgt)
      norb(kap) = ngrp; iorbk(1:ngrp, kap) = rep(1:ngrp); lorbk(1:ngrp, kap) = lorb9(rep(1:ngrp))
      worbk(1:ngrp, kap) = real(wgt(1:ngrp), dp)
      call reduction_ops(dgo, mono, rep(1:ngrp), red)
      cop(kap)%d = dgo(:, :, red)
      call build_augcls_operator(sgr(kap)%cg%g, mdl, cls, sgr(kap)%perm(:, gidx(red)), cop(kap)%d, 1, ngen, cop(kap)%op, &
                                 inp%maxrank, sgr(kap)%flip, stype=sgr(kap)%stype)
      res%nodes(kap) = cop(kap)%op%nn
      if (vb) print '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)', '  sublattice ', kap, ': cluster sites ', sgr(kap)%cg%n, &
        ', site group ', size(gidx), ', start orbitals ', ngrp, ', reduction group ', size(red), ', components ', inp%nsc(kap), &
        ', augmented-space nodes ', cop(kap)%op%nn
    end do
    nomax = maxval(norb)
    res%cluster_radius = sgr(1)%cg%radius
    res%cluster_sites = [(sgr(kap)%cg%n, kap=1, ncls)]
    if (vb) print '(a,i0,a,i0,a,i0,a,i0,a,i0)', '  class ASR: components=', ncomp, '  levels=', ngen, '  recursion steps=', nlev, &
      '  rank<=', inp%maxrank, '  nspin=', nsp
    ! ---- screened fluctuation electrostatics: F = (G00 + G1 Ahat)^-1 (W00 + W1 Ahat)
    allocate (ahat(ncomp, ncomp), res%fscr(ncomp, ncomp)); ahat = 0.0_dp
    if (inp%sro) then
      do p = 1, ncomp
        do p2 = 1, ncomp
          ahat(p, p2) = zsh(pcl(p), pcl(p2)) * (inp%sro_cond(pk(p), pk(p2), pcl(p), pcl(p2)) - inp%conc(pk(p2), pcl(p2)))
        end do
      end do
    end if
    allocate (res%gscr(ncomp, ncomp), res%wscr(ncomp, ncomp), res%mode_lam(ncomp), res%mode_gam(ncomp))
    res%uform = inp%sro .and. .not. res%screen%unscreened
    call screen_matrix(res%screen, ahat, w, gm=res%gscr, wm=res%wscr, lam=res%mode_lam, gam=res%mode_gam)
    if (minval(res%mode_gam) > 1.0e-9_dp * res%screen%g00) then
      call screen_matrix(res%screen, ahat, w, fm=res%fscr)
    else
      res%fscr = huge(1.0_dp)                          ! F infinite on a mode: only the drive form is defined
    end if
    if (vb .and. inp%sro) then
      print '(a)', '  SRO modes: eigenvalue lambda of Ahat, Gamma = dq/db (closure), G00 + G1 lambda (first-shell form), F (Ry/e)'
      do p = 1, ncomp
        if (res%mode_gam(p) > 1.0e-9_dp * res%screen%g00) then
          print '(4x,f9.4,2f11.5,f12.4)', res%mode_lam(p), res%mode_gam(p), res%screen%g00 + res%screen%g1 * res%mode_lam(p), &
            mode_f(res%screen, res%mode_lam(p), res%mode_gam(p))
        else
          print '(4x,f9.4,2f11.5,a)', res%mode_lam(p), res%mode_gam(p), res%screen%g00 + res%screen%g1 * res%mode_lam(p), &
            '    infinite (full clustering: this charge deviation vanishes)'
        end if
      end do
    end if
    ! ---- mixing vector: per component and spin the smooth potential (n) and E_nu (3); with SRO then the drive u (ncomp)
    allocate (offx(ncomp + 1)); offx(1) = 0
    do p = 1, ncomp; offx(p + 1) = offx(p) + (cr(p)%mesh%n + 3) * nsp; end do
    nx = offx(ncomp + 1); offu = nx
    if (res%uform) nx = nx + ncomp
    allocate (xin(nx), xout(nx), xnext(nx))
    call anderson_init(mixer, nx, inp%mix_history)
    allocate (mom(ncomp, 2, 0:2, 0:2), dosef(ncomp, 2, 0:2), qch(ncomp, 2, 0:2), dsm(ncomp, 2, 0:2, 0:2))
    ed_depth_tol = inp%depth_tol; ed_tail = inp%ed_tail
    ef = maxval(enu(2, 1:nsp, :)) + 0.35_dp
    if (ef_in < huge(1.0_dp)) ef = ef_in                 ! E_F of the start potential's run
    allocate (res%hist_dv(inp%maxit), res%hist_ef(inp%maxit), res%hist_moment(inp%maxit, ncomp), res%hist_qnet(inp%maxit, ncomp), &
              res%hist_etot(inp%maxit), esum(ncomp))
    allocate (res%hist_npass(inp%maxit), res%hist_wall(inp%maxit), res%hist_denu(inp%maxit), ecp(n_et, ncomp), &
              res%dosef(ncomp, nsp, 0:2))
    res%hist_npass = 0; res%hist_wall = 0.0_dp; res%hist_denu = 0.0_dp; ecp = 0.0_dp; res%dosef = 0.0_dp
    br = inp%beta; if (inp%beta_restart > 0.0_dp) br = inp%beta_restart
    call timer_stop(res%timer, tb_setup)

    do it = 1, inp%maxit
      ed_nchain = 0; ed_nlev = 0; ed_nfall = 0; res%nclamp = 0
      tw_it = wall_seconds()
      do p = 1, ncomp
        n = cr(p)%mesh%n
        do sg = 1, nsp
          i = offx(p) + (sg - 1) * (n + 3)
          xin(i + 2:i + n) = cr(p)%v(2:n, sg) + 2.0_dp * inp%spec(psp(p))%z / cr(p)%mesh%r(2:n)
          xin(i + 1) = xin(i + 2)
          xin(i + n + 1:i + n + 3) = enu(:, sg, p)
        end do
      end do
      if (res%uform) xin(offu + 1:offu + ncomp) = ucur
      call all_params()
      if (vb .and. it == 1) then
        call ed_moments(6.0_dp, mom, dosef, qch, dsm)
        do p = 1, ncomp
          print '(a,i0,a,3f9.4)', '  sum rule (E_F=6 Ry), component ', p, ': M0(s,p,d)/(1,3,5) =', (qch(p, 1, l) / real(2 * l + 1, dp), l=0, 2)
          if (ed3) print '(a,i0,a,3f10.5)', '  ed3 component ', p, ': p^gamma(s,p,d) of this potential =', pps(p, 1)%p
        end do
        if (ed3) print '(a,f9.4,a)', '  ed3 sum rule: weighted channel charge per orbital =', &
          sum([(w(p) * sum(qch(p, 1, :)) / 9.0_dp, p=1, ncomp)]), ' (1 only if all pencil eigenvalues lie below 6 Ry; single channels differ from 1)'
        print '(a)', '  (this sum rule runs the production contour up to 6 Ry: a rough check, a few per cent off at 32 contour points)'
      end if
      if (it == 1 .and. inp%ef_scan == 0) then
        call weight_below(res%wbelow1)
        if (vb .or. abs(res%wbelow1) >= wbelow_tol) print '(a,f6.2,a,es10.2,a)', '  valence weight within 3 Ry below the contour bottom (', &
          inp%ebot_offset, ' Ry below the lowest E_nu):', res%wbelow1, ' electrons per site (must vanish)'
        if (abs(res%wbelow1) >= wbelow_tol) print '(a)', '  WARNING: states below the contour bottom at the start potential'
      end if
      if (inp%ef_scan > 0) then
        call ef_scan_run()
        stop 'class_scf: ef_scan done'
      end if
      ! E_F: Newton steps until N - N_elec changes sign, then Illinois false position on the bracket (as asr_alloy_scf), stopped at
      ! |N - N_elec| < 2e-4 or when the Newton step is below 1e-3 Ry (one contour per iteration close to self-consistency: the
      ! first-order step below is exact at the fixed point, where E_F no longer moves between iterations).  The electron count is then made exact (see below): the residual left by the search went straight into
      ! the sphere charges (self-Hartree 2 q/S: up to ~1.5e-4 Ry) and set an SCF noise floor above tol whenever E_F lies in a gap
      ! (Co2CrAl: dV ~1e-4 for 120 iterations at converged moments).
      nfe = 16
      hlo = .false.; hhi = .false.; side = 0
      do jf = 1, nfe
        call ed_moments(ef, mom, dosef, qch, dsm)
        nn = fs * sum([(w(p) * sum(qch(p, 1:nsp, :)), p=1, ncomp)])
        dvc = fs * sum([(w(p) * sum(dosef(p, 1:nsp, :)), p=1, ncomp)])
        if (vb .and. it == 1) print '(a,i2,a,f10.6,a,f9.5,a,f9.4)', '    E_F search ', jf, ': E_F=', ef, '  N-nelec=', nn - nelec, '  DOS(E_F)=', dvc
        if (nn < nelec) then
          if (side == -1 .and. hhi) fhi = 0.5_dp * fhi
          elo = ef; flo = nn - nelec; hlo = .true.; side = -1
          rlo = nn - nelec; molo = mom; dslo = dosef; qclo = qch; dmlo = dsm
        else
          if (side == 1 .and. hlo) flo = 0.5_dp * flo
          ehi = ef; fhi = nn - nelec; hhi = .true.; side = 1
          rhi = nn - nelec; mohi = mom; dshi = dosef; qchi = qch; dmhi = dsm
        end if
        npass = jf
        if (abs(nn - nelec) < max(2.0e-4_dp, 1.0e-3_dp * dvc) .or. jf == nfe) exit
        ! one contour per iteration (ef_linear_max > 0): a first-order step up to ef_linear_max replaces the next contour.  The
        ! electron count stays exact (below).  The moments carry an error of first order in the step: the step uses the DOS at
        ! E_F + i ed_eta, while the derivative of the contour moments is the thermally broadened DOS (fermi_kt), and the two
        ! differ per component and spin.  The error vanishes at self-consistency, where E_F no longer moves between iterations
        ! (the carried E_F is a quasi-Newton iterate of N(E_F) = N_elec).  So the fixed point is unchanged; only the path
        ! towards it differs.  The SCF stops on dV only: check the printed shift of the last iteration.
        if (jf == 1 .and. inp%ef_linear_max > 0.0_dp) then
          if (abs(nelec - nn) / max(dvc, 0.1_dp) <= inp%ef_linear_max) exit
        end if
        if (hlo .and. hhi) then
          e = (elo * fhi - ehi * flo) / (fhi - flo)
          if (.not. (e > min(elo, ehi) + 1.0e-3_dp * abs(ehi - elo) .and. e < max(elo, ehi) - 1.0e-3_dp * abs(ehi - elo))) e = 0.5_dp * (elo + ehi)
          e = e - ef
        else
          e = max(-0.15_dp, min(0.15_dp, (nelec - nn) / max(dvc, 1.0_dp)))
        end if
        ef = ef + e
      end do
      if (abs(nn - nelec) >= max(2.0e-4_dp, 1.0e-3_dp * dvc) .and. vb) print '(a,i0,a,es10.3)', '    warning: E_F search not converged after ', nfe, &
        ' passes: N - N_elec =', nn - nelec
      efmode = 0; efshift = 0.0_dp
      ! exact electron count: N is linear in the moments, so on a narrow bracket the moments interpolated linearly between its two
      ! evaluations give N = N_elec exactly (error O(width^2) in the moments, no extra recursion); otherwise a first-order step
      ! with the DOS at E_F
      if (hlo .and. hhi .and. abs(ehi - elo) < 0.02_dp) then
        e = rlo / (rlo - rhi)
        efshift = elo + e * (ehi - elo) - ef; efmode = 1
        ef = elo + e * (ehi - elo)
        mom = molo + e * (mohi - molo); dosef = dslo + e * (dshi - dslo)
        qch = qclo + e * (qchi - qclo); dsm = dmlo + e * (dmhi - dmlo)
      else if (abs(nn - nelec) > 0.0_dp) then
        e = max(-0.15_dp, min(0.15_dp, (nelec - nn) / max(dvc, 0.1_dp)))
        if (ed3) then                                  ! the q_k are matrix elements, not energy moments: each has its own DOS
          mom = mom + e * dsm; qch = qch + e * dosef
        else
          do p = 1, ncomp; do sg = 1, nsp; do l = 0, 2; do k = 0, 2
            mom(p, sg, k, l) = mom(p, sg, k, l) + e * (ef + 0.5_dp * e - enu(l, sg, p))**k * dosef(p, sg, l)
          end do; end do; end do; end do
          qch = mom(:, :, 0, :)
        end if
        ef = ef + e; efshift = e; efmode = 2
      end if
      ! densities, core, charges; esum = band + core eigenvalue sums - int n v~_in (kinetic energy + int n 2Z/r)
      call timer_start(res%timer, tb_density)
      dv = 0.0_dp; pdv = 0; sgdv = 0; kdv = 0; denu = 0.0_dp
      do p = 1, ncomp
        n = cr(p)%mesh%n
        allocate (nval(n, nsp), ncore(n, nsp), ntot(n), tmp(n), cum(n), yc(2, 0:0, n))
        nval = 0.0_dp; ncore = 0.0_dp; esum(p) = 0.0_dp; ecp(:, p) = 0.0_dp
        do sg = 1, nsp
          do l = 0, 2
            ! band energy of the channel: E_nu times the sphere charge (edc: M0; ed3: q0 + p^gamma q2, since
            ! e z^T O z = z^T H z = sum_RL [q1 + E_nu (q0 + p^gamma q2)]; test_band_energy) plus the first moment
            esum(p) = esum(p) + fs * (mom(p, sg, 1, l) + enu(l, sg, p) * qch(p, sg, l))
            ecp(et_band, p) = ecp(et_band, p) + fs * (mom(p, sg, 1, l) + enu(l, sg, p) * qch(p, sg, l))
            nval(:, sg) = nval(:, sg) + fs * ( &
                mom(p, sg, 0, l) * (y(1, 0, 1:n, l, p, sg)**2 + y(2, 0, 1:n, l, p, sg)**2 / c**2) &
              + 2.0_dp * mom(p, sg, 1, l) * (y(1, 0, 1:n, l, p, sg) * y(1, 1, 1:n, l, p, sg) + y(2, 0, 1:n, l, p, sg) * y(2, 1, 1:n, l, p, sg) / c**2) &
              + mom(p, sg, 2, l) * (y(1, 1, 1:n, l, p, sg)**2 + y(2, 1, 1:n, l, p, sg)**2 / c**2 &
                                    + phi2 * (y(1, 0, 1:n, l, p, sg) * y(1, 2, 1:n, l, p, sg) + y(2, 0, 1:n, l, p, sg) * y(2, 2, 1:n, l, p, sg) / c**2)))
          end do
          do is = 1, size(inp%spec(psp(p))%core)
            nod = inp%spec(psp(p))%core(is)%n - inp%spec(psp(p))%core(is)%l - 1
            call bound_state(cr(p)%mesh, cr(p)%v(:, sg), inp%spec(psp(p))%core(is)%l, nod, inp%spec(psp(p))%z, c, n, &
                             -1.6_dp * inp%spec(psp(p))%z**2 - 20.0_dp, cr(p)%v(n, sg) - 1.0e-6_dp, e, yc, core_ok)
            ncore(:, sg) = ncore(:, sg) + inp%spec(psp(p))%core(is)%occ / real(nsp, dp) * (yc(1, 0, :)**2 + yc(2, 0, :)**2 / c**2)
            esum(p) = esum(p) + inp%spec(psp(p))%core(is)%occ / real(nsp, dp) * e
            ecp(et_core, p) = ecp(et_core, p) + inp%spec(psp(p))%core(is)%occ / real(nsp, dp) * e
          end do
        end do
        cr(p)%nrad = nval + ncore
        do sg = 1, nsp
          tmp(1) = 0.0_dp
          tmp(2:n) = cr(p)%nrad(2:n, sg) * (cr(p)%v(2:n, sg) + 2.0_dp * inp%spec(psp(p))%z / cr(p)%mesh%r(2:n))
          call cumint(cr(p)%mesh, tmp, cum); esum(p) = esum(p) - cum(n); ecp(et_dc, p) = ecp(et_dc, p) - cum(n)
        end do
        ntot = sum(cr(p)%nrad, dim=2)
        tmp = ntot; call cumint(cr(p)%mesh, tmp, cum)
        cr(p)%qnet = cum(n) - inp%spec(psp(p))%z
        cr(p)%moment = 0.0_dp
        if (nsp == 2) then
          tmp = cr(p)%nrad(:, 1) - cr(p)%nrad(:, 2); call cumint(cr(p)%mesh, tmp, cum); cr(p)%moment = cum(n)
        end if
        deallocate (nval, ncore, ntot, tmp, cum, yc)
      end do
      qbar = 0.0_dp
      do p = 1, ncomp; qbar(pcl(p)) = qbar(pcl(p)) + xp(p) * cr(p)%qnet; end do
      qcell = sum(sfrac * qbar)            ! residual of the E_F search: the Madelung sum needs neutral charges
      ! SRO: potential of the charge deviations from the drive u (screen_update).  u_out = dq + (I - Gamma) u_in is the bare
      ! charge deviation of the model, dq + v_in / a, if the input potential holds the SRO potential v_in = W u_in.  A start
      ! without u (free atoms, a potential of a run without SRO) holds none, so u starts at zero.  The first SRO potential is
      ! then W dq, at most a |dq|, and it follows the charges.  (Until 2026-11 the first u was Gamma^+ dq.  That stands for
      ! an input potential that holds F dq.  From free atoms the first charges are several electrons, the drive moves by
      ! 0.15 / a per iteration only, and W (I - Gamma) u kept several Ry of them in every later output potential.)
      if (res%uform) then
        do p = 1, ncomp; dqv(p) = cr(p)%qnet - qbar(pcl(p)); end do
        call screen_update(res%gscr, res%wscr, xp, pcl, dqv, ucur, uout, vfl)
        xout(offu + 1:offu + ncomp) = uout
      end if
      ! new potentials: sphere Hartree + xc + Madelung (sublattice averages) + SIM (deviation from the sublattice average)
      etot = 0.0_dp
      do p = 1, ncomp
        n = cr(p)%mesh%n
        allocate (ntot(n), vh(n), vxc(n, 2), rho(n, 2), exc(n))
        ntot = sum(cr(p)%nrad, dim=2)
        call poisson_hartree(cr(p)%mesh, ntot, vh)
        do k = 2, n; rho(k, 1:nsp) = cr(p)%nrad(k, 1:nsp) / (4.0_dp * pi * cr(p)%mesh%r(k)**2); end do
        rho(1, :) = rho(2, :)
        if (nsp == 1) then
          call xc_eval(xf, rho(:, 1), exc, vxc(:, 1), cr(p)%mesh%r); vxc(:, 2) = vxc(:, 1)
        else
          call xc_eval_pol(xf, rho(:, 1), rho(:, 2), exc, vxc(:, 1), vxc(:, 2), cr(p)%mesh%r)
        end if
        vm = dot_product(phimad(pcl(p), :), qbar - qcell)
        if (res%uform) then
          vm = vm + vfl(p)
        else
          do p2 = 1, ncomp; vm = vm + res%fscr(p, p2) * (cr(p2)%qnet - qbar(pcl(p2))); end do
        end if
        ecp(et_mad, p) = 0.5_dp * cr(p)%qnet * dot_product(phimad(pcl(p), :), qbar - qcell)
        ecp(et_sim, p) = 0.5_dp * cr(p)%qnet * vm - ecp(et_mad, p)
        allocate (tmp(n), cum(n))
        tmp = ntot * (0.5_dp * vh + exc); tmp(1) = 0.0_dp
        call cumint(cr(p)%mesh, tmp, cum)
        ep = esum(p) + cum(n) + 0.5_dp * cr(p)%qnet * vm
        etot = etot + w(p) * ep
        tmp = 0.5_dp * ntot * vh; tmp(1) = 0.0_dp
        call cumint(cr(p)%mesh, tmp, cum); ecp(et_hartree, p) = cum(n)
        tmp = ntot * exc; tmp(1) = 0.0_dp
        call cumint(cr(p)%mesh, tmp, cum); ecp(et_xc, p) = cum(n)
        deallocate (tmp, cum)
        do sg = 1, nsp
          i = offx(p) + (sg - 1) * (n + 3)
          xout(i + 2:i + n) = vh(2:n) + vxc(2:n, sg) + vm
          xout(i + 1) = xout(i + 2)
          do k = 2, n
            rdv = abs(xout(i + k) - xin(i + k)) * cr(p)%mesh%r(k) / (1.0_dp + cr(p)%mesh%r(k))
            if (rdv > dv) then; dv = rdv; pdv = p; sgdv = sg; kdv = k; end if
          end do
          do l = 0, 2
            xout(i + n + 1 + l) = enu(l, sg, p) + 0.7_dp * mom(p, sg, 1, l) / max(mom(p, sg, 0, l), 1.0e-6_dp)
            denu = max(denu, abs(xout(i + n + 1 + l) - xin(i + n + 1 + l)))
          end do
        end do
        deallocate (ntot, vh, vxc, rho, exc)
      end do
      call timer_stop(res%timer, tb_density)
      call timer_start(res%timer, tb_mix)
      if (.not. allocated(wmix)) then           ! Anderson metric = the dV measure (r/(1+r))^2 on the potentials, 1 on E_nu
        allocate (wmix(nx)); wmix = 1.0_dp
        do p = 1, ncomp
          n = cr(p)%mesh%n
          do sg = 1, nsp
            i = offx(p) + (sg - 1) * (n + 3)
            wmix(i + 1:i + n) = (cr(p)%mesh%r(1:n) / (1.0_dp + cr(p)%mesh%r(1:n)))**2
          end do
        end do
        if (res%uform) wmix(offu + 1:offu + ncomp) = res%screen%a**2     ! u measured as the potential a u (Ry)
      end if
      ! the drive is part of the fixed point: its residual a |u_out - u_in| (Ry) enters dV
      if (res%uform) dv = max(dv, res%screen%a * maxval(abs(uout - ucur)))
      ! restart safeguards only far from self-consistency (dV > 1e-2); closer the map is linear (Co2Fe0.4Cr0.6Al: the Cr response
      ! is linear to 2e-5 e for 5e-4 Ry shifts) and restarts every few iterations stalled the SCF at dV ~ 1e-4; Fe66Ni14Cr20 at
      ! levels 4 restarted at dV ~ 1.5e-3 (threshold 1e-3) 15 times in 35 iterations and never converged
      r2n = sqrt(sum(wmix * (xout - xin)**2))
      call anderson_step(mixer, xin, xout, inp%beta, xnext, wmix, safeguard=dv > 1.0e-2_dp, beta_restart=br)
      do p = 1, ncomp
        n = cr(p)%mesh%n
        do sg = 1, nsp
          i = offx(p) + (sg - 1) * (n + 3)
          do k = 1, n
            xnext(i + k) = xin(i + k) + max(-0.15_dp, min(0.15_dp, xnext(i + k) - xin(i + k)))
          end do
          cr(p)%v(2:n, sg) = xnext(i + 2:i + n) - 2.0_dp * inp%spec(psp(p))%z / cr(p)%mesh%r(2:n)
          cr(p)%v(1, sg) = cr(p)%v(2, sg)
          do l = 0, 2
            enu(l, sg, p) = xin(i + n + 1 + l) + max(-0.08_dp, min(0.08_dp, xnext(i + n + 1 + l) - xin(i + n + 1 + l)))
            if (clamp_enu(cr(p)%mesh, cr(p)%v(:, sg), l, inp%spec(psp(p))%nvl(l) - l - 1, inp%spec(psp(p))%z, c, enu(l, sg, p))) then
              res%nclamp = res%nclamp + 1; res%nclamp_run = res%nclamp_run + 1
              if (vb) print '(a,i0,a,i0,a,i0,a,f9.5)', '    E_nu outside its band window: component ', p, ' spin ', sg, ' l ', l, ' -> ', enu(l, sg, p)
            end if
          end do
        end do
      end do
      if (res%uform) then                          ! drive: the step limit of the potentials, 0.15 Ry, as a u step 0.15 / a
        do p = 1, ncomp
          ucur(p) = xin(offu + p) + max(-0.15_dp, min(0.15_dp, res%screen%a * (xnext(offu + p) - xin(offu + p)))) / res%screen%a
        end do
        call screen_project(xp, pcl, ucur)
      end if
      call timer_stop(res%timer, tb_mix)
      res%iter = it; res%dv = dv; res%ef = ef; res%etot = etot
      res%nnode = nnode; res%nnoncausal = npsd
      res%hist_dv(it) = dv; res%hist_ef(it) = ef; res%hist_etot(it) = etot
      res%hist_moment(it, :) = cr(:)%moment; res%hist_qnet(it, :) = cr(:)%qnet
      res%hist_npass(it) = npass; res%hist_denu(it) = denu
      res%ecomp = ecp; res%eterm = matmul(ecp, w)
      res%dosef = dosef(:, 1:nsp, :); res%qcell = qcell
      res%mix_nrestart = mixer%nrestart; res%mix_nh = mixer%nh
      res%depth_mean = 0.0_dp; if (ed_nchain > 0) res%depth_mean = real(ed_nlev, dp) / real(ed_nchain, dp)
      res%nfall = int(ed_nfall)
      if (vb) then
        print '(a,i3,a,es10.3,a,f10.6,a,*(f8.4))', 'it', it, '  dV=', dv, '  EF=', ef, '  qnet=', cr(:)%qnet
        if (nsp == 2) print '(a,*(f8.4))', '                                      m   =', cr(:)%moment
        print '(a,f18.8,a)', '      E_total =', etot, ' Ry/site'
        if (ed_nchain > 0) print '(a,f6.1,a,i0,a)', '      mean recursion depth used ', real(ed_nlev, dp) / real(ed_nchain, dp), ' of ', nlev, ' levels'
        if (ed_nfall > 0) print '(a,i0,a,i0,a)', '      bp tail fixed point failed (hull used) in ', ed_nfall, ' of ', ed_nchain, ' chains'
        if (ed3) then
          do p = 1, ncomp
            print '(a,i0,a,*(3f10.5,2x))', '      channel charge q0 + p q2, component ', p, ' (s p d per spin):', (qch(p, sg, :), sg=1, nsp)
            print '(a,i0,a,*(3f10.5,2x))', '      q0 q1 q2 of the d channel, component ', p, ':', (mom(p, sg, :, 2), sg=1, nsp)
          end do
          print '(a,es9.2,a,i0,a,i0,a,es9.2)', '      ed3: max |dq1| / sqrt(|q0| |q2|) ', asym_max, ';  block -Im G not positive at ', npsd, ' of ', nnode, &
            ' node chains, max b^2/(a c) ', psd_max
        else
          print '(a,i0,a,i0,a)', '      edc: -Im G not positive at ', npsd, ' of ', nnode, ' node values'
        end if
        print '(a,i0,a,i0)', '      Anderson history ', mixer%nh, '  restarts so far ', mixer%nrestart
        print '(a,i0,a,a,a,es9.2,a,es10.3,a,i0,a,i0,a,f8.4,a,es10.3)', '      E_F passes ', npass, ' final ', &
          trim(merge('bracket', 'dos    ', efmode == 1)), ' shift ', efshift, '  |r|_2=', r2n, '  max dV: comp ', pdv, ' spin ', sgdv, &
          ' r=', cr(max(pdv, 1))%mesh%r(max(kdv, 1)), '  max dE_nu=', denu
      end if
      if (len_trim(inp%pot_out) > 0) then                          ! every iteration: a stopped run restarts from here
        call timer_start(res%timer, tb_io); call write_potentials(); call timer_stop(res%timer, tb_io)
      end if
      res%hist_wall(it) = wall_seconds() - tw_it
      ! dV below the tolerance with E_nu held at its window is not a fixed point (the result would depend on the window margin)
      if (dv < inp%tol .and. it > 3) then; res%converged = (res%nclamp == 0); exit; end if
    end do
    if (res%nclamp > 0) print '(a,i0,a)', '  WARNING: E_nu held at its band window in the last iteration (', res%nclamp, &
      ' channels): not a fixed point, the result depends on the window margin'
    if (res%nnoncausal > 0) print '(a,i0,a,i0,a)', '  WARNING: Green function not causal at ', res%nnoncausal, ' of ', res%nnode, &
      ' contour node values in the last iteration (truncated energy-dependent chain; manual, theory/energy_dependent)'
    ! valence weight below the contour bottom at the final potential (the occupied moments would miss it)
    call all_params()
    call weight_below(res%wbelow)
    res%below_ok = abs(res%wbelow) < wbelow_tol .and. abs(res%wbelow1) < wbelow_tol
    if (vb .or. abs(res%wbelow) >= wbelow_tol) print '(a,es10.2,a)', '  valence weight within 3 Ry below the contour bottom, final potential:', &
      res%wbelow, ' electrons per site (must vanish)'
    if (abs(res%wbelow) >= wbelow_tol) print '(a)', '  WARNING: states below the contour bottom: the charges and the energy miss them'
    if (res%uform) res%udrive = ucur
    call xc_finalize(xf)
    do p = 1, ncomp
      do sg = 1, nsp
        do l = 0, 2
          call potpar_at_enu(cr(p)%mesh, cr(p)%v(:, sg), l, enu(l, sg, p), inp%spec(psp(p))%z, c, avw, cr(p)%pp(l, sg))
          cr(p)%q(0:2, l, sg) = mom(p, sg, 0:2, l)
        end do
        cr(p)%enu(:, sg) = enu(:, sg, p)
      end do
    end do
    if (inp%ndos > 0) call dos_final()
    if (inp%spectral%on) then
      call timer_start(res%timer, tb_spectral); call spectral_final(); call timer_stop(res%timer, tb_spectral)
    end if
    if (inp%epi_shells > 0) then
      call timer_start(res%timer, tb_epi); call epi_run(); call timer_stop(res%timer, tb_epi)
    end if
    res%ncomp = ncomp; res%nspin = nsp
    res%comp_cls = pcl; res%comp_spec = psp; res%comp_w = w; res%comp_x = xp
    res%comp = cr
    res%moment_site = sum(w * cr(:)%moment)
    res%qch = qch
    res%timer%wall_total = wall_seconds() - tw0; res%timer%cpu_total = cpu_seconds() - tc0
    do kap = 1, size(cop)
      call augcls_free(cop(kap)%op)                      ! the node tables are pointers (asr_augcls): release them here
    end do

  contains

    subroutine all_params()
      integer :: pp_, gg, ll
      call timer_start(res%timer, tb_params)
      do pp_ = 1, ncomp
        do gg = 1, nsp
          do ll = 0, 2
            call potpar_at_enu(cr(pp_)%mesh, cr(pp_)%v(:, gg), ll, enu(ll, gg, pp_), inp%spec(psp(pp_))%z, c, avw, ppl, &
                               y(:, :, 1:cr(pp_)%mesh%n, ll, pp_, gg))
            cr(pp_)%pp(ll, gg) = ppl
            pps(pp_, gg)%enu(ll) = ppl%enu; pps(pp_, gg)%c(ll) = ppl%c; pps(pp_, gg)%sqrdel(ll) = ppl%sqrdel
            pps(pp_, gg)%p(ll) = ppl%p; pps(pp_, gg)%gamma(ll) = ppl%gamma; pps(pp_, gg)%alpha_rep(ll) = ppl%gamma
          end do
        end do
        if (nsp == 1) pps(pp_, 2) = pps(pp_, 1)
      end do
      call timer_stop(res%timer, tb_params)
    end subroutine all_params

    !> occupied moments mo(component, spin, k, l) on the contour to efv, the DOS at E_F (node nq + 1), the channel charges qc
    !> (edc: M0; ed3: q0 + p^gamma q2, dosv is their DOS) and with ed3 the DOS dm(component, spin, k, l) of each q_k
    subroutine ed_moments(efv, mo, dosv, qc, dm)
      real(dp), intent(in) :: efv
      real(dp), intent(out) :: mo(:, :, 0:, 0:), dosv(:, :, 0:), qc(:, :, 0:), dm(:, :, 0:, 0:)
      complex(dp), allocatable :: zq(:), wq(:), gz(:, :, :, :), gb(:, :, :, :, :)
      real(dp) :: xm(0:2), xd(0:2), oa(0:2), a_, b_, c_, t_
      complex(dp) :: q0z
      type(potpar_t) :: qa
      integer :: j, o, gg, kk, ll, nq, pp_
      call contour(efv, zq, wq, nq)
      allocate (gz(nq + 1, nomax, ncomp, nsp))
      mo = 0.0_dp; dosv = 0.0_dp; dm = 0.0_dp
      if (.not. ed3) then
        call chain_greens(zq, gz)
        npsd = 0; nnode = 0
        do gg = 1, nsp
          do pp_ = 1, ncomp
            do o = 1, norb(pcl(pp_))
              ll = lorbk(o, pcl(pp_))
              ! causality at every node in the upper half plane: -Im G > 0 for a Herglotz function (manual/src/theory/energy_dependent.md, causality)
              do j = 1, nq + 1
                if (aimag(zq(j)) <= 0.0_dp) cycle
                nnode = nnode + 1
                if (aimag(gz(j, o, pp_, gg)) >= 0.0_dp) npsd = npsd + 1
              end do
              dosv(pp_, gg, ll) = dosv(pp_, gg, ll) - worbk(o, pcl(pp_)) * aimag(gz(nq + 1, o, pp_, gg)) / pi
              do j = 1, nq
                do kk = 0, 2
                  mo(pp_, gg, kk, ll) = mo(pp_, gg, kk, ll) - worbk(o, pcl(pp_)) * aimag(wq(j) * (zq(j) - enu(ll, gg, pp_))**kk &
                                                                                       * gz(j, o, pp_, gg)) / pi
                end do
              end do
            end do
          end do
        end do
        qc = mo(:, :, 0, :)
        return
      end if
      allocate (gb(nq + 1, nomax, ncomp, nsp, 0:3))
      call chain_greens(zq, gz, gb)
      asym_max = 0.0_dp; psd_max = 0.0_dp; npsd = 0; nnode = 0
      qc = 0.0_dp
      do gg = 1, nsp
        do pp_ = 1, ncomp
          qa = to_alpha(pps(pp_, gg), alpha); oa = o_alpha(qa, alpha)
          do o = 1, norb(pcl(pp_))
            ll = lorbk(o, pcl(pp_))
            xm = 0.0_dp
            do kk = 0, 2
              do j = 1, nq
                xm(kk) = xm(kk) - worbk(o, pcl(pp_)) * aimag(wq(j) * gb(j, o, pp_, gg, kk)) / pi
              end do
              xd(kk) = -worbk(o, pcl(pp_)) * aimag(gb(nq + 1, o, pp_, gg, kk)) / pi
            end do
            mo(pp_, gg, 0, ll) = mo(pp_, gg, 0, ll) + xm(0) + 2.0_dp * oa(ll) * xm(1) + oa(ll)**2 * xm(2)
            mo(pp_, gg, 1, ll) = mo(pp_, gg, 1, ll) + xm(1) + oa(ll) * xm(2)
            mo(pp_, gg, 2, ll) = mo(pp_, gg, 2, ll) + xm(2)
            dm(pp_, gg, 0, ll) = dm(pp_, gg, 0, ll) + xd(0) + 2.0_dp * oa(ll) * xd(1) + oa(ll)**2 * xd(2)
            dm(pp_, gg, 1, ll) = dm(pp_, gg, 1, ll) + xd(1) + oa(ll) * xd(2)
            dm(pp_, gg, 2, ll) = dm(pp_, gg, 2, ll) + xd(2)
            qc(pp_, gg, ll) = qc(pp_, gg, ll) + xm(0) + 2.0_dp * oa(ll) * xm(1) + qa%p(ll) * xm(2)
            dosv(pp_, gg, ll) = dosv(pp_, gg, ll) + xd(0) + 2.0_dp * oa(ll) * xd(1) + qa%p(ll) * xd(2)
            ! physics bound at every node: -Im of the block matrix (q0, q1; q1, q2) is positive semidefinite in the upper half plane
            do j = 1, nq + 1
              if (aimag(zq(j)) <= 0.0_dp) cycle
              q0z = gb(j, o, pp_, gg, 0) + 2.0_dp * oa(ll) * gb(j, o, pp_, gg, 1) + oa(ll)**2 * gb(j, o, pp_, gg, 2)
              a_ = -aimag(q0z); c_ = -aimag(gb(j, o, pp_, gg, 2))
              b_ = -aimag(gb(j, o, pp_, gg, 1) + oa(ll) * gb(j, o, pp_, gg, 2))
              t_ = abs(gb(j, o, pp_, gg, 3)) / max(sqrt(abs(q0z) * abs(gb(j, o, pp_, gg, 2))), tiny(1.0_dp))
              asym_max = max(asym_max, t_)
              nnode = nnode + 1
              if (a_ <= 0.0_dp .or. c_ <= 0.0_dp) then
                npsd = npsd + 1
              else
                psd_max = max(psd_max, b_**2 / (a_ * c_))
                if (b_**2 > a_ * c_) npsd = npsd + 1
              end if
            end do
          end do
        end do
      end do
    end subroutine ed_moments

    !> diagnostic ([asr] ef_scan = n): N(E) - N_elec at the start potential on 2n + 1 energies around ef_scan_center, split into
    !> the parts of the Fermi contour (vertical leg, main line, Fermi step, tail line, Matsubara poles), the component charges
    !> (valence electrons per atom) and the bp fallback count.  N(E) of an analytic G is smooth: a jump between neighbouring
    !> energies locates a non-analytic step in G(z) (terminator branch or eigenvalue selection) in that part of the contour.
    subroutine ef_scan_run()
      complex(dp), allocatable :: zq(:), wq(:), gz(:, :, :, :)
      real(dp) :: e0, ev, x, dsum, part(5), qc(ncomp)
      integer :: js, j, o, gg, pp_, nq, ip, bnd(4)
      integer(8) :: nf0
      if (inp%fermi_kt <= 0.0_dp) error stop 'class_scf: ef_scan needs fermi_kt > 0'
      e0 = ef; if (inp%ef_scan_center < huge(1.0_dp)) e0 = inp%ef_scan_center
      print '(a,i0,a,es9.2,a,f10.6,a,a)', '  ef_scan: ', 2 * inp%ef_scan + 1, ' energies, step ', inp%ef_scan_step, ' Ry, centre ', e0, &
        ', tail ', trim(inp%ed_tail)
      print '(a)', '  ef_scan    E (Ry)     N-nelec       DOS     vert         main         step         tail         poles     &
                   & | q per component | bp fallbacks'
      do js = -inp%ef_scan, inp%ef_scan
        ev = e0 + js * inp%ef_scan_step
        nf0 = ed_nfall
        call contour(ev, zq, wq, nq)
        bnd = [4, 4 + inp%ncont, 4 + inp%ncont + 16, 4 + inp%ncont + 24]   ! last node of each part (contour_nodes_fermi order)
        allocate (gz(nq + 1, nomax, ncomp, nsp))
        if (len_trim(inp%ef_scan_dump) > 0) then
          allocate (dmp_a(nlev, nq + 1, nomax, ncomp, nsp), dmp_b(nlev, nq + 1, nomax, ncomp, nsp), dmp_m(nq + 1, nomax, ncomp, nsp))
          dmp_a = (0.0_dp, 0.0_dp); dmp_b = (0.0_dp, 0.0_dp); dmp_m = 0
        end if
        call chain_greens(zq, gz)
        if (allocated(dmp_a)) then
          call dump_chains(js, ev, zq, wq, nq)
          deallocate (dmp_a, dmp_b, dmp_m)
        end if
        part = 0.0_dp; qc = 0.0_dp; dsum = 0.0_dp
        do gg = 1, nsp
          do pp_ = 1, ncomp
            do o = 1, norb(pcl(pp_))
              dsum = dsum - fs * w(pp_) * worbk(o, pcl(pp_)) * aimag(gz(nq + 1, o, pp_, gg)) / pi
              do j = 1, nq
                ip = count(j > bnd) + 1
                x = -fs * worbk(o, pcl(pp_)) * aimag(wq(j) * gz(j, o, pp_, gg)) / pi
                part(ip) = part(ip) + w(pp_) * x; qc(pp_) = qc(pp_) + x
              end do
            end do
          end do
        end do
        print '(a,f12.7,es14.6,f9.3,5es13.5,a,*(f11.7))', '  ef_scan ', ev, sum(part) - nelec, dsum, part, ' |', qc
        print '(a,i0)', '  ef_scan   bp fallbacks: ', ed_nfall - nf0
        deallocate (gz)
      end do
    end subroutine ef_scan_run

    !> ef_scan_dump: the contour (node, z, w) and every chain (a_k, b_k, k = 1..m) at one scan energy, appended to the dump file;
    !> header: nelec fs, component weights w(p), orbital weights worb(o)
    subroutine dump_chains(js, ev, zq, wq, nq)
      integer, intent(in) :: js, nq
      real(dp), intent(in) :: ev
      complex(dp), intent(in) :: zq(:), wq(:)
      integer :: u, j, o, pp_, gg, k
      if (.not. par_root()) return
      if (js == -inp%ef_scan) then
        if (any(norb /= norb(1))) error stop 'class_scf: ef_scan_dump needs the same start orbitals on every sublattice'
        do k = 2, ncls
          if (any(worbk(:, k) /= worbk(:, 1))) error stop 'class_scf: ef_scan_dump needs the same start orbitals on every sublattice'
        end do
        open (newunit=u, file=trim(inp%ef_scan_dump), status='replace', action='write')
        write (u, '(a,es25.16,f6.1,3i6)') 'H ', nelec, fs, ncomp, nsp, nlev
        write (u, '(a,*(es25.16))') 'W ', w
        write (u, '(a,*(f6.1))') 'O ', worbk(1:norb(1), 1)
      else
        open (newunit=u, file=trim(inp%ef_scan_dump), status='old', position='append', action='write')
      end if
      write (u, '(a,i4,es25.16,i6)') 'E ', js, ev, nq
      do j = 1, nq + 1
        write (u, '(a,i5,4es25.16)') 'Z ', j, zq(j), wq(j)
      end do
      do gg = 1, nsp; do pp_ = 1, ncomp; do o = 1, norb(pcl(pp_)); do j = 1, nq + 1
        write (u, '(a,4i5)') 'C ', j, o, pp_, gg
        write (u, '(i4)') dmp_m(j, o, pp_, gg)
        do k = 1, dmp_m(j, o, pp_, gg)
          write (u, '(4es25.16)') dmp_a(k, j, o, pp_, gg), dmp_b(k, j, o, pp_, gg)
        end do
      end do; end do; end do; end do
      close (u)
    end subroutine dump_chains

    !> G(z) of every chain gz(energy, start orbital, component, spin) at the energies zs, the chains in parallel.  ed3: gz is
    !> G + 2 o hG + p hGh, and gb(.., 0:3) gets the blocks of ed3_blocks.  With MPI each rank computes its share of the chain
    !> tasks (chain_worker) and holds zeros elsewhere; one sum over the ranks then gives every rank the full arrays.  Each entry
    !> has one non-zero term, so the sum is exact and the same on every rank in any reduction order.
    subroutine chain_greens(zs, gz, gb, tbk)
      complex(dp), intent(in) :: zs(:)
      complex(dp), intent(out) :: gz(:, :, :, :)
      complex(dp), intent(out), optional :: gb(:, :, :, :, 0:)
      integer, intent(in), optional :: tbk           !! timer bucket (asr_runstats); default tb_chains
      type(class_alloy_t), allocatable :: ms(:, :)
      type(edn_t), allocatable :: es(:, :)
      type(ed3_t), allocatable :: e3(:, :)
      integer :: gg, pp_, nbe, kp, ib, nth, nout, nin, nmy, ntg_, t_, lvl
      real(dp) :: tw0
      complex(dp), allocatable :: trk(:)
      integer(8) :: cnt0(3), cnt(3)
      ib = tb_chains; if (present(tbk)) ib = tbk
      call timer_start(res%timer, ib)
      allocate (ms(ncomp, nsp), es(ncomp, nsp), e3(ncomp, nsp))
      do gg = 1, nsp
        do pp_ = 1, ncomp
          if (ed3) then
            call ed3_setup(mdl, pps(:, gg), alpha, pp_, e3(pp_, gg))
          else
            call edc_setup(mdl, pps(:, gg), alpha, pp_, es(pp_, gg))
          end if
          ms(pp_, gg) = mdl
        end do
      end do
      ! block size: energies per chain block, the same on every thread, capped so that one thread's block work space (10 vectors
      ! per member) stays below 400 MB on the largest class operator; the ef_scan dump keeps one energy per chain
      nbe = 0
      if (inp%eblock > 0 .and. .not. allocated(dmp_a)) then
        nbe = min(inp%eblock, size(zs))
        do kp = 1, ncls
          nbe = min(nbe, blk_nbmax(inp%eblock, cop(kp)%op%dim(), 10, 400.0_dp))
        end do
      end if
      if (par_size() > 1) then
        gz = (0.0_dp, 0.0_dp)
        if (present(gb)) gb = (0.0_dp, 0.0_dp)
        cnt0 = [ed_nchain, ed_nlev, ed_nfall]
      end if
      ! threads: when this rank holds fewer chains than threads, teams of nin threads share the matrix-vector products of one
      ! chain (asr_chain::chain_teams; nested OpenMP, the inner loops of asr_augcls).  Otherwise one thread per chain, as before.
      ! The inner loops write disjoint parts of the vectors, so the results do not depend on the split.
      ntg_ = size(zs); if (nbe > 0) ntg_ = (size(zs) + nbe - 1) / nbe
      nmy = 0
      do t_ = 1, ntg_ * nomax * ncomp * nsp
        if (mod((t_ - 1) / ntg_, nomax) + 1 > norb(pcl(mod((t_ - 1) / (ntg_ * nomax), ncomp) + 1))) cycle   ! no such start orbital
        if (mod(t_ - 1 + (t_ - 1) / ntg_, par_size()) == par_rank()) nmy = nmy + 1
      end do
      nth = 1
      !$ nth = omp_get_max_threads()
      call chain_teams(nmy, nth, nout, nin)
      call forced_teams(nout, nin)                           ! diagnostic: ASR_CHAIN_TEAMS=nout,nin
      lvl = 1
      !$ lvl = omp_get_max_active_levels()
      !$ if (nin > 1) call omp_set_max_active_levels(max(lvl, 2))
      if (vb .and. (teams_said(1) /= nout .or. teams_said(2) /= nin)) then
        if (nin > 1) then
          print '(a,i0,a,i0,a,i0,a,i0,a)', '  chains: ', nmy, ' per rank on ', nth, ' threads: ', nout, ' teams of ', nin, ' threads'
        else if (teams_said(1) /= 0) then
          print '(a,i0,a,i0,a)', '  chains: ', nmy, ' per rank on ', nth, ' threads: one thread per chain'
        end if
        teams_said = [nout, nin]
      end if
      tw0 = wall_seconds()
      chain_vec_par = nin > 1                                 ! the element-wise chain loops: parallel only inside a team
      !$omp parallel num_threads(nout)
      !$ if (nin > 1) call omp_set_num_threads(nin)
      call chain_worker(zs, ms, es, e3, nbe, gz, gb)
      !$omp end parallel
      chain_vec_par = .true.
      !$ if (nin > 1) call omp_set_max_active_levels(lvl)
      ! each rank's own chain time (before the sum over ranks waits for the slowest): the load balance between ranks
      allocate (trk(par_size())); trk = (0.0_dp, 0.0_dp); trk(par_rank() + 1) = cmplx(wall_seconds() - tw0, 0.0_dp, dp)
      if (par_size() > 1) call par_sum_z(trk, size(trk))
      if (.not. allocated(res%chain_rank_wall)) then; allocate (res%chain_rank_wall(par_size())); res%chain_rank_wall = 0.0_dp; end if
      res%chain_rank_wall = res%chain_rank_wall + real(trk)
      res%chain_split = [nout, nin]
      if (par_size() > 1) then
        call par_sum_z(gz, size(gz))
        if (present(gb)) call par_sum_z(gb, size(gb))
        if (allocated(dmp_a)) then
          call par_sum_z(dmp_a, size(dmp_a)); call par_sum_z(dmp_b, size(dmp_b)); call par_sum_i(dmp_m, size(dmp_m))
        end if
        cnt = [ed_nchain, ed_nlev, ed_nfall] - cnt0
        call par_sum_i8(cnt, 3)
        ed_nchain = cnt0(1) + cnt(1); ed_nlev = cnt0(2) + cnt(2); ed_nfall = cnt0(3) + cnt(3)
      end if
      call timer_stop(res%timer, ib, int(size(zs), i8) * int(sum(norb(pcl)) * nsp, i8))
    end subroutine chain_greens

    !> Diagnostic: the environment variable ASR_CHAIN_TEAMS = "nout,nin" forces the thread split of chain_greens (timing of the
    !> nested teams on one machine, where the contour always has more chains than threads).  The results do not depend on it.
    subroutine forced_teams(nout, nin)
      integer, intent(inout) :: nout, nin
      character(len=32) :: v
      integer :: ln, st_, ios, a, b
      call get_environment_variable('ASR_CHAIN_TEAMS', v, ln, st_)
      if (st_ /= 0 .or. ln == 0) return
      read (v(1:ln), *, iostat=ios) a, b
      if (ios /= 0 .or. a < 1 .or. b < 1) error stop 'ASR_CHAIN_TEAMS must be "nout,nin" with two positive integers'
      nout = a; nin = b
    end subroutine forced_teams

    !> one thread's share of the chains (called inside the parallel region; locals are private): a private copy of the class
    !> operator per sublattice with its own work space, re-parametrised / re-frozen / re-started only when the chain needs it
    subroutine chain_worker(zs, ms, es, e3, nbe, gz, gb)
      complex(dp), intent(in) :: zs(:)
      type(class_alloy_t), intent(in) :: ms(:, :)
      type(edn_t), intent(in) :: es(:, :)
      type(ed3_t), intent(in) :: e3(:, :)
      integer, intent(in) :: nbe                     !! energies per block; 0: one energy per chain (real kernels)
      complex(dp), intent(inout) :: gz(:, :, :, :)
      complex(dp), intent(inout), optional :: gb(:, :, :, :, 0:)
      type(augcls_op), allocatable :: lop(:)
      type(edn_work_t), allocatable :: lwk(:)
      type(ed3_work_t), allocatable :: lw3(:)
      type(blk_work_t), allocatable :: lbw(:)
      logical :: newstart(ncls)
      complex(dp) :: bk(0:3)
      complex(dp), allocatable :: bkb(:, :)
      integer :: cur(3, ncls), t, nz, j, o, gg, pp_, kp, ks, nbk, jb, j1, j2, iv, ntg, rk, nr
      nz = size(zs)
      nbk = 1; if (nbe > 0) nbk = (nz + nbe - 1) / nbe        ! blocks per (start orbital, component, spin)
      ntg = merge(nz, nbk, nbe == 0)                          ! tasks per (start orbital, component, spin)
      rk = par_rank(); nr = par_size()
      allocate (lop(ncls), lwk(ncls), lw3(ncls), lbw(ncls)); cur = -1; newstart = .true.
      if (nbe > 0) allocate (bkb(0:3, nbe))
      !$omp do schedule(dynamic)
      do t = 1, ntg * nomax * ncomp * nsp
        ! MPI: this rank's share, cyclic in t; the shift by the group index rotates the energies over the ranks, so that no
        ! rank keeps the same energies (whose cost differs) in every group.  One rank: every task.
        if (mod(t - 1 + (t - 1) / ntg, nr) /= rk) cycle
        if (nbe == 0) then
          j = mod(t - 1, nz) + 1; o = mod((t - 1) / nz, nomax) + 1
          pp_ = mod((t - 1) / (nomax * nz), ncomp) + 1; gg = (t - 1) / (nomax * nz * ncomp) + 1
          j1 = j; j2 = j
        else
          jb = mod(t - 1, nbk) + 1; o = mod((t - 1) / nbk, nomax) + 1
          pp_ = mod((t - 1) / (nomax * nbk), ncomp) + 1; gg = (t - 1) / (nomax * nbk * ncomp) + 1
          j1 = (jb - 1) * nbe + 1; j2 = min(nz, jb * nbe)
        end if
        kp = pcl(pp_)
        if (o > norb(kp)) cycle                                 ! this sublattice has fewer start orbitals
        if (cur(1, kp) < 0) then
          lop(kp) = cop(kp)%op
          nullify (lop(kp)%wk); allocate (lop(kp)%wk(size(cop(kp)%op%wk, 1), size(cop(kp)%op%wk, 2)))
          if (ed3) then
            call ed3_work_init(lw3(kp), lop(kp)%dim())
          else if (nbe == 0) then
            call edn_work_init(lwk(kp), lop(kp)%dim())
          end if
          if (nbe > 0) call blk_work_init(lbw(kp), nbe, lop(kp)%dim())
          cur(:, kp) = 0
        end if
        if (cur(1, kp) /= pp_ .or. cur(2, kp) /= gg) then
          call augcls_set_params(lop(kp), ms(pp_, gg)); call augcls_freeze_site(lop(kp), 1, pp_)
          if (inp%sro) then
            do ks = 1, ncls
              call augcls_set_class_conc(lop(kp), ncls + ks, inp%sro_cond(pk(pp_), 1:inp%nsc(ks), kp, ks))
            end do
          end if
          cur(1, kp) = pp_; cur(2, kp) = gg; newstart(kp) = .true.
        end if
        if (cur(3, kp) /= o) then
          call augcls_set_start(lop(kp), iorbk(o, kp), cop(kp)%d); cur(3, kp) = o; newstart(kp) = .true.
        end if
        iv = iorbk(o, kp)
        if (ed3) then
          if (newstart(kp)) then; call ed3_start(lop(kp), lw3(kp), e3(pp_, gg), iv); newstart(kp) = .false.; end if
          if (nbe == 0) then
            call ed3_blocks(lop(kp), lw3(kp), e3(pp_, gg), zs(j1), iv, nlev, bk)
            gz(j1, o, pp_, gg) = bk(0) + 2.0_dp * e3(pp_, gg)%o(iv, pp_) * bk(1) + e3(pp_, gg)%p(iv, pp_) * bk(2)
            if (present(gb)) gb(j1, o, pp_, gg, :) = bk
          else
            call ed3_blocks_blk(lop(kp), lw3(kp), lbw(kp), e3(pp_, gg), zs(j1:j2), iv, nlev, bkb(:, 1:j2 - j1 + 1))
            gz(j1:j2, o, pp_, gg) = bkb(0, 1:j2 - j1 + 1) + 2.0_dp * e3(pp_, gg)%o(iv, pp_) * bkb(1, 1:j2 - j1 + 1) &
                                    + e3(pp_, gg)%p(iv, pp_) * bkb(2, 1:j2 - j1 + 1)
            if (present(gb)) then
              do j = j1, j2; gb(j, o, pp_, gg, :) = bkb(:, j - j1 + 1); end do
            end if
          end if
        else if (allocated(dmp_a)) then
          call edc_green(lop(kp), lwk(kp), es(pp_, gg), zs(j1), iv, nlev, gz(j1, o, pp_, gg), &
                         acoef=dmp_a(:, j1, o, pp_, gg), bcoef=dmp_b(:, j1, o, pp_, gg), mused=dmp_m(j1, o, pp_, gg))
        else if (nbe == 0) then
          call edc_green(lop(kp), lwk(kp), es(pp_, gg), zs(j1), iv, nlev, gz(j1, o, pp_, gg))
        else
          call edc_green_blk(lop(kp), lbw(kp), es(pp_, gg), zs(j1:j2), iv, nlev, gz(j1:j2, o, pp_, gg))
        end if
      end do
      !$omp end do
      do kp = 1, ncls
        if (cur(1, kp) >= 0) deallocate (lop(kp)%wk)
      end do
    end subroutine chain_worker

    !> Effective pair interactions of the final potentials (asr_epi) on the first inp%epi_shells shells of every sublattice.
    !> One operator per pair type: origin on sublattice kp, the partner site frozen, every other site averaged with the
    !> concentrations of its sublattice (no SRO shell: the pair table is that of the uncorrelated medium around the pair).
    !> The chains (path species x path point x partner species x spin x start vector x energy) run in parallel as in
    !> chain_greens.  A pair type between two sublattices is computed once, from the lower sublattice.
    subroutine epi_run()
      type(epi_pair_t), allocatable :: prs(:)
      type(epi_table_t), allocatable :: tab(:)
      type(class_alloy_t) :: mp
      type(augcls_op) :: pop
      type(epi_par_t), allocatable :: par(:, :)
      integer, allocatable :: cl0(:), clp(:), pidx(:)
      logical, allocatable :: flip(:)
      real(dp), allocatable :: pd(:, :, :), chi(:, :), xl(:), wl(:), db(:), rvs(:, :), gr(:), wr(:)
      complex(dp), allocatable :: zq(:), wq(:), tz(:, :, :, :, :, :), zb(:), wb(:), tb(:, :, :, :, :, :)
      real(dp) :: vec(9, 9), tailmax, ebot, dls, dl0, acc, dtol
      integer :: lo(9), kp, k2, ip, nq, n1, n2, nlam, ntab, c, s, il, gg, ia, j, jsite, nfr, nbe
      character(len=:), allocatable :: emsg
      integer, parameter :: nba = 16, nbv = 4
      nlam = inp%epi_nlam
      dtol = ed_depth_tol; ed_depth_tol = 0.0_dp      ! full depth: U is a small difference of large numbers (route ed3 adapts)
      allocate (par(ncomp, nsp), xl(nlam), wl(nlam), db(ncomp), tab(0), zb(nba + nbv), wb(nba + nbv))
      do gg = 1, nsp
        do c = 1, ncomp; par(c, gg) = epi_species(pps(c, gg), alpha); end do
      end do
      call gauss_legendre(nlam, xl, wl)
      xl = 0.5_dp * (1.0_dp + xl); wl = 0.5_dp * wl
      call contour(ef, zq, wq, nq)
      ebot = minval(enu(:, 1:nsp, :)) - inp%ebot_offset
      call contour_nodes_closed(ebot - 3.0_dp, ebot, 0.05_dp, nba, nbv, zb, wb)
      ! drive deviations of the charge-excess model: the SRO drive itself, or dq / G00 in the random alloy
      if (res%uform) then
        db = ucur
      else
        do c = 1, ncomp; db(c) = (cr(c)%qnet - qbar(pcl(c))) / res%screen%g00; end do
      end if
      res%epi_wbelow = 0.0_dp; ntab = 0
      if (vb) print '(a,i0,a,i0,a,a,a,i0)', '  effective pair interactions: ', inp%epi_shells, ' shells, ', nlam, ' path points', &
        trim(merge(', pair-exchange symmetry used', '                             ', inp%epi_sym)), ', energy_block ', inp%eblock
      do kp = 1, ncls
        cl0 = sgr(kp)%cls
        call epi_pair_list_c(sgr(kp)%cg%disp, sgr(kp)%perm, sgr(kp)%rot, cl0, inp%epi_shells, prs, emsg)
        if (len(emsg) > 0) error stop emsg
        if (size(prs) == 0) cycle
        allocate (rvs(3, size(prs)), gr(size(prs)), wr(size(prs)), pidx(size(prs)))
        do ip = 1, size(prs)
          rvs(:, ip) = inp%alat * prs(ip)%rvec
          pidx(ip) = sgr(kp)%cg%basis(prs(ip)%site)
        end do
        call screen_structure_shells(st, inp%alat, res%screen, rvs, gr, wr, b0=sgr(kp)%cg%b0, b2=pidx)
        deallocate (pidx)
        do ip = 1, size(prs)
          k2 = prs(ip)%cls2; jsite = prs(ip)%site
          n1 = inp%nsc(kp); n2 = inp%nsc(k2)
          if (k2 < kp .or. n1 < 2 .or. n2 < 2) cycle        ! the other sublattice computes it; one species: nothing to exchange
          ! pair operator: the sublattice classes, the partner in a class of its own, one more species (the path site)
          clp = cl0; clp(jsite) = ncls + 1
          mp%nl = 9; mp%nsp = ncomp + 1; mp%ncls = ncls + 1
          if (allocated(mp%nsc)) deallocate (mp%nsc, mp%spc, mp%conc)
          allocate (mp%nsc(ncls + 1), mp%spc(size(mdl%spc, 1), ncls + 1), mp%conc(size(mdl%conc, 1), ncls + 1))
          mp%nsc(1:ncls) = mdl%nsc(1:ncls); mp%spc(:, 1:ncls) = mdl%spc(:, 1:ncls); mp%conc(:, 1:ncls) = mdl%conc(:, 1:ncls)
          mp%nsc(ncls + 1) = mdl%nsc(k2); mp%spc(:, ncls + 1) = mdl%spc(:, k2); mp%conc(:, ncls + 1) = mdl%conc(:, k2)
          mp%block = mdl%block; mp%block0 = mdl%block0
          flip = sgr(kp)%flip; flip(jsite) = .false.
          call epi_pair_ops_c(sgr(kp)%perm, sgr(kp)%rot, clp, pidx)
          if (allocated(pd)) deallocate (pd)
          allocate (pd(9, 9, size(pidx)))
          do j = 1, size(pidx); call orbital_matrix(sgr(kp)%rot(:, :, pidx(j)), pd(:, :, j)); end do
          call build_augcls_operator(sgr(kp)%cg%g, mp, clp, sgr(kp)%perm(:, pidx), pd, 1, ngen, pop, inp%maxrank, flip, &
                                     stype=sgr(kp)%stype)
          if (allocated(chi)) deallocate (chi)
          allocate (chi(size(pidx), 9))
          call epi_adapted(pd, vec, chi, lo)
          call augcls_freeze_any(pop, jsite, mdl%spc(1, k2), nfr)
          if (nfr == 0) print '(a,i0,a)', '  WARNING: effective pair interactions, shell ', prs(ip)%shell, &
            ': the partner site is outside the node set (increase [asr] levels); this pair table is zero'
          allocate (tz(nq, 9, nsp, n2, nlam, n1), tb(nba + nbv, 9, nsp, 1, nlam, n1))
          tz = (0.0_dp, 0.0_dp); tb = (0.0_dp, 0.0_dp); tailmax = 0.0_dp
          ! energies per chain block ([asr] energy_block, as in chain_greens): the same on every thread, one thread's block
          ! work space (10 vectors per member) below 400 MB
          nbe = 0
          if (inp%eblock > 0) nbe = blk_nbmax(inp%eblock, pop%dim(), 10, 400.0_dp)
          chain_vec_par = .false.                             ! one thread per chain (asr_chain::chain_vec_par)
          !$omp parallel
          call epi_worker(pop, mp, par, kp, k2, jsite, vec, chi, lo, zq(1:nq), xl, n2, .false., min(nbe, nq), tz, tailmax)
          call epi_worker(pop, mp, par, kp, k2, jsite, vec, chi, lo, zb, xl, 1, .true., min(nbe, nba + nbv), tb, tailmax)
          !$omp end parallel
          chain_vec_par = .true.
          if (par_size() > 1) then
            call par_sum_z(tz, size(tz)); call par_sum_z(tb, size(tb))
          end if
          ntab = ntab + 1
          tab = [tab, epi_table_t()]
          associate (t => tab(ntab))
            t%kap = kp; t%kap2 = k2; t%shell = prs(ip)%shell; t%mult = prs(ip)%mult; t%nop = size(pidx); t%nodes = pop%nn
            t%vec = prs(ip)%vec; t%orbit = prs(ip)%orbit; t%rvec = prs(ip)%rvec; t%rorbit = prs(ip)%rorbit
            t%dist = norm2(rvs(:, ip)); t%wr = wr(ip); t%tail = tailmax
            allocate (t%uband(n1, n2), t%ues(n1, n2)); t%uband = 0.0_dp; t%ues = 0.0_dp
            do c = 2, n1
              do il = 1, nlam
                dl0 = 0.0_dp
                do s = 1, n2
                  if (epi_skip(inp%epi_sym, kp == k2, c, s)) cycle
                  dls = 0.0_dp
                  do gg = 1, nsp; do ia = 1, 9; do j = 1, nq
                    dls = dls - fs * aimag(wq(j) * tz(j, ia, gg, s, il, c)) / pi
                  end do; end do; end do
                  if (s == 1) dl0 = dls
                  t%uband(c, s) = t%uband(c, s) + wl(il) * (dls - dl0)
                end do
                acc = 0.0_dp                                 ! weight of the path site below the contour bottom
                do gg = 1, nsp; do ia = 1, 9; do j = 1, nba + nbv
                  acc = acc - fs * aimag(wb(j) * tb(j, ia, gg, 1, il, c)) / pi
                end do; end do; end do
                res%epi_wbelow = max(res%epi_wbelow, abs(acc))
              end do
              do s = 2, n2
                t%ues(c, s) = wr(ip) * (db(mdl%spc(c, kp)) - db(mdl%spc(1, kp))) * (db(mdl%spc(s, k2)) - db(mdl%spc(1, k2)))
              end do
            end do
            if (inp%epi_sym .and. kp == k2) call epi_mirror(t%uband)
            if (vb) then
              if (maxval(abs(2.0_dp * t%rvec - real(t%vec, dp))) < 1.0e-6_dp) then
                print '(a,i0,a,i0,a,i0,a,3i3,a,f8.4,a,i0,a,i0,a,i0,a,i0,a,es9.2)', &
                  '  pair: sublattices ', kp, ' - ', k2, '  shell ', t%shell, '  vector (a/2)', t%vec, '  R=', t%dist, &
                  ' Bohr  sites ', t%mult, '  group order ', t%nop, '  nodes ', t%nodes, &
                  '  energies per block ', max(1, min(nbe, nq)), '  X1 mismatch ', t%tail
              else
                print '(a,i0,a,i0,a,i0,a,3f7.3,a,f8.4,a,i0,a,i0,a,i0,a,i0,a,es9.2)', &
                  '  pair: sublattices ', kp, ' - ', k2, '  shell ', t%shell, '  vector (a)', t%rvec, '  R=', t%dist, &
                  ' Bohr  sites ', t%mult, '  group order ', t%nop, '  nodes ', t%nodes, &
                  '  energies per block ', max(1, min(nbe, nq)), '  X1 mismatch ', t%tail
              end if
              do c = 2, n1
                print '(a,i0,a,*(f10.4))', '    U(', c, ', s) band (mRy), s = 2..:', 1.0e3_dp * t%uband(c, 2:n2)
                print '(a,i0,a,*(f10.4))', '    U(', c, ', s) screened Coulomb (mRy):', 1.0e3_dp * t%ues(c, 2:n2)
              end do
            end if
          end associate
          deallocate (tz, tb); call augcls_free(pop)              ! the node tables are pointers (asr_augcls)
        end do
        deallocate (rvs, gr, wr)
      end do
      res%epi = tab(1:ntab)
      ed_depth_tol = dtol
      if (vb .or. res%epi_wbelow >= wbelow_tol) print '(a,es10.2,a)', &
        '  effective pair interactions: largest weight of a path site below the contour bottom:', res%epi_wbelow, ' electrons (must vanish)'
      if (res%epi_wbelow >= wbelow_tol) print '(a)', '  WARNING: a path species has states below the contour bottom: the pair tables miss them'
    end subroutine epi_run

    !> one thread's share of the chains of one pair type (called inside the parallel region; locals are private).  Tasks:
    !> (path species c = 2..n1) x (path point) x (partner species 1..ns2) x spin x (start vector) x (energy zs, or block of nbe
    !> energies; nbe = 0: one energy per chain, real kernels).  tz gets the integrand of asr_epi, or with g0 the block X0 alone
    !> (weight below the contour bottom).  tl: largest X1 mismatch.
    subroutine epi_worker(pop, mp, par, kp, k2, jsite, vec, chi, lo, zs, xl, ns2, g0, nbe, tz, tl)
      type(augcls_op), intent(in) :: pop
      type(class_alloy_t), intent(in) :: mp
      type(epi_par_t), intent(in) :: par(:, :)
      integer, intent(in) :: kp, k2, jsite, lo(9), ns2, nbe
      real(dp), intent(in) :: vec(9, 9), chi(:, :), xl(:)
      complex(dp), intent(in) :: zs(:)
      logical, intent(in) :: g0
      complex(dp), intent(inout) :: tz(:, :, :, :, :, :)
      real(dp), intent(inout) :: tl
      type(augcls_op) :: lop
      type(epi_work_t) :: wk
      type(ed3_work_t) :: w3
      type(blk_work_t) :: bw
      type(class_alloy_t) :: ml
      type(edn_t) :: edl
      type(ed3_t) :: e3l
      type(epi_path_t) :: pt
      complex(dp) :: xb(0:3, max(nbe, 1))
      real(dp) :: tloc
      integer :: cur(4), cura, t, r, j, jb, j1, j2, ia, gg, s, il, c, nz, nl, n1, rk, nr, ntg
      nz = size(zs); nl = size(xl); n1 = inp%nsc(kp)
      ntg = nz; if (nbe > 0) ntg = (nz + nbe - 1) / nbe          ! tasks per (start vector, spin, partner, path point, path species)
      rk = par_rank(); nr = par_size()
      lop = pop
      nullify (lop%wk); allocate (lop%wk(size(pop%wk, 1), size(pop%wk, 2)))
      if (ed3) then
        call ed3_work_init(w3, lop%dim())
      else
        call epi_work_init(wk, lop%dim())
      end if
      if (nbe > 0) call blk_work_init(bw, nbe, lop%dim())
      cur = -1; cura = -1; tloc = 0.0_dp
      !$omp do schedule(dynamic)
      do t = 1, ntg * 9 * nsp * ns2 * nl * (n1 - 1)
        if (mod(t - 1 + (t - 1) / ntg, nr) /= rk) cycle
        jb = mod(t - 1, ntg) + 1; r = (t - 1) / ntg
        ia = mod(r, 9) + 1; r = r / 9
        gg = mod(r, nsp) + 1; r = r / nsp
        s = mod(r, ns2) + 1; r = r / ns2
        il = mod(r, nl) + 1; c = r / nl + 2
        if (.not. g0 .and. epi_skip(inp%epi_sym, kp == k2, c, s)) cycle
        if (any(cur /= [c, il, s, gg])) then
          pt = epi_path(par(mp%spc(1, kp), gg), par(mp%spc(c, kp), gg), xl(il))
          ml = mp
          if (ed3) then
            call epi3_setup(ml, par(:, gg), pt, e3l)
          else
            call epi_setup(ml, par(:, gg), pt, edl)
          end if
          call augcls_set_params(lop, ml)
          call augcls_freeze_any(lop, 1, ncomp + 1); call augcls_freeze_any(lop, jsite, mp%spc(s, k2))
          cur = [c, il, s, gg]; cura = -1
        end if
        if (cura /= ia) then
          call augcls_set_chi(lop, chi(1:lop%nop, ia))
          if (ed3) then
            call epi3_start(lop, w3, e3l, vec(:, ia), lo(ia))
          else
            call epi_start(lop, wk, vec(:, ia))
          end if
          cura = ia
        end if
        if (nbe == 0) then
          j1 = jb; j2 = jb
          if (ed3) then
            call ed3_blocks(lop, w3, e3l, zs(j1), lo(ia), nlev, xb(:, 1))
          else
            call epi_blocks(lop, wk, edl, zs(j1), nlev, xb(:, 1))
          end if
        else
          j1 = (jb - 1) * nbe + 1; j2 = min(nz, jb * nbe)
          if (ed3) then
            call ed3_blocks_blk(lop, w3, bw, e3l, zs(j1:j2), lo(ia), nlev, xb(:, 1:j2 - j1 + 1))
          else
            call epi_blocks_blk(lop, wk, bw, edl, zs(j1:j2), nlev, xb(:, 1:j2 - j1 + 1))
          end if
        end if
        do j = j1, j2
          if (g0) then
            tz(j, ia, gg, s, il, c) = xb(0, j - j1 + 1)
          else
            if (ed3) then
              tz(j, ia, gg, s, il, c) = epi3_integrand(pt, lo(ia), zs(j), xb(:, j - j1 + 1))
            else
              tz(j, ia, gg, s, il, c) = epi_integrand(pt, lo(ia), zs(j), xb(:, j - j1 + 1))
            end if
            tloc = max(tloc, abs(xb(3, j - j1 + 1)) / max(abs(xb(1, j - j1 + 1)), tiny(1.0_dp)))
          end if
        end do
      end do
      !$omp end do
      !$omp critical (epi_tail)
      tl = max(tl, tloc)
      !$omp end critical (epi_tail)
      deallocate (lop%wk)
    end subroutine epi_worker

    subroutine contour(efv, zq, wq, nq)
      real(dp), intent(in) :: efv
      complex(dp), allocatable, intent(out) :: zq(:), wq(:)
      integer, intent(out) :: nq
      real(dp) :: ebot
      integer :: npl
      ebot = minval(enu(:, 1:nsp, :)) - inp%ebot_offset
      if (inp%fermi_kt > 0.0_dp) then
        npl = fermi_npole(inp%fermi_kt); if (inp%fermi_poles > 0) npl = inp%fermi_poles
        nq = fermi_nodes_count(npl, inp%ncont, 16)
        allocate (zq(nq + 1), wq(nq + 1))
        call contour_nodes_fermi(ebot, efv, inp%fermi_kt, npl, inp%ncont, 16, zq(1:nq), wq(1:nq))
      else
        nq = inp%ncont + nvert
        allocate (zq(nq + 1), wq(nq + 1))
        call contour_nodes_closed(ebot, efv, inp%ed_eta, inp%ncont, nvert, zq(1:nq), wq(1:nq))
      end if
      zq(nq + 1) = cmplx(efv, inp%ed_eta, dp); wq(nq + 1) = (0.0_dp, 0.0_dp)
    end subroutine contour

    !> wb = valence electrons per site with energies in [ebot - 3 Ry, ebot], ebot = the bottom of the occupation contour: the
    !> integral of G over a closed path around that interval (arc through the upper half plane, 16 + 4 Gauss-Legendre nodes).
    !> G is analytic there when no state lies below ebot; then wb vanishes to quadrature accuracy.  A finite value is the weight
    !> that the occupied moments miss.  Cost: 20 energies, about a third of one E_F pass.
    subroutine weight_below(wb)
      real(dp), intent(out) :: wb
      integer, parameter :: na = 16, nv = 4
      complex(dp) :: zq(na + nv), wq(na + nv)
      complex(dp), allocatable :: gz(:, :, :, :)
      real(dp) :: ebot
      integer :: j, o, gg, pp_
      ebot = minval(enu(:, 1:nsp, :)) - inp%ebot_offset
      call contour_nodes_closed(ebot - 3.0_dp, ebot, 0.05_dp, na, nv, zq, wq)
      allocate (gz(na + nv, nomax, ncomp, nsp))
      call chain_greens(zq, gz, tbk=tb_below)
      wb = 0.0_dp
      do gg = 1, nsp; do pp_ = 1, ncomp; do o = 1, norb(pcl(pp_)); do j = 1, na + nv
        wb = wb - fs * w(pp_) * worbk(o, pcl(pp_)) * aimag(wq(j) * gz(j, o, pp_, gg)) / pi
      end do; end do; end do; end do
    end subroutine weight_below

    !> SRO: the first shell of the start site of the cluster of sublattice kp goes to the classes ncls + kappa'.
    !> If SRO is given for specific sublattice pairs (kp, k), the nearest shell of sublattice k around kp is selected; otherwise
    !> the nearest shell of any sublattice around kp is selected.
    subroutine shell_classes(kp, cl)
      integer, intent(in) :: kp
      integer, intent(inout) :: cl(:)
      call shell_classes_dsp(kp, sgr(kp)%cg%disp, sgr(kp)%cls, cl)
    end subroutine shell_classes

    subroutine shell_classes_dsp(kp, dsp, scls, cl)
      integer, intent(in) :: kp, scls(:)
      real(dp), intent(in) :: dsp(:, :)
      integer, intent(inout) :: cl(:)
      integer :: i_, k_
      real(dp) :: r2m
      logical :: has_k(ncls)
      do k_ = 1, ncls
        has_k(k_) = any(abs(inp%sro_cond(1:inp%nsc(kp), 1:inp%nsc(k_), kp, k_) - &
                            spread(inp%conc(1:inp%nsc(k_), k_), 1, inp%nsc(kp))) > 1.0e-12_dp)
      end do
      if (any(has_k)) then
        do k_ = 1, ncls
          if (.not. has_k(k_)) cycle
          r2m = huge(1.0_dp)
          do i_ = 2, size(dsp, 2)
            if (scls(i_) == k_) r2m = min(r2m, sum(dsp(:, i_)**2))
          end do
          do i_ = 2, size(dsp, 2)
            if (scls(i_) == k_ .and. sum(dsp(:, i_)**2) <= r2m + 1.0e-9_dp) cl(i_) = ncls + k_
          end do
        end do
      else
        r2m = huge(1.0_dp)
        do i_ = 2, size(dsp, 2); r2m = min(r2m, sum(dsp(:, i_)**2)); end do
        do i_ = 2, size(dsp, 2)
          if (sum(dsp(:, i_)**2) <= r2m + 1.0e-9_dp) cl(i_) = ncls + cl(i_)
        end do
      end if
    end subroutine shell_classes_dsp

    !> The site operations of the cluster of sublattice kp that keep the classes cl (gidx: indices into sgr(kp)%perm, identity
    !> first), their orbital matrices d (asr_crystal_symmetry::orbital_matrix) and which of them are monomial
    subroutine class_group(kp, cl, gidx_, d, mono_)
      integer, intent(in) :: kp, cl(:)
      integer, allocatable, intent(out) :: gidx_(:)
      real(dp), allocatable, intent(out) :: d(:, :, :)
      logical, allocatable, intent(out) :: mono_(:)
      call class_group_ops(sgr(kp)%perm, sgr(kp)%rot, cl, gidx_, d, mono_)
    end subroutine class_group

    subroutine class_group_ops(pm, rt, cl, gidx_, d, mono_)
      integer, intent(in) :: pm(:, :), cl(:)
      real(dp), intent(in) :: rt(:, :, :)
      integer, allocatable, intent(out) :: gidx_(:)
      real(dp), allocatable, intent(out) :: d(:, :, :)
      logical, allocatable, intent(out) :: mono_(:)
      integer :: k_
      gidx_ = pack([(k_, k_=1, size(pm, 2))], [(all(cl(pm(:, k_)) == cl), k_=1, size(pm, 2))])
      if (gidx_(1) /= 1) error stop 'class_scf: the identity does not keep the classes'
      allocate (d(9, 9, size(gidx_)), mono_(size(gidx_)))
      do k_ = 1, size(gidx_)
        call orbital_matrix(rt(:, :, gidx_(k_)), d(:, :, k_)); mono_(k_) = is_monomial(rt(:, :, gidx_(k_)))
      end do
    end subroutine class_group_ops

    !> Bloch spectral function A(k, E) at the final potential (asr_spectral, route edc).  Per translation-inversion origin class:
    !> the operator reduced with D2h, the target sites within rmax and their nodes; per spin, origin-component pair, start orbital
    !> and energy one projected chain (edc_green_proj).  The diagonal elements <G_{0L, RL}> are averaged over the centre species,
    !> symmetrised between partner origins, and summed over R with the Bloch phases and origin site weights.
    subroutine spectral_final()
      type(sub_graph_t), allocatable :: ogr(:)
      type(class_op_t), allocatable :: sop(:)
      type(class_alloy_t), allocatable :: ms(:, :)
      type(edn_t), allocatable :: es(:, :)
      integer, allocatable :: bclass(:), orig_b(:), orig_sign(:), orig_sub(:), oc_orig(:), oc_comp(:)
      integer, allocatable :: clsb(:), clsr(:), gix(:), ntk(:), tsite(:, :), tmir(:, :), tbase(:, :), tsign(:, :), tcls(:, :)
      integer, allocatable :: tnode(:, :, :), tgop(:, :, :), tn_(:, :), tg_(:, :), shell(:, :), kop(:), t0(:)
      real(dp), allocatable :: orig_w(:), tdisp(:, :, :), tw(:, :, :, :, :), ee(:), v(:), r2sh(:), wwin(:, :)
      complex(dp), allocatable :: zs(:), gdg(:, :, :, :, :), gdg4(:, :, :, :, :)
      real(dp) :: rmax2, gmax, gmax0, a_, de, amax, rmin_s
      integer :: norig, noc, oc_, ne, nst, nt, ntmax, nsm, ks, cl_, t_, uo, u2, kp, pp_, gg, i, j, l_, ik, nkp, nsh, nopd, it_
      real(dp), allocatable :: dsp(:, :, :)
      logical, allocatable :: mn_(:)
      if (ed3) error stop 'class_scf: the spectral function needs route = edc'
      if (.not. spectral_structure_ok(st)) error stop &
        'class_scf: spectral function: a sublattice holds symmetry-inequivalent sites; split them into separate sublattices'
      nst = inp%spectral%nsteps; if (nst <= 0) nst = nlev
      ne = inp%spectral%ne; nsm = maxval(inp%nsc)
      call all_params()
      call spectral_origin_classes(st, sops, bclass, orig_b, orig_sign, orig_w)
      norig = size(orig_b)
      rmin_s = inp%spectral%rmax
      allocate (orig_sub(norig), ogr(norig))
      do uo = 1, norig
        kp = st%label(orig_b(uo)); orig_sub(uo) = kp
        if (orig_b(uo) == sgr(kp)%cg%b0) then
          ogr(uo) = sgr(kp)
        else
          call crystal_graph(st, sops, inp%rcl, ngen, ogr(uo)%cg, smsg, start=orig_b(uo), rmin=rmin_s)
          if (len(smsg) > 0) error stop smsg
          call graph_distances(ogr(uo)%cg%g, 1, ogr(uo)%dist)
          ogr(uo)%flip = ogr(uo)%dist <= inp%flip_radius; ogr(uo)%flip(1) = .false.
          ogr(uo)%cls = st%label(ogr(uo)%cg%basis)
          ogr(uo)%stype = styp(ogr(uo)%cg%basis)
          call graph_site_ops(ogr(uo)%cg, st, sops, ogr(uo)%perm, ogr(uo)%rot, gix)
        end if
      end do
      ! ---- targets per origin class: cluster sites within rmax of the origin, and their partner in the target's origin class
      rmax2 = 4.0_dp * inp%spectral%rmax**2
      allocate (ntk(norig), t0(norig)); ntk = 0; t0 = 0
      do uo = 1, norig
        ntk(uo) = count([(4.0_dp * sum(ogr(uo)%cg%disp(:, i)**2) <= rmax2 + 1.0e-9_dp, i=1, ogr(uo)%cg%n)])
      end do
      ntmax = maxval(ntk)
      allocate (tsite(ntmax, norig), tdisp(3, ntmax, norig), tmir(ntmax, norig), tbase(ntmax, norig), tsign(ntmax, norig), &
                tcls(ntmax, norig))
      tsite = 1; tdisp = 0.0_dp; tmir = 0; tbase = 1; tsign = 1; tcls = 1
      res%spec_maxhop = 0
      do uo = 1, norig
        associate (cgk => ogr(uo)%cg)
          nt = 0
          do i = 1, cgk%n
            if (4.0_dp * sum(cgk%disp(:, i)**2) <= rmax2 + 1.0e-9_dp) then
              nt = nt + 1; tsite(nt, uo) = i; tdisp(:, nt, uo) = cgk%disp(:, i)
              if (i == 1) t0(uo) = nt
              res%spec_maxhop = max(res%spec_maxhop, ogr(uo)%dist(i))
            end if
          end do
          tbase(1:nt, uo) = bclass(cgk%basis(tsite(1:nt, uo)))
          tsign(1:nt, uo) = orig_sign(cgk%basis(tsite(1:nt, uo)))
        end associate
      end do
      if (any(t0 == 0)) error stop 'class_scf: the origin is not a spectral target'
      do uo = 1, norig
        do t_ = 1, ntk(uo)
          u2 = tbase(t_, uo)
          tmir(t_, uo) = 0
          do j = 1, ntk(u2)
            if (tbase(j, u2) == uo .and. all(abs(tdisp(:, j, u2) + tsign(t_, uo) * tdisp(:, t_, uo)) < 1.0e-9_dp)) then
              tmir(t_, uo) = j; exit
            end if
          end do
          ! Fallback when no partner at -R exists: for two sites in one origin class related by inversion, reciprocity maps
          ! G(b, b + R) to G(b', b' - R) and inversion maps (b', -R) back to (b, +R) with D_LL(-I)^2 = 1, while spectral_unfold
          ! sums cos(2 pi k.R) which is even in R, so the element is its own partner (tmir = t_).
          if (tmir(t_, uo) == 0) tmir(t_, uo) = t_
        end do
      end do
      res%spec_nt = ntk(1); res%spec_hoplim = 2 * (ngen - 1); res%spec_nsteps = nst
      ! ---- per origin class: D2h-reduced operator with the target nodes
      allocate (sop(norig), tnode(nsm, ntmax, norig), tgop(nsm, ntmax, norig), res%spec_nodes(norig))
      tnode = 0; tgop = 1
      do uo = 1, norig
        nt = ntk(uo); kp = orig_sub(uo)
        clsb = ogr(uo)%cls
        clsr = clsb; if (inp%sro) call shell_classes_dsp(kp, ogr(uo)%cg%disp, ogr(uo)%cls, clsr)
        call class_group_ops(ogr(uo)%perm, ogr(uo)%rot, clsr, gix, dsp, mn_)
        if (allocated(kop)) deallocate (kop)
        allocate (kop(size(gix))); nopd = 0                  ! D2h: the operations with a diagonal matrix (identity first)
        do j = 1, size(gix)
          if (maxval(abs(ogr(uo)%rot(:, :, gix(j)))) <= 1.0_dp + 1.0e-9_dp .and. &
              all([(abs(abs(ogr(uo)%rot(i, i, gix(j))) - 1.0_dp) < 1.0e-9_dp, i=1, 3)])) then
            nopd = nopd + 1; kop(nopd) = j
          end if
        end do
        sop(uo)%d = dsp(:, :, kop(1:nopd))
        call build_augcls_operator(ogr(uo)%cg%g, mdl, clsr, ogr(uo)%perm(:, gix(kop(1:nopd))), sop(uo)%d, 1, ngen, sop(uo)%op, &
                                   inp%maxrank, ogr(uo)%flip, tsite=tsite(1:nt, uo), tnode=tn_, tgop=tg_, stype=ogr(uo)%stype)
        if (any(tn_(1, 1:nt) == 0)) error stop 'class_scf: a spectral target lies outside the node set (lower rmax or raise levels)'
        tnode(:, 1:nt, uo) = tn_; tgop(:, 1:nt, uo) = tg_; tcls(1:nt, uo) = clsr(tsite(1:nt, uo))
        res%spec_nodes(uo) = sop(uo)%op%nn
        if (vb) print '(a,i0,a,i0,a,i0,a,i0)', '  spectral: origin class ', uo, ', D2h operations ', nopd, ', nodes ', &
          sop(uo)%op%nn, ', target sites ', nt
      end do
      if (vb .and. res%spec_maxhop > res%spec_hoplim) print '(a,i0,a,i0,a)', '  WARNING: spectral targets reach ', &
        res%spec_maxhop, ' graph hops; the interior of the node set ends at ', res%spec_hoplim, ' (lower rmax or raise levels)'
      ! ---- species parameters per (component, spin) and dressing weights per (origin, component) task
      noc = sum([(inp%nsc(orig_sub(uo)), uo=1, norig)])
      allocate (oc_orig(noc), oc_comp(noc), ms(ncomp, nsp), es(ncomp, nsp), tw(9, nsm, ntmax, noc, nsp))
      tw = 0.0_dp; oc_ = 0
      do uo = 1, norig
        do pp_ = 1, ncomp
          if (pcl(pp_) == orig_sub(uo)) then
            oc_ = oc_ + 1; oc_orig(oc_) = uo; oc_comp(oc_) = pp_
          end if
        end do
      end do
      do gg = 1, nsp
        do pp_ = 1, ncomp
          call edc_setup(mdl, pps(:, gg), alpha, pp_, es(pp_, gg)); ms(pp_, gg) = mdl
        end do
        do oc_ = 1, noc
          uo = oc_orig(oc_); pp_ = oc_comp(oc_); kp = orig_sub(uo)
          do t_ = 1, ntk(uo)
            cl_ = tcls(t_, uo)
            if (cl_ <= ncls) then
              call spectral_dressing(es(pp_, gg)%dc, es(pp_, gg)%del, mdl%spc(1:inp%nsc(cl_), cl_), &
                                     inp%conc(1:inp%nsc(cl_), cl_), tsite(t_, uo) == 1, tw(:, :, t_, oc_, gg))
            else
              ks = cl_ - ncls
              call spectral_dressing(es(pp_, gg)%dc, es(pp_, gg)%del, mdl%spc(1:inp%nsc(ks), cl_), &
                                     inp%sro_cond(pk(pp_), 1:inp%nsc(ks), kp, ks), tsite(t_, uo) == 1, tw(:, :, t_, oc_, gg))
            end if
          end do
        end do
      end do
      ! ---- energies, chains
      allocate (ee(ne), zs(ne))
      do j = 1, ne
        ee(j) = inp%spectral%emin + (inp%spectral%emax - inp%spectral%emin) * real(j - 1, dp) / real(max(ne - 1, 1), dp)
      end do
      zs = cmplx(ee, inp%spectral%eta, dp)
      allocate (gdg(9, ntmax, ne, nsp, norig), gdg4(9, ntmax, ne, nsp, norig)); gdg = (0.0_dp, 0.0_dp); gdg4 = (0.0_dp, 0.0_dp)
      chain_vec_par = .false.                                 ! one thread per chain (asr_chain::chain_vec_par)
      !$omp parallel
      call spectral_worker(sop, ms, es, zs, nst, oc_orig, oc_comp, tnode, tgop, tw, gdg, gdg4)
      !$omp end parallel
      chain_vec_par = .true.
      if (par_size() > 1) then
        call par_sum_z(gdg, size(gdg)); call par_sum_z(gdg4, size(gdg4))
      end if
      ! ---- diagnostics on the averaged elements
      gmax = maxval(abs(gdg))
      gmax0 = 0.0_dp
      do uo = 1, norig; gmax0 = max(gmax0, maxval(abs(gdg(:, t0(uo), :, :, uo)))); end do
      res%spec_rdef = 0.0_dp; res%spec_spread0 = 0.0_dp; res%spec_spread = 0.0_dp
      do uo = 1, norig
        res%spec_spread0 = max(res%spec_spread0, maxval(abs(gdg(:, t0(uo), :, :, uo) - gdg4(:, t0(uo), :, :, uo))))
        do t_ = 1, ntk(uo)
          if (tmir(t_, uo) /= t_ .or. tbase(t_, uo) /= uo) then
            res%spec_rdef = max(res%spec_rdef, maxval(abs(gdg(:, t_, :, :, uo) - gdg(:, tmir(t_, uo), :, :, tbase(t_, uo)))))
          end if
          if (t_ /= t0(uo)) res%spec_spread = max(res%spec_spread, maxval(abs(gdg(:, t_, :, :, uo) - gdg4(:, t_, :, :, uo))))
        end do
      end do
      res%spec_rdef = res%spec_rdef / max(gmax, tiny(1.0_dp))
      res%spec_spread0 = res%spec_spread0 / max(gmax0, tiny(1.0_dp))
      res%spec_spread = res%spec_spread / max(gmax, tiny(1.0_dp))
      ! ---- symmetrise between partner origins
      gdg4 = gdg
      do uo = 1, norig
        do t_ = 1, ntk(uo)
          if (tmir(t_, uo) /= t_ .or. tbase(t_, uo) /= uo) then
            gdg(:, t_, :, :, uo) = 0.5_dp * (gdg4(:, t_, :, :, uo) + gdg4(:, tmir(t_, uo), :, :, tbase(t_, uo)))
          end if
        end do
      end do
      ! shell decay of Im G relative to the on-site element
      allocate (r2sh(ntmax * norig), shell(ntmax, norig)); nsh = 0; shell = 0
      do uo = 1, norig
        do t_ = 1, ntk(uo)
          a_ = sum(tdisp(:, t_, uo)**2)
          do j = 1, nsh
            if (abs(r2sh(j) - a_) < 1.0e-9_dp) then; shell(t_, uo) = j; exit; end if
          end do
          if (shell(t_, uo) == 0) then; nsh = nsh + 1; r2sh(nsh) = a_; shell(t_, uo) = nsh; end if
        end do
      end do
      allocate (res%spec_rsh(nsh), res%spec_decay(nsh)); res%spec_rsh = sqrt(r2sh(1:nsh)); res%spec_decay = 0.0_dp
      gmax0 = 0.0_dp
      do uo = 1, norig; gmax0 = max(gmax0, maxval(abs(aimag(gdg(:, t0(uo), :, :, uo))))); end do
      do uo = 1, norig
        do t_ = 1, ntk(uo)
          res%spec_decay(shell(t_, uo)) = max(res%spec_decay(shell(t_, uo)), &
                                              maxval(abs(aimag(gdg(:, t_, :, :, uo)))) / max(gmax0, tiny(1.0_dp)))
        end do
      end do
      ! ---- Bloch sums along the path
      call spectral_kpath(inp%spectral%kpts, inp%spectral%nk, res%spec_k, res%spec_dist, res%spec_icorner)
      nkp = size(res%spec_dist)
      allocate (res%spec_e(ne), res%spec_a(ne, nkp, 0:2, nsp), res%spec_ac(ne, nkp, ncls, nsp), v(ntmax), wwin(ntmax, norig))
      wwin = 0.0_dp
      do uo = 1, norig
        do t_ = 1, ntk(uo); wwin(t_, uo) = spectral_window(inp%spectral%window, norm2(tdisp(:, t_, uo)), inp%spectral%rmax); end do
      end do
      res%spec_e = ee; res%spec_a = 0.0_dp; res%spec_ac = 0.0_dp
      do gg = 1, nsp
        do uo = 1, norig
          nt = ntk(uo); kp = orig_sub(uo)
          do j = 1, ne
            do l_ = 0, 2
              v(1:nt) = 0.0_dp
              do it_ = 1, 9
                if (lorb9(it_) == l_) v(1:nt) = v(1:nt) - aimag(gdg(it_, 1:nt, j, gg, uo)) / pi
              end do
              v(1:nt) = v(1:nt) * wwin(1:nt, uo)
              do ik = 1, nkp
                a_ = orig_w(uo) * spectral_unfold(tdisp(:, 1:nt, uo), v(1:nt), res%spec_k(:, ik))
                res%spec_a(j, ik, l_, gg) = res%spec_a(j, ik, l_, gg) + a_
                res%spec_ac(j, ik, kp, gg) = res%spec_ac(j, ik, kp, gg) + a_
              end do
            end do
          end do
        end do
      end do
      if (all(inp%nsc == 1)) call ordered_bands(ms, es)
      amax = maxval(sum(res%spec_a, dim=3))
      res%spec_nneg = count(sum(res%spec_a, dim=3) < -1.0e-9_dp * amax)
      res%spec_amin = minval(sum(res%spec_a, dim=3)) / max(amax, tiny(1.0_dp))
      allocate (res%spec_sum(2, nsp)); res%spec_sum = 0.0_dp
      de = (inp%spectral%emax - inp%spectral%emin) / real(max(ne - 1, 1), dp)
      do gg = 1, nsp
        do j = 1, ne
          a_ = merge(0.5_dp, 1.0_dp, j == 1 .or. j == ne) * de
          res%spec_sum(1, gg) = res%spec_sum(1, gg) + a_ * sum(res%spec_a(j, 1, :, gg))
          res%spec_sum(2, gg) = res%spec_sum(2, gg) + a_ * sum(res%spec_a(j, nkp, :, gg))
        end do
      end do
      if (vb) then
        print '(a,i0,a,f7.4,a,i0,a,i0,a,i0,a)', '  spectral: ', ne, ' energies, eta ', inp%spectral%eta, ' Ry, ', nst, &
          ' steps, ', nkp, ' k points, ', nt, ' target sites'
        print '(a,es9.2,a,es9.2,a,es9.2)', '  spectral: reciprocity defect ', res%spec_rdef, ', step spread (4 steps) off-site ', &
          res%spec_spread, ' on-site ', res%spec_spread0
        print '(a,a,a,i0,a,es9.2)', '  spectral: window ', trim(inp%spectral%window), ';  points with A < 0: ', res%spec_nneg, &
          ', min A / max A ', res%spec_amin
        print '(a,*(f7.3))', '  spectral: shell radius (alat)        ', res%spec_rsh
        print '(a,*(f7.3))', '  spectral: max |Im G| / on-site value ', res%spec_decay
        print '(a,*(f8.3))', '  spectral: energy integral of A at the first and last k (per spin; the number of bands inside the window, at most 9) ', res%spec_sum
      end if
      do uo = 1, norig
        call augcls_free(sop(uo)%op)                     ! shared node tables of the spectral operators (asr_augcls)
      end do
    end subroutine spectral_final

    !> Eigenvalues and unfolding weights of H^(2)(k) = E_nu + h - h o h on any ordered structure (one species per sublattice).
    !> Basis sites related by pure centering translations of sops are reduced to nb_p primitive representatives, giving 9 * nb_p
    !> bands in the tau-gauge Bloch basis with unfolding weights W_m(k) = (1 / nb_p) sum_L |sum_r psi_m(r, L; k)|^2.
    subroutine ordered_bands(ms, es)
      type(class_alloy_t), intent(in) :: ms(:, :)
      type(edn_t), intent(in) :: es(:, :)
      complex(dp), allocatable :: sk(:, :), hk(:, :), cw(:)
      real(dp), allocatable :: rw(:), ev(:), cv(:), dv(:), ov(:), ev0(:)
      integer :: tclass(st%nb), rep_b(st%nb)
      complex(dp) :: csum
      real(dp) :: dsl(3), eye(3, 3), w_
      integer :: nb_, ndim, lwk_, r1_, r2_, b1_, b2_, g_, t1_, p_, i0_, j0_, ks_, ns_, info_, m_, n_, l_, gg, ik, nkp
      eye = 0.0_dp; do l_ = 1, 3; eye(l_, l_) = 1.0_dp; end do
      tclass = 0; nb_ = 0
      do b1_ = 1, st%nb
        if (tclass(b1_) /= 0) cycle
        nb_ = nb_ + 1; rep_b(nb_) = b1_; tclass(b1_) = nb_
        do b2_ = b1_ + 1, st%nb
          if (st%label(b2_) /= st%label(b1_) .or. tclass(b2_) /= 0) cycle
          do g_ = 1, size(sops)
            if (sops(g_)%perm(b1_) == b2_ .and. maxval(abs(sops(g_)%rot - eye)) < 1.0e-7_dp) then
              tclass(b2_) = nb_; exit
            end if
          end do
        end do
      end do
      ndim = 9 * nb_; lwk_ = max(36, 4 * ndim)
      nkp = size(res%spec_dist)
      ns_ = sgr(1)%cg%nslot
      allocate (sk(ndim, ndim), hk(ndim, ndim), cw(lwk_), rw(3 * ndim), ev(ndim))
      allocate (cv(ndim), dv(ndim), ov(ndim), ev0(ndim))
      allocate (res%spec_bands(ndim, nkp, nsp), res%spec_bwgt(ndim, nkp, nsp))
      do gg = 1, nsp
        do r1_ = 1, nb_
          b1_ = rep_b(r1_); p_ = st%label(b1_); i0_ = 9 * (r1_ - 1)
          cv(i0_ + 1:i0_ + 9) = ms(p_, gg)%c(:, p_)
          dv(i0_ + 1:i0_ + 9) = ms(p_, gg)%d(:, p_)
          ov(i0_ + 1:i0_ + 9) = ms(p_, gg)%o(:, p_)
          ev0(i0_ + 1:i0_ + 9) = es(p_, gg)%enu(:, p_)
        end do
        do ik = 1, nkp
          sk = (0.0_dp, 0.0_dp)
          do r1_ = 1, nb_
            b1_ = rep_b(r1_); t1_ = styp(b1_); i0_ = 9 * (r1_ - 1)
            sk(i0_ + 1:i0_ + 9, i0_ + 1:i0_ + 9) = ms(1, gg)%block0(:, :, t1_)
            do ks_ = 1, ns_
              b2_ = sgr(1)%cg%btgt(ks_, b1_)
              if (b2_ == 0) cycle
              r2_ = tclass(b2_); j0_ = 9 * (r2_ - 1)
              dsl = sgr(1)%cg%bvec(:, ks_)
              sk(i0_ + 1:i0_ + 9, j0_ + 1:j0_ + 9) = sk(i0_ + 1:i0_ + 9, j0_ + 1:j0_ + 9) + &
                ms(1, gg)%block(:, :, ks_, t1_) * exp(cmplx(0.0_dp, 2.0_dp * pi * dot_product(res%spec_k(:, ik), dsl), dp))
            end do
          end do
          sk = 0.5_dp * (sk + conjg(transpose(sk)))
          do n_ = 1, ndim
            do m_ = 1, ndim; hk(m_, n_) = dv(m_) * sk(m_, n_) * dv(n_); end do
            hk(n_, n_) = hk(n_, n_) + cv(n_)
          end do
          sk = hk
          do n_ = 1, ndim
            do m_ = 1, ndim; hk(m_, n_) = sk(m_, n_) - sum(sk(m_, :) * ov * sk(:, n_)); end do
            hk(n_, n_) = hk(n_, n_) + ev0(n_)
          end do
          hk = 0.5_dp * (hk + conjg(transpose(hk)))
          call zheev('V', 'U', ndim, hk, ndim, ev, cw, lwk_, rw, info_)
          if (info_ /= 0) error stop 'class_scf: zheev failed (spectral bands)'
          res%spec_bands(:, ik, gg) = ev
          do m_ = 1, ndim
            w_ = 0.0_dp
            do l_ = 1, 9
              csum = (0.0_dp, 0.0_dp)
              do r1_ = 1, nb_; csum = csum + hk(9 * (r1_ - 1) + l_, m_); end do
              w_ = w_ + abs(csum)**2
            end do
            res%spec_bwgt(m_, ik, gg) = w_ / real(nb_, dp)
          end do
        end do
      end do
    end subroutine ordered_bands

    !> one thread's share of the spectral chains (inside the parallel region of spectral_final): tasks (energy, start orbital,
    !> origin-component pair, spin), a private copy of the origin-class operator per thread
    subroutine spectral_worker(sop, ms, es, zs, nst, oc_orig, oc_comp, tnode, tgop, tw, gdg, gdg4)
      type(class_op_t), intent(in) :: sop(:)
      type(class_alloy_t), intent(in) :: ms(:, :)
      type(edn_t), intent(in) :: es(:, :)
      complex(dp), intent(in) :: zs(:)
      integer, intent(in) :: nst, oc_orig(:), oc_comp(:), tnode(:, :, :), tgop(:, :, :)
      real(dp), intent(in) :: tw(:, :, :, :, :)
      complex(dp), intent(inout) :: gdg(:, :, :, :, :), gdg4(:, :, :, :, :)
      type(augcls_op), allocatable :: lop(:)
      type(edn_work_t), allocatable :: lwk(:)
      complex(dp), allocatable :: gp(:, :), gp4(:, :)
      complex(dp) :: gc
      integer :: cur(3, size(sop)), tk, j_, o_, oc_, u_, p_, k_, q_, s_, rk, nr, ne, noc, norig
      rk = par_rank(); nr = par_size(); ne = size(zs); noc = size(oc_orig); norig = size(sop)
      allocate (lop(norig), lwk(norig), gp(9, size(tnode, 2)), gp4(9, size(tnode, 2))); cur = -1
      !$omp do schedule(dynamic)
      do tk = 1, ne * 9 * noc * nsp
        if (mod(tk - 1 + (tk - 1) / ne, nr) /= rk) cycle
        j_ = mod(tk - 1, ne) + 1; o_ = mod((tk - 1) / ne, 9) + 1
        oc_ = mod((tk - 1) / (9 * ne), noc) + 1; s_ = (tk - 1) / (9 * ne * noc) + 1
        u_ = oc_orig(oc_); p_ = oc_comp(oc_); k_ = pcl(p_)
        if (cur(1, u_) < 0) then
          lop(u_) = sop(u_)%op
          nullify (lop(u_)%wk); allocate (lop(u_)%wk(size(sop(u_)%op%wk, 1), size(sop(u_)%op%wk, 2)))
          call edn_work_init(lwk(u_), lop(u_)%dim())
          cur(:, u_) = 0
        end if
        if (cur(1, u_) /= p_ .or. cur(2, u_) /= s_) then
          call augcls_set_params(lop(u_), ms(p_, s_)); call augcls_freeze_site(lop(u_), 1, p_)
          if (inp%sro) then
            do q_ = 1, ncls
              call augcls_set_class_conc(lop(u_), ncls + q_, inp%sro_cond(pk(p_), 1:inp%nsc(q_), k_, q_))
            end do
          end if
          cur(1, u_) = p_; cur(2, u_) = s_; cur(3, u_) = -1
        end if
        if (cur(3, u_) /= o_) then; call augcls_set_start(lop(u_), o_, sop(u_)%d); cur(3, u_) = o_; end if
        call edc_green_proj(lop(u_), lwk(u_), es(p_, s_), zs(j_), o_, nst, tnode(:, :, u_), tgop(:, :, u_), &
                            tw(:, :, :, oc_, s_), gc, gp, gp4)
        !$omp critical (spectral_acc)
        gdg(o_, :, j_, s_, u_) = gdg(o_, :, j_, s_, u_) + xp(p_) * gp(o_, :)
        gdg4(o_, :, j_, s_, u_) = gdg4(o_, :, j_, s_, u_) + xp(p_) * gp4(o_, :)
        !$omp end critical (spectral_acc)
      end do
      !$omp end do
      do u_ = 1, norig
        if (cur(1, u_) >= 0) deallocate (lop(u_)%wk)
      end do
    end subroutine spectral_worker

    subroutine dos_final()
      complex(dp), allocatable :: zs(:), gz(:, :, :, :)
      integer :: ie, o, gg, pp_
      allocate (res%edos_e(inp%ndos), res%edos(inp%ndos, 0:2, ncomp, nsp)); res%edos = 0.0_dp
      do ie = 1, inp%ndos
        res%edos_e(ie) = inp%dos_emin + (inp%dos_emax - inp%dos_emin) * (ie - 1) / real(max(inp%ndos - 1, 1), dp)
      end do
      call all_params()
      zs = cmplx(res%edos_e, inp%dos_eta, dp)
      allocate (gz(inp%ndos, nomax, ncomp, nsp))
      call chain_greens(zs, gz, tbk=tb_dos)
      do gg = 1, nsp
        do pp_ = 1, ncomp
          do o = 1, norb(pcl(pp_))
            res%edos(:, lorbk(o, pcl(pp_)), pp_, gg) = res%edos(:, lorbk(o, pcl(pp_)), pp_, gg) &
                                                       - worbk(o, pcl(pp_)) * aimag(gz(:, o, pp_, gg)) / pi
          end do
        end do
      end do
    end subroutine dos_final

    !> component potentials (with the nucleus) and E_nu: nspin, ncomp; per component z n r(n), E_nu(0:2) per spin, v(r_k, spin)
    subroutine write_potentials()
      integer :: u, pp_, gg, kk
      if (.not. par_root()) return
      open (newunit=u, file=trim(inp%pot_out), status='replace', action='write')
      write (u, '(a,es25.16)') '# asr sublattice component potentials (Ry): nspin ncomp; per component: z n r(n); E_nu(0:2) per spin; &
                             &v(r_k, spin);  E_F= ', ef
      write (u, '(2i5)') nsp, ncomp
      do pp_ = 1, ncomp
        write (u, '(f8.2,i8,es25.16)') inp%spec(psp(pp_))%z, cr(pp_)%mesh%n, cr(pp_)%mesh%r(cr(pp_)%mesh%n)
        do gg = 1, nsp; write (u, '(3es25.16)') enu(:, gg, pp_); end do
        do kk = 1, cr(pp_)%mesh%n; write (u, '(2es25.16)') cr(pp_)%v(kk, :); end do
      end do
      if (res%uform) then                          ! older readers stop before this block
        write (u, '(a)') '# u: drive of the SRO charge deviations (e), one value per component'
        write (u, '(*(es25.16))') ucur
      end if
      close (u)
    end subroutine write_potentials

    subroutine read_potentials()
      integer :: u, pp_, gg, kk, nf, nspf, ncf, ios
      real(dp) :: zf, rf, vv2(2), ee(0:2, 2)
      character(len=512) :: line
      open (newunit=u, file=trim(inp%pot_in), status='old', action='read', iostat=ios)
      if (ios /= 0) error stop 'class_scf: cannot open start_potential'
      read (u, '(a)') line
      kk = index(line, 'E_F=')
      if (kk > 0) then
        read (line(kk + 4:), *, iostat=ios) ef_in
        if (ios /= 0) ef_in = huge(1.0_dp)
      end if
      read (u, *) nspf, ncf
      if (nspf > nsp) error stop 'class_scf: start_potential is spin polarised, the run is not'
      if (ncf /= ncomp) error stop 'class_scf: start_potential has a different number of components'
      do pp_ = 1, ncomp
        read (u, *) zf, nf, rf
        if (nf /= cr(pp_)%mesh%n .or. abs(zf - inp%spec(psp(pp_))%z) > 1.0e-6_dp .or. abs(rf - cr(pp_)%mesh%r(nf)) > 1.0e-8_dp * rf) &
          error stop 'class_scf: start_potential components / radial meshes differ from this input'
        do gg = 1, nspf; read (u, *) ee(:, gg); end do
        if (nspf == 1) ee(:, 2) = ee(:, 1)
        do kk = 1, nf
          read (u, *) vv2(1:nspf)
          if (nspf == 1) vv2(2) = vv2(1)
          cr(pp_)%v(kk, :) = vv2(1:nsp)
        end do
        enu(:, 1:nsp, pp_) = ee(:, 1:nsp)
      end do
      read (u, '(a)', iostat=ios) line              ! optional: the SRO drive u (written by runs with SRO)
      if (ios == 0 .and. index(line, '# u:') == 1 .and. inp%sro) then
        read (u, *, iostat=ios) vfl
        if (ios == 0) then; ucur = vfl; uload = .true.; end if
        vfl = 0.0_dp
      end if
      close (u)
      if (vb) print '(3a)', '  start potentials from ', trim(inp%pot_in)
      if (vb .and. uload) print '(a)', '  SRO drive u from the start potential'
    end subroutine read_potentials
  end subroutine class_scf

end module asr_class_scf
