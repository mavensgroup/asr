!> Self-consistent binary alloy A_x B_(1-x) from the augmented-space recursion (single-site charge level,
!> Screened Impurity Model for the electrostatics).  Collinear spin polarisation (nspin = 1 or 2).  Every
!> parameter is generated here.
!>
!> Per iteration:  species/spin potentials -> potential parameters (phi, phidot, phiddot) -> alloy operator
!> (H = E + h - h o h, common screening, spin-dependent parameters on one node set) -> four recursions per spin
!> (s, p_z, t2g, eg start orbitals) with stored Lanczos bases -> species-resolved Green functions
!> G_A = <<n G>>/x, G_B = (G - <<n G>>)/(1-x) -> contour moments M_k = int (E - E_nu)^k N(E) dE (common E_F from
!> the electron count) -> orthogonal-representation density
!>    rho_l = M0 phi^2 + 2 M1 phi phidot + M2 (phidot^2 + phi phiddot)
!> -> relaxed core, Hartree, xc (spin-polarised libxc), SIM Madelung -> Anderson mixing of the smooth
!> potentials and E_nu.
module asr_alloy_scf
  use asr_kinds, only: dp, pi
  use asr_radial
  use asr_xc
  use asr_atom, only: cumint, poisson_hartree, shell_t
  use asr_potpar, only: potpar_l_t, potpar_at_enu, clamp_enu
  use asr_lmto_ham, only: potpar_t
  use asr_asa, only: asa_class_t, asa_input_t, asa_class_result_t, atomic_start, core_shells_for
  use asr_lmto_cell, only: alpha_from_sigma, cell_t, cluster_t, build_cluster, setup_cell, recip_vectors
  use asr_strconst, only: avg_ws_radius, standard_alpha
  use asr_augspace, only: graph_t, graph_distances
  use asr_augorb
  use asr_augsym
  use asr_endep
  use asr_ed
  use asr_pfun, only: pfun_exact
  use asr_symmetry, only: d4h_ops, orbital_rep
  use asr_lmto_alloy, only: build_lmto_alloy, set_alloy_params
  use asr_partial, only: lanczos_basis, transition_green
  use asr_contour, only: contour_nodes, contour_nodes_closed
  use asr_greens, only: cf_green, sqrt_tail, tail_params_avg2
  use asr_terminator, only: bp_asymptote
  use asr_mixing, only: anderson_t, anderson_init, anderson_step
  use asr_lmto_realization, only: lmto_real_op, build_lmto_ordered, lmto_hph
  use asr_recursion, only: lanczos
  use asr_ewald, only: madelung_phi
  implicit none
  private
  public :: alloy_input_t, alloy_result_t, alloy_scf, alloy_ldos, alloy_alpha

  type :: alloy_input_t
    type(asa_class_t) :: spec(2)                 !! 1 = A (concentration x), 2 = B
    real(dp) :: x = 0.5_dp
    real(dp) :: alat = 0.0_dp                    !! Bohr
    character(len=8) :: lattice = 'fcc', xc = 'vbh'
    real(dp) :: sigma = 0.7_dp, rcl = 1.2_dp     !! hard-core screening (screening = 'hardcore'), kept S^alpha range (alat)
    character(len=10) :: screening = 'auto'      !! 'standard' TB-LMTO alpha (0.3485, 0.05303, 0.010714) | 'hardcore' (sigma) |
                                                 !! 'auto' = standard for the energy-independent H^(2), hardcore for the energy-dependent P(z) - S
    real(dp) :: rinv = 0.0_dp                     !! inversion cluster radius for S^alpha (alat); <= rcl: invert on the rcl cluster
    integer :: lsuper = 24                       !! supercell (conventional cells) of the ASR graph
    integer :: nlevels = 8, maxrank = 2, flip_radius = 2
    integer :: nrec = 0                          !! Lanczos steps (0: = nlevels); > nlevels runs on the generated node set as a hard wall (AMscf: map
                                                 !! exact to 6 hops, 12 steps)
    integer :: ncont = 60                        !! contour points
    real(dp) :: beta = 0.3_dp, tol = 1.0e-5_dp
    !> Screened Impurity Model (Ruban & Skriver, PRB 66, 024201 (2002); thesis sec. 2.6.6): the net charge qnet = N - Z of a
    !> sphere is screened by the first coordination shell (radius R1), V = -sim_beta 2 qnet / R1 (Ry), added to the
    !> vacuum-boundary sphere Hartree potential (which already holds +2 qnet / S).  sim_beta = 1: full first-shell screening.
    real(dp) :: sim_beta = 1.0_dp
    integer :: maxit = 60, mix_history = 6, nspin = 1
    real(dp) :: zeta0 = 0.1_dp                   !! initial spin polarisation of the atomic start (nspin = 2); 0.4 runs away for Fe
    logical :: relativistic = .true.
    logical :: symmetry = .true.                !! D4h point-group reduction of the augmented space (exact)
    logical :: ordered = .false.                !! ordered B2 on the bcc graph (A even / B odd sublattice): real-space recursion, no augmented space, exact Madelung
    logical :: energy_dependent = .false.       !! recursion on P(z) - S per complex energy (needs symmetry = true)
    logical :: ed_asr = .false.               !! with energy_dependent: AMscf's literal form H(E) = E - T(E - E_nu - h)T - T hoh T (asr_ed), one
                                                 !! recursion per central species, start site frozen; default screening = standard
    real(dp) :: ed_eta = 0.01_dp                  !! energy-dependent route: the contour arc ends at Im z = ed_eta
    logical :: ed_scan = .false.                !! diagnostic: print N(E_F) on a grid in iteration 1, then stop
    integer :: ndos = 0                          !! energy-dependent route: DOS grid computed inside the SCF (final iteration)
    character(len=8) :: terminator = 'avg2'      !! energy-independent route tail: 'avg2' (last two coefficients) | 'bp' (Beer-Pettifor)
    integer :: order = 2                         !! energy-independent Hamiltonian: 2 = E_nu + h - h o h; 3 = + h o h o h - (h p h h + h h p h)/2
    real(dp) :: order_switch = 1.0e-3_dp          !! order 3: iterate with order 2 until dV < order_switch (the cubic terms are unreliable
                                                 !! far from self-consistency, E_nu away from the band centres), then with order 3
    logical :: ed_exact = .false.                !! energy-dependent route: exact P_l(z) from the radial equation at each contour energy
    integer :: kkr_check = 0                     !! ordered ED: after the SCF, same-model k-space KKR-ASA check on an nk^3 mesh (0 = off)
    logical :: ed_exact_density = .true.         !! with ed_exact: exact density -1/pi Im int [G phi(z)^2 + T(z)] (T: single-site term, asr_pfun);
                                                 !! .false.: density from the energy moments of the exact G (phi, phidot, phiddot at E_nu)
    logical :: ed_third_order = .false.          !! energy-dependent route: Taylor P = (z - C + p eps^3)/Delta (experimental; SCF diverges)
    real(dp) :: depth_tol = 1.0e-5_dp             !! energy-dependent route: stop a chain when G changes < depth_tol (levels = maximum)
    character(len=8) :: ed_tail = 'hull'          !! energy-dependent route: tail of the chain, 'hull' | 'bp' | 'avg2' (asr_endep%ed_tail)
    character(len=8) :: charges = 'lowdin'       !! order 3: 'lowdin' (orthonormal-orbital weights) | 'sphere' (start vectors (1 - hph/2)|L>,
                                                 !! ASA sphere charges; loses ~0.03 e per cell at this order -> E_F shifts ~10 mRy in CuZn)
    character(len=10) :: density = 'expand'      !! M2 term: 'expand' = M2 (phidot^2 + phi phiddot) | 'stuttgart' = M2 (phidot^2 - p phi^2) (LMTO47 newrho)
    real(dp) :: dos_emin = -1.0_dp, dos_emax = 1.0_dp, dos_eta = 0.02_dp
    character(len=256) :: pot_in = ''            !! start from the species potentials and E_nu of this file (written by pot_out) instead of atoms
    character(len=256) :: pot_out = ''           !! write the final species potentials and E_nu (restart / fixed-potential diagnostics)
    logical :: pot_flip(2) = .false.             !! with pot_in, nspin = 2: swap the spin channels of species s (antiparallel starts)
  end type alloy_input_t

  type :: alloy_result_t
    type(asa_class_result_t) :: cls(2)
    real(dp) :: ef = 0.0_dp, dv = 0.0_dp, moment_cell = 0.0_dp, x = 0.5_dp
    integer :: iter = 0, nodes = 0, nspin = 1
    logical :: converged = .false.
    !! recursion coefficients of the last iteration (4 start orbitals x spin): for the LDOS at any energy
    real(dp), allocatable :: acoef(:, :, :), bcoef(:, :, :), ccoef(:, :, :)
    integer :: mm(4, 2) = 0
    logical :: ed = .false.
    logical :: ordered = .false.
    character(len=8) :: terminator = 'avg2'
    real(dp), allocatable :: oa(:, :, :, :), ob(:, :, :, :)       !! ordered mode: recursion coefficients (level, orbital, spin, sublattice)
    integer :: omm(4, 2, 2) = 0                                 !! ordered mode: number of coefficients (orbital, spin, sublattice)
    real(dp), allocatable :: hist_dv(:), hist_ef(:), hist_moment(:, :), hist_qnet(:, :)   !! SCF history, rows 1:iter
    real(dp), allocatable :: edos_e(:), edos(:, :, :, :)       !! energy-dependent route: DOS(e, l, species, spin)
    real(dp) :: cnrm(4, 2) = 1.0_dp, onrm(4, 2, 2) = 1.0_dp    !! |v|^2 of the start vectors (1 - hph/2)|L> (order 3; else 1)
  end type alloy_result_t

  integer, parameter :: iorb(4) = [1, 3, 5, 7]       ! s, p_z, t2g (xy), eg (3z^2-r^2) start orbitals
  integer, parameter :: lorb(4) = [0, 1, 2, 2]
  integer, parameter :: nvert = 8                     ! vertical-segment nodes of the closed ED contour
  integer, parameter :: lorb9(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]
  real(dp), parameter :: worb(4) = [1.0_dp, 3.0_dp, 3.0_dp, 2.0_dp]

contains

  !> Species- and l-resolved local DOS per atom and per spin channel (states/Ry) at E + i eta, from the
  !> stored recursion coefficients.  dos(ie, l, species, spin);  species 1 = A (Cu in Cu-x...), 2 = B.
  subroutine alloy_ldos(res, egrid, eta, dos)
    type(alloy_result_t), intent(in) :: res
    real(dp), intent(in) :: egrid(:), eta
    real(dp), intent(out) :: dos(size(egrid), 0:2, 2, res%nspin)
    real(dp) :: ainf, binf
    complex(dp) :: z, gt, gavg, fa
    integer :: ie, o, sg, ll, ip
    dos = 0.0_dp
    if (res%ordered) then
      do sg = 1, res%nspin
        do ip = 1, 2
          do o = 1, 4
            call tail_params(res%terminator, res%oa(:, o, sg, ip), res%ob(:, o, sg, ip), res%omm(o, sg, ip), ainf, binf)
            ll = lorb(o)
            do ie = 1, size(egrid)
              z = cmplx(egrid(ie), eta, dp)
              gavg = res%onrm(o, sg, ip) * cf_green(res%oa(:, o, sg, ip), res%ob(:, o, sg, ip), res%omm(o, sg, ip), z, sqrt_tail(z, ainf, binf))
              dos(ie, ll, ip, sg) = dos(ie, ll, ip, sg) - worb(o) * aimag(gavg) / pi
            end do
          end do
        end do
      end do
      return
    end if
    do sg = 1, res%nspin
      do o = 1, 4
        call tail_params(res%terminator, res%acoef(:, o, sg), res%bcoef(:, o, sg), res%mm(o, sg), ainf, binf)
        ll = lorb(o)
        do ie = 1, size(egrid)
          z = cmplx(egrid(ie), eta, dp)
          gt = sqrt_tail(z, ainf, binf)
          gavg = res%cnrm(o, sg) * cf_green(res%acoef(:, o, sg), res%bcoef(:, o, sg), res%mm(o, sg), z, gt)
          fa = sqrt(res%cnrm(o, sg)) * transition_green(res%acoef(:, o, sg), res%bcoef(:, o, sg), res%mm(o, sg), res%ccoef(:, o, sg), z, gt)
          dos(ie, ll, 1, sg) = dos(ie, ll, 1, sg) - worb(o) * aimag(fa / res%x) / pi
          dos(ie, ll, 2, sg) = dos(ie, ll, 2, sg) - worb(o) * aimag((gavg - fa) / (1.0_dp - res%x)) / pi
        end do
      end do
    end do
  end subroutine alloy_ldos

  !> nearest-neighbour distance in units of alat
  function nn_distance(kind) result(d)
    character(len=*), intent(in) :: kind
    real(dp) :: d
    select case (kind)
    case ('fcc'); d = sqrt(0.5_dp)
    case ('bcc'); d = 0.5_dp * sqrt(3.0_dp)
    case default; d = 1.0_dp
    end select
  end function nn_distance

  subroutine alloy_scf(inp, res, verbose)
    type(alloy_input_t), intent(in) :: inp
    type(alloy_result_t), intent(out) :: res
    logical, intent(in), optional :: verbose
    type(xc_lda_t) :: xf
    type(asa_input_t) :: ai
    type(potpar_t) :: pps(2, 2)                    ! (species, spin)
    type(potpar_l_t) :: ppl
    type(graph_t) :: g
    type(orb_alloy_t) :: mdl
    type(augorb_op) :: op
    type(augsym_op) :: sop, sopc(2)
    type(asr_ed_t) :: edp_c(2)
    logical :: asr_ed
    integer, allocatable :: perm(:, :), mats(:, :, :)
    real(dp), allocatable :: dmat(:, :, :)
    integer :: nd, nnodes, flip1, nfe, jf
    type(ed_params_t) :: eda(2), edb(2)
    type(ed_work_t) :: edw
    real(dp) :: alpha9(9), dosef(2, 2, 0:2), ef_old, nn_old, slope
    logical :: ed
    real(dp) :: md1, sflip
    type(anderson_t) :: mixer
    type(asa_class_result_t) :: cr(2)
    integer, allocatable :: coord(:, :), dist(:), offx(:)
    logical, allocatable :: flippable(:)
    real(dp), allocatable :: ub(:, :), y(:, :, :, :, :, :), U(:, :, :, :), acoef(:, :, :), bcoef(:, :, :), ccoef(:, :, :), v0(:), vv(:)
    real(dp), allocatable :: wmix(:)
    real(dp), allocatable :: xin(:), xout(:), xnext(:), nval(:, :), ncore(:, :), ntot(:), vh(:), vxc(:, :), rho(:, :), exc(:), tmp(:), cum(:)
    real(dp), allocatable :: yc(:, :, :), vtmp(:, :), etmp(:, :)
    real(dp) :: alpha(0:2), c, avw, enu(0:2, 2, 2), xs(2), ef, nelec, r1, lo, hi, nn, vm, dvc, dv, e, fs, s_s
    real(dp) :: mom(2, 2, 0:2, 0:2)                  ! (species, spin, k, l)
    integer :: nlev, ngen, mm(4, 2), s, sg, nsp, l, oc, n, it, i, k, nx, nmax, is, nod
    logical :: vb, core_ok
    logical :: exact_now
    real(dp) :: ebot_exact
    integer :: order_now
    real(dp) :: dv_last, cnrm(4, 2), onrm(4, 2, 2)
    real(dp), allocatable :: hvec(:), rhoex(:, :, :)
    logical :: sphq
    complex(dp), allocatable :: zc(:), wc(:)
    logical :: ordered
    logical, allocatable :: isa(:)
    type(lmto_real_op) :: rop
    real(dp), allocatable :: oa(:, :, :, :), ob(:, :, :, :), ta(:), tb(:), v0o(:), phimad(:, :)
    real(dp) :: vmo(2), qat2(2), bvec3(3, 3), tau2(3, 2)
    integer :: omm(4, 2, 2), isite_b, ip

    vb = .false.; if (present(verbose)) vb = verbose
    nsp = inp%nspin
    fs = 2.0_dp / real(nsp, dp)
    c = c_nonrel; if (inp%relativistic) c = c_ry
    ordered = inp%ordered
    xs = [inp%x, 1.0_dp - inp%x]
    if (ordered) xs = [0.5_dp, 0.5_dp]
    avw = avg_ws_radius(trim(inp%lattice)) * inp%alat
    alpha = alloy_alpha(inp)
    r1 = nn_distance(trim(inp%lattice)) * inp%alat
    nelec = xs(1) * real(inp%spec(1)%nval, dp) + xs(2) * real(inp%spec(2)%nval, dp)
    ngen = inp%nlevels                               ! augmented-space generation depth
    nlev = inp%nlevels; if (inp%nrec > 0) nlev = inp%nrec   ! recursion steps
    ed = inp%energy_dependent
    asr_ed = ed .and. inp%ed_asr
    omm = 0
    if (ed .and. .not. ordered .and. .not. inp%symmetry) error stop 'alloy_scf: energy_dependent requires symmetry = .true.'
    if (asr_ed .and. (ordered .or. inp%ed_exact)) error stop 'alloy_scf: the AMscf form is for the disordered alloy with parametrised potential functions'
    do i = 1, 9; alpha9(i) = alpha(lorb9(i)); end do
    call xc_init(xf, inp%xc)
    ! ---- meshes and atomic starts
    nmax = 0; enu = 0.0_dp
    do s = 1, 2
      s_s = merge(inp%spec(s)%sws, avw, inp%spec(s)%sws > 0.0_dp)
      call make_radmesh(inp%spec(s)%z, s_s, 0.012_dp, cr(s)%mesh, c)
      n = cr(s)%mesh%n; nmax = max(nmax, n)
      allocate (cr(s)%v(n, nsp), cr(s)%nrad(n, nsp))
      ai%alat = inp%alat; ai%xc = inp%xc; ai%relativistic = inp%relativistic; ai%nspin = nsp; ai%zeta0 = inp%zeta0
      if (allocated(ai%cls)) deallocate (ai%cls)
      allocate (ai%cls(1)); ai%cls(1) = inp%spec(s)
      allocate (vtmp(n, nsp), etmp(0:2, nsp))
      call atomic_start(ai, 1, cr(s)%mesh, c, xf, nsp, vtmp, etmp)
      cr(s)%v = vtmp
      do sg = 1, nsp; enu(:, sg, s) = etmp(:, sg); end do
      deallocate (vtmp, etmp)
    end do
    if (len_trim(inp%pot_in) > 0) call read_potentials()
    allocate (y(2, 0:2, nmax, 0:2, 2, 2)); y = 0.0_dp
    ! ---- first parameters, structure (graph, blocks), operator
    call all_params()
    call build_lmto_alloy(trim(inp%lattice), inp%lsuper, inp%rcl, inp%sigma, pps(1, 1), pps(2, 1), inp%x, g, coord, mdl, inp%rinv, alpha)
    if (ordered) then
      if (trim(inp%lattice) /= 'bcc') error stop 'alloy_scf: ordered mode is B2 on lattice bcc'
      allocate (isa(g%n)); isa = modulo(coord(1, :), 2) == 0         ! A on the even (corner) sublattice, B on the odd (body-centre) one
      isite_b = 0
      do i = 1, g%n
        if (modulo(coord(1, i), 2) == 1) then; isite_b = i; exit; end if
      end do
      nd = 9 * g%n; nnodes = g%n; flip1 = 0; md1 = 0.0_dp; sflip = 0.0_dp
      allocate (oa(nlev, 4, nsp, 2), ob(nlev, 4, nsp, 2), ta(nlev), tb(nlev), v0o(nd), phimad(2, 2))
      oa = 0.0_dp; ob = 0.0_dp
      bvec3 = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
      tau2 = 0.0_dp; tau2(:, 2) = 0.5_dp
      call madelung_phi(bvec3, tau2, inp%alat, phimad)
      if (ed) then                       ! ED: S, the site species and the neighbour table are all that is needed (spin independent)
        call set_alloy_params(mdl, pps(1, 1), pps(2, 1), alpha)
        call build_lmto_ordered(g, mdl, isa, rop); rop%order = 2
        call ed_work_init(edw, nd)
      end if
    else
    call graph_distances(g, 1, dist)
    allocate (flippable(g%n)); flippable = dist <= inp%flip_radius
    if (inp%symmetry) then
      call d4h_ops(coord, 2 * inp%lsuper, perm, mats)
      call orbital_rep(mats, dmat)
      if (asr_ed) then                                  ! one operator per central species, start site frozen (not flippable)
        flippable(1) = .false.
        do i = 1, 2
          call build_augsym_operator(g, mdl, perm, dmat, 1, ngen, sopc(i), inp%maxrank, flippable, order=2)
          call augsym_freeze_site(sopc(i), 1, merge(1.0_dp, 0.0_dp, i == 1))
        end do
        nd = sopc(1)%dim(); nnodes = sopc(1)%nn; flip1 = 0; md1 = 0.0_dp; sflip = 0.0_dp
      else
        call build_augsym_operator(g, mdl, perm, dmat, 1, ngen, sop, inp%maxrank, flippable, order=merge(inp%order, 2, .not. ed))
        nd = sop%dim(); nnodes = sop%nn; flip1 = sop%flipr(1); md1 = sop%md(1); sflip = sop%s
      end if
    else
      call build_augorb_operator(g, mdl, 1, ngen, op, inp%maxrank, flippable, order=merge(inp%order, 2, .not. ed))
      nd = op%dim(); nnodes = op%nn; flip1 = op%flip(1); md1 = op%md(1); sflip = op%s
    end if
    end if
    res%nodes = nnodes
    if (vb) print '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,l1,a,l1,a,l1)', '  alloy ASR: nodes=', nnodes, '  levels=', ngen, '  recursion steps=', nlev, '  rank<=', inp%maxrank, '  flip radius=', inp%flip_radius, '  nspin=', nsp, '  sym=', inp%symmetry, '  energy-dependent=', inp%energy_dependent, '  ordered=', ordered
    allocate (offx(3)); offx(1) = 0
    do s = 1, 2; offx(s + 1) = offx(s) + (cr(s)%mesh%n + 3) * nsp; end do
    nx = offx(3)
    allocate (xin(nx), xout(nx), xnext(nx))
    call anderson_init(mixer, nx, inp%mix_history)
    allocate (zc(inp%ncont), wc(inp%ncont))
    if (ordered) then
      if (ed) ef = maxval(enu(2, 1:nsp, :)) + 0.35_dp
      allocate (U(1, 1, 1, 1), acoef(nlev, 4, nsp), bcoef(nlev, 4, nsp), ccoef(nlev, 4, nsp), v0(1), vv(1))
      acoef = 0.0_dp; bcoef = 0.0_dp; ccoef = 0.0_dp; mm = 0
    else if (ed) then
      allocate (U(1, 1, 1, 1), acoef(nlev, 4, nsp), bcoef(nlev, 4, nsp), ccoef(nlev, 4, nsp), v0(1), vv(1))
      acoef = 0.0_dp; bcoef = 0.0_dp; ccoef = 0.0_dp; mm = 0
      call ed_work_init(edw, nd)
      ef = maxval(enu(2, 1:nsp, :)) + 0.35_dp
    else
      allocate (U(nd, nlev, 4, nsp), acoef(nlev, 4, nsp), bcoef(nlev, 4, nsp), ccoef(nlev, 4, nsp), v0(nd), vv(nd))
    end if

    order_now = merge(2, min(inp%order, 2), ed .and. .not. inp%ed_third_order); dv_last = huge(1.0_dp); cnrm = 1.0_dp; onrm = 1.0_dp
    ed_depth_tol = 0.0_dp; if (ed) ed_depth_tol = inp%depth_tol
    ed_tail = inp%ed_tail
    allocate (hvec(nd), rhoex(nmax, 2, 2)); rhoex = 0.0_dp
    exact_now = inp%ed_exact
    do it = 1, inp%maxit
      ed_nchain = 0; ed_nlev = 0
      if (inp%order == 3 .and. order_now == 2 .and. dv_last < inp%order_switch) then
        order_now = 3
        call anderson_init(mixer, nx, inp%mix_history)          ! the map changes: drop the mixing history
        if (vb) print '(a,i0)', '  switching to the third-order Hamiltonian at iteration ', it
      end if
      sphq = (order_now == 3 .and. .not. ed .and. trim(inp%charges) == 'sphere')   ! sphere-charge moments: O^-1/2 |L> ~ (1 - hph/2)|L>
      do s = 1, 2
        n = cr(s)%mesh%n
        do sg = 1, nsp
          i = offx(s) + (sg - 1) * (n + 3)
          xin(i + 2:i + n) = cr(s)%v(2:n, sg) + 2.0_dp * inp%spec(s)%z / cr(s)%mesh%r(2:n)
          xin(i + 1) = xin(i + 2)
          xin(i + n + 1:i + n + 3) = enu(:, sg, s)
        end do
      end do
      call all_params()
      if (ordered .and. .not. ed) then
        do sg = 1, nsp
          call set_alloy_params(mdl, pps(1, sg), pps(2, sg), alpha)
          call build_lmto_ordered(g, mdl, isa, rop); rop%order = order_now
          do ip = 1, 2
            do oc = 1, 4
              v0o = 0.0_dp
              v0o(9 * (merge(1, isite_b, ip == 1) - 1) + iorb(oc)) = 1.0_dp
              if (sphq) then
                call lmto_hph(rop, v0o, hvec); v0o = v0o - 0.5_dp * hvec
              end if
              onrm(oc, sg, ip) = dot_product(v0o, v0o)
              call lanczos(rop, v0o, nlev, ta, tb, omm(oc, sg, ip))
              oa(:, oc, sg, ip) = ta; ob(:, oc, sg, ip) = tb
            end do
          end do
        end do
        lo = minval(enu(:, 1:nsp, :)) - 0.5_dp; hi = maxval(enu(:, 1:nsp, :)) + 1.5_dp
        do i = 1, 60
          ef = 0.5_dp * (lo + hi)
          call moments_ord(ef, mom)
          nn = electrons(mom)
          if (nn > nelec) then; hi = ef; else; lo = ef; end if
          if (hi - lo < 1.0e-11_dp) exit
        end do
        ef = 0.5_dp * (lo + hi)
        call moments_ord(ef, mom)
      else if (ed) then
        do sg = 1, nsp
          call ed_load(pps(1, sg), eda(sg)); call ed_load(pps(2, sg), edb(sg))
        end do
        if (vb .and. it == 1) then
          call ed_moments(6.0_dp, mom, dosef)
          print '(a,3f9.4,a,3f9.4)', '  sum rule (E_F=6 Ry)  A: M0(s,p,d)/(1,3,5) =', (mom(1, 1, 0, l) / real(2 * l + 1, dp), l=0, 2), '   B:', (mom(2, 1, 0, l) / real(2 * l + 1, dp), l=0, 2)
        end if
        if (inp%ed_scan .and. it == 1) then
          do jf = 0, 40
            e = -0.35_dp + 0.01_dp * jf
            call ed_moments(e, mom, dosef)
            nn = fs * (xs(1) * sum(mom(1, 1:nsp, 0, :)) + xs(2) * sum(mom(2, 1:nsp, 0, :)))
            print '(a,f8.4,a,f10.5,a,f9.4)', '  scan E=', e, '  N-nelec=', nn - nelec, '  DOS=', fs * (xs(1) * sum(dosef(1, 1:nsp, :)) + xs(2) * sum(dosef(2, 1:nsp, :)))
          end do
          stop 0
        end if
        if (exact_now) call choose_ebot()
        nfe = 12                       ! contour passes; a converged SCF needs 1-2.  With 6, wild iterations ran out and the linear
        ef_old = 0.0_dp; nn_old = 0.0_dp ! correction left N /= N_elec (species charges unbalanced by up to 0.3 e)
        do jf = 1, nfe
          call ed_moments(ef, mom, dosef, rhoex)
          nn = fs * (xs(1) * sum(mom(1, 1:nsp, 0, :)) + xs(2) * sum(mom(2, 1:nsp, 0, :)))
          dvc = fs * (xs(1) * sum(dosef(1, 1:nsp, :)) + xs(2) * sum(dosef(2, 1:nsp, :)))
          if (vb .and. it == 1) print '(a,i2,a,f10.6,a,f9.5,a,f9.4)', '    E_F search ', jf, ': E_F=', ef, '  N-nelec=', nn - nelec, '  DOS(E_F)=', dvc
          if (abs(nn - nelec) < 2.0e-4_dp) exit
          slope = max(dvc, 1.0_dp)
          if (jf > 1) then
            if ((nn - nn_old) * (ef - ef_old) > 0.0_dp .and. abs(ef - ef_old) > 1.0e-9_dp) slope = (nn - nn_old) / (ef - ef_old)
          end if
          e = max(-0.15_dp, min(0.15_dp, (nelec - nn) / slope))
          ef_old = ef; nn_old = nn
          if (jf == nfe) then          ! out of evaluations: first-order correction of the moments with the DOS at E_F
            e = max(-0.15_dp, min(0.15_dp, (nelec - nn) / max(dvc, 0.1_dp)))
            if (vb) print '(a,i0,a,es10.3,a,f9.5,a,es10.3)', '    warning: E_F search not converged after ', nfe, ' passes: N - N_elec =', &
              nn - nelec, ';  first-order E_F shift', e, ', remaining', nn + e * dvc - nelec
            do s = 1, 2; do sg = 1, nsp; do l = 0, 2; do k = 0, 2
              mom(s, sg, k, l) = mom(s, sg, k, l) + e * (ef + 0.5_dp * e - enu(l, sg, s))**k * dosef(s, sg, l)
            end do; end do; end do; end do
          end if
          ef = ef + e
        end do
        if (inp%ed_exact .and. jf > nfe) call ed_moments(ef, mom, dosef, rhoex)   ! density at the final E_F
      else
      ! recursions: four start orbitals per spin
      do sg = 1, nsp
        call set_alloy_params(mdl, pps(1, sg), pps(2, sg), alpha)
        if (inp%symmetry) then; call augsym_set_params(sop, mdl); sop%order = order_now
        else; call augorb_set_params(op, mdl); op%order = order_now; end if
        do oc = 1, 4
          v0 = 0.0_dp; v0(iorb(oc)) = 1.0_dp
          if (inp%symmetry) call augsym_set_start(sop, iorb(oc), dmat)
          if (sphq) then
            if (inp%symmetry) then; call augsym_hph(sop, v0, hvec); else; call augorb_hph(op, v0, hvec); end if
            v0 = v0 - 0.5_dp * hvec
          end if
          cnrm(oc, sg) = dot_product(v0, v0)
          if (inp%symmetry) then
            call lanczos_basis(sop, v0, nlev, acoef(:, oc, sg), bcoef(:, oc, sg), ub, mm(oc, sg))
          else
            call lanczos_basis(op, v0, nlev, acoef(:, oc, sg), bcoef(:, oc, sg), ub, mm(oc, sg))
          end if
          U(:, :, oc, sg) = ub
          vv = 0.0_dp
          vv(iorb(oc)) = md1
          if (flip1 > 0) vv(9 * (flip1 - 1) + iorb(oc)) = vv(9 * (flip1 - 1) + iorb(oc)) + sflip
          if (sphq) then                                  ! species part of O^-1/2 |L>:  (1 - hph/2) N_R |L>
            if (inp%symmetry) then; call augsym_hph(sop, vv, hvec); else; call augorb_hph(op, vv, hvec); end if
            vv = vv - 0.5_dp * hvec
          end if
          do k = 1, mm(oc, sg); ccoef(k, oc, sg) = dot_product(U(:, k, oc, sg), vv); end do
        end do
      end do
      if (vb .and. it == 1) then          ! sum rule: with E_F far above the bands every orbital carries unit weight
        call moments(6.0_dp, mom)
        print '(a,3f9.4,a,3f9.4)', '  sum rule (E_F=6 Ry)  A: M0(s,p,d)/(1,3,5) =', (mom(1, 1, 0, l) / real(2 * l + 1, dp), l=0, 2), '   B:', (mom(2, 1, 0, l) / real(2 * l + 1, dp), l=0, 2)
      end if
      ! Fermi level from the electron count
      lo = minval(enu(:, 1:nsp, :)) - 0.5_dp; hi = maxval(enu(:, 1:nsp, :)) + 1.5_dp
      do i = 1, 60
        ef = 0.5_dp * (lo + hi)
        call moments(ef, mom)
        nn = electrons(mom)
        if (nn > nelec) then; hi = ef; else; lo = ef; end if
        if (hi - lo < 1.0e-11_dp) exit
      end do
      ef = 0.5_dp * (lo + hi)
      call moments(ef, mom)
      end if
      ! densities, core, electrostatics, new potentials
      dv = 0.0_dp
      do s = 1, 2
        n = cr(s)%mesh%n
        allocate (nval(n, nsp), ncore(n, nsp), ntot(n), vh(n), vxc(n, 2), rho(n, 2), exc(n), tmp(n), cum(n), yc(2, 0:0, n))
        nval = 0.0_dp; ncore = 0.0_dp
        do sg = 1, nsp
          do l = 0, 2
            nval(:, sg) = nval(:, sg) + fs * ( &
                mom(s, sg, 0, l) * (y(1, 0, 1:n, l, s, sg)**2 + y(2, 0, 1:n, l, s, sg)**2 / c**2) &
              + 2.0_dp * mom(s, sg, 1, l) * (y(1, 0, 1:n, l, s, sg) * y(1, 1, 1:n, l, s, sg) + y(2, 0, 1:n, l, s, sg) * y(2, 1, 1:n, l, s, sg) / c**2) &
              + mom(s, sg, 2, l) * (y(1, 1, 1:n, l, s, sg)**2 + y(2, 1, 1:n, l, s, sg)**2 / c**2))
            if (sphq) then                                   ! sphere-charge moments m_k (as the k-space ASA): n = m0 phi^2 + 2 m1 phi phidot + m2 phidot^2
            else if (trim(inp%density) == 'stuttgart') then     ! q0 = m0 + p m2 (LMTO47): n = (q0 - p q2) phi^2 + 2 q1 phi phidot + q2 phidot^2
              nval(:, sg) = nval(:, sg) - fs * mom(s, sg, 2, l) * pps(s, sg)%p(l) &
                                          * (y(1, 0, 1:n, l, s, sg)**2 + y(2, 0, 1:n, l, s, sg)**2 / c**2)
            else                                            ! normalised phi(E) expanded to second order
              nval(:, sg) = nval(:, sg) + fs * mom(s, sg, 2, l) &
                  * (y(1, 0, 1:n, l, s, sg) * y(1, 2, 1:n, l, s, sg) + y(2, 0, 1:n, l, s, sg) * y(2, 2, 1:n, l, s, sg) / c**2)
            end if
          end do
          if (ed .and. inp%ed_exact .and. inp%ed_exact_density) nval(:, sg) = fs * rhoex(1:n, sg, s)    ! exact KKR-ASA density -1/pi Im int G phi(z)^2
          do is = 1, size(inp%spec(s)%core)
            nod = inp%spec(s)%core(is)%n - inp%spec(s)%core(is)%l - 1
            call bound_state(cr(s)%mesh, cr(s)%v(:, sg), inp%spec(s)%core(is)%l, nod, inp%spec(s)%z, c, n, &
                             -1.6_dp * inp%spec(s)%z**2 - 20.0_dp, cr(s)%v(n, sg) - 1.0e-6_dp, e, yc, core_ok)
            ncore(:, sg) = ncore(:, sg) + inp%spec(s)%core(is)%occ / real(nsp, dp) * (yc(1, 0, :)**2 + yc(2, 0, :)**2 / c**2)
          end do
        end do
        cr(s)%nrad = nval + ncore
        ntot = sum(cr(s)%nrad, dim=2)
        tmp = ntot; call cumint(cr(s)%mesh, tmp, cum)
        cr(s)%qnet = cum(n) - inp%spec(s)%z
        if (nsp == 2) then
          tmp = cr(s)%nrad(:, 1) - cr(s)%nrad(:, 2); call cumint(cr(s)%mesh, tmp, cum); cr(s)%moment = cum(n)
        end if
        deallocate (nval, ncore, ntot, vh, vxc, rho, exc, tmp, cum, yc)
      end do
      if (ordered) then                    ! exact Madelung shift of the two-site cell, V_i = 2 sum_j Phi_ij qnet_j
        qat2 = [cr(1)%qnet, cr(2)%qnet]
        vmo = 2.0_dp * matmul(phimad, qat2)
      end if
      do s = 1, 2
        n = cr(s)%mesh%n
        allocate (ntot(n), vh(n), vxc(n, 2), rho(n, 2), exc(n))
        ntot = sum(cr(s)%nrad, dim=2)
        call poisson_hartree(cr(s)%mesh, ntot, vh)
        do k = 2, n; rho(k, :) = cr(s)%nrad(k, :) / (4.0_dp * pi * cr(s)%mesh%r(k)**2); end do
        rho(1, :) = rho(2, :)
        if (nsp == 1) then
          call xc_eval(xf, rho(:, 1), exc, vxc(:, 1), cr(s)%mesh%r); vxc(:, 2) = vxc(:, 1)
        else
          call xc_eval_pol(xf, rho(:, 1), rho(:, 2), exc, vxc(:, 1), vxc(:, 2), cr(s)%mesh%r)
        end if
        vm = -inp%sim_beta * 2.0_dp * cr(s)%qnet / r1                               ! SIM
        if (ordered) vm = vmo(s)
        do sg = 1, nsp
          i = offx(s) + (sg - 1) * (n + 3)
          xout(i + 2:i + n) = vh(2:n) + vxc(2:n, sg) + vm
          xout(i + 1) = xout(i + 2)
          dvc = maxval(abs(xout(i + 2:i + n) - xin(i + 2:i + n)) * cr(s)%mesh%r(2:n) / (1.0_dp + cr(s)%mesh%r(2:n)))
          dv = max(dv, dvc)
          do l = 0, 2
            xout(i + n + 1 + l) = enu(l, sg, s) + 0.7_dp * mom(s, sg, 1, l) / max(mom(s, sg, 0, l), 1.0e-6_dp)
          end do
        end do
        deallocate (ntot, vh, vxc, rho, exc)
      end do
      if (.not. allocated(wmix)) then           ! Anderson metric = the dV measure (r/(1+r))^2 on the potentials, 1 on E_nu (see asr_asa)
        allocate (wmix(nx)); wmix = 1.0_dp
        do s = 1, 2
          n = cr(s)%mesh%n
          do sg = 1, nsp
            i = offx(s) + (sg - 1) * (n + 3)
            wmix(i + 1:i + n) = (cr(s)%mesh%r(1:n) / (1.0_dp + cr(s)%mesh%r(1:n)))**2
          end do
        end do
      end if
      call anderson_step(mixer, xin, xout, inp%beta, xnext, wmix)
      do s = 1, 2
        n = cr(s)%mesh%n
        do sg = 1, nsp
          i = offx(s) + (sg - 1) * (n + 3)
          do k = 1, n
            xnext(i + k) = xin(i + k) + max(-0.15_dp, min(0.15_dp, xnext(i + k) - xin(i + k)))
          end do
          cr(s)%v(2:n, sg) = xnext(i + 2:i + n) - 2.0_dp * inp%spec(s)%z / cr(s)%mesh%r(2:n)
          cr(s)%v(1, sg) = cr(s)%v(2, sg)
          do l = 0, 2
            enu(l, sg, s) = xin(i + n + 1 + l) + max(-0.08_dp, min(0.08_dp, xnext(i + n + 1 + l) - xin(i + n + 1 + l)))
            if (clamp_enu(cr(s)%mesh, cr(s)%v(:, sg), l, inp%spec(s)%nvl(l) - l - 1, inp%spec(s)%z, c, enu(l, sg, s))) then
              if (vb) print '(a,i0,a,i0,a,i0,a,f9.5)', '    E_nu outside its band window: species ', s, ' spin ', sg, ' l ', l, ' -> ', enu(l, sg, s)
            end if
          end do
        end do
      end do
      res%iter = it; res%dv = dv; res%ef = ef; dv_last = dv
      if (.not. allocated(res%hist_dv)) then
        allocate (res%hist_dv(inp%maxit), res%hist_ef(inp%maxit), res%hist_moment(inp%maxit, 2), res%hist_qnet(inp%maxit, 2))
      end if
      res%hist_dv(it) = dv; res%hist_ef(it) = ef
      res%hist_moment(it, :) = [cr(1)%moment, cr(2)%moment]; res%hist_qnet(it, :) = [cr(1)%qnet, cr(2)%qnet]
      if (vb) print '(a,i3,a,es10.3,a,f10.6,a,2f9.5,a,2f8.4)', 'it', it, '  dV=', dv, '  EF=', ef, '  qnet(A,B)=', cr(1)%qnet, cr(2)%qnet, '  m(A,B)=', cr(1)%moment, cr(2)%moment
      if (vb .and. ed .and. ed_nchain > 0) print '(a,f6.1,a,i0,a)', '      mean recursion depth used ', real(ed_nlev, dp) / real(ed_nchain, dp), ' of ', nlev, ' levels'
      if (dv < inp%tol .and. it > 3 .and. order_now == inp%order) then; res%converged = .true.; exit; end if
    end do
    call xc_finalize(xf)
    if (len_trim(inp%pot_out) > 0) call write_potentials()
    do s = 1, 2
      do sg = 1, nsp
        do l = 0, 2
          call potpar_at_enu(cr(s)%mesh, cr(s)%v(:, sg), l, enu(l, sg, s), inp%spec(s)%z, c, avw, cr(s)%pp(l, sg))
          cr(s)%q(0, l, sg) = mom(s, sg, 0, l); cr(s)%q(1, l, sg) = mom(s, sg, 1, l); cr(s)%q(2, l, sg) = mom(s, sg, 2, l)
        end do
        cr(s)%enu(:, sg) = enu(:, sg, s)
      end do
    end do
    res%x = xs(1); res%nspin = nsp
    res%ordered = ordered
    if (ordered) then
      res%oa = oa; res%ob = ob; res%omm = omm
    end if
    res%acoef = acoef; res%bcoef = bcoef; res%ccoef = ccoef; res%mm = mm
    res%ed = ed
    res%cnrm = cnrm; res%onrm = onrm
    res%terminator = inp%terminator
    if (ed .and. inp%ndos > 0) call ed_dos_final()
    if (ed .and. ordered .and. inp%kkr_check > 0 .and. exact_now) call kkr_check(inp%kkr_check)
    res%cls = cr
    res%moment_cell = xs(1) * cr(1)%moment + xs(2) * cr(2)%moment
    if (ordered) res%moment_cell = cr(1)%moment + cr(2)%moment      ! per two-atom cell
  contains
    !> species potentials (with the nucleus) and E_nu, both spins: file format of read_potentials
    subroutine write_potentials()
      integer :: u, ss, gg, kk
      open (newunit=u, file=trim(inp%pot_out), status='replace', action='write')
      write (u, '(a)') '# asr alloy species potentials (Ry): nspin; per species: z n r(n); E_nu(0:2) per spin; v(r_k, spin), k = 1..n'
      write (u, '(i3)') nsp
      do ss = 1, 2
        write (u, '(f8.2,i8,es25.16)') inp%spec(ss)%z, cr(ss)%mesh%n, cr(ss)%mesh%r(cr(ss)%mesh%n)
        do gg = 1, nsp; write (u, '(3es25.16)') enu(:, gg, ss); end do
        do kk = 1, cr(ss)%mesh%n; write (u, '(2es25.16)') cr(ss)%v(kk, :); end do
      end do
      close (u)
    end subroutine write_potentials

    !> replace the atomic starts by a written potential; nspin 1 -> 2 copies the potential to both spins (no splitting);
    !> pot_flip(s) swaps the spin channels of species s
    subroutine read_potentials()
      integer :: u, ss, gg, kk, nf, nspf, ios
      real(dp) :: zf, rf, vv2(2), ee(0:2, 2)
      open (newunit=u, file=trim(inp%pot_in), status='old', action='read', iostat=ios)
      if (ios /= 0) error stop 'alloy_scf: cannot open pot_in'
      read (u, *); read (u, *) nspf
      if (nspf > nsp) error stop 'alloy_scf: pot_in is spin polarised, the run is not'
      do ss = 1, 2
        read (u, *) zf, nf, rf
        if (nf /= cr(ss)%mesh%n .or. abs(zf - inp%spec(ss)%z) > 1.0e-6_dp .or. abs(rf - cr(ss)%mesh%r(nf)) > 1.0e-8_dp * rf) &
          error stop 'alloy_scf: pot_in species / radial mesh differ from this input'
        do gg = 1, nspf; read (u, *) ee(:, gg); end do
        if (nspf == 1) ee(:, 2) = ee(:, 1)
        do kk = 1, nf
          read (u, *) vv2(1:nspf)
          if (nspf == 1) vv2(2) = vv2(1)
          cr(ss)%v(kk, :) = vv2(1:nsp)
        end do
        enu(:, 1:nsp, ss) = ee(:, 1:nsp)
        if (nsp == 2 .and. inp%pot_flip(ss)) then
          cr(ss)%v(:, [1, 2]) = cr(ss)%v(:, [2, 1]); enu(:, [1, 2], ss) = enu(:, [2, 1], ss)
        end if
      end do
      close (u)
      if (vb) print '(3a,i0,a,2l2)', '  start potentials from ', trim(inp%pot_in), ' (nspin ', nspf, ')  spin flip(A,B):', inp%pot_flip
    end subroutine read_potentials

    !> electrons per cell from the occupied moments: sum of m0 (+ p m2 for sphere-charge moments, order 3)
    real(dp) function electrons(mo)
      real(dp), intent(in) :: mo(2, 2, 0:2, 0:2)
      integer :: ss, gg, ll
      real(dp) :: q(2)
      q = 0.0_dp
      do ss = 1, 2
        do gg = 1, nsp
          do ll = 0, 2
            q(ss) = q(ss) + mo(ss, gg, 0, ll)
            if (sphq) q(ss) = q(ss) + pps(ss, gg)%p(ll) * mo(ss, gg, 2, ll)
          end do
        end do
      end do
      electrons = fs * (xs(1) * q(1) + xs(2) * q(2))
    end function electrons

    subroutine ed_load(pp, sp)
      type(potpar_t), intent(in) :: pp
      type(ed_params_t), intent(out) :: sp
      integer :: ii
      do ii = 1, 9
        sp%c(ii) = pp%c(lorb9(ii)); sp%del(ii) = pp%sqrdel(lorb9(ii))**2; sp%gam(ii) = pp%gamma(lorb9(ii))
        sp%enu(ii) = pp%enu(lorb9(ii))
        if (order_now == 3 .and. inp%ed_third_order) sp%p(ii) = pp%p(lorb9(ii))   ! third-order P(z): diverges in SCF (STATUS 2026-10-07), off by default
      end do
    end subroutine ed_load

    !> occupied moments from the energy-dependent recursion on the contour to efv, DOS at E_F (eta = ed_eta) and, with
    !> exact potential functions (inp%ed_exact), the valence radial density rhov(r, spin, species) = -1/pi Im int G phi^2.
    subroutine ed_moments(efv, mo, dosv, rhov)
      real(dp), intent(in) :: efv
      real(dp), intent(out) :: mo(2, 2, 0:2, 0:2), dosv(2, 2, 0:2)
      real(dp), intent(out), optional :: rhov(:, :, :)
      complex(dp), allocatable :: zq(:), wq(:), ph(:, :, :), ts(:, :, :)
      complex(dp) :: gs(2), z
      real(dp) :: ebot
      integer :: j, o, gg, kk, ll, nq, pp, site, ss
      logical :: dens
      nq = inp%ncont + nvert
      allocate (zq(nq + 1), wq(nq + 1))
      ebot = minval(enu(:, 1:nsp, :)) - 1.2_dp
      if (exact_now) ebot = ebot_exact                    ! chosen once per iteration away from the poles of the exact P
      call contour_nodes_closed(ebot, efv, inp%ed_eta, inp%ncont, nvert, zq(1:nq), wq(1:nq))   ! arc to E_F + i eta, closed by the vertical
      zq(nq + 1) = cmplx(efv, inp%ed_eta, dp); wq(nq + 1) = (0.0_dp, 0.0_dp)                  ! last node: DOS at E_F only
      dens = present(rhov) .and. exact_now .and. inp%ed_exact_density
      if (dens) then; allocate (ph(nmax, 0:2, 2), ts(nmax, 0:2, 2)); rhov = 0.0_dp; end if
      mo = 0.0_dp; dosv = 0.0_dp
      do gg = 1, nsp
        if (asr_ed) call asr_ed_prepare(gg)
        do j = 1, nq + 1
          z = zq(j)
          if (exact_now) call set_exact(z, gg, dens .and. j <= nq, ph, ts)
          do o = 1, 4
            ll = lorb(o)
            if (asr_ed) then
              do pp = 1, 2
                call augsym_set_start(sopc(pp), iorb(o), dmat)
                call asr_ed_green(sopc(pp), edw, edp_c(pp), z, iorb(o), nlev, gs(pp))
              end do
            else if (ordered) then
              do pp = 1, 2
                site = merge(1, isite_b, pp == 1)
                call ed_ordered_green(rop, edw, z, eda(gg), edb(gg), alpha9, site, iorb(o), nlev, gs(pp))
              end do
            else
              call augsym_set_start(sop, iorb(o), dmat)
              call ed_species_green(sop, edw, z, eda(gg), edb(gg), alpha9, inp%x, iorb(o), nlev, flip1, md1, sflip, gs(1), gs(2))
            end if
            if (j == nq + 1) then
              dosv(:, gg, ll) = dosv(:, gg, ll) - worb(o) * aimag(gs) / pi
              cycle
            end if
            do ss = 1, 2
              do kk = 0, 2
                mo(ss, gg, kk, ll) = mo(ss, gg, kk, ll) - worb(o) * aimag(wq(j) * (z - enu(ll, gg, ss))**kk * gs(ss)) / pi
              end do
              if (dens) rhov(1:cr(ss)%mesh%n, gg, ss) = rhov(1:cr(ss)%mesh%n, gg, ss) &
                  - worb(o) * aimag(wq(j) * (gs(ss) * ph(1:cr(ss)%mesh%n, ll, ss) + ts(1:cr(ss)%mesh%n, ll, ss))) / pi
            end do
          end do
        end do
      end do
    end subroutine ed_moments

    !> Same-model reference for the ordered ED route at the final potentials: k-space KKR-ASA, G_L(z) = lambda + mu^2
    !> <[P(z) - S(k)]^-1>_k with the same exact P, S^alpha (clusters of radius rcl, same alpha) and contour.  Prints E_F and the
    !> per-(species, spin, l) occupations from the BZ sum and from the recursion at the same E_F.
    subroutine kkr_check(nkk)
      integer, intent(in) :: nkk
      type(cell_t) :: kc
      type(cluster_t) :: kcl(2)
      complex(dp), allocatable :: sk(:, :, :), zq(:), wq(:), gk(:, :, :)
      real(dp) :: als(0:2, 2), bv(3, 3), kv(3), lo, hi, efk, nk3, qk(0:2, 2, 2), mrc(2, 2, 0:2, 0:2), dd(2, 2, 0:2), nnk
      integer :: ik, ix, iy, iz, ib, i, a, nq, it2, ss, gg, ll
      complex(dp), parameter :: ci = (0.0_dp, 1.0_dp)
      kc%nat = 2; kc%avec = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
      allocate (kc%tau(3, 2), kc%isp(2)); kc%tau = 0.0_dp; kc%tau(:, 2) = 0.5_dp; kc%isp = [1, 2]
      call setup_cell(kc)
      als(:, 1) = alpha; als(:, 2) = alpha
      do ib = 1, 2; call build_cluster(kc, ib, inp%rcl, als, kcl(ib)); end do
      bv = recip_vectors(kc%avec)
      nk3 = real(nkk, dp)**3
      allocate (sk(18, 18, nkk**3))
      do iz = 0, nkk - 1; do iy = 0, nkk - 1; do ix = 0, nkk - 1
        ik = 1 + ix + nkk * (iy + nkk * iz)
        kv = (ix + 0.5_dp) / nkk * bv(:, 1) + (iy + 0.5_dp) / nkk * bv(:, 2) + (iz + 0.5_dp) / nkk * bv(:, 3)
        sk(:, :, ik) = (0.0_dp, 0.0_dp)
        do ib = 1, 2
          do i = 1, kcl(ib)%n
            a = kcl(ib)%iat(i)
            sk(9 * (ib - 1) + 1:9 * ib, 9 * (a - 1) + 1:9 * a, ik) = sk(9 * (ib - 1) + 1:9 * ib, 9 * (a - 1) + 1:9 * a, ik) &
              + transpose(kcl(ib)%s(:, :, i)) * exp(ci * 2.0_dp * pi * dot_product(kv, kcl(ib)%d(:, i)))
          end do
        end do
        sk(:, :, ik) = 0.5_dp * (sk(:, :, ik) + conjg(transpose(sk(:, :, ik))))
      end do; end do; end do
      nq = inp%ncont + nvert
      allocate (zq(nq), wq(nq), gk(18, nq, 2))
      lo = res%ef - 0.1_dp; hi = res%ef + 0.1_dp
      do it2 = 1, 40
        efk = 0.5_dp * (lo + hi)
        call kocc(efk, qk, sk, nkk, nq, zq, wq)
        nnk = fs * sum(qk)
        if (nnk > 2.0_dp * nelec) then; hi = efk; else; lo = efk; end if
        if (hi - lo < 2.0e-6_dp) exit
      end do
      efk = 0.5_dp * (lo + hi)
      call kocc(efk, qk, sk, nkk, nq, zq, wq)
      call ed_moments(efk, mrc, dd)
      print '(a,f10.6,a,f10.6,a,i0,a)', '  KKR k-space check: E_F(k-space) = ', efk, '   E_F(recursion) = ', res%ef, '   (nk = ', nkk, ')'
      do ss = 1, 2
        do gg = 1, nsp
          print '(a,i1,a,i1,a,3f9.5,a,3f9.5)', '    species ', ss, ' spin ', gg, '  occupations s p d: k-space', qk(:, gg, ss), &
              '   recursion', mrc(ss, gg, 0, :)
        end do
      end do
    end subroutine kkr_check

    !> k-space KKR-ASA occupied moments m0(l, spin, species) to efv (used by kkr_check)
      subroutine kocc(efv, q, sk, nkk, nq, zq, wq)
        real(dp), intent(in) :: efv
        real(dp), intent(out) :: q(0:2, 2, 2)
        integer, intent(in) :: nkk, nq
        complex(dp), intent(in) :: sk(:, :, :)
        complex(dp), intent(inout) :: zq(nq), wq(nq)
        real(dp) :: nk3
        integer :: gg, ss
        complex(dp) :: pdiag(18), mu18(18), lm18(18), mm(18, 18), gsum(18), wk(18 * 64)
        integer :: j, kk2, ipv(18), info, ll2
        complex(dp), allocatable :: dummy(:, :, :)
        nk3 = real(nkk, dp)**3
        call contour_nodes_closed(ebot_exact, efv, inp%ed_eta, inp%ncont, nvert, zq, wq)
        allocate (dummy(1, 0:2, 2))
        q = 0.0_dp
        do gg = 1, nsp
          do j = 1, nq
            call set_exact(zq(j), gg, .false., dummy)
            pdiag(1:9) = eda(gg)%pfx; pdiag(10:18) = edb(gg)%pfx
            mu18(1:9) = eda(gg)%mux; mu18(10:18) = edb(gg)%mux
            lm18(1:9) = eda(gg)%lmx; lm18(10:18) = edb(gg)%lmx
            gsum = 0.0_dp
            !$omp parallel do private(mm, ipv, info, wk, ll2) reduction(+:gsum)
            do kk2 = 1, nkk**3
              mm = -sk(:, :, kk2)
              do ll2 = 1, 18; mm(ll2, ll2) = mm(ll2, ll2) + pdiag(ll2); end do
              call zgetrf(18, 18, mm, 18, ipv, info)
              call zgetri(18, mm, 18, ipv, wk, size(wk), info)
              do ll2 = 1, 18; gsum(ll2) = gsum(ll2) + mm(ll2, ll2); end do
            end do
            !$omp end parallel do
            gsum = lm18 + mu18 * gsum / nk3
            do ll2 = 1, 18
              ss = (ll2 - 1) / 9 + 1
              q(lorb9(ll2 - 9 * (ss - 1)), gg, ss) = q(lorb9(ll2 - 9 * (ss - 1)), gg, ss) - aimag(wq(j) * gsum(ll2)) / pi
            end do
          end do
        end do
      end subroutine kocc

    !> bottom of the exact-P contour: on the real axis between min C_s - 0.6 and min C_s - 0.15 Ry (below the valence band),
    !> at the point farthest from the poles of every P^alpha_l (smallest max |lambda| = |C'/C| over species, l, spin)
    subroutine choose_ebot()
      complex(dp) :: pf1, mu1, lm1
      real(dp) :: lo, e0, best, worst
      integer :: k, ss, gg, ll
      lo = minval([(minval(pps(ss, 1:nsp)%c(0)), ss=1, 2)]) - 0.6_dp
      best = huge(1.0_dp); ebot_exact = lo
      do k = 0, 45
        e0 = lo + 0.01_dp * k
        worst = 0.0_dp
        do ss = 1, 2
          do gg = 1, nsp
            do ll = 0, 2
              call pfun_exact(cr(ss)%mesh, cr(ss)%v(:, gg), ll, cmplx(e0, 0.0_dp, dp), inp%spec(ss)%z, c, avw, alpha(ll), pf1, mu1, lm1)
              worst = max(worst, abs(lm1))
            end do
          end do
        end do
        if (worst < best) then; best = worst; ebot_exact = e0; end if
      end do
    end subroutine choose_ebot

    !> exact P^alpha, mu^2, lambda (and phi^2 if wanted) of both species for spin gg at energy z
    subroutine set_exact(z, gg, want_phi, ph, ts)
      complex(dp), intent(in) :: z
      integer, intent(in) :: gg
      logical, intent(in) :: want_phi
      complex(dp), intent(inout) :: ph(:, 0:, :)
      complex(dp), intent(inout), optional :: ts(:, 0:, :)
      complex(dp) :: pf3(0:2), mu3(0:2), lm3(0:2)
      integer :: ss, ll, ii, nn
      do ss = 1, 2
        nn = cr(ss)%mesh%n
        do ll = 0, 2
          if (want_phi) then
            call pfun_exact(cr(ss)%mesh, cr(ss)%v(:, gg), ll, z, inp%spec(ss)%z, c, avw, alpha(ll), pf3(ll), mu3(ll), lm3(ll), ph(1:nn, ll, ss), ts(1:nn, ll, ss))
          else
            call pfun_exact(cr(ss)%mesh, cr(ss)%v(:, gg), ll, z, inp%spec(ss)%z, c, avw, alpha(ll), pf3(ll), mu3(ll), lm3(ll))
          end if
        end do
        do ii = 1, 9
          if (ss == 1) then
            eda(gg)%exact = .true.; eda(gg)%pfx(ii) = pf3(lorb9(ii)); eda(gg)%mux(ii) = mu3(lorb9(ii)); eda(gg)%lmx(ii) = lm3(lorb9(ii))
          else
            edb(gg)%exact = .true.; edb(gg)%pfx(ii) = pf3(lorb9(ii)); edb(gg)%mux(ii) = mu3(lorb9(ii)); edb(gg)%lmx(ii) = lm3(lorb9(ii))
          end if
        end do
      end do
    end subroutine set_exact

    subroutine ed_dos_final()
      complex(dp) :: ga1, gb1, phd(1, 0:2, 2)
      integer :: ie, o, gg, ll, pp
      allocate (res%edos_e(inp%ndos), res%edos(inp%ndos, 0:2, 2, nsp)); res%edos = 0.0_dp
      do ie = 1, inp%ndos
        res%edos_e(ie) = inp%dos_emin + (inp%dos_emax - inp%dos_emin) * (ie - 1) / real(max(inp%ndos - 1, 1), dp)
      end do
      if (ordered) then
        do gg = 1, nsp
          do ip = 1, 2
            do o = 1, 4
              ll = lorb(o)
              do ie = 1, inp%ndos
                if (inp%ed_exact) call set_exact(cmplx(res%edos_e(ie), inp%dos_eta, dp), gg, .false., phd, phd)
                call ed_ordered_green(rop, edw, cmplx(res%edos_e(ie), inp%dos_eta, dp), eda(gg), edb(gg), alpha9, merge(1, isite_b, ip == 1), iorb(o), nlev, ga1)
                res%edos(ie, ll, ip, gg) = res%edos(ie, ll, ip, gg) - worb(o) * aimag(ga1) / pi
              end do
            end do
          end do
        end do
        return
      end if
      do gg = 1, nsp
        if (asr_ed) call asr_ed_prepare(gg)
        do o = 1, 4
          if (.not. asr_ed) call augsym_set_start(sop, iorb(o), dmat)
          if (asr_ed) then
            ll = lorb(o)
            do ie = 1, inp%ndos
              do pp = 1, 2
                call augsym_set_start(sopc(pp), iorb(o), dmat)
                call asr_ed_green(sopc(pp), edw, edp_c(pp), cmplx(res%edos_e(ie), inp%dos_eta, dp), iorb(o), nlev, ga1)
                res%edos(ie, ll, pp, gg) = res%edos(ie, ll, pp, gg) - worb(o) * aimag(ga1) / pi
              end do
            end do
            cycle
          end if
          ll = lorb(o)
          do ie = 1, inp%ndos
            if (inp%ed_exact) call set_exact(cmplx(res%edos_e(ie), inp%dos_eta, dp), gg, .false., phd, phd)
            call ed_species_green(sop, edw, cmplx(res%edos_e(ie), inp%dos_eta, dp), eda(gg), edb(gg), alpha9, inp%x, iorb(o), nlev, flip1, md1, sflip, ga1, gb1)
            res%edos(ie, ll, 1, gg) = res%edos(ie, ll, 1, gg) - worb(o) * aimag(ga1) / pi
            res%edos(ie, ll, 2, gg) = res%edos(ie, ll, 2, gg) - worb(o) * aimag(gb1) / pi
          end do
        end do
      end do
    end subroutine ed_dos_final

    !> AMscf form: operator parameters (central species A, B) and the energy-dependent diagonals for spin gg
    subroutine asr_ed_prepare(gg)
      integer, intent(in) :: gg
      integer :: cc
      do cc = 1, 2
        call asr_ed_setup(mdl, pps(1, gg), pps(2, gg), alpha, cc == 1, edp_c(cc))
        call augsym_set_params(sopc(cc), mdl)
      end do
    end subroutine asr_ed_prepare

    subroutine all_params()
      integer :: ss, gg, ll
      do ss = 1, 2
        do gg = 1, nsp
          do ll = 0, 2
            call potpar_at_enu(cr(ss)%mesh, cr(ss)%v(:, gg), ll, enu(ll, gg, ss), inp%spec(ss)%z, c, avw, ppl, &
                               y(:, :, 1:cr(ss)%mesh%n, ll, ss, gg))
            cr(ss)%pp(ll, gg) = ppl
            pps(ss, gg)%enu(ll) = ppl%enu; pps(ss, gg)%c(ll) = ppl%c; pps(ss, gg)%sqrdel(ll) = ppl%sqrdel
            pps(ss, gg)%p(ll) = ppl%p; pps(ss, gg)%gamma(ll) = ppl%gamma; pps(ss, gg)%alpha_rep(ll) = ppl%gamma
          end do
        end do
      end do
    end subroutine all_params

    !> ordered mode: sublattice/spin-resolved occupied moments from the real-space recursion started on an A and a B site
    subroutine moments_ord(efv, mo)
      real(dp), intent(in) :: efv
      real(dp), intent(out) :: mo(2, 2, 0:2, 0:2)
      real(dp) :: ebot, ainf, binf
      complex(dp), allocatable :: g1(:)
      integer :: j, o, kk, ll, gg, pp
      allocate (g1(inp%ncont))
      ebot = minval(enu(:, 1:nsp, :)) - 1.2_dp
      call contour_nodes(ebot, efv, inp%ncont, zc, wc)
      mo = 0.0_dp
      do gg = 1, nsp
        do pp = 1, 2
          do o = 1, 4
            call tail_params(inp%terminator, oa(:, o, gg, pp), ob(:, o, gg, pp), omm(o, gg, pp), ainf, binf)
            do j = 1, inp%ncont
              g1(j) = onrm(o, gg, pp) * cf_green(oa(:, o, gg, pp), ob(:, o, gg, pp), omm(o, gg, pp), zc(j), sqrt_tail(zc(j), ainf, binf))
            end do
            ll = lorb(o)
            do kk = 0, 2
              mo(pp, gg, kk, ll) = mo(pp, gg, kk, ll) + worb(o) * (-aimag(sum(wc * (zc - enu(ll, gg, pp))**kk * g1)) / pi)
            end do
          end do
        end do
      end do
    end subroutine moments_ord

    !> species/spin-resolved occupied moments  mo(species, spin, k, l)  up to E_F
    subroutine moments(efv, mo)
      real(dp), intent(in) :: efv
      real(dp), intent(out) :: mo(2, 2, 0:2, 0:2)
      real(dp) :: ebot, ainf, binf
      complex(dp) :: gavg, fa, gtail
      complex(dp), allocatable :: ga(:), gb(:)
      integer :: j, o, kk, ll, gg
      allocate (ga(inp%ncont), gb(inp%ncont))
      ebot = minval(enu(:, 1:nsp, :)) - 1.2_dp
      call contour_nodes(ebot, efv, inp%ncont, zc, wc)
      mo = 0.0_dp
      do gg = 1, nsp
        do o = 1, 4
          call tail_params(inp%terminator, acoef(:, o, gg), bcoef(:, o, gg), mm(o, gg), ainf, binf)
          do j = 1, inp%ncont
            gtail = sqrt_tail(zc(j), ainf, binf)
            gavg = cnrm(o, gg) * cf_green(acoef(:, o, gg), bcoef(:, o, gg), mm(o, gg), zc(j), gtail)
            fa = sqrt(cnrm(o, gg)) * transition_green(acoef(:, o, gg), bcoef(:, o, gg), mm(o, gg), ccoef(:, o, gg), zc(j), gtail)
            ga(j) = fa / xs(1)
            gb(j) = (gavg - fa) / xs(2)
          end do
          ll = lorb(o)
          do kk = 0, 2
            mo(1, gg, kk, ll) = mo(1, gg, kk, ll) + worb(o) * (-aimag(sum(wc * (zc - enu(ll, gg, 1))**kk * ga)) / pi)
            mo(2, gg, kk, ll) = mo(2, gg, kk, ll) + worb(o) * (-aimag(sum(wc * (zc - enu(ll, gg, 2))**kk * gb)) / pi)
          end do
        end do
      end do
    end subroutine moments
  end subroutine alloy_scf

  !> screening constants (s, p, d) of the recursion Hamiltonian; H^(2) is representation dependent (STATUS.md 2026-10-06)
  pure function alloy_alpha(inp) result(al)
    type(alloy_input_t), intent(in) :: inp
    real(dp) :: al(0:2)
    ! H^(2) needs the standard constants (Fe moment 2.26 vs 2.34 with sigma 0.7); the energy-dependent route is formally representation
    ! independent but its SCF did not converge for CuZn with the standard alpha (P(z) pole of Zn d just below E_F), so it keeps sigma.
    if (trim(inp%screening) == 'hardcore' .or. (trim(inp%screening) == 'auto' .and. inp%energy_dependent .and. .not. inp%ed_asr)) then
      al = alpha_from_sigma(inp%sigma, 1.0_dp)
    else
      al = standard_alpha()
    end if
  end function alloy_alpha

  !> asymptotic (a_inf, b_inf) of the square-root tail: 'bp' = Beer-Pettifor (band edges = extreme eigenvalues of the
  !> terminated chain; falls back to avg2 if the fixed point fails), otherwise the average of the last two coefficients
  subroutine tail_params(method, a, b, n, ainf, binf)
    character(len=*), intent(in) :: method
    real(dp), intent(in) :: a(:), b(:)
    integer, intent(in) :: n
    real(dp), intent(out) :: ainf, binf
    logical :: ok
    if (trim(method) == 'bp' .and. n >= 2) then
      call bp_asymptote(a, b, n, ainf, binf, converged=ok)
      if (ok .and. binf > 0.0_dp) return
    end if
    call tail_params_avg2(a, b, n, ainf, binf)
  end subroutine tail_params
end module asr_alloy_scf
