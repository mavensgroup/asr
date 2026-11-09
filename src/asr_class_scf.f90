!> Self-consistent alloy with sublattices (site classes) from the energy-dependent augmented-space recursion (asr_ed form,
!> asr_augcls + asr_edc): ordered compounds, partial (sublattice) disorder and multicomponent alloys.  Collinear spin.
!>
!> A component p = (sublattice kappa, chemical species s) has its own potential, E_nu and potential parameters; its weight per
!> site is w_p = f_kappa x_{s,kappa} (f_kappa: fraction of the sites on sublattice kappa).  Equal ASA spheres (one common S).
!> Per iteration:  component potentials -> potential parameters -> for every sublattice one class operator with the origin on
!> that sublattice; for every component p of it the start site is frozen to p and the operator re-parametrised with p's Delta
!> (asr_edc_setup) -> G_p(z) on the Fermi contour, four start orbitals (s, p_z, t2g, e_g; cubic site symmetry is checked)
!> -> moments M_k (common E_F from sum_p w_p N_p = sum_p w_p Z_val,p) -> density M0 phi^2 + 2 M1 phi phidot + M2 (phidot^2 +
!> phi phiddot), relaxed core, Hartree, xc ->
!>     V_p = sum_kappa' Phi(kappa, kappa') qbar_kappa'  +  sum_p' F(p, p') (q_p' - qbar_kappa'),   qbar_kappa = sum_{p in kappa} x_p q_p:
!> Madelung of the sublattice-average charges over the decoration cell (Ewald) plus the screened potential of the deviations
!> from the sublattice averages, F from asr_screening (charge-excess model of the graph lattice, hardness solved from sim_beta).
!> Without SRO F = -2 sim_beta / R1 (the screened-impurity model; ordered sublattice: no deviation).  With first-shell SRO,
!> F = (G00 + G1 Ahat)^-1 (W00 + W1 Ahat), Ahat(p, p') = z(kappa, kappa') (X_s'^(p) - x_s'): the conditional shell charge and the
!> part of the neighbours' screening clouds on the centre (gate B3, FINDINGS_LOG "B3 re-derived").
!> -> Anderson mixing of the smooth potentials and E_nu.
!> Total energy per site (Ry, Kohn-Sham form: output density, input potential), E = sum_p w_p E_p with
!>     E_p = E_band + sum eps_core - int n v~_in + 1/2 int n V_H[n] + int n eps_xc[n] + 1/2 q_p vm_p,
!> E_band = sum_l (M1_l + E_nu,l M0_l), v~ = v + 2Z/r (the nuclear terms -2Z int n/r of T and of the sphere electrostatics cancel),
!> vm_p = Madelung + screened fluctuation shift from the output charges (1/2 q vm = Madelung energy of the sublattice charges +
!> 1/2 sum w dq F dq; diag(w) F is symmetric, so this is variational).  The Fermi-smearing entropy (k T ~ 1 mRy) is not included.
!> The chains (spin x component x start orbital x energy) are independent (H depends on E) and run in parallel (OpenMP), one
!> per task; every thread keeps its own copy of the class operator.
!>
!> route = 'ed3' (docs/DERIVATION_ed3.md, asr_ed3): the full generalised ASA problem instead of H^(2).  The chains give the central
!> blocks xm0 = G_11, xm1 = (hG)_11, xm2 = (hGh)_11 and the moments are those of k-space (bands_moments): q0 = xm0 + 2 o xm1 + o^2 xm2,
!> q1 = xm1 + o xm2, q2 = xm2, density q0 phi^2 + 2 q1 phi phidot + q2 phidot^2, channel charge q0 + p^gamma q2 = xm0 + 2 o xm1 + p xm2.
!> E_band = sum_l (q1 + E_nu,l q0) as in asa_scf.  E_F comes from the channel charges.  gz(...) then holds G + 2 o hG + p hGh, whose
!> -Im / pi is the DOS of the sphere charge, so dos_final and ef_scan need no other change.
module asr_class_scf
  use asr_kinds, only: dp, pi
  use asr_radial
  use asr_xc
  use asr_atom, only: cumint, poisson_hartree
  use asr_potpar, only: potpar_l_t, potpar_at_enu, clamp_enu
  use asr_lmto_ham, only: potpar_t, to_alpha, o_alpha
  use asr_asa, only: asa_class_t, asa_input_t, asa_class_result_t, atomic_start
  use asr_strconst, only: avg_ws_radius, standard_alpha
  use asr_augspace, only: graph_t, graph_distances
  use asr_augorb, only: orb_alloy_t
  use asr_augcls
  use asr_edc
  use asr_ed3
  use asr_augblk, only: blk_work_t, blk_work_init, blk_nbmax
  use asr_endep, only: ed_depth_tol, ed_nchain, ed_nlev, ed_nfall, ed_tail
  use asr_symmetry, only: class_ops, orbital_rep
  use asr_lmto_alloy, only: build_lmto_alloy
  use asr_sublattice
  use asr_contour, only: contour_nodes_closed, contour_nodes_fermi, fermi_npole, fermi_nodes_count
  use asr_mixing, only: anderson_t, anderson_init, anderson_step
  use asr_screening, only: screen_t, screen_lattice, screen_matrix
  use asr_mpi, only: par_rank, par_size, par_root, par_sum_z, par_sum_i8, par_sum_i
  implicit none
  private
  public :: class_scf_input_t, class_scf_result_t, class_scf

  type :: class_scf_input_t
    type(asa_class_t), allocatable :: spec(:)          !! chemical species
    type(decoration_t) :: dec                          !! sublattices on the graph lattice
    integer, allocatable :: nsc(:), spc(:, :)           !! per sublattice: number of species, chemical species spc(k, sublattice)
    real(dp), allocatable :: conc(:, :)                 !! conc(k, sublattice)
    logical :: sro = .false.                            !! first-shell short-range order (asr_sro)
    real(dp), allocatable :: sro_cond(:, :, :, :)       !! sro_cond(k, k', kappa, kappa') = P(species k' on a first-shell site of
                                                        !! sublattice kappa' | species k at a site of kappa)
    real(dp) :: alat = 0.0_dp                          !! graph lattice constant (Bohr)
    character(len=8) :: lattice = 'bcc', xc = 'vbh'
    real(dp) :: rcl = 1.2_dp, rinv = 0.0_dp            !! kept S^alpha range, inversion cluster radius (alat)
    integer :: lsuper = 8, nlevels = 4, nrec = 0, maxrank = 1, flip_radius = 99, ncont = 32
    real(dp) :: beta = 0.3_dp, tol = 1.0e-5_dp, sim_beta = 1.0_dp, zeta0 = 0.1_dp
    real(dp) :: beta_restart = 0.0_dp                  !! linear step after an Anderson restart (0: beta)
    integer :: maxit = 60, mix_history = 6, nspin = 1
    logical :: relativistic = .true.
    real(dp) :: ed_eta = 0.01_dp, fermi_kt = 0.001_dp, depth_tol = 1.0e-5_dp
    integer :: fermi_poles = 0
    character(len=8) :: ed_tail = 'hull'                !! asr_endep%ed_tail; 'bp' is not safe on complex chains
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
    type(screen_t) :: screen                               !! screened fluctuation electrostatics (asr_screening)
    real(dp), allocatable :: fscr(:, :)                    !! F(p, p'): potential of component p per unit deviation of p' (Ry/e)
    real(dp), allocatable :: qch(:, :, :)                  !! qch(component, spin, l): sphere charge per channel, one spin channel
  end type class_scf_result_t

  integer, parameter :: iorb(4) = [1, 3, 5, 7], lorb(4) = [0, 1, 2, 2], nvert = 8
  real(dp), parameter :: worb(4) = [1.0_dp, 3.0_dp, 3.0_dp, 2.0_dp]

  type :: class_op_t
    type(augcls_op) :: op
    real(dp), allocatable :: d(:, :, :)
  end type class_op_t

contains

  subroutine class_scf(inp, res, verbose)
    type(class_scf_input_t), intent(in) :: inp
    type(class_scf_result_t), intent(out) :: res
    logical, intent(in), optional :: verbose
    type(xc_lda_t) :: xf
    type(asa_input_t) :: ai
    type(potpar_t), allocatable :: pps(:, :)                ! (component, spin)
    type(potpar_l_t) :: ppl
    type(graph_t) :: g
    type(orb_alloy_t) :: m0
    type(class_alloy_t) :: mdl
    type(class_op_t), allocatable :: cop(:)
    type(edn_t) :: ed
    type(anderson_t) :: mixer
    type(asa_class_result_t), allocatable :: cr(:)
    integer, allocatable :: coord(:, :), dist(:), offx(:), cls(:), perm(:, :), mats(:, :, :), pcl(:), psp(:)
    logical, allocatable :: flippable(:)
    real(dp), allocatable :: y(:, :, :, :, :, :), enu(:, :, :), mom(:, :, :, :), dosef(:, :, :), w(:), xp(:), phimad(:, :), qbar(:)
    complex(dp), allocatable :: dmp_a(:, :, :, :, :), dmp_b(:, :, :, :, :)   ! ef_scan_dump: chain coefficients (level, node, orb, comp, spin)
    integer, allocatable :: dmp_m(:, :, :, :)
    real(dp), allocatable :: xin(:), xout(:), xnext(:), wmix(:), nval(:, :), ncore(:, :), ntot(:), vh(:), vxc(:, :), rho(:, :), exc(:)
    real(dp), allocatable :: tmp(:), cum(:), yc(:, :, :), vtmp(:, :), etmp(:, :), esum(:)
    real(dp) :: alpha(0:2), c, avw, r1, fs, nelec, ef, nn, dvc, dv, e, elo, ehi, flo, fhi, vm, qcell, rlo, rhi, etot, ep, br
    real(dp) :: ef_in, r2n, denu, rdv, efshift
    integer :: npass, efmode, pdv, sgdv, kdv
    real(dp), allocatable :: molo(:, :, :, :), mohi(:, :, :, :), dslo(:, :, :), dshi(:, :, :)
    integer, allocatable :: pk(:), zsh(:, :)
    real(dp), allocatable :: ahat(:, :)
    character(len=:), allocatable :: smsg
    integer :: p2
    integer :: nsp, ncomp, ncls, p, kap, k, n, nmax, it, i, sg, l, nx, jf, nfe, side, ngen, nlev, is, nod, mg, nmc, r2min
    logical :: vb, cubic, hlo, hhi, core_ok, ed3
    real(dp), allocatable :: qch(:, :, :), dsm(:, :, :, :), qclo(:, :, :), qchi(:, :, :), dmlo(:, :, :, :), dmhi(:, :, :, :)
    real(dp) :: phi2, asym_max, psd_max
    integer :: npsd, nnode

    ed3 = trim(inp%route) == 'ed3'
    if (ed3 .and. len_trim(inp%ef_scan_dump) > 0) error stop 'class_scf: ef_scan_dump is not available with route = ed3'
    phi2 = merge(0.0_dp, 1.0_dp, ed3)                  ! weight of phi phiddot in the density
    vb = .false.; if (present(verbose)) vb = verbose
    nsp = inp%nspin; fs = 2.0_dp / real(nsp, dp)
    c = c_nonrel; if (inp%relativistic) c = c_ry
    ncls = inp%dec%ncls
    avw = avg_ws_radius(trim(inp%lattice)) * inp%alat
    alpha = standard_alpha() * inp%alpha_scale
    r1 = nn_distance(trim(inp%lattice)) * inp%alat
    call screen_lattice(trim(inp%lattice), inp%alat, inp%sim_beta, res%screen, smsg)
    if (len(smsg) > 0) error stop smsg
    ngen = inp%nlevels
    nlev = inp%nlevels; if (inp%nrec > 0) nlev = inp%nrec
    ! ---- components
    ncomp = sum(inp%nsc)
    allocate (pcl(ncomp), psp(ncomp), pk(ncomp), w(ncomp), xp(ncomp))
    p = 0
    do kap = 1, ncls
      do k = 1, inp%nsc(kap)
        p = p + 1; pcl(p) = kap; pk(p) = k; psp(p) = inp%spc(k, kap); xp(p) = inp%conc(k, kap)
        w(p) = inp%dec%frac(kap) * inp%conc(k, kap)
      end do
    end do
    nelec = sum([(w(p) * real(inp%spec(psp(p))%nval, dp), p=1, ncomp)])
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
    call decoration_madelung(inp%dec, inp%alat, phimad)
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
      call atomic_start(ai, 1, cr(p)%mesh, c, xf, nsp, vtmp, etmp)
      cr(p)%v = vtmp
      do sg = 1, nsp; enu(:, sg, p) = etmp(:, sg); end do
      deallocate (vtmp, etmp)
    end do
    ef_in = huge(1.0_dp)
    if (len_trim(inp%pot_in) > 0) call read_potentials()
    allocate (y(2, 0:2, nmax, 0:2, ncomp, 2)); y = 0.0_dp
    call all_params()
    ! ---- graph, S blocks, one class operator per sublattice (origin on it)
    call build_lmto_alloy(trim(inp%lattice), inp%lsuper, inp%rcl, 0.7_dp, pps(1, 1), pps(1, 1), 0.5_dp, g, coord, m0, &
                          max(inp%rcl, inp%rinv), alpha)
    mdl%block = m0%block; mdl%block0 = m0%block0
    allocate (mdl%c(9, ncomp), mdl%d(9, ncomp), mdl%o(9, ncomp), mdl%e(9, ncomp))
    call edc_setup(mdl, pps(:, 1), alpha, 1, ed)
    mg = 2 * inp%lsuper
    call graph_distances(g, 1, dist)
    allocate (flippable(g%n)); flippable = dist <= inp%flip_radius; flippable(1) = .false.
    allocate (cop(ncls), res%nodes(ncls), zsh(ncls, ncls)); zsh = 0
    do kap = 1, ncls
      call decoration_classes(inp%dec, coord, mg, kap, cls)
      if (inp%sro) then                                ! first shell of the origin (site 1 at coord 0; minimum image)
        r2min = huge(1)
        do i = 2, g%n; r2min = min(r2min, sum((modulo(coord(:, i) + mg / 2, mg) - mg / 2)**2)); end do
        do i = 2, g%n
          if (sum((modulo(coord(:, i) + mg / 2, mg) - mg / 2)**2) == r2min) cls(i) = ncls + cls(i)
        end do
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
      call class_ops(coord, mg, cls, perm, mats, cubic)
      if (.not. cubic) error stop 'class_scf: a sublattice site is not cubic (needs the p_x, p_y / t2g, e_g partners separately)'
      call orbital_rep(mats, cop(kap)%d)
      call build_augcls_operator(g, mdl, cls, perm, cop(kap)%d, 1, ngen, cop(kap)%op, inp%maxrank, flippable)
      res%nodes(kap) = cop(kap)%op%nn
      if (vb) print '(a,i0,a,i0,a,i0,a,i0)', '  sublattice ', kap, ': site group order (in D4h) ', size(perm, 2), ', components ', &
        inp%nsc(kap), ', augmented-space nodes ', cop(kap)%op%nn
    end do
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
    call screen_matrix(res%screen, ahat, res%fscr)
    ! ---- mixing vector: per component and spin the smooth potential (n) and E_nu (3)
    allocate (offx(ncomp + 1)); offx(1) = 0
    do p = 1, ncomp; offx(p + 1) = offx(p) + (cr(p)%mesh%n + 3) * nsp; end do
    nx = offx(ncomp + 1)
    allocate (xin(nx), xout(nx), xnext(nx))
    call anderson_init(mixer, nx, inp%mix_history)
    allocate (mom(ncomp, 2, 0:2, 0:2), dosef(ncomp, 2, 0:2), qch(ncomp, 2, 0:2), dsm(ncomp, 2, 0:2, 0:2))
    ed_depth_tol = inp%depth_tol; ed_tail = inp%ed_tail
    ef = maxval(enu(2, 1:nsp, :)) + 0.35_dp
    if (ef_in < huge(1.0_dp)) ef = ef_in                 ! E_F of the start potential's run
    allocate (res%hist_dv(inp%maxit), res%hist_ef(inp%maxit), res%hist_moment(inp%maxit, ncomp), res%hist_qnet(inp%maxit, ncomp), &
              res%hist_etot(inp%maxit), esum(ncomp))
    br = inp%beta; if (inp%beta_restart > 0.0_dp) br = inp%beta_restart

    do it = 1, inp%maxit
      ed_nchain = 0; ed_nlev = 0; ed_nfall = 0
      do p = 1, ncomp
        n = cr(p)%mesh%n
        do sg = 1, nsp
          i = offx(p) + (sg - 1) * (n + 3)
          xin(i + 2:i + n) = cr(p)%v(2:n, sg) + 2.0_dp * inp%spec(psp(p))%z / cr(p)%mesh%r(2:n)
          xin(i + 1) = xin(i + 2)
          xin(i + n + 1:i + n + 3) = enu(:, sg, p)
        end do
      end do
      call all_params()
      if (vb .and. it == 1) then
        call ed_moments(6.0_dp, mom, dosef, qch, dsm)
        do p = 1, ncomp
          print '(a,i0,a,3f9.4)', '  sum rule (E_F=6 Ry), component ', p, ': M0(s,p,d)/(1,3,5) =', (qch(p, 1, l) / real(2 * l + 1, dp), l=0, 2)
          if (ed3) print '(a,i0,a,3f10.5)', '  ed3 component ', p, ': p^gamma(s,p,d) of this potential =', pps(p, 1)%p
        end do
        if (ed3) print '(a,f9.4,a)', '  ed3 sum rule: weighted channel charge per orbital =', &
          sum([(w(p) * sum(qch(p, 1, :)) / 9.0_dp, p=1, ncomp)]), ' (1 only if all pencil eigenvalues lie below 6 Ry; single channels differ from 1)'
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
        ! electron count stays exact (below); the moments carry an error O(step^2), which vanishes at self-consistency, where E_F
        ! no longer moves between iterations.  So the fixed point is unchanged; only the path towards it differs.
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
      dv = 0.0_dp; pdv = 0; sgdv = 0; kdv = 0; denu = 0.0_dp
      do p = 1, ncomp
        n = cr(p)%mesh%n
        allocate (nval(n, nsp), ncore(n, nsp), ntot(n), tmp(n), cum(n), yc(2, 0:0, n))
        nval = 0.0_dp; ncore = 0.0_dp; esum(p) = 0.0_dp
        do sg = 1, nsp
          do l = 0, 2
            esum(p) = esum(p) + fs * (mom(p, sg, 1, l) + enu(l, sg, p) * mom(p, sg, 0, l))
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
          end do
        end do
        cr(p)%nrad = nval + ncore
        do sg = 1, nsp
          tmp(1) = 0.0_dp
          tmp(2:n) = cr(p)%nrad(2:n, sg) * (cr(p)%v(2:n, sg) + 2.0_dp * inp%spec(psp(p))%z / cr(p)%mesh%r(2:n))
          call cumint(cr(p)%mesh, tmp, cum); esum(p) = esum(p) - cum(n)
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
      qcell = sum(inp%dec%frac * qbar)            ! residual of the E_F search: the Madelung sum needs neutral charges
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
        do p2 = 1, ncomp; vm = vm + res%fscr(p, p2) * (cr(p2)%qnet - qbar(pcl(p2))); end do
        allocate (tmp(n), cum(n))
        tmp = ntot * (0.5_dp * vh + exc); tmp(1) = 0.0_dp
        call cumint(cr(p)%mesh, tmp, cum)
        ep = esum(p) + cum(n) + 0.5_dp * cr(p)%qnet * vm
        etot = etot + w(p) * ep
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
      if (.not. allocated(wmix)) then           ! Anderson metric = the dV measure (r/(1+r))^2 on the potentials, 1 on E_nu
        allocate (wmix(nx)); wmix = 1.0_dp
        do p = 1, ncomp
          n = cr(p)%mesh%n
          do sg = 1, nsp
            i = offx(p) + (sg - 1) * (n + 3)
            wmix(i + 1:i + n) = (cr(p)%mesh%r(1:n) / (1.0_dp + cr(p)%mesh%r(1:n)))**2
          end do
        end do
      end if
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
              if (vb) print '(a,i0,a,i0,a,i0,a,f9.5)', '    E_nu outside its band window: component ', p, ' spin ', sg, ' l ', l, ' -> ', enu(l, sg, p)
            end if
          end do
        end do
      end do
      res%iter = it; res%dv = dv; res%ef = ef; res%etot = etot
      res%hist_dv(it) = dv; res%hist_ef(it) = ef; res%hist_etot(it) = etot
      res%hist_moment(it, :) = cr(:)%moment; res%hist_qnet(it, :) = cr(:)%qnet
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
        end if
        print '(a,i0,a,i0)', '      Anderson history ', mixer%nh, '  restarts so far ', mixer%nrestart
        print '(a,i0,a,a,a,es9.2,a,es10.3,a,i0,a,i0,a,f8.4,a,es10.3)', '      E_F passes ', npass, ' final ', &
          trim(merge('bracket', 'dos    ', efmode == 1)), ' shift ', efshift, '  |r|_2=', r2n, '  max dV: comp ', pdv, ' spin ', sgdv, &
          ' r=', cr(max(pdv, 1))%mesh%r(max(kdv, 1)), '  max dE_nu=', denu
      end if
      if (len_trim(inp%pot_out) > 0) call write_potentials()       ! every iteration: a stopped run restarts from here
      if (dv < inp%tol .and. it > 3) then; res%converged = .true.; exit; end if
    end do
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
    res%ncomp = ncomp; res%nspin = nsp
    res%comp_cls = pcl; res%comp_spec = psp; res%comp_w = w; res%comp_x = xp
    res%comp = cr
    res%moment_site = sum(w * cr(:)%moment)
    res%qch = qch

  contains

    subroutine all_params()
      integer :: pp_, gg, ll
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
      allocate (gz(nq + 1, 4, ncomp, nsp))
      mo = 0.0_dp; dosv = 0.0_dp; dm = 0.0_dp
      if (.not. ed3) then
        call chain_greens(zq, gz)
        do gg = 1, nsp
          do pp_ = 1, ncomp
            do o = 1, 4
              ll = lorb(o)
              dosv(pp_, gg, ll) = dosv(pp_, gg, ll) - worb(o) * aimag(gz(nq + 1, o, pp_, gg)) / pi
              do j = 1, nq
                do kk = 0, 2
                  mo(pp_, gg, kk, ll) = mo(pp_, gg, kk, ll) - worb(o) * aimag(wq(j) * (zq(j) - enu(ll, gg, pp_))**kk * gz(j, o, pp_, gg)) / pi
                end do
              end do
            end do
          end do
        end do
        qc = mo(:, :, 0, :)
        return
      end if
      allocate (gb(nq + 1, 4, ncomp, nsp, 0:3))
      call chain_greens(zq, gz, gb)
      asym_max = 0.0_dp; psd_max = 0.0_dp; npsd = 0; nnode = 0
      qc = 0.0_dp
      do gg = 1, nsp
        do pp_ = 1, ncomp
          qa = to_alpha(pps(pp_, gg), alpha); oa = o_alpha(qa, alpha)
          do o = 1, 4
            ll = lorb(o)
            xm = 0.0_dp
            do kk = 0, 2
              do j = 1, nq
                xm(kk) = xm(kk) - worb(o) * aimag(wq(j) * gb(j, o, pp_, gg, kk)) / pi
              end do
              xd(kk) = -worb(o) * aimag(gb(nq + 1, o, pp_, gg, kk)) / pi
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
        allocate (gz(nq + 1, 4, ncomp, nsp))
        if (len_trim(inp%ef_scan_dump) > 0) then
          allocate (dmp_a(nlev, nq + 1, 4, ncomp, nsp), dmp_b(nlev, nq + 1, 4, ncomp, nsp), dmp_m(nq + 1, 4, ncomp, nsp))
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
            do o = 1, 4
              dsum = dsum - fs * w(pp_) * worb(o) * aimag(gz(nq + 1, o, pp_, gg)) / pi
              do j = 1, nq
                ip = count(j > bnd) + 1
                x = -fs * worb(o) * aimag(wq(j) * gz(j, o, pp_, gg)) / pi
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
        open (newunit=u, file=trim(inp%ef_scan_dump), status='replace', action='write')
        write (u, '(a,es25.16,f6.1,3i6)') 'H ', nelec, fs, ncomp, nsp, nlev
        write (u, '(a,*(es25.16))') 'W ', w
        write (u, '(a,*(f6.1))') 'O ', worb
      else
        open (newunit=u, file=trim(inp%ef_scan_dump), status='old', position='append', action='write')
      end if
      write (u, '(a,i4,es25.16,i6)') 'E ', js, ev, nq
      do j = 1, nq + 1
        write (u, '(a,i5,4es25.16)') 'Z ', j, zq(j), wq(j)
      end do
      do gg = 1, nsp; do pp_ = 1, ncomp; do o = 1, 4; do j = 1, nq + 1
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
    subroutine chain_greens(zs, gz, gb)
      complex(dp), intent(in) :: zs(:)
      complex(dp), intent(out) :: gz(:, :, :, :)
      complex(dp), intent(out), optional :: gb(:, :, :, :, 0:)
      type(class_alloy_t), allocatable :: ms(:, :)
      type(edn_t), allocatable :: es(:, :)
      type(ed3_t), allocatable :: e3(:, :)
      integer :: gg, pp_, nbe, kp
      integer(8) :: cnt0(3), cnt(3)
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
      !$omp parallel
      call chain_worker(zs, ms, es, e3, nbe, gz, gb)
      !$omp end parallel
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
    end subroutine chain_greens

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
      do t = 1, ntg * 4 * ncomp * nsp
        ! MPI: this rank's share, cyclic in t; the shift by the group index rotates the energies over the ranks, so that no
        ! rank keeps the same energies (whose cost differs) in every group.  One rank: every task.
        if (mod(t - 1 + (t - 1) / ntg, nr) /= rk) cycle
        if (nbe == 0) then
          j = mod(t - 1, nz) + 1; o = mod((t - 1) / nz, 4) + 1
          pp_ = mod((t - 1) / (4 * nz), ncomp) + 1; gg = (t - 1) / (4 * nz * ncomp) + 1
          j1 = j; j2 = j
        else
          jb = mod(t - 1, nbk) + 1; o = mod((t - 1) / nbk, 4) + 1
          pp_ = mod((t - 1) / (4 * nbk), ncomp) + 1; gg = (t - 1) / (4 * nbk * ncomp) + 1
          j1 = (jb - 1) * nbe + 1; j2 = min(nz, jb * nbe)
        end if
        kp = pcl(pp_)
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
        if (cur(3, kp) /= o) then; call augcls_set_start(lop(kp), iorb(o), cop(kp)%d); cur(3, kp) = o; newstart(kp) = .true.; end if
        iv = iorb(o)
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

    subroutine contour(efv, zq, wq, nq)
      real(dp), intent(in) :: efv
      complex(dp), allocatable, intent(out) :: zq(:), wq(:)
      integer, intent(out) :: nq
      real(dp) :: ebot
      integer :: npl
      ebot = minval(enu(:, 1:nsp, :)) - 1.2_dp
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

    subroutine dos_final()
      complex(dp), allocatable :: zs(:), gz(:, :, :, :)
      integer :: ie, o, gg, pp_
      allocate (res%edos_e(inp%ndos), res%edos(inp%ndos, 0:2, ncomp, nsp)); res%edos = 0.0_dp
      do ie = 1, inp%ndos
        res%edos_e(ie) = inp%dos_emin + (inp%dos_emax - inp%dos_emin) * (ie - 1) / real(max(inp%ndos - 1, 1), dp)
      end do
      call all_params()
      zs = cmplx(res%edos_e, inp%dos_eta, dp)
      allocate (gz(inp%ndos, 4, ncomp, nsp))
      call chain_greens(zs, gz)
      do gg = 1, nsp
        do pp_ = 1, ncomp
          do o = 1, 4
            res%edos(:, lorb(o), pp_, gg) = res%edos(:, lorb(o), pp_, gg) - worb(o) * aimag(gz(:, o, pp_, gg)) / pi
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
      close (u)
      if (vb) print '(3a)', '  start potentials from ', trim(inp%pot_in)
    end subroutine read_potentials
  end subroutine class_scf

  function nn_distance(kind) result(d)
    character(len=*), intent(in) :: kind
    real(dp) :: d
    select case (kind)
    case ('sc'); d = 1.0_dp
    case ('bcc'); d = sqrt(3.0_dp) / 2.0_dp
    case ('fcc'); d = 1.0_dp / sqrt(2.0_dp)
    case default; error stop 'nn_distance: unknown lattice'
    end select
  end function nn_distance
end module asr_class_scf
