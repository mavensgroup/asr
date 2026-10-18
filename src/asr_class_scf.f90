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
!>     V_p = sum_kappa' Phi(kappa, kappa') qbar_kappa'  -  sim_beta 2 (q_p - qbar_kappa) / R1,   qbar_kappa = sum_{p in kappa} x_p q_p:
!> Madelung of the sublattice-average charges over the decoration cell (Ewald) plus the screened-impurity term for the deviation
!> from the sublattice average (binary alloy: qbar = 0 -> the SIM of asr_alloy_scf; ordered sublattice: no SIM term)
!> -> Anderson mixing of the smooth potentials and E_nu.
module asr_class_scf
  use asr_kinds, only: dp, pi
  use asr_radial
  use asr_xc
  use asr_atom, only: cumint, poisson_hartree
  use asr_potpar, only: potpar_l_t, potpar_at_enu, clamp_enu
  use asr_lmto_ham, only: potpar_t
  use asr_asa, only: asa_class_t, asa_input_t, asa_class_result_t, atomic_start
  use asr_strconst, only: avg_ws_radius, standard_alpha
  use asr_augspace, only: graph_t, graph_distances
  use asr_augorb, only: orb_alloy_t
  use asr_augcls
  use asr_edc
  use asr_endep, only: ed_depth_tol, ed_nchain, ed_nlev, ed_tail
  use asr_symmetry, only: class_ops, orbital_rep
  use asr_lmto_alloy, only: build_lmto_alloy
  use asr_sublattice
  use asr_contour, only: contour_nodes_closed, contour_nodes_fermi, fermi_npole, fermi_nodes_count
  use asr_mixing, only: anderson_t, anderson_init, anderson_step
  implicit none
  private
  public :: class_scf_input_t, class_scf_result_t, class_scf

  type :: class_scf_input_t
    type(asa_class_t), allocatable :: spec(:)          !! chemical species
    type(decoration_t) :: dec                          !! sublattices on the graph lattice
    integer, allocatable :: nsc(:), spc(:, :)           !! per sublattice: number of species, chemical species spc(k, sublattice)
    real(dp), allocatable :: conc(:, :)                 !! conc(k, sublattice)
    real(dp) :: alat = 0.0_dp                          !! graph lattice constant (Bohr)
    character(len=8) :: lattice = 'bcc', xc = 'vbh'
    real(dp) :: rcl = 1.2_dp, rinv = 0.0_dp            !! kept S^alpha range, inversion cluster radius (alat)
    integer :: lsuper = 8, nlevels = 4, nrec = 0, maxrank = 1, flip_radius = 99, ncont = 32
    real(dp) :: beta = 0.3_dp, tol = 1.0e-5_dp, sim_beta = 1.0_dp, zeta0 = 0.1_dp
    integer :: maxit = 60, mix_history = 6, nspin = 1
    logical :: relativistic = .true.
    real(dp) :: ed_eta = 0.01_dp, fermi_kt = 0.001_dp, depth_tol = 1.0e-5_dp
    integer :: fermi_poles = 0
    character(len=8) :: ed_tail = 'bp'
    integer :: ndos = 0
    real(dp) :: dos_emin = -1.0_dp, dos_emax = 1.0_dp, dos_eta = 0.02_dp
    character(len=256) :: pot_in = '', pot_out = ''
  end type class_scf_input_t

  type :: class_scf_result_t
    integer :: ncomp = 0, nspin = 1, iter = 0
    integer, allocatable :: comp_cls(:), comp_spec(:), nodes(:)
    real(dp), allocatable :: comp_w(:), comp_x(:)       !! weight per site f_kappa x; concentration on the sublattice x
    type(asa_class_result_t), allocatable :: comp(:)
    real(dp) :: ef = 0.0_dp, dv = 0.0_dp, moment_site = 0.0_dp
    logical :: converged = .false.
    real(dp), allocatable :: hist_dv(:), hist_ef(:), hist_moment(:, :), hist_qnet(:, :)
    real(dp), allocatable :: edos_e(:), edos(:, :, :, :)   !! DOS(e, l, component, spin) per atom
  end type class_scf_result_t

  integer, parameter :: iorb(4) = [1, 3, 5, 7], lorb(4) = [0, 1, 2, 2], nvert = 8
  real(dp), parameter :: worb(4) = [1.0_dp, 3.0_dp, 3.0_dp, 2.0_dp]

  type :: class_op_t
    type(augcls_op) :: op
    real(dp), allocatable :: d(:, :, :)
    type(edn_work_t) :: wk
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
    real(dp), allocatable :: xin(:), xout(:), xnext(:), wmix(:), nval(:, :), ncore(:, :), ntot(:), vh(:), vxc(:, :), rho(:, :), exc(:)
    real(dp), allocatable :: tmp(:), cum(:), yc(:, :, :), vtmp(:, :), etmp(:, :)
    real(dp) :: alpha(0:2), c, avw, r1, fs, nelec, ef, nn, dvc, dv, e, elo, ehi, flo, fhi, vm, qcell
    integer :: nsp, ncomp, ncls, p, kap, k, n, nmax, it, i, sg, l, nx, jf, nfe, side, ngen, nlev, is, nod, mg
    logical :: vb, cubic, hlo, hhi, core_ok

    vb = .false.; if (present(verbose)) vb = verbose
    nsp = inp%nspin; fs = 2.0_dp / real(nsp, dp)
    c = c_nonrel; if (inp%relativistic) c = c_ry
    ncls = inp%dec%ncls
    avw = avg_ws_radius(trim(inp%lattice)) * inp%alat
    alpha = standard_alpha()
    r1 = nn_distance(trim(inp%lattice)) * inp%alat
    ngen = inp%nlevels
    nlev = inp%nlevels; if (inp%nrec > 0) nlev = inp%nrec
    ! ---- components
    ncomp = sum(inp%nsc)
    allocate (pcl(ncomp), psp(ncomp), w(ncomp), xp(ncomp))
    p = 0
    do kap = 1, ncls
      do k = 1, inp%nsc(kap)
        p = p + 1; pcl(p) = kap; psp(p) = inp%spc(k, kap); xp(p) = inp%conc(k, kap); w(p) = inp%dec%frac(kap) * inp%conc(k, kap)
      end do
    end do
    nelec = sum([(w(p) * real(inp%spec(psp(p))%nval, dp), p=1, ncomp)])
    mdl%nl = 9; mdl%nsp = ncomp; mdl%ncls = ncls
    allocate (mdl%nsc(ncls), mdl%spc(maxval(inp%nsc), ncls), mdl%conc(maxval(inp%nsc), ncls))
    mdl%nsc = inp%nsc; mdl%spc = 0; mdl%conc = 0.0_dp
    p = 0
    do kap = 1, ncls
      do k = 1, inp%nsc(kap)
        p = p + 1; mdl%spc(k, kap) = p; mdl%conc(k, kap) = inp%conc(k, kap)
      end do
    end do
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
    allocate (cop(ncls), res%nodes(ncls))
    do kap = 1, ncls
      call decoration_classes(inp%dec, coord, mg, kap, cls)
      call class_ops(coord, mg, cls, perm, mats, cubic)
      if (.not. cubic) error stop 'class_scf: a sublattice site is not cubic (needs the p_x, p_y / t2g, e_g partners separately)'
      call orbital_rep(mats, cop(kap)%d)
      call build_augcls_operator(g, mdl, cls, perm, cop(kap)%d, 1, ngen, cop(kap)%op, inp%maxrank, flippable)
      call edn_work_init(cop(kap)%wk, cop(kap)%op%dim())
      res%nodes(kap) = cop(kap)%op%nn
      if (vb) print '(a,i0,a,i0,a,i0,a,i0)', '  sublattice ', kap, ': site group order (in D4h) ', size(perm, 2), ', components ', &
        inp%nsc(kap), ', augmented-space nodes ', cop(kap)%op%nn
    end do
    if (vb) print '(a,i0,a,i0,a,i0,a,i0,a,i0)', '  class ASR: components=', ncomp, '  levels=', ngen, '  recursion steps=', nlev, &
      '  rank<=', inp%maxrank, '  nspin=', nsp
    ! ---- mixing vector: per component and spin the smooth potential (n) and E_nu (3)
    allocate (offx(ncomp + 1)); offx(1) = 0
    do p = 1, ncomp; offx(p + 1) = offx(p) + (cr(p)%mesh%n + 3) * nsp; end do
    nx = offx(ncomp + 1)
    allocate (xin(nx), xout(nx), xnext(nx))
    call anderson_init(mixer, nx, inp%mix_history)
    allocate (mom(ncomp, 2, 0:2, 0:2), dosef(ncomp, 2, 0:2))
    ed_depth_tol = inp%depth_tol; ed_tail = inp%ed_tail
    ef = maxval(enu(2, 1:nsp, :)) + 0.35_dp
    allocate (res%hist_dv(inp%maxit), res%hist_ef(inp%maxit), res%hist_moment(inp%maxit, ncomp), res%hist_qnet(inp%maxit, ncomp))

    do it = 1, inp%maxit
      ed_nchain = 0; ed_nlev = 0
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
        call ed_moments(6.0_dp, mom, dosef)
        do p = 1, ncomp
          print '(a,i0,a,3f9.4)', '  sum rule (E_F=6 Ry), component ', p, ': M0(s,p,d)/(1,3,5) =', (mom(p, 1, 0, l) / real(2 * l + 1, dp), l=0, 2)
        end do
      end if
      ! E_F: Newton steps until N - N_elec changes sign, then Illinois false position on the bracket (as asr_alloy_scf)
      nfe = 16
      hlo = .false.; hhi = .false.; side = 0
      do jf = 1, nfe
        call ed_moments(ef, mom, dosef)
        nn = fs * sum([(w(p) * sum(mom(p, 1:nsp, 0, :)), p=1, ncomp)])
        dvc = fs * sum([(w(p) * sum(dosef(p, 1:nsp, :)), p=1, ncomp)])
        if (vb .and. it == 1) print '(a,i2,a,f10.6,a,f9.5,a,f9.4)', '    E_F search ', jf, ': E_F=', ef, '  N-nelec=', nn - nelec, '  DOS(E_F)=', dvc
        if (abs(nn - nelec) < 2.0e-4_dp) exit
        if (nn < nelec) then
          if (side == -1 .and. hhi) fhi = 0.5_dp * fhi
          elo = ef; flo = nn - nelec; hlo = .true.; side = -1
        else
          if (side == 1 .and. hlo) flo = 0.5_dp * flo
          ehi = ef; fhi = nn - nelec; hhi = .true.; side = 1
        end if
        if (hlo .and. hhi) then
          e = (elo * fhi - ehi * flo) / (fhi - flo)
          if (.not. (e > min(elo, ehi) + 1.0e-3_dp * abs(ehi - elo) .and. e < max(elo, ehi) - 1.0e-3_dp * abs(ehi - elo))) e = 0.5_dp * (elo + ehi)
          e = e - ef
        else
          e = max(-0.15_dp, min(0.15_dp, (nelec - nn) / max(dvc, 1.0_dp)))
        end if
        if (jf == nfe) then
          e = max(-0.15_dp, min(0.15_dp, (nelec - nn) / max(dvc, 0.1_dp)))
          if (vb) print '(a,i0,a,es10.3,a,f9.5)', '    warning: E_F search not converged after ', nfe, ' passes: N - N_elec =', &
            nn - nelec, ';  first-order E_F shift', e
          do p = 1, ncomp; do sg = 1, nsp; do l = 0, 2; do k = 0, 2
            mom(p, sg, k, l) = mom(p, sg, k, l) + e * (ef + 0.5_dp * e - enu(l, sg, p))**k * dosef(p, sg, l)
          end do; end do; end do; end do
        end if
        ef = ef + e
      end do
      ! densities, core, charges
      dv = 0.0_dp
      do p = 1, ncomp
        n = cr(p)%mesh%n
        allocate (nval(n, nsp), ncore(n, nsp), ntot(n), tmp(n), cum(n), yc(2, 0:0, n))
        nval = 0.0_dp; ncore = 0.0_dp
        do sg = 1, nsp
          do l = 0, 2
            nval(:, sg) = nval(:, sg) + fs * ( &
                mom(p, sg, 0, l) * (y(1, 0, 1:n, l, p, sg)**2 + y(2, 0, 1:n, l, p, sg)**2 / c**2) &
              + 2.0_dp * mom(p, sg, 1, l) * (y(1, 0, 1:n, l, p, sg) * y(1, 1, 1:n, l, p, sg) + y(2, 0, 1:n, l, p, sg) * y(2, 1, 1:n, l, p, sg) / c**2) &
              + mom(p, sg, 2, l) * (y(1, 1, 1:n, l, p, sg)**2 + y(2, 1, 1:n, l, p, sg)**2 / c**2 &
                                    + y(1, 0, 1:n, l, p, sg) * y(1, 2, 1:n, l, p, sg) + y(2, 0, 1:n, l, p, sg) * y(2, 2, 1:n, l, p, sg) / c**2))
          end do
          do is = 1, size(inp%spec(psp(p))%core)
            nod = inp%spec(psp(p))%core(is)%n - inp%spec(psp(p))%core(is)%l - 1
            call bound_state(cr(p)%mesh, cr(p)%v(:, sg), inp%spec(psp(p))%core(is)%l, nod, inp%spec(psp(p))%z, c, n, &
                             -1.6_dp * inp%spec(psp(p))%z**2 - 20.0_dp, cr(p)%v(n, sg) - 1.0e-6_dp, e, yc, core_ok)
            ncore(:, sg) = ncore(:, sg) + inp%spec(psp(p))%core(is)%occ / real(nsp, dp) * (yc(1, 0, :)**2 + yc(2, 0, :)**2 / c**2)
          end do
        end do
        cr(p)%nrad = nval + ncore
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
        vm = dot_product(phimad(pcl(p), :), qbar - qcell) - inp%sim_beta * 2.0_dp * (cr(p)%qnet - qbar(pcl(p))) / r1
        do sg = 1, nsp
          i = offx(p) + (sg - 1) * (n + 3)
          xout(i + 2:i + n) = vh(2:n) + vxc(2:n, sg) + vm
          xout(i + 1) = xout(i + 2)
          dvc = maxval(abs(xout(i + 2:i + n) - xin(i + 2:i + n)) * cr(p)%mesh%r(2:n) / (1.0_dp + cr(p)%mesh%r(2:n)))
          dv = max(dv, dvc)
          do l = 0, 2
            xout(i + n + 1 + l) = enu(l, sg, p) + 0.7_dp * mom(p, sg, 1, l) / max(mom(p, sg, 0, l), 1.0e-6_dp)
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
      call anderson_step(mixer, xin, xout, inp%beta, xnext, wmix)
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
      res%iter = it; res%dv = dv; res%ef = ef
      res%hist_dv(it) = dv; res%hist_ef(it) = ef
      res%hist_moment(it, :) = cr(:)%moment; res%hist_qnet(it, :) = cr(:)%qnet
      if (vb) then
        print '(a,i3,a,es10.3,a,f10.6,a,*(f8.4))', 'it', it, '  dV=', dv, '  EF=', ef, '  qnet=', cr(:)%qnet
        if (nsp == 2) print '(a,*(f8.4))', '                                      m   =', cr(:)%moment
        if (ed_nchain > 0) print '(a,f6.1,a,i0,a)', '      mean recursion depth used ', real(ed_nlev, dp) / real(ed_nchain, dp), ' of ', nlev, ' levels'
      end if
      if (dv < inp%tol .and. it > 3) then; res%converged = .true.; exit; end if
    end do
    call xc_finalize(xf)
    if (len_trim(inp%pot_out) > 0) call write_potentials()
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

    !> occupied moments mo(component, spin, k, l) on the contour to efv and the DOS at E_F (node nq + 1)
    subroutine ed_moments(efv, mo, dosv)
      real(dp), intent(in) :: efv
      real(dp), intent(out) :: mo(:, :, 0:, 0:), dosv(:, :, 0:)
      complex(dp), allocatable :: zq(:), wq(:)
      complex(dp) :: gz
      integer :: j, o, gg, kk, ll, nq, pp_, kp
      call contour(efv, zq, wq, nq)
      mo = 0.0_dp; dosv = 0.0_dp
      do gg = 1, nsp
        do pp_ = 1, ncomp
          kp = pcl(pp_)
          call edc_setup(mdl, pps(:, gg), alpha, pp_, ed)
          call augcls_set_params(cop(kp)%op, mdl)
          call augcls_freeze_site(cop(kp)%op, 1, pp_)
          do o = 1, 4
            ll = lorb(o)
            call augcls_set_start(cop(kp)%op, iorb(o), cop(kp)%d)
            do j = 1, nq + 1
              call edc_green(cop(kp)%op, cop(kp)%wk, ed, zq(j), iorb(o), nlev, gz)
              if (j == nq + 1) then
                dosv(pp_, gg, ll) = dosv(pp_, gg, ll) - worb(o) * aimag(gz) / pi
                cycle
              end if
              do kk = 0, 2
                mo(pp_, gg, kk, ll) = mo(pp_, gg, kk, ll) - worb(o) * aimag(wq(j) * (zq(j) - enu(ll, gg, pp_))**kk * gz) / pi
              end do
            end do
          end do
        end do
      end do
    end subroutine ed_moments

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
      complex(dp) :: gz
      integer :: ie, o, gg, ll, pp_, kp
      allocate (res%edos_e(inp%ndos), res%edos(inp%ndos, 0:2, ncomp, nsp)); res%edos = 0.0_dp
      do ie = 1, inp%ndos
        res%edos_e(ie) = inp%dos_emin + (inp%dos_emax - inp%dos_emin) * (ie - 1) / real(max(inp%ndos - 1, 1), dp)
      end do
      call all_params()
      do gg = 1, nsp
        do pp_ = 1, ncomp
          kp = pcl(pp_)
          call edc_setup(mdl, pps(:, gg), alpha, pp_, ed)
          call augcls_set_params(cop(kp)%op, mdl)
          call augcls_freeze_site(cop(kp)%op, 1, pp_)
          do o = 1, 4
            ll = lorb(o)
            call augcls_set_start(cop(kp)%op, iorb(o), cop(kp)%d)
            do ie = 1, inp%ndos
              call edc_green(cop(kp)%op, cop(kp)%wk, ed, cmplx(res%edos_e(ie), inp%dos_eta, dp), iorb(o), nlev, gz)
              res%edos(ie, ll, pp_, gg) = res%edos(ie, ll, pp_, gg) - worb(o) * aimag(gz) / pi
            end do
          end do
        end do
      end do
    end subroutine dos_final

    !> component potentials (with the nucleus) and E_nu: nspin, ncomp; per component z n r(n), E_nu(0:2) per spin, v(r_k, spin)
    subroutine write_potentials()
      integer :: u, pp_, gg, kk
      open (newunit=u, file=trim(inp%pot_out), status='replace', action='write')
      write (u, '(a)') '# asr sublattice component potentials (Ry): nspin ncomp; per component: z n r(n); E_nu(0:2) per spin; v(r_k, spin)'
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
      open (newunit=u, file=trim(inp%pot_in), status='old', action='read', iostat=ios)
      if (ios /= 0) error stop 'class_scf: cannot open start_potential'
      read (u, *); read (u, *) nspf, ncf
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
