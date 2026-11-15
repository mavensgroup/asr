!> probe_cuni [x_Cu=0.5] [maxdepth=3] [maxrank=2] [L_truth=8] [eta=0.02]
!> fcc Cu_x Ni_(1-x) from our own ASA SCF of the two elements (non-magnetic, common lattice constant a = 3.57 A,
!> common screening), augmented-space recursion (depth, rank cut) vs KPM on random realisations.
program probe_cuni
  use asr_kinds, only: dp, pi
  use asr_asa
  use asr_atom, only: shell_t
  use asr_lmto_cell
  use asr_lmto_ham, only: potpar_t
  use asr_augspace, only: graph_t
  use asr_augorb
  use asr_lmto_alloy
  use asr_lmto_realization
  use asr_recursion, only: lanczos
  use asr_terminator
  use asr_kpm
  use asr_kpm, only: spectral_bounds
  use asr_contour, only: contour_nodes
  use asr_greens, only: cf_green, sqrt_tail, tail_params_avg2
  implicit none
  type(asa_input_t) :: inp
  type(asa_result_t) :: res
  type(cell_t) :: cell
  type(potpar_t) :: ppa, ppb
  type(graph_t) :: g, gt
  type(orb_alloy_t) :: mdl, mdlt
  type(augorb_op) :: op
  type(lmto_real_op) :: rop
  integer, allocatable :: coord(:, :)
  real(dp), allocatable :: v0(:), a(:), b(:), mu(:), vv(:)
  real(dp) :: x, eta, scale, err, c0, s0, emin, emax
  integer :: maxdepth, maxrank, lt, depth, iorb, io, ie, im, nst, m, nreal, nvec, nmom, ir, nd
  integer, parameter :: ne = 801
  real(dp) :: egrid(ne), rho(ne), rhosum(ne), rhob(ne)
  complex(dp) :: z, gg
  type(terminator_t) :: t
  complex(dp) :: zc(80), wc(80), gtl
  real(dp) :: ef_t, mt(0:2), ma(0:2), enu_s, ainf, binf, de_t
  integer :: kk, jj
  character(len=16) :: arg
  integer :: iarg
  real(dp) :: ss_tmp
  x = 0.5_dp; maxdepth = 3; maxrank = 2; lt = 8; eta = 0.02_dp
  iarg = command_argument_count()
  if (iarg >= 1) then; call get_command_argument(1, arg); read (arg, *) x; end if
  if (iarg >= 2) then; call get_command_argument(2, arg); read (arg, *) maxdepth; end if
  if (iarg >= 3) then; call get_command_argument(3, arg); read (arg, *) maxrank; end if
  if (iarg >= 4) then; call get_command_argument(4, arg); read (arg, *) lt; end if
  if (iarg >= 5) then; call get_command_argument(5, arg); read (arg, *) eta; end if

  ! ---- species parameters from our own SCF (non-magnetic)
  allocate (inp%cls(1)); cell%nat = 1
  cell%avec = reshape([0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp], [3, 3])
  allocate (cell%tau(3, 1), cell%isp(1)); cell%tau = 0.0_dp; cell%isp = 1
  inp%alat = 3.57_dp / 0.529177_dp; inp%nk = 10; inp%maxit = 80
  call setup_cell(cell); inp%cls(1)%sws = cell%w * inp%alat
  inp%cls(1)%core = core_shells_for([4, 4, 3, 4]); inp%cls(1)%nvl = [4, 4, 3]
  inp%cls(1)%z = 29.0_dp; inp%cls(1)%nval = 11
  inp%cls(1)%atom_shells = [inp%cls(1)%core, shell_t(4, 0, 1.0_dp), shell_t(3, 2, 10.0_dp)]
  call asa_scf(inp, cell, res); call to_pp(res, ppa); print '(a,f10.5)', 'Cu E_F=', res%ef
  inp%cls(1)%z = 28.0_dp; inp%cls(1)%nval = 10
  inp%cls(1)%atom_shells = [inp%cls(1)%core, shell_t(4, 0, 2.0_dp), shell_t(3, 2, 8.0_dp)]
  call asa_scf(inp, cell, res); call to_pp(res, ppb); print '(a,f10.5)', 'Ni E_F=', res%ef

  ! ---- truth: KPM on random realisations (supercell lt^3 cubic cells)
  call build_lmto_alloy('fcc', lt, 1.2_dp, 0.7_dp, ppa, ppb, x, gt, coord, mdlt)
  iorb = 1                                           ! s orbital
  emin = -1.0_dp; emax = 1.2_dp
  do ie = 1, ne; egrid(ie) = emin + (emax - emin) * (ie - 1) / real(ne - 1, dp); end do
  nmom = 500; nvec = 40; nreal = 4
  rhosum = 0.0_dp
  allocate (mu(0:nmom - 1))
  call build_lmto_real(gt, mdlt, 299, rop)
  call spectral_bounds(rop, 80, 0.15_dp, emin, emax)
  c0 = 0.5_dp * (emin + emax); s0 = 0.5_dp * (emax - emin) * 1.02_dp
  print '(a,2f10.4)', 'spectral bounds (Ry): ', emin, emax
  emin = -1.0_dp; emax = 1.2_dp
  do ie = 1, ne; egrid(ie) = emin + (emax - emin) * (ie - 1) / real(ne - 1, dp); end do
  do ir = 1, nreal
    call build_lmto_real(gt, mdlt, 300 + ir, rop)
    call kpm_orbital(rop, iorb, c0, s0, nmom, nvec, 11 + ir, mu)
    call kpm_dos(mu, c0, s0, egrid, rho)
    rhosum = rhosum + rho
  end do
  rho = rhosum / nreal
  call lorentz_broaden(egrid, rho, eta, rhob)
  scale = maxval(rhob)
  print '(a,i0,a,f5.2,a,f6.3,a,i0)', 'truth: ', gt%n * nreal, ' sites, x_Cu=', x, '  eta=', eta, '   s-orbital LDOS peak = ', nint(scale)

  ! ---- truth moments of the s-orbital LDOS up to E_F (about E_nu of the species-average s level)
  ef_t = -0.10_dp; enu_s = -0.5_dp
  de_t = egrid(2) - egrid(1)
  do kk = 0, 2
    mt(kk) = sum(rho * (egrid - enu_s)**kk, mask=egrid <= ef_t) * de_t
  end do
  print '(a,f7.3,a,3f10.5)', 'truth s-orbital moments up to E_F=', ef_t, ' (M0,M1,M2) =', mt
  ! ---- ASR
  call build_lmto_alloy('fcc', 2 * (4 * maxdepth) / 2 + 3, 1.2_dp, 0.7_dp, ppa, ppb, x, g, coord, mdl)
  print '(a)', 'depth   nodes(orb 9)  | max|dLDOS|/peak:  sqrt-last / sqrt-avg2 / BP / BP+LN'
  do depth = 2, maxdepth
    call build_augorb_operator(g, mdl, 1, depth, op, maxrank)
    nd = op%dim()
    allocate (a(depth), b(depth), v0(nd)); v0 = 0.0_dp; v0(iorb) = 1.0_dp
    call lanczos(op, v0, depth, a, b, m, reorth=.true.)
    write (*, '(i4,i12,3x)', advance='no') depth, op%nn
    do im = 1, 4
      select case (im)
      case (1); t = make_terminator('sqrt-last', 'butt', a, b, depth, 4)
      case (2); t = make_terminator('sqrt-avg2', 'butt', a, b, depth, 4)
      case (3); t = make_terminator('beer-pettifor', 'butt', a, b, depth, 4)
      case (4); t = make_terminator('beer-pettifor', 'luchini-nex', a, b, depth, 4)
      end select
      err = 0.0_dp
      do ie = 1, ne
        z = cmplx(egrid(ie), eta, dp); gg = terminated_green(t, z)
        err = max(err, abs(-aimag(gg) / pi - rhob(ie)))
      end do
      write (*, '(es11.3)', advance='no') err / scale
    end do
    write (*, *)
    ! integrated moments (what the SCF uses) from the contour
    call contour_nodes(-1.1_dp, ef_t, 80, zc, wc)
    call tail_params_avg2(a, b, depth, ainf, binf)
    do kk = 0, 2
      ma(kk) = 0.0_dp
      do jj = 1, 80
        gtl = sqrt_tail(zc(jj), ainf, binf)
        ma(kk) = ma(kk) - aimag(wc(jj) * (zc(jj) - enu_s)**kk * cf_green(a, b, depth, zc(jj), gtl)) / pi
      end do
    end do
    print '(a,3f10.5,a,3f8.2,a)', '        ASR s moments (M0,M1,M2) =', ma, '   error vs truth (%):', 100.0_dp * (ma - mt) / abs(mt), '  (sqrt-avg2 tail)'
    deallocate (a, b, v0)
  end do
contains
  subroutine to_pp(r, pp)
    type(asa_result_t), intent(in) :: r
    type(potpar_t), intent(out) :: pp
    integer :: l
    do l = 0, 2
      pp%enu(l) = r%cls(1)%pp(l, 1)%enu; pp%c(l) = r%cls(1)%pp(l, 1)%c; pp%sqrdel(l) = r%cls(1)%pp(l, 1)%sqrdel
      pp%p(l) = r%cls(1)%pp(l, 1)%p; pp%gamma(l) = r%cls(1)%pp(l, 1)%gamma; pp%alpha_rep(l) = r%cls(1)%pp(l, 1)%gamma
    end do
  end subroutine to_pp

  !> KPM moments with random vectors supported on one orbital of every site (trace over that orbital).
  subroutine kpm_orbital(opr, orb, cc, ss, nm, nv, seed, mu_)
    use asr_kpm, only: kpm_moments
    type(lmto_real_op), intent(in) :: opr
    integer, intent(in) :: orb, nm, nv, seed
    real(dp), intent(in) :: cc, ss
    real(dp), intent(out) :: mu_(0:nm - 1)
    real(dp), allocatable :: r(:), t0(:), t1(:), t2(:), w(:)
    integer, allocatable :: sd(:)
    integer :: n, iv, k, dim, i
    real(dp) :: ss_tmp
    dim = opr%dim()
    allocate (r(dim), t0(dim), t1(dim), t2(dim), w(dim))
    call random_seed(size=k); allocate (sd(k)); sd = seed + 101 * [(iv, iv=1, k)]; call random_seed(put=sd)
    mu_ = 0.0_dp
    do iv = 1, nv
      r = 0.0_dp
      do i = 1, opr%n
        call random_number(ss_tmp)
        r(opr%nl * (i - 1) + orb) = merge(1.0_dp, -1.0_dp, ss_tmp > 0.5_dp)
      end do
      t0 = r
      call opr%apply(r, w); t1 = (w - cc * r) / ss
      mu_(0) = mu_(0) + dot_product(r, t0) / opr%n
      mu_(1) = mu_(1) + dot_product(r, t1) / opr%n
      do n = 2, nm - 1
        call opr%apply(t1, w)
        t2 = 2.0_dp * (w - cc * t1) / ss - t0
        mu_(n) = mu_(n) + dot_product(r, t2) / opr%n
        t0 = t1; t1 = t2
      end do
    end do
    mu_ = mu_ / real(nv, dp)
  end subroutine kpm_orbital
end program probe_cuni
