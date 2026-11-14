!> asr_screening: lattice Coulomb kernel and the single-site screened electrostatics with first-shell SRO (gate B3).
!>  Ewald M(k): CsCl (bcc H) and rock-salt (sc R) Madelung constants.  screen_lattice: beta round trip, the unscreened limit,
!>  the reference values of the derivation (bcc a = 2.85 A: a = 1.263 Ry gives beta 1.045, lambda 0.212), the lattice maximum
!>  of beta, mesh convergence.  screen_matrix: no SRO = SIM (linear form, bit for bit), diag(w) F symmetric, sublattice sum
!>  rule, the unscreened limit = bare first-shell charge.  Closure (screen_gamma): reference values of the numpy closure that was
!>  tested against exact supercells (tools/sro_closure/sro_closure.py, same nk = 28 mesh), slope G1 at lambda = 0,
!>  Gamma = 0 at full clustering, Gamma = G(H) for perfect B2 order (exact), positive and decreasing on the realizable range,
!>  continuous at both limits of the spherical model, finite where the first-shell form has its pole (fcc beta 1.2).
!>  screen_update: the drive iteration converges to V = F dq, and to dq = 0 on a mode with Gamma = 0.
program test_screening
  use asr_kinds, only: dp, pi
  use asr_screening
  use asr_sro, only: sro_conditional
  use asr_testing
  implicit none
  type(screen_t) :: s, s2
  character(len=:), allocatable :: msg
  real(dp) :: p(3, 3), kv(3, 1), m(1), alat, afcc, al, x(3), at(3, 3), ah(3, 3), f(3, 3), w(3), f2(2, 2), ah2(2, 2), ev
  real(dp) :: g3(3, 3), w3(3, 3), gi3(3, 3), g2(2, 2), w2(2, 2), gi2(2, 2), lam2(2), gam2(2), lam3(3), gam3(3)
  real(dp) :: h, gp, gm_, gprev, gnow, lam, dmax
  real(dp), allocatable :: cond(:, :)
  integer :: i, j, k
  logical :: ok

  print '(a)', 'test_screening'
  alat = 2.85_dp / 0.529177210903_dp
  afcc = 3.59_dp / 0.529177210903_dp
  p = reshape([-1, 1, 1, 1, -1, 1, 1, 1, -1], [3, 3]) * 0.5_dp * alat
  kv(:, 1) = [2.0_dp * pi / alat, 0.0_dp, 0.0_dp]
  call coulomb_k(p, kv, m)
  call check_close(-m(1) * sqrt(3.0_dp) / 2.0_dp * alat / 2.0_dp, 1.7626748_dp, 1.0e-6_dp, 'bcc H point: CsCl Madelung constant')
  p = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3]) * alat
  kv(:, 1) = [1.0_dp, 1.0_dp, 1.0_dp] * pi / alat
  call coulomb_k(p, kv, m)
  call check_close(-m(1) * alat / 2.0_dp, 1.7475646_dp, 1.0e-6_dp, 'sc R point: rock-salt Madelung constant')

  call screen_lattice('bcc', alat, 1.0_dp, s, msg)
  call check(len(msg) == 0 .and. .not. s%unscreened, 'bcc beta = 1 accepted')
  call check_close(-s%w00 * s%r1 / (2.0_dp * s%g00), 1.0_dp, 1.0e-10_dp, 'beta round trip: -W00 R1 / (2 G00) = beta')
  call check_close(s%w00, s%a * (1.0_dp - s%g00), 1.0e-12_dp, 'W00 = a (1 - G00)')
  call check(s%z == 8 .and. abs(s%r1 - sqrt(3.0_dp) / 2.0_dp * alat) < 1.0e-10_dp, 'bcc first shell: z = 8 at sqrt(3)/2 a')
  call check_close(s%a_min * s%r1 / 2.0_dp, 1.7626748_dp, 1.0e-6_dp, 'bcc stability bound: a_min R1 / 2 = CsCl Madelung constant (H)')
  call check(s%beta_max > 1.135_dp .and. s%beta_max < 1.7626748_dp, &
             'bcc beta maximum below a_min R1/2 (G00 stays finite at the instability)')
  call screen_lattice('bcc', alat, 1.045_dp, s2, msg)
  call check_close(s2%a, 1.263_dp, 0.01_dp, 'derivation reference: beta 1.045 <-> a 1.263 Ry')
  call check_close(s2%lambda, 0.212_dp, 0.003_dp, 'derivation reference: lambda 0.212')
  call screen_lattice('bcc', alat, 1.045_dp, s, msg, nk=36)
  call check_close(s%lambda, s2%lambda, 2.0e-3_dp, 'mesh 28 vs 36: lambda converged')
  call screen_lattice('bcc', alat, 1.8_dp, s, msg)
  call check(index(msg, 'maximum') > 0, 'beta above the lattice maximum rejected')
  call screen_lattice('fcc', afcc, 1.0_dp, s, msg)
  call check(len(msg) == 0 .and. s%z == 12, 'fcc beta = 1 accepted, z = 12')
  call check(s%lambda > 0.0_dp .and. s%lambda < 1.0_dp, 'fcc: 0 < lambda < 1')
  call check_close(s%c1_ord, -1.0_dp / 3.0_dp, 1.0e-14_dp, 'fcc ordering point X: c1 = -1/3')

  ! ---- single-site matrix, one sublattice, ternary with SRO
  call screen_lattice('bcc', alat, 1.0_dp, s, msg)
  call check_close(s%c1_ord, -1.0_dp, 1.0e-14_dp, 'bcc ordering point H: c1 = -1')
  x = [0.66_dp, 0.14_dp, 0.20_dp]; w = x
  at = 0.0_dp; at(1, 2) = -0.1_dp; at(2, 1) = -0.1_dp; at(1, 3) = 0.2_dp; at(3, 1) = 0.2_dp; at(2, 3) = 0.3_dp; at(3, 2) = 0.3_dp
  call sro_conditional(x, x, at, .true., cond, msg)
  do i = 1, 3; do j = 1, 3; ah(i, j) = s%z * (cond(i, j) - x(j)); end do; end do
  call screen_matrix(s, 0.0_dp * ah, w, fm=f, gm=g3, wm=w3)
  call check_close(maxval(abs(f - (s%w00 / s%g00) * eye3())), 0.0_dp, 0.0_dp, 'no SRO: F = (W00/G00) I = -2 beta / R1 (the SIM), bit for bit')
  call check_close(maxval(abs(g3 - s%g00 * eye3())) + maxval(abs(w3 - s%w00 * eye3())), 0.0_dp, 0.0_dp, &
                   'no SRO: Gamma = G00 I, W = W00 I exactly')
  call check_close(s%w00 / s%g00, -2.0_dp * s%beta / s%r1, 1.0e-10_dp, 'no SRO: W00/G00 = -2 beta / R1')
  call screen_matrix(s, ah, w, fm=f, gm=g3, wm=w3, lam=lam3, gam=gam3)
  call check_close(maxval(abs(spread(w, 2, 3) * f - transpose(spread(w, 2, 3) * f))), 0.0_dp, 1.0e-13_dp, &
                   'ternary SRO: diag(w) F symmetric (variational energy)')
  call check_close(maxval(abs(matmul(x, f) - (s%w00 / s%g00) * x)), 0.0_dp, 1.0e-13_dp, &
                   'ternary SRO: x^T F = (W00/G00) x^T (no potential from the sublattice mean)')
  call check_close(maxval(abs(matmul(g3, [1.0_dp, 1.0_dp, 1.0_dp]) - s%g00)), 0.0_dp, 1.0e-13_dp, &
                   'ternary SRO: Gamma 1 = G00 1 (the sublattice mean is not correlated)')
  call check_close(maxval(abs(matmul(g3, ah) - matmul(ah, g3))), 0.0_dp, 1.0e-12_dp, 'ternary SRO: Gamma is a function of Ahat (commutes)')
  call check_close(maxval(abs(w3 - s%a * (eye3() - g3))), 0.0_dp, 1.0e-12_dp, 'ternary SRO: W = a (I - Gamma)')
  call check_close(maxval(abs(f - matmul(inv3(g3), w3))), 0.0_dp, 1.0e-12_dp, 'ternary SRO: F = Gamma^-1 W')
  dmax = 0.0_dp
  do i = 1, 3; dmax = max(dmax, abs(gam3(i) - screen_gamma(s, lam3(i)))); end do
  call check_close(dmax, 0.0_dp, 1.0e-14_dp, 'ternary SRO: mode Gammas are the scalar closure at the eigenvalues of Ahat')
  ! binary 50:50: eigenvalue on (1, -1) is a [1/Gamma(z alpha) - 1]
  al = -0.15_dp
  ah2 = reshape([0.5_dp * al, -0.5_dp * al, -0.5_dp * al, 0.5_dp * al], [2, 2]) * s%z
  call screen_matrix(s, ah2, [0.5_dp, 0.5_dp], fm=f2)
  ev = f2(1, 1) - f2(1, 2)
  call check_close(ev, s%a * (1.0_dp / screen_gamma(s, s%z * al) - 1.0_dp), 1.0e-12_dp, &
                   'binary 50:50: F(1,-1) = a [1/Gamma(z alpha) - 1]')

  ! ---- closure: reference values (numpy closure, nk = 28; tested against exact supercells)
  call screen_lattice('fcc', afcc, 1.0_dp, s, msg)
  call check_close(screen_gamma(s, 7.92_dp), 0.461601648026_dp, 2.0e-8_dp, 'closure fcc beta 1.0, lambda 7.92 (condensate) = numpy')
  call check_close(screen_gamma(s, 2.0_dp), 1.131680702912_dp, 2.0e-8_dp, 'closure fcc beta 1.0, lambda 2 = numpy')
  call check_close(screen_gamma(s, -1.5_dp), 1.514195756237_dp, 2.0e-8_dp, 'closure fcc beta 1.0, lambda -1.5 = numpy')
  call check_close(s%g00 + s%g1 * 7.92_dp, 0.486728364529_dp, 1.0e-9_dp, 'same lattice functions as numpy (G00 + G1 7.92)')
  ! slope G1 at lambda = 0 (the linear form is the first order of the closure)
  h = 0.01_dp
  gp = screen_gamma(s, h); gm_ = screen_gamma(s, -h)
  call check_close((gp - gm_) / (2.0_dp * h), s%g1, 1.0e-3_dp * abs(s%g1), 'closure: dGamma/dlambda = G1 at lambda = 0')
  call check_close(screen_gamma(s, 0.0_dp), s%g00, 0.0_dp, 'closure: Gamma(0) = G00')
  ! full clustering: Gamma = 0
  call check_close(screen_gamma(s, real(s%z, dp)), 0.0_dp, 1.0e-15_dp, 'closure: Gamma = 0 at full clustering (lambda = z)')
  ! positive and decreasing on the realizable range z c1(X) .. z, continuous at both limits
  call check_positive_decreasing(s, 'fcc beta 1.0')
  call check_close(screen_gamma(s, s%z * s%alpha_cl * (1.0_dp + 1.0e-7_dp)), screen_gamma(s, s%z * s%alpha_cl * (1.0_dp - 1.0e-7_dp)), &
                   1.0e-6_dp, 'closure: continuous at the clustering limit alpha_c')
  call check_close(screen_gamma(s, s%z * s%alpha_ord * (1.0_dp + 1.0e-7_dp)), screen_gamma(s, s%z * s%alpha_ord * (1.0_dp - 1.0e-7_dp)), &
                   1.0e-6_dp, 'closure: continuous at the ordering limit')
  call check(s%alpha_cl > 0.24_dp .and. s%alpha_cl < 0.26_dp, 'fcc clustering limit alpha_c near 0.256 (mesh value)')
  ! fcc beta 1.2: the first-shell form has its pole at lambda = 7.67; the closure is finite and positive
  call screen_lattice('fcc', afcc, 1.2_dp, s, msg)
  call check(s%g00 + s%g1 * 7.92_dp < 0.0_dp, 'fcc beta 1.2: first-shell form negative at lambda 7.92 (past its pole)')
  call check_close(screen_gamma(s, 7.92_dp), 0.689508478842_dp, 2.0e-8_dp, 'closure fcc beta 1.2, lambda 7.92: positive, = numpy')
  call check_positive_decreasing(s, 'fcc beta 1.2')
  ! bcc: ordering side (CuZn), and perfect B2 order = G(H) exactly
  call screen_lattice('bcc', alat, 1.0_dp, s, msg)
  call check_close(screen_gamma(s, -1.184_dp), 1.552197257666_dp, 2.0e-8_dp, 'closure bcc beta 1.0, lambda -1.184 (CuZn alpha -0.148) = numpy')
  call check_close(screen_gamma(s, 2.0_dp), 1.079149185461_dp, 2.0e-8_dp, 'closure bcc beta 1.0, lambda 2 = numpy')
  call check_close(screen_gamma(s, -real(s%z, dp)), s%a / (s%a - 2.0_dp * 1.7626748_dp / s%r1), 1.0e-6_dp, &
                   'closure: perfect B2 order (alpha_1 = -1) gives Gamma = G(H) = a / (a + M(H)) exactly')
  call check_positive_decreasing(s, 'bcc beta 1.0')

  ! ---- full clustering in screen_matrix: binary 50:50, alpha = 1
  call screen_lattice('fcc', afcc, 1.0_dp, s, msg)
  al = 1.0_dp
  ah2 = reshape([0.5_dp * al, -0.5_dp * al, -0.5_dp * al, 0.5_dp * al], [2, 2]) * s%z
  call screen_matrix(s, ah2, [0.5_dp, 0.5_dp], gm=g2, wm=w2, ginv=gi2, lam=lam2, gam=gam2)
  call check_close(g2(1, 1) - g2(1, 2), 0.0_dp, 1.0e-12_dp, 'alpha = 1: Gamma = 0 on the mode (1, -1)')
  call check_close(w2(1, 1) - w2(1, 2), s%a, 1.0e-12_dp, 'alpha = 1: W = a on the mode (1, -1) (V = a db, finite)')
  call check_close(gi2(1, 1) - gi2(1, 2), 0.0_dp, 1.0e-12_dp, 'alpha = 1: pseudo-inverse leaves the Gamma = 0 mode out')
  call check_close(g2(1, 1) + g2(1, 2), s%g00, 1.0e-12_dp, 'alpha = 1: Gamma = G00 on the uniform mode')

  ! ---- unscreened limit: F = bare first-shell charge, unchanged
  call screen_lattice('bcc', alat, 0.0_dp, s, msg)
  call check(len(msg) == 0 .and. s%unscreened, 'beta = 0: unscreened limit')
  call screen_matrix(s, ah, w, fm=f)
  call check_close(maxval(abs(f - (2.0_dp / s%r1) * ah)), 0.0_dp, 0.0_dp, 'beta = 0: F = (2/R1) Ahat (bare shell charge), as before')

  ! ---- screen_update: drive iteration with a model response dq = P (c - chi v)
  ! (a) moderate SRO (bcc ternary above): the fixed point is v = F dq
  call screen_lattice('bcc', alat, 1.0_dp, s, msg)
  call screen_matrix(s, ah, w, fm=f, gm=g3, wm=w3, ginv=gi3)
  call drive_fixed_point(g3, w3, x, [0.10_dp, -0.30_dp, -0.06_dp], 2.0_dp, ok, dmax, f)
  call check(ok, 'drive iteration converges (moderate SRO)')
  call check_close(dmax, 0.0_dp, 1.0e-10_dp, 'drive fixed point: v = F dq (the old form) for moderate SRO')
  ! (b) fcc (1,1,0) triad: Fe has only Fe neighbours, Ahat eigenvalue 12, Gamma = 0 on that mode
  call screen_lattice('fcc', afcc, 1.0_dp, s, msg)
  at = 0.0_dp; at(1, 2) = 1.0_dp; at(2, 1) = 1.0_dp; at(1, 3) = 1.0_dp; at(3, 1) = 1.0_dp
  call sro_conditional(x, x, at, .true., cond, msg)
  call check(abs(cond(1, 1) - 1.0_dp) < 1.0e-12_dp, 'triad (1,1,0): P(Fe | Fe) = 1')
  do i = 1, 3; do j = 1, 3; ah(i, j) = s%z * (cond(i, j) - x(j)); end do; end do
  call screen_matrix(s, ah, w, gm=g3, wm=w3, ginv=gi3, lam=lam3, gam=gam3)
  call check(abs(maxval(lam3) - 12.0_dp) < 1.0e-10_dp .and. minval(gam3) < 1.0e-12_dp, 'triad (1,1,0): mode lambda = 12 has Gamma = 0')
  call drive_fixed_point(g3, w3, x, [0.10_dp, -0.30_dp, -0.06_dp], 2.0_dp, ok, dmax)
  call check(ok, 'drive iteration converges with a Gamma = 0 mode (F infinite)')
  call check_close(dmax, 0.0_dp, 1.0e-10_dp, 'drive fixed point: dq has no component on the Gamma = 0 mode')
  call report()

contains

  function eye3() result(e)
    real(dp) :: e(3, 3)
    integer :: ii
    e = 0.0_dp; do ii = 1, 3; e(ii, ii) = 1.0_dp; end do
  end function eye3

  function inv3(a) result(b)
    real(dp), intent(in) :: a(3, 3)
    real(dp) :: b(3, 3), dt
    b(1, 1) = a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2); b(1, 2) = a(1, 3) * a(3, 2) - a(1, 2) * a(3, 3)
    b(1, 3) = a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2); b(2, 1) = a(2, 3) * a(3, 1) - a(2, 1) * a(3, 3)
    b(2, 2) = a(1, 1) * a(3, 3) - a(1, 3) * a(3, 1); b(2, 3) = a(1, 3) * a(2, 1) - a(1, 1) * a(2, 3)
    b(3, 1) = a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1); b(3, 2) = a(1, 2) * a(3, 1) - a(1, 1) * a(3, 2)
    b(3, 3) = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)
    dt = a(1, 1) * b(1, 1) + a(1, 2) * b(2, 1) + a(1, 3) * b(3, 1)
    b = b / dt
  end function inv3

  !> Gamma(lambda) > 0 and non-increasing on lambda in [z c1(k_o), z) (the realizable range of one sublattice), 0 at z.
  subroutine check_positive_decreasing(sc, label)
    type(screen_t), intent(in) :: sc
    character(len=*), intent(in) :: label
    logical :: pos, dec
    pos = .true.; dec = .true.; gprev = huge(1.0_dp)
    do k = 0, 400
      lam = sc%z * (sc%c1_ord + (1.0_dp - sc%c1_ord) * k / 400.0_dp)
      gnow = screen_gamma(sc, lam)
      if (k < 400 .and. gnow <= 0.0_dp) pos = .false.
      if (gnow > gprev + 1.0e-12_dp) dec = .false.
      gprev = gnow
    end do
    call check(pos, 'closure '//label//': Gamma > 0 on the realizable range (no pole)')
    call check(dec, 'closure '//label//': Gamma decreases with lambda')
  end subroutine check_positive_decreasing

  !> Damped iteration of u with dq = P (c - chi W u) (the potential implied by u acts on a linear model response), then
  !> u_out = screen_update.  ok: converged.  dmax: |v - F dq| if fref is given, else |Gamma Gamma^+ dq - dq| (dq in the range
  !> of Gamma, i.e. no charge on a Gamma = 0 mode).
  subroutine drive_fixed_point(gm, wm, xc, c0, chi, ok, dmax, fref)
    real(dp), intent(in) :: gm(3, 3), wm(3, 3), xc(3), c0(3), chi
    logical, intent(out) :: ok
    real(dp), intent(out) :: dmax
    real(dp), intent(in), optional :: fref(3, 3)
    real(dp) :: u(3), uo(3), v(3), c(3), dq(3), gi(3, 3)
    integer :: it, grp(3)
    grp = 1; c = c0
    call screen_project(xc, grp, c)
    u = 0.0_dp; ok = .false.
    do it = 1, 2000
      dq = c - chi * matmul(wm, u)
      call screen_project(xc, grp, dq)
      call screen_update(gm, wm, xc, grp, dq, u, uo, v)
      if (maxval(abs(uo - u)) < 1.0e-14_dp) then; ok = .true.; exit; end if
      u = u + 0.3_dp * (uo - u)
    end do
    if (present(fref)) then
      dmax = maxval(abs(v - matmul(fref, dq)))
    else
      call screen_matrix(s, ah, w, ginv=gi)
      dmax = maxval(abs(matmul(gm, matmul(gi, dq)) - dq))
    end if
  end subroutine drive_fixed_point
end program test_screening
