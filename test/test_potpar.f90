!> Potential parameters: C must agree with the exact energy where D_phi(E) = -(l+1)
!> (zero of W{K,phi}); the error is second order in (E_nu - C), so it must collapse
!> when E_nu is moved toward C.  Tested for a model well with a Coulomb core, l = 0,1,2.
program test_potpar
  use asr_kinds, only: dp
  use asr_radial
  use asr_potpar
  use asr_testing
  implicit none
  type(radmesh_t) :: m
  real(dp), allocatable :: v(:), y(:, :, :)
  type(potpar_l_t) :: pp
  real(dp) :: z, s, avw, e, ec_exact, err(3), enu(3), cgam
  integer :: l, k, i
  real(dp) :: dlog, ddot, mm

  print '(a)', 'test_potpar'
  z = 29.0_dp; s = 2.7_dp; avw = 2.7_dp
  call make_radmesh(z, s, 0.012_dp, m)
  allocate (v(m%n), y(2, 0:2, m%n))
  ! model potential: screened nucleus + constant shift
  do i = 2, m%n
    v(i) = -2.0_dp * z / m%r(i) * exp(-m%r(i) * 1.4_dp) - 2.0_dp * 1.0_dp / max(m%r(i), 0.3_dp) * (1.0_dp - exp(-m%r(i))) - 0.5_dp
  end do
  v(1) = v(2)
  do l = 0, 2
    ! exact energy where D_phi(E) = -(l+1): Newton on f(E) = D + l + 1 using dD/dE from phi, phidot
    e = -0.3_dp
    do k = 1, 60
      call potpar_at_enu(m, v, l, e, z, c_ry, avw, pp)
      mm = 1.0_dp + (e - v(m%n)) / c_ry**2
      dlog = pp%dphi
      ! dD/dE = (D_phidot - D_phi) phidot/phi  (log-derivative energy dependence of the normalised family)
      ddot = (pp%dphidot - pp%dphi) * pp%phidot / pp%phi
      e = e - (dlog + l + 1) / ddot
      if (abs(dlog + l + 1) < 1.0e-13_dp) exit
    end do
    ec_exact = e
    ! parameter C at three E_nu moving toward C_exact
    enu = [ec_exact + 0.4_dp, ec_exact + 0.2_dp, ec_exact + 0.1_dp]
    do k = 1, 3
      call potpar_at_enu(m, v, l, enu(k), z, c_ry, avw, pp)
      err(k) = abs(pp%c - ec_exact)
    end do
    print '(a,i1,a,f10.5,a,3es11.3)', '    l=', l, '  exact C=', ec_exact, '  |C(Enu)-C_exact| at dE=0.4,0.2,0.1: ', err
    call check(err(1) > err(2) .and. err(2) > err(3), 'C error decreases as E_nu -> C')
    ! second order: halving (E_nu - C) should reduce the error ~4x (allow 2.5..8)
    call check(err(1) / err(2) > 2.5_dp .and. err(1) / err(2) < 9.0_dp, 'C error ~ quadratic in (E_nu - C)')
    call potpar_at_enu(m, v, l, ec_exact, z, c_ry, avw, pp)
    call check_close(pp%c, ec_exact, 1.0e-9_dp, 'C reproduces itself at E_nu = C_exact')
    call check(pp%sqrdel /= 0.0_dp .and. pp%p > 0.0_dp, 'sqrt(Delta) nonzero, p > 0')
    cgam = pp%gamma
  end do
  call report()
end program test_potpar
