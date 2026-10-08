!> KPM on a periodic sc tight-binding lattice reproduces the exact k-space DOS (independent truth for disordered tests).
program test_kpm
  use asr_kinds, only: dp, pi
  use asr_operator, only: csr_op
  use asr_augspace, only: graph_t, binary_model_t
  use asr_lattices, only: periodic_graph, random_alloy_op
  use asr_kpm
  use asr_testing
  implicit none
  type(graph_t) :: g
  type(binary_model_t) :: mdl
  type(csr_op) :: op
  integer, parameter :: ne = 1261, nmom = 400, ns = 41
  real(dp) :: mu(0:nmom - 1), e(ne), rho(ne), rhob(ne), ref(ns), es(ns), kx, ky, kz, eps, err
  integer :: ie, ix, iy, iz, nk, is
  print '(a)', 'test_kpm'
  call periodic_graph('sc', 12, g)
  call check(g%n == 12**3, 'sc supercell has L^3 sites')
  mdl%x = 1.0_dp; mdl%epsa = 0.0_dp; mdl%epsb = 0.0_dp; mdl%tAA = -1.0_dp; mdl%tAB = -1.0_dp; mdl%tBB = -1.0_dp
  call random_alloy_op(g, mdl, 1, op)
  call kpm_moments(op, 0.0_dp, 6.3_dp, nmom, 60, 5, mu)
  do ie = 1, ne; e(ie) = -6.3_dp + 12.6_dp * (ie - 1) / (ne - 1); end do
  do is = 1, ns; es(is) = -6.0_dp + 12.0_dp * (is - 1) / (ns - 1); end do
  call kpm_dos(mu, 0.0_dp, 6.3_dp, e, rho)
  call lorentz_broaden(e, rho, 0.3_dp, rhob)
  ! exact: <1/(E+i eta - eps(k))> imaginary part on a fine k grid
  nk = 12; ref = 0.0_dp                    ! the periodic 12^3 supercell's own exact spectrum
  do iz = 0, nk - 1; do iy = 0, nk - 1; do ix = 0, nk - 1
        kx = 2 * pi * ix / nk; ky = 2 * pi * iy / nk; kz = 2 * pi * iz / nk
        eps = -2.0_dp * (cos(kx) + cos(ky) + cos(kz))
        do is = 1, ns
          ref(is) = ref(is) + (0.3_dp / pi) / ((es(is) - eps)**2 + 0.09_dp)
        end do
  end do; end do; end do
  ref = ref / real(nk, dp)**3
  err = 0.0_dp
  do is = 1, ns
    ie = 1 + nint((es(is) + 6.3_dp) / 12.6_dp * (ne - 1))
    err = max(err, abs(rhob(ie) - ref(is)))
  end do
  err = err / maxval(ref)
  print '(a,es10.3)', '    max relative error of broadened KPM DOS vs exact k-space: ', err
  call check(err < 2.5e-2_dp, 'KPM (12^3 sites, 60 random vectors) within 2.5% of the exact spectrum of the same supercell')
  call report()
end program test_kpm
