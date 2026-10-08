!> Contour moments from a terminated continued fraction vs the tetrahedron method (independent), for the
!> ordered sc tight-binding band:  M0 = occupied states, M1 = int E N dE up to E_F.
program test_contour
  use asr_kinds, only: dp, pi
  use asr_operator, only: csr_op
  use asr_lattice_toy, only: cubic_cluster
  use asr_recursion, only: lanczos
  use asr_greens, only: cf_green, sqrt_tail, tail_params_avg2
  use asr_contour
  use asr_tetra, only: tetra_weights
  use asr_testing
  implicit none
  integer, parameter :: nmax = 60, nk = 40, nz = 80
  type(csr_op) :: op
  real(dp) :: a(nmax), b(nmax), ainf, binf, ef, m0, m1, t0, t1
  real(dp), allocatable :: v(:), e(:, :), w(:, :)
  complex(dp) :: zc(nz), wc(nz), g(nz)
  integer :: centre, m, ix, iy, iz, ik, j, ie
  real(dp), parameter :: pi2 = 2.0_dp * pi
  real(dp) :: efs(3) = [-3.0_dp, -1.0_dp, 0.0_dp]
  print '(a)', 'test_contour'
  call cubic_cluster('sc', 2 * nmax + 6, 1.0_dp, op, centre)
  allocate (v(op%n)); v = 0.0_dp; v(centre) = 1.0_dp
  call lanczos(op, v, nmax, a, b, m)
  call tail_params_avg2(a, b, nmax, ainf, binf)
  allocate (e(1, nk**3), w(1, nk**3))
  do iz = 0, nk - 1; do iy = 0, nk - 1; do ix = 0, nk - 1
        ik = 1 + ix + nk * (iy + nk * iz)
        e(1, ik) = -2.0_dp * (cos(pi2 * ix / nk) + cos(pi2 * iy / nk) + cos(pi2 * iz / nk))
  end do; end do; end do
  do ie = 1, 3
    ef = efs(ie)
    call contour_nodes(-7.5_dp, ef, nz, zc, wc)
    do j = 1, nz
      g(j) = cf_green(a, b, nmax, zc(j), sqrt_tail(zc(j), ainf, binf))
    end do
    m0 = contour_moment(zc, wc, g, 0.0_dp, 0)
    m1 = contour_moment(zc, wc, g, 0.0_dp, 1)
    ! tetrahedron reference: occupied fraction and sum of w*E
    call tetra_weights(nk, e, ef, w)
    t0 = sum(w); t1 = sum(w(1, :) * e(1, :))
    print '(a,f6.2,a,2f10.5,a,2f10.5)', '    E_F=', ef, '   contour M0,M1 =', m0, m1, '   tetra =', t0, t1
    call check_close(m0, t0, 4.0e-3_dp, 'M0 (occupied states) contour vs tetra')
    call check_close(m1, t1, 6.0e-3_dp, 'M1 (energy moment) contour vs tetra')
  end do
  call report()
end program test_contour
