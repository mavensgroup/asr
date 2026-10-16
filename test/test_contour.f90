!> Contour moments from a terminated continued fraction vs the tetrahedron method (independent), for the
!> ordered sc tight-binding band:  M0 = occupied states, M1 = int E N dE up to E_F.
!> Fermi-smeared contour: a single level h = 1/(z - eps) is occupied by f(eps); the smeared M0 of the sc band equals the
!> tetrahedron weights convolved with f (Fermi-Dirac occupation of the k-mesh levels).
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
  ! Fermi-smeared contour
  block
    real(dp) :: kt, mu, eps, occ, f0, t0f
    complex(dp), allocatable :: zf(:), wf(:), gf(:)
    integer :: np, nf, k2
    kt = 0.003_dp; mu = 0.0_dp
    np = fermi_npole(kt); nf = fermi_nodes_count(np, 40, 24)
    allocate (zf(nf), wf(nf), gf(nf))
    call contour_nodes_fermi(-1.5_dp, mu, kt, np, 40, 24, zf, wf)
    print '(a,i0,a,i0)', '    Fermi contour: Matsubara poles ', np, ', nodes ', nf
    do k2 = -3, 3
      eps = mu + 2.0_dp * kt * k2
      occ = -aimag(sum(wf / (zf - eps))) / pi
      f0 = 1.0_dp / (exp((eps - mu) / kt) + 1.0_dp)
      call check_close(occ, f0, 1.0e-6_dp, 'single level occupied by f(eps)')
    end do
    deallocate (zf, wf, gf)
    do ie = 1, 3
      mu = efs(ie)
      np = fermi_npole(kt); nf = fermi_nodes_count(np, 80, 24)
      allocate (zf(nf), wf(nf), gf(nf))
      call contour_nodes_fermi(-7.5_dp, mu, kt, np, 80, 24, zf, wf)
      do j = 1, nf
        gf(j) = cf_green(a, b, nmax, zf(j), sqrt_tail(zf(j), ainf, binf))
      end do
      m0 = contour_moment(zf, wf, gf, 0.0_dp, 0)
      call tetra_weights(nk, e, mu, w); t0 = sum(w)
      t0f = sum(1.0_dp / (exp(max(-600.0_dp, min(600.0_dp, (e(1, :) - mu) / kt))) + 1.0_dp)) / real(nk**3, dp)
      print '(a,f6.2,a,f10.5,a,2f10.5)', '    mu=', mu, '   Fermi-contour M0 =', m0, '   tetra (T=0), k-mesh Fermi-Dirac =', t0, t0f
      call check_close(m0, t0, 4.0e-3_dp, 'smeared M0 vs tetra (kT = 3 mRy)')
      deallocate (zf, wf, gf)
    end do
  end block
  call report()
end program test_contour
