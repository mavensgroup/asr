!> Tetrahedron integration: free-electron sphere (exact occupied volume), tight-binding half filling,
!> weights sum rule, and convergence with the grid.
program test_tetra
  use asr_kinds, only: dp, pi
  use asr_tetra
  use asr_testing
  implicit none
  integer :: nk, ix, iy, iz, ik, it
  integer, parameter :: nks(2) = [16, 32]
  real(dp), allocatable :: e(:, :), w(:, :)
  real(dp) :: kx, ky, kz, ef, exact, nn, err(2)
  print '(a)', 'test_tetra'
  ! free electrons e = k^2 inside [-pi,pi)^3 (non-periodic band, fine for E << pi^2)
  do it = 1, 2
    nk = nks(it)
    allocate (e(1, nk**3), w(1, nk**3))
    do iz = 0, nk - 1; do iy = 0, nk - 1; do ix = 0, nk - 1
          ik = 1 + ix + nk * (iy + nk * iz)
          kx = 2 * pi * (modulo(ix + nk / 2, nk) - nk / 2) / nk
          ky = 2 * pi * (modulo(iy + nk / 2, nk) - nk / 2) / nk
          kz = 2 * pi * (modulo(iz + nk / 2, nk) - nk / 2) / nk
          e(1, ik) = kx * kx + ky * ky + kz * kz
    end do; end do; end do
    ef = 2.0_dp
    call tetra_weights(nk, e, ef, w)
    exact = (4.0_dp * pi / 3.0_dp) * ef**1.5_dp / (2.0_dp * pi)**3       ! occupied fraction of the BZ (volume (2 pi)^3)
    nn = sum(w)
    err(it) = abs(nn - exact) / exact
    print '(a,i3,a,f10.6,a,f10.6,a,es9.2)', '    free electrons nk=', nk, '  N=', nn, '  exact=', exact, '  rel.err=', err(it)
    deallocate (e, w)
  end do
  ! linear tetrahedra on a strongly curved small sphere: second-order convergence (error x4 per doubling)
  call check(err(1) / err(2) > 3.0_dp .and. err(1) / err(2) < 5.5_dp, 'second-order convergence of the occupied volume')
  call check(err(2) < 2.0e-2_dp, 'free-electron occupied volume within 2% at nk=32 (sphere only ~4 spacings wide)')
  ! tight-binding sc band at half filling: exactly half the states
  nk = 12
  allocate (e(1, nk**3), w(1, nk**3))
  do iz = 0, nk - 1; do iy = 0, nk - 1; do ix = 0, nk - 1
        ik = 1 + ix + nk * (iy + nk * iz)
        e(1, ik) = -2.0_dp * (cos(2 * pi * ix / nk + 0.0_dp) + cos(2 * pi * iy / nk) + cos(2 * pi * iz / nk))
  end do; end do; end do
  call tetra_weights(nk, e, 0.0_dp, w)
  call check_close(sum(w), 0.5_dp, 1.0e-12_dp, 'sc tight-binding band: half filling -> 0.5')
  call tetra_weights(nk, e, 100.0_dp, w)
  call check_close(sum(w), 1.0_dp, 1.0e-12_dp, 'E above the band -> full')
  call tetra_weights(nk, e, -100.0_dp, w)
  call check_close(sum(w), 0.0_dp, 1.0e-12_dp, 'E below the band -> empty')
  ! Fermi level from the electron count: quarter filling of two spin-degenerate bands
  block
    real(dp) :: efq, f(1, 1)
    call tetra_fermi(nk, e, 0.5_dp, 1.0_dp, efq)
    call check_close(efq, 0.0_dp, 1.0e-9_dp, 'tetra_fermi: half filling gives E_F = 0')
    f = 0.0_dp
  end block
  ! tetra_fermi (closed-form counting) must agree with summing the Bloechl weights
  block
    real(dp) :: efc, ec(2, nk**3), wc(2, nk**3)
    ec(1, :) = e(1, :); ec(2, :) = 0.7_dp * e(1, :) + 0.3_dp
    call tetra_fermi(nk, ec, 0.8_dp, 1.0_dp, efc)
    call tetra_weights(nk, ec, efc, wc)
    call check_close(sum(wc), 0.8_dp, 1.0e-10_dp, 'tetra_fermi consistent with the weights sum (two bands)')
  end block
  call report()
end program test_tetra
