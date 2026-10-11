!> Terminator unit tests + the gapped-band stress test.
program test_terminators
  use asr_kinds, only: dp, pi
  use asr_operator, only: csr_op
  use asr_lattice_toy, only: cubic_cluster
  use asr_recursion, only: lanczos
  use asr_greens, only: sqrt_tail, periodic_g
  use asr_kspace, only: kspace_green_b2
  use asr_terminator
  use asr_testing, only: check
  use asr_endep, only: cf_species, ed_tail
  implicit none
  integer, parameter :: ne = 81
  integer :: n, k, centre, m, ie
  real(dp) :: a(60), b(60), ainf, binf, egrid(ne), e_avg, e_per, d
  real(dp), allocatable :: v(:)
  complex(dp) :: z, g1, g2, gref(ne), g
  type(csr_op) :: op
  type(terminator_t) :: t
  integer :: ni
  logical :: ok
  complex(dp) :: ca(6), cb(6), cc(6), fdum

  print '(a)', 'test_terminators'
  ! 1. p=1 periodic tail == sqrt terminator
  z = cmplx(1.3_dp, 0.2_dp, dp)
  g1 = periodic_g(z, [0.4_dp], [1.7_dp]); g2 = sqrt_tail(z, 0.4_dp, 1.7_dp) / 1.7_dp**2
  call check_c(abs(g1 - g2), 1.0e-13_dp, 'periodic_g(p=1) == sqrt terminator')
  ! 2. period-2 tail is its own fixed point: g = 1/(z-a1-b1^2/(z-a2-b2^2 g))
  g1 = periodic_g(z, [0.5_dp, -0.5_dp], [1.2_dp, 0.9_dp])
  g2 = 1.0_dp / (z - 0.5_dp - 1.2_dp**2 / (z + 0.5_dp - 0.9_dp**2 * g1))
  call check_c(abs(g1 - g2), 1.0e-13_dp, 'period-2 tail is a fixed point')
  call check_c(max(0.0_dp, aimag(g1)), 1.0e-14_dp, 'period-2 tail on physical sheet (Im g<0)')
  ! 3. Beer-Pettifor on a uniform chain recovers (0,1) as depth grows
  a = 0.0_dp; b = 1.0_dp
  do n = 8, 56, 16
    call bp_asymptote(a, b, n, ainf, binf)
    print '(a,i3,a,es10.3,a,es10.3)', '    BP uniform chain n=', n, '  |A|=', abs(ainf), '  |B-1|=', abs(binf - 1.0_dp)
  end do
  call check(abs(binf - 1.0_dp) < 1.0e-2_dp .and. abs(ainf) < 1.0e-12_dp, 'BP -> (0,1) on uniform chain, n=56')
  ! 4. Luchini-Nex joins linearly between last computed and asymptote
  a = 0.0_dp; b = 1.0_dp; a(8) = 0.4_dp; b(8) = 1.4_dp
  ni = 4
  t = make_terminator('sqrt-avg2', 'luchini-nex', a, b, 8, ni)
  call check(t%n == 12, 'LN extends the chain by ninterp levels')
  call check_c(abs(t%a(12) - (a(8) + 4.0_dp / 5.0_dp * (t%ainf - a(8)))), 1.0e-14_dp, 'LN a interpolation')
  call check_c(abs(t%b(10) - (b(8) + 2.0_dp / 5.0_dp * (t%binf - b(8)))), 1.0e-14_dp, 'LN b interpolation')
  ! 5. gapped B2 compound: period-2 tail must beat every one-band tail
  do ie = 1, ne
    egrid(ie) = -9.0_dp + 18.0_dp * (ie - 1) / real(ne - 1, dp)
  end do
  d = 3.0_dp
  gref = kspace_green_b2(1.0_dp, d, egrid, 0.1_dp, 96)
  call cubic_cluster('bcc', 2 * 32 + 6, 1.0_dp, op, centre, b2_delta=d)
  allocate (v(op%n)); v = 0.0_dp; v(centre) = 1.0_dp
  call lanczos(op, v, 32, a, b, m)
  e_avg = 0.0_dp; e_per = 0.0_dp
  do k = 1, 2
    if (k == 1) then; t = make_terminator('sqrt-avg2', 'butt', a, b, 32, 0)
    else; t = make_terminator('periodic2', 'butt', a, b, 32, 0); end if
    do ie = 1, ne
      z = cmplx(egrid(ie), 0.1_dp, dp)
      g = terminated_green(t, z)
      if (k == 1) then; e_avg = max(e_avg, abs(aimag(g - gref(ie))) / pi)
      else; e_per = max(e_per, abs(aimag(g - gref(ie))) / pi); end if
    end do
  end do
  print '(a,es10.3,a,es10.3)', '    B2 gapped, depth 32: sqrt-avg2 err=', e_avg, '   periodic2 err=', e_per
  ok = e_per < 1.0e-3_dp .and. e_per * 20.0_dp < e_avg
  call check(ok, 'B2: periodic2 < 1e-3 and >=20x better than one-band tail')
  ! 6. complex-symmetric chain tail (asr_endep%cf_species, ed_tail = 'bp') on a real chain == real Beer-Pettifor (bp_asymptote)
  a(1:6) = [0.10_dp, -0.05_dp, 0.20_dp, 0.00_dp, 0.15_dp, 0.05_dp]; b(1:6) = [0.90_dp, 1.10_dp, 0.95_dp, 1.05_dp, 1.00_dp, 1.00_dp]
  call bp_asymptote(a, b, 6, ainf, binf)
  ca = a(1:6); cb = b(1:6); cc = (0.0_dp, 0.0_dp)
  ed_tail = 'bp'
  d = 0.0_dp
  do ie = 1, 9
    z = cmplx(-2.4_dp + 0.6_dp * (ie - 1), 0.05_dp, dp)
    call cf_species(ca, cb, cc, 6, z, g1, fdum)
    g2 = 1.0_dp / (z - a(6) - sqrt_tail(z, ainf, binf))
    do k = 5, 1, -1
      g2 = 1.0_dp / (z - a(k) - b(k)**2 * g2)
    end do
    d = max(d, abs(g1 - g2))
  end do
  ed_tail = 'hull'
  call check_c(d, 1.0e-10_dp, 'cf_species bp tail == bp_asymptote on a real chain')
  call report_all()
contains
  subroutine check_c(val, tol, name)
    real(dp), intent(in) :: val, tol
    character(*), intent(in) :: name
    call check(val <= tol, name)
  end subroutine check_c
  subroutine report_all()
    use asr_testing, only: report
    call report()
  end subroutine report_all
end program test_terminators
