!> Structure constants: analytic limits, symmetry, and a gate against the
!> legacy INFO_FCC on-site block (Stuttgart TBSTR, sigma=0.7, from the 2009 code).
program test_strconst
  use asr_kinds, only: dp
  use asr_strconst
  use asr_testing
  implicit none
  real(dp) :: b(9, 9), bt(9, 9), w, d(3), s, al(0:3)
  real(dp), allocatable :: vec(:, :), blk(:, :, :)
  real(dp), parameter :: ref(9) = [3.2326_dp, 3.0223_dp, 3.0223_dp, 3.0223_dp, 2.7422_dp, &
                                   2.7422_dp, 1.9492_dp, 2.7422_dp, 1.9492_dp]
  integer :: i

  print '(a)', 'test_strconst'
  w = avg_ws_radius('fcc')
  d = [0.5_dp, 0.5_dp, 0.0_dp]
  call bare_s0_block(d, w, b)
  call check_close(b(1, 1), -2.0_dp * w / sqrt(0.5_dp), 1.0e-13_dp, 'S0_ss = -2 w/d')
  call bare_s0_block(-d, w, bt)
  s = maxval(abs(b - transpose(bt)))
  call check_close(s, 0.0_dp, 1.0e-12_dp, 'S0_{LL''}(d) = S0_{L''L}(-d)')
  al = hardcore_alpha(0.7_dp)
  call check_close(al(0) + al(1) + al(2), 0.35_dp + 0.7_dp**3 / 6 + 0.7_dp**5 / 10, 1.0e-14_dp, 'alpha_l formula')

  call lattice_cluster('fcc', 1.55_dp, vec)
  print '(a,i0)', '    fcc cluster sites within 1.55 alat: ', size(vec, 2)
  call screened_cluster(vec, w, hardcore_alpha(0.7_dp), blk)
  s = 0.0_dp
  do i = 1, 9
    print '(a,i2,a,f9.4,a,f9.4)', '    onsite diag ', i, '  ours=', blk(i, i, 1), '   legacy=', ref(i)
    s = max(s, abs(blk(i, i, 1) - ref(i)))
  end do
  call check_close(s, 0.0_dp, 2.0e-3_dp, 'fcc screened on-site block == legacy INFO_FCC (sigma=0.7)')
  call report()
end program test_strconst
