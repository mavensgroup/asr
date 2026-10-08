program test_harmonics
  use asr_kinds, only: dp, pi
  use asr_harmonics
  use asr_testing
  implicit none
  real(dp), allocatable :: g(:, :, :)
  real(dp) :: y(25), s
  integer :: i, j
  print '(a)', 'test_harmonics'
  call gaunt_table(2, g)
  call check_close(g(1, 1, 1), 1.0_dp / sqrt(4.0_dp * pi), 1.0e-13_dp, 'G(s,s,s) = 1/sqrt(4pi)')
  ! orthonormality: G(L,L',s) * sqrt(4pi) = delta_LL'
  s = 0.0_dp
  do j = 1, 9; do i = 1, 9
      s = max(s, abs(g(i, j, 1) * sqrt(4.0_dp * pi) - merge(1.0_dp, 0.0_dp, i == j)))
  end do; end do
  call check_close(s, 0.0_dp, 1.0e-13_dp, 'Y_L orthonormal (spd)')
  ! addition theorem  sum_m Y_lm(r)Y_lm(r) = (2l+1)/4pi
  call ylm_real(4, 0.3_dp, -0.5_dp, 0.8_dp, y)
  call check_close(sum(y(lm_idx(2, -2):lm_idx(2, 2))**2), 5.0_dp / (4.0_dp * pi), 1.0e-13_dp, 'addition theorem l=2')
  call check_close(sum(y(lm_idx(4, -4):lm_idx(4, 4))**2), 9.0_dp / (4.0_dp * pi), 1.0e-13_dp, 'addition theorem l=4')
  ! z-axis: only m=0 nonzero;  y(3)=z-like
  call ylm_real(2, 0.0_dp, 0.0_dp, 1.0_dp, y(1:9))
  call check_close(abs(y(lm_idx(1, 0))), sqrt(3.0_dp / (4.0_dp * pi)), 1.0e-13_dp, 'p_z on z axis')
  call check_close(abs(y(lm_idx(1, 1))) + abs(y(lm_idx(1, -1))), 0.0_dp, 1.0e-13_dp, 'p_x,p_y vanish on z axis')
  call report()
end program test_harmonics
