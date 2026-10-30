!> asr_sro: Warren-Cowley first-shell tables -> conditional probabilities P(s | c) = y_s (1 - alpha_cs).
!>  alpha = 0 -> P = y; one-sublattice ternary: derived diagonal, rows sum to 1, homogeneity sum_c x_c P(s|c) = x_s, detailed
!>  balance x_c P(s|c) = x_s P(c|s); a pair with alpha = 1 never touches (P = 0); input errors (asymmetric, diagonal given,
!>  alpha > 1, collective budget, two-sublattice marginals).
program test_sro
  use asr_kinds, only: dp
  use asr_sro, only: sro_conditional
  use asr_testing
  implicit none
  real(dp) :: x(3), y(2), a(3, 3), b(2, 2), c2(3, 2), h(3)
  real(dp), allocatable :: p(:, :)
  character(len=:), allocatable :: msg
  integer :: s

  print '(a)', 'test_sro'
  x = [0.66_dp, 0.14_dp, 0.20_dp]
  a = 0.0_dp
  call sro_conditional(x, x, a, .true., p, msg)
  call check(len(msg) == 0 .and. maxval(abs(p - spread(x, 1, 3))) < 1.0e-14_dp, 'alpha = 0: P(s | c) = x_s')

  a = 0.0_dp                                                   ! Fe-Ni -0.1, Fe-Cr 0.2, Ni-Cr 0.3
  a(1, 2) = -0.1_dp; a(2, 1) = -0.1_dp; a(1, 3) = 0.2_dp; a(3, 1) = 0.2_dp; a(2, 3) = 0.3_dp; a(3, 2) = 0.3_dp
  call sro_conditional(x, x, a, .true., p, msg)
  call check(len(msg) == 0, 'ternary: accepted')
  call check_close(maxval(abs(sum(p, 2) - 1.0_dp)), 0.0_dp, 1.0e-14_dp, 'ternary: rows sum to 1')
  do s = 1, 3; h(s) = sum(x * p(:, s)); end do
  call check_close(maxval(abs(h - x)), 0.0_dp, 1.0e-14_dp, 'ternary: homogeneity sum_c x_c P(s|c) = x_s')
  call check_close(maxval(abs(spread(x, 2, 3) * p - transpose(spread(x, 2, 3) * p))), 0.0_dp, 1.0e-14_dp, &
                   'ternary: detailed balance x_c P(s|c) = x_s P(c|s)')
  call check_close(p(1, 1), x(1) + x(2) * (-0.1_dp) + x(3) * 0.2_dp, 1.0e-14_dp, 'ternary: X_c^(c) = x_c + sum_s x_s alpha_cs')
  call check_close(a(2, 2), -(x(1) * (-0.1_dp) + x(3) * 0.3_dp) / x(2), 1.0e-14_dp, 'ternary: derived diagonal alpha_cc')

  a = 0.0_dp; a(2, 3) = 1.0_dp; a(3, 2) = 1.0_dp               ! Ni-Cr never neighbours (Alam triad 0-0-1)
  call sro_conditional(x, x, a, .true., p, msg)
  call check(len(msg) == 0 .and. p(2, 3) == 0.0_dp .and. p(3, 2) == 0.0_dp, 'alpha = 1: the pair never touches (P = 0)')

  a = 0.0_dp; a(1, 2) = 0.1_dp; a(2, 1) = 0.2_dp
  call sro_conditional(x, x, a, .true., p, msg)
  call check(len(msg) > 0, 'asymmetric alpha rejected')
  a = 0.0_dp; a(1, 1) = 0.1_dp
  call sro_conditional(x, x, a, .true., p, msg)
  call check(len(msg) > 0, 'diagonal alpha rejected (derived)')
  a = 0.0_dp; a(1, 2) = 1.2_dp; a(2, 1) = 1.2_dp
  call sro_conditional(x, x, a, .true., p, msg)
  call check(len(msg) > 0, 'alpha > 1 rejected (P < 0)')
  a = 0.0_dp; a(2, 3) = -6.0_dp; a(3, 2) = -6.0_dp             ! X_Ni^(Ni) = 0.14 - 0.2 * 6 < 0
  call sro_conditional(x, x, a, .true., p, msg)
  call check(len(msg) > 0, 'ordering beyond the collective budget rejected')

  y = [0.4_dp, 0.6_dp]                                         ! two sublattices: binary x binary, one free parameter
  ! rows: sum_s y_s b(c, s) = 0 and columns: sum_c x_c b(c, s) = 0 with x = y = (0.4, 0.6): b = t [[-3, 2], [2, -4/3]] / 10
  b = 0.1_dp * reshape([-3.0_dp, 2.0_dp, 2.0_dp, -4.0_dp / 3.0_dp], [2, 2])
  call sro_conditional(y, y, b, .false., p, msg)
  call check(len(msg) == 0, 'two sublattices: consistent table accepted')
  if (len(msg) == 0) then
    call check_close(maxval(abs(sum(p, 2) - 1.0_dp)), 0.0_dp, 1.0e-14_dp, 'two sublattices: rows sum to 1')
    call check_close(maxval(abs(matmul(y, p) - y)), 0.0_dp, 1.0e-14_dp, 'two sublattices: homogeneity')
  end if
  c2 = 0.0_dp; c2(1, 2) = 0.2_dp
  call sro_conditional(x, y, c2, .false., p, msg)
  call check(len(msg) > 0, 'two sublattices: table violating the marginals rejected')
  call report()
end program test_sro
