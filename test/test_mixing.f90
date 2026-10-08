!> Anderson mixing on a stiff nonlinear fixed point: x = F(x) with a large-gain (unstable under linear mixing) map.
program test_mixing
  use asr_kinds, only: dp
  use asr_mixing
  use asr_testing
  implicit none
  type(anderson_t) :: a
  real(dp) :: x(4), xo(4), xn(4), xlin(4), g(4)
  integer :: it, itl
  real(dp), parameter :: tgt(4) = [1.0_dp, -2.0_dp, 0.5_dp, 3.0_dp]
  print '(a)', 'test_mixing'
  ! F(x) = tgt - G (x - tgt) with gains 12, 8, 1.5, 0.3 (|1 - beta*(1+g)| > 1 for beta=0.3 on the first two)
  g = [12.0_dp, 8.0_dp, 1.5_dp, 0.3_dp]
  xlin = 0.0_dp
  do itl = 1, 200
    xlin = xlin + 0.3_dp * (F(xlin) - xlin)
    if (maxval(abs(xlin - tgt)) > 1.0e10_dp) exit
  end do
  call check(maxval(abs(xlin - tgt)) > 1.0e-3_dp, 'plain linear mixing does NOT converge here (stiff map)')
  call anderson_init(a, 4, 6)
  x = 0.0_dp
  do it = 1, 60
    xo = F(x)
    call anderson_step(a, x, xo, 0.3_dp, xn)
    x = xn
    if (maxval(abs(x - tgt)) < 1.0e-10_dp) exit
  end do
  print '(a,i0,a,es9.2)', '    Anderson iterations: ', it, '  error: ', maxval(abs(x - tgt))
  call check(maxval(abs(x - tgt)) < 1.0e-9_dp, 'Anderson mixing converges the stiff map')
  call check(it < 40, 'in fewer than 40 iterations')
  call report()
contains
  function F(y) result(f_)
    real(dp), intent(in) :: y(4)
    real(dp) :: f_(4)
    f_ = tgt - g * (y - tgt) + 0.1_dp * sin(y - tgt)
  end function F
end program test_mixing
