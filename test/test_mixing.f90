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

  ! after a restart the newest secant pair is the history: a kick far from the fixed point of a non-linear 1-d map restarts it
  ! (|r| grows), and the next step must be the secant step through the last two iterates (old ring bug: a stale pair was fitted)
  block
    type(anderson_t) :: b
    real(dp) :: y(1), yo(1), yn(1), yprev(1), rprev(1), sec
    integer :: k
    call anderson_init(b, 1, 3)
    y = 0.0_dp
    do k = 1, 8                                         ! fill the ring past its end so that the write position is not column 1
      yo = G1(y); call anderson_step(b, y, yo, 0.3_dp, yn); y = yn
    end do
    yprev = y + 0.01_dp; rprev = G1(yprev) - yprev
    yo = G1(yprev); call anderson_step(b, yprev, yo, 0.3_dp, yn)
    y = yprev + 3.0_dp                                  ! kick: the residual grows -> restart
    yo = G1(y); call anderson_step(b, y, yo, 0.3_dp, yn)
    sec = y(1) - (yo(1) - y(1)) * (y(1) - yprev(1)) / ((yo(1) - y(1)) - rprev(1))
    call check(b%nrestart >= 1 .and. abs(yn(1) - sec) < 1.0e-6_dp * (1.0_dp + abs(sec)), &
               'after a restart the step uses the newest pair (secant through the last two iterates)')
  end block
  call report()
contains
  function G1(y) result(f_)
    real(dp), intent(in) :: y(1)
    real(dp) :: f_(1)
    f_ = 1.0_dp - 4.0_dp * (y - 1.0_dp) - 0.5_dp * (y - 1.0_dp)**3
  end function G1

  function F(y) result(f_)
    real(dp), intent(in) :: y(4)
    real(dp) :: f_(4)
    f_ = tgt - g * (y - tgt) + 0.1_dp * sin(y - tgt)
  end function F
end program test_mixing
