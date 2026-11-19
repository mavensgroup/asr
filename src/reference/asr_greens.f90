!> Continued-fraction Green function and tails (terminators).
module asr_greens
  use asr_kinds, only: dp, pi
  implicit none
  private
  public :: periodic_g, sqrt_tail, tail_params_last, tail_params_avg2, cf_green, ldos_from_g

contains

  !> Square-root tail  t(z) = binf^2 * g(z),  g = 1/(z-ainf-binf^2 g).
  !> Branch: the root with |binf*g| <= 1 (decaying, physical sheet).
  pure function sqrt_tail(z, ainf, binf) result(t)
    complex(dp), intent(in) :: z
    real(dp), intent(in) :: ainf, binf
    complex(dp) :: t, w, s, g1, g2, g
    w = z - ainf
    s = sqrt(w * w - 4.0_dp * binf * binf)
    g1 = (w - s) / (2.0_dp * binf * binf)
    g2 = (w + s) / (2.0_dp * binf * binf)
    if (abs(g1) <= abs(g2)) then; g = g1; else; g = g2; end if
    t = binf * binf * g
  end function sqrt_tail

  !> Green function g_{n+1} of a semi-infinite chain whose levels repeat with
  !> period p:  levels a_p(1..p), couplings b_p(k) between level k and k+1 (b_p(p)
  !> links the last level of one period to the first of the next).  Fixed point
  !> of the composed Moebius maps g_k = 1/(z - a_k - b_k^2 g_{k+1}); the physical
  !> solution is the attracting one (eigenvector of the larger |eigenvalue|).
  !> p = 1 reduces to the square-root terminator.
  pure function periodic_g(z, ap, bp) result(g)
    complex(dp), intent(in) :: z
    real(dp), intent(in) :: ap(:), bp(:)
    complex(dp) :: g
    complex(dp) :: m(2, 2), mk(2, 2), tr, det, l1, l2, lam, disc
    integer :: k, p
    p = size(ap)
    m = reshape([(1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), (1.0_dp, 0.0_dp)], [2, 2])
    do k = 1, p
      mk(1, 1) = 0.0_dp; mk(1, 2) = 1.0_dp
      mk(2, 1) = -bp(k)**2; mk(2, 2) = z - ap(k)
      m = matmul(m, mk)            ! (x_1,y_1) = M_1 M_2 ... M_p (x_{p+1},y_{p+1})
    end do
    tr = m(1, 1) + m(2, 2); det = m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)
    disc = sqrt(tr * tr - 4.0_dp * det)
    l1 = 0.5_dp * (tr + disc); l2 = 0.5_dp * (tr - disc)
    if (abs(l1) >= abs(l2)) then; lam = l1; else; lam = l2; end if
    ! eigenvector (x,y) of m for lam:  (m11 - lam) x + m12 y = 0
    if (abs(m(1, 2)) > 0.0_dp .or. abs(m(1, 1) - lam) > 0.0_dp) then
      ! x/y = m12 / (lam - m11)  (if lam == m11, use second row)
      if (abs(lam - m(1, 1)) > abs(m(2, 1))) then
        g = m(1, 2) / (lam - m(1, 1))
      else
        g = (lam - m(2, 2)) / m(2, 1)
      end if
    else
      g = (lam - m(2, 2)) / m(2, 1)
    end if
  end function periodic_g

  !> Baseline tail parameters: last coefficient.
  pure subroutine tail_params_last(a, b, n, ainf, binf)
    real(dp), intent(in) :: a(:), b(:)
    integer, intent(in) :: n
    real(dp), intent(out) :: ainf, binf
    ainf = a(n); binf = b(n)
  end subroutine tail_params_last

  !> Baseline tail parameters: mean of the last two levels (kills the
  !> even/odd oscillation of bipartite lattices such as bcc).
  pure subroutine tail_params_avg2(a, b, n, ainf, binf)
    real(dp), intent(in) :: a(:), b(:)
    integer, intent(in) :: n
    real(dp), intent(out) :: ainf, binf
    ainf = 0.5_dp * (a(n) + a(n - 1))
    binf = sqrt(0.5_dp * (b(n)**2 + b(n - 1)**2))
  end subroutine tail_params_avg2

  !> G(z) = 1/(z-a1-b1^2/(z-a2- ... -b_{n-1}^2/(z-a_n - tail)))
  !> tail given as a complex value (use 0 for an exhausted Krylov space).
  pure function cf_green(a, b, n, z, tail) result(g)
    real(dp), intent(in) :: a(:), b(:)
    integer, intent(in) :: n
    complex(dp), intent(in) :: z, tail
    complex(dp) :: g
    integer :: i
    g = 1.0_dp / (z - a(n) - tail)
    do i = n - 1, 1, -1
      g = 1.0_dp / (z - a(i) - b(i) * b(i) * g)
    end do
  end function cf_green

  elemental function ldos_from_g(g) result(rho)
    complex(dp), intent(in) :: g
    real(dp) :: rho
    rho = -aimag(g) / pi
  end function ldos_from_g
end module asr_greens
