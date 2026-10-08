!> Real spherical harmonics and Gaunt coefficients.
!> Index order within an l-shell is m = -l..l (so s | y z x | xy yz 3z2-r2 xz x2-y2),
!> the Stuttgart/INFO order.  Gaunt coefficients are computed by exact
!> quadrature, so they are independent of any tabulated Clebsch-Gordan data.
module asr_harmonics
  use asr_kinds, only: dp, pi
  implicit none
  private
  public :: ylm_real, gaunt_table, lm_idx

contains

  pure integer function lm_idx(l, m)
    integer, intent(in) :: l, m
    lm_idx = l * l + l + m + 1        ! 1-based; s=1, y=2, z=3, x=4, ...
  end function lm_idx

  !> Real orthonormal Y_lm(rhat) for l = 0..lmax; y(lm_idx(l,m)).
  pure subroutine ylm_real(lmax, x, y, z, yl)
    integer, intent(in) :: lmax
    real(dp), intent(in) :: x, y, z          !! need not be normalised
    real(dp), intent(out) :: yl((lmax + 1)**2)
    real(dp) :: r, ct, st, phi, p(0:lmax, 0:lmax), nrm, cm, sm
    integer :: l, m
    r = sqrt(x * x + y * y + z * z)
    ct = z / r
    st = sqrt(max(0.0_dp, 1.0_dp - ct * ct))
    phi = atan2(y, x)
    ! associated Legendre P_l^m(ct) without Condon-Shortley phase
    p = 0.0_dp
    p(0, 0) = 1.0_dp
    do m = 1, lmax
      p(m, m) = p(m - 1, m - 1) * (2 * m - 1) * st
    end do
    do m = 0, lmax - 1
      p(m + 1, m) = ct * (2 * m + 1) * p(m, m)
    end do
    do m = 0, lmax
      do l = m + 2, lmax
        p(l, m) = (ct * (2 * l - 1) * p(l - 1, m) - (l + m - 1) * p(l - 2, m)) / real(l - m, dp)
      end do
    end do
    do l = 0, lmax
      do m = 0, l
        nrm = sqrt((2 * l + 1) / (4.0_dp * pi) * fact_ratio(l, m))
        if (m == 0) then
          yl(lm_idx(l, 0)) = nrm * p(l, 0)
        else
          cm = cos(m * phi); sm = sin(m * phi)
          yl(lm_idx(l, m)) = sqrt(2.0_dp) * nrm * p(l, m) * cm
          yl(lm_idx(l, -m)) = sqrt(2.0_dp) * nrm * p(l, m) * sm
        end if
      end do
    end do
  end subroutine ylm_real

  !> (l-m)!/(l+m)!
  pure real(dp) function fact_ratio(l, m)
    integer, intent(in) :: l, m
    integer :: k
    fact_ratio = 1.0_dp
    do k = l - m + 1, l + m
      fact_ratio = fact_ratio / real(k, dp)
    end do
  end function fact_ratio

  !> G(L, L', L'') = int Y_L Y_L' Y_L'' dOmega for all L, L' (l<=lmx) and
  !> L'' (l''<=2 lmx), exact Gauss-Legendre x trapezoid quadrature.
  subroutine gaunt_table(lmx, g)
    integer, intent(in) :: lmx
    real(dp), allocatable, intent(out) :: g(:, :, :)
    integer :: nth, nph, it, ip, n1, n3, i, j, k
    real(dp), allocatable :: xg(:), wg(:)
    real(dp) :: ct, st, ph, w, yv((2 * lmx + 1)**2)
    n1 = (lmx + 1)**2; n3 = (2 * lmx + 1)**2
    allocate (g(n1, n1, n3)); g = 0.0_dp
    nth = 3 * lmx + 3; nph = 6 * lmx + 6
    allocate (xg(nth), wg(nth))
    call gauss_legendre(nth, xg, wg)
    do it = 1, nth
      ct = xg(it); st = sqrt(1.0_dp - ct * ct)
      do ip = 0, nph - 1
        ph = 2.0_dp * pi * ip / nph
        w = wg(it) * 2.0_dp * pi / nph
        call ylm_real(2 * lmx, st * cos(ph), st * sin(ph), ct, yv)
        do k = 1, n3
          do j = 1, n1
            do i = 1, n1
              g(i, j, k) = g(i, j, k) + w * yv(i) * yv(j) * yv(k)
            end do
          end do
        end do
      end do
    end do
    where (abs(g) < 1.0e-13_dp) g = 0.0_dp
  end subroutine gaunt_table

  subroutine gauss_legendre(n, x, w)
    integer, intent(in) :: n
    real(dp), intent(out) :: x(n), w(n)
    integer :: i, k, it
    real(dp) :: z, z1, p1, p2, p3, pp
    do i = 1, (n + 1) / 2
      z = cos(pi * (i - 0.25_dp) / (n + 0.5_dp))
      do it = 1, 100
        p1 = 1.0_dp; p2 = 0.0_dp
        do k = 1, n
          p3 = p2; p2 = p1
          p1 = ((2 * k - 1) * z * p2 - (k - 1) * p3) / k
        end do
        pp = n * (z * p1 - p2) / (z * z - 1.0_dp)
        z1 = z; z = z1 - p1 / pp
        if (abs(z - z1) < 1.0e-15_dp) exit
      end do
      x(i) = -z; x(n + 1 - i) = z
      w(i) = 2.0_dp / ((1.0_dp - z * z) * pp * pp); w(n + 1 - i) = w(i)
    end do
  end subroutine gauss_legendre
end module asr_harmonics
