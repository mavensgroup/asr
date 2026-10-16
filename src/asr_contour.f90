!> Complex-energy contour integration of Green functions: occupied moments
!>   M_k = -(1/pi) Im  int_C (z - E_nu)^k G(z) dz ,   C = semicircle in the upper half plane from E_bot to E_F,
!> which equals  int_{E_bot}^{E_F} (E - E_nu)^k N(E) dE  when the Green function is analytic there (as the
!> terminated continued fraction is).  Gauss-Legendre in the angle.
module asr_contour
  use asr_kinds, only: dp, pi
  implicit none
  private
  public :: contour_nodes, contour_nodes_closed, contour_moment, contour_nodes_fermi, fermi_npole, fermi_nodes_count

contains

  !> z(j), w(j): int_C f(z) dz ~ sum_j w(j) f(z(j)), C: E_bot -> E_F through the upper half plane.
  subroutine contour_nodes(ebot, ef, n, z, w, eta_end)
    real(dp), intent(in) :: ebot, ef
    real(dp), intent(in), optional :: eta_end         !! stop the arc at Im z = eta_end (smooth, monotone N(E_F)); default: real axis
    integer, intent(in) :: n
    complex(dp), intent(out) :: z(n), w(n)
    real(dp) :: x(n), wg(n), c, r, th, th0
    integer :: j
    call gauss_legendre(n, x, wg)
    c = 0.5_dp * (ebot + ef); r = 0.5_dp * (ef - ebot)
    th0 = 0.0_dp
    if (present(eta_end)) th0 = asin(min(eta_end / r, 0.5_dp))
    do j = 1, n
      th = th0 + 0.5_dp * (pi - th0) * (1.0_dp - x(j))   ! theta from pi (at E_bot) to th0 (near E_F)
      z(j) = cmplx(c + r * cos(th), r * sin(th), dp)
      ! dz = i r e^{i th} dth ; dth = -(pi/2) dx  (th decreases with x)
      w(j) = wg(j) * (-0.5_dp * (pi - th0)) * cmplx(0.0_dp, r, dp) * exp(cmplx(0.0_dp, th, dp))
    end do
  end subroutine contour_nodes

  !> Path-independent version: circular arc from E_bot to E_F + i eta (n_arc nodes) closed by the vertical segment
  !> E_F + i eta -> E_F (n_vert nodes, all strictly in the upper half plane).  The sum equals the real-axis integral
  !> to E_F by Cauchy for any analytic integrand, so N(E_F) does not depend on eta.
  subroutine contour_nodes_closed(ebot, ef, eta, n_arc, n_vert, z, w)
    real(dp), intent(in) :: ebot, ef, eta
    integer, intent(in) :: n_arc, n_vert
    complex(dp), intent(out) :: z(n_arc + n_vert), w(n_arc + n_vert)
    real(dp) :: xa(n_arc), wa(n_arc), xv(n_vert), wv(n_vert), c, r, th, th0
    integer :: j
    call gauss_legendre(n_arc, xa, wa)
    call gauss_legendre(n_vert, xv, wv)
    c = (ef * ef + eta * eta - ebot * ebot) / (2.0_dp * (ef - ebot)); r = c - ebot      ! circle through E_bot and E_F + i eta
    th0 = atan2(eta, ef - c)
    do j = 1, n_arc
      th = th0 + 0.5_dp * (pi - th0) * (1.0_dp - xa(j))
      z(j) = cmplx(c + r * cos(th), r * sin(th), dp)
      w(j) = wa(j) * (-0.5_dp * (pi - th0)) * cmplx(0.0_dp, r, dp) * exp(cmplx(0.0_dp, th, dp))
    end do
    do j = 1, n_vert
      z(n_arc + j) = cmplx(ef, 0.5_dp * eta * (1.0_dp + xv(j)), dp)
      w(n_arc + j) = cmplx(0.0_dp, -0.5_dp * eta * wv(j), dp)
    end do
  end subroutine contour_nodes_closed

  !> Matsubara poles below the Fermi contour: Gamma = 2 npole pi kT ~ 0.1 Ry (G is smooth there), 2 <= npole <= 40.
  pure integer function fermi_npole(kt)
    real(dp), intent(in) :: kt
    fermi_npole = max(2, min(40, nint(0.1_dp / (2.0_dp * pi * kt))))
  end function fermi_npole

  pure integer function fermi_nodes_count(npole, n_main, n_step)
    integer, intent(in) :: npole, n_main, n_step
    fermi_nodes_count = 4 + n_main + n_step + npole
  end function fermi_nodes_count

  !> Fermi-smeared contour (Wildberger, Lang, Zeller, Dederichs, PRB 52, 11502 (1995)):
  !>   int f(E) [-(1/pi) Im h(E + i0)] dE  =  -(1/pi) Im sum_j w(j) h(z(j))      (h analytic in the upper half plane)
  !> f(E) = 1 / (exp((E - mu)/kT) + 1).  Path: E_bot -> E_bot + i Gamma (4 nodes) -> mu - 10 kT + i Gamma (n_main) -> mu + 40 kT + i Gamma
  !> (n_step, where the occupation falls); Gamma = 2 npole pi kT, on which f(x + i Gamma) = f(x) is real.  The weights carry f; the
  !> leg back to the real axis beyond mu + 40 kT is dropped (f < 5e-18).  The npole Matsubara poles mu + i (2n+1) pi kT enclosed
  !> between path and real axis enter as nodes with w = -2 pi i kT.  N(mu) is smooth and monotone (no discrete steps at E_F).
  subroutine contour_nodes_fermi(ebot, mu, kt, npole, n_main, n_step, z, w)
    real(dp), intent(in) :: ebot, mu, kt
    integer, intent(in) :: npole, n_main, n_step
    complex(dp), intent(out) :: z(:), w(:)
    real(dp) :: xv(4), wv(4), xm(n_main), wm(n_main), xs(n_step), ws(n_step), gam, lo, hi, x
    integer :: j, k
    gam = 2.0_dp * npole * pi * kt
    call gauss_legendre(4, xv, wv); call gauss_legendre(n_main, xm, wm); call gauss_legendre(n_step, xs, ws)
    k = 0
    do j = 1, 4                                        ! vertical leg at E_bot (far below mu: f = 1)
      k = k + 1; z(k) = cmplx(ebot, 0.5_dp * gam * (1.0_dp + xv(j)), dp)
      w(k) = cmplx(0.0_dp, 0.5_dp * gam * wv(j), dp) * fermi_c(z(k))
    end do
    lo = ebot; hi = mu - 10.0_dp * kt
    do j = 1, n_main
      k = k + 1; x = 0.5_dp * (lo + hi) + 0.5_dp * (hi - lo) * xm(j); z(k) = cmplx(x, gam, dp)
      w(k) = 0.5_dp * (hi - lo) * wm(j) * fermi_r(x)
    end do
    lo = mu - 10.0_dp * kt; hi = mu + 40.0_dp * kt
    do j = 1, n_step
      k = k + 1; x = 0.5_dp * (lo + hi) + 0.5_dp * (hi - lo) * xs(j); z(k) = cmplx(x, gam, dp)
      w(k) = 0.5_dp * (hi - lo) * ws(j) * fermi_r(x)
    end do
    do j = 0, npole - 1
      k = k + 1; z(k) = cmplx(mu, (2 * j + 1) * pi * kt, dp); w(k) = cmplx(0.0_dp, -2.0_dp * pi * kt, dp)
    end do
  contains
    real(dp) function fermi_r(e)
      real(dp), intent(in) :: e
      real(dp) :: t
      t = (e - mu) / kt
      if (t > 0.0_dp) then; fermi_r = exp(-t) / (1.0_dp + exp(-t)); else; fermi_r = 1.0_dp / (1.0_dp + exp(t)); end if
    end function fermi_r
    complex(dp) function fermi_c(zz)                   ! complex Fermi function (only used far below mu, where it is ~1)
      complex(dp), intent(in) :: zz
      complex(dp) :: t
      t = (zz - mu) / kt
      if (real(t) < -600.0_dp) then; fermi_c = (1.0_dp, 0.0_dp); else; fermi_c = 1.0_dp / (exp(t) + 1.0_dp); end if
    end function fermi_c
  end subroutine contour_nodes_fermi

  !> M_k from G sampled on the contour nodes.
  pure function contour_moment(z, w, g, enu, k) result(m)
    complex(dp), intent(in) :: z(:), w(:), g(:)
    real(dp), intent(in) :: enu
    integer, intent(in) :: k
    real(dp) :: m
    m = -aimag(sum(w * (z - enu)**k * g)) / pi
  end function contour_moment

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
          p3 = p2; p2 = p1; p1 = ((2 * k - 1) * z * p2 - (k - 1) * p3) / k
        end do
        pp = n * (z * p1 - p2) / (z * z - 1.0_dp)
        z1 = z; z = z1 - p1 / pp
        if (abs(z - z1) < 1.0e-15_dp) exit
      end do
      x(i) = -z; x(n + 1 - i) = z
      w(i) = 2.0_dp / ((1.0_dp - z * z) * pp * pp); w(n + 1 - i) = w(i)
    end do
  end subroutine gauss_legendre
end module asr_contour
