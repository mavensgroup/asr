!> ASA potential parameters (kappa^2 = 0) from the radial solutions at the sphere radius S:
!> energy E_nu, C, sqrt(Delta), p, gamma in the gamma representation (alpha = gamma).
!> Formulas follow the TB-LMTO notation (Wronskians of the bare K, J functions with the
!> normalised phi and its energy derivative phidot, <phi|phidot> = 0):
!>   W{K,phi}    = (D_phi    + l + 1) phi(S)
!>   W{K,phidot} = (D_phidot + l + 1) phidot(S)
!>   W{J,phidot} = (D_phidot - l) phidot(S) / (2 (2l+1))
!>   C = E_nu - W{K,phi}/W{K,phidot}
!>   sqrt(Delta) = - sqrt( (S/w)^(2l+1) / (2 S) ) / W{K,phidot}
!>   gamma = (S/w)^(2l+1) W{J,phidot} / W{K,phidot},   p = <phidot|phidot>
module asr_potpar
  use asr_kinds, only: dp
  use asr_radial
  implicit none
  private
  public :: potpar_l_t, potpar_at_enu

  type :: potpar_l_t
    real(dp) :: enu = 0.0_dp, c = 0.0_dp, sqrdel = 0.0_dp, p = 0.0_dp, gamma = 0.0_dp
    real(dp) :: dphi = 0.0_dp, dphidot = 0.0_dp, phi = 0.0_dp, phidot = 0.0_dp   !! values/log-derivs at S
  end type potpar_l_t

contains

  !> m: mesh ending exactly at S (odd number of points); v: true potential on the mesh;
  !> avw: average WS radius (Bohr).  The radial equation is integrated outward at E = enu.
  subroutine potpar_at_enu(m, v, l, enu, z, c, avw, pp, y)
    type(radmesh_t), intent(in) :: m
    real(dp), intent(in) :: v(:), enu, z, c, avw
    integer, intent(in) :: l
    type(potpar_l_t), intent(out) :: pp
    real(dp), intent(out), optional :: y(2, 0:2, m%n)
    real(dp), allocatable :: yy(:, :, :), f(:)
    real(dp) :: s, mm, wk, wkdot, wjdot, sdivw, pc, nrm
    integer :: n, i
    n = m%n
    allocate (yy(2, 0:2, n), f(n))
    s = m%r(n)
    call integrate_regular(m, v, l, enu, z, c, 2, n, yy)
    call norm_pq(m, c, 2, n, yy, nrm)
    mm = 1.0_dp + (enu - v(n)) / (c * c)
    pp%phi = yy(1, 0, n) / s
    pp%phidot = yy(1, 1, n) / s
    pp%dphi = mm * yy(2, 0, n) * s / yy(1, 0, n)
    pp%dphidot = s * (mm * yy(2, 1, n) + yy(2, 0, n) / c**2) / yy(1, 1, n)
    ! p = <phidot|phidot> over the sphere
    f = [(yy(1, 1, i)**2 + yy(2, 1, i)**2 / c**2, i=1, n)] * m%drdt
    pc = f(1) + f(n)
    do i = 2, n - 1
      if (mod(i, 2) == 0) then; pc = pc + 4.0_dp * f(i); else; pc = pc + 2.0_dp * f(i); end if
    end do
    pp%p = pc / 3.0_dp
    sdivw = (s / avw)**(2 * l + 1)
    wk = (pp%dphi + l + 1) * pp%phi
    wkdot = (pp%dphidot + l + 1) * pp%phidot
    wjdot = (pp%dphidot - l) / (2.0_dp * (2 * l + 1)) * pp%phidot
    pp%enu = enu
    pp%c = enu - wk / wkdot
    pp%sqrdel = -sqrt(0.5_dp * sdivw / s) / wkdot
    pp%gamma = sdivw * wjdot / wkdot
    if (present(y)) y = yy
  end subroutine potpar_at_enu
end module asr_potpar
