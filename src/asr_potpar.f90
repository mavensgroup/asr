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
  public :: potpar_l_t, potpar_at_enu, band_window, clamp_enu

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

  !> Wigner-Seitz window of the l band on the radial branch with `nodes` nodes inside the sphere:
  !>   ebot: D_l(E) = 0 (band bottom);  etop: phi_l(S) = 0 (band top, the next node enters at S).
  !> Past etop the parametrisation jumps to the next branch (sqrt(Delta) changes sign, C runs away); Questaal avoids this by
  !> floating pnu with a fixed principal quantum number.
  subroutine band_window(m, v, l, nodes, z, c, ebot, etop)
    type(radmesh_t), intent(in) :: m
    real(dp), intent(in) :: v(:), z, c
    integer, intent(in) :: l, nodes
    real(dp), intent(out) :: ebot, etop
    real(dp) :: elo, ehi, emid, eprev
    integer :: it
    elo = v(m%n) - 20.0_dp                      ! below the valence branches (semicore tops are above ~ -10 Ry)
    ehi = v(m%n) + 1.0_dp
    do while (count_nodes(ehi) <= nodes .and. ehi < 400.0_dp)
      ehi = v(m%n) + 2.0_dp * (ehi - v(m%n))
    end do
    etop = bisect_nodes(elo, ehi, nodes)
    eprev = elo
    if (nodes > 0) eprev = bisect_nodes(elo, etop, nodes - 1)
    ! on (eprev, etop) the branch has `nodes` nodes and D_l falls from +inf to -inf: bisect on its sign
    elo = eprev; ehi = etop
    do it = 1, 60
      emid = 0.5_dp * (elo + ehi)
      if (log_deriv(emid) > 0.0_dp) then; elo = emid; else; ehi = emid; end if
    end do
    ebot = 0.5_dp * (elo + ehi)
  contains
    integer function count_nodes(e)
      real(dp), intent(in) :: e
      real(dp) :: y(2, 0:0, m%n)
      integer :: i
      call integrate_regular(m, v, l, e, z, c, 0, m%n, y)
      count_nodes = 0
      do i = 3, m%n
        if (y(1, 0, i) * y(1, 0, i - 1) < 0.0_dp) count_nodes = count_nodes + 1
      end do
    end function count_nodes

    real(dp) function log_deriv(e)
      real(dp), intent(in) :: e
      real(dp) :: y(2, 0:0, m%n)
      call integrate_regular(m, v, l, e, z, c, 0, m%n, y)
      log_deriv = (1.0_dp + (e - v(m%n)) / (c * c)) * y(2, 0, m%n) * m%r(m%n) / y(1, 0, m%n)
    end function log_deriv

    !> lowest energy at which the branch has more than k nodes (between a, which has <= k, and b, which has > k)
    real(dp) function bisect_nodes(a, b, k)
      real(dp), intent(in) :: a, b
      integer, intent(in) :: k
      real(dp) :: lo, hi, mid
      integer :: j
      lo = a; hi = b
      do j = 1, 60
        mid = 0.5_dp * (lo + hi)
        if (count_nodes(mid) > k) then; hi = mid; else; lo = mid; end if
      end do
      bisect_nodes = 0.5_dp * (lo + hi)
    end function bisect_nodes
  end subroutine band_window

  !> E_nu kept inside [ebot, etop - 5% of the band window] of its own branch (see band_window); returns .true. if it was moved
  logical function clamp_enu(m, v, l, nodes, z, c, enu) result(moved)
    type(radmesh_t), intent(in) :: m
    real(dp), intent(in) :: v(:), z, c
    integer, intent(in) :: l, nodes
    real(dp), intent(inout) :: enu
    real(dp) :: eb, et, e0
    call band_window(m, v, l, nodes, z, c, eb, et)
    e0 = enu
    moved = .false.
    if (.not. (eb < et .and. abs(eb) < 1.0e3_dp .and. abs(et) < 1.0e3_dp)) return     ! window not found: leave E_nu
    enu = max(eb, min(et - 0.05_dp * (et - eb), enu))
    moved = (enu /= e0)
  end function clamp_enu
end module asr_potpar
