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
  public :: potpar_l_t, potpar_at_enu, band_window, clamp_enu, enu_window_margin

  !> clamp_enu keeps E_nu below etop - enu_window_margin * (etop - ebot).  A numerical safeguard of the SCF path: a result is
  !> valid only if the clamp is not active in the last iteration (the SCF drivers count it; then the result does not depend on
  !> this number).  Changed by tests only.
  real(dp) :: enu_window_margin = 0.05_dp

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

  !> Window of the l channel on the radial branch with `nodes` nodes inside the sphere:
  !>   ebot: D_l(E) = l (free-electron log derivative: Questaal's floor on the fractional part of pnu; below the Wigner-Seitz
  !>   band bottom D_l = 0);  etop: phi_l(S) = 0 (band top, the next node enters at S).
  !> Past etop the parametrisation jumps to the next branch (sqrt(Delta) changes sign, C runs away); Questaal avoids this by
  !> floating pnu with a fixed principal quantum number.
  !>
  !> Root search (one outward integration per energy, nshoot of them on output if present).  A node enters the sphere
  !> exactly where phi_l(S; E) changes sign, and nodes never leave it as E rises (Sturm).  So once the node count is k at one
  !> end of a bracket and k + 1 at the other, phi_l(S; E) has exactly one zero inside: it is found by root_t on phi_l(S).  On
  !> (eprev, etop) phi_l(S) has no zero and D_l falls monotonically, so g = M S Q - l P = phi(S) (D_l - l) is continuous with
  !> one sign change there: ebot is its root.  Both to 1e-13 relative (the bisections they replace ran 60 steps, to rounding).
  !> If a bracket does not show the expected sign change, the plain bisection of before runs instead.
  subroutine band_window(m, v, l, nodes, z, c, ebot, etop, nshoot)
    type(radmesh_t), intent(in) :: m
    real(dp), intent(in) :: v(:), z, c
    integer, intent(in) :: l, nodes
    real(dp), intent(out) :: ebot, etop
    integer, intent(out), optional :: nshoot
    real(dp) :: elo, ehi, emid, eprev, pv, gv, glo, ghi
    real(dp) :: y(2, 0:0, m%n)
    integer :: it, nc, nev
    type(root_t) :: rt
    nev = 0
    elo = v(m%n) - 20.0_dp                      ! below the valence branches (semicore tops are above ~ -10 Ry)
    ehi = v(m%n) + 1.0_dp
    call shoot(ehi, nc, pv, gv)
    do while (nc <= nodes .and. ehi < 400.0_dp)
      ehi = v(m%n) + 2.0_dp * (ehi - v(m%n))
      call shoot(ehi, nc, pv, gv)
    end do
    etop = transition(elo, ehi, nodes)
    eprev = elo
    if (nodes > 0) eprev = transition(elo, etop, nodes - 1)
    ! on (eprev, etop) the branch has `nodes` nodes and D_l falls from +inf to -inf: the root of g, D_l = l
    call shoot(eprev, nc, pv, glo); call shoot(etop, nc, pv, ghi)
    if (glo * ghi < 0.0_dp) then
      call root_init(rt, eprev, glo, etop, ghi, 0.0_dp, 1.0e-13_dp)
      do it = 1, 300
        if (rt%done) exit
        emid = root_next(rt)
        call shoot(emid, nc, pv, gv)
        call root_update(rt, emid, gv)
      end do
      ebot = rt%x
    else                                        ! as before: bisection on the sign of D_l - l
      elo = eprev; ehi = etop
      do it = 1, 60
        emid = 0.5_dp * (elo + ehi)
        call shoot(emid, nc, pv, gv)
        if (gv / pv > 0.0_dp) then; elo = emid; else; ehi = emid; end if
      end do
      ebot = 0.5_dp * (elo + ehi)
    end if
    if (present(nshoot)) nshoot = nev
  contains
    !> outward integration at e: node count in the sphere, P(S), and g = M S Q(S) - l P(S) (D_l - l = g / P), both divided by
    !> the Pruefer radius sqrt(P^2 + (M S Q)^2) at S
    subroutine shoot(e, nn, ps, gs)
      real(dp), intent(in) :: e
      integer, intent(out) :: nn
      real(dp), intent(out) :: ps, gs
      real(dp) :: qs, rs
      integer :: i
      nev = nev + 1
      call integrate_regular(m, v, l, e, z, c, 0, m%n, y)
      nn = 0
      do i = 3, m%n
        if (y(1, 0, i) * y(1, 0, i - 1) < 0.0_dp) nn = nn + 1
      end do
      ps = y(1, 0, m%n)
      qs = (1.0_dp + (e - v(m%n)) / (c * c)) * y(2, 0, m%n) * m%r(m%n)
      gs = qs - real(l, dp) * ps
      ! scaled by |(P, M S Q)| (the Pruefer radius at S, never zero): the same signs and zeros, but bounded and smooth in E,
      ! where P(S) itself grows exponentially with E
      rs = hypot(ps, qs)
      if (rs > 0.0_dp) then; ps = ps / rs; gs = gs / rs; end if
    end subroutine shoot

    !> lowest energy at which the branch has more than k nodes (between a, which has <= k, and b, which has > k).  If a has more
    !> than k nodes already (a semicore branch whose lower neighbour lies below a), a: what the bisection converged to.
    real(dp) function transition(a, b, k)
      real(dp), intent(in) :: a, b
      integer, intent(in) :: k
      real(dp) :: lo, hi, mid, plo, phi, pm, gm
      integer :: j, nlo, nhi, nm
      lo = a; hi = b
      call shoot(lo, nlo, plo, gm)
      if (nlo > k) then                         ! already more than k nodes at a: the bisection would converge onto a
        transition = a; return
      end if
      call shoot(hi, nhi, phi, gm)
      do j = 1, 60
        if (nlo == k .and. nhi == k + 1 .and. plo * phi < 0.0_dp) then   ! one node enters in (lo, hi): the zero of P(S)
          call root_init(rt, lo, plo, hi, phi, 0.0_dp, 1.0e-13_dp)
          do nm = 1, 300
            if (rt%done) exit
            mid = root_next(rt)
            call shoot(mid, nc, pm, gm)
            call root_update(rt, mid, pm)
          end do
          transition = rt%x
          return
        end if
        mid = 0.5_dp * (lo + hi)
        call shoot(mid, nm, pm, gm)
        if (nm > k) then; hi = mid; nhi = nm; phi = pm; else; lo = mid; nlo = nm; plo = pm; end if
      end do
      transition = 0.5_dp * (lo + hi)
    end function transition
  end subroutine band_window

  !> E_nu kept inside [ebot, etop - enu_window_margin (5 %) of the band window] of its own branch (see band_window); returns
  !> .true. if it was moved
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
    enu = max(eb, min(et - enu_window_margin * (et - eb), enu))
    moved = (enu /= e0)
  end function clamp_enu
end module asr_potpar
