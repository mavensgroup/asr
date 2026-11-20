!> The root searches of bound_state (core states) and band_window (E_nu windows) against the plain bisections they replaced
!> (copied below from main 5d42a7f as old_bound_state / old_band_window): the same eigenvalues, wave functions and windows to
!> 1e-11 relative, with fewer radial integrations.  Free Fe and Al atoms (scalar relativistic, and Fe non-relativistic): the
!> core shells in the Wigner-Seitz sphere (as asr_class_scf: elo = -1.6 Z^2 - 20, ehi = v(S) - 1e-6), every shell of the free
!> atom on its own mesh (as asr_atom: ehi = -1e-4), and the windows of the valence and semicore branches.
!> The algorithms of main 5d42a7f (bisection), with a count of the radial integrations.
module test_root_search_old
  use asr_kinds, only: dp
  use asr_radial
  implicit none
  private
  public :: old_bound_state, old_band_window
contains
  subroutine old_bound_state(m, v, l, nod, z, c, isph, elo, ehi, e, y, nshoot)
    type(radmesh_t), intent(in) :: m
    real(dp), intent(in) :: v(:), z, c, elo, ehi
    integer, intent(in) :: l, nod, isph
    real(dp), intent(out) :: e
    real(dp), intent(out) :: y(2, 0:0, isph)
    integer, intent(out) :: nshoot
    logical :: found
    real(dp) :: e1, e2, ea, eb, fa, fb, fm, em, vout, s, ehi_eff, ll
    integer :: it, ntot, nc
    real(dp), allocatable :: yo(:, :, :), yi(:, :)
    allocate (yo(2, 0:0, isph), yi(2, 2:isph))
    found = .true.; nshoot = 0
    s = m%r(isph); vout = v(isph); ll = real(l * (l + 1), dp)
    ehi_eff = min(ehi, vout - 1.0e-8_dp)
    e1 = bisect_nodes(nod)
    e2 = bisect_nodes(nod + 1)
    ea = e1; eb = e2
    call eval(ea, fa, ntot, nc); call eval(eb, fb, ntot, nc)
    if (fa * fb > 0.0_dp) call scan_bracket(ea, eb, fa, fb)
    do it = 1, 200
      em = 0.5_dp * (ea + eb)
      call eval(em, fm, ntot, nc)
      if (fa * fm <= 0.0_dp) then; eb = em; fb = fm; else; ea = em; fa = fm; end if
      if (abs(eb - ea) < 1.0e-13_dp * max(1.0_dp, abs(em))) exit
    end do
    e = 0.5_dp * (ea + eb)
    ! final wavefunction: outward to nc, inward beyond, matched in P
    call eval(e, fm, ntot, nc)
    y(:, 0, 1:nc) = yo(:, 0, 1:nc)
    block
      real(dp) :: sc, nrm
      integer :: j
      sc = yo(1, 0, nc) / yi(1, nc)
      if (abs(yi(1, nc)) < 1.0e-300_dp) sc = 1.0_dp
      do j = nc + 1, isph
        y(1, 0, j) = sc * yi(1, j); y(2, 0, j) = sc * yi(2, j)
      end do
      call norm_pq(m, c, 0, isph, y, nrm)
    end block
  contains
    subroutine eval(en, f, nt, nmatch)
      real(dp), intent(in) :: en
      real(dp), intent(out) :: f
      integer, intent(out) :: nt, nmatch
      real(dp) :: kap, qp, fo, fi, cumk, kloc
      integer :: j, it2, istart
      logical :: usewkb
      nshoot = nshoot + 1
      it2 = 0
      do j = 3, isph
        if (en - v(j) - ll / m%r(j)**2 > 0.0_dp) it2 = j
      end do
      nmatch = min(it2 + 10, isph - 1)
      if (it2 == 0) nmatch = min(isph - 1, 5)
      call integrate_regular(m, v, l, en, z, c, 0, nmatch, yo(:, :, 1:nmatch))
      ! where to start the inward integration: Hankel boundary condition at the sphere edge if the
      ! decay is weak; otherwise a local WKB decaying start where the accumulated decay is ~ e^-35
      cumk = 0.0_dp; istart = isph; usewkb = .false.
      do j = nmatch + 1, isph
        kloc = sqrt(max(v(j) + ll / m%r(j)**2 - en, 0.0_dp))
        cumk = cumk + kloc * (m%r(j) - m%r(j - 1))
        if (cumk > 35.0_dp) then; istart = j; usewkb = .true.; exit; end if
      end do
      if (usewkb) then
        kap = sqrt(max(v(istart) + ll / m%r(istart)**2 - en, 1.0e-12_dp))
        qp = -(kap + 1.0_dp / m%r(istart)) / (1.0_dp + (en - v(istart)) / c**2)
      else
        kap = sqrt(max(vout - en, 1.0e-12_dp))
        qp = hankel_logder(l, kap, s, en, vout, c)
      end if
      yi = 0.0_dp
      call integrate_inward(m, v, l, en, c, istart, nmatch, qp, yi(:, nmatch:istart))
      nt = 0
      do j = 3, nmatch
        if (yo(1, 0, j) * yo(1, 0, j - 1) < 0.0_dp) nt = nt + 1
      end do
      do j = nmatch + 1, isph
        if (yi(1, j) * yi(1, j - 1) < 0.0_dp) nt = nt + 1
      end do
      ! scale-safe: normalise each (P,Q) pair by its own magnitude before forming the cross product
      fo = max(abs(yo(1, 0, nmatch)), abs(yo(2, 0, nmatch)))
      fi = max(abs(yi(1, nmatch)), abs(yi(2, nmatch)))
      f = (yo(2, 0, nmatch) / fo) * (yi(1, nmatch) / fi) - (yi(2, nmatch) / fi) * (yo(1, 0, nmatch) / fo)
    end subroutine eval
    function bisect_nodes(target) result(eb_)
      integer, intent(in) :: target
      real(dp) :: eb_, lo, hi, mid, f
      integer :: k, nt, nmatch
      lo = elo; hi = ehi_eff
      do k = 1, 200
        mid = 0.5_dp * (lo + hi)
        call eval(mid, f, nt, nmatch)
        if (nt >= target) then; hi = mid; else; lo = mid; end if
        if (hi - lo < 1.0e-12_dp) exit
      end do
      eb_ = hi
    end function bisect_nodes
    subroutine scan_bracket(ea_, eb_, fa_, fb_)
      real(dp), intent(inout) :: ea_, eb_, fa_, fb_
      real(dp) :: x0, x1, f0, f1
      integer :: k, nt, nmatch
      x0 = ea_; call eval(x0, f0, nt, nmatch)
      do k = 1, 400
        x1 = ea_ + (eb_ - ea_) * k / 400.0_dp
        call eval(x1, f1, nt, nmatch)
        if (f0 * f1 <= 0.0_dp) then; ea_ = x0; eb_ = x1; fa_ = f0; fb_ = f1; return; end if
        x0 = x1; f0 = f1
      end do
      found = .false.
      ea_ = eb_; fa_ = fb_                       ! best estimate: upper end of the node window
    end subroutine scan_bracket
  end subroutine old_bound_state

  subroutine old_band_window(m, v, l, nodes, z, c, ebot, etop, nshoot)
    type(radmesh_t), intent(in) :: m
    real(dp), intent(in) :: v(:), z, c
    integer, intent(in) :: l, nodes
    real(dp), intent(out) :: ebot, etop
    integer, intent(out) :: nshoot
    real(dp) :: elo, ehi, emid, eprev
    integer :: it
    nshoot = 0
    elo = v(m%n) - 20.0_dp                      ! below the valence branches (semicore tops are above ~ -10 Ry)
    ehi = v(m%n) + 1.0_dp
    do while (count_nodes(ehi) <= nodes .and. ehi < 400.0_dp)
      ehi = v(m%n) + 2.0_dp * (ehi - v(m%n))
    end do
    etop = bisect_nodes(elo, ehi, nodes)
    eprev = elo
    if (nodes > 0) eprev = bisect_nodes(elo, etop, nodes - 1)
    ! on (eprev, etop) the branch has `nodes` nodes and D_l falls from +inf to -inf: bisect on the sign of D_l - l
    elo = eprev; ehi = etop
    do it = 1, 60
      emid = 0.5_dp * (elo + ehi)
      if (log_deriv(emid) > real(l, dp)) then; elo = emid; else; ehi = emid; end if
    end do
    ebot = 0.5_dp * (elo + ehi)
  contains
    integer function count_nodes(e)
      real(dp), intent(in) :: e
      real(dp) :: y(2, 0:0, m%n)
      integer :: i
      nshoot = nshoot + 1
      call integrate_regular(m, v, l, e, z, c, 0, m%n, y)
      count_nodes = 0
      do i = 3, m%n
        if (y(1, 0, i) * y(1, 0, i - 1) < 0.0_dp) count_nodes = count_nodes + 1
      end do
    end function count_nodes

    real(dp) function log_deriv(e)
      real(dp), intent(in) :: e
      real(dp) :: y(2, 0:0, m%n)
      nshoot = nshoot + 1
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
  end subroutine old_band_window
end module test_root_search_old

program test_root_search
  use asr_kinds, only: dp
  use asr_testing
  use asr_radial
  use asr_atom
  use asr_potpar
  use test_root_search_old
  implicit none
  integer :: nnew, nold
  print '(a)', 'test_root_search'
  nnew = 0; nold = 0
  call one_atom(26.0_dp, fe_shells(), 2.6671_dp, c_ry, 'Fe rel')
  call one_atom(26.0_dp, fe_shells(), 2.6671_dp, c_nonrel, 'Fe nonrel')
  call one_atom(13.0_dp, al_shells(), 2.99_dp, c_ry, 'Al rel')
  print '(a,i0,a,i0)', '    radial integrations: new ', nnew, ', old ', nold
  call check(real(nnew, dp) < 0.3_dp * real(nold, dp), 'under 30 % of the radial integrations of the bisections')
  call report()
contains
  subroutine one_atom(z, sh, s, c, tag)
    real(dp), intent(in) :: z, s, c
    type(shell_t), intent(in) :: sh(:)
    character(len=*), intent(in) :: tag
    type(atom_result_t) :: at
    type(radmesh_t) :: m
    real(dp), allocatable :: v(:), yn(:, :, :), yo(:, :, :)
    real(dp) :: en, eo, bn, tn, bo, to
    integer :: k, is, nod, ns1, ns2, l, nd, iw
    integer, parameter :: lw(5) = [0, 0, 1, 1, 2], nw(5) = [3, 2, 2, 1, 0]
    logical :: ok
    call solve_atom(z, sh, 'vbh', 30.0_dp, 0.01_dp, c < 1.0e5_dp, at)
    ! every shell of the free atom on its own mesh
    allocate (yn(2, 0:0, at%mesh%n), yo(2, 0:0, at%mesh%n))
    do is = 1, size(sh)
      nod = sh(is)%n - sh(is)%l - 1
      call bound_state(at%mesh, at%v, sh(is)%l, nod, z, c, at%mesh%n, -1.6_dp * z * z - 20.0_dp, -1.0e-4_dp, en, yn, ok, ns1)
      call old_bound_state(at%mesh, at%v, sh(is)%l, nod, z, c, at%mesh%n, -1.6_dp * z * z - 20.0_dp, -1.0e-4_dp, eo, yo, ns2)
      call compare(tag // ' free atom', sh(is)%n, sh(is)%l, en, eo, yn, yo, ns1, ns2)
    end do
    deallocate (yn, yo)
    ! core shells in the sphere
    call make_radmesh(z, s, 0.012_dp, m, c)
    allocate (v(m%n))
    do k = 2, m%n
      v(k) = lin(at%mesh%r, at%v * at%mesh%r, m%r(k)) / m%r(k)
    end do
    v(1) = v(2)
    allocate (yn(2, 0:0, m%n), yo(2, 0:0, m%n))
    do is = 1, size(sh)
      if (sh(is)%n >= 4 .or. sh(is)%l == 2) cycle                       ! the valence shells are not core states
      if (z < 20.0_dp .and. sh(is)%n >= 3) cycle
      nod = sh(is)%n - sh(is)%l - 1
      call bound_state(m, v, sh(is)%l, nod, z, c, m%n, -1.6_dp * z * z - 20.0_dp, v(m%n) - 1.0e-6_dp, en, yn, ok, ns1)
      call old_bound_state(m, v, sh(is)%l, nod, z, c, m%n, -1.6_dp * z * z - 20.0_dp, v(m%n) - 1.0e-6_dp, eo, yo, ns2)
      call compare(tag // ' sphere core', sh(is)%n, sh(is)%l, en, eo, yn, yo, ns1, ns2)
    end do
    ! windows: the valence branches and, below them, the semicore branches
    do iw = 1, 5                                ! (l, nodes): Fe 4s 3s 4p 3p 3d, Al 3s 3p 3d (and 2s 2p below)
      l = lw(iw); nd = nw(iw); if (z < 20.0_dp) nd = nd - 1
      if (nd < 0) cycle
      block
        call band_window(m, v, l, nd, z, c, bn, tn, ns1)
        call old_band_window(m, v, l, nd, z, c, bo, to, ns2)
        nnew = nnew + ns1; nold = nold + ns2
        print '(4x,a,a,i2,a,i2,a,2f13.7,a,2es9.1,a,i4,a,i4)', tag, ' window l', l, ' nodes', nd, ':', bn, tn, '  diff', &
          bn - bo, tn - to, '  shots', ns1, ' /', ns2
        call check(bn < tn, tag // ' window ordered')
        call check(abs(bn - bo) < 1.0e-11_dp * max(1.0_dp, abs(bo)), tag // ' ebot as the bisection')
        call check(abs(tn - to) < 1.0e-11_dp * max(1.0_dp, abs(to)), tag // ' etop as the bisection')
      end block
    end do
  end subroutine one_atom

  subroutine compare(tag, n, l, en, eo, yn, yo, ns1, ns2)
    character(len=*), intent(in) :: tag
    integer, intent(in) :: n, l, ns1, ns2
    real(dp), intent(in) :: en, eo, yn(:, 0:, :), yo(:, 0:, :)
    real(dp) :: dy
    nnew = nnew + ns1; nold = nold + ns2
    dy = maxval(abs(yn(1, 0, :) - yo(1, 0, :))) / maxval(abs(yo(1, 0, :)))
    print '(4x,a,a,i1,a,i1,a,f16.8,a,es9.1,a,es9.1,a,i4,a,i4)', tag, ' n', n, ' l', l, ': E', en, '  dE', en - eo, '  dP', dy, &
      '  shots', ns1, ' /', ns2
    call check(abs(en - eo) < 1.0e-11_dp * max(1.0_dp, abs(eo)), tag // ' eigenvalue as the bisection')
    call check(dy < 1.0e-8_dp, tag // ' wave function as the bisection')
  end subroutine compare

  real(dp) function lin(r, f, r0)
    real(dp), intent(in) :: r(:), f(:), r0
    integer :: j
    j = max(1, min(size(r) - 1, count(r <= r0)))
    lin = f(j) + (f(j + 1) - f(j)) * (r0 - r(j)) / (r(j + 1) - r(j))
  end function lin

  function fe_shells() result(sh)
    type(shell_t), allocatable :: sh(:)
    sh = [shell_t(1, 0, 2.0_dp), shell_t(2, 0, 2.0_dp), shell_t(2, 1, 6.0_dp), shell_t(3, 0, 2.0_dp), shell_t(3, 1, 6.0_dp), &
          shell_t(3, 2, 6.0_dp), shell_t(4, 0, 2.0_dp)]
  end function fe_shells

  function al_shells() result(sh)
    type(shell_t), allocatable :: sh(:)
    sh = [shell_t(1, 0, 2.0_dp), shell_t(2, 0, 2.0_dp), shell_t(2, 1, 6.0_dp), shell_t(3, 0, 2.0_dp), shell_t(3, 1, 1.0_dp)]
  end function al_shells

end program test_root_search
