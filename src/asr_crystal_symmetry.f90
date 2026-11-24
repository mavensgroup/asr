!> Symmetry of a crystal structure (asr_crystal), found from the structure alone.  No lattice names.
!>
!> Space group.  A rotation is an integer matrix m in lattice coordinates (f -> m f + t) with det m = +-1 that keeps the metric,
!> m^T G m = G with G = avec^T avec; its Cartesian form is rot = avec m avec^-1.  Entries in -2..2 are searched, which holds the
!> point group of every cell that is not strongly skewed (the prototypes have entries in -1..1).  For each rotation the
!> translations follow from mapping one basis site onto every site with its label.  An operation is kept when every site maps onto
!> a site with the same label modulo a lattice vector.  Its action on the sites (n, b) of the infinite crystal is exact and integer:
!> (n, b) -> (m n + shift(:, b), perm(b)).
!>
!> Site group of basis site b: the operations with perm(b) = b, as rotations about the site.
!>
!> Orbital matrices.  d(k, j) with Y_k(rot p) = sum_j d(k, j) Y_j(p) for the nine real harmonics (s | y z x | xy yz 3z2-r2 xz
!> x2-y2), as asr_symmetry::orbital_rep but for any rotation (a threefold axis gives entries +-1/2, +-sqrt(3)/2).
!>
!> Orbital groups.  P = (1/|G|) sum_g D(g) (x) D(g) projects onto the matrices that commute with the site group.  The diagonal
!> elements of every such matrix (a site Green function block) satisfy M_LL = M_L'L' exactly when P vec(E_LL) = P vec(E_L'L'),
!> that is when (1/|G|) sum_g d(:, L, g) d(:, L, g)^T is the same matrix for L and L'.  Orbitals with the same matrix form a group.
!> One representative per group runs a chain; the group size is its weight.  Cubic sites give s | p (3) | t2g (3) | eg (2).
!>
!> Representatives and the reduction group.  The reduced augmented space needs monomial operations (signed permutations: every
!> real harmonic maps onto +- one harmonic) under which the start orbital is one dimensional.  The representative of a group is
!> its orbital with the largest little group {g monomial : |d(L, L, g)| = 1}; ties follow the order s, p_z, d_3z2, d_xy, d_x2-y2,
!> p_y, p_x, d_yz, d_xz, which gives today's cubic choice s, p_z, d_xy, d_3z2.  The reduction group is the intersection of the
!> little groups of all representatives: the D4h subset (16) on cubic sites, as today, 8 on simple hexagonal, 4 on hcp and NiAs.
module asr_crystal_symmetry
  use asr_kinds, only: dp
  use asr_crystal, only: structure_t, geo_tol, same_frac, inv3
  use asr_harmonics, only: ylm_real
  use asr_lapack, only: dgesv
  implicit none
  private
  public :: symop_t, space_group, site_ops, basis_orbits, is_monomial, orbital_matrix, orbital_groups, reduction_ops

  type :: symop_t
    integer :: m(3, 3) = 0                !! rotation in lattice coordinates
    real(dp) :: rot(3, 3) = 0.0_dp        !! the rotation in Cartesian coordinates (orthogonal)
    real(dp) :: t(3) = 0.0_dp             !! translation, fractional, in [0, 1)
    integer, allocatable :: perm(:)       !! basis site b goes to perm(b)
    integer, allocatable :: shift(:, :)   !! (3, nb): m frac_b + t = frac_perm(b) + shift(:, b)
  end type symop_t

  !> Tie-break order of the representatives (orbital indices).
  integer, parameter :: pref(9) = [1, 3, 7, 5, 9, 2, 4, 6, 8]
  real(dp), parameter :: otol = 1.0e-9_dp

contains

  !> All operations {m | t} that map the structure onto itself, labels included.  The identity is operation 1.
  subroutine space_group(s, ops)
    type(structure_t), intent(in) :: s
    type(symop_t), allocatable, intent(out) :: ops(:)
    type(symop_t), allocatable :: tmp(:)
    type(symop_t) :: swp
    real(dp) :: g(3, 3), ainv(3, 3), rm(3, 3), t(3), fb(3)
    integer :: m(3, 3), code, k, i0, j, b, c, n, nlab, lab, perm(s%nb), shift(3, s%nb)
    logical :: ok
    g = matmul(transpose(s%avec), s%avec)
    ainv = inv3(s%avec)
    ! anchor: the first site of the label with the fewest sites
    i0 = 1; nlab = s%nb + 1
    do lab = 1, s%nlabel
      n = count(s%label == lab)
      if (n > 0 .and. n < nlab) then; nlab = n; i0 = findloc(s%label, lab, 1); end if
    end do
    allocate (tmp(48 * s%nb)); n = 0
    do code = 0, 5**9 - 1
      k = code
      do j = 1, 9
        m(mod(j - 1, 3) + 1, (j - 1) / 3 + 1) = mod(k, 5) - 2; k = k / 5
      end do
      if (abs(idet3(m)) /= 1) cycle
      rm = real(m, dp)
      if (maxval(abs(matmul(transpose(rm), matmul(g, rm)) - g)) > 1.0e-8_dp) cycle
      do j = 1, s%nb
        if (s%label(j) /= s%label(i0)) cycle
        t = s%frac(:, j) - matmul(rm, s%frac(:, i0)); t = t - floor(t + geo_tol)
        ok = .true.
        do b = 1, s%nb
          fb = matmul(rm, s%frac(:, b)) + t
          perm(b) = 0
          do c = 1, s%nb
            if (s%label(c) == s%label(b) .and. same_frac(fb, s%frac(:, c))) then; perm(b) = c; exit; end if
          end do
          if (perm(b) == 0) then; ok = .false.; exit; end if
          shift(:, b) = nint(fb - s%frac(:, perm(b)))
        end do
        if (.not. ok) cycle
        if (known(m, t)) cycle
        n = n + 1
        tmp(n)%m = m; tmp(n)%t = t; tmp(n)%rot = matmul(s%avec, matmul(rm, ainv))
        tmp(n)%perm = perm; tmp(n)%shift = shift
      end do
    end do
    ! identity first
    do k = 1, n
      if (all(tmp(k)%m == reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])) .and. all(abs(tmp(k)%t) < geo_tol)) exit
    end do
    if (k > n) error stop 'space_group: the identity was not found'
    if (k /= 1) then; swp = tmp(1); tmp(1) = tmp(k); tmp(k) = swp; end if
    allocate (ops(n)); ops = tmp(1:n)
  contains
    logical function known(mm, tt)
      integer, intent(in) :: mm(3, 3)
      real(dp), intent(in) :: tt(3)
      integer :: q
      known = .false.
      do q = 1, n
        if (all(tmp(q)%m == mm) .and. same_frac(tmp(q)%t, tt)) then; known = .true.; return; end if
      end do
    end function known
  end subroutine space_group

  !> Indices (into ops) of the site group of basis site b: the operations that map b onto itself modulo a lattice vector.
  subroutine site_ops(ops, b, idx)
    type(symop_t), intent(in) :: ops(:)
    integer, intent(in) :: b
    integer, allocatable, intent(out) :: idx(:)
    integer :: k
    idx = pack([(k, k=1, size(ops))], [(ops(k)%perm(b) == b, k=1, size(ops))])
  end subroutine site_ops

  !> Space-group orbits of the nb basis sites (symmetry-distinct sites): orbit(b) in 1 .. norbit, numbered by their first site.
  subroutine basis_orbits(ops, nb, orbit, norbit)
    type(symop_t), intent(in) :: ops(:)
    integer, intent(in) :: nb
    integer, allocatable, intent(out) :: orbit(:)
    integer, intent(out) :: norbit
    integer :: b, k
    allocate (orbit(nb)); orbit = 0; norbit = 0
    do b = 1, nb
      if (orbit(b) > 0) cycle
      norbit = norbit + 1
      do k = 1, size(ops)
        orbit(ops(k)%perm(b)) = norbit
      end do
    end do
  end subroutine basis_orbits

  !> Is a Cartesian rotation a signed permutation (one entry +-1 per row and column, zeros elsewhere)?
  pure logical function is_monomial(rot)
    real(dp), intent(in) :: rot(3, 3)
    integer :: i
    is_monomial = all(abs(rot) < geo_tol .or. abs(abs(rot) - 1.0_dp) < geo_tol)
    if (.not. is_monomial) return
    do i = 1, 3
      if (count(abs(rot(i, :)) > 0.5_dp) /= 1 .or. count(abs(rot(:, i)) > 0.5_dp) /= 1) is_monomial = .false.
    end do
  end function is_monomial

  !> Orbital matrix of a Cartesian rotation on the nine real harmonics: Y_k(rot p) = sum_j d(k, j) Y_j(p).  Least squares on fixed
  !> points (a Weyl sequence), so the result is the same in every call; entries within 1e-9 of an integer are set to it.
  subroutine orbital_matrix(rot, d)
    real(dp), intent(in) :: rot(3, 3)
    real(dp), intent(out) :: d(9, 9)
    integer, parameter :: np = 24
    real(dp), parameter :: wq(3) = [0.414213562373095_dp, 0.732050807568877_dp, 0.236067977499790_dp]
    real(dp) :: pts(3, np), a(np, 9), bm(np, 9), yv(9), mp(3), ata(9, 9), x(9, 9)
    integer :: ip, k, ipiv(9), info
    do ip = 1, np
      pts(:, ip) = real(ip, dp) * wq; pts(:, ip) = pts(:, ip) - aint(pts(:, ip)) - 0.5_dp
      call ylm_real(2, pts(1, ip), pts(2, ip), pts(3, ip), yv); a(ip, :) = yv
      mp = matmul(rot, pts(:, ip))
      call ylm_real(2, mp(1), mp(2), mp(3), yv); bm(ip, :) = yv
    end do
    ata = matmul(transpose(a), a)
    x = matmul(transpose(a), bm)
    call dgesv(9, 9, ata, 9, ipiv, x, 9, info)
    if (info /= 0) error stop 'orbital_matrix: singular least-squares system'
    d = transpose(x)
    do k = 1, 9
      where (abs(d(:, k) - anint(d(:, k))) < otol) d(:, k) = anint(d(:, k))
    end do
  end subroutine orbital_matrix

  !> Orbital groups of a site group.  d(9, 9, nop): orbital matrices of all site operations (monomial or not); mono(nop): which
  !> of them are monomial.  Out: grp(9) = group of each orbital (1..ngrp), rep(1:ngrp) = representative orbital,
  !> weight(1:ngrp) = group size.  Groups are numbered by their first orbital.
  subroutine orbital_groups(d, mono, grp, ngrp, rep, weight)
    real(dp), intent(in) :: d(:, :, :)
    logical, intent(in) :: mono(:)
    integer, intent(out) :: grp(9), ngrp, rep(9), weight(9)
    real(dp) :: u(9, 9, 9)
    integer :: l, l2, g, nlit, best, k
    u = 0.0_dp
    do l = 1, 9
      do g = 1, size(d, 3)
        u(:, :, l) = u(:, :, l) + spread(d(:, l, g), 2, 9) * spread(d(:, l, g), 1, 9)
      end do
    end do
    u = u / real(size(d, 3), dp)
    grp = 0; ngrp = 0
    do l = 1, 9
      do l2 = 1, l - 1
        if (maxval(abs(u(:, :, l) - u(:, :, l2))) < 1.0e-8_dp) then; grp(l) = grp(l2); exit; end if
      end do
      if (grp(l) == 0) then; ngrp = ngrp + 1; grp(l) = ngrp; end if
    end do
    rep = 0; weight = 0
    do k = 1, ngrp
      weight(k) = count(grp == k)
      best = -1
      do l = 1, 9                                   ! in preference order
        if (grp(pref(l)) /= k) cycle
        nlit = count([(mono(g) .and. abs(abs(d(pref(l), pref(l), g)) - 1.0_dp) < otol, g=1, size(d, 3))])
        if (nlit > best) then; best = nlit; rep(k) = pref(l); end if
      end do
    end do
  end subroutine orbital_groups

  !> Indices (1..nop) of the reduction group: the monomial operations under which every orbital in reps is one dimensional.
  subroutine reduction_ops(d, mono, reps, idx)
    real(dp), intent(in) :: d(:, :, :)
    logical, intent(in) :: mono(:)
    integer, intent(in) :: reps(:)
    integer, allocatable, intent(out) :: idx(:)
    integer :: g, k
    logical :: keep(size(d, 3))
    do g = 1, size(d, 3)
      keep(g) = mono(g)
      do k = 1, size(reps)
        if (abs(abs(d(reps(k), reps(k), g)) - 1.0_dp) >= otol) keep(g) = .false.
      end do
    end do
    idx = pack([(g, g=1, size(d, 3))], keep)
  end subroutine reduction_ops

  pure integer function idet3(m)
    integer, intent(in) :: m(3, 3)
    idet3 = m(1, 1) * (m(2, 2) * m(3, 3) - m(2, 3) * m(3, 2)) - m(1, 2) * (m(2, 1) * m(3, 3) - m(2, 3) * m(3, 1)) &
            + m(1, 3) * (m(2, 1) * m(3, 2) - m(2, 2) * m(3, 1))
  end function idet3
end module asr_crystal_symmetry
