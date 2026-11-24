!> Real-space graph of a crystal structure (asr_crystal) for the recursion: a finite cluster of sites around the start site, with
!> exact site labels and one union list of bond slots.  No lattice names: everything follows from the structure, the bond radius
!> and the recursion depth.
!>
!> Why a cluster.  The recursion computes G_00(z) of the infinite crystal.  It needs a list of real sites with neighbour tables to
!> hop on, and it only ever visits the sites that the closure of asr_augcls reaches from the start site: two bonds per level.
!> So the graph holds every crystal site within the radius 2 levels bmax of the start site (bmax: the longest bond), or rmin if
!> that is larger (spectral targets, pair shells).  It is a ball about the start site, so every operation of the site group maps
!> it onto itself, and the displacements are the true vectors.  There is no periodic box: nothing wraps around, and nothing
!> depends on a box size.
!>
!> Sites.  A site is the exact label (cell n, basis site b), n relative to the cell of the start site.  Site 1 is the start site
!> (basis site b0 in cell 0).  The sites are numbered by distance from site 1, ties by z, then y, then x (structure_cluster).
!>
!> Slots.  The bond vectors of all basis sites up to the radius rcl form one list, sorted by length, then x, y, z (the order of
!> asr_lattices).  Slot k of site i is g%nbr(g%nbr_ptr(i) + k - 1): the site at r_i + bvec(k); 0 when the basis site of i has
!> no such bond (sites of different basis sites have different neighbour sets: hcp, NiAs); -1 when the neighbour lies outside
!> the cluster.  Every site within 2 levels - 1 bonds of site 1 has no entry -1 (asr_augcls stops if the closure meets one).
!> btgt(k, b) and bdn(:, k, b) give the neighbour exactly: basis site btgt, cell n + bdn.
!>
!> Site operations (graph_site_ops).  The space-group operations that fix the start site act on the sites as (n, b) ->
!> (m n + shift(:, b) - shift(:, b0), perm(b)).  They are ordered canonically: monomial operations first, by the axis permutation
!> (123, 132, 213, 231, 312, 321) and then by the signs of the rows (bit i - 1 set: row i negative), as
!> asr_symmetry::cubic_site_perms; the other operations follow in space-group order.  The identity is operation 1.
module asr_crystal_graph
  use asr_kinds, only: dp
  use asr_crystal, only: structure_t, geo_tol
  use asr_crystal_symmetry, only: symop_t, is_monomial, basis_orbits
  use asr_augspace, only: graph_t
  implicit none
  private
  public :: crystal_graph_t, crystal_graph, graph_site, graph_site_ops, structure_cluster

  !> Tolerance of distances (alat units): a site at distance r + rtol is still inside r (asr_lattices, asr_strconst).
  real(dp), parameter, public :: rtol = 1.0e-9_dp

  type :: crystal_graph_t
    type(graph_t) :: g                      !! neighbour lists, nslot entries per site (0: no such bond, -1: outside the cluster)
    integer :: n = 0, nslot = 0, nb = 0, b0 = 1, norbit = 0, levels = 0
    real(dp) :: radius = 0.0_dp             !! cluster radius (alat units)
    real(dp) :: bmax = 0.0_dp               !! longest bond (alat units)
    integer :: nbox(3) = 0                  !! cells -nbox .. nbox are indexed
    integer, allocatable :: cell(:, :)      !! (3, n) cell of each site, relative to the cell of site 1
    integer, allocatable :: basis(:)        !! (n) basis site of each site; the structure constants are tabulated per basis site
    integer, allocatable :: orbit(:)        !! (n) space-group orbit of the basis site (symmetry-distinct site), 1 .. norbit
    integer, allocatable :: index(:, :, :, :)   !! (-nbox(1):nbox(1), ..., nb) site of (cell, basis site); 0 outside the cluster
    real(dp), allocatable :: bvec(:, :)     !! (3, nslot) bond vectors, Cartesian, alat units
    integer, allocatable :: btgt(:, :)      !! (nslot, nb) basis site of the neighbour on slot k of basis site b (0: none)
    integer, allocatable :: bdn(:, :, :)    !! (3, nslot, nb) cell offset of that neighbour
    real(dp), allocatable :: disp(:, :)     !! (3, n) displacement of every site from site 1, Cartesian, alat units
  end type crystal_graph_t

  integer, parameter :: axperm(3, 6) = reshape([1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, 3, 2, 1], [3, 6])

contains

  !> Graph of structure s with all bonds up to rcl (alat units) on the cluster for a recursion of `levels` levels about basis site
  !> start (default 1).  ops: the space group of s (orbits).  rmin (optional, alat units): the cluster holds at least the ball of
  !> this radius.  msg is empty on success.
  subroutine crystal_graph(s, ops, rcl, levels, cg, msg, start, rmin)
    type(structure_t), intent(in) :: s
    type(symop_t), intent(in) :: ops(:)
    real(dp), intent(in) :: rcl
    integer, intent(in) :: levels
    type(crystal_graph_t), intent(out) :: cg
    character(len=:), allocatable, intent(out) :: msg
    integer, intent(in), optional :: start
    real(dp), intent(in), optional :: rmin
    real(dp), allocatable :: tv(:, :), key(:, :)
    integer, allocatable :: tb(:), tt(:), tdn(:, :), slot_of(:), ord(:), lab(:), orb(:)
    integer :: nmax(3), b, c, i, j, k, n1, n2, n3, nt, ns, cnt, t
    real(dp) :: f(3), v(3)
    msg = ''
    cg%b0 = 1; if (present(start)) cg%b0 = start
    if (cg%b0 < 1 .or. cg%b0 > s%nb) then; msg = 'crystal_graph: start is not a basis site'; return; end if
    if (levels < 0) then; msg = 'crystal_graph: levels must not be negative'; return; end if
    cg%nb = s%nb; cg%levels = levels
    ! ---- bonds of every basis site up to rcl
    do i = 1, 3; nmax(i) = ceiling(rcl * norm2(s%bvec(:, i))) + 1; end do
    allocate (tv(3, 64), tb(64), tt(64), tdn(3, 64)); nt = 0
    do b = 1, s%nb
      do c = 1, s%nb
        do n3 = -nmax(3), nmax(3); do n2 = -nmax(2), nmax(2); do n1 = -nmax(1), nmax(1)
          if (c == b .and. n1 == 0 .and. n2 == 0 .and. n3 == 0) cycle
          f = real([n1, n2, n3], dp) + s%frac(:, c) - s%frac(:, b)
          v = matmul(s%avec, f)
          if (norm2(v) > rcl + rtol) cycle
          if (nt == size(tb)) call grow()
          nt = nt + 1; tv(:, nt) = v; tb(nt) = b; tt(nt) = c; tdn(:, nt) = [n1, n2, n3]
        end do; end do; end do
      end do
    end do
    if (nt == 0) then; msg = 'crystal_graph: no bond within rcl'; return; end if
    ! ---- union of the bond vectors, sorted by length, then x, y, z
    allocate (slot_of(nt), key(4, nt)); ns = 0
    do i = 1, nt
      slot_of(i) = 0
      do j = 1, i - 1
        if (maxval(abs(tv(:, j) - tv(:, i))) < geo_tol) then; slot_of(i) = slot_of(j); exit; end if
      end do
      if (slot_of(i) == 0) then
        ns = ns + 1; slot_of(i) = ns; key(:, ns) = [sum(tv(:, i)**2), tv(:, i)]
      end if
    end do
    call sort_keys(key(:, 1:ns), ord)
    cg%nslot = ns
    allocate (cg%bvec(3, ns), cg%btgt(ns, s%nb), cg%bdn(3, ns, s%nb), lab(ns)); cg%btgt = 0; cg%bdn = 0
    do k = 1, ns; lab(ord(k)) = k; cg%bvec(:, k) = key(2:4, ord(k)); end do
    do i = 1, nt
      k = lab(slot_of(i))
      cg%btgt(k, tb(i)) = tt(i); cg%bdn(:, k, tb(i)) = tdn(:, i)
    end do
    cg%bmax = maxval(norm2(cg%bvec, 1))
    ! ---- sites: the ball about the start site
    cg%radius = 2.0_dp * real(levels, dp) * cg%bmax
    if (present(rmin)) cg%radius = max(cg%radius, rmin)
    call structure_cluster(s, cg%b0, cg%radius, cg%disp, cg%basis, cg%cell)
    cg%n = size(cg%basis)
    do i = 1, 3; cg%nbox(i) = maxval(abs(cg%cell(i, :))) + 1; end do
    allocate (cg%index(-cg%nbox(1):cg%nbox(1), -cg%nbox(2):cg%nbox(2), -cg%nbox(3):cg%nbox(3), s%nb)); cg%index = 0
    do i = 1, cg%n
      cg%index(cg%cell(1, i), cg%cell(2, i), cg%cell(3, i), cg%basis(i)) = i
    end do
    ! ---- neighbour lists
    cg%g%n = cg%n
    allocate (cg%g%nbr_ptr(cg%n + 1), cg%g%nbr(cg%n * ns))
    cnt = 0; cg%g%nbr_ptr(1) = 1
    do i = 1, cg%n
      b = cg%basis(i)
      do k = 1, ns
        cnt = cnt + 1; cg%g%nbr(cnt) = 0
        if (cg%btgt(k, b) == 0) cycle
        t = graph_site(cg, cg%cell(:, i) + cg%bdn(:, k, b), cg%btgt(k, b))
        cg%g%nbr(cnt) = merge(t, -1, t > 0)
      end do
      cg%g%nbr_ptr(i + 1) = cg%g%nbr_ptr(i) + ns
    end do
    ! ---- symmetry-distinct sites
    call basis_orbits(ops, s%nb, orb, cg%norbit)
    allocate (cg%orbit(cg%n))
    cg%orbit = orb(cg%basis)
  contains
    subroutine grow()
      real(dp), allocatable :: t2(:, :)
      integer, allocatable :: i1(:), i2(:), i3(:, :)
      allocate (t2(3, 2 * size(tb)), i1(2 * size(tb)), i2(2 * size(tb)), i3(3, 2 * size(tb)))
      t2(:, 1:nt) = tv(:, 1:nt); i1(1:nt) = tb(1:nt); i2(1:nt) = tt(1:nt); i3(:, 1:nt) = tdn(:, 1:nt)
      call move_alloc(t2, tv); call move_alloc(i1, tb); call move_alloc(i2, tt); call move_alloc(i3, tdn)
    end subroutine grow
  end subroutine crystal_graph

  !> Site of the crystal site (n, b) (n relative to the cell of site 1); 0 when it lies outside the cluster.
  pure integer function graph_site(cg, n, b)
    type(crystal_graph_t), intent(in) :: cg
    integer, intent(in) :: n(3), b
    graph_site = 0
    if (any(abs(n) > cg%nbox)) return
    graph_site = cg%index(n(1), n(2), n(3), b)
  end function graph_site

  !> The space-group operations (ops of s) that fix site 1, in canonical order (see the header).  perm(site, k): image of every
  !> site under operation k; rot(:, :, k): its Cartesian rotation; idx(k): its index in ops.  Operation 1 is the identity.  The
  !> cluster is a ball about site 1, so every such operation maps it onto itself.
  subroutine graph_site_ops(cg, s, ops, perm, rot, idx)
    type(crystal_graph_t), intent(in) :: cg
    type(structure_t), intent(in) :: s
    type(symop_t), intent(in) :: ops(:)
    integer, allocatable, intent(out) :: perm(:, :), idx(:)
    real(dp), allocatable, intent(out) :: rot(:, :, :)
    integer, allocatable :: cand(:), ord(:)
    real(dp), allocatable :: key(:, :)
    integer :: k, nk, i, b, b0, n2(3)
    if (s%nb /= cg%nb) error stop 'graph_site_ops: the graph belongs to another structure'
    b0 = cg%b0
    allocate (cand(size(ops))); nk = 0
    do k = 1, size(ops)
      if (ops(k)%perm(b0) /= b0) cycle
      nk = nk + 1; cand(nk) = k
    end do
    allocate (key(1, nk))
    do i = 1, nk
      key(1, i) = real(op_key(ops(cand(i))%rot, cand(i)), dp)
    end do
    call sort_keys(key, ord)
    allocate (perm(cg%n, nk), rot(3, 3, nk), idx(nk))
    do i = 1, nk
      k = cand(ord(i)); idx(i) = k; rot(:, :, i) = ops(k)%rot
      do b = 1, cg%n
        n2 = matmul(ops(k)%m, cg%cell(:, b)) + ops(k)%shift(:, cg%basis(b)) - ops(k)%shift(:, b0)
        perm(b, i) = graph_site(cg, n2, ops(k)%perm(cg%basis(b)))
      end do
    end do
    if (any(perm == 0)) error stop 'graph_site_ops: an operation maps a site out of the cluster'
    if (any(perm(:, 1) /= [(b, b=1, cg%n)])) error stop 'graph_site_ops: operation 1 is not the identity'
    if (any(perm(1, :) /= 1)) error stop 'graph_site_ops: an operation moves the start site'
  end subroutine graph_site_ops

  !> Cluster of all crystal sites within rmax (alat units) of basis site b: vec(:, i) relative to b, Cartesian.  The origin is
  !> first, then the sites by distance, ties by z, then y, then x (the order of asr_strconst::lattice_cluster).
  !> cb, cn (optional): basis site and cell (relative to the cell of b) of every cluster site.
  subroutine structure_cluster(s, b, rmax, vec, cb, cn)
    type(structure_t), intent(in) :: s
    integer, intent(in) :: b
    real(dp), intent(in) :: rmax
    real(dp), allocatable, intent(out) :: vec(:, :)
    integer, allocatable, intent(out), optional :: cb(:), cn(:, :)
    real(dp), allocatable :: tv(:, :), key(:, :)
    integer, allocatable :: tc(:), tn(:, :), ord(:)
    integer :: nmax(3), c, n1, n2, n3, n, i
    real(dp) :: v(3)
    do i = 1, 3; nmax(i) = ceiling(rmax * norm2(s%bvec(:, i))) + 1; end do
    n = s%nb * product(2 * nmax + 1)
    allocate (tv(3, n), tc(n), tn(3, n)); n = 0
    do c = 1, s%nb
      do n3 = -nmax(3), nmax(3); do n2 = -nmax(2), nmax(2); do n1 = -nmax(1), nmax(1)
        v = matmul(s%avec, real([n1, n2, n3], dp) + s%frac(:, c) - s%frac(:, b))
        if (sqrt(sum(v * v)) > rmax + rtol) cycle
        n = n + 1; tv(:, n) = v; tc(n) = c; tn(:, n) = [n1, n2, n3]
      end do; end do; end do
    end do
    allocate (key(4, n))
    do i = 1, n; key(:, i) = [sum(tv(:, i)**2), tv(3, i), tv(2, i), tv(1, i)]; end do
    call sort_keys(key, ord)
    allocate (vec(3, n)); vec = tv(:, ord)
    if (present(cb)) cb = tc(ord)
    if (present(cn)) cn = tn(:, ord)
    if (sum(vec(:, 1)**2) > 0.0_dp) error stop 'structure_cluster: the origin is not first'
  end subroutine structure_cluster

  ! ------------------------------------------------------------------ helpers
  !> Canonical key of a site operation: monomial operations by axis permutation and row signs, then the rest by index.
  pure integer function op_key(rot, k)
    real(dp), intent(in) :: rot(3, 3)
    integer, intent(in) :: k
    integer :: pp(3), i, ip, sgn
    if (.not. is_monomial(rot)) then; op_key = 1000 + k; return; end if
    sgn = 0
    do i = 1, 3
      pp(i) = maxloc(abs(rot(i, :)), 1)
      if (rot(i, pp(i)) < 0.0_dp) sgn = sgn + 2**(i - 1)
    end do
    do ip = 1, 6
      if (all(axperm(:, ip) == pp)) exit
    end do
    op_key = 8 * (ip - 1) + sgn
  end function op_key

  !> Stable sort of columns by their keys (first row first); entries within rtol are equal.
  subroutine sort_keys(key, ord)
    real(dp), intent(in) :: key(:, :)
    integer, allocatable, intent(out) :: ord(:)
    integer, allocatable :: tmp(:)
    integer :: n, w, lo, mid, hi, i, j, k
    n = size(key, 2)
    allocate (ord(n), tmp(n)); ord = [(i, i=1, n)]
    w = 1
    do while (w < n)                                 ! bottom-up merge sort
      lo = 1
      do while (lo <= n)
        mid = min(lo + w - 1, n); hi = min(lo + 2 * w - 1, n)
        i = lo; j = mid + 1; k = lo
        do while (i <= mid .and. j <= hi)
          if (less(ord(j), ord(i))) then; tmp(k) = ord(j); j = j + 1
          else; tmp(k) = ord(i); i = i + 1; end if
          k = k + 1
        end do
        do while (i <= mid); tmp(k) = ord(i); i = i + 1; k = k + 1; end do
        do while (j <= hi); tmp(k) = ord(j); j = j + 1; k = k + 1; end do
        lo = lo + 2 * w
      end do
      ord = tmp
      w = 2 * w
    end do
  contains
    logical function less(a, b)
      integer, intent(in) :: a, b
      integer :: r
      less = .false.
      do r = 1, size(key, 1)
        if (key(r, a) < key(r, b) - rtol) then; less = .true.; return; end if
        if (key(r, a) > key(r, b) + rtol) return
      end do
    end function less
  end subroutine sort_keys
end module asr_crystal_graph
