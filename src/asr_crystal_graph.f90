!> Periodic graph of a crystal structure (asr_crystal) on a supercell, with exact site labels and one union list of bond slots.
!> No lattice names: everything follows from the structure, the supercell and the cluster radius.
!>
!> Supercell.  Its vectors are the columns of an integer matrix msc in primitive lattice coordinates (avec msc in Cartesian), so
!> any supercell can be given: L1 x L2 x L3 primitive cells (msc diagonal) or L conventional cubes of a centred cubic lattice
!> (supercell_from_cart).  A cell n (integer, primitive coordinates) is reduced modulo the supercell lattice by its upper
!> triangular basis h (column operations on msc): 0 <= n_i < h(i, i).  So a site of the graph is the exact label (cell, basis).
!>
!> Sites.  Site 1 is basis site b0 in cell 0 (the start site).  The sites are numbered by their position relative to site 1 in
!> supercell coordinates wrapped into [0, 1), third coordinate slowest.  For a cubic supercell of a cubic lattice this is the
!> half-lattice order of asr_lattices (z, then y, then x), so the graph is the same graph, site by site.
!>
!> Slots.  The bond vectors of all basis sites up to the radius rcl form one list, sorted by length, then x, y, z (the order of
!> asr_lattices).  Slot k of site i is the neighbour at r_i + bvec(k): g%nbr(g%nbr_ptr(i) + k - 1), or 0 when that point is not
!> a site of the crystal (sites of different basis sites have different neighbour sets: hcp, NiAs).  Every site has nslot
!> entries.  btgt(k, b) and bdn(:, k, b) give the neighbour exactly: basis site btgt, cell n + bdn.
!>
!> Displacements.  disp(:, i) is the shortest Cartesian vector from site 1 to an image of site i.  The supercell can be skewed
!> (hexagonal cells), so the images in the neighbouring supercells are searched after wrapping the supercell coordinates into
!> [-1/2, 1/2); a tie keeps the wrapped image (for an orthogonal supercell: the half-open box of asr_lattices).
!>
!> Site operations (graph_site_ops).  The space-group operations that fix site 1 act on the graph as (n, b) ->
!> (m n + shift(:, b) - shift(:, b0), perm(b)) when they map the supercell lattice onto itself (the others are dropped).  They
!> are ordered canonically: monomial operations first, by the axis permutation (123, 132, 213, 231, 312, 321) and then by the
!> signs of the rows (bit i - 1 set: row i negative), as asr_symmetry::cubic_site_perms; the other operations follow in
!> space-group order.  The identity is operation 1.
module asr_crystal_graph
  use asr_kinds, only: dp
  use asr_crystal, only: structure_t, geo_tol, inv3
  use asr_crystal_symmetry, only: symop_t, is_monomial, basis_orbits
  use asr_augspace, only: graph_t
  implicit none
  private
  public :: crystal_graph_t, crystal_graph, graph_site, graph_reduce_cell, graph_site_ops, structure_cluster, &
            supercell_from_cart

  !> Tolerance of distances (alat units): a site at distance rcl + rtol is still inside rcl (asr_lattices, asr_strconst).
  real(dp), parameter, public :: rtol = 1.0e-9_dp

  type :: crystal_graph_t
    type(graph_t) :: g                      !! neighbour lists, nslot entries per site (0: no neighbour on that slot)
    integer :: n = 0, nslot = 0, nb = 0, b0 = 1, ncell = 0, norbit = 0
    integer :: msc(3, 3) = 0                !! supercell vectors (columns), primitive lattice coordinates
    integer :: h(3, 3) = 0                  !! upper triangular basis of the supercell lattice (cell labels)
    logical :: aliased = .false.            !! the supercell is too small: a slot of some site hits the site itself or a site
                                            !! that another slot of the same site hits already
    integer, allocatable :: cell(:, :)      !! (3, n) reduced cell of each site
    integer, allocatable :: basis(:)        !! (n) basis site of each site; the structure constants are tabulated per basis site
    integer, allocatable :: orbit(:)        !! (n) space-group orbit of the basis site (symmetry-distinct site), 1 .. norbit
    integer, allocatable :: index(:, :)     !! (0:ncell - 1, nb) graph site of (cell label, basis site)
    real(dp), allocatable :: bvec(:, :)     !! (3, nslot) bond vectors, Cartesian, alat units
    integer, allocatable :: btgt(:, :)      !! (nslot, nb) basis site of the neighbour on slot k of basis site b (0: none)
    integer, allocatable :: bdn(:, :, :)    !! (3, nslot, nb) cell offset of that neighbour
    real(dp), allocatable :: disp(:, :)     !! (3, n) minimum-image displacement of every site from site 1, Cartesian, alat units
  end type crystal_graph_t

  integer, parameter :: axperm(3, 6) = reshape([1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, 3, 2, 1], [3, 6])

contains

  !> Graph of structure s on the supercell msc with all bonds up to rcl (alat units).  ops: the space group of s (orbits).
  !> start (optional, default 1): basis site of site 1.  msg is empty on success.
  subroutine crystal_graph(s, ops, msc, rcl, cg, msg, start)
    type(structure_t), intent(in) :: s
    type(symop_t), intent(in) :: ops(:)
    integer, intent(in) :: msc(3, 3)
    real(dp), intent(in) :: rcl
    type(crystal_graph_t), intent(out) :: cg
    character(len=:), allocatable, intent(out) :: msg
    integer, intent(in), optional :: start
    real(dp), allocatable :: tv(:, :), key(:, :), skey(:, :)
    integer, allocatable :: tb(:), tt(:), tdn(:, :), slot_of(:), ord(:), lab(:), orb(:), scell(:, :), sbas(:)
    integer :: nmax(3), b, c, i, j, k, n1, n2, n3, nt, ns, det, adj(3, 3), nc(3), cnt, m(3), e1, e2, e3
    real(dp) :: f(3), v(3), p(3), q(3), d(3), dbest(3), best, dd
    msg = ''
    cg%b0 = 1; if (present(start)) cg%b0 = start
    if (cg%b0 < 1 .or. cg%b0 > s%nb) then; msg = 'crystal_graph: start is not a basis site'; return; end if
    det = idet3(msc)
    if (det == 0) then; msg = 'crystal_graph: the supercell matrix is singular'; return; end if
    cg%msc = msc; cg%nb = s%nb
    call upper_basis(msc, cg%h)
    cg%ncell = cg%h(1, 1) * cg%h(2, 2) * cg%h(3, 3)
    adj = iadj3(msc)
    ! ---- bonds of every basis site up to rcl
    do i = 1, 3; nmax(i) = ceiling(rcl * sqrt(sum(s%bvec(:, i)**2))) + 1; end do
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
    ! ---- sites, numbered by supercell coordinates relative to site 1, wrapped into [0, 1), third slowest
    cg%n = cg%ncell * s%nb
    allocate (cg%cell(3, cg%n), cg%basis(cg%n), cg%index(0:cg%ncell - 1, s%nb), skey(3, cg%n), scell(3, cg%n), sbas(cg%n))
    i = 0
    do n3 = 0, cg%h(3, 3) - 1; do n2 = 0, cg%h(2, 2) - 1; do n1 = 0, cg%h(1, 1) - 1
      do b = 1, s%nb
        i = i + 1; scell(:, i) = [n1, n2, n3]; sbas(i) = b
        p = real([n1, n2, n3], dp) + s%frac(:, b) - s%frac(:, cg%b0)
        q = matmul(real(adj, dp), p) / real(det, dp)
        q = q - floor(q + geo_tol)
        skey(:, i) = [q(3), q(2), q(1)]
      end do
    end do; end do; end do
    call sort_keys(skey, ord)
    do i = 1, cg%n
      cg%cell(:, i) = scell(:, ord(i)); cg%basis(i) = sbas(ord(i))
      cg%index(cell_label(cg%h, cg%cell(:, i)), cg%basis(i)) = i
    end do
    if (cg%basis(1) /= cg%b0 .or. any(cg%cell(:, 1) /= 0)) error stop 'crystal_graph: site 1 is not the start site'
    ! ---- neighbour lists
    cg%g%n = cg%n
    allocate (cg%g%nbr_ptr(cg%n + 1), cg%g%nbr(cg%n * ns))
    cnt = 0; cg%g%nbr_ptr(1) = 1
    do i = 1, cg%n
      b = cg%basis(i)
      do k = 1, ns
        cnt = cnt + 1; cg%g%nbr(cnt) = 0
        if (cg%btgt(k, b) == 0) cycle
        nc = graph_reduce_cell(cg, cg%cell(:, i) + cg%bdn(:, k, b))
        cg%g%nbr(cnt) = cg%index(cell_label(cg%h, nc), cg%btgt(k, b))
      end do
      cg%g%nbr_ptr(i + 1) = cg%g%nbr_ptr(i) + ns
    end do
    do i = 1, cg%n                                ! aliasing: the supercell is too small for rcl
      do k = 1, ns
        j = cg%g%nbr(cg%g%nbr_ptr(i) + k - 1)
        if (j == 0) cycle
        if (j == i .or. any(cg%g%nbr(cg%g%nbr_ptr(i):cg%g%nbr_ptr(i) + k - 2) == j)) cg%aliased = .true.
      end do
    end do
    ! ---- minimum-image displacements from site 1
    allocate (cg%disp(3, cg%n))
    do i = 1, cg%n
      p = real(cg%cell(:, i), dp) + s%frac(:, cg%basis(i)) - s%frac(:, cg%b0)
      q = matmul(real(adj, dp), p) / real(det, dp)
      m = floor(q + 0.5_dp + geo_tol)
      p = p - real(matmul(msc, m), dp)
      d = matmul(s%avec, p); dbest = d; best = sum(d * d)
      do e3 = -2, 2; do e2 = -2, 2; do e1 = -2, 2
        if (e1 == 0 .and. e2 == 0 .and. e3 == 0) cycle
        v = matmul(s%avec, p + real(matmul(msc, [e1, e2, e3]), dp))
        dd = sum(v * v)
        if (dd < best - rtol) then; best = dd; dbest = v; end if
      end do; end do; end do
      cg%disp(:, i) = dbest
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

  !> Graph site of the crystal site (n, b): cell n (any integers) and basis site b.
  pure integer function graph_site(cg, n, b)
    type(crystal_graph_t), intent(in) :: cg
    integer, intent(in) :: n(3), b
    graph_site = cg%index(cell_label(cg%h, graph_reduce_cell(cg, n)), b)
  end function graph_site

  !> The cell n reduced modulo the supercell lattice: 0 <= n_i < h(i, i).
  pure function graph_reduce_cell(cg, n) result(r)
    type(crystal_graph_t), intent(in) :: cg
    integer, intent(in) :: n(3)
    integer :: r(3), i
    r = n
    do i = 3, 1, -1
      r = r - fdiv(r(i), cg%h(i, i)) * cg%h(:, i)
    end do
  end function graph_reduce_cell

  !> The space-group operations (ops of s) that fix site 1 and map the supercell lattice onto themselves, in canonical order (see
  !> the header).  perm(site, k): image of every graph site under operation k; rot(:, :, k): its Cartesian rotation;
  !> idx(k): its index in ops.  Operation 1 is the identity.
  subroutine graph_site_ops(cg, s, ops, perm, rot, idx)
    type(crystal_graph_t), intent(in) :: cg
    type(structure_t), intent(in) :: s
    type(symop_t), intent(in) :: ops(:)
    integer, allocatable, intent(out) :: perm(:, :), idx(:)
    real(dp), allocatable, intent(out) :: rot(:, :, :)
    integer, allocatable :: cand(:), ord(:)
    real(dp), allocatable :: key(:, :)
    integer :: k, nk, i, b, b0, adj(3, 3), det, t(3, 3), n2(3)
    b0 = cg%b0
    det = idet3(cg%msc); adj = iadj3(cg%msc)
    allocate (cand(size(ops))); nk = 0
    do k = 1, size(ops)
      if (ops(k)%perm(b0) /= b0) cycle
      t = matmul(adj, matmul(ops(k)%m, cg%msc))                ! msc^-1 m msc integer: the supercell lattice is kept
      if (any(modulo(t, det) /= 0)) cycle
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
    if (any(perm(:, 1) /= [(b, b=1, cg%n)])) error stop 'graph_site_ops: operation 1 is not the identity'
    if (any(perm(1, :) /= 1)) error stop 'graph_site_ops: an operation moves the start site'
  end subroutine graph_site_ops

  !> Cluster of all crystal sites within rmax (alat units) of basis site b: vec(:, i) relative to b, Cartesian.  The origin is
  !> first, then the sites by distance, ties by z, then y, then x (the order of asr_strconst::lattice_cluster).
  !> cb (optional): basis site of every cluster site.
  subroutine structure_cluster(s, b, rmax, vec, cb)
    type(structure_t), intent(in) :: s
    integer, intent(in) :: b
    real(dp), intent(in) :: rmax
    real(dp), allocatable, intent(out) :: vec(:, :)
    integer, allocatable, intent(out), optional :: cb(:)
    real(dp), allocatable :: tv(:, :), key(:, :)
    integer, allocatable :: tc(:), ord(:)
    integer :: nmax(3), c, n1, n2, n3, n, i
    real(dp) :: v(3)
    do i = 1, 3; nmax(i) = ceiling(rmax * sqrt(sum(s%bvec(:, i)**2))) + 1; end do
    n = s%nb * product(2 * nmax + 1)
    allocate (tv(3, n), tc(n)); n = 0
    do c = 1, s%nb
      do n3 = -nmax(3), nmax(3); do n2 = -nmax(2), nmax(2); do n1 = -nmax(1), nmax(1)
        v = matmul(s%avec, real([n1, n2, n3], dp) + s%frac(:, c) - s%frac(:, b))
        if (sqrt(sum(v * v)) > rmax + rtol) cycle
        n = n + 1; tv(:, n) = v; tc(n) = c
      end do; end do; end do
    end do
    allocate (key(4, n))
    do i = 1, n; key(:, i) = [sum(tv(:, i)**2), tv(3, i), tv(2, i), tv(1, i)]; end do
    call sort_keys(key, ord)
    allocate (vec(3, n)); vec = tv(:, ord)
    if (present(cb)) cb = tc(ord)
    if (sum(vec(:, 1)**2) > 0.0_dp) error stop 'structure_cluster: the origin is not first'
  end subroutine structure_cluster

  !> Supercell matrix (primitive lattice coordinates) of the supercell with Cartesian vectors c (columns, alat units): for example
  !> L conventional cubes of bcc or fcc, c = L times the identity.  msg is empty on success.
  subroutine supercell_from_cart(s, c, msc, msg)
    type(structure_t), intent(in) :: s
    real(dp), intent(in) :: c(3, 3)
    integer, intent(out) :: msc(3, 3)
    character(len=:), allocatable, intent(out) :: msg
    real(dp) :: f(3, 3)
    msg = ''
    f = matmul(inv3(s%avec), c)
    msc = nint(f)
    if (maxval(abs(f - real(msc, dp))) > geo_tol) msg = 'supercell_from_cart: the vectors are not lattice vectors'
    if (len(msg) == 0 .and. idet3(msc) == 0) msg = 'supercell_from_cart: the vectors are linearly dependent'
  end subroutine supercell_from_cart

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

  !> Upper triangular basis h of the lattice spanned by the columns of m (integer column operations): h(:, 1) = (h11, 0, 0),
  !> h(:, 2) = (h12, h22, 0), h(:, 3) = (h13, h23, h33), with positive diagonal.
  pure subroutine upper_basis(m, h)
    integer, intent(in) :: m(3, 3)
    integer, intent(out) :: h(3, 3)
    integer :: t(3), q, j
    h = m
    do j = 1, 2
      do while (h(3, j) /= 0)
        q = h(3, 3) / h(3, j)
        h(:, 3) = h(:, 3) - q * h(:, j)
        t = h(:, 3); h(:, 3) = h(:, j); h(:, j) = t
      end do
    end do
    if (h(3, 3) < 0) h(:, 3) = -h(:, 3)
    do while (h(2, 1) /= 0)
      q = h(2, 2) / h(2, 1)
      h(:, 2) = h(:, 2) - q * h(:, 1)
      t = h(:, 2); h(:, 2) = h(:, 1); h(:, 1) = t
    end do
    if (h(2, 2) < 0) h(:, 2) = -h(:, 2)
    if (h(1, 1) < 0) h(:, 1) = -h(:, 1)
  end subroutine upper_basis

  !> Label 0 .. ncell - 1 of a reduced cell.
  pure integer function cell_label(h, n)
    integer, intent(in) :: h(3, 3), n(3)
    cell_label = n(1) + h(1, 1) * (n(2) + h(2, 2) * n(3))
  end function cell_label

  !> floor(a / b) for b > 0
  pure integer function fdiv(a, b)
    integer, intent(in) :: a, b
    fdiv = (a - modulo(a, b)) / b
  end function fdiv

  pure integer function idet3(m)
    integer, intent(in) :: m(3, 3)
    idet3 = m(1, 1) * (m(2, 2) * m(3, 3) - m(2, 3) * m(3, 2)) - m(1, 2) * (m(2, 1) * m(3, 3) - m(2, 3) * m(3, 1)) &
            + m(1, 3) * (m(2, 1) * m(3, 2) - m(2, 2) * m(3, 1))
  end function idet3

  !> Adjugate: m adj(m) = det(m) times the identity.
  pure function iadj3(m) result(a)
    integer, intent(in) :: m(3, 3)
    integer :: a(3, 3)
    a(1, 1) = m(2, 2) * m(3, 3) - m(2, 3) * m(3, 2); a(1, 2) = m(1, 3) * m(3, 2) - m(1, 2) * m(3, 3)
    a(1, 3) = m(1, 2) * m(2, 3) - m(1, 3) * m(2, 2); a(2, 1) = m(2, 3) * m(3, 1) - m(2, 1) * m(3, 3)
    a(2, 2) = m(1, 1) * m(3, 3) - m(1, 3) * m(3, 1); a(2, 3) = m(1, 3) * m(2, 1) - m(1, 1) * m(2, 3)
    a(3, 1) = m(2, 1) * m(3, 2) - m(2, 2) * m(3, 1); a(3, 2) = m(1, 2) * m(3, 1) - m(1, 1) * m(3, 2)
    a(3, 3) = m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)
  end function iadj3
end module asr_crystal_graph
