!> Periodic nearest-neighbour graphs (sc, bcc, fcc) and random binary-alloy realisations.
module asr_lattices
  use asr_kinds, only: dp
  use asr_augspace, only: graph_t, binary_model_t
  use asr_operator, only: csr_op, csr_from_triplets
  implicit none
  private
  public :: periodic_graph, random_alloy_op, periodic_graph_shells

contains

  !> L^3 conventional cubic cells; sites are lattice points (sc: 1/cell, bcc: 2, fcc: 4), periodic.
  !> Site index map is returned through `coord` (integer half-lattice coordinates 0..2L-1).
  subroutine periodic_graph(kind, l, g, coord)
    character(len=*), intent(in) :: kind
    integer, intent(in) :: l
    type(graph_t), intent(out) :: g
    integer, allocatable, intent(out), optional :: coord(:, :)
    integer, allocatable :: idx(:, :, :), vec(:, :), cc(:, :)
    integer :: i, j, k, n, c, nv, p(3), q(3), m, cnt, nn
    m = 2 * l
    allocate (idx(0:m - 1, 0:m - 1, 0:m - 1)); idx = 0
    n = 0
    do k = 0, m - 1; do j = 0, m - 1; do i = 0, m - 1
          if (on_lat(kind, i, j, k)) then; n = n + 1; idx(i, j, k) = n; end if
    end do; end do; end do
    call nn_vec(kind, vec); nv = size(vec, 2)
    allocate (cc(3, n))
    g%n = n
    allocate (g%nbr_ptr(n + 1), g%nbr(n * nv))
    g%nbr_ptr(1) = 1; cnt = 0
    do k = 0, m - 1; do j = 0, m - 1; do i = 0, m - 1
          if (idx(i, j, k) == 0) cycle
          cc(:, idx(i, j, k)) = [i, j, k]
          do c = 1, nv
            q = [modulo(i + vec(1, c), m), modulo(j + vec(2, c), m), modulo(k + vec(3, c), m)]
            cnt = cnt + 1; g%nbr(cnt) = idx(q(1), q(2), q(3))
          end do
    end do; end do; end do
    ! indices increase in loop order, so pointer = nv * (idx-1) + 1
    do i = 1, n; g%nbr_ptr(i + 1) = g%nbr_ptr(i) + nv; end do
    nn = 0
    if (present(coord)) then; allocate (coord(3, n)); coord = cc; end if
    p = 0
  end subroutine periodic_graph

  !> Periodic graph with every neighbour within r2max (half-lattice units squared: e.g. fcc NN = 2,
  !> 2nd = 4, 3rd = 6, 4th = 8).  Slot k of every site is the bond vector bvec(:, k) (same order everywhere).
  subroutine periodic_graph_shells(kind, l, r2max, g, coord, bvec)
    character(len=*), intent(in) :: kind
    integer, intent(in) :: l, r2max
    type(graph_t), intent(out) :: g
    integer, allocatable, intent(out) :: coord(:, :), bvec(:, :)
    integer, allocatable :: idx(:, :, :), tmp(:, :), cc(:, :), key(:)
    integer :: i, j, k, n, m, nv, v(3), c, q(3), cnt, t(3), kk, jj, rmax
    m = 2 * l
    allocate (idx(0:m - 1, 0:m - 1, 0:m - 1)); idx = 0
    n = 0
    do k = 0, m - 1; do j = 0, m - 1; do i = 0, m - 1
          if (on_lat(kind, i, j, k)) then; n = n + 1; idx(i, j, k) = n; end if
    end do; end do; end do
    rmax = int(sqrt(real(r2max))) + 1
    allocate (tmp(3, (2 * rmax + 1)**3)); nv = 0
    do k = -rmax, rmax; do j = -rmax, rmax; do i = -rmax, rmax
          if (i == 0 .and. j == 0 .and. k == 0) cycle
          if (.not. on_lat(kind, modulo(i, 2), modulo(j, 2), modulo(k, 2)) .and. kind /= 'sc') cycle
          if (kind == 'sc' .and. any(modulo([i, j, k], 2) /= 0)) cycle
          if (i * i + j * j + k * k > r2max) cycle
          nv = nv + 1; tmp(:, nv) = [i, j, k]
    end do; end do; end do
    ! sort by (length^2, then components)
    allocate (key(nv))
    do i = 1, nv; key(i) = (tmp(1, i)**2 + tmp(2, i)**2 + tmp(3, i)**2) * 100000 + (tmp(1, i) + 50) * 10000 + (tmp(2, i) + 50) * 100 + (tmp(3, i) + 50); end do
    do i = 2, nv
      kk = key(i); t = tmp(:, i); jj = i - 1
      do while (jj >= 1)
        if (key(jj) <= kk) exit
        key(jj + 1) = key(jj); tmp(:, jj + 1) = tmp(:, jj); jj = jj - 1
      end do
      key(jj + 1) = kk; tmp(:, jj + 1) = t
    end do
    allocate (bvec(3, nv)); bvec = tmp(:, 1:nv)
    g%n = n
    allocate (g%nbr_ptr(n + 1), g%nbr(n * nv), cc(3, n))
    cnt = 0
    do k = 0, m - 1; do j = 0, m - 1; do i = 0, m - 1
          if (idx(i, j, k) == 0) cycle
          cc(:, idx(i, j, k)) = [i, j, k]
          do c = 1, nv
            q = [modulo(i + bvec(1, c), m), modulo(j + bvec(2, c), m), modulo(k + bvec(3, c), m)]
            cnt = cnt + 1; g%nbr(cnt) = idx(q(1), q(2), q(3))
          end do
    end do; end do; end do
    g%nbr_ptr(1) = 1
    do i = 1, n; g%nbr_ptr(i + 1) = g%nbr_ptr(i) + nv; end do
    allocate (coord(3, n)); coord = cc
    v = 0
  end subroutine periodic_graph_shells

  pure logical function on_lat(kind, i, j, k)
    character(len=*), intent(in) :: kind
    integer, intent(in) :: i, j, k
    select case (kind)
    case ('sc'); on_lat = modulo(i, 2) == 0 .and. modulo(j, 2) == 0 .and. modulo(k, 2) == 0
    case ('bcc'); on_lat = modulo(i, 2) == modulo(j, 2) .and. modulo(j, 2) == modulo(k, 2)
    case ('fcc'); on_lat = modulo(i + j + k, 2) == 0
    case default; on_lat = .false.
    end select
  end function on_lat

  subroutine nn_vec(kind, vec)
    character(len=*), intent(in) :: kind
    integer, allocatable, intent(out) :: vec(:, :)
    integer :: i, j, k, c
    select case (kind)
    case ('sc')
      allocate (vec(3, 6)); c = 0
      do i = 1, 3; do j = -1, 1, 2
          c = c + 1; vec(:, c) = 0; vec(i, c) = 2 * j
      end do; end do
    case ('bcc')
      allocate (vec(3, 8)); c = 0
      do i = -1, 1, 2; do j = -1, 1, 2; do k = -1, 1, 2
            c = c + 1; vec(:, c) = [i, j, k]
      end do; end do; end do
    case ('fcc')
      allocate (vec(3, 12)); c = 0
      do i = -1, 1, 2; do j = -1, 1, 2
          c = c + 1; vec(:, c) = [i, j, 0]
          c = c + 1; vec(:, c) = [i, 0, j]
          c = c + 1; vec(:, c) = [0, i, j]
      end do; end do
    case default
      error stop 'nn_vec: unknown lattice'
    end select
  end subroutine nn_vec

  !> One random realisation (A with probability x) of the binary model as a CSR operator.
  subroutine random_alloy_op(g, mdl, seed, op)
    type(graph_t), intent(in) :: g
    type(binary_model_t), intent(in) :: mdl
    integer, intent(in) :: seed
    type(csr_op), intent(out) :: op
    logical, allocatable :: isa(:)
    integer, allocatable :: ti(:), tj(:), sd(:)
    real(dp), allocatable :: tv(:), r(:)
    integer :: i, j, k, nt
    allocate (isa(g%n), r(g%n))
    call random_seed(size=k); allocate (sd(k)); sd = seed + 37 * [(i, i=1, k)]; call random_seed(put=sd)
    call random_number(r)
    isa = r < mdl%x
    allocate (ti(g%n + size(g%nbr)), tj(g%n + size(g%nbr)), tv(g%n + size(g%nbr))); nt = 0
    do i = 1, g%n
      nt = nt + 1; ti(nt) = i; tj(nt) = i
      if (isa(i)) then; tv(nt) = mdl%epsa; else; tv(nt) = mdl%epsb; end if
      do k = g%nbr_ptr(i), g%nbr_ptr(i + 1) - 1
        j = g%nbr(k)
        nt = nt + 1; ti(nt) = i; tj(nt) = j
        if (isa(i) .and. isa(j)) then; tv(nt) = mdl%tAA
        else if (isa(i) .neqv. isa(j)) then; tv(nt) = mdl%tAB
        else; tv(nt) = mdl%tBB; end if
      end do
    end do
    call csr_from_triplets(g%n, ti(1:nt), tj(1:nt), tv(1:nt), op)
  end subroutine random_alloy_op
end module asr_lattices
