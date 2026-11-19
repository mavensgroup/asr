!> Nearest-neighbour single-orbital cubic-family clusters (sc/bcc/fcc) for
!> tests and the first ordered-DOS driver.  Integer coordinates in units of
!> a/2:  sc: all even; bcc: all same parity; fcc: sum even.
module asr_lattice_toy
  use asr_kinds, only: dp
  use asr_operator, only: csr_op, csr_from_triplets
  implicit none
  private
  public :: cubic_cluster, nn_vectors

contains

  subroutine nn_vectors(kind, vec)
    character(*), intent(in) :: kind
    integer, allocatable, intent(out) :: vec(:, :)
    integer :: i, j, k, c
    select case (kind)
    case ('sc')
      allocate (vec(3, 6)); c = 0
      do i = 1, 3
        do j = -1, 1, 2
          c = c + 1; vec(:, c) = 0; vec(i, c) = 2 * j
        end do
      end do
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
      error stop 'unknown lattice'
    end select
  end subroutine nn_vectors

  logical function on_lattice(kind, p)
    character(*), intent(in) :: kind
    integer, intent(in) :: p(3)
    select case (kind)
    case ('sc');  on_lattice = all(modulo(p, 2) == 0)
    case ('bcc'); on_lattice = all(modulo(p, 2) == modulo(p(1), 2))
    case ('fcc'); on_lattice = modulo(sum(p), 2) == 0
    case default; on_lattice = .false.
    end select
  end function on_lattice

  !> Open cluster of all lattice points with max|p_i| <= r (half-units),
  !> hopping -t between nearest neighbours (t=1 default sign: H = -t sum).
  !> `centre` returns the index of the origin site.
  subroutine cubic_cluster(kind, r, t, op, centre, b2_delta)
    character(*), intent(in) :: kind
    integer, intent(in) :: r
    real(dp), intent(in) :: t
    real(dp), intent(in), optional :: b2_delta   !! bcc only: on-site +delta (all-even sites) / -delta (all-odd)
    type(csr_op), intent(out) :: op
    integer, intent(out) :: centre
    integer, allocatable :: idx(:, :, :), vec(:, :), ti(:), tj(:)
    real(dp), allocatable :: tv(:)
    integer :: i, j, k, n, c, nv, q(3), nt
    logical :: b2

    b2 = present(b2_delta)
    if (b2 .and. kind /= 'bcc') error stop 'b2_delta requires bcc'
    call nn_vectors(kind, vec); nv = size(vec, 2)
    allocate (idx(-r:r, -r:r, -r:r)); idx = 0
    n = 0
    do k = -r, r; do j = -r, r; do i = -r, r
          if (on_lattice(kind, [i, j, k])) then; n = n + 1; idx(i, j, k) = n; end if
    end do; end do; end do
    centre = idx(0, 0, 0)
    allocate (ti(n * (nv + 1)), tj(n * (nv + 1)), tv(n * (nv + 1))); nt = 0
    do k = -r, r; do j = -r, r; do i = -r, r
          if (idx(i, j, k) == 0) cycle
          if (b2) then
            nt = nt + 1; ti(nt) = idx(i, j, k); tj(nt) = idx(i, j, k)
            if (modulo(i, 2) == 0) then; tv(nt) = b2_delta; else; tv(nt) = -b2_delta; end if
          end if
          do c = 1, nv
            q = [i, j, k] + vec(:, c)
            if (any(abs(q) > r)) cycle
            nt = nt + 1; ti(nt) = idx(i, j, k); tj(nt) = idx(q(1), q(2), q(3)); tv(nt) = -t
          end do
    end do; end do; end do
    call csr_from_triplets(n, ti(1:nt), tj(1:nt), tv(1:nt), op)
  end subroutine cubic_cluster
end module asr_lattice_toy
