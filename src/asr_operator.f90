!> Abstract Hermitian (real symmetric for now) linear operator.
!> Everything the recursion needs from a Hamiltonian is `apply`.
!> LMTO, Wannier, plain TB and the augmented-space Hamiltonian all extend it.
module asr_operator
  use asr_kinds, only: dp
  implicit none
  private
  public :: linop, csr_op, csr_from_triplets

  type, abstract :: linop
  contains
    procedure(apply_if), deferred :: apply   ! y = H x
    procedure(dim_if),   deferred :: dim
    procedure :: mapply => linop_mapply          ! y = M x, the configuration projector (augmented-space operators only)
  end type linop

  abstract interface
    subroutine apply_if(self, x, y)
      import :: linop, dp
      class(linop), intent(in) :: self
      real(dp), intent(in)  :: x(:)
      real(dp), intent(out) :: y(:)
    end subroutine apply_if
    pure integer function dim_if(self)
      import :: linop
      class(linop), intent(in) :: self
    end function dim_if
  end interface

  !> Sparse real symmetric/general matrix in CSR (stores all entries).
  type, extends(linop) :: csr_op
    integer :: n = 0
    integer, allocatable :: ia(:), ja(:)
    real(dp), allocatable :: a(:)
  contains
    procedure :: apply => csr_apply
    procedure :: dim => csr_dim
  end type csr_op

contains

  subroutine linop_mapply(self, x, y)
    class(linop), intent(in) :: self
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    y = 0.0_dp
    error stop 'linop%mapply: not an augmented-space operator'
  end subroutine linop_mapply

  subroutine csr_apply(self, x, y)
    class(csr_op), intent(in) :: self
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    integer :: i, k
    real(dp) :: s
    !$omp parallel do private(k, s) schedule(static)
    do i = 1, self%n
      s = 0.0_dp
      do k = self%ia(i), self%ia(i + 1) - 1
        s = s + self%a(k) * x(self%ja(k))
      end do
      y(i) = s
    end do
    !$omp end parallel do
  end subroutine csr_apply

  pure integer function csr_dim(self)
    class(csr_op), intent(in) :: self
    csr_dim = self%n
  end function csr_dim

  !> Build CSR from (i,j,v) triplets; duplicates are summed.
  subroutine csr_from_triplets(n, ti, tj, tv, op)
    integer, intent(in) :: n, ti(:), tj(:)
    real(dp), intent(in) :: tv(:)
    type(csr_op), intent(out) :: op
    integer :: nt, k, i, p
    integer, allocatable :: cnt(:), pos(:)
    nt = size(ti)
    allocate (cnt(n + 1), source=0)
    do k = 1, nt
      cnt(ti(k) + 1) = cnt(ti(k) + 1) + 1
    end do
    op%n = n
    allocate (op%ia(n + 1), op%ja(nt), op%a(nt))
    op%ia(1) = 1
    do i = 1, n
      op%ia(i + 1) = op%ia(i) + cnt(i + 1)
    end do
    allocate (pos(n)); pos = op%ia(1:n)
    do k = 1, nt
      p = pos(ti(k)); op%ja(p) = tj(k); op%a(p) = tv(k); pos(ti(k)) = p + 1
    end do
  end subroutine csr_from_triplets
end module asr_operator
