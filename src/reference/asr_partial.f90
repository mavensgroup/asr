!> Species-resolved (partial) Green functions from the augmented-space recursion.
!>
!> <<n_R G_{RL,RL}>> = <u1| M_R G |u1>  with u1 = |R,{},L>  (the occupation operator n_R -> M_R).  With the
!> Lanczos basis {u_k} of u1 and v = M_R u1 this is  sum_k c_k [G]_{k1},  c_k = <u_k|v>,  where the first
!> column of the resolvent of the tridiagonal matrix follows from the continued-fraction tails:
!>   g_j = 1/(z - a_j - b_j^2 g_{j+1}),   [G]_{k1} = g_1 prod_{j<k} (b_j g_{j+1}).
!> The B-species function is  (G_avg - x G_A)/(1 - x)  with G_A = <<n G>>/x.
module asr_partial
  use asr_kinds, only: dp
  use asr_operator, only: linop
  implicit none
  private
  public :: lanczos_basis, transition_green

contains

  !> Lanczos with the n basis vectors kept (columns of u), full reorthogonalisation.
  subroutine lanczos_basis(op, v0, n, a, b, u, m)
    class(linop), intent(in) :: op
    real(dp), intent(in) :: v0(:)
    integer, intent(in) :: n
    real(dp), intent(out) :: a(n), b(n)
    real(dp), allocatable, intent(out) :: u(:, :)
    integer, intent(out) :: m
    real(dp), allocatable :: w(:)
    real(dp) :: bn
    integer :: i, dim
    dim = size(v0)
    allocate (u(dim, n), w(dim)); u = 0.0_dp; a = 0.0_dp; b = 0.0_dp
    u(:, 1) = v0 / sqrt(dot_product(v0, v0))
    m = 0
    do i = 1, n
      call op%apply(u(:, i), w)
      a(i) = dot_product(u(:, i), w)
      w = w - a(i) * u(:, i)
      if (i > 1) w = w - b(i - 1) * u(:, i - 1)
      w = w - matmul(u(:, 1:i), matmul(transpose(u(:, 1:i)), w))
      w = w - matmul(u(:, 1:i), matmul(transpose(u(:, 1:i)), w))
      m = i
      bn = sqrt(dot_product(w, w))
      b(i) = bn
      if (i == n .or. bn < 1.0e-12_dp) exit
      u(:, i + 1) = w / bn
    end do
  end subroutine lanczos_basis

  !> F(z) = sum_k c_k [G]_{k1}  for c_k = <u_k|v>  (k = 1..n), tail = terminator tail for level n.
  pure function transition_green(a, b, n, c, z, tail) result(f)
    real(dp), intent(in) :: a(:), b(:), c(:)
    integer, intent(in) :: n
    complex(dp), intent(in) :: z, tail
    complex(dp) :: f
    complex(dp) :: gt(n + 1), prod
    integer :: j, k
    gt(n + 1) = (0.0_dp, 0.0_dp)
    gt(n) = 1.0_dp / (z - a(n) - tail)
    do j = n - 1, 1, -1
      gt(j) = 1.0_dp / (z - a(j) - b(j)**2 * gt(j + 1))
    end do
    f = (0.0_dp, 0.0_dp); prod = gt(1)
    do k = 1, n
      f = f + c(k) * prod
      if (k < n) prod = prod * b(k) * gt(k + 1)
    end do
  end function transition_green
end module asr_partial
