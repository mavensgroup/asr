!> Anderson (Pulay/DIIS-type) mixing for fixed-point iterations x = F(x).
!>   r_k = F(x_k) - x_k ;  x_{k+1} = x_k + beta r_k - sum_j gamma_j (dX_j + beta dR_j)
!> with gamma from min | r_k - sum_j gamma_j dR_j |.  History length m.
module asr_mixing
  use asr_kinds, only: dp
  implicit none
  private
  public :: anderson_t, anderson_init, anderson_step

  type :: anderson_t
    integer :: n = 0, m = 0, nh = 0, k = 0
    real(dp), allocatable :: dx(:, :), dr(:, :), xold(:), rold(:)
  end type anderson_t

contains

  subroutine anderson_init(a, n, m)
    type(anderson_t), intent(out) :: a
    integer, intent(in) :: n, m
    a%n = n; a%m = m; a%nh = 0; a%k = 0
    allocate (a%dx(n, m), a%dr(n, m), a%xold(n), a%rold(n))
  end subroutine anderson_init

  !> xin: the input of this iteration, xout = F(xin); returns xnext.
  subroutine anderson_step(a, xin, xout, beta, xnext)
    type(anderson_t), intent(inout) :: a
    real(dp), intent(in) :: xin(:), xout(:), beta
    real(dp), intent(out) :: xnext(size(xin))
    real(dp) :: r(size(xin)), amat(a%m, a%m), rhs(a%m), gam(a%m), piv(a%m)
    integer :: j, i, nh, info, slot, ipiv(a%m)
    external :: dgesv
    r = xout - xin
    if (a%k > 0) then
      slot = mod(a%k - 1, a%m) + 1
      a%dx(:, slot) = xin - a%xold
      a%dr(:, slot) = r - a%rold
      a%nh = min(a%nh + 1, a%m)
    end if
    a%xold = xin; a%rold = r; a%k = a%k + 1
    nh = a%nh
    xnext = xin + beta * r
    if (nh == 0) return
    ! least squares  min | r - dR gamma |   via normal equations (ridge-regularised)
    do j = 1, nh
      do i = 1, nh
        amat(i, j) = dot_product(a%dr(:, i), a%dr(:, j))
      end do
      rhs(j) = dot_product(a%dr(:, j), r)
    end do
    do j = 1, nh; amat(j, j) = amat(j, j) * (1.0_dp + 1.0e-8_dp) + 1.0e-30_dp; end do
    gam = 0.0_dp
    call dgesv(nh, 1, amat(1:nh, 1:nh), nh, ipiv, rhs(1:nh), nh, info)
    if (info /= 0) return
    gam(1:nh) = rhs(1:nh)
    do j = 1, nh
      xnext = xnext - gam(j) * (a%dx(:, j) + beta * a%dr(:, j))
    end do
    piv = 0.0_dp
  end subroutine anderson_step
end module asr_mixing
