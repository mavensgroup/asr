!> Anderson (Pulay/DIIS-type) mixing for fixed-point iterations x = F(x).
!>   r_k = F(x_k) - x_k ;  x_{k+1} = x_k + beta r_k - sum_j gamma_j (dX_j + beta dR_j)
!> with gamma from min | r_k - sum_j gamma_j dR_j |.  History length m.
!> Safeguard: the history is dropped (plain linear step x_k + beta r_k) when the Anderson step points against the current
!> residual, (x_{k+1} - x_k).r_k <= 0, or when |r_k| grew by more than a factor 2: in strongly non-linear early iterations
!> (magnetic starts) the extrapolation otherwise marched the potential away from F(x) at the safety-clamp rate.
!> safeguard = .false. (the caller's choice, e.g. once close to self-consistency where the map is linear) keeps the history: near
!> convergence the two tests fired every few iterations and reduced Anderson to linear mixing, which crawls along soft modes.
!> After a restart the history is refilled from column 1 (nfill counts the pairs since the restart): the fit uses columns 1..nh,
!> so writing on at the old ring position left the newest pair out and fitted stale pairs from before the restart.
module asr_mixing
  use asr_kinds, only: dp
  implicit none
  private
  public :: anderson_t, anderson_init, anderson_step

  type :: anderson_t
    integer :: n = 0, m = 0, nh = 0, k = 0, nfill = 0
    integer :: nrestart = 0                     !! restarts so far (diagnostic)
    real(dp), allocatable :: dx(:, :), dr(:, :), xold(:), rold(:)
  end type anderson_t

contains

  subroutine anderson_init(a, n, m)
    type(anderson_t), intent(out) :: a
    integer, intent(in) :: n, m
    a%n = n; a%m = m; a%nh = 0; a%k = 0; a%nfill = 0; a%nrestart = 0
    allocate (a%dx(n, m), a%dr(n, m), a%xold(n), a%rold(n))
  end subroutine anderson_init

  !> xin: the input of this iteration, xout = F(xin); returns xnext.
  !> w (optional): weights of the inner products (the least-squares metric); default 1.
  subroutine anderson_step(a, xin, xout, beta, xnext, w, safeguard)
    type(anderson_t), intent(inout) :: a
    real(dp), intent(in) :: xin(:), xout(:), beta
    real(dp), intent(out) :: xnext(size(xin))
    real(dp), intent(in), optional :: w(:)
    logical, intent(in), optional :: safeguard
    real(dp) :: r(size(xin)), amat(a%m, a%m), rhs(a%m), gam(a%m), piv(a%m), rprev, wt(size(xin))
    integer :: j, i, nh, info, slot, ipiv(a%m)
    logical :: guard
    external :: dgesv
    r = xout - xin
    wt = 1.0_dp; if (present(w)) wt = w
    guard = .true.; if (present(safeguard)) guard = safeguard
    rprev = huge(1.0_dp); if (a%k > 0) rprev = sqrt(sum(wt * a%rold**2))
    if (guard .and. a%k > 0 .and. sqrt(sum(wt * r**2)) > 2.0_dp * rprev) call restart()    ! residual blew up: restart the history
    if (a%k > 0) then
      slot = mod(a%nfill, a%m) + 1
      a%dx(:, slot) = xin - a%xold
      a%dr(:, slot) = r - a%rold
      a%nfill = a%nfill + 1
      a%nh = min(a%nfill, a%m)
    end if
    a%xold = xin; a%rold = r; a%k = a%k + 1
    nh = a%nh
    xnext = xin + beta * r
    if (nh == 0) return
    ! least squares  min | r - dR gamma |   via normal equations (ridge-regularised)
    do j = 1, nh
      do i = 1, nh
        amat(i, j) = sum(wt * a%dr(:, i) * a%dr(:, j))
      end do
      rhs(j) = sum(wt * a%dr(:, j) * r)
    end do
    do j = 1, nh; amat(j, j) = amat(j, j) * (1.0_dp + 1.0e-8_dp) + 1.0e-30_dp; end do
    gam = 0.0_dp
    call dgesv(nh, 1, amat(1:nh, 1:nh), nh, ipiv, rhs(1:nh), nh, info)
    if (info /= 0) return
    gam(1:nh) = rhs(1:nh)
    do j = 1, nh
      xnext = xnext - gam(j) * (a%dx(:, j) + beta * a%dr(:, j))
    end do
    if (guard .and. sum(wt * (xnext - xin) * r) <= 0.0_dp) then                      ! step against the residual: linear step, new history
      xnext = xin + beta * r
      call restart()
    end if
    piv = 0.0_dp
  contains
    subroutine restart()
      a%nh = 0; a%nfill = 0; a%nrestart = a%nrestart + 1
    end subroutine restart
  end subroutine anderson_step
end module asr_mixing
