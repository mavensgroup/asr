!> Anderson (Pulay/DIIS-type) mixing for fixed-point iterations x = F(x).
!>   r_k = F(x_k) - x_k ;  x_{k+1} = x_k + beta r_k - sum_j gamma_j (dX_j + beta dR_j)
!> with gamma from min | r_k - sum_j gamma_j dR_j |.  History length m.
!> Safeguard: the history is dropped when |r_k| grew by more than a factor 2 (the newest pair is kept), and when the Anderson step
!> points against the current residual, (x_{k+1} - x_k).r_k <= 0: then the step is first retried with the newest pair alone, and
!> only if that also points against r_k it is the linear step x_k + beta_restart r_k with an empty history.  In strongly
!> non-linear early iterations (magnetic starts) the extrapolation otherwise marched the potential away from F(x) at the
!> safety-clamp rate.  beta_restart (default beta) may be smaller than beta: a stiff map (charge sloshing between alloy
!> components, Fe66Ni14Cr20 levels 4: gain ~45) is unstable under the linear step at beta = 0.2 (dV x8 after every restart, a
!> limit cycle), while Anderson with a history controls it.
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
  subroutine anderson_step(a, xin, xout, beta, xnext, w, safeguard, beta_restart)
    type(anderson_t), intent(inout) :: a
    real(dp), intent(in) :: xin(:), xout(:), beta
    real(dp), intent(out) :: xnext(size(xin))
    real(dp), intent(in), optional :: w(:)
    logical, intent(in), optional :: safeguard
    real(dp), intent(in), optional :: beta_restart
    real(dp) :: r(size(xin)), rprev, wt(size(xin)), br
    integer :: slot
    logical :: guard
    r = xout - xin
    wt = 1.0_dp; if (present(w)) wt = w
    guard = .true.; if (present(safeguard)) guard = safeguard
    br = beta; if (present(beta_restart)) br = beta_restart
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
    call fit_step()
    if (.not. guard .or. a%nh == 0) return
    if (sum(wt * (xnext - xin) * r) > 0.0_dp) return
    if (a%nh > 1) then                                                 ! step against the residual: the newest pair alone
      slot = mod(a%nfill - 1, a%m) + 1
      a%dx(:, 1) = a%dx(:, slot); a%dr(:, 1) = a%dr(:, slot)
      call restart(); a%nfill = 1; a%nh = 1
      call fit_step()
      if (sum(wt * (xnext - xin) * r) > 0.0_dp) return
    end if
    xnext = xin + br * r                                               ! still against it: linear step, empty history
    call restart()
  contains
    subroutine restart()
      a%nh = 0; a%nfill = 0; a%nrestart = a%nrestart + 1
    end subroutine restart

    !> xnext from the pairs 1..nh (least squares min | r - dR gamma | via ridge-regularised normal equations)
    subroutine fit_step()
      real(dp) :: amat(a%m, a%m), rhs(a%m)
      integer :: i, j, nh, info, ipiv(a%m)
      external :: dgesv
      nh = a%nh
      xnext = xin + beta * r
      if (nh == 0) return
      do j = 1, nh
        do i = 1, nh
          amat(i, j) = sum(wt * a%dr(:, i) * a%dr(:, j))
        end do
        rhs(j) = sum(wt * a%dr(:, j) * r)
      end do
      do j = 1, nh; amat(j, j) = amat(j, j) * (1.0_dp + 1.0e-8_dp) + 1.0e-30_dp; end do
      call dgesv(nh, 1, amat(1:nh, 1:nh), nh, ipiv, rhs(1:nh), nh, info)
      if (info /= 0) return
      do j = 1, nh
        xnext = xnext - rhs(j) * (a%dx(:, j) + beta * a%dr(:, j))
      end do
    end subroutine fit_step
  end subroutine anderson_step
end module asr_mixing
