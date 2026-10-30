!> Warren-Cowley short-range order between a site and its first coordination shell (Mookerjee & Prasad PRB 48, 17724 (1993);
!> Alam & Mookerjee JPCM 21, 195503 (2009)).  The joint probability of species c on a site of sublattice kappa and species s on a
!> first-shell neighbour on sublattice kappa' is
!>     P(c, s) = x_c y_s (1 - alpha_cs),        x, y: concentrations on kappa, kappa',
!> so the conditional probability of the neighbour given the centre is  X_s^(c) = P(s | c) = y_s (1 - alpha_cs).
!> Marginals: sum_s y_s alpha_cs = 0 for every c and sum_c x_c alpha_cs = 0 for every s; the table for (kappa', kappa) is the
!> transpose.  Unlike pair: alpha > 0 = clustering (the pair avoids nearest-neighbour contact; alpha = 1: never neighbours),
!> alpha < 0 = ordering (the pair is favoured).
!>   - one sublattice (kappa = kappa'): alpha is symmetric, the user gives the off-diagonal entries and the diagonal is derived,
!>       alpha_cc = -sum_{s /= c} x_s alpha_cs / x_c,   i.e.  X_c^(c) = x_c + sum_{s /= c} x_s alpha_cs;
!>   - two sublattices: every entry is given (unset entries are 0) and both marginal conditions are checked.
!> Every X_s^(c) must be >= 0.  Only centre-shell correlations are modelled: shell sites are independent of each other given the
!> centre, and sites beyond the first shell are uncorrelated (the model of the papers).
module asr_sro
  use asr_kinds, only: dp
  implicit none
  private
  public :: sro_conditional

  real(dp), parameter :: tol = 1.0e-10_dp

contains

  !> Conditional probabilities cond(k, k') = P(species k' on the neighbour | species k at the centre) from the concentrations of
  !> the centre's sublattice x(1:n) and of the neighbour's sublattice y(1:m) and the Warren-Cowley table alpha(n, m).  same =
  !> .true.: one sublattice (x = y), alpha must be symmetric with a zero diagonal on input; the full table (derived diagonal) is
  !> returned in alpha.  msg is empty on success.
  subroutine sro_conditional(x, y, alpha, same, cond, msg)
    real(dp), intent(in) :: x(:), y(:)
    real(dp), intent(inout) :: alpha(:, :)
    logical, intent(in) :: same
    real(dp), allocatable, intent(out) :: cond(:, :)
    character(len=:), allocatable, intent(out) :: msg
    character(len=160) :: buf
    integer :: n, m, i, j
    n = size(x); m = size(y); msg = ''
    allocate (cond(n, m)); cond = 0.0_dp
    if (size(alpha, 1) /= n .or. size(alpha, 2) /= m) then; msg = 'sro: alpha must be n(centre) x n(neighbour)'; return; end if
    if (any(x <= 0.0_dp) .or. any(y <= 0.0_dp) .or. abs(sum(x) - 1.0_dp) > 1.0e-8_dp .or. abs(sum(y) - 1.0_dp) > 1.0e-8_dp) then
      msg = 'sro: concentrations must be > 0 and sum to 1'; return
    end if
    if (same) then
      if (n /= m .or. maxval(abs(x - y)) > tol) then; msg = 'sro: one sublattice needs x = y'; return; end if
      do i = 1, n
        if (abs(alpha(i, i)) > tol) then; msg = 'sro: the diagonal alpha_cc is derived, not input'; return; end if
        do j = i + 1, n
          if (abs(alpha(i, j) - alpha(j, i)) > tol) then; msg = 'sro: alpha must be symmetric on one sublattice'; return; end if
        end do
      end do
      do i = 1, n
        alpha(i, i) = -sum(x * alpha(i, :)) / x(i)
      end do
    else
      do i = 1, n
        if (abs(sum(y * alpha(i, :))) > 1.0e-8_dp) then
          write (buf, '(a,i0,a,es10.3)') 'sro: sum_s y_s alpha(c,s) must vanish (centre species ', i, '): ', sum(y * alpha(i, :))
          msg = trim(buf); return
        end if
      end do
      do j = 1, m
        if (abs(sum(x * alpha(:, j))) > 1.0e-8_dp) then
          write (buf, '(a,i0,a,es10.3)') 'sro: sum_c x_c alpha(c,s) must vanish (neighbour species ', j, '): ', sum(x * alpha(:, j))
          msg = trim(buf); return
        end if
      end do
    end if
    do i = 1, n
      do j = 1, m
        cond(i, j) = y(j) * (1.0_dp - alpha(i, j))
        if (cond(i, j) < -tol) then
          write (buf, '(a,i0,a,i0,a,es10.3,a)') 'sro: P(', j, ' | ', i, ') = ', cond(i, j), &
            ' < 0 (alpha > 1, or the diagonal of a strongly ordering triad)'
          msg = trim(buf); return
        end if
        if (abs(cond(i, j)) <= tol) cond(i, j) = 0.0_dp
      end do
      if (abs(sum(cond(i, :)) - 1.0_dp) > 1.0e-8_dp) then; msg = 'sro: internal error, conditional sum /= 1'; return; end if
    end do
  end subroutine sro_conditional
end module asr_sro
