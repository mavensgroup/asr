!> Terminators for the continued fraction.
!>
!> A terminator is (explicit extended coefficients, asymptotic tail).
!>   asymptote estimators : last | avg2 | beer-pettifor
!>   joiners              : butt (none) | luchini-nex (linear interpolation)
!> so every combination is available and can be compared on the same
!> recursion coefficients.
!>
!> Luchini-Nex here follows the published *description* (Luchini & Nex,
!> J. Phys. C 20, 3125 (1987), via secondary sources): after n computed levels
!> the coefficients are interpolated linearly from (a_n, b_n) to the
!> asymptotic (a_inf, b_inf) over ninterp further levels, then held constant.
!> NOT yet checked against the original paper.
module asr_terminator
  use asr_kinds, only: dp
  use asr_greens, only: cf_green, sqrt_tail, periodic_g
  implicit none
  private
  public :: terminator_t, make_terminator, terminated_green, bp_asymptote

  type :: terminator_t
    integer :: n = 0                       !! explicit levels incl. interpolated
    real(dp), allocatable :: a(:), b(:)
    real(dp) :: ainf = 0.0_dp, binf = 0.0_dp
    integer :: period = 1                  !! 1: sqrt tail; 2: period-2 tail
    real(dp) :: ap(2) = 0.0_dp, bp(2) = 0.0_dp, blink = 0.0_dp
    character(len=:), allocatable :: label
  end type terminator_t

contains

  !> Beer-Pettifor asymptote: choose (A,B) such that the continued fraction,
  !> terminated by the square-root tail of (A,B), has its denominator vanish
  !> exactly at both band edges A +- 2B (minimum band, no loss of weight):
  !>   A+2B = max eig of T_n with a_n -> a_n + B
  !>   A-2B = min eig of T_n with a_n -> a_n - B
  !> solved by fixed-point iteration from (A,B) = avg2 estimate.
  subroutine bp_asymptote(a, b, n, ainf, binf, niter, converged)
    real(dp), intent(in) :: a(:), b(:)
    integer, intent(in) :: n
    real(dp), intent(out) :: ainf, binf
    integer, intent(out), optional :: niter
    logical, intent(out), optional :: converged
    real(dp) :: ep, em, an, bn, d(n), e(max(n - 1, 1)), w(max(n, 1))
    integer :: it, info
    logical :: ok
    ainf = 0.5_dp * (a(n) + a(n - 1))
    binf = sqrt(0.5_dp * (b(n)**2 + b(n - 1)**2))
    ok = .false.
    do it = 1, 500
      d = a(1:n); e(1:n - 1) = b(1:n - 1)
      d(n) = a(n) + binf
      call extreme(d, e, n, ep, .true.)
      d = a(1:n); e(1:n - 1) = b(1:n - 1)
      d(n) = a(n) - binf
      call extreme(d, e, n, em, .false.)
      an = 0.5_dp * (ep + em); bn = 0.25_dp * (ep - em)
      if (abs(an - ainf) + abs(bn - binf) < 1.0e-13_dp) then
        ainf = an; binf = bn; ok = .true.; exit
      end if
      ainf = an; binf = bn
    end do
    if (present(niter)) niter = it
    if (present(converged)) converged = ok
    info = 0; w = 0.0_dp
  end subroutine bp_asymptote

  !> Largest (upper=.true.) or smallest eigenvalue of a symmetric tridiagonal.
  subroutine extreme(d, e, n, val, upper)
    real(dp), intent(inout) :: d(:), e(:)
    integer, intent(in) :: n
    real(dp), intent(out) :: val
    logical, intent(in) :: upper
    integer :: info
    external :: dsterf
    call dsterf(n, d, e, info)           ! eigenvalues ascending in d
    if (info /= 0) error stop 'bp_asymptote: dsterf failed'
    if (upper) then; val = d(n); else; val = d(1); end if
  end subroutine extreme

  !> method: 'sqrt-last' | 'sqrt-avg2' | 'beer-pettifor'
  !> join:   'butt' | 'luchini-nex' (needs ninterp > 0)
  function make_terminator(method, join, a, b, n, ninterp) result(t)
    character(len=*), intent(in) :: method, join
    real(dp), intent(in) :: a(:), b(:)
    integer, intent(in) :: n, ninterp
    type(terminator_t) :: t
    integer :: k, ni
    real(dp) :: f

    select case (method)
    case ('sqrt-last'); t%ainf = a(n); t%binf = b(n)
    case ('sqrt-avg2')
      t%ainf = 0.5_dp * (a(n) + a(n - 1)); t%binf = sqrt(0.5_dp * (b(n)**2 + b(n - 1)**2))
    case ('beer-pettifor'); call bp_asymptote(a, b, n, t%ainf, t%binf)
    case ('periodic2')
      ! levels repeat with period 2: level n+1 ~ level n-1, link b_n is exact
      t%period = 2
      t%ap = [a(n - 1), a(n)]; t%bp = [b(n - 1), b(n)]; t%blink = b(n)
      t%ainf = 0.5_dp * (a(n) + a(n - 1)); t%binf = sqrt(0.5_dp * (b(n)**2 + b(n - 1)**2))
    case default; error stop 'make_terminator: unknown method'
    end select

    ni = 0
    if (join == 'luchini-nex' .and. t%period == 1) ni = ninterp
    t%n = n + ni
    allocate (t%a(t%n), t%b(t%n))
    t%a(1:n) = a(1:n); t%b(1:n) = b(1:n)
    do k = 1, ni
      f = real(k, dp) / real(ni + 1, dp)
      t%a(n + k) = a(n) + f * (t%ainf - a(n))
      t%b(n + k) = b(n) + f * (t%binf - b(n))
    end do
    t%label = trim(method)//'+'//trim(join)
  end function make_terminator

  function terminated_green(t, z) result(g)
    type(terminator_t), intent(in) :: t
    complex(dp), intent(in) :: z
    complex(dp) :: g
    if (t%period == 1) then
      g = cf_green(t%a, t%b, t%n, z, sqrt_tail(z, t%ainf, t%binf))
    else
      g = cf_green(t%a, t%b, t%n, z, t%blink**2 * periodic_g(z, t%ap, t%bp))
    end if
  end function terminated_green
end module asr_terminator
