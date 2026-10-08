!> Minimal test harness (no external dependency).
module asr_testing
  use asr_kinds, only: dp
  implicit none
  private
  integer :: nfail = 0, npass = 0
  public :: check, check_close, report

contains

  subroutine check(cond, name)
    logical, intent(in) :: cond
    character(*), intent(in) :: name
    if (cond) then
      npass = npass + 1
      print '(a,a)', '  ok    ', name
    else
      nfail = nfail + 1
      print '(a,a)', '  FAIL  ', name
    end if
  end subroutine check

  subroutine check_close(val, ref, tol, name)
    real(dp), intent(in) :: val, ref, tol
    character(*), intent(in) :: name
    logical :: ok
    ok = abs(val - ref) <= tol
    if (ok) then
      npass = npass + 1
      print '(a,a,a,es9.2)', '  ok    ', name, '  err=', abs(val - ref)
    else
      nfail = nfail + 1
      print '(a,a,a,es9.2,a,es9.2)', '  FAIL  ', name, '  err=', abs(val - ref), '  tol=', tol
    end if
  end subroutine check_close

  subroutine report()
    print '(a,i0,a,i0,a)', 'passed ', npass, ', failed ', nfail, '.'
    if (nfail > 0) error stop 1
  end subroutine report
end module asr_testing
