!> The resolved input (asr_inputdump): every key the reader used, defaults marked, method-foreign sections left out,
!> reals in their shortest form.
program test_inputdump
  use asr_driver, only: run_input_t, read_run_input
  use asr_testing
  implicit none
  type(run_input_t) :: inp
  character(len=:), allocatable :: msg, t

  print '(a)', 'test_inputdump'
  call read_run_input('example/cuzn_alloy.toml', inp, msg)
  call check(len(msg) == 0 .and. allocated(inp%resolved), 'cuzn_alloy.toml read, resolved input made')
  t = inp%resolved
  print '(a)', t
  call check(index(t, '[asr]') > 0 .and. index(t, '[scf]') > 0 .and. index(t, '[[species]]') > 0 &
             .and. index(t, '[[sublattice]]') > 0, 'sections present')
  call check(index(t, '[kspace]') == 0 .and. index(t, 'params_file') == 0, 'k-space keys left out for method = asr')
  call check(index(t, 'levels = 4') > 0 .and. index(t, 'mixing = 0.3') > 0, 'given values, shortest reals')
  call check(has_default(t, 'tolerance = 1.0E-05') .and. has_default(t, 'fermi_kt = 0.001') .and. has_default(t, 'route = "edc"'), &
             'defaults filled and marked')
  call check(.not. has_default(t, 'levels = 4'), 'a given key is not marked')
  call check(has_default(t, 'konfig = [4, 4, 3, 4]') .or. index(t, 'konfig = [4, 4, 3, 4]') > 0, 'species konfig present')

  call read_run_input('data/kspace/ag_fcc.toml', inp, msg)
  call check(len(msg) == 0 .and. index(inp%resolved, '[kspace]') > 0 .and. index(inp%resolved, '[asr]') == 0, &
             'method = kspace: [kspace] shown, [asr] left out')
  call report()
contains
  !> the line that starts with key = value carries the default mark
  logical function has_default(text, kv)
    character(len=*), intent(in) :: text, kv
    integer :: i, j
    has_default = .false.
    i = index(text, new_line('a') // kv)
    if (i == 0) return
    j = index(text(i + 1:), new_line('a'))
    if (j == 0) j = len(text) - i + 1
    has_default = index(text(i + 1:i + j - 1), '# default') > 0
  end function has_default
end program test_inputdump
