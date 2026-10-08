!> TOML input: defaults, unknown-key rejection, value validation.
program test_input
  use asr_kinds, only: dp
  use asr_input, only: input_t, read_input
  use asr_testing
  implicit none
  type(input_t) :: inp
  character(len=:), allocatable :: msg
  integer :: u
  character(len=*), parameter :: f = 'build/test_input_tmp.toml'

  print '(a)', 'test_input'
  call put('[structure]'//new_line('a')//'lattice = "bcc"'//new_line('a')//'[recursion]'//new_line('a')//'nlevels = 12')
  call read_input(f, inp, msg)
  call check(len(msg) == 0, 'valid input accepted')
  call check(inp%lattice == 'bcc' .and. inp%nlevels == 12, 'values read')
  call check(inp%eta == 0.05_dp .and. inp%npoints == 801, 'defaults filled')

  call put('[recursion]'//new_line('a')//'nlevles = 12')
  call read_input(f, inp, msg)
  call check(index(msg, 'unknown key') > 0, 'typo in key rejected')

  call put('[strukture]'//new_line('a')//'lattice = "bcc"')
  call read_input(f, inp, msg)
  call check(index(msg, 'unknown key') > 0, 'typo in section rejected')

  call put('[structure]'//new_line('a')//'lattice = "hcp"')
  call read_input(f, inp, msg)
  call check(index(msg, 'lattice') > 0, 'bad lattice rejected')

  call read_input('build/does_not_exist.toml', inp, msg)
  call check(len(msg) > 0, 'missing file reported')
  call report()
contains
  subroutine put(text)
    character(len=*), intent(in) :: text
    open (newunit=u, file=f, status='replace'); write (u, '(a)') text; close (u)
  end subroutine put
end program test_input
