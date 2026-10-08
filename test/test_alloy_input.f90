!> [calculation] kind = "alloy" input: parsing, defaults, strict key checking, validation.
program test_alloy_input
  use asr_kinds, only: dp
  use asr_input_lmto, only: peek_kind
  use asr_alloy_driver
  use asr_testing
  implicit none
  type(alloy_run_input_t) :: inp
  character(len=:), allocatable :: msg
  integer :: u
  character(len=*), parameter :: f = 'build/test_alloy_input_tmp.toml'
  character(len=1), parameter :: nl = new_line('a')
  character(len=*), parameter :: sp1 = '[[species]]'//nl//'symbol = "Cu"'//nl//'z = 29'//nl//'nval = 11'//nl//'concentration = 0.3'//nl
  character(len=*), parameter :: sp2 = '[[species]]'//nl//'symbol = "Zn"'//nl//'z = 30'//nl//'nval = 12'//nl
  character(len=*), parameter :: base = '[calculation]'//nl//'kind = "alloy"'//nl//'[structure]'//nl//'lattice = "bcc"'//nl//'alat = 2.85'//nl//sp1//sp2

  print '(a)', 'test_alloy_input'
  call put(base)
  call check(peek_kind(f) == 'alloy', 'peek_kind sees alloy')
  call read_alloy_input(f, inp, msg)
  call check(len(msg) == 0, 'valid alloy input accepted')
  call check(abs(inp%ai%x - 0.3_dp) < 1.0e-12_dp .and. inp%species(2)%nval == 12, 'concentration and species read')
  call check(inp%ai%nlevels == 4 .and. inp%ai%symmetry .and. inp%ai%nspin == 1 .and. .not. inp%ai%energy_dependent, 'defaults filled (energy-independent recursion by default until the ED route is validated for alloys)')
  call put(base//'[asr]'//nl//'levls = 3'//nl)
  call read_alloy_input(f, inp, msg)
  call check(index(msg, 'unknown key') > 0, 'typo in [asr] rejected')
  call put('[calculation]'//nl//'kind = "alloy"'//nl//'[structure]'//nl//'lattice = "bcc"'//nl//'alat = 2.85'//nl//sp1)
  call read_alloy_input(f, inp, msg)
  call check(index(msg, 'two') > 0, 'one species rejected')
  call put('[calculation]'//nl//'kind = "alloy"'//nl//'[structure]'//nl//'lattice = "bcc"'//nl//'alat = 2.85'//nl//sp2//sp2)
  call read_alloy_input(f, inp, msg)
  call check(index(msg, 'concentration') > 0, 'missing concentration rejected')
  call report()
contains
  subroutine put(text)
    character(len=*), intent(in) :: text
    open (newunit=u, file=f, status='replace'); write (u, '(a)') text; close (u)
  end subroutine put
end program test_alloy_input
