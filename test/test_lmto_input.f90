!> [calculation] kind = "lmto" input: parsing, defaults, strict key checking, validation.
program test_lmto_input
  use asr_kinds, only: dp
  use asr_input_lmto
  use asr_testing
  implicit none
  type(lmto_input_t) :: inp
  character(len=:), allocatable :: msg
  integer :: u
  character(len=*), parameter :: f = 'build/test_lmto_input_tmp.toml'
  character(len=1), parameter :: nl = new_line('a')
  character(len=*), parameter :: base = &
    '[calculation]'//nl//'kind = "lmto"'//nl//'[structure]'//nl//'lattice = "fcc"'//nl//'alat = 4.09'//nl// &
    '[[species]]'//nl//'symbol = "Ag"'//nl//'z = 47'//nl//'nval = 11'//nl//'konfig = [5,5,4,4]'//nl//'occupations = [1.0, 10.0]'//nl

  print '(a)', 'test_lmto_input'
  call put(base)
  call check(peek_kind(f) == 'lmto', 'peek_kind sees lmto')
  call read_lmto_input(f, inp, msg)
  call check(len(msg) == 0, 'valid lmto input accepted')
  call check(inp%species(1)%z == 47.0_dp .and. inp%species(1)%nval == 11 .and. inp%species(1)%konfig(2) == 4, 'species read')
  call check(inp%spin == 1 .and. inp%nk == 12 .and. inp%bz == 'tetra', 'scf defaults filled')
  call put(base//'[scf]'//nl//'nkk = 8'//nl)
  call read_lmto_input(f, inp, msg)
  call check(index(msg, 'unknown key') > 0, 'typo in [scf] rejected')
  call put(base//'[structure]'//nl)
  call read_lmto_input(f, inp, msg)
  call check(len(msg) > 0, 'duplicate/odd structure handled without a crash')
  call put(replace(base, 'fcc', 'hcp'))
  call read_lmto_input(f, inp, msg)
  call check(index(msg, 'lattice') > 0, 'bad lattice rejected')
  call put('[calculation]'//nl//'kind = "lmto"'//nl//'[structure]'//nl//'lattice = "b2"'//nl//'alat = 3.0'//nl)
  call read_lmto_input(f, inp, msg)
  call check(index(msg, 'species') > 0, 'missing species rejected')
  call report()
contains
  subroutine put(text)
    character(len=*), intent(in) :: text
    open (newunit=u, file=f, status='replace'); write (u, '(a)') text; close (u)
  end subroutine put
  function replace(s, a, b) result(r)
    character(len=*), intent(in) :: s, a, b
    character(len=:), allocatable :: r
    integer :: i
    i = index(s, a)
    r = s(1:i - 1)//b//s(i + len(a):)
  end function replace
end program test_lmto_input
