!> Runtime input: one validated record read from a TOML file.
!> Same section layout as asr2 ([calculation] [structure] [model]
!> [recursion] [dos] [output]); more sections ([[species]], [scf],
!> [augmented_space], [sro] ...) are added as the code grows.
!> Unknown sections/keys are an error, so a typo can never silently fall
!> back to a default.
module asr_input
  use asr_kinds, only: dp
  use tomlf, only: toml_table, toml_error, toml_load, get_value, toml_key
  implicit none
  private
  public :: input_t, read_input

  type :: input_t
    character(len=:), allocatable :: name, kind
    character(len=:), allocatable :: lattice
    real(dp) :: hopping = 1.0_dp
    integer :: nlevels = 40
    character(len=:), allocatable :: terminator
    real(dp) :: emin = -8.0_dp, emax = 8.0_dp, eta = 0.05_dp
    integer :: npoints = 801
    character(len=:), allocatable :: dos_file
    character(len=:), allocatable :: workdir   !! directory of the input file
  end type input_t

contains

  !> Read and validate; on failure `msg` is non-empty.
  subroutine read_input(filename, inp, msg)
    character(len=*), intent(in) :: filename
    type(input_t), intent(out) :: inp
    character(len=:), allocatable, intent(out) :: msg
    type(toml_table), allocatable :: table
    type(toml_table), pointer :: child
    type(toml_error), allocatable :: perr
    integer :: slash

    msg = ''
    call toml_load(table, filename, error=perr)
    if (allocated(perr)) then
      msg = "cannot read '"//filename//"': "//perr%message; return
    end if
    slash = index(filename, '/', back=.true.)
    if (slash > 0) then; inp%workdir = filename(1:slash); else; inp%workdir = ''; end if

    call check_keys(table, ['calculation', 'structure  ', 'model      ', 'recursion  ', 'dos        ', 'output     '], 'top level', msg)
    if (len(msg) > 0) return

    call section(table, 'calculation', child)
    call check_keys(child, ['name', 'kind'], '[calculation]', msg); if (len(msg) > 0) return
    call get_value(child, 'name', inp%name, 'asr')
    call get_value(child, 'kind', inp%kind, 'model')

    call section(table, 'structure', child)
    call check_keys(child, ['lattice'], '[structure]', msg); if (len(msg) > 0) return
    call get_value(child, 'lattice', inp%lattice, 'sc')

    call section(table, 'model', child)
    call check_keys(child, ['hopping'], '[model]', msg); if (len(msg) > 0) return
    call get_value(child, 'hopping', inp%hopping, 1.0_dp)

    call section(table, 'recursion', child)
    call check_keys(child, ['nlevels   ', 'terminator'], '[recursion]', msg); if (len(msg) > 0) return
    call get_value(child, 'nlevels', inp%nlevels, 40)
    call get_value(child, 'terminator', inp%terminator, 'sqrt-avg2')

    call section(table, 'dos', child)
    call check_keys(child, ['emin   ', 'emax   ', 'eta    ', 'npoints'], '[dos]', msg); if (len(msg) > 0) return
    call get_value(child, 'emin', inp%emin, -8.0_dp)
    call get_value(child, 'emax', inp%emax, 8.0_dp)
    call get_value(child, 'eta', inp%eta, 0.05_dp)
    call get_value(child, 'npoints', inp%npoints, 801)

    call section(table, 'output', child)
    call check_keys(child, ['dos_file'], '[output]', msg); if (len(msg) > 0) return
    call get_value(child, 'dos_file', inp%dos_file, 'dos.dat')

    ! ---- value validation
    if (inp%kind /= 'model') msg = "[calculation] kind '"//inp%kind//"' not supported yet (only 'model')"
    if (inp%lattice /= 'sc' .and. inp%lattice /= 'bcc' .and. inp%lattice /= 'fcc') msg = "[structure] lattice must be sc|bcc|fcc"
    if (inp%nlevels < 2) msg = '[recursion] nlevels must be >= 2'
    if (inp%eta <= 0.0_dp) msg = '[dos] eta must be > 0'
    if (inp%emax <= inp%emin) msg = '[dos] emax must exceed emin'
    if (inp%npoints < 2) msg = '[dos] npoints must be >= 2'
    if (inp%terminator /= 'sqrt-last' .and. inp%terminator /= 'sqrt-avg2') msg = "[recursion] terminator must be sqrt-last|sqrt-avg2"
  end subroutine read_input

  !> Get a sub-table; an absent section is an empty one (all defaults).
  subroutine section(table, name, child)
    type(toml_table), intent(inout) :: table
    character(len=*), intent(in) :: name
    type(toml_table), pointer, intent(out) :: child
    nullify (child)
    call get_value(table, name, child, requested=.false.)
    if (.not. associated(child)) call get_value(table, name, child)  ! create empty
  end subroutine section

  subroutine check_keys(table, allowed, where, msg)
    type(toml_table), intent(inout) :: table
    character(len=*), intent(in) :: allowed(:), where
    character(len=:), allocatable, intent(inout) :: msg
    type(toml_key), allocatable :: keys(:)
    integer :: i, j
    logical :: known
    call table%get_keys(keys)
    do i = 1, size(keys)
      known = .false.
      do j = 1, size(allowed)
        if (trim(keys(i)%key) == trim(allowed(j))) known = .true.
      end do
      if (.not. known) then
        msg = 'unknown key "'//trim(keys(i)%key)//'" in '//where; return
      end if
    end do
  end subroutine check_keys
end module asr_input
