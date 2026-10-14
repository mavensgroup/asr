!> TOML input for  [calculation] kind = "lmto":  structure + species -> SCF -> parameters -> DOS.
!>
!>   [calculation]  name, kind = "lmto"
!>   [structure]    lattice = "fcc"|"bcc"|"sc"|"b2", alat (Angstrom)  [, ratio = r_A/r_B for b2]
!>   [[species]]    symbol, z, nval, konfig = [ns, np, nd, nf] (valence principal quantum numbers),
!>                  occupations = [n_s, n_d]  (neutral-atom valence s and d occupation, starting potential)
!>                  spin_sign = 1 | -1  (sign of this species' starting spin polarisation zeta0; -1: antiparallel start)
!>   [scf]          spin = 1|2, xc = "vbh"|"pw"|"vwn", nk, bz = "tetra"|"gauss", smear, mixing,
!>                  tolerance, max_iterations, zeta0, sigma, relativistic
!>   [dos]          emin, emax, npoints, smear, nk
!>   [output]       dos_file, params_file, out_file
!> Unknown sections/keys are errors.
module asr_input_lmto
  use asr_kinds, only: dp
  use tomlf, only: toml_table, toml_array, toml_error, toml_load, get_value, len
  use asr_tomlutil, only: section, check_keys
  implicit none
  private
  public :: lmto_species_t, lmto_input_t, read_lmto_input, peek_kind

  type :: lmto_species_t
    character(len=:), allocatable :: symbol
    real(dp) :: z = 0.0_dp
    integer :: nval = 0
    integer :: konfig(0:3) = [4, 4, 3, 4]
    real(dp) :: occ_s = 2.0_dp, occ_d = 0.0_dp
    real(dp) :: spin_sign = 1.0_dp              !! sign of the starting spin polarisation (zeta0) of this species
  end type lmto_species_t

  type :: lmto_input_t
    character(len=:), allocatable :: name, lattice, xc, bz, dos_file, params_file, out_file, workdir, infile
    real(dp) :: alat_angstrom = 0.0_dp, ratio = 1.0_dp
    type(lmto_species_t), allocatable :: species(:)
    integer :: spin = 1, nk = 12, maxit = 100
    real(dp) :: smear = 0.01_dp, mixing = 0.3_dp, tolerance = 1.0e-6_dp, zeta0 = 0.4_dp, sigma = 0.7_dp
    logical :: relativistic = .true.
    real(dp) :: emin = -1.0_dp, emax = 1.0_dp, dos_smear = 0.01_dp
    integer :: npoints = 801, dos_nk = 24
  end type lmto_input_t

contains

  !> The [calculation] kind of a TOML file ('' if unreadable).
  function peek_kind(filename) result(kind)
    character(len=*), intent(in) :: filename
    character(len=:), allocatable :: kind
    type(toml_table), allocatable :: table
    type(toml_table), pointer :: child
    type(toml_error), allocatable :: perr
    kind = ''
    call toml_load(table, filename, error=perr)
    if (allocated(perr)) return
    call section(table, 'calculation', child)
    call get_value(child, 'kind', kind, 'model')
  end function peek_kind

  subroutine read_lmto_input(filename, inp, msg)
    character(len=*), intent(in) :: filename
    type(lmto_input_t), intent(out) :: inp
    character(len=:), allocatable, intent(out) :: msg
    type(toml_table), allocatable :: table
    type(toml_table), pointer :: child, sp
    type(toml_array), pointer :: arr, ka
    type(toml_error), allocatable :: perr
    integer :: i, j, slash, nsp
    real(dp) :: x

    msg = ''
    call toml_load(table, filename, error=perr)
    if (allocated(perr)) then; msg = "cannot read '"//filename//"': "//perr%message; return; end if
    inp%infile = filename
    slash = index(filename, '/', back=.true.)
    if (slash > 0) then; inp%workdir = filename(1:slash); else; inp%workdir = ''; end if
    call check_keys(table, ['calculation', 'structure  ', 'species    ', 'scf        ', 'dos        ', 'output     '], 'top level', msg)
    if (len(msg) > 0) return

    call section(table, 'calculation', child)
    call check_keys(child, ['name', 'kind'], '[calculation]', msg); if (len(msg) > 0) return
    call get_value(child, 'name', inp%name, 'asr-lmto')

    call section(table, 'structure', child)
    call check_keys(child, ['lattice', 'alat   ', 'ratio  '], '[structure]', msg); if (len(msg) > 0) return
    call get_value(child, 'lattice', inp%lattice, 'fcc')
    call get_value(child, 'alat', inp%alat_angstrom, 0.0_dp)
    call get_value(child, 'ratio', inp%ratio, 1.0_dp)

    nullify (arr)
    call get_value(table, 'species', arr, requested=.false.)
    if (.not. associated(arr)) then; msg = 'no [[species]] given'; return; end if
    nsp = len(arr)
    allocate (inp%species(nsp))
    do i = 1, nsp
      call get_value(arr, i, sp)
      call check_keys(sp, ['symbol     ', 'z          ', 'nval       ', 'konfig     ', 'occupations', 'spin_sign  '], '[[species]]', msg)
      if (len(msg) > 0) return
      call get_value(sp, 'symbol', inp%species(i)%symbol, 'X')
      call get_value(sp, 'z', inp%species(i)%z, 0.0_dp)
      call get_value(sp, 'nval', inp%species(i)%nval, 0)
      call get_value(sp, 'spin_sign', x, 1.0_dp); inp%species(i)%spin_sign = sign(1.0_dp, x)
      nullify (ka)
      call get_value(sp, 'konfig', ka, requested=.false.)
      if (associated(ka)) then
        if (len(ka) /= 4) then; msg = '[[species]] konfig needs 4 integers (s,p,d,f)'; return; end if
        do j = 1, 4; call get_value(ka, j, inp%species(i)%konfig(j - 1)); end do
      end if
      nullify (ka)
      call get_value(sp, 'occupations', ka, requested=.false.)
      if (associated(ka)) then
        if (len(ka) /= 2) then; msg = '[[species]] occupations needs 2 numbers (s,d)'; return; end if
        call get_value(ka, 1, x); inp%species(i)%occ_s = x
        call get_value(ka, 2, x); inp%species(i)%occ_d = x
      end if
    end do

    call section(table, 'scf', child)
    call check_keys(child, ['spin           ', 'xc             ', 'nk             ', 'bz             ', 'smear          ', &
                            'mixing         ', 'tolerance      ', 'max_iterations ', 'zeta0          ', 'sigma          ', &
                            'relativistic   '], '[scf]', msg)
    if (len(msg) > 0) return
    call get_value(child, 'spin', inp%spin, 1)
    call get_value(child, 'xc', inp%xc, 'vbh')
    call get_value(child, 'nk', inp%nk, 12)
    call get_value(child, 'bz', inp%bz, 'tetra')
    call get_value(child, 'smear', inp%smear, 0.01_dp)
    call get_value(child, 'mixing', inp%mixing, 0.3_dp)
    call get_value(child, 'tolerance', inp%tolerance, 1.0e-6_dp)
    call get_value(child, 'max_iterations', inp%maxit, 100)
    call get_value(child, 'zeta0', inp%zeta0, 0.4_dp)
    call get_value(child, 'sigma', inp%sigma, 0.7_dp)
    call get_value(child, 'relativistic', inp%relativistic, .true.)

    call section(table, 'dos', child)
    call check_keys(child, ['emin   ', 'emax   ', 'npoints', 'smear  ', 'nk     '], '[dos]', msg); if (len(msg) > 0) return
    call get_value(child, 'emin', inp%emin, -1.0_dp)
    call get_value(child, 'emax', inp%emax, 1.0_dp)
    call get_value(child, 'npoints', inp%npoints, 801)
    call get_value(child, 'smear', inp%dos_smear, 0.01_dp)
    call get_value(child, 'nk', inp%dos_nk, 24)

    call section(table, 'output', child)
    call check_keys(child, ['dos_file   ', 'params_file', 'out_file   '], '[output]', msg); if (len(msg) > 0) return
    call get_value(child, 'dos_file', inp%dos_file, 'dos.dat')
    call get_value(child, 'params_file', inp%params_file, 'params.dat')
    call get_value(child, 'out_file', inp%out_file, 'asr.out')

    ! ---- validation
    if (inp%lattice /= 'fcc' .and. inp%lattice /= 'bcc' .and. inp%lattice /= 'sc' .and. inp%lattice /= 'b2') &
      msg = '[structure] lattice must be fcc|bcc|sc|b2'
    if (inp%alat_angstrom <= 0.0_dp) msg = '[structure] alat (Angstrom) must be > 0'
    if (inp%lattice == 'b2' .and. nsp /= 2) msg = 'lattice b2 needs exactly two [[species]]'
    if (inp%lattice /= 'b2' .and. nsp /= 1) msg = 'lattice '//inp%lattice//' needs exactly one [[species]]'
    do i = 1, nsp
      if (inp%species(i)%z <= 0.0_dp .or. inp%species(i)%nval <= 0) msg = '[[species]] needs z and nval'
    end do
    if (inp%spin /= 1 .and. inp%spin /= 2) msg = '[scf] spin must be 1 or 2'
    if (inp%bz /= 'tetra' .and. inp%bz /= 'gauss') msg = "[scf] bz must be 'tetra' or 'gauss'"
    if (inp%nk < 4) msg = '[scf] nk must be >= 4'
    if (inp%emax <= inp%emin) msg = '[dos] emax must exceed emin'
  end subroutine read_lmto_input
end module asr_input_lmto
