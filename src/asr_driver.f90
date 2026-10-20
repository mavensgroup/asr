!> The one input format and driver of `asr <input.toml>`: a structure + species, run self-consistently by
!>   method = "asr"    (default) energy-dependent augmented-space recursion (asr_class_scf, asr_ed form): any decorated structure,
!>                     ordered or with disordered sublattices (pure metal, B2, Heusler, partial disorder, binary .. HEA alloys);
!>   method = "kspace" TB-LMTO-ASA in k-space (asr_lmto_driver), ordered structures only: the exact reference of the same model.
!>
!>   [calculation]  name, method = "asr" | "kspace"
!>   [structure]    lattice = "sc"|"bcc"|"fcc" (the graph lattice: every site of it belongs to one sublattice), alat (Angstrom, of the
!>                  graph lattice), period (decoration cell = period^3 conventional cubes; default 1), centering = "P"|"I"|"F"
!>                  (translations of the decoration cell applied to every position; default "P")
!>                  or lattice = "cell" (method = "kspace" only; e.g. SQS supercells, tools/sqs.py): alat, vectors = [[a1], [a2], [a3]]
!>                  (alat units), site_classes = true|false (one potential per site instead of per species; default false)
!>   [[species]]    symbol, z, nval, konfig = [ns, np, nd, nf] (valence principal quantum numbers), occupations = [n_s, n_d]
!>                  (neutral-atom start), spin_sign = 1 | -1 (sign of this species' start polarisation)
!>   [[sublattice]] name, positions = [[x, y, z], ...] (fractional coordinates of the decoration cell), species = ["A", ...],
!>                  concentrations = [x_A, ...] (sum 1; may be omitted for one species)
!>   [[site]]       (lattice = "cell") position = [x, y, z] (Cartesian, alat units), species = "<symbol>"
!>   [scf]          spin, xc, mixing, tolerance, max_iterations, zeta0, relativistic; method = "asr" also: sim_beta,
!>                  save_potential / start_potential
!>   [asr]          (method = "asr") levels, recursion_steps, rank, supercell (multiple of period), flip_radius, contour_points,
!>                  cluster_radius, inversion_radius, ed_terminator (hull|bp|avg2), fermi_kt, fermi_poles, depth_tol
!>   [kspace]       (method = "kspace") nk, bz = "tetra"|"gauss", smear, sigma, trace, dos_nk, sphere_ratio (r_A / r_B of a B2
!>                  structure; default 1 = equal spheres)
!>   [dos]          emin, emax, npoints, eta (recursion: imaginary part; k-space: smearing)
!>   [output]       dos_file, out_file, params_file (method = "kspace")
!> Examples: pure bcc Fe = lattice "bcc", centering "I", one sublattice [[0,0,0]] with species ["Fe"]; B2 FeCo = "bcc", period 1,
!> centering "P", Fe (0,0,0) and Co (1/2,1/2,1/2); L2_1 Co2Fe0.4Cr0.6Al on the bcc graph (a_bcc = a/2): period = 2, centering = "F";
!> Al (0,0,0); Fe/Cr (1/2,1/2,1/2); Co (1/4,1/4,1/4) and (3/4,3/4,3/4).  Unknown sections/keys are errors.
module asr_driver
  use asr_kinds, only: dp
  use tomlf, only: toml_table, toml_array, toml_error, toml_load, get_value, len
  use asr_tomlutil, only: section, check_keys
  use asr_input_lmto, only: lmto_species_t, lmto_input_t
  use asr_asa, only: core_shells_for, asa_class_result_t
  use asr_lmto_driver, only: run_lmto, atom_shells
  use asr_sublattice, only: decoration_t, decoration_build
  use asr_class_scf
  use asr_report
  !$ use omp_lib, only: omp_get_max_threads
  implicit none
  private
  public :: run_input_t, read_run_input, run_input, kspace_cell

  type :: run_input_t
    character(len=:), allocatable :: name, method, dos_file, out_file, workdir, infile
    type(class_scf_input_t) :: ci                      !! method = "asr"
    type(lmto_input_t) :: li                           !! method = "kspace"
    type(lmto_species_t), allocatable :: species(:)
    character(len=16), allocatable :: subname(:)
    real(dp) :: emin = -1.0_dp, emax = 1.0_dp, eta = 0.02_dp
    integer :: npoints = 801
  end type run_input_t

  real(dp), parameter :: bohr_per_angstrom = 1.0_dp / 0.529177210903_dp

contains

  subroutine read_run_input(filename, inp, msg)
    character(len=*), intent(in) :: filename
    type(run_input_t), intent(out) :: inp
    character(len=:), allocatable, intent(out) :: msg
    type(toml_table), allocatable :: table
    type(toml_table), pointer :: child, sp
    type(toml_array), pointer :: arr, ka, kb
    type(toml_error), allocatable :: perr
    integer :: i, j, k, slash, nsp, nsub, npos, nsc, period, nmax
    real(dp) :: x, alat, ratio
    real(dp), allocatable :: pos(:, :), tpos(:, :)
    integer, allocatable :: owner(:), towner(:)
    character(len=:), allocatable :: lat, xcs, edt, cen, str
    logical :: kspace, cell, decor

    msg = ''
    call toml_load(table, filename, error=perr)
    if (allocated(perr)) then; msg = "cannot read '"//filename//"': "//perr%message; return; end if
    inp%infile = filename
    slash = index(filename, '/', back=.true.)
    if (slash > 0) then; inp%workdir = filename(1:slash); else; inp%workdir = ''; end if
    call check_keys(table, ['calculation', 'structure  ', 'species    ', 'sublattice ', 'site       ', 'scf        ', 'asr        ', &
                            'kspace     ', 'dos        ', 'output     '], 'top level', msg)
    if (len(msg) > 0) return

    call section(table, 'calculation', child)
    if (child%has_key('kind')) then
      msg = '[calculation] kind is no longer read: every system uses one input format, run by method = "asr" (default) or, '// &
            'for an ordered structure, method = "kspace" (see example/)'
      return
    end if
    call check_keys(child, ['name  ', 'method'], '[calculation]', msg); if (len(msg) > 0) return
    call get_value(child, 'name', inp%name, 'asr')
    call get_value(child, 'method', str, 'asr'); inp%method = str
    if (inp%method /= 'asr' .and. inp%method /= 'kspace') then; msg = '[calculation] method must be "asr" or "kspace"'; return; end if
    kspace = inp%method == 'kspace'

    call section(table, 'structure', child)
    call check_keys(child, ['lattice     ', 'alat        ', 'period      ', 'centering   ', 'vectors     ', 'site_classes'], &
                    '[structure]', msg)
    if (len(msg) > 0) return
    decor = child%has_key('period') .or. child%has_key('centering')      ! (get_value with a default adds the key)
    call get_value(child, 'lattice', lat, 'bcc')
    call get_value(child, 'alat', alat, 0.0_dp)
    call get_value(child, 'period', period, 1)
    call get_value(child, 'centering', cen, 'P')
    cell = lat == 'cell'
    if (lat /= 'sc' .and. lat /= 'bcc' .and. lat /= 'fcc' .and. .not. cell) then
      msg = '[structure] lattice must be sc|bcc|fcc (decorated, [[sublattice]]) or cell ([[site]])'; return
    end if
    if (alat <= 0.0_dp) then; msg = '[structure] alat (Angstrom) must be > 0'; return; end if
    if (period < 1) then; msg = '[structure] period must be >= 1'; return; end if
    if (cell) then
      if (decor) then
        msg = '[structure] period / centering belong to a decorated lattice, not to lattice = "cell"'; return
      end if
      if (.not. kspace) then
        msg = 'lattice = "cell" (explicit sites) runs with [calculation] method = "kspace" only: the recursion needs a decorated '// &
              'sc|bcc|fcc lattice with [[sublattice]] entries'
        return
      end if
      nullify (arr)
      call get_value(child, 'vectors', arr, requested=.false.)
      if (.not. associated(arr)) then; msg = 'lattice "cell" needs [structure] vectors = [[a1], [a2], [a3]]'; return; end if
      if (len(arr) /= 3) then; msg = '[structure] vectors needs three vectors'; return; end if
      do i = 1, 3
        nullify (ka); call get_value(arr, i, ka)
        if (.not. associated(ka)) then; msg = '[structure] vectors: each vector is an array of 3 numbers'; return; end if
        if (len(ka) /= 3) then; msg = '[structure] vectors: each vector is an array of 3 numbers'; return; end if
        do j = 1, 3; call get_value(ka, j, x); inp%li%cell_vec(j, i) = x; end do
      end do
      call get_value(child, 'site_classes', inp%li%site_classes, .false.)
    else
      if (child%has_key('vectors') .or. child%has_key('site_classes')) then
        msg = '[structure] vectors / site_classes need lattice = "cell"'; return
      end if
    end if
    inp%ci%lattice = lat; inp%ci%alat = alat * bohr_per_angstrom

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
      if (inp%species(i)%z <= 0.0_dp .or. inp%species(i)%nval <= 0) then; msg = '[[species]] needs z and nval'; return; end if
      do j = 1, i - 1
        if (inp%species(j)%symbol == inp%species(i)%symbol) then; msg = 'duplicate [[species]] symbol '//inp%species(i)%symbol; return; end if
      end do
    end do

    if (cell) then
      if (table%has_key('sublattice')) then; msg = 'lattice = "cell" takes [[site]] entries, not [[sublattice]]'; return; end if
      nullify (arr)
      call get_value(table, 'site', arr, requested=.false.)
      if (.not. associated(arr)) then; msg = 'lattice "cell" needs [[site]] entries'; return; end if
      allocate (inp%li%site_pos(3, len(arr)), inp%li%site_sp(len(arr)))
      do i = 1, len(arr)
        call get_value(arr, i, sp)
        call check_keys(sp, ['position', 'species '], '[[site]]', msg); if (len(msg) > 0) return
        nullify (ka); call get_value(sp, 'position', ka, requested=.false.)
        if (.not. associated(ka)) then; msg = '[[site]] needs position = [x, y, z]'; return; end if
        if (len(ka) /= 3) then; msg = '[[site]] position needs 3 numbers'; return; end if
        do j = 1, 3; call get_value(ka, j, x); inp%li%site_pos(j, i) = x; end do
        call get_value(sp, 'species', str, '')
        inp%li%site_sp(i) = 0
        do k = 1, nsp
          if (inp%species(k)%symbol == str) inp%li%site_sp(i) = k
        end do
        if (inp%li%site_sp(i) == 0) then; msg = '[[site]] species "'//str//'" is not a [[species]] symbol'; return; end if
      end do
    else
      if (table%has_key('site')) then; msg = '[[site]] entries need lattice = "cell"; a decorated lattice takes [[sublattice]]'; return; end if
      nullify (arr)
      call get_value(table, 'sublattice', arr, requested=.false.)
      if (.not. associated(arr)) then; msg = 'no [[sublattice]] given'; return; end if
      nsub = len(arr)
      allocate (inp%subname(nsub), inp%ci%nsc(nsub), inp%ci%spc(nsp, nsub), inp%ci%conc(nsp, nsub), pos(3, 0), owner(0))
      inp%ci%spc = 0; inp%ci%conc = 0.0_dp
      do i = 1, nsub
        call get_value(arr, i, sp)
        call check_keys(sp, ['name          ', 'positions     ', 'species       ', 'concentrations'], '[[sublattice]]', msg)
        if (len(msg) > 0) return
        call get_value(sp, 'name', str, 'S'); inp%subname(i) = str
        nullify (ka)
        call get_value(sp, 'positions', ka, requested=.false.)
        if (.not. associated(ka)) then; msg = '[[sublattice]] needs positions'; return; end if
        npos = len(ka)
        do j = 1, npos
          nullify (kb)
          call get_value(ka, j, kb)
          if (.not. associated(kb)) then; msg = '[[sublattice]] positions must be [[x, y, z], ...]'; return; end if
          if (len(kb) /= 3) then; msg = '[[sublattice]] positions must be [[x, y, z], ...]'; return; end if
          allocate (tpos(3, size(pos, 2) + 1), towner(size(owner) + 1))
          tpos(:, 1:size(pos, 2)) = pos; towner(1:size(owner)) = owner
          do k = 1, 3; call get_value(kb, k, x); tpos(k, size(tpos, 2)) = x; end do
          towner(size(towner)) = i
          call move_alloc(tpos, pos); call move_alloc(towner, owner)
        end do
        nullify (ka)
        call get_value(sp, 'species', ka, requested=.false.)
        if (.not. associated(ka)) then; msg = '[[sublattice]] needs species'; return; end if
        nsc = len(ka)
        if (nsc < 1 .or. nsc > nsp) then; msg = '[[sublattice]] species: 1 .. number of [[species]] entries'; return; end if
        inp%ci%nsc(i) = nsc
        do j = 1, nsc
          call get_value(ka, j, str)
          do k = 1, nsp
            if (inp%species(k)%symbol == str) inp%ci%spc(j, i) = k
          end do
          if (inp%ci%spc(j, i) == 0) then; msg = '[[sublattice]] species '//str//' is not a [[species]] symbol'; return; end if
          if (any(inp%ci%spc(1:j - 1, i) == inp%ci%spc(j, i))) then; msg = '[[sublattice]] lists a species twice'; return; end if
        end do
        nullify (ka)
        call get_value(sp, 'concentrations', ka, requested=.false.)
        if (associated(ka)) then
          if (len(ka) /= nsc) then; msg = '[[sublattice]] concentrations must match species'; return; end if
          do j = 1, nsc; call get_value(ka, j, x); inp%ci%conc(j, i) = x; end do
        else
          if (nsc /= 1) then; msg = '[[sublattice]] with several species needs concentrations'; return; end if
          inp%ci%conc(1, i) = 1.0_dp
        end if
        if (any(inp%ci%conc(1:nsc, i) <= 0.0_dp) .or. abs(sum(inp%ci%conc(1:nsc, i)) - 1.0_dp) > 1.0e-8_dp) then
          msg = '[[sublattice]] concentrations must be > 0 and sum to 1'; return
        end if
      end do
      nmax = maxval(inp%ci%nsc)
      inp%ci%spc = inp%ci%spc(1:nmax, :); inp%ci%conc = inp%ci%conc(1:nmax, :)
      call decoration_build(lat, period, cen, pos, owner, nsub, inp%ci%dec, msg)
      if (len(msg) > 0) then; msg = '[[sublattice]]: '//msg; return; end if
      if (kspace .and. nmax > 1) then
        msg = 'method = "kspace" needs an ordered structure (one species per sublattice); disorder runs with method = "asr"'; return
      end if
    end if

    call section(table, 'scf', child)
    call check_keys(child, ['spin           ', 'xc             ', 'sim_beta       ', 'mixing         ', 'tolerance      ', &
                            'max_iterations ', 'zeta0          ', 'relativistic   ', 'start_potential', 'save_potential '], '[scf]', msg)
    if (len(msg) > 0) return
    if (kspace .and. (child%has_key('sim_beta') .or. child%has_key('start_potential') .or. child%has_key('save_potential'))) then
      msg = '[scf] sim_beta / start_potential / save_potential are read with method = "asr" only'; return
    end if
    call get_value(child, 'spin', inp%ci%nspin, 1)
    call get_value(child, 'xc', xcs, 'vbh'); inp%ci%xc = xcs
    call get_value(child, 'sim_beta', inp%ci%sim_beta, 1.0_dp)
    call get_value(child, 'mixing', inp%ci%beta, 0.3_dp)
    call get_value(child, 'tolerance', inp%ci%tol, merge(1.0e-6_dp, 1.0e-5_dp, kspace))
    call get_value(child, 'max_iterations', inp%ci%maxit, merge(100, 60, kspace))
    call get_value(child, 'zeta0', inp%ci%zeta0, merge(0.4_dp, 0.1_dp, kspace))
    call get_value(child, 'relativistic', inp%ci%relativistic, .true.)
    call get_value(child, 'start_potential', xcs, ''); inp%ci%pot_in = xcs
    call get_value(child, 'save_potential', xcs, ''); inp%ci%pot_out = xcs

    if (kspace .and. table%has_key('asr')) then; msg = '[asr] is read with method = "asr" only'; return; end if
    if (.not. kspace .and. table%has_key('kspace')) then; msg = '[kspace] is read with method = "kspace" only'; return; end if
    call section(table, 'asr', child)
    call check_keys(child, ['levels          ', 'recursion_steps ', 'rank            ', 'supercell       ', 'flip_radius     ', &
                            'contour_points  ', 'cluster_radius  ', 'inversion_radius', 'ed_terminator   ', 'fermi_kt        ', &
                            'fermi_poles     ', 'depth_tol       '], '[asr]', msg)
    if (len(msg) > 0) return
    call get_value(child, 'levels', inp%ci%nlevels, 4)
    call get_value(child, 'recursion_steps', inp%ci%nrec, 0)
    call get_value(child, 'rank', inp%ci%maxrank, 1)
    call get_value(child, 'supercell', inp%ci%lsuper, 8)
    call get_value(child, 'flip_radius', inp%ci%flip_radius, 99)
    call get_value(child, 'contour_points', inp%ci%ncont, 32)
    call get_value(child, 'cluster_radius', inp%ci%rcl, inp%ci%rcl)
    call get_value(child, 'inversion_radius', inp%ci%rinv, inp%ci%rinv)
    call get_value(child, 'ed_terminator', edt, 'bp'); inp%ci%ed_tail = edt
    call get_value(child, 'fermi_kt', inp%ci%fermi_kt, 0.001_dp)
    call get_value(child, 'fermi_poles', inp%ci%fermi_poles, 0)
    call get_value(child, 'depth_tol', inp%ci%depth_tol, inp%ci%depth_tol)

    call section(table, 'kspace', child)
    call check_keys(child, ['nk          ', 'bz          ', 'smear       ', 'sigma       ', 'trace       ', 'dos_nk      ', &
                            'sphere_ratio'], '[kspace]', msg)
    if (len(msg) > 0) return
    call get_value(child, 'nk', inp%li%nk, 12)
    call get_value(child, 'bz', str, 'tetra'); inp%li%bz = str
    call get_value(child, 'smear', inp%li%smear, 0.01_dp)
    call get_value(child, 'sigma', inp%li%sigma, 0.7_dp)
    call get_value(child, 'trace', inp%li%trace, .false.)
    call get_value(child, 'dos_nk', inp%li%dos_nk, 24)
    call get_value(child, 'sphere_ratio', ratio, 1.0_dp)

    call section(table, 'dos', child)
    call check_keys(child, ['emin   ', 'emax   ', 'npoints', 'eta    '], '[dos]', msg); if (len(msg) > 0) return
    call get_value(child, 'emin', inp%emin, -1.0_dp)
    call get_value(child, 'emax', inp%emax, 1.0_dp)
    call get_value(child, 'npoints', inp%npoints, merge(801, 401, kspace))
    call get_value(child, 'eta', inp%eta, merge(0.01_dp, 0.02_dp, kspace))

    call section(table, 'output', child)
    call check_keys(child, ['dos_file   ', 'out_file   ', 'params_file'], '[output]', msg); if (len(msg) > 0) return
    if (.not. kspace .and. child%has_key('params_file')) then; msg = '[output] params_file is written by method = "kspace" only'; return; end if
    call get_value(child, 'dos_file', inp%dos_file, 'dos.dat')
    call get_value(child, 'out_file', inp%out_file, 'asr.out')
    call get_value(child, 'params_file', inp%li%params_file, 'params.dat')

    ! ---- validation
    if (inp%ci%nspin /= 1 .and. inp%ci%nspin /= 2) msg = '[scf] spin must be 1 or 2'
    if (inp%ci%nlevels < 2) msg = '[asr] levels must be >= 2'
    if (inp%ci%nrec /= 0 .and. inp%ci%nrec < inp%ci%nlevels) msg = '[asr] recursion_steps must be 0 or >= levels'
    if (inp%ci%maxrank < 0) msg = '[asr] rank must be >= 0'
    if (.not. cell .and. modulo(inp%ci%lsuper, period) /= 0) msg = '[asr] supercell must be a multiple of [structure] period'
    if (inp%ci%fermi_poles < 0 .or. inp%ci%fermi_poles > 200) msg = '[asr] fermi_poles must be 0 (automatic) .. 200'
    if (edt /= 'hull' .and. edt /= 'bp' .and. edt /= 'avg2') msg = "[asr] ed_terminator must be 'hull', 'bp' or 'avg2'"
    if (inp%li%bz /= 'tetra' .and. inp%li%bz /= 'gauss') msg = "[kspace] bz must be 'tetra' or 'gauss'"
    if (inp%li%nk < 4) msg = '[kspace] nk must be >= 4'
    if (ratio <= 0.0_dp) msg = '[kspace] sphere_ratio must be > 0'
    if (inp%emax <= inp%emin) msg = '[dos] emax must exceed emin'
    if (len(msg) > 0) return

    if (kspace) then
      ! the k-space input of the same structure
      inp%li%name = inp%name; inp%li%infile = inp%infile; inp%li%workdir = inp%workdir
      inp%li%dos_file = inp%dos_file; inp%li%out_file = inp%out_file
      inp%li%alat_angstrom = alat; inp%li%xc = trim(inp%ci%xc); inp%li%spin = inp%ci%nspin
      inp%li%mixing = inp%ci%beta; inp%li%tolerance = inp%ci%tol; inp%li%maxit = inp%ci%maxit; inp%li%zeta0 = inp%ci%zeta0
      inp%li%relativistic = inp%ci%relativistic
      inp%li%emin = inp%emin; inp%li%emax = inp%emax; inp%li%npoints = inp%npoints; inp%li%dos_smear = inp%eta
      if (cell) then
        if (abs(ratio - 1.0_dp) > 1.0e-12_dp) then; msg = '[kspace] sphere_ratio needs a decorated B2 structure'; return; end if
        inp%li%lattice = 'cell'; inp%li%species = inp%species
      else
        call kspace_cell(inp%ci%dec, inp%species, inp%ci%spc(1, :), inp%subname, ratio, inp%li, msg)
      end if
    end if
  end subroutine read_run_input

  !> The ordered decorated structure as a k-space cell (li%lattice, cell_vec, site_pos, site_sp, species): the decoration cell is
  !> reduced to its primitive cell by the centering translations (I or F) that map every sublattice onto itself; one ASA class per
  !> sublattice (a species on several sublattices is labelled species@sublattice).  sphere_ratio /= 1: unequal spheres, B2 only
  !> (two sublattices, one site each, in the cubic cell; the sublattice at the origin is A).
  subroutine kspace_cell(dec, species, spc1, subname, ratio, li, msg)
    type(decoration_t), intent(in) :: dec
    type(lmto_species_t), intent(in) :: species(:)
    integer, intent(in) :: spc1(:)
    character(len=*), intent(in) :: subname(:)
    real(dp), intent(in) :: ratio
    type(lmto_input_t), intent(inout) :: li
    character(len=:), allocatable, intent(inout) :: msg
    integer :: tr(3, 4), ntr, i, j, t, k, p, nkeep, h(3)
    integer, allocatable :: keep(:)
    logical :: dup
    p = dec%period
    ! centering translations (half-lattice units) that are symmetries of the decoration
    tr(:, 1) = 0; ntr = 1
    if (invariant([0, p, p]) .and. invariant([p, 0, p]) .and. invariant([p, p, 0])) then
      tr(:, 2) = [0, p, p]; tr(:, 3) = [p, 0, p]; tr(:, 4) = [p, p, 0]; ntr = 4
      li%cell_vec = 0.5_dp * real(p, dp) * reshape([0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp], [3, 3])
    else if (invariant([p, p, p])) then
      tr(:, 2) = [p, p, p]; ntr = 2
      li%cell_vec = 0.5_dp * real(p, dp) * reshape([-1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp], [3, 3])
    else
      li%cell_vec = 0.0_dp
      do k = 1, 3; li%cell_vec(k, k) = real(p, dp); end do
    end if
    allocate (keep(dec%nsite)); nkeep = 0
    do i = 1, dec%nsite
      dup = .false.
      do j = 1, nkeep
        do t = 1, ntr
          if (all(modulo(dec%pts(:, i) - dec%pts(:, keep(j)) - tr(:, t), dec%m) == 0)) dup = .true.
        end do
      end do
      if (.not. dup) then; nkeep = nkeep + 1; keep(nkeep) = i; end if
    end do
    li%site_pos = 0.5_dp * real(dec%pts(:, keep(1:nkeep)), dp)
    li%site_sp = dec%pcls(keep(1:nkeep))
    li%site_classes = .false.
    allocate (li%species(dec%ncls))
    do k = 1, dec%ncls
      li%species(k) = species(spc1(k))
      if (count(spc1 == spc1(k)) > 1) li%species(k)%symbol = trim(species(spc1(k))%symbol)//'@'//trim(subname(k))
    end do
    li%lattice = 'cell'
    if (abs(ratio - 1.0_dp) > 1.0e-12_dp) then
      if (dec%ncls /= 2 .or. nkeep /= 2 .or. ntr /= 1 .or. p /= 1) then
        msg = '[kspace] sphere_ratio is implemented for B2 (two sublattices on the simple-cubic cell of the bcc graph) only'; return
      end if
      ! asr_lmto_driver's B2: species 1 at the origin, species 2 at the cube centre, r_1 / r_2 = ratio
      li%lattice = 'b2'; li%ratio = ratio
      li%species = [li%species(li%site_sp(1)), li%species(li%site_sp(2))]
    end if
  contains
    logical function invariant(v)
      integer, intent(in) :: v(3)
      integer :: n
      invariant = .true.
      do n = 1, dec%nsite
        h = modulo(dec%pts(:, n) + v, dec%m)
        if (dec%tab(h(1), h(2), h(3)) /= dec%pcls(n)) invariant = .false.
      end do
    end function invariant
  end subroutine kspace_cell

  subroutine run_input(inp, ok)
    type(run_input_t), intent(inout) :: inp
    logical, intent(out) :: ok
    if (inp%method == 'kspace') then
      call run_lmto(inp%li, ok)
    else
      call run_asr(inp, ok)
    end if
  end subroutine run_input

  subroutine run_asr(inp, ok)
    type(run_input_t), intent(inout) :: inp
    logical, intent(out) :: ok
    type(class_scf_result_t) :: res
    integer :: i, ie, l, p, sg, u, nsp, ncomp, ru, nthr, kap
    character(len=:), allocatable :: dfile
    character(len=8), allocatable :: labels(:)
    real(dp), allocatable :: zc(:), swsc(:), tsp(:)
    integer, allocatable :: nvc(:), natc(:)
    real(dp) :: t0, t_scf, avw, tot, qcell, mcell, fsp

    ok = .false.
    t0 = rep_wall()
    nsp = size(inp%species)
    allocate (inp%ci%spec(nsp))
    do i = 1, nsp
      inp%ci%spec(i)%z = inp%species(i)%z; inp%ci%spec(i)%nval = inp%species(i)%nval; inp%ci%spec(i)%spin_sign = inp%species(i)%spin_sign
      inp%ci%spec(i)%nvl = inp%species(i)%konfig(0:2)
      inp%ci%spec(i)%core = core_shells_for(inp%species(i)%konfig)
      inp%ci%spec(i)%atom_shells = atom_shells(inp%species(i))
    end do
    print '(a,a)', 'asr (augmented-space recursion)  ', inp%name
    print '(a,a,a,f9.5,a,i0,a,i0,a,i0)', '  graph lattice ', trim(inp%ci%lattice), '  a(Bohr)=', inp%ci%alat, '  period=', inp%ci%dec%period, &
      '  sites per decoration cell=', inp%ci%dec%nsite, '  spin=', inp%ci%nspin
    do kap = 1, inp%ci%dec%ncls
      print '(a,a8,a,f7.4,a,*(a))', '  sublattice ', inp%subname(kap), '  site fraction=', inp%ci%dec%frac(kap), '  species: ', &
        (trim(inp%species(inp%ci%spc(i, kap))%symbol)//' ', i=1, inp%ci%nsc(kap))
    end do
    inp%ci%ndos = inp%npoints; inp%ci%dos_emin = inp%emin; inp%ci%dos_emax = inp%emax; inp%ci%dos_eta = inp%eta
    t_scf = rep_wall()
    call class_scf(inp%ci, res, verbose=.true.)
    t_scf = rep_wall() - t_scf
    ncomp = res%ncomp
    allocate (labels(ncomp), zc(ncomp), swsc(ncomp), nvc(ncomp), natc(ncomp), tsp(res%nspin))
    avw = avg_ws(trim(inp%ci%lattice)) * inp%ci%alat
    qcell = 0.0_dp; mcell = 0.0_dp
    do p = 1, ncomp
      labels(p) = trim(inp%species(res%comp_spec(p))%symbol)//'@'//trim(inp%subname(res%comp_cls(p)))
      zc(p) = inp%species(res%comp_spec(p))%z; nvc(p) = inp%species(res%comp_spec(p))%nval; swsc(p) = avw; natc(p) = 1
      qcell = qcell + res%comp_w(p) * res%comp(p)%qnet * inp%ci%dec%nsite
      mcell = mcell + res%comp_w(p) * res%comp(p)%moment * inp%ci%dec%nsite
    end do
    print '(a,l1,a,i0,a,f10.6,a)', '  SCF converged=', res%converged, '  iterations=', res%iter, '  E_F=', res%ef, ' Ry'
    if (.not. res%converged) print '(a)', '  WARNING: SCF did not converge'
    do p = 1, ncomp
      print '(a,a12,a,f7.4,a,f9.5,a,f9.5)', '  component ', labels(p), '  x=', res%comp_x(p), '  qnet=', res%comp(p)%qnet, '  moment=', res%comp(p)%moment
    end do
    print '(a,i0,a,f9.5,a,es10.2)', '  per decoration cell (', inp%ci%dec%nsite, ' sites): moment=', mcell, '  net charge=', qcell

    dfile = inp%workdir//inp%dos_file
    open (newunit=u, file=dfile, status='replace', action='write')
    write (u, '(a)') '# '//inp%name//'   states/Ry per atom of each component, summed over spin; total = per site (weighted)'
    write (u, '(a,f12.8)') '# E_F = ', res%ef
    write (u, '(a,*(a))') '# E(Ry)  total  ', ('  '//trim(labels(p))//'_s _p _d', p=1, ncomp), '  total_spin1 total_spin2'
    fsp = 2.0_dp / real(res%nspin, dp)                 ! res%edos is per spin channel; spin = 1 holds both spins in one
    do ie = 1, inp%npoints
      tsp = 0.0_dp
      do sg = 1, res%nspin
        tsp(sg) = fsp * sum([(res%comp_w(p) * sum(res%edos(ie, :, p, sg)), p=1, ncomp)])
      end do
      tot = sum(tsp)
      write (u, '(f12.6,*(es14.6))') res%edos_e(ie), tot, ((fsp * sum(res%edos(ie, l, p, :)), l=0, 2), p=1, ncomp), tsp
    end do
    close (u)
    print '(a,a)', '  wrote ', dfile

    nthr = 1
    !$ nthr = omp_get_max_threads()
    call rep_open(ru, inp%workdir//inp%out_file, 'Sublattice alloy, augmented-space recursion + TB-LMTO-ASA: '//inp%name, inp%infile, nthr)
    call rep_echo_input(ru, inp%infile)
    call rep_section(ru, 'method')
    call rep_line(ru, 'Hamiltonian   TB-LMTO-ASA, energy-dependent ASR (asr_ed form), one potential per (sublattice, species) component')
    write (ru, '(a,i0,a,i0,a,i0,a,i0)') '  levels=', inp%ci%nlevels, '  recursion steps=', merge(inp%ci%nrec, inp%ci%nlevels, inp%ci%nrec > 0), &
      '  max rank=', inp%ci%maxrank, '  supercell=', inp%ci%lsuper
    write (ru, '(a,*(i0,1x))') '  augmented-space nodes per sublattice: ', res%nodes
    write (ru, '(a,f6.3,a,f6.3,a,f7.4,a)') '  SIM beta=', inp%ci%sim_beta, '   mixing=', inp%ci%beta, '   fermi_kt=', inp%ci%fermi_kt, ' Ry'
    call rep_line(ru, 'electrostatics  Madelung of the sublattice-average charges (Ewald, decoration cell) + SIM of the deviation')
    call rep_section(ru, 'structure')
    write (ru, '(a,a,a,f12.6,a,f10.6,a)') ' graph lattice ', trim(inp%ci%lattice), '  a = ', inp%ci%alat, ' Bohr = ', inp%ci%alat * 0.529177210903_dp, ' Angstrom'
    write (ru, '(a,i0,a,i0,a,f9.5,a)') ' decoration cell ', inp%ci%dec%period, '^3 conventional cubes, ', inp%ci%dec%nsite, ' sites;  S_ws = ', avw, ' Bohr (equal spheres)'
    do kap = 1, inp%ci%dec%ncls
      write (ru, '(a,a8,a,f7.4,a,*(a))') ' sublattice ', inp%subname(kap), '  site fraction ', inp%ci%dec%frac(kap), '   ', &
        (trim(inp%species(inp%ci%spc(i, kap))%symbol)//' ', i=1, inp%ci%nsc(kap))
      write (ru, '(a,*(f8.4))') '      concentrations ', inp%ci%conc(1:inp%ci%nsc(kap), kap)
    end do
    call rep_scf_history(ru, labels, res%iter, res%hist_dv, res%hist_ef, res%hist_moment, res%hist_qnet, inp%ci%tol, res%converged)
    block
      type(asa_class_result_t), allocatable :: rc(:)
      rc = res%comp
      if (res%nspin == 1) then
        do p = 1, ncomp; rc(p)%q = 2.0_dp * rc(p)%q; end do
      end if
      call rep_classes(ru, labels, zc, nvc, swsc, natc, rc, res%nspin, res%ef, mcell, ' (per component atom; cell = decoration cell)')
    end block
    write (ru, '(a)') ' components are (species@sublattice); qnet/moments per atom of that component; cell moment per decoration cell'
    call rep_dos_at_ef(ru, labels, res%edos_e, res%edos, res%ef, inp%eta, res%nspin, 'per-component LDOS from the recursion (eta broadening)')
    call rep_section(ru, 'timing')
    write (ru, '(a,f10.2,a)') ' SCF + DOS  ', t_scf, ' s'
    write (ru, '(a)') ' files written:'
    write (ru, '(a)') '   '//dfile
    call rep_close(ru, res%converged, rep_wall() - t0)
    print '(a,a)', '  wrote ', inp%workdir//inp%out_file
    ok = res%converged
  contains
    real(dp) function avg_ws(kind)
      use asr_strconst, only: avg_ws_radius
      character(len=*), intent(in) :: kind
      avg_ws = avg_ws_radius(kind)
    end function avg_ws
  end subroutine run_asr
end module asr_driver
