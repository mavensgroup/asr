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
!>   [[sro]]        (method = "asr") first-shell Warren-Cowley short-range order: sublattices = ["A", "B"] (optional with one
!>                  sublattice), shell = 1, pairs = [["Fe", "Ni", alpha], ...] (one sublattice: unlike pairs, symmetric, diagonal
!>                  derived; alpha > 0 = the pair avoids contact, 1 = never neighbours); see asr_sro
!>   [scf]          spin, xc, mixing, tolerance, max_iterations, zeta0, relativistic; method = "asr" also: sim_beta,
!>                  save_potential (written every iteration) / start_potential, mixing_restart (linear step after an Anderson
!>                  restart; default = mixing; smaller for stiff charge-transfer maps)
!>                  sim_beta (0 <= beta < lattice maximum, about 1.39 on bcc) fixes the hardness of the charge-excess screening
!>                  model; with [[sro]] it also sets the screened potential of the correlated first shell (asr_screening)
!>   [asr]          (method = "asr") levels, recursion_steps, rank, supercell (multiple of period), flip_radius, contour_points,
!>                  cluster_radius, inversion_radius, ed_terminator (hull (default) | bp | avg2; bp can put a pole of the
!>                  terminated fraction next to the contour on complex chains, see docs/FINDINGS_LOG.md), fermi_kt, fermi_poles, depth_tol,
!>                  ef_scan, ef_scan_step, ef_scan_center, ef_scan_dump (diagnostic: N(E) at the start potential, then stop),
!>                  route = "edc" (default: H^(2), asr_edc) | "ed3" (full generalised ASA problem, docs/DERIVATION_ed3.md),
!>                  alpha_scale (diagnostic: screening constants times this factor; default 1),
!>                  energy_block (contour energies per chain block, default 8, capped by memory; 0 = one energy per chain with the
!>                  older real kernels, for checks), ef_linear_max (Ry; E_F steps up to this size are taken to first order without
!>                  a new contour, exact at self-consistency; default 0 = off)
!>   [kspace]       (method = "kspace") nk, bz = "tetra"|"gauss", smear, sigma, trace, dos_nk, sphere_ratio (r_A / r_B of a B2
!>                  structure; default 1 = equal spheres), alpha_scale (diagnostic: screening constants times this factor),
!>                  std_alpha (diagnostic: the recursion's standard alpha instead of the hard-core alpha of sigma),
!>                  second_order (diagnostic, default false: the orthogonal second-order
!>                  H^(2) = E_nu + h - h o h of the recursion route, no h p h overlap, density M0 phi^2 + 2 M1 phi phidot +
!>                  M2 (phidot^2 + phi phiddot), so the sphere charge is M0 as in the recursion)
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
  use asr_sro, only: sro_conditional
  use asr_screening, only: screen_t, screen_lattice
  use asr_class_scf
  use asr_report
  use asr_mpi, only: par_root, par_size
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
                            'kspace     ', 'dos        ', 'output     ', 'sro        '], 'top level', msg)
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
    if (table%has_key('sro')) then
      if (kspace .or. cell) then; msg = '[[sro]] needs method = "asr" on a decorated lattice'; return; end if
      call read_sro(table, inp, msg); if (len(msg) > 0) return
    end if

    call section(table, 'scf', child)
    call check_keys(child, ['spin           ', 'xc             ', 'sim_beta       ', 'mixing         ', 'tolerance      ', &
                            'max_iterations ', 'zeta0          ', 'relativistic   ', 'start_potential', 'save_potential ', &
                            'mixing_restart '], '[scf]', msg)
    if (len(msg) > 0) return
    if (kspace .and. (child%has_key('sim_beta') .or. child%has_key('start_potential') .or. child%has_key('save_potential') &
                      .or. child%has_key('mixing_restart'))) then
      msg = '[scf] sim_beta / start_potential / save_potential / mixing_restart are read with method = "asr" only'; return
    end if
    call get_value(child, 'spin', inp%ci%nspin, 1)
    call get_value(child, 'xc', xcs, 'vbh'); inp%ci%xc = xcs
    call get_value(child, 'sim_beta', inp%ci%sim_beta, 1.0_dp)
    if (.not. kspace) then                         ! sim_beta fixes the hardness of the screening model; it has a lattice maximum
      block
        type(screen_t) :: scr
        call screen_lattice(trim(inp%ci%lattice), inp%ci%alat, inp%ci%sim_beta, scr, msg)
        if (len(msg) > 0) then; msg = '[scf] sim_beta: '//msg; return; end if
      end block
    end if
    call get_value(child, 'mixing', inp%ci%beta, 0.3_dp)
    call get_value(child, 'mixing_restart', inp%ci%beta_restart, 0.0_dp)
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
                            'fermi_poles     ', 'depth_tol       ', 'ef_scan         ', 'ef_scan_step    ', 'ef_scan_center  ', &
                            'ef_scan_dump    ', 'route           ', 'alpha_scale     ', 'energy_block    ', 'ef_linear_max   '], &
                   '[asr]', msg)
    if (len(msg) > 0) return
    call get_value(child, 'levels', inp%ci%nlevels, 4)
    call get_value(child, 'recursion_steps', inp%ci%nrec, 0)
    call get_value(child, 'rank', inp%ci%maxrank, 1)
    call get_value(child, 'supercell', inp%ci%lsuper, 8)
    call get_value(child, 'flip_radius', inp%ci%flip_radius, 99)
    call get_value(child, 'contour_points', inp%ci%ncont, 32)
    call get_value(child, 'cluster_radius', inp%ci%rcl, inp%ci%rcl)
    call get_value(child, 'inversion_radius', inp%ci%rinv, inp%ci%rinv)
    call get_value(child, 'ed_terminator', edt, 'hull'); inp%ci%ed_tail = edt
    call get_value(child, 'fermi_kt', inp%ci%fermi_kt, 0.001_dp)
    call get_value(child, 'fermi_poles', inp%ci%fermi_poles, 0)
    call get_value(child, 'depth_tol', inp%ci%depth_tol, inp%ci%depth_tol)
    call get_value(child, 'ef_scan', inp%ci%ef_scan, 0)
    call get_value(child, 'ef_scan_step', inp%ci%ef_scan_step, inp%ci%ef_scan_step)
    call get_value(child, 'ef_scan_center', inp%ci%ef_scan_center, inp%ci%ef_scan_center)
    call get_value(child, 'ef_scan_dump', xcs, ''); inp%ci%ef_scan_dump = xcs
    call get_value(child, 'route', xcs, 'edc'); inp%ci%route = xcs
    if (xcs /= 'edc' .and. xcs /= 'ed3') then; msg = "[asr] route must be 'edc' or 'ed3'"; return; end if
    call get_value(child, 'alpha_scale', inp%ci%alpha_scale, 1.0_dp)
    if (inp%ci%alpha_scale <= 0.0_dp) then; msg = '[asr] alpha_scale must be positive'; return; end if
    call get_value(child, 'energy_block', inp%ci%eblock, 8)
    if (inp%ci%eblock < 0) then; msg = '[asr] energy_block must be >= 0'; return; end if
    call get_value(child, 'ef_linear_max', inp%ci%ef_linear_max, 0.0_dp)
    if (inp%ci%ef_linear_max < 0.0_dp) then; msg = '[asr] ef_linear_max must be >= 0'; return; end if

    call section(table, 'kspace', child)
    call check_keys(child, ['nk          ', 'bz          ', 'smear       ', 'sigma       ', 'trace       ', 'dos_nk      ', &
                            'sphere_ratio', 'second_order', 'alpha_scale ', 'std_alpha   '], '[kspace]', msg)
    if (len(msg) > 0) return
    call get_value(child, 'nk', inp%li%nk, 12)
    call get_value(child, 'bz', str, 'tetra'); inp%li%bz = str
    call get_value(child, 'smear', inp%li%smear, 0.01_dp)
    call get_value(child, 'sigma', inp%li%sigma, 0.7_dp)
    call get_value(child, 'trace', inp%li%trace, .false.)
    call get_value(child, 'dos_nk', inp%li%dos_nk, 24)
    call get_value(child, 'second_order', inp%li%second_order, .false.)
    call get_value(child, 'alpha_scale', inp%li%alpha_scale, 1.0_dp)
    call get_value(child, 'std_alpha', inp%li%std_alpha, .false.)
    if (inp%li%alpha_scale <= 0.0_dp) then; msg = '[kspace] alpha_scale must be positive'; return; end if
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

  !> [[sro]] entries: sublattices = ["A", "B"] (optional with one sublattice), shell = 1, pairs = [["Fe", "Ni", alpha], ...].
  !> One sublattice: symmetric off-diagonal Warren-Cowley parameters, the diagonal is derived.  Two sublattices: pair [a, b, alpha]
  !> is species a on the first, b on the second; unset pairs are 0 and both marginals must vanish (asr_sro).  Fills
  !> inp%ci%sro_cond for every sublattice pair (alpha = 0 where no entry).
  subroutine read_sro(table, inp, msg)
    type(toml_table), intent(inout) :: table
    type(run_input_t), intent(inout) :: inp
    character(len=:), allocatable, intent(out) :: msg
    type(toml_array), pointer :: arr, ka, kb
    type(toml_table), pointer :: sp
    real(dp), allocatable :: alpha(:, :, :, :), cond(:, :)
    logical, allocatable :: given(:, :), seen(:, :, :, :)
    character(len=:), allocatable :: str
    integer :: nsub, nmax, i, j, k, ish, ks(2), kk(2), n1, n2
    real(dp) :: x
    msg = ''
    nsub = size(inp%ci%nsc); nmax = maxval(inp%ci%nsc)
    allocate (alpha(nmax, nmax, nsub, nsub), given(nsub, nsub), seen(nmax, nmax, nsub, nsub))
    alpha = 0.0_dp; given = .false.; seen = .false.
    call get_value(table, 'sro', arr, requested=.false.)
    if (.not. associated(arr)) then; msg = '[[sro]] must be an array of tables'; return; end if
    do i = 1, len(arr)
      call get_value(arr, i, sp)
      if (.not. associated(sp)) then; msg = '[[sro]] must be an array of tables'; return; end if
      call check_keys(sp, ['sublattices', 'shell      ', 'pairs      '], '[[sro]]', msg); if (len(msg) > 0) return
      call get_value(sp, 'shell', ish, 1)
      if (ish /= 1) then; msg = '[[sro]] shell: only the first shell (1) is implemented'; return; end if
      nullify (ka); call get_value(sp, 'sublattices', ka, requested=.false.)
      if (associated(ka)) then
        if (len(ka) /= 2) then; msg = '[[sro]] sublattices = ["A", "B"] (two names, may be equal)'; return; end if
        do j = 1, 2
          call get_value(ka, j, str); ks(j) = 0
          do k = 1, nsub
            if (trim(inp%subname(k)) == str) ks(j) = k
          end do
          if (ks(j) == 0) then; msg = '[[sro]] sublattices: '//str//' is not a [[sublattice]] name'; return; end if
        end do
      else
        if (nsub /= 1) then; msg = '[[sro]] needs sublattices = ["A", "B"] when there are several sublattices'; return; end if
        ks = 1
      end if
      if (given(ks(1), ks(2)) .or. given(ks(2), ks(1))) then; msg = '[[sro]] sublattice pair given twice'; return; end if
      given(ks(1), ks(2)) = .true.
      nullify (ka); call get_value(sp, 'pairs', ka, requested=.false.)
      if (.not. associated(ka)) then; msg = '[[sro]] needs pairs = [["A", "B", alpha], ...]'; return; end if
      do j = 1, len(ka)
        nullify (kb); call get_value(ka, j, kb)
        if (.not. associated(kb)) then; msg = '[[sro]] pairs must be [["A", "B", alpha], ...]'; return; end if
        if (len(kb) /= 3) then; msg = '[[sro]] pairs must be [["A", "B", alpha], ...]'; return; end if
        do k = 1, 2
          call get_value(kb, k, str); kk(k) = 0
          if (.not. allocated(str)) then; msg = '[[sro]] pairs must be [["A", "B", alpha], ...]'; return; end if
          do n1 = 1, inp%ci%nsc(ks(k))
            if (inp%species(inp%ci%spc(n1, ks(k)))%symbol == str) kk(k) = n1
          end do
          if (kk(k) == 0) then; msg = '[[sro]] pair species '//str//' is not on sublattice '//trim(inp%subname(ks(k))); return; end if
        end do
        call get_value(kb, 3, x)
        if (ks(1) == ks(2) .and. kk(1) == kk(2)) then
          msg = '[[sro]] alpha of a species with itself is derived (one sublattice); give the unlike pairs only'; return
        end if
        if (seen(kk(1), kk(2), ks(1), ks(2))) then; msg = '[[sro]] pair given twice'; return; end if
        alpha(kk(1), kk(2), ks(1), ks(2)) = x; alpha(kk(2), kk(1), ks(2), ks(1)) = x
        seen(kk(1), kk(2), ks(1), ks(2)) = .true.; seen(kk(2), kk(1), ks(2), ks(1)) = .true.
      end do
    end do
    inp%ci%sro = .true.
    allocate (inp%ci%sro_cond(nmax, nmax, nsub, nsub)); inp%ci%sro_cond = 0.0_dp
    do i = 1, nsub
      do j = 1, nsub
        n1 = inp%ci%nsc(i); n2 = inp%ci%nsc(j)
        call sro_conditional(inp%ci%conc(1:n1, i), inp%ci%conc(1:n2, j), alpha(1:n1, 1:n2, i, j), i == j, cond, msg)
        if (len(msg) > 0) then
          msg = '[[sro]] '//trim(inp%subname(i))//' - '//trim(inp%subname(j))//': '//msg; return
        end if
        inp%ci%sro_cond(1:n1, 1:n2, i, j) = cond
      end do
    end do
  end subroutine read_sro

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
      ! k space has no MPI split: rank 0 runs it, the other ranks wait at the end
      ok = .true.
      if (par_size() > 1) print '(a,i0,a)', '  kspace: no MPI split, rank 0 of ', par_size(), ' runs the calculation'
      if (par_root()) call run_lmto(inp%li, ok)
    else
      call run_asr(inp, ok)
    end if
  end subroutine run_input

  subroutine run_asr(inp, ok)
    type(run_input_t), intent(inout) :: inp
    logical, intent(out) :: ok
    type(class_scf_result_t) :: res
    integer :: i, ie, l, p, sg, u, nsp, ncomp, ru, nthr, kap, k
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
    print '(a,f18.8,a)', '  total energy =', res%etot, ' Ry per site'
    if (.not. res%converged) print '(a)', '  WARNING: SCF did not converge'
    do p = 1, ncomp
      print '(a,a12,a,f7.4,a,f9.5,a,f9.5)', '  component ', labels(p), '  x=', res%comp_x(p), '  qnet=', res%comp(p)%qnet, '  moment=', res%comp(p)%moment
    end do
    print '(a,i0,a,f9.5,a,es10.2)', '  per decoration cell (', inp%ci%dec%nsite, ' sites): moment=', mcell, '  net charge=', qcell

    ! every rank holds the full result (asr_class_scf); rank 0 alone writes the files
    if (.not. par_root()) then
      ok = res%converged; return
    end if
    dfile = inp%workdir//inp%dos_file
    open (newunit=u, file=dfile, status='replace', action='write')
    write (u, '(a)') '# '//inp%name//'   states/Ry per atom of each component, summed over spin; total = per site (weighted)'
    write (u, '(a,f12.8)') '# E_F = ', res%ef
    if (res%nspin == 2) then                           ! spin-resolved columns per component follow (per atom, per spin)
      write (u, '(a,*(a))') '# E(Ry)  total  ', ('  '//trim(labels(p))//'_s _p _d', p=1, ncomp), '  total_spin1 total_spin2', &
        ('  '//trim(labels(p))//'_up_s _p _d _dn_s _p _d', p=1, ncomp)
    else
      write (u, '(a,*(a))') '# E(Ry)  total  ', ('  '//trim(labels(p))//'_s _p _d', p=1, ncomp), '  total_spin1 total_spin2'
    end if
    fsp = 2.0_dp / real(res%nspin, dp)                 ! res%edos is per spin channel; spin = 1 holds both spins in one
    do ie = 1, inp%npoints
      tsp = 0.0_dp
      do sg = 1, res%nspin
        tsp(sg) = fsp * sum([(res%comp_w(p) * sum(res%edos(ie, :, p, sg)), p=1, ncomp)])
      end do
      tot = sum(tsp)
      if (res%nspin == 2) then
        write (u, '(f12.6,*(es14.6))') res%edos_e(ie), tot, ((fsp * sum(res%edos(ie, l, p, :)), l=0, 2), p=1, ncomp), tsp, &
          (((res%edos(ie, l, p, sg), l=0, 2), sg=1, 2), p=1, ncomp)
      else
        write (u, '(f12.6,*(es14.6))') res%edos_e(ie), tot, ((fsp * sum(res%edos(ie, l, p, :)), l=0, 2), p=1, ncomp), tsp
      end if
    end do
    close (u)
    print '(a,a)', '  wrote ', dfile

    nthr = 1
    !$ nthr = omp_get_max_threads()
    call rep_open(ru, inp%workdir//inp%out_file, 'Sublattice alloy, augmented-space recursion + TB-LMTO-ASA: '//inp%name, inp%infile, nthr)
    if (par_size() > 1) write (ru, '(a,i0,a)') ' MPI ranks: ', par_size(), ' (the recursion chains of each energy pass are split over the ranks)'
    call rep_echo_input(ru, inp%infile)
    call rep_section(ru, 'method')
    if (trim(inp%ci%route) == 'ed3') then
      call rep_line(ru, 'Hamiltonian   TB-LMTO-ASA, full generalised ASA problem (route ed3, two chains per energy), one potential per (sublattice, species) component')
    else
      call rep_line(ru, 'Hamiltonian   TB-LMTO-ASA, energy-dependent ASR (asr_ed form), one potential per (sublattice, species) component')
    end if
    if (inp%ci%alpha_scale /= 1.0_dp) write (ru, '(a,f7.4)') '  screening constants scaled by ', inp%ci%alpha_scale
    write (ru, '(a,i0,a,i0,a,i0,a,i0)') '  levels=', inp%ci%nlevels, '  recursion steps=', merge(inp%ci%nrec, inp%ci%nlevels, inp%ci%nrec > 0), &
      '  max rank=', inp%ci%maxrank, '  supercell=', inp%ci%lsuper
    write (ru, '(a,*(i0,1x))') '  augmented-space nodes per sublattice: ', res%nodes
    write (ru, '(a,f6.3,a,f6.3,a,f7.4,a)') '  SIM beta=', inp%ci%sim_beta, '   mixing=', inp%ci%beta, '   fermi_kt=', inp%ci%fermi_kt, ' Ry'
    if (res%screen%unscreened) then
      write (ru, '(a)') '  screening: none (sim_beta = 0; fluctuations interact as bare point charges)'
    else
      write (ru, '(a,f8.4,a,f7.4,a,f8.4,a,f8.4,a,f7.4,a,f7.4)') '  screening (charge-excess model of the graph lattice): hardness a=', &
        res%screen%a, ' Ry/e  G00=', res%screen%g00, '  G1=', res%screen%g1, '  W00/G00=', res%screen%w00 / res%screen%g00, &
        '  first-shell weight lambda=', res%screen%lambda, '  beta max=', res%screen%beta_max
    end if
    call rep_line(ru, 'electrostatics  Madelung of the sublattice-average charges (Ewald, decoration cell) + screened deviations (SIM; with SRO'// &
                  ' the conditional first shell, asr_screening)')
    call rep_section(ru, 'structure')
    write (ru, '(a,a,a,f12.6,a,f10.6,a)') ' graph lattice ', trim(inp%ci%lattice), '  a = ', inp%ci%alat, ' Bohr = ', inp%ci%alat * 0.529177210903_dp, ' Angstrom'
    write (ru, '(a,i0,a,i0,a,f9.5,a)') ' decoration cell ', inp%ci%dec%period, '^3 conventional cubes, ', inp%ci%dec%nsite, ' sites;  S_ws = ', avw, ' Bohr (equal spheres)'
    do kap = 1, inp%ci%dec%ncls
      write (ru, '(a,a8,a,f7.4,a,*(a))') ' sublattice ', inp%subname(kap), '  site fraction ', inp%ci%dec%frac(kap), '   ', &
        (trim(inp%species(inp%ci%spc(i, kap))%symbol)//' ', i=1, inp%ci%nsc(kap))
      write (ru, '(a,*(f8.4))') '      concentrations ', inp%ci%conc(1:inp%ci%nsc(kap), kap)
    end do
    if (inp%ci%sro) then
      write (ru, '(a)') ' first-shell short-range order (Warren-Cowley): P(neighbour species | centre species); shell sites'// &
        ' independent given the centre, uncorrelated beyond the first shell; electrostatics: screened conditional shell (asr_screening)'
      do kap = 1, inp%ci%dec%ncls
        do k = 1, inp%ci%dec%ncls
          if (all(abs(inp%ci%sro_cond(1:inp%ci%nsc(kap), 1:inp%ci%nsc(k), kap, k) - &
                      spread(inp%ci%conc(1:inp%ci%nsc(k), k), 1, inp%ci%nsc(kap))) < 1.0e-12_dp)) cycle
          do i = 1, inp%ci%nsc(kap)
            write (ru, '(a,a8,a,a8,a,a4,a,*(f8.4))') '   centre on ', inp%subname(kap), ' neighbour on ', inp%subname(k), &
              ' centre ', trim(inp%species(inp%ci%spc(i, kap))%symbol), ': ', inp%ci%sro_cond(i, 1:inp%ci%nsc(k), kap, k)
          end do
        end do
      end do
    end if
    call rep_scf_history(ru, labels, res%iter, res%hist_dv, res%hist_ef, res%hist_moment, res%hist_qnet, inp%ci%tol, res%converged)
    write (ru, '(a)') ' total energy per site (Ry; Kohn-Sham form, output density in the input potential), by iteration:'
    write (ru, '(i6,f20.8)') (i, res%hist_etot(i), i=1, res%iter)
    write (ru, '(a,f20.8,a)') ' total energy =', res%etot, ' Ry per site'
    block
      type(asa_class_result_t), allocatable :: rc(:)
      rc = res%comp
      if (res%nspin == 1) then
        do p = 1, ncomp; rc(p)%q = 2.0_dp * rc(p)%q; end do
      end if
      call rep_classes(ru, labels, zc, nvc, swsc, natc, rc, res%nspin, res%ef, mcell, ' (per component atom; cell = decoration cell)', &
                       natw=res%comp_w * inp%ci%dec%nsite)
    end block
    if (trim(inp%ci%route) == 'ed3') then
      write (ru, '(a)') ' route ed3: sphere charge per channel, s p d, per spin channel (q0 + p q2; the table above lists q0 alone)'
      do p = 1, ncomp
        write (ru, '(3x,a12,*(3f10.5,2x))') labels(p), (2.0_dp / res%nspin * res%qch(p, sg, :), sg=1, res%nspin)
      end do
    end if
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
