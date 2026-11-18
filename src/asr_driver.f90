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
!>                  energy_block (contour energies per chain block, capped by memory; default 0 = one energy per chain with the
!>                  real kernels; 8 was 1.35x faster for CoCrFeNi at about 6x the memory; the 1e-4 block difference in B2 ed3 is
!>                  rounding scatter inside the length error of the chain, docs/REVIEW_fable.md 3.8), ef_linear_max (Ry; E_F
!>                  steps up to this size are taken to first order without a new contour, exact at self-consistency;
!>                  default 0 = off)
!>   [kspace]       (method = "kspace") nk, bz = "tetra"|"gauss", smear, sigma, trace, dos_nk, sphere_ratio (r_A / r_B of a B2
!>                  structure; default 1 = equal spheres), alpha_scale (diagnostic: screening constants times this factor),
!>                  std_alpha (diagnostic: the recursion's standard alpha instead of the hard-core alpha of sigma),
!>                  second_order (diagnostic, default false: the orthogonal second-order
!>                  H^(2) = E_nu + h - h o h of the recursion route, no h p h overlap, density M0 phi^2 + 2 M1 phi phidot +
!>                  M2 (phidot^2 + phi phiddot), so the sphere charge is M0 as in the recursion.  H^(2) depends on the
!>                  screening constants at third order: it is the recursion's H^(2) only with std_alpha = true)
!>   [dos]          emin, emax, npoints, eta (recursion: imaginary part; k-space: smearing)
!>   [output]       dos_file, out_file, params_file (method = "kspace")
!> Examples: pure bcc Fe = lattice "bcc", centering "I", one sublattice [[0,0,0]] with species ["Fe"]; B2 FeCo = "bcc", period 1,
!> centering "P", Fe (0,0,0) and Co (1/2,1/2,1/2); L2_1 Co2Fe0.4Cr0.6Al on the bcc graph (a_bcc = a/2): period = 2, centering = "F";
!> Al (0,0,0); Fe/Cr (1/2,1/2,1/2); Co (1/4,1/4,1/4) and (3/4,3/4,3/4).  Unknown sections/keys are errors.
module asr_driver
  use asr_kinds, only: dp
  use tomlf, only: toml_table, toml_array, toml_error, toml_load, get_value, set_value, len
  use asr_inputdump, only: resolved_input
  use asr_tomlutil, only: section, check_keys
  use asr_input_lmto, only: lmto_species_t, lmto_input_t
  use asr_asa, only: core_shells_for, asa_class_result_t
  use asr_lmto_driver, only: run_lmto, atom_shells
  use asr_sublattice, only: decoration_t, decoration_build
  use asr_sro, only: sro_conditional
  use asr_screening, only: screen_t, screen_lattice, mode_f
  use asr_class_scf
  use asr_report
  use asr_mpi, only: par_root, par_size
  !$ use omp_lib, only: omp_get_max_threads
  implicit none
  private
  public :: run_input_t, read_run_input, run_input, kspace_cell, class_species, rep_from_class_scf

  type :: run_input_t
    character(len=:), allocatable :: name, method, dos_file, out_file, workdir, infile
    character(len=:), allocatable :: spec_file          !! [spectral] file (method = "asr")
    character(len=:), allocatable :: resolved           !! every input key with the value used (asr_inputdump); for the report
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
                            'kspace     ', 'dos        ', 'output     ', 'sro        ', 'spectral   '], 'top level', msg)
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
      else                                             ! the default, written into the table for the resolved input
        call get_value(sp, 'konfig', ka); call set_value(ka, inp%species(i)%konfig)
      end if
      nullify (ka)
      call get_value(sp, 'occupations', ka, requested=.false.)
      if (associated(ka)) then
        if (len(ka) /= 2) then; msg = '[[species]] occupations needs 2 numbers (s,d)'; return; end if
        call get_value(ka, 1, x); inp%species(i)%occ_s = x
        call get_value(ka, 2, x); inp%species(i)%occ_d = x
      else
        call get_value(sp, 'occupations', ka); call set_value(ka, [inp%species(i)%occ_s, inp%species(i)%occ_d])
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
          call get_value(sp, 'concentrations', ka); call set_value(ka, [1.0_dp])
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
    call get_value(child, 'energy_block', inp%ci%eblock, 0)
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

    if (table%has_key('spectral')) then                ! Bloch spectral function after the SCF (asr_spectral)
      if (kspace) then; msg = '[spectral] is read with method = "asr" only'; return; end if
      call section(table, 'spectral', child)
      call check_keys(child, ['kpath  ', 'labels ', 'nk     ', 'emin   ', 'emax   ', 'npoints', 'eta    ', 'rmax   ', 'steps  ', &
                              'file   ', 'window '], '[spectral]', msg)
      if (len(msg) > 0) return
      if (trim(inp%ci%route) /= 'edc') then; msg = '[spectral] needs route = "edc"'; return; end if
      nullify (arr); call get_value(child, 'kpath', arr, requested=.false.)
      if (.not. associated(arr)) then; msg = '[spectral] needs kpath = [[kx, ky, kz], ...] (Cartesian, units 2 pi / a)'; return; end if
      if (len(arr) < 2) then; msg = '[spectral] kpath needs at least two points'; return; end if
      allocate (inp%ci%spectral%kpts(3, len(arr)), inp%ci%spectral%labels(len(arr))); inp%ci%spectral%labels = ''
      do i = 1, len(arr)
        nullify (ka); call get_value(arr, i, ka)
        if (.not. associated(ka)) then; msg = '[spectral] kpath: each point is an array of 3 numbers'; return; end if
        if (len(ka) /= 3) then; msg = '[spectral] kpath: each point is an array of 3 numbers'; return; end if
        do j = 1, 3; call get_value(ka, j, x); inp%ci%spectral%kpts(j, i) = x; end do
      end do
      nullify (arr); call get_value(child, 'labels', arr, requested=.false.)
      if (associated(arr)) then
        if (len(arr) /= size(inp%ci%spectral%kpts, 2)) then; msg = '[spectral] labels needs one label per kpath point'; return; end if
        do i = 1, len(arr); call get_value(arr, i, str); inp%ci%spectral%labels(i) = str; end do
      end if
      call get_value(child, 'nk', inp%ci%spectral%nk, 30)
      call get_value(child, 'emin', inp%ci%spectral%emin, inp%emin)
      call get_value(child, 'emax', inp%ci%spectral%emax, inp%emax)
      call get_value(child, 'npoints', inp%ci%spectral%ne, 201)
      call get_value(child, 'eta', inp%ci%spectral%eta, 0.01_dp)
      call get_value(child, 'rmax', inp%ci%spectral%rmax, 2.0_dp)
      call get_value(child, 'window', str, 'ball')
      if (str /= 'ball' .and. str /= 'sharp') then; msg = '[spectral] window must be "ball" or "sharp"'; return; end if
      inp%ci%spectral%window = str
      call get_value(child, 'steps', inp%ci%spectral%nsteps, 0)
      call get_value(child, 'file', inp%spec_file, 'spectral.dat')
      if (inp%ci%spectral%nk < 1) then; msg = '[spectral] nk must be >= 1'; return; end if
      if (inp%ci%spectral%ne < 2) then; msg = '[spectral] npoints must be >= 2'; return; end if
      if (inp%ci%spectral%eta <= 0.0_dp) then; msg = '[spectral] eta must be > 0'; return; end if
      if (inp%ci%spectral%rmax <= 0.0_dp) then; msg = '[spectral] rmax must be > 0'; return; end if
      if (inp%ci%spectral%emax <= inp%ci%spectral%emin) then; msg = '[spectral] emax must be above emin'; return; end if
      if (inp%ci%spectral%nsteps < 0) then; msg = '[spectral] steps must be >= 0'; return; end if
      inp%ci%spectral%on = .true.
      ! the same bound as class_scf (spectral_final), checked here so that the run fails before the SCF, not after it
      if (2.0_dp * inp%ci%spectral%rmax >= real(inp%ci%lsuper, dp)) then
        msg = '[spectral] rmax must be below [asr] supercell / 2'; return
      end if
    end if

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
    call resolved_input(table, filename, inp%method, inp%resolved)

    if (kspace) then
      ! the k-space input of the same structure
      inp%li%resolved = inp%resolved
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

  !> the species of the class SCF input inp%ci (valence, core shells, atomic start) from the [[species]] of the input file
  subroutine class_species(inp)
    type(run_input_t), intent(inout) :: inp
    integer :: i, nsp
    nsp = size(inp%species)
    if (allocated(inp%ci%spec)) deallocate (inp%ci%spec)
    allocate (inp%ci%spec(nsp))
    do i = 1, nsp
      inp%ci%spec(i)%z = inp%species(i)%z; inp%ci%spec(i)%nval = inp%species(i)%nval; inp%ci%spec(i)%spin_sign = inp%species(i)%spin_sign
      inp%ci%spec(i)%nvl = inp%species(i)%konfig(0:2)
      inp%ci%spec(i)%core = core_shells_for(inp%species(i)%konfig)
      inp%ci%spec(i)%atom_shells = atom_shells(inp%species(i))
    end do
  end subroutine class_species

  subroutine run_asr(inp, ok)
    type(run_input_t), intent(inout) :: inp
    logical, intent(out) :: ok
    type(class_scf_result_t) :: res
    type(rep_run_t) :: rep
    character(len=10) :: start_d, start_t
    integer :: ie, l, p, sg, u, nsp, ncomp, nthr
    character(len=:), allocatable :: dfile
    character(len=16), allocatable :: labels(:)
    real(dp), allocatable :: tsp(:)
    real(dp) :: t0, t_scf, avw, tot, qcell, mcell, fsp

    ok = .false.
    call date_and_time(start_d, start_t)
    t0 = rep_wall()
    nsp = size(inp%species)
    call class_species(inp)
    print '(a,a)', 'asr (augmented-space recursion)  ', inp%name
    print '(a,a,a,f9.5,a,i0,a,i0,a,i0)', '  graph lattice ', trim(inp%ci%lattice), '  a(Bohr)=', inp%ci%alat, '  period=', inp%ci%dec%period, &
      '  sites per decoration cell=', inp%ci%dec%nsite, '  spin=', inp%ci%nspin
    do p = 1, inp%ci%dec%ncls
      print '(a,a8,a,f7.4,a,*(a))', '  sublattice ', inp%subname(p), '  site fraction=', inp%ci%dec%frac(p), '  species: ', &
        (trim(inp%species(inp%ci%spc(ie, p))%symbol)//' ', ie=1, inp%ci%nsc(p))
    end do
    inp%ci%ndos = inp%npoints; inp%ci%dos_emin = inp%emin; inp%ci%dos_emax = inp%emax; inp%ci%dos_eta = inp%eta
    t_scf = rep_wall()
    call class_scf(inp%ci, res, verbose=.true.)
    t_scf = rep_wall() - t_scf
    ncomp = res%ncomp
    allocate (labels(ncomp), tsp(res%nspin))
    avw = avg_ws(trim(inp%ci%lattice)) * inp%ci%alat
    qcell = 0.0_dp; mcell = 0.0_dp
    do p = 1, ncomp
      labels(p) = trim(inp%species(res%comp_spec(p))%symbol)//'@'//trim(inp%subname(res%comp_cls(p)))
      qcell = qcell + res%comp_w(p) * res%comp(p)%qnet * inp%ci%dec%nsite
      mcell = mcell + res%comp_w(p) * res%comp(p)%moment * inp%ci%dec%nsite
    end do
    print '(a,l1,a,i0,a,f10.6,a)', '  SCF converged=', res%converged, '  iterations=', res%iter, '  E_F=', res%ef, ' Ry'
    print '(a,f18.8,a)', '  total energy =', res%etot, ' Ry per site'
    if (.not. res%converged) print '(a)', '  WARNING: SCF did not converge'
    do p = 1, ncomp
      print '(a,a16,a,f7.4,a,f9.5,a,f9.5)', '  component ', labels(p), '  x=', res%comp_x(p), '  qnet=', res%comp(p)%qnet, '  moment=', res%comp(p)%moment
    end do
    print '(a,i0,a,f9.5,a,es10.2)', '  per decoration cell (', inp%ci%dec%nsite, ' sites): moment=', mcell, '  net charge=', qcell

    ! every rank holds the full result (asr_class_scf); rank 0 alone writes the files
    if (.not. par_root()) then
      ok = res%converged; return
    end if

    nthr = 1
    !$ nthr = omp_get_max_threads()
    call rep_from_class_scf(inp, res, rep, t0, nthr, start_d, start_t)

    dfile = inp%workdir//inp%dos_file
    open (newunit=u, file=dfile, status='replace', action='write')
    call rep_write_dos_header(u, rep, inp%npoints, inp%emin, inp%emax)
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
        write (u, '(f12.6,*(es14.6))') res%edos_e(ie), tot, ((fsp * sum(res%edos(ie, l, p, :)), l=0, 2), p=1, ncomp)
      end if
    end do
    close (u)
    print '(a,a)', '  wrote ', dfile
    if (allocated(res%spec_a)) call write_spectral()

    call rep_write_report(rep, inp%workdir//inp%out_file)
    print '(a,a)', '  wrote ', inp%workdir//inp%out_file

    call rep_write_json(rep, inp%workdir//rep%json_file)
    print '(a,a)', '  wrote ', inp%workdir//rep%json_file

    ok = res%converged
  contains
    real(dp) function avg_ws(kind)
      use asr_strconst, only: avg_ws_radius, standard_alpha
      character(len=*), intent(in) :: kind
      avg_ws = avg_ws_radius(kind)
    end function avg_ws

    !> [spectral] file: one block per k point (blank line between blocks), columns path distance, k (2 pi / a), E (Ry), A summed
    !> over spin, A_s A_p A_d, per spin channel (spin 2), per origin sublattice (states/Ry per graph site)
    subroutine write_spectral()
      character(len=:), allocatable :: sfile
      integer :: us, ik, je, kp, nkp, gs, ns, nl, i0, ie
      real(dp) :: fs_
      character(len=13), allocatable :: lab(:)
      integer, allocatable :: wid(:)
      sfile = inp%workdir//inp%spec_file
      fs_ = 2.0_dp / real(res%nspin, dp)
      nkp = size(res%spec_dist)
      open (newunit=us, file=sfile, status='replace', action='write')
      call spectral_header(us, 'Bloch spectral function A(k, E)', .false.)
      write (us, '(a)') '# blocks     : one block of energies per k point, blocks separated by one blank line (gnuplot index)'
      write (us, '(a)') '# columns    : A summed over spins and over s p d unless the label says;  A <sublattice>: part of each'
      write (us, '(a)') '#              origin sublattice (the parts add up to A)'
      ns = size(res%spec_ac, 3)
      nl = 9 + ns + merge(2, 0, res%nspin == 2)
      allocate (lab(nl), wid(nl))
      lab(1:9) = [character(len=13) :: 'dist', 'kx', 'ky', 'kz', 'E [Ry]', 'A', 'A s', 'A p', 'A d']
      wid(1:9) = [10, 9, 9, 9, 11, 13, 13, 13, 13]
      i0 = 9
      if (res%nspin == 2) then
        lab(10:11) = [character(len=13) :: 'A up', 'A down']; wid(10:11) = 13; i0 = 11
      end if
      do kp = 1, ns
        lab(i0 + kp) = 'A ' // trim(inp%subname(kp)); wid(i0 + kp) = 13
      end do
      call rep_label_rows(us, lab, wid)
      do ik = 1, nkp
        do je = 1, size(res%spec_e)
          if (res%nspin == 2) then
            write (us, '(f10.5,3f9.5,f11.6,*(es13.5))') res%spec_dist(ik), res%spec_k(:, ik), res%spec_e(je), &
              fs_ * sum(res%spec_a(je, ik, :, :)), (fs_ * sum(res%spec_a(je, ik, kp, :)), kp=0, 2), &
              (sum(res%spec_a(je, ik, :, gs)), gs=1, 2), (fs_ * sum(res%spec_ac(je, ik, kp, :)), kp=1, size(res%spec_ac, 3))
          else
            write (us, '(f10.5,3f9.5,f11.6,*(es13.5))') res%spec_dist(ik), res%spec_k(:, ik), res%spec_e(je), &
              fs_ * sum(res%spec_a(je, ik, :, :)), (fs_ * sum(res%spec_a(je, ik, kp, :)), kp=0, 2), &
              (fs_ * sum(res%spec_ac(je, ik, kp, :)), kp=1, size(res%spec_ac, 3))
          end if
        end do
        write (us, '(a)') ''
      end do
      close (us)
      print '(a,a)', '  wrote ', sfile
      if (allocated(res%spec_bands)) then              ! ordered, one species: the exact bands of the same H^(2)
        open (newunit=us, file=sfile//'.bands', status='replace', action='write')
        call spectral_header(us, 'eigenvalues of H^(2)(k), from the structure constant blocks and parameters of the recursion', .true.)
        write (us, '(a)') '# columns    : 9 eigenvalues per spin channel, ascending (Ry)'
        nl = 4 + 9 * size(res%spec_bands, 3)
        deallocate (lab, wid); allocate (lab(nl), wid(nl))
        lab(1:4) = [character(len=13) :: 'dist', 'kx', 'ky', 'kz']; wid(1:4) = [10, 9, 9, 9]
        do gs = 1, size(res%spec_bands, 3)
          do ie = 1, 9
            lab(4 + (gs - 1) * 9 + ie) = 'e' // trim(adjustl(itoa_(ie)))
            if (size(res%spec_bands, 3) == 2) lab(4 + (gs - 1) * 9 + ie) = trim(lab(4 + (gs - 1) * 9 + ie)) // merge(' up  ', ' down', gs == 1)
            wid(4 + (gs - 1) * 9 + ie) = 11
          end do
        end do
        call rep_label_rows(us, lab, wid)
        do ik = 1, nkp
          write (us, '(f10.5,3f9.5,*(f11.6))') res%spec_dist(ik), res%spec_k(:, ik), res%spec_bands(:, ik, :)
        end do
        close (us)
        print '(a,a)', '  wrote ', sfile//'.bands'
      end if
    end subroutine write_spectral

    !> common header of the spectral files: run, route, E_F, path corners, units
    subroutine spectral_header(us, what, bands)
      integer, intent(in) :: us
      character(len=*), intent(in) :: what
      logical, intent(in) :: bands                    !! the .bands file: no energy mesh, no A
      integer :: ic
      write (us, '(a)') '# ASR ' // what
      write (us, '(a)') '# run        : ' // trim(rep%title) // '   (' // trim(rep%infile) // ', git ' // trim(rep%git_rev) // ')'
      write (us, '(a,i0,a,f8.5,a,i0,a,f6.3,a)') '# route      : recursion edc, levels ', inp%ci%nlevels, ';  eta ', inp%ci%spectral%eta, &
        ' Ry;  chain steps ', res%spec_nsteps, ';  rmax ', inp%ci%spectral%rmax, ' alat;  window ' // trim(inp%ci%spectral%window)
      write (us, '(a,f12.8,a)') '# E_F        : ', res%ef, ' Ry'
      if (bands) then
        write (us, '(a,i0,a)') '# mesh       : ', size(res%spec_dist), ' k points'
        write (us, '(a)') '# units      : k Cartesian in 2 pi / a; dist along the path in 2 pi / a; eigenvalues in Ry'
      else
        write (us, '(a,i0,a,i0,a)') '# mesh       : ', size(res%spec_dist), ' k points x ', size(res%spec_e), ' energies'
        write (us, '(a)') '# units      : k Cartesian in 2 pi / a; dist along the path in 2 pi / a; E in Ry; A in states/Ry per graph site'
      end if
      write (us, '(a)') '# corners    :  label   block   dist        kx       ky       kz'
      do ic = 1, size(res%spec_icorner)
        write (us, '(a,a8,i6,f10.5,3f9.5)') '#              ', inp%ci%spectral%labels(ic), res%spec_icorner(ic), &
          res%spec_dist(res%spec_icorner(ic)), inp%ci%spectral%kpts(:, ic)
      end do
    end subroutine spectral_header

    function itoa_(i) result(t)
      integer, intent(in) :: i
      character(len=12) :: t
      write (t, '(i0)') i
    end function itoa_
  end subroutine run_asr

  !> Fill rep_run_t from recursion class_scf_result_t
  subroutine rep_from_class_scf(inp, res, rep, t0, nthr, start_d, start_t)
    use iso_fortran_env, only: compiler_version, compiler_options
    use asr_build_info, only: git_rev, git_dirty
    use asr_kinds, only: dp, pi
    use asr_mpi, only: par_size
    use asr_strconst, only: avg_ws_radius, standard_alpha
    use asr_contour, only: fermi_npole, fermi_nodes_count
    use asr_runstats, only: cpu_seconds
    type(run_input_t), intent(in) :: inp
    type(class_scf_result_t), intent(in) :: res
    type(rep_run_t), intent(out) :: rep
    real(dp), intent(in) :: t0
    integer, intent(in) :: nthr
    character(len=10), intent(in) :: start_d, start_t
    real(dp) :: fsp, avw
    integer :: p, sg, l, kap, k, i, j, n1, n2, nmax, dot

    rep%prog_name = 'asr'
    rep%prog_version = '0.1.0'
    rep%git_rev = git_rev
    rep%git_dirty = git_dirty
    rep%compiler_ver = compiler_version()
    rep%compiler_opts = compiler_options()
    rep%host = get_host()
    rep%pid = get_pid()
    rep%start_date = start_d
    rep%start_time = start_t
    call date_and_time(rep%finish_date, rep%finish_time)
    rep%wall_total = rep_wall() - t0
    rep%mpi_ranks = par_size()
    rep%omp_threads = nthr
    rep%infile = inp%infile
    rep%fingerprint = fnv1a64_file(inp%infile)
    rep%workdir = inp%workdir
    rep%out_file = inp%out_file
    dot = index(inp%out_file, '.', back=.true.)
    if (dot > 0) then
      rep%json_file = inp%out_file(1:dot - 1)//'.json'
    else
      rep%json_file = trim(inp%out_file)//'.json'
    end if
    rep%dos_file = inp%dos_file
    if (len_trim(inp%ci%pot_out) > 0) rep%pot_file = inp%ci%pot_out
    rep%title = inp%name
    if (allocated(inp%resolved)) rep%resolved = inp%resolved

    ! Method
    rep%method = 'asr'
    rep%route = trim(inp%ci%route)
    rep%xc = trim(inp%ci%xc)
    rep%relativistic = inp%ci%relativistic
    rep%nspin = inp%ci%nspin
    rep%ncomp = res%ncomp
    rep%nlevels = inp%ci%nlevels
    rep%nrec = merge(inp%ci%nrec, inp%ci%nlevels, inp%ci%nrec > 0)
    rep%maxrank = inp%ci%maxrank
    rep%lsuper = inp%ci%lsuper
    if (allocated(res%nodes)) rep%nodes = res%nodes
    allocate (rep%subname(inp%ci%dec%ncls))
    rep%subname = inp%subname
    rep%s_range = inp%ci%rcl
    rep%inv_radius = max(inp%ci%rcl, inp%ci%rinv)            ! rinv = 0: the inversion cluster is the kept range
    rep%alpha = standard_alpha() * inp%ci%alpha_scale
    rep%alpha_scale = inp%ci%alpha_scale
    rep%terminator = inp%ci%ed_tail
    rep%depth_tol = inp%ci%depth_tol
    ! occupation contour as asr_class_scf builds it (contour): Fermi function with poles, or a closed contour without smearing
    rep%contour_main = inp%ci%ncont
    if (inp%ci%fermi_kt > 0.0_dp) then
      rep%contour_poles = fermi_npole(inp%ci%fermi_kt)
      if (inp%ci%fermi_poles > 0) rep%contour_poles = inp%ci%fermi_poles
      rep%contour_nodes = fermi_nodes_count(rep%contour_poles, inp%ci%ncont, 16)
    else
      rep%contour_poles = 0
      rep%contour_nodes = inp%ci%ncont + 8
    end if
    rep%fermi_kt = inp%ci%fermi_kt
    rep%fermi_poles = inp%ci%fermi_poles
    rep%ebot_offset = inp%ci%ebot_offset
    rep%eta = inp%eta
    rep%ed_eta = inp%ci%ed_eta
    rep%sim_beta = inp%ci%sim_beta
    rep%has_screen = .true.
    rep%unscreened = res%screen%unscreened
    if (.not. res%screen%unscreened) then
      rep%screen_a = res%screen%a
      rep%screen_g00 = res%screen%g00
      rep%screen_g1 = res%screen%g1
      rep%screen_w00_g00 = res%screen%w00 / res%screen%g00
      rep%screen_lambda = res%screen%lambda
      rep%screen_beta_max = res%screen%beta_max
    end if

    ! SRO
    rep%sro%has_sro = inp%ci%sro
    if (inp%ci%sro) then
      if (inp%ci%dec%ncls == 1) then
        rep%sro%sub_centre = inp%subname(1)
        rep%sro%sub_neigh = inp%subname(1)
        rep%sro%coord_num = res%screen%z
        nmax = inp%ci%nsc(1)
        rep%sro%nspec = nmax
        allocate (rep%sro%species(nmax), rep%sro%cond_p(nmax, nmax), rep%sro%alpha(nmax, nmax))
        do i = 1, nmax
          rep%sro%species(i) = inp%species(inp%ci%spc(i, 1))%symbol
          do j = 1, nmax
            rep%sro%cond_p(i, j) = inp%ci%sro_cond(i, j, 1, 1)
            rep%sro%alpha(i, j) = 1.0_dp - inp%ci%sro_cond(i, j, 1, 1) / max(1.0e-12_dp, inp%ci%conc(j, 1))
          end do
        end do
      else
        find_pair: do kap = 1, inp%ci%dec%ncls
          do k = 1, inp%ci%dec%ncls
            if (any(abs(inp%ci%sro_cond(1:inp%ci%nsc(kap), 1:inp%ci%nsc(k), kap, k) - &
                        spread(inp%ci%conc(1:inp%ci%nsc(k), k), 1, inp%ci%nsc(kap))) > 1.0e-12_dp)) then
              rep%sro%sub_centre = inp%subname(kap)
              rep%sro%sub_neigh = inp%subname(k)
              rep%sro%coord_num = res%screen%z
              n1 = inp%ci%nsc(kap); n2 = inp%ci%nsc(k)
              rep%sro%nspec = n1
              allocate (rep%sro%species(n1), rep%sro%cond_p(n1, n2), rep%sro%alpha(n1, n2))
              do i = 1, n1
                rep%sro%species(i) = inp%species(inp%ci%spc(i, kap))%symbol
                do j = 1, n2
                  rep%sro%cond_p(i, j) = inp%ci%sro_cond(i, j, kap, k)
                  rep%sro%alpha(i, j) = 1.0_dp - inp%ci%sro_cond(i, j, kap, k) / max(1.0e-12_dp, inp%ci%conc(j, k))
                end do
              end do
              exit find_pair
            end if
          end do
        end do find_pair
      end if
      if (allocated(res%mode_lam)) then
        rep%sro%nmodes = size(res%mode_lam)
        allocate (rep%sro%mode_lam(rep%sro%nmodes), rep%sro%mode_gam(rep%sro%nmodes), &
                  rep%sro%mode_gform(rep%sro%nmodes), rep%sro%mode_f(rep%sro%nmodes), &
                  rep%sro%mode_inf(rep%sro%nmodes))
        rep%sro%mode_lam = res%mode_lam
        rep%sro%mode_gam = res%mode_gam
        rep%sro%mode_gform = res%screen%g00 + res%screen%g1 * res%mode_lam
        do i = 1, rep%sro%nmodes
          rep%sro%mode_inf(i) = (res%mode_gam(i) <= 1.0e-9_dp * res%screen%g00)
          if (.not. rep%sro%mode_inf(i)) then
            rep%sro%mode_f(i) = mode_f(res%screen, res%mode_lam(i), res%mode_gam(i))
          else
            rep%sro%mode_f(i) = 0.0_dp
          end if
        end do
      end if
      rep%sro%uform = res%uform
      if (allocated(res%udrive)) rep%sro%udrive = res%udrive
    end if

    ! Mixing
    rep%mix_type = 'Anderson'
    rep%mix_beta = inp%ci%beta
    rep%mix_history = inp%ci%mix_history
    rep%mix_restart = inp%ci%beta_restart
    rep%dvmax = 0.15_dp
    rep%denumax = 0.08_dp
    rep%tol = inp%ci%tol

    ! Structure
    rep%lattice = trim(inp%ci%lattice)
    rep%alat = inp%ci%alat
    rep%period = inp%ci%dec%period
    rep%nsite = inp%ci%dec%nsite
    avw = avg_ws_radius(trim(inp%ci%lattice)) * inp%ci%alat
    rep%sws_avg = avw
    rep%cell_vol = (4.0_dp * pi / 3.0_dp) * avw**3 * inp%ci%dec%nsite
    rep%nsub = inp%ci%dec%ncls
    allocate (rep%sub(rep%nsub))
    do kap = 1, rep%nsub
      rep%sub(kap)%name = inp%subname(kap)
      rep%sub(kap)%frac = inp%ci%dec%frac(kap)
      rep%sub(kap)%nspec = inp%ci%nsc(kap)
      allocate (rep%sub(kap)%species(rep%sub(kap)%nspec), rep%sub(kap)%conc(rep%sub(kap)%nspec))
      do i = 1, rep%sub(kap)%nspec
        rep%sub(kap)%species(i) = inp%species(inp%ci%spc(i, kap))%symbol
        rep%sub(kap)%conc(i) = inp%ci%conc(i, kap)
      end do
    end do

    ! Components
    allocate (rep%comp(res%ncomp))
    fsp = 2.0_dp / real(res%nspin, dp)
    do p = 1, res%ncomp
      rep%comp(p)%label = trim(inp%species(res%comp_spec(p))%symbol)//'@'//trim(inp%subname(res%comp_cls(p)))
      rep%comp(p)%species = inp%species(res%comp_spec(p))%symbol
      rep%comp(p)%subname = inp%subname(res%comp_cls(p))
      rep%comp(p)%z = inp%species(res%comp_spec(p))%z
      rep%comp(p)%valence = inp%species(res%comp_spec(p))%nval
      rep%comp(p)%concentration = res%comp_x(p)
      rep%comp(p)%weight_per_site = res%comp_w(p)
      rep%comp(p)%sws = avw
      rep%comp(p)%qnet = res%comp(p)%qnet
      rep%comp(p)%moment = res%comp(p)%moment
      rep%comp(p)%enu = res%comp(p)%enu
      do sg = 1, res%nspin
        do l = 0, 2
          rep%comp(p)%c(l, sg) = res%comp(p)%pp(l, sg)%c
          rep%comp(p)%sqrdel(l, sg) = res%comp(p)%pp(l, sg)%sqrdel
          rep%comp(p)%p(l, sg) = res%comp(p)%pp(l, sg)%p
          rep%comp(p)%gamma(l, sg) = res%comp(p)%pp(l, sg)%gamma
        end do
      end do
      rep%comp(p)%q = fsp * res%comp(p)%q
      do sg = 1, res%nspin
        do l = 0, 2
          if (allocated(res%qch)) then
            rep%comp(p)%qch(l, sg) = fsp * res%qch(p, sg, l)
          else
            rep%comp(p)%qch(l, sg) = rep%comp(p)%q(0, l, sg)
          end if
          if (allocated(res%dosef)) then
            rep%comp(p)%dosef(l, sg) = fsp * res%dosef(p, sg, l)
          end if
        end do
      end do
      if (allocated(res%ecomp)) then
        rep%comp(p)%ecomp = res%ecomp(:, p)
      end if
    end do

    ! SCF history and observables
    rep%iter = res%iter
    rep%converged = res%converged
    rep%ef = res%ef
    rep%etot = res%etot
    rep%eterm = res%eterm
    if (allocated(res%ecomp)) rep%ecomp = res%ecomp
    if (allocated(res%hist_etot)) rep%hist_etot = res%hist_etot
    if (allocated(res%hist_dv)) rep%hist_dv = res%hist_dv
    if (allocated(res%hist_denu)) rep%hist_denu = res%hist_denu
    if (allocated(res%hist_ef)) rep%hist_ef = res%hist_ef
    if (allocated(res%hist_npass)) rep%hist_npass = res%hist_npass
    if (allocated(res%hist_wall)) rep%hist_wall = res%hist_wall
    if (allocated(res%hist_moment)) rep%hist_moment = res%hist_moment
    if (allocated(res%hist_qnet)) rep%hist_qnet = res%hist_qnet
    rep%qcell = res%qcell
    rep%moment_site = res%moment_site
    rep%dosef_site = sum([(res%comp_w(p) * sum(rep%comp(p)%dosef), p=1, res%ncomp)])

    ! Diagnostics
    rep%nclamp = res%nclamp
    rep%nclamp_run = res%nclamp_run
    rep%nnoncausal = res%nnoncausal
    rep%nnode = res%nnode
    rep%wbelow1 = res%wbelow1
    rep%wbelow = res%wbelow
    rep%below_ok = res%below_ok
    rep%depth_mean = res%depth_mean
    rep%nfall = res%nfall
    rep%mix_nrestart = res%mix_nrestart
    rep%mix_nh = res%mix_nh

    ! Bloch spectral function
    if (allocated(res%spec_a)) then
      rep%spec%on = .true.
      rep%spec%file = inp%spec_file
      rep%spec%has_bands = allocated(res%spec_bands)
      rep%spec%nk = size(res%spec_dist); rep%spec%ne = size(res%spec_e)
      rep%spec%emin = res%spec_e(1); rep%spec%emax = res%spec_e(size(res%spec_e))
      rep%spec%eta = inp%ci%spectral%eta; rep%spec%rmax = inp%ci%spectral%rmax; rep%spec%window = inp%ci%spectral%window
      rep%spec%nsteps = res%spec_nsteps; rep%spec%nt = res%spec_nt
      rep%spec%maxhop = res%spec_maxhop; rep%spec%hoplim = res%spec_hoplim; rep%spec%nneg = res%spec_nneg
      rep%spec%rdef = res%spec_rdef; rep%spec%spread = res%spec_spread; rep%spec%spread0 = res%spec_spread0
      rep%spec%amin = res%spec_amin
      if (allocated(res%spec_nodes)) rep%spec%nodes = res%spec_nodes
      if (allocated(res%spec_rsh)) rep%spec%rsh = res%spec_rsh
      if (allocated(res%spec_decay)) rep%spec%decay = res%spec_decay
      if (allocated(res%spec_sum)) rep%spec%sumrule = transpose(res%spec_sum)
      rep%spec%clabel = inp%ci%spectral%labels
      rep%spec%cdist = res%spec_dist(res%spec_icorner)
      rep%spec%ck = inp%ci%spectral%kpts
    end if

    ! Timing: the parts from the SCF; the total is the whole run so far (input, SCF, DOS), as the end-of-run line
    rep%timer = res%timer
    rep%timer%wall_total = rep%wall_total
    rep%timer%cpu_total = cpu_seconds()
  end subroutine rep_from_class_scf
end module asr_driver
