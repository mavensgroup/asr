!> Driver for [calculation] kind = "sublattice": ordered compounds, partial (sublattice) disorder and multicomponent alloys on a
!> decorated sc / bcc / fcc lattice, self-consistent from the energy-dependent augmented-space recursion (asr_class_scf), then the
!> component- and l-resolved DOS.
!>
!>   [calculation]  name, kind = "sublattice"
!>   [structure]    lattice = "sc"|"bcc"|"fcc" (the graph lattice: every site of it belongs to one sublattice), alat (Angstrom, of the
!>                  graph lattice), period (decoration cell = period^3 conventional cubes; default 1), centering = "P"|"I"|"F"
!>                  (translations of the decoration cell applied to every position; default "P")
!>   [[species]]    symbol, z, nval, konfig, occupations, spin_sign (as kind = "alloy")
!>   [[sublattice]] name, positions = [[x, y, z], ...] (fractional coordinates of the decoration cell), species = ["A", ...],
!>                  concentrations = [x_A, ...] (sum 1; may be omitted for one species)
!>   [scf]          spin, xc, sim_beta, mixing, tolerance, max_iterations, zeta0, relativistic, save_potential / start_potential
!>   [asr]          levels, recursion_steps, rank, supercell (multiple of period), flip_radius, contour_points, cluster_radius,
!>                  inversion_radius, ed_terminator (hull|bp|avg2), fermi_kt, fermi_poles, depth_tol
!>   [dos]          emin, emax, npoints, eta
!>   [output]       dos_file, out_file
!> Example (L2_1 Co2Fe0.4Cr0.6Al on the bcc graph, a_bcc = a/2): period = 2, centering = "F"; Al (0,0,0); Fe/Cr (1/2,1/2,1/2);
!> Co (1/4,1/4,1/4) and (3/4,3/4,3/4).  Unknown sections/keys are errors.
module asr_class_driver
  use asr_kinds, only: dp
  use tomlf, only: toml_table, toml_array, toml_error, toml_load, get_value, len
  use asr_tomlutil, only: section, check_keys
  use asr_input_lmto, only: lmto_species_t
  use asr_asa, only: core_shells_for, asa_class_result_t
  use asr_lmto_driver, only: atom_shells
  use asr_sublattice, only: decoration_build
  use asr_class_scf
  use asr_report
  !$ use omp_lib, only: omp_get_max_threads
  implicit none
  private
  public :: class_run_input_t, read_class_input, run_class

  type :: class_run_input_t
    character(len=:), allocatable :: name, dos_file, out_file, workdir, infile
    type(class_scf_input_t) :: ci
    type(lmto_species_t), allocatable :: species(:)
    character(len=16), allocatable :: subname(:)
    real(dp) :: emin = -1.0_dp, emax = 1.0_dp, eta = 0.02_dp
    integer :: npoints = 801
  end type class_run_input_t

  real(dp), parameter :: bohr_per_angstrom = 1.0_dp / 0.529177210903_dp

contains

  subroutine read_class_input(filename, inp, msg)
    character(len=*), intent(in) :: filename
    type(class_run_input_t), intent(out) :: inp
    character(len=:), allocatable, intent(out) :: msg
    type(toml_table), allocatable :: table
    type(toml_table), pointer :: child, sp
    type(toml_array), pointer :: arr, ka, kb
    type(toml_error), allocatable :: perr
    integer :: i, j, k, slash, nsp, nsub, npos, nsc, period, nmax
    real(dp) :: x, alat
    real(dp), allocatable :: pos(:, :), tpos(:, :)
    integer, allocatable :: owner(:), towner(:)
    character(len=:), allocatable :: lat, xcs, edt, cen, str

    msg = ''
    call toml_load(table, filename, error=perr)
    if (allocated(perr)) then; msg = "cannot read '"//filename//"': "//perr%message; return; end if
    inp%infile = filename
    slash = index(filename, '/', back=.true.)
    if (slash > 0) then; inp%workdir = filename(1:slash); else; inp%workdir = ''; end if
    call check_keys(table, ['calculation', 'structure  ', 'species    ', 'sublattice ', 'scf        ', 'asr        ', 'dos        ', &
                            'output     '], 'top level', msg)
    if (len(msg) > 0) return

    call section(table, 'calculation', child)
    call check_keys(child, ['name', 'kind'], '[calculation]', msg); if (len(msg) > 0) return
    call get_value(child, 'name', inp%name, 'asr-sublattice')

    call section(table, 'structure', child)
    call check_keys(child, ['lattice  ', 'alat     ', 'period   ', 'centering'], '[structure]', msg); if (len(msg) > 0) return
    call get_value(child, 'lattice', lat, 'bcc')
    call get_value(child, 'alat', alat, 0.0_dp)
    call get_value(child, 'period', period, 1)
    call get_value(child, 'centering', cen, 'P')
    inp%ci%lattice = lat; inp%ci%alat = alat * bohr_per_angstrom
    if (lat /= 'sc' .and. lat /= 'bcc' .and. lat /= 'fcc') then; msg = '[structure] lattice must be sc|bcc|fcc'; return; end if
    if (alat <= 0.0_dp) then; msg = '[structure] alat (Angstrom) must be > 0'; return; end if
    if (period < 1) then; msg = '[structure] period must be >= 1'; return; end if

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

    call section(table, 'scf', child)
    call check_keys(child, ['spin           ', 'xc             ', 'sim_beta       ', 'mixing         ', 'tolerance      ', &
                            'max_iterations ', 'zeta0          ', 'relativistic   ', 'start_potential', 'save_potential '], '[scf]', msg)
    if (len(msg) > 0) return
    call get_value(child, 'spin', inp%ci%nspin, 1)
    call get_value(child, 'xc', xcs, 'vbh'); inp%ci%xc = xcs
    call get_value(child, 'sim_beta', inp%ci%sim_beta, 1.0_dp)
    call get_value(child, 'mixing', inp%ci%beta, 0.3_dp)
    call get_value(child, 'tolerance', inp%ci%tol, 1.0e-5_dp)
    call get_value(child, 'max_iterations', inp%ci%maxit, 60)
    call get_value(child, 'zeta0', inp%ci%zeta0, 0.1_dp)
    call get_value(child, 'relativistic', inp%ci%relativistic, .true.)
    call get_value(child, 'start_potential', xcs, ''); inp%ci%pot_in = xcs
    call get_value(child, 'save_potential', xcs, ''); inp%ci%pot_out = xcs

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

    call section(table, 'dos', child)
    call check_keys(child, ['emin   ', 'emax   ', 'npoints', 'eta    '], '[dos]', msg); if (len(msg) > 0) return
    call get_value(child, 'emin', inp%emin, -1.0_dp)
    call get_value(child, 'emax', inp%emax, 1.0_dp)
    call get_value(child, 'npoints', inp%npoints, 401)
    call get_value(child, 'eta', inp%eta, 0.02_dp)

    call section(table, 'output', child)
    call check_keys(child, ['dos_file', 'out_file'], '[output]', msg); if (len(msg) > 0) return
    call get_value(child, 'dos_file', inp%dos_file, 'dos.dat')
    call get_value(child, 'out_file', inp%out_file, 'asr.out')

    ! ---- validation
    if (inp%ci%nspin /= 1 .and. inp%ci%nspin /= 2) msg = '[scf] spin must be 1 or 2'
    if (inp%ci%nlevels < 2) msg = '[asr] levels must be >= 2'
    if (inp%ci%nrec /= 0 .and. inp%ci%nrec < inp%ci%nlevels) msg = '[asr] recursion_steps must be 0 or >= levels'
    if (inp%ci%maxrank < 0) msg = '[asr] rank must be >= 0'
    if (modulo(inp%ci%lsuper, period) /= 0) msg = '[asr] supercell must be a multiple of [structure] period'
    if (inp%ci%fermi_poles < 0 .or. inp%ci%fermi_poles > 200) msg = '[asr] fermi_poles must be 0 (automatic) .. 200'
    if (edt /= 'hull' .and. edt /= 'bp' .and. edt /= 'avg2') msg = "[asr] ed_terminator must be 'hull', 'bp' or 'avg2'"
    if (inp%emax <= inp%emin) msg = '[dos] emax must exceed emin'
  end subroutine read_class_input

  subroutine run_class(inp, ok)
    type(class_run_input_t), intent(inout) :: inp
    logical, intent(out) :: ok
    type(class_scf_result_t) :: res
    integer :: i, ie, l, p, sg, u, nsp, ncomp, ru, nthr, kap
    character(len=:), allocatable :: dfile
    character(len=8), allocatable :: labels(:)
    real(dp), allocatable :: zc(:), swsc(:), tsp(:)
    integer, allocatable :: nvc(:), natc(:)
    real(dp) :: t0, t_scf, avw, tot, qcell, mcell

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
    print '(a,a)', 'asr sublattice alloy  ', inp%name
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
    do ie = 1, inp%npoints
      tsp = 0.0_dp
      do sg = 1, res%nspin
        tsp(sg) = sum([(res%comp_w(p) * sum(res%edos(ie, :, p, sg)), p=1, ncomp)])
      end do
      tot = sum(tsp)
      write (u, '(f12.6,*(es14.6))') res%edos_e(ie), tot, ((sum(res%edos(ie, l, p, :)), l=0, 2), p=1, ncomp), tsp
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
  end subroutine run_class
end module asr_class_driver
