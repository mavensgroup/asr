!> Driver for [calculation] kind = "alloy": binary substitutional alloy A_x B_(1-x) on fcc/bcc, self-consistent
!> from the augmented-space recursion, then species- and l-resolved DOS.
!>
!>   [calculation]  name, kind = "alloy"
!>   [structure]    lattice = "fcc"|"bcc" (random alloy) | "b2" (ordered CsCl, 50:50; no concentration), alat (Angstrom)
!>   [[species]]    (two) symbol, z, nval, konfig, occupations, concentration (A only; B = 1 - x_A), spin_sign (+-1: start polarisation)
!>   [scf]          spin = 1|2, xc, sigma, sim_beta, mixing, tolerance, max_iterations, zeta0, relativistic,
!>                  save_potential / start_potential (file), flip_spin ("A"|"B"|"AB": swap that species' spins in the start)
!>   [asr]          levels, rank, supercell, flip_radius, contour_points, cluster_radius, inversion_radius, screening, terminator,
!>                  density, order (2|3: second/third-order H), symmetry, hamiltonian, ed_terminator (hull|bp|avg2),
!>                  recursion_steps (Lanczos steps on the generated node set, hard wall; default = levels)
!>   [dos]          emin, emax, npoints, eta
!>   [output]       dos_file, out_file
!> Unknown sections/keys are errors.
module asr_alloy_driver
  use asr_kinds, only: dp
  use tomlf, only: toml_table, toml_array, toml_error, toml_load, get_value, len
  use asr_tomlutil, only: section, check_keys
  use asr_input_lmto, only: lmto_species_t
  use asr_asa, only: core_shells_for, asa_class_result_t
  use asr_lmto_driver, only: atom_shells
  use asr_alloy_scf
  use asr_report
  !$ use omp_lib, only: omp_get_max_threads
  implicit none
  private
  public :: alloy_run_input_t, read_alloy_input, run_alloy

  type :: alloy_run_input_t
    character(len=:), allocatable :: name, dos_file, out_file, workdir, infile
    type(alloy_input_t) :: ai
    type(lmto_species_t) :: species(2)
    real(dp) :: emin = -1.0_dp, emax = 1.0_dp, eta = 0.02_dp
    integer :: npoints = 801
  end type alloy_run_input_t

  real(dp), parameter :: bohr_per_angstrom = 1.0_dp / 0.529177210903_dp

contains

  subroutine read_alloy_input(filename, inp, msg)
    character(len=*), intent(in) :: filename
    type(alloy_run_input_t), intent(out) :: inp
    character(len=:), allocatable, intent(out) :: msg
    type(toml_table), allocatable :: table
    type(toml_table), pointer :: child, sp
    type(toml_array), pointer :: arr, ka
    type(toml_error), allocatable :: perr
    integer :: i, j, slash, nsp
    real(dp) :: x, alat, xs(2)
    character(len=:), allocatable :: lat, xcs, hamil, scr, trm, dns, chg, pfs, edt

    msg = ''
    call toml_load(table, filename, error=perr)
    if (allocated(perr)) then; msg = "cannot read '"//filename//"': "//perr%message; return; end if
    inp%infile = filename
    slash = index(filename, '/', back=.true.)
    if (slash > 0) then; inp%workdir = filename(1:slash); else; inp%workdir = ''; end if
    call check_keys(table, ['calculation', 'structure  ', 'species    ', 'scf        ', 'asr        ', 'dos        ', 'output     '], 'top level', msg)
    if (len(msg) > 0) return

    call section(table, 'calculation', child)
    call check_keys(child, ['name', 'kind'], '[calculation]', msg); if (len(msg) > 0) return
    call get_value(child, 'name', inp%name, 'asr-alloy')

    call section(table, 'structure', child)
    call check_keys(child, ['lattice', 'alat   '], '[structure]', msg); if (len(msg) > 0) return
    call get_value(child, 'lattice', lat, 'fcc')
    call get_value(child, 'alat', alat, 0.0_dp)
    inp%ai%lattice = lat; inp%ai%alat = alat * bohr_per_angstrom
    inp%ai%ordered = (lat == 'b2')
    if (inp%ai%ordered) inp%ai%lattice = 'bcc'      ! B2 = ordered binary on the bcc graph (A on the corner, B on the body centre)

    nullify (arr)
    call get_value(table, 'species', arr, requested=.false.)
    if (.not. associated(arr)) then; msg = 'no [[species]] given'; return; end if
    nsp = len(arr)
    if (nsp /= 2) then; msg = 'kind = "alloy" needs exactly two [[species]]'; return; end if
    xs = -1.0_dp
    do i = 1, 2
      call get_value(arr, i, sp)
      call check_keys(sp, ['symbol       ', 'z            ', 'nval         ', 'konfig       ', 'occupations  ', 'concentration', &
                            'spin_sign    '], '[[species]]', msg)
      if (len(msg) > 0) return
      call get_value(sp, 'symbol', inp%species(i)%symbol, 'X')
      call get_value(sp, 'z', inp%species(i)%z, 0.0_dp)
      call get_value(sp, 'nval', inp%species(i)%nval, 0)
      call get_value(sp, 'spin_sign', x, 1.0_dp); inp%species(i)%spin_sign = sign(1.0_dp, x)
      call get_value(sp, 'concentration', xs(i), -1.0_dp)
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
    if (inp%ai%ordered) then
      if (xs(1) >= 0.0_dp .and. abs(xs(1) - 0.5_dp) > 1.0e-8_dp) then; msg = 'lattice b2 is the 50:50 ordered crystal; omit concentration'; return; end if
      xs(1) = 0.5_dp
    else
      if (xs(1) < 0.0_dp) then; msg = 'first [[species]] needs concentration'; return; end if
      if (xs(2) >= 0.0_dp .and. abs(xs(1) + xs(2) - 1.0_dp) > 1.0e-8_dp) then; msg = 'concentrations must sum to 1'; return; end if
    end if
    inp%ai%x = xs(1)

    call section(table, 'scf', child)
    call check_keys(child, ['spin           ', 'xc             ', 'sigma          ', 'sim_beta       ', 'mixing         ', &
                            'tolerance      ', 'max_iterations ', 'zeta0          ', 'relativistic   ', &
                            'start_potential', 'save_potential ', 'flip_spin      '], '[scf]', msg)
    if (len(msg) > 0) return
    call get_value(child, 'spin', inp%ai%nspin, 1)
    call get_value(child, 'xc', xcs, 'vbh'); inp%ai%xc = xcs
    call get_value(child, 'sigma', inp%ai%sigma, 0.7_dp)
    call get_value(child, 'sim_beta', inp%ai%sim_beta, 1.0_dp)
    call get_value(child, 'mixing', inp%ai%beta, 0.3_dp)
    call get_value(child, 'tolerance', inp%ai%tol, 1.0e-5_dp)
    call get_value(child, 'max_iterations', inp%ai%maxit, 60)
    call get_value(child, 'zeta0', inp%ai%zeta0, 0.1_dp)
    call get_value(child, 'relativistic', inp%ai%relativistic, .true.)
    call get_value(child, 'start_potential', xcs, ''); inp%ai%pot_in = xcs
    call get_value(child, 'save_potential', xcs, ''); inp%ai%pot_out = xcs
    call get_value(child, 'flip_spin', xcs, '')             ! species whose spin channels are swapped in the start potential: "A", "B", "AB"
    inp%ai%pot_flip = [index(xcs, 'A') > 0, index(xcs, 'B') > 0]
    if (any(inp%ai%pot_flip) .and. len_trim(inp%ai%pot_in) == 0) then; msg = 'flip_spin needs start_potential'; return; end if

    call section(table, 'asr', child)
    call check_keys(child, ['levels            ', 'rank              ', 'supercell         ', 'flip_radius       ', 'contour_points    ', &
                            'cluster_radius    ', 'symmetry          ', 'hamiltonian       ', 'inversion_radius  ', 'screening         ', 'terminator        ', 'density           ', 'order             ', 'order_switch      ', 'charges           ', 'depth_tol         ', 'potential_function', 'exact_density     ', 'kkr_check         ', 'ed_terminator     ', 'recursion_steps   '], '[asr]', msg)
    if (len(msg) > 0) return
    call get_value(child, 'levels', inp%ai%nlevels, 4)
    call get_value(child, 'recursion_steps', inp%ai%nrec, 0)
    call get_value(child, 'rank', inp%ai%maxrank, 2)
    call get_value(child, 'supercell', inp%ai%lsuper, 8)
    call get_value(child, 'flip_radius', inp%ai%flip_radius, 99)
    call get_value(child, 'contour_points', inp%ai%ncont, 60)
    call get_value(child, 'cluster_radius', inp%ai%rcl, inp%ai%rcl)
    call get_value(child, 'inversion_radius', inp%ai%rinv, inp%ai%rinv)
    call get_value(child, 'screening', scr, 'auto'); inp%ai%screening = scr
    call get_value(child, 'terminator', trm, 'avg2'); inp%ai%terminator = trm
    call get_value(child, 'order', inp%ai%order, inp%ai%order)
    call get_value(child, 'order_switch', inp%ai%order_switch, inp%ai%order_switch)
    call get_value(child, 'charges', chg, 'lowdin'); inp%ai%charges = chg
    call get_value(child, 'depth_tol', inp%ai%depth_tol, inp%ai%depth_tol)
    call get_value(child, 'ed_terminator', edt, 'hull'); inp%ai%ed_tail = edt
    call get_value(child, 'potential_function', pfs, 'parametrised'); inp%ai%ed_exact = (pfs == 'exact')
    call get_value(child, 'exact_density', inp%ai%ed_exact_density, inp%ai%ed_exact_density)
    call get_value(child, 'kkr_check', inp%ai%kkr_check, inp%ai%kkr_check)
    call get_value(child, 'density', dns, 'expand'); inp%ai%density = dns
    call get_value(child, 'symmetry', inp%ai%symmetry, .true.)
    call get_value(child, 'hamiltonian', hamil, 'energy-independent')
    inp%ai%energy_dependent = (hamil == 'energy-dependent' .or. hamil == 'asr_ed')
    inp%ai%ed_asr = (hamil == 'asr_ed')

    call section(table, 'dos', child)
    call check_keys(child, ['emin   ', 'emax   ', 'npoints', 'eta    '], '[dos]', msg); if (len(msg) > 0) return
    call get_value(child, 'emin', inp%emin, -1.0_dp)
    call get_value(child, 'emax', inp%emax, 1.0_dp)
    call get_value(child, 'npoints', inp%npoints, 801)
    call get_value(child, 'eta', inp%eta, 0.02_dp)

    call section(table, 'output', child)
    call check_keys(child, ['dos_file', 'out_file'], '[output]', msg); if (len(msg) > 0) return
    call get_value(child, 'dos_file', inp%dos_file, 'dos.dat')
    call get_value(child, 'out_file', inp%out_file, 'asr.out')

    ! ---- validation
    if (lat /= 'fcc' .and. lat /= 'bcc' .and. lat /= 'b2') msg = '[structure] lattice must be fcc|bcc|b2'
    if (alat <= 0.0_dp) msg = '[structure] alat (Angstrom) must be > 0'
    if (.not. inp%ai%ordered .and. (inp%ai%x <= 0.0_dp .or. inp%ai%x >= 1.0_dp)) msg = 'concentration must be strictly between 0 and 1'
    do i = 1, 2
      if (inp%species(i)%z <= 0.0_dp .or. inp%species(i)%nval <= 0) msg = '[[species]] needs z and nval'
    end do
    if (inp%ai%nspin /= 1 .and. inp%ai%nspin /= 2) msg = '[scf] spin must be 1 or 2'
    if (inp%ai%nlevels < 2) msg = '[asr] levels must be >= 2'
    if (inp%ai%nrec /= 0 .and. inp%ai%nrec < inp%ai%nlevels) msg = '[asr] recursion_steps must be 0 or >= levels'
    if (inp%ai%maxrank < 0) msg = '[asr] rank must be >= 0'
    if (inp%ai%symmetry .and. inp%ai%lsuper < 2) msg = '[asr] supercell must be >= 2'
    if (inp%emax <= inp%emin) msg = '[dos] emax must exceed emin'
    if (inp%ai%order /= 2 .and. inp%ai%order /= 3) msg = '[asr] order must be 2 or 3'
    if (pfs /= 'exact' .and. pfs /= 'parametrised') msg = "[asr] potential_function must be 'exact' or 'parametrised'"
    if (chg /= 'lowdin' .and. chg /= 'sphere') msg = "[asr] charges must be 'lowdin' or 'sphere'"
    if (trm /= 'avg2' .and. trm /= 'bp') msg = "[asr] terminator must be 'avg2' or 'bp'"
    if (edt /= 'hull' .and. edt /= 'bp' .and. edt /= 'avg2') msg = "[asr] ed_terminator must be 'hull', 'bp' or 'avg2'"
    if (dns /= 'expand' .and. dns /= 'stuttgart') msg = "[asr] density must be 'expand' or 'stuttgart'"
    if (scr /= 'standard' .and. scr /= 'hardcore' .and. scr /= 'auto') msg = "[asr] screening must be 'auto', 'standard' or 'hardcore'"
    if (hamil /= 'energy-dependent' .and. hamil /= 'energy-independent' .and. hamil /= 'asr_ed') msg = "[asr] hamiltonian must be 'energy-dependent', 'energy-independent' or 'asr_ed'"
    if (inp%ai%energy_dependent .and. .not. inp%ai%ordered .and. .not. inp%ai%symmetry) msg = '[asr] energy-dependent recursion needs symmetry = true'
  end subroutine read_alloy_input

  subroutine run_alloy(inp, ok)
    type(alloy_run_input_t), intent(inout) :: inp
    logical, intent(out) :: ok
    type(alloy_result_t) :: res
    real(dp), allocatable :: egrid(:), dos(:, :, :, :)
    real(dp) :: tsp(2), tot, x
    integer :: i, ie, l, s, sg, u, np
    character(len=:), allocatable :: dfile
    character(len=8) :: syms(2)
    real(dp) :: zc(2), swsc(2), t0, t_scf, t_dos, vat, avw
    integer :: nvc(2), natc(2), ru, nthr
    type(asa_class_result_t) :: rcls(2)

    ok = .false.
    t0 = rep_wall()
    do i = 1, 2
      inp%ai%spec(i)%z = inp%species(i)%z; inp%ai%spec(i)%nval = inp%species(i)%nval; inp%ai%spec(i)%spin_sign = inp%species(i)%spin_sign
      inp%ai%spec(i)%nvl = inp%species(i)%konfig(0:2)
      inp%ai%spec(i)%core = core_shells_for(inp%species(i)%konfig)
      inp%ai%spec(i)%atom_shells = atom_shells(inp%species(i))
    end do
    x = inp%ai%x
    if (inp%ai%ordered) then
      print '(a,a)', 'asr alloy (ordered B2, real-space recursion)  ', inp%name
    else
      print '(a,a)', 'asr alloy  ', inp%name
    end if
    print '(a,a,a,f9.5,a,a,a,f6.3,a,a,f6.3,a,i0)', '  lattice ', trim(inp%ai%lattice), '  a(Bohr)=', inp%ai%alat, '  ', inp%species(1)%symbol, ' x=', x, '  ', &
      inp%species(2)%symbol//' x=', 1.0_dp - x, '  spin=', inp%ai%nspin
    inp%ai%ndos = inp%npoints; inp%ai%dos_emin = inp%emin; inp%ai%dos_emax = inp%emax; inp%ai%dos_eta = inp%eta
    t_scf = rep_wall()
    call alloy_scf(inp%ai, res, verbose=.true.)
    t_scf = rep_wall() - t_scf
    print '(a,l1,a,i0,a,f10.6,a)', '  SCF converged=', res%converged, '  iterations=', res%iter, '  E_F=', res%ef, ' Ry'
    if (.not. res%converged) print '(a)', '  WARNING: SCF did not converge'
    do s = 1, 2
      print '(a,a,a,f9.5,a,f9.5)', '  species ', inp%species(s)%symbol, '  qnet=', res%cls(s)%qnet, '  moment=', res%cls(s)%moment
    end do
    print '(a,i0,a,f9.5)', '  augmented-space nodes: ', res%nodes, '   cell moment=', res%moment_cell

    np = inp%npoints
    allocate (egrid(np), dos(np, 0:2, 2, inp%ai%nspin))
    do ie = 1, np; egrid(ie) = inp%emin + (inp%emax - inp%emin) * (ie - 1) / real(np - 1, dp); end do
    t_dos = rep_wall()
    if (res%ed) then
      dos = res%edos
    else
      call alloy_ldos(res, egrid, inp%eta, dos)
    end if
    dfile = inp%workdir//inp%dos_file
    open (newunit=u, file=dfile, status='replace', action='write')
    write (u, '(a)') '# '//inp%name//'   states/Ry/atom, summed over spin'
    write (u, '(a,f12.8)') '# E_F = ', res%ef
    write (u, '(a)') '# E(Ry)  total  A_s A_p A_d  B_s B_p B_d  total_spin1 total_spin2  (A = '//inp%species(1)%symbol//', B = '// &
                     inp%species(2)%symbol//'; per-spin totals per atom, concentration weighted)'
    do ie = 1, np
      tot = 0.0_dp; tsp = 0.0_dp
      do sg = 1, inp%ai%nspin
        tsp(sg) = x * sum(dos(ie, :, 1, sg)) + (1.0_dp - x) * sum(dos(ie, :, 2, sg))
        tot = tot + tsp(sg)
      end do
      write (u, '(f12.6,9es14.6)') egrid(ie), tot, ((sum(dos(ie, l, s, :)), l=0, 2), s=1, 2), tsp
    end do
    close (u)
    print '(a,a)', '  wrote ', dfile
    t_dos = rep_wall() - t_dos

    ! ---- report (.out)
    vat = inp%ai%alat**3 / merge(4.0_dp, 2.0_dp, trim(inp%ai%lattice) == 'fcc')
    avw = (3.0_dp * vat / (4.0_dp * 3.14159265358979323846_dp))**(1.0_dp / 3.0_dp)
    do i = 1, 2
      syms(i) = inp%species(i)%symbol; zc(i) = inp%species(i)%z; nvc(i) = inp%species(i)%nval
      swsc(i) = merge(inp%ai%spec(i)%sws, avw, inp%ai%spec(i)%sws > 0.0_dp)
    end do
    natc = 1
    nthr = 1
    !$ nthr = omp_get_max_threads()
    call rep_open(ru, inp%workdir//inp%out_file, 'Disordered alloy, augmented-space recursion + TB-LMTO-ASA: '//inp%name, inp%infile, nthr)
    call rep_echo_input(ru, inp%infile)
    call rep_section(ru, 'method')
    call rep_line(ru, 'Hamiltonian   TB-LMTO-ASA, random binary alloy (A = '//trim(syms(1))//', B = '//trim(syms(2))//'), single-site classes')
    if (inp%ai%ordered) then
      call rep_line(ru, 'ordered B2 (CsCl): A on the corner sublattice, B on the body centre; real-space recursion on a periodic bcc supercell,')
      call rep_line(ru, 'no augmented space (nothing is random); exact two-site Madelung shifts')
    else
      write (ru, '(a,f7.4,a,f7.4)') '  concentration  x_A = ', x, '   x_B = ', 1.0_dp - x
    end if
    write (ru, '(a,a,a,i0,a,i0,a,i0)') '  recursion      ', merge(merge('energy-dependent (asr_ed)', 'energy-dependent         ', inp%ai%ed_asr), 'energy-independent       ', res%ed), '   levels=', inp%ai%nlevels, &
      '  max rank=', inp%ai%maxrank, '  supercell=', inp%ai%lsuper
    if (.not. res%ed) write (ru, '(a,i0,a)') '  Hamiltonian    orthonormal ASA expanded to order ', inp%ai%order, ' in h (2: E_nu + h - hoh; 3: + hohoh - (hphh + hhph)/2)'
    write (ru, '(a,a,a,3f9.5,a,f6.3,a,f6.3)') '  screening      ', trim(inp%ai%screening), '  alpha(s,p,d)=', alloy_alpha(inp%ai), &
      '   kept S range=', inp%ai%rcl, '  inversion radius=', max(inp%ai%rcl, inp%ai%rinv), ' (alat)'
    if (inp%ai%ordered) then
      write (ru, '(a,i0,a,i0,a,f6.3)') '  real-space sites ', res%nodes, '   (orbitals per site: 9)   supercell=', inp%ai%lsuper, '   mixing=', inp%ai%beta
    else
      write (ru, '(a,i0,a,l1,a,f6.3,a,f6.3)') '  augmented-space nodes ', res%nodes, '   D4h symmetry=', inp%ai%symmetry, '   SIM beta=', inp%ai%sim_beta, '   mixing=', inp%ai%beta
    end if
    call rep_section(ru, 'structure')
    write (ru, '(a,a)') ' lattice               ', trim(inp%ai%lattice)
    write (ru, '(a,f12.6,a,f12.6,a)') ' lattice constant       ', inp%ai%alat, ' Bohr  = ', inp%ai%alat * 0.529177210903_dp, ' Angstrom'
    write (ru, '(a,f12.5,a,f9.5,a)') ' volume per atom        ', vat * 0.529177210903_dp**3, ' Angstrom^3   S_ws = ', avw, ' Bohr'
    if (inp%ai%ordered) then
      write (ru, '(a)') ' two-atom cell: '//trim(syms(1))//' at (0,0,0), '//trim(syms(2))//' at (1/2,1/2,1/2) (units of the cubic lattice constant)'
    else
      write (ru, '(a)') ' disordered substitutional alloy: no explicit positions (virtual-crystal cell, one site per cell)'
    end if
    call rep_scf_history(ru, syms, res%iter, res%hist_dv, res%hist_ef, res%hist_moment, res%hist_qnet, inp%ai%tol, res%converged)
    ! alloy_scf stores the moments q per spin channel; for spin = 1 the report needs the full (two-spin) occupation
    rcls = res%cls
    if (inp%ai%nspin == 1) then
      do i = 1, 2; rcls(i)%q = 2.0_dp * rcls(i)%q; end do
    end if
    call rep_classes(ru, syms, zc, nvc, swsc, natc, rcls, inp%ai%nspin, res%ef, res%moment_cell, ' (per species, not concentration-weighted)')
    write (ru, '(a)') ' species are alloy components; qnet/moments above are per atom of that species'
    call rep_dos_at_ef(ru, syms, egrid, dos, res%ef, inp%eta, inp%ai%nspin, &
                       'per-species LDOS from the recursion (complex-energy broadening eta shown as broadening); weight by concentration for the alloy DOS')
    call rep_section(ru, 'timing')
    write (ru, '(a,f10.2,a)') ' SCF  ', t_scf, ' s'
    write (ru, '(a,f10.2,a)') ' DOS  ', t_dos, ' s'
    write (ru, '(a)') ' files written:'
    write (ru, '(a)') '   '//dfile
    call rep_close(ru, res%converged, rep_wall() - t0)
    print '(a,a)', '  wrote ', inp%workdir//inp%out_file
    ok = res%converged
  end subroutine run_alloy
end module asr_alloy_driver
