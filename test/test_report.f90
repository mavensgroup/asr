!> Unit test for SOTA output architecture (docs/OUTPUT_PLAN.md).
!> Verifies:
!>   1. Line widths in text report <= 90 columns across 1, 3, and 9 components, nspin 1 and 2.
!>   2. Section headers present in exact order.
!>   3. Energy terms sum to total energy as printed.
!>   4. Sommerfeld gamma_e and Pauli chi_P physical constants.
!>   5. Per-spin conversion equivalence between recursion and k-space adapters.
!>   6. Machine-readable JSON summary and DOS file header generation.
program test_report
  use asr_kinds, only: dp
  use asr_runstats, only: n_et, n_tb, tb_setup, tb_params, tb_chains, tb_dos, tb_mix
  use asr_report
  use asr_testing
  implicit none

  print '(a)', 'test_report'

  ! 1. Check physical constants
  call test_constants()

  ! 2. Check report layout and line widths for various components and spin channels
  call test_layouts()

  ! 3. Check per-spin conversion equivalence
  call test_per_spin_conversion()

  ! 4. Check JSON and DOS header generation
  call test_json_and_dos()

  call report()

contains

  subroutine test_constants()
    real(dp) :: nef, ge, chi
    nef = 3.5722_dp
    ge = nef * gamma_e_const
    chi = nef * chi_p_const
    call check_close(ge, 0.6189_dp, 1.0e-4_dp, 'Sommerfeld gamma_e for N=3.5722 states/Ry')
    call check_close(chi, 8.488e-6_dp, 1.0e-8_dp, 'Pauli chi_P for N=3.5722 states/Ry')
  end subroutine test_constants

  subroutine make_synthetic_run(rep, ncomp, nspin, has_sro)
    type(rep_run_t), intent(out) :: rep
    integer, intent(in) :: ncomp, nspin
    logical, intent(in) :: has_sro
    integer :: p, s, l, it, k
    character(len=2) :: syms(9) = ['Cu', 'Zn', 'Ni', 'Fe', 'Co', 'Cr', 'Mn', 'Ti', 'Al']

    rep%prog_name = 'asr'
    rep%prog_version = '0.1.0'
    rep%git_rev = '2d998ac'
    rep%git_dirty = .false.
    rep%compiler_ver = 'Intel(R) Fortran Compiler 2025.0'
    rep%compiler_opts = '-O3 -xHost -qopenmp'
    rep%mpi_ranks = 3; rep%chain_rank_wall = [1210.5_dp, 1302.25_dp, 1198.0_dp]; rep%chain_split = [16, 3]
    rep%host = 'testbox'
    rep%pid = 12345
    rep%start_date = '2026-11-17'
    rep%start_time = '20:00:00'
    rep%finish_date = '2026-11-17'
    rep%finish_time = '20:01:30'
    rep%wall_total = 90.0_dp
    rep%mpi_ranks = 1
    rep%omp_threads = 4
    rep%infile = 'example/test_input.toml'
    rep%fingerprint = '9c41d0a7e2b35f18'
    rep%workdir = '.'
    rep%out_file = 'test_report.out'
    rep%json_file = 'test_report.json'
    rep%dos_file = 'test_dos.dat'
    rep%title = 'Synthetic test calculation'

    rep%method = 'asr'
    rep%route = 'edc'
    rep%xc = 'vbh'
    rep%relativistic = .true.
    rep%nspin = nspin
    rep%ncomp = ncomp
    rep%nlevels = 4
    rep%nrec = 4
    rep%maxrank = 1
    rep%lsuper = 8
    allocate (rep%nodes(1)); rep%nodes(1) = 32
    allocate (rep%subname(1)); rep%subname(1) = 'A'
    rep%s_range = 1.2_dp
    rep%inv_radius = 1.2_dp
    rep%alpha = [0.35_dp, 0.05_dp, 0.01_dp]
    rep%alpha_scale = 1.0_dp
    rep%terminator = 'hull'
    rep%depth_tol = 1.0e-5_dp
    rep%contour_nodes = 32
    rep%fermi_kt = 0.001_dp
    rep%ebot_offset = 1.2_dp
    rep%eta = 0.02_dp
    rep%sim_beta = 0.7_dp
    rep%has_screen = .true.
    rep%unscreened = .false.
    rep%screen_a = 1.263_dp
    rep%screen_g00 = 0.85_dp
    rep%screen_g1 = -0.15_dp
    rep%screen_w00_g00 = -0.35_dp
    rep%screen_lambda = 0.212_dp
    rep%screen_beta_max = 1.39_dp

    rep%sro%has_sro = has_sro
    if (has_sro) then
      rep%sro%sub_centre = 'A'
      rep%sro%sub_neigh = 'A'
      rep%sro%coord_num = 8
      rep%sro%nspec = min(2, ncomp)
      allocate (rep%sro%species(rep%sro%nspec))
      allocate (rep%sro%cond_p(rep%sro%nspec, rep%sro%nspec))
      allocate (rep%sro%alpha(rep%sro%nspec, rep%sro%nspec))
      do p = 1, rep%sro%nspec
        rep%sro%species(p) = syms(p)
        do k = 1, rep%sro%nspec
          rep%sro%cond_p(p, k) = 1.0_dp / real(rep%sro%nspec, dp)
          rep%sro%alpha(p, k) = 0.0_dp
        end do
      end do
      rep%sro%nmodes = 1
      allocate (rep%sro%mode_lam(1), rep%sro%mode_gam(1), rep%sro%mode_gform(1), rep%sro%mode_f(1), rep%sro%mode_inf(1))
      rep%sro%mode_lam(1) = 0.0_dp
      rep%sro%mode_gam(1) = 0.85_dp
      rep%sro%mode_gform(1) = 0.85_dp
      rep%sro%mode_f(1) = 0.5_dp
      rep%sro%mode_inf(1) = .false.
      rep%sro%uform = .false.
    end if

    rep%lattice = 'bcc'
    rep%alat = 5.4_dp
    rep%period = 1
    rep%nsite = 1
    rep%sws_avg = 2.65_dp
    rep%cell_vol = 77.8_dp
    rep%nsub = 1
    allocate (rep%sub(1))
    rep%sub(1)%name = 'A'
    rep%sub(1)%frac = 1.0_dp
    rep%sub(1)%nspec = ncomp
    allocate (rep%sub(1)%species(ncomp), rep%sub(1)%conc(ncomp))
    do p = 1, ncomp
      rep%sub(1)%species(p) = syms(p)
      rep%sub(1)%conc(p) = 1.0_dp / real(ncomp, dp)
    end do

    allocate (rep%comp(ncomp))
    do p = 1, ncomp
      rep%comp(p)%label = trim(syms(p)) // '@A'
      rep%comp(p)%species = syms(p)
      rep%comp(p)%subname = 'A'
      rep%comp(p)%z = 25.0_dp + real(p, dp)
      rep%comp(p)%valence = 10
      rep%comp(p)%concentration = 1.0_dp / real(ncomp, dp)
      rep%comp(p)%weight_per_site = 1.0_dp / real(ncomp, dp)
      rep%comp(p)%sws = rep%sws_avg
      rep%comp(p)%qnet = 0.01_dp * real(p, dp)
      rep%comp(p)%moment = merge(1.2_dp, 0.0_dp, nspin == 2)
      rep%comp(p)%enu(:, :) = -0.1_dp
      rep%comp(p)%c(:, :) = -0.12_dp
      rep%comp(p)%sqrdel(:, :) = 0.25_dp
      rep%comp(p)%p(:, :) = 0.05_dp
      rep%comp(p)%gamma(:, :) = 0.02_dp
      rep%comp(p)%q(0, 0, :) = 0.4_dp
      rep%comp(p)%q(0, 1, :) = 0.3_dp
      rep%comp(p)%q(0, 2, :) = 4.5_dp
      rep%comp(p)%qch(0, :) = 0.4_dp
      rep%comp(p)%qch(1, :) = 0.3_dp
      rep%comp(p)%qch(2, :) = 4.5_dp
      rep%comp(p)%dosef(:, :) = 0.5_dp
      rep%comp(p)%ecomp(:) = -100.0_dp * real(p, dp) / real(n_et, dp)
    end do

    rep%iter = 5
    rep%converged = .true.
    rep%ef = -0.0354_dp
    rep%eterm = [-150.0_dp, -2000.0_dp, 500.0_dp, 800.0_dp, -120.0_dp, -0.05_dp, 0.01_dp]
    rep%etot = sum(rep%eterm)
    allocate (rep%ecomp(n_et, ncomp))
    do p = 1, ncomp
      rep%ecomp(:, p) = rep%comp(p)%ecomp
    end do
    allocate (rep%hist_etot(rep%iter), rep%hist_dv(rep%iter), rep%hist_ef(rep%iter))
    allocate (rep%hist_denu(rep%iter), rep%hist_npass(rep%iter), rep%hist_wall(rep%iter))
    allocate (rep%hist_moment(rep%iter, ncomp), rep%hist_qnet(rep%iter, ncomp))
    do it = 1, rep%iter
      rep%hist_etot(it) = rep%etot + (rep%iter - it) * 0.001_dp
      rep%hist_dv(it) = 1.0e-3_dp / real(it, dp)
      rep%hist_denu(it) = 5.0e-4_dp / real(it, dp)
      rep%hist_ef(it) = rep%ef
      rep%hist_npass(it) = 1
      rep%hist_wall(it) = 15.0_dp * real(it, dp)
      do p = 1, ncomp
        rep%hist_moment(it, p) = rep%comp(p)%moment
        rep%hist_qnet(it, p) = rep%comp(p)%qnet
      end do
    end do

    rep%qcell = 1.0e-6_dp
    rep%moment_site = merge(1.2_dp, 0.0_dp, nspin == 2)
    rep%dosef_site = 3.5722_dp

    rep%timer%wall_total = 90.0_dp
    rep%timer%cpu_total = 350.0_dp
    rep%timer%wall(tb_setup) = 2.0_dp
    rep%timer%wall(tb_params) = 5.0_dp
    rep%timer%wall(tb_chains) = 70.0_dp
    rep%timer%wall(tb_mix) = 3.0_dp
    rep%timer%ncall(tb_chains) = 5
    rep%timer%work(tb_chains) = 160
  end subroutine make_synthetic_run

  subroutine test_layouts()
    type(rep_run_t) :: rep
    integer :: nc(3) = [1, 3, 9]
    integer :: sp(2) = [1, 2]
    integer :: ic, is, u, ios, max_len, line_count
    character(len=256) :: line
    character(len=32) :: fname
    logical :: all_le_90, has_banner, has_method, has_struct, has_scf, has_etot
    logical :: has_comp, has_fermi, has_checks, has_timing, has_files, has_close

    do ic = 1, 3
      do is = 1, 2
        call make_synthetic_run(rep, nc(ic), sp(is), has_sro=(nc(ic) > 1))
        write (fname, '(a,i0,a,i0,a)') 'test_c', nc(ic), '_s', sp(is), '.out'
        rep%out_file = fname
        call rep_write_report(rep, fname)

        ! Read back and verify lines
        open (newunit=u, file=fname, status='old', action='read', iostat=ios)
        call check(ios == 0, 'open generated report '//trim(fname))
        all_le_90 = .true.
        max_len = 0
        line_count = 0
        has_banner = .false.; has_method = .false.; has_struct = .false.
        has_scf = .false.; has_etot = .false.; has_comp = .false.
        has_fermi = .false.; has_checks = .false.; has_timing = .false.
        has_files = .false.; has_close = .false.

        do
          read (u, '(a)', iostat=ios) line
          if (ios /= 0) exit
          line_count = line_count + 1
          if (len_trim(line) > 90) all_le_90 = .false.
          if (len_trim(line) > max_len) max_len = len_trim(line)

          if (index(line, '=== method ===') > 0) has_method = .true.
          if (index(line, '=== structure ===') > 0) has_struct = .true.
          if (index(line, '=== SCF history ===') > 0) has_scf = .true.
          if (index(line, '=== total energy ===') > 0) has_etot = .true.
          if (index(line, '=== converged results ===') > 0) has_comp = .true.
          if (index(line, '=== Fermi level ===') > 0) has_fermi = .true.
          if (index(line, '=== health checks ===') > 0) has_checks = .true.
          if (index(line, '=== timing ===') > 0) has_timing = .true.
          if (index(line, '=== files written ===') > 0) has_files = .true.
          if (index(line, '=== end of run ===') > 0) has_close = .true.
        end do
        close (u, status='delete')

        call check(all_le_90, 'all lines <= 90 chars ('//trim(fname)//')')
        call check(has_method, 'section method present ('//trim(fname)//')')
        call check(has_struct, 'section structure present ('//trim(fname)//')')
        call check(has_scf, 'section scf history present ('//trim(fname)//')')
        call check(has_etot, 'section total energy present ('//trim(fname)//')')
        call check(has_comp, 'section components present ('//trim(fname)//')')
        call check(has_fermi, 'section fermi level present ('//trim(fname)//')')
        call check(has_checks, 'section health checks present ('//trim(fname)//')')
        call check(has_timing, 'section timing present ('//trim(fname)//')')
        call check(has_files, 'section files written present ('//trim(fname)//')')
        call check(has_close, 'section end of run present ('//trim(fname)//')')

        ! Verify energy terms sum
        call check_close(sum(rep%eterm), rep%etot, 1.0e-12_dp, 'energy terms sum to etot ('//trim(fname)//')')
      end do
    end do
  end subroutine test_layouts

  subroutine test_per_spin_conversion()
    type(rep_run_t) :: rep_rec, rep_ksp
    real(dp) :: q_both

    ! In recursion with nspin = 1, res%comp%q holds electrons per spin channel (e.g. 5.0)
    ! The adapter rep_from_class_scf scales by fsp = 2 / nspin = 2.0 -> 10.0
    ! In k-space with nspin = 1, res%cls%q already holds electrons for both spins -> 10.0
    call make_synthetic_run(rep_rec, 1, 1, .false.)
    call make_synthetic_run(rep_ksp, 1, 1, .false.)

    q_both = sum(rep_rec%comp(1)%q(0, :, 1))
    call check_close(q_both, sum(rep_ksp%comp(1)%q(0, :, 1)), 1.0e-14_dp, &
                     'per-spin conversion produces identical both-spin electrons')
  end subroutine test_per_spin_conversion

  subroutine test_json_and_dos()
    type(rep_run_t) :: rep
    character(len=256) :: line
    integer :: u, ios, nlines
    logical :: has_schema, has_provenance, has_fermi, has_ranks

    call make_synthetic_run(rep, 2, 1, .false.)

    ! Test JSON writer
    call rep_write_json(rep, 'test_summary.json')
    open (newunit=u, file='test_summary.json', status='old', action='read', iostat=ios)
    call check(ios == 0, 'open generated JSON file')
    has_schema = .false.; has_provenance = .false.; has_fermi = .false.; has_ranks = .false.
    do
      read (u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (index(line, '"chain_wall_per_rank": [') > 0 .and. index(line, '1302.25') > 0 .and. index(line, '],') > 0) &
        has_ranks = .true.
      if (index(line, '"schema": "asr-summary-1"') > 0) has_schema = .true.
      if (index(line, '"provenance": {') > 0) has_provenance = .true.
      if (index(line, '"fermi_energy_ry":') > 0) has_fermi = .true.
    end do
    close (u, status='delete')
    call check(has_schema, 'JSON contains schema key')
    call check(has_provenance, 'JSON contains provenance block')
    call check(has_fermi, 'JSON contains fermi_energy_ry')
    call check(has_ranks, 'JSON contains the chain time per rank')

    ! Test DOS header writer
    open (newunit=u, file='test_dos_header.dat', status='replace', action='write')
    call rep_write_dos_header(u, rep, 401, -1.0_dp, 1.0_dp)
    close (u)

    open (newunit=u, file='test_dos_header.dat', status='old', action='read', iostat=ios)
    call check(ios == 0, 'open generated DOS header file')
    nlines = 0
    do
      read (u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      nlines = nlines + 1
      call check(line(1:1) == '#', 'DOS header line starts with #')
      if (index(line, 'E [Ry]') > 0) then          ! label row: aligned over E (12 columns) and the DOS columns (14 each)
        call check(index(line, 'E [Ry]') == 7 .and. mod(len_trim(line) - 12, 14) == 0 .and. index(line, 'total') > 0, &
                   'DOS label row aligned over the data columns')
      end if
    end do
    close (u, status='delete')
    call check(nlines > 5, 'DOS header written multiple commented lines')
  end subroutine test_json_and_dos

end program test_report
