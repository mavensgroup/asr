!> Human-readable run report (asr.out) and machine-readable summary (asr.json).
!> Implements the SOTA output architecture described in docs/OUTPUT_PLAN.md:
!>   - Single canonical unit system (Ry, Bohr, e, mu_B, states/Ry), units in every header
!>   - Maximum 90 columns width, guaranteed no number collisions or wrapping
!>   - Fixed greppable section headers ("=== name ===")
!>   - Unified SCF history table, converged energetics decomposition, component tables
!>   - Fermi level card with Sommerfeld gamma_e and Pauli chi_P
!>   - Health checks table with status tags (OK, WARN, FAIL, info)
!>   - Hierarchical timing profile
!>   - Plain JSON summary (<out>.json)
!>   - Standardized self-describing DOS file headers
module asr_report
  use asr_kinds, only: dp, i8, pi
  use asr_runstats, only: run_timer_t, n_tb, tb_name, tb_setup, tb_params, tb_chains, tb_below, tb_dos, &
                          tb_bands, tb_density, tb_mix, tb_io, n_et, et_name, et_band, et_core, et_dc, &
                          et_hartree, et_xc, et_mad, et_sim, wall_seconds
  use asr_build_info, only: git_rev, git_dirty
  use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options
  use, intrinsic :: iso_c_binding, only: c_int
  implicit none
  private

  public :: rep_run_t, rep_component_t, rep_sublattice_t, rep_sro_t, rep_site_t
  public :: rep_write_report, rep_write_json, rep_write_dos_header, rep_wall, fnv1a64_file
  public :: rep_banner, rep_echo_input, rep_method, rep_structure, rep_sro, rep_scf_history
  public :: rep_energy, rep_components, rep_fermi, rep_checks, rep_timing, rep_warnings, rep_files, rep_close
  public :: ry_to_ev, bohr_to_ang, gamma_e_const, chi_p_const, get_host, get_pid, collect_warnings, clean_flags, fmt_dt

  character(len=*), parameter :: program_version = 'asr 0.1.0'
  real(dp), parameter :: ry_to_ev = 13.605693123_dp
  real(dp), parameter :: bohr_to_ang = 0.529177211_dp
  real(dp), parameter :: gamma_e_const = 0.173247_dp      !! mJ/(mol K^2) per (state/Ry, both spins)
  real(dp), parameter :: chi_p_const   = 2.37605e-6_dp    !! emu/mol per (state/Ry, both spins)

  interface
    function c_getpid() bind(c, name='getpid') result(res)
      import :: c_int
      integer(c_int) :: res
    end function c_getpid
  end interface

  type :: rep_component_t
    character(len=16) :: label = ''        !! e.g. 'Cu:A' or 'Ag'
    character(len=8)  :: species = ''      !! e.g. 'Cu'
    character(len=8)  :: subname = ''      !! e.g. 'A'
    real(dp) :: z = 0.0_dp
    integer  :: valence = 0
    real(dp) :: concentration = 1.0_dp
    real(dp) :: weight_per_site = 1.0_dp
    real(dp) :: sws = 0.0_dp               !! sphere radius in Bohr
    real(dp) :: qnet = 0.0_dp              !! electrons in sphere - Z
    real(dp) :: moment = 0.0_dp            !! spin moment in sphere (mu_B)
    real(dp) :: enu(0:2, 2) = 0.0_dp       !! linearisation energies (Ry)
    real(dp) :: c(0:2, 2) = 0.0_dp         !! potential parameters (Ry)
    real(dp) :: sqrdel(0:2, 2) = 0.0_dp
    real(dp) :: p(0:2, 2) = 0.0_dp
    real(dp) :: gamma(0:2, 2) = 0.0_dp
    real(dp) :: q(0:2, 0:2, 2) = 0.0_dp    !! moments q0, q1, q2 (both-spin electrons)
    real(dp) :: qch(0:2, 2) = 0.0_dp       !! channel sphere charge Q_l (both-spin electrons)
    real(dp) :: dosef(0:2, 2) = 0.0_dp     !! DOS at E_F (states/Ry, both spins)
    real(dp) :: ecomp(n_et) = 0.0_dp       !! component energy terms (Ry/atom)
  end type rep_component_t

  type :: rep_sublattice_t
    character(len=16) :: name = ''
    real(dp) :: frac = 1.0_dp
    integer :: nspec = 0
    character(len=8), allocatable :: species(:)
    real(dp), allocatable :: conc(:)
  end type rep_sublattice_t

  type :: rep_site_t
    integer :: ia = 0
    character(len=8) :: species = ''
    integer :: iclass = 0
    real(dp) :: pos(3) = 0.0_dp
    real(dp) :: sws = 0.0_dp
  end type rep_site_t

  type :: rep_sro_t
    logical :: has_sro = .false.
    character(len=16) :: sub_centre = '', sub_neigh = ''
    integer :: coord_num = 0               !! e.g. 8 for bcc, 12 for fcc
    integer :: nspec = 0
    character(len=8), allocatable :: species(:)
    real(dp), allocatable :: cond_p(:, :)  !! P(j|i)
    real(dp), allocatable :: alpha(:, :)   !! alpha_ij = 1 - P(j|i)/c_j
    integer :: nmodes = 0
    real(dp), allocatable :: mode_lam(:)   !! eigenvalue of Ahat
    real(dp), allocatable :: mode_gam(:)   !! Gamma = dq/db
    real(dp), allocatable :: mode_gform(:) !! G00 + G1 lambda
    real(dp), allocatable :: mode_f(:)     !! F (Ry/e)
    logical, allocatable :: mode_inf(:)    !! infinite flag (Gamma = 0)
    logical :: uform = .false.
    real(dp), allocatable :: udrive(:)     !! drive u per component (e)
  end type rep_sro_t

  type :: rep_run_t
    ! Provenance and execution metadata
    character(len=64) :: prog_name = 'asr'
    character(len=32) :: prog_version = '0.1.0'
    character(len=64) :: git_rev = 'unknown'
    logical :: git_dirty = .false.
    character(len=128) :: compiler_ver = ''
    character(len=4096) :: compiler_opts = ''
    character(len=64) :: host = ''
    integer :: pid = 0
    character(len=10) :: start_date = '', start_time = ''
    character(len=10) :: finish_date = '', finish_time = ''
    real(dp) :: wall_total = 0.0_dp
    integer :: mpi_ranks = 1
    integer :: omp_threads = 1
    character(len=256) :: infile = ''
    character(len=16) :: fingerprint = '0000000000000000'
    character(len=256) :: workdir = ''
    character(len=256) :: out_file = 'asr.out'
    character(len=256) :: json_file = 'asr.json'
    character(len=256) :: dos_file = ''
    character(len=256) :: pot_file = ''
    character(len=256) :: params_file = ''
    character(len=128) :: title = ''

    ! Method and model
    character(len=16) :: method = 'asr'       !! 'asr' or 'kspace'
    character(len=16) :: route = 'edc'        !! 'edc', 'ed3', or 'kspace'
    character(len=16) :: xc = 'vbh'
    logical :: relativistic = .true.
    integer :: nspin = 1
    integer :: ncomp = 0
    type(rep_component_t), allocatable :: comp(:)

    ! Recursion route specific knobs
    integer :: nlevels = 0
    integer :: nrec = 0
    integer :: maxrank = 0
    integer :: lsuper = 0
    integer, allocatable :: nodes(:)
    character(len=16), allocatable :: subname(:)
    real(dp) :: s_range = 1.2_dp, inv_radius = 1.2_dp
    real(dp) :: alpha(0:2) = 0.0_dp
    real(dp) :: alpha_scale = 1.0_dp
    character(len=16) :: terminator = 'hull'
    real(dp) :: depth_tol = 1.0e-5_dp
    integer :: contour_nodes = 0            !! all nodes of the occupation contour
    integer :: contour_main = 0             !! nodes on the main line
    integer :: contour_poles = 0            !! Matsubara poles (0: closed contour without Fermi smearing)
    real(dp) :: fermi_kt = 0.0_dp
    integer :: fermi_poles = 0
    real(dp) :: ebot_offset = 1.2_dp
    real(dp) :: eta = 0.02_dp               !! [dos] eta: broadening of the DOS curve
    real(dp) :: ed_eta = 0.01_dp            !! Im z of the E_F probe: broadening of N(E_F)
    real(dp) :: sim_beta = 0.7_dp
    logical :: has_screen = .false.
    logical :: unscreened = .false.
    real(dp) :: screen_a = 0.0_dp, screen_g00 = 0.0_dp, screen_g1 = 0.0_dp
    real(dp) :: screen_w00_g00 = 0.0_dp, screen_lambda = 0.0_dp, screen_beta_max = 0.0_dp
    type(rep_sro_t) :: sro

    ! k-space route specific knobs
    integer :: k_nk = 0
    character(len=16) :: k_bz = 'tetra'
    real(dp) :: k_smear = 0.01_dp
    integer :: k_dos_nk = 0
    real(dp) :: k_sigma = 0.7_dp

    ! Mixing and convergence criteria
    character(len=16) :: mix_type = 'Anderson'
    real(dp) :: mix_beta = 0.3_dp
    integer :: mix_history = 6
    real(dp) :: mix_restart = 0.3_dp
    real(dp) :: dvmax = 0.15_dp
    real(dp) :: denumax = 0.08_dp
    real(dp) :: tol = 1.0e-5_dp

    ! Structure
    character(len=16) :: lattice = ''
    real(dp) :: alat = 0.0_dp             !! in Bohr
    real(dp) :: avec(3, 3) = 0.0_dp
    real(dp) :: cell_vol = 0.0_dp         !! Bohr^3
    integer :: period = 1
    integer :: nsite = 1
    real(dp) :: sws_avg = 0.0_dp
    logical :: is_space_filling = .true.
    integer :: nsub = 0
    type(rep_sublattice_t), allocatable :: sub(:)
    integer :: nat = 0
    type(rep_site_t), allocatable :: sites(:)

    ! SCF history
    integer :: iter = 0
    logical :: converged = .false.
    real(dp), allocatable :: hist_etot(:)
    real(dp), allocatable :: hist_dv(:)
    real(dp), allocatable :: hist_denu(:)
    real(dp), allocatable :: hist_ef(:)
    integer, allocatable :: hist_npass(:)
    real(dp), allocatable :: hist_wall(:)
    real(dp), allocatable :: hist_moment(:, :)   !! (iter, ncomp)
    real(dp), allocatable :: hist_qnet(:, :)     !! (iter, ncomp)

    ! Energy breakdown
    real(dp) :: etot = 0.0_dp
    real(dp) :: eterm(n_et) = 0.0_dp
    real(dp), allocatable :: ecomp(:, :)        !! (n_et, ncomp)

    ! Converged observables
    real(dp) :: ef = 0.0_dp
    real(dp) :: qcell = 0.0_dp
    real(dp) :: moment_site = 0.0_dp
    real(dp) :: dosef_site = 0.0_dp

    ! Diagnostics
    integer :: nclamp = 0, nclamp_run = 0
    integer :: nnoncausal = 0, nnode = 0
    real(dp) :: wbelow1 = 0.0_dp, wbelow = 0.0_dp, wbelow_lim = 1.0e-3_dp
    logical :: below_ok = .true.
    real(dp) :: depth_mean = 0.0_dp
    integer :: nfall = 0
    integer :: mix_nrestart = 0, mix_nh = 0
    real(dp) :: eband = 0.0_dp, sumev = 0.0_dp   !! k-space band energy check

    ! Timing
    type(run_timer_t) :: timer

    ! Collected warnings
    character(len=128), allocatable :: warnings(:)
  end type rep_run_t

contains

  !> Wall-clock seconds since arbitrary origin (differences only).
  function rep_wall() result(t)
    real(dp) :: t
    t = wall_seconds()
  end function rep_wall

  !> Read host name from /proc/sys/kernel/hostname or environment
  function get_host() result(h)
    character(len=64) :: h
    integer :: u, ios
    h = ''
    open (newunit=u, file='/proc/sys/kernel/hostname', status='old', action='read', iostat=ios)
    if (ios == 0) then
      read (u, '(a)', iostat=ios) h
      close (u)
    end if
    if (len_trim(h) == 0) call get_environment_variable('HOSTNAME', h)
    if (len_trim(h) == 0) h = 'localhost'
    h = trim(adjustl(h))
  end function get_host

  !> Process ID via POSIX C getpid()
  function get_pid() result(p)
    integer :: p
    p = int(c_getpid())
  end function get_pid

  !> 64-bit FNV-1a hash of a file as a 16-character lowercase hex string.
  function fnv1a64_file(path) result(hex)
    character(len=*), intent(in) :: path
    character(len=16) :: hex
    integer(i8), parameter :: fnv_offset = -3750763034362895579_i8  !! z'CBF29CE484222325'
    integer(i8), parameter :: fnv_prime  = 1099511628211_i8          !! z'00000100000001B3'
    integer(i8) :: h
    integer :: u, ios, i
    character(len=1) :: b
    h = fnv_offset
    open (newunit=u, file=path, status='old', access='stream', form='unformatted', action='read', iostat=ios)
    if (ios /= 0) then
      hex = '0000000000000000'
      return
    end if
    do
      read (u, iostat=ios) b
      if (ios /= 0) exit
      h = ieor(h, int(iachar(b), i8))
      h = h * fnv_prime
    end do
    close (u)
    write (hex, '(z16.16)') h
    do i = 1, 16
      if (hex(i:i) >= 'A' .and. hex(i:i) <= 'Z') hex(i:i) = char(iachar(hex(i:i)) + 32)
    end do
  end function fnv1a64_file

  !> Write a line trimmed, ensuring length at most 90 columns.
  subroutine wline(u, text)
    integer, intent(in) :: u
    character(len=*), intent(in) :: text
    if (len_trim(text) > 90) then
      write (u, '(a)') text(1:90)
    else
      write (u, '(a)') trim(text)
    end if
  end subroutine wline

  !> Write fixed section header with width exactly 90.
  subroutine rep_section(u, name)
    integer, intent(in) :: u
    character(len=*), intent(in) :: name
    character(len=90) :: s
    integer :: n
    s = '=== ' // trim(name) // ' ==='
    n = len_trim(s)
    if (n < 90) s(n+1:90) = repeat('=', 90 - n)
    call wline(u, s)
  end subroutine rep_section

  !> Wrap tokens to fit within 90 columns with a given indentation.
  subroutine write_wrapped(u, prefix, text, indent)
    integer, intent(in) :: u
    character(len=*), intent(in) :: prefix, text, indent
    character(len=512) :: line, token
    integer :: p, q, n
    line = trim(prefix)
    n = len_trim(text)
    if (n == 0) then
      call wline(u, trim(line) // ' (none)')
      return
    end if
    p = 1
    do while (p <= n)
      do while (p <= n .and. text(p:p) == ' '); p = p + 1; end do
      if (p > n) exit
      q = p
      do while (q <= n .and. text(q:q) /= ' '); q = q + 1; end do
      token = text(p:q-1)
      p = q
      if (len_trim(line) > len_trim(prefix) .and. len_trim(line) + 1 + len_trim(token) > 90) then
        call wline(u, line)
        line = indent // trim(token)
      else
        if (len_trim(line) > 0) then
          line = trim(line) // ' ' // trim(token)
        else
          line = trim(token)
        end if
      end if
    end do
    if (len_trim(line) > 0) call wline(u, line)
  end subroutine write_wrapped


  !> Section 4.1: Banner and provenance
  subroutine rep_banner(u, r)
    integer, intent(in) :: u
    type(rep_run_t), intent(in) :: r
    character(len=90) :: line
    character(len=32) :: git_status
    character(len=24), parameter :: ind = '                      '

    if (r%git_dirty) then; git_status = '(working tree dirty)'; else; git_status = '(working tree clean)'; end if
    call wline(u, repeat('=', 90))
    call wline(u, ' ASR  augmented-space recursion for disordered alloys  *  TB-LMTO-ASA')
    call wline(u, repeat('=', 90))
    call write_wrapped(u, ' program            : ', trim(program_version) // '   git ' // trim(r%git_rev) // '   ' // &
                       trim(git_status), ind)
    call write_wrapped(u, ' compiler           : ', trim(r%compiler_ver), ind)
    call write_wrapped(u, ' flags              : ', trim(clean_flags(r%compiler_opts)), ind)
    call write_wrapped(u, ' started            : ', trim(fmt_dt(r%start_date, r%start_time)) // '   host ' // trim(r%host) // &
                       '   pid ' // trim(itoa(r%pid)), ind)
    write (line, '(a,i0,a,i0,a)') ' parallel           : MPI ranks ', r%mpi_ranks, ', OpenMP threads ', r%omp_threads, ' per rank'
    call wline(u, line)
    call write_wrapped(u, ' input              : ', trim(r%infile) // '   fingerprint ' // trim(r%fingerprint), ind)
    call write_wrapped(u, ' working directory  : ', trim(abs_dir(r%workdir)), ind)
    call wline(u, ' units              : energy Ry, length Bohr, charge e, moment mu_B, DOS states/Ry')
    call wline(u, ' conversions        : 1 Ry = 13.605693123 eV = 0.5 Ha;  1 Bohr = 0.529177211 Angstrom')
    call wline(u, repeat('=', 90))
  end subroutine rep_banner

  !> Section 4.2: Echo of verbatim input lines with prefix '   | '
  subroutine rep_echo_input(u, infile)
    integer, intent(in) :: u
    character(len=*), intent(in) :: infile
    character(len=512) :: line
    integer :: v, ios
    call rep_section(u, 'input (verbatim)')
    open (newunit=v, file=infile, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      call wline(u, '   | (input file not readable)')
      return
    end if
    do
      read (v, '(a)', iostat=ios) line
      if (ios /= 0) exit
      call wline(u, '   | '//trim(line))
    end do
    close (v)
  end subroutine rep_echo_input


  !> Section 4.3: Method specifications
  subroutine rep_method(u, r)
    integer, intent(in) :: u
    type(rep_run_t), intent(in) :: r
    character(len=90) :: line
    integer :: i

    call rep_section(u, 'method')
    if (trim(r%method) == 'kspace') then
      call wline(u, ' Hamiltonian        : TB-LMTO-ASA, k-space route (internal check)')
      call wline(u, ' basis              : s p d, scalar relativistic; parameters in the gamma representation')
      write (line, '(a,a,a,a)') ' XC functional      : ', trim(r%xc), ' (libxc), LDA;  spin polarised: ', merge('yes', 'no ', r%nspin == 2)
      call wline(u, line)
      write (line, '(a,f6.3)') ' screening sigma    : ', r%k_sigma
      call wline(u, line)
      write (line, '(a,a,a,i0,a,f7.4,a,i0)') ' BZ integration     : ', trim(r%k_bz), ', nk ', r%k_nk, ', smearing ', r%k_smear, &
        ' Ry, DOS nk ', r%k_dos_nk
      call wline(u, line)
    else
      if (trim(r%route) == 'ed3') then
        call wline(u, ' Hamiltonian        : TB-LMTO-ASA, full generalised ASA problem (route ed3)')
      else
        call wline(u, ' Hamiltonian        : TB-LMTO-ASA, energy-dependent ASR, route edc (H^(2) form)')
      end if
      call wline(u, ' basis              : s p d, scalar relativistic; parameters in the gamma representation')
      write (line, '(a,a,a,a)') ' XC functional      : ', trim(r%xc), ' (libxc), LDA;  spin polarised: ', merge('yes', 'no ', r%nspin == 2)
      call wline(u, line)
      write (line, '(a,i0,a)') ' components         : ', r%ncomp, ' (sublattice, species) pairs, one potential each'
      call wline(u, line)
      write (line, '(a,i0,a,i0,a,i0,a,i0)') ' augmented space    : levels ', r%nlevels, ', recursion steps ', r%nrec, &
        ', max rank ', r%maxrank, ', supercell ', r%lsuper
      call wline(u, line)
      if (allocated(r%nodes) .and. allocated(r%subname)) then
        line = ' nodes / sublattice :'
        do i = 1, size(r%nodes)
          line = trim(line) // ' ' // trim(r%subname(i)) // ' ' // trim(itoa(r%nodes(i)))
          if (i < size(r%nodes)) line = trim(line) // ','
        end do
        call wline(u, line)
      end if
      write (line, '(a,a,a,a,a)') ' S^alpha            : kept range ', trim(ftoa(r%s_range, 3)), ' alat, inversion radius ', &
        trim(ftoa(r%inv_radius, 3)), ' alat'
      call wline(u, line)
      write (line, '(a,f8.5,a,f8.5,a,f8.5)') ' screening alpha    : s', r%alpha(0), '   p', r%alpha(1), '   d', r%alpha(2)
      call wline(u, line)
      write (line, '(a,a,a,a,a)') ' terminator         : ', trim(r%terminator), ';  depth tolerance ', trim(esf(r%depth_tol, 1)), &
        ' (relative)'
      call wline(u, line)
      if (r%contour_poles > 0) then
        write (line, '(a,i0,a,i0,a,i0,a,a,a)') ' occupation contour : ', r%contour_nodes, ' nodes: main line ', r%contour_main, &
          ', Matsubara poles ', r%contour_poles, ' (kT ', trim(ftoa(r%fermi_kt * 1000.0_dp, 3)), ' mRy)'
      else
        write (line, '(a,i0,a,i0,a)') ' occupation contour : ', r%contour_nodes, ' nodes: closed contour, main line ', &
          r%contour_main, ' (no Fermi smearing)'
      end if
      call wline(u, line)
      write (line, '(a,a,a,a,a)') '                      bottom ', trim(ftoa(r%ebot_offset, 3)), &
        ' Ry below the lowest E_nu;  E_F probe eta ', trim(ftoa(r%ed_eta, 3)), ' Ry'
      call wline(u, line)
      call wline(u, ' electrostatics     : Madelung of the sublattice averages (Ewald, decoration cell)')
      if (r%unscreened) then
        call wline(u, '                      screening: none (sim_beta = 0; bare point charges)')
      else
        write (line, '(a,a,a)') '                      + screened charge deviations (SIM, sim_beta ', trim(ftoa(r%sim_beta, 3)), ')'
        call wline(u, line)
        write (line, '(a,a,a,a,a,a,a,a)') ' screening model    : hardness a ', trim(ftoa(r%screen_a, 4)), ' Ry/e   G00 ', &
          trim(ftoa(r%screen_g00, 4)), '   G1 ', trim(ftoa(r%screen_g1, 4)), '   W00/G00 ', trim(ftoa(r%screen_w00_g00, 4))
        call wline(u, line)
      end if
    end if
    write (line, '(a,a,a,i0)') ' mixing             : Anderson, beta ', trim(ftoa(r%mix_beta, 3)), ', history ', r%mix_history
    call wline(u, line)
    write (line, '(a,a,a,a,a)') '                      step limits ', trim(ftoa(r%dvmax, 2)), ' Ry (V), ', trim(ftoa(r%denumax, 2)), &
      ' Ry (E_nu)'
    call wline(u, line)
    write (line, '(a,a,a)') ' convergence        : max over r of r/(1+r) |V_out - V_in| < ', trim(esf(r%tol, 1)), ' Ry'
    call wline(u, line)
    if (r%sro%uform) then
      call wline(u, '                      and a |u_out - u_in| of the SRO drive;  E_nu inside its band window')
    else
      call wline(u, '                      and E_nu inside its band window')
    end if
  end subroutine rep_method


  !> Section 4.4: Structure
  subroutine rep_structure(u, r)
    integer, intent(in) :: u
    type(rep_run_t), intent(in) :: r
    character(len=90) :: line
    character(len=64) :: item
    integer :: k, i
    real(dp) :: vsite, vsph

    call rep_section(u, 'structure')
    if (trim(r%method) == 'kspace') then
      if (trim(r%lattice) == 'cell') then
        call wline(u, ' lattice                 : explicit cell (primitive vectors below)')
      else
        call wline(u, ' lattice                 : ' // trim(r%lattice))
      end if
    else
      call wline(u, ' graph lattice           : ' // trim(r%lattice))
    end if
    write (line, '(a,f11.6,a,f11.6,a)') ' lattice constant        : ', r%alat, ' Bohr  (', r%alat * bohr_to_ang, ' Angstrom)'
    call wline(u, line)
    if (trim(r%method) == 'kspace') then
      call wline(u, ' primitive vectors (alat units):')
      do k = 1, 3
        write (line, '(a,i1,a,3f11.6)') '   a', k, '                    :', r%avec(:, k)
        call wline(u, line)
      end do
      vsite = r%cell_vol / real(max(1, r%nat), dp)
      write (line, '(a,i0)') ' atoms per cell          : ', r%nat
      call wline(u, line)
      write (line, '(a,f11.5,a,f11.5,a)') ' volume per atom         : ', vsite, ' Bohr^3  (', vsite * bohr_to_ang**3, ' Angstrom^3)'
      call wline(u, line)
      vsph = 0.0_dp
      if (allocated(r%sites)) vsph = sum([(4.0_dp * pi / 3.0_dp * r%sites(i)%sws**3, i=1, r%nat)])
      write (line, '(a,f9.5)') ' sphere / cell volume    : ', vsph / max(r%cell_vol, tiny(1.0_dp))
      call wline(u, line)
    else
      write (line, '(a,i0,a,i0,a)') ' decoration cell         : ', r%period, '^3 conventional cubes, ', r%nsite, ' sites'
      call wline(u, line)
      write (line, '(a,f11.6,a)') ' sphere radius S         : ', r%sws_avg, ' Bohr (equal spheres, space filling)'
      call wline(u, line)
      vsite = r%cell_vol / real(max(1, r%nsite), dp)
      write (line, '(a,f11.5,a,f11.5,a)') ' volume per site         : ', vsite, ' Bohr^3  (', vsite * bohr_to_ang**3, ' Angstrom^3)'
      call wline(u, line)
    end if

    if (allocated(r%sub)) then
      call wline(u, '')
      call wline(u, ' sublattice  site fraction  species and concentrations')
      call wline(u, ' ' // repeat('-', 69))
      do k = 1, size(r%sub)
        line = ' ' // lj(r%sub(k)%name, 10) // rj(ftoa(r%sub(k)%frac, 4), 15) // ' '
        if (allocated(r%sub(k)%species) .and. allocated(r%sub(k)%conc)) then
          do i = 1, r%sub(k)%nspec
            item = trim(r%sub(k)%species(i)) // ' ' // trim(ftoa(r%sub(k)%conc(i), 4))
            if (len_trim(line) + len_trim(item) + 3 > 88) then
              call wline(u, line)
              line = repeat(' ', 27) // trim(item)
            else
              line = line(1:max(27, len_trim(line))) // '  ' // trim(item)
            end if
          end do
        end if
        call wline(u, line)
      end do
    end if

    if (allocated(r%sites) .and. r%nat > 0) then
      call wline(u, '')
      call wline(u, ' site  species  class      x (alat)    y (alat)    z (alat)      S [Bohr]')
      call wline(u, ' ' // repeat('-', 75))
      do i = 1, r%nat
        write (line, '(i5,2x,a8,i6,1x,3f12.6,f14.5)') r%sites(i)%ia, r%sites(i)%species, r%sites(i)%iclass, &
          r%sites(i)%pos, r%sites(i)%sws
        call wline(u, line)
      end do
    end if
  end subroutine rep_structure

  !> Section 4.5: Short-range order (only with [[sro]])
  subroutine rep_sro(u, r)
    integer, intent(in) :: u
    type(rep_run_t), intent(in) :: r
    character(len=90) :: line
    character(len=32) :: buf
    character(len=64) :: ustr
    integer :: i, j

    if (.not. r%sro%has_sro) return

    call rep_section(u, 'short-range order')
    write (line, '(a,a,a,a,a,i0,a)') ' first shell: centre on ', trim(r%sro%sub_centre), ', neighbours on ', &
      trim(r%sro%sub_neigh), ', ', r%sro%coord_num, ' sites'
    call wline(u, line)

    if (allocated(r%sro%species) .and. allocated(r%sro%cond_p)) then
      line = '   centre '
      do j = 1, r%sro%nspec
        line = trim(line) // '   P(' // trim(r%sro%species(j)) // '|centre)'
      end do
      do j = 1, r%sro%nspec
        line = trim(line) // '  alpha(' // trim(r%sro%species(j)) // ')'
      end do
      call wline(u, line)
      call wline(u, '   --------------------------------------------------------------')
      do i = 1, r%sro%nspec
        write (line, '(3x,a8)') r%sro%species(i)
        do j = 1, r%sro%nspec
          write (buf, '(f14.4)') r%sro%cond_p(i, j)
          line = trim(line) // buf
        end do
        if (allocated(r%sro%alpha)) then
          do j = 1, r%sro%nspec
            write (buf, '(f11.4)') r%sro%alpha(i, j)
            line = trim(line) // buf
          end do
        end if
        call wline(u, line)
      end do
    end if

    call wline(u, ' alpha_ij = 1 - P(j|i) / c_j.  For i /= j: alpha > 0 means fewer unlike neighbours')
    call wline(u, ' (clustering), alpha < 0 more (ordering).')
    call wline(u, '')

    if (allocated(r%sro%mode_lam) .and. r%sro%nmodes > 0) then
      call wline(u, ' screening modes of the SRO closure')
      call wline(u, '   mode  lambda(Ahat)    Gamma=dq/db    G00+G1 lam    F [Ry/e]')
      call wline(u, '   -------------------------------------------------------------')
      do i = 1, r%sro%nmodes
        if (r%sro%mode_inf(i)) then
          write (line, '(3x,i4,f14.4,2f15.5,a)') i, r%sro%mode_lam(i), r%sro%mode_gam(i), &
            r%sro%mode_gform(i), '   infinite (clustered)'
        else
          write (line, '(3x,i4,f14.4,2f15.5,f12.4)') i, r%sro%mode_lam(i), r%sro%mode_gam(i), &
            r%sro%mode_gform(i), r%sro%mode_f(i)
        end if
        call wline(u, line)
      end do
    end if

    if (allocated(r%sro%udrive) .and. allocated(r%comp)) then
      ustr = ''
      do i = 1, size(r%sro%udrive)
        ustr = trim(ustr) // ' ' // trim(r%comp(i)%label) // ' ' // trim(ftoa(r%sro%udrive(i), 5))
      end do
      call wline(u, ' drive u [e] (last iteration): ' // trim(ustr))
    end if
  end subroutine rep_sro


  !> Section 4.6: SCF history (main table and per-component blocks)
  subroutine rep_scf_history(u, r)
    integer, intent(in) :: u
    type(rep_run_t), intent(in) :: r
    character(len=90) :: line
    character(len=32) :: buf
    integer :: it, p, pstart, pend, bsize
    real(dp) :: sitem

    call rep_section(u, 'SCF history')
    call wline(u, ' iter' // rj('E_tot [Ry]', 16) // rj('dE [Ry]', 12) // rj('dV [Ry]', 11) // rj('dE_nu [Ry]', 11) // &
                  rj('E_F [Ry]', 14) // rj('pass', 6) // rj('wall [s]', 9))
    call wline(u, ' ' // repeat('-', 83))
    do it = 1, r%iter
      write (line, '(i5,f16.8)') it, r%hist_etot(it)
      if (it == 1) then
        line = line(1:21) // rj('---', 12)
      else
        write (buf, '(es12.3)') r%hist_etot(it) - r%hist_etot(it - 1)
        line = line(1:21) // buf(1:12)
      end if
      write (buf, '(es11.3)') r%hist_dv(it); line = line(1:33) // buf(1:11)
      if (allocated(r%hist_denu)) then
        write (buf, '(es11.3)') r%hist_denu(it)
      else
        buf = rj('---', 11)
      end if
      line = line(1:44) // buf(1:11)
      write (buf, '(f14.8)') r%hist_ef(it); line = line(1:55) // buf(1:14)
      if (allocated(r%hist_npass)) then
        write (buf, '(i6)') r%hist_npass(it)
      else
        buf = rj('---', 6)
      end if
      line = line(1:69) // buf(1:6)
      if (allocated(r%hist_wall)) then
        write (buf, '(f9.2)') r%hist_wall(it)
      else
        buf = rj('---', 9)
      end if
      line = line(1:75) // buf(1:9)
      call wline(u, line)
    end do
    call wline(u, ' ' // repeat('-', 83))
    if (r%converged) then
      call wline(u, ' converged after ' // trim(itoa(r%iter)) // ' iterations (dV < ' // trim(esf(r%tol, 1)) // &
                 ' Ry, no E_nu held at its window)')
    else
      call wline(u, ' WARNING: not converged after ' // trim(itoa(r%iter)) // ' iterations (dV < ' // trim(esf(r%tol, 1)) // ' Ry)')
    end if
    call wline(u, '')

    if (.not. (allocated(r%hist_moment) .and. allocated(r%hist_qnet) .and. r%ncomp > 0)) return
    if (r%nspin == 2) then
      call wline(u, ' per component:  m [mu_B per atom]  qnet [e per atom];  site m = sum_p w_p m_p')
      bsize = 3
    else
      call wline(u, ' per component:  qnet [e per atom]')
      bsize = 7
    end if
    pstart = 1
    do while (pstart <= r%ncomp)
      pend = min(r%ncomp, pstart + bsize - 1)
      line = ' iter'
      do p = pstart, pend
        if (r%nspin == 2) line = line(1:len_trim(line)) // rj(trim(r%comp(p)%label) // ' m', 11)
        line = line(1:len_trim(line)) // rj(trim(r%comp(p)%label) // ' qnet', 13)
      end do
      if (r%nspin == 2 .and. pstart == 1) line = line(1:len_trim(line)) // rj('site m', 10)
      call wline(u, line)
      call wline(u, ' ' // repeat('-', len_trim(line) - 1))
      do it = 1, r%iter
        write (line, '(i5)') it
        do p = pstart, pend
          if (r%nspin == 2) then
            write (buf, '(f11.5)') nz(r%hist_moment(it, p), 5); line = line(1:len_trim(line)) // buf(1:11)
          end if
          write (buf, '(f13.5)') nz(r%hist_qnet(it, p), 5); line = line(1:len_trim(line)) // buf(1:13)
        end do
        if (r%nspin == 2 .and. pstart == 1) then
          sitem = sum(r%hist_moment(it, 1:r%ncomp) * r%comp(1:r%ncomp)%weight_per_site)
          write (buf, '(f10.5)') nz(sitem, 5); line = line(1:len_trim(line)) // buf(1:10)
        end if
        call wline(u, line)
      end do
      call wline(u, '')
      pstart = pend + 1
    end do
  end subroutine rep_scf_history


  !> Section 4.7: Total energy decomposition
  subroutine rep_energy(u, r)
    integer, intent(in) :: u
    type(rep_run_t), intent(in) :: r
    character(len=90) :: line
    character(len=32) :: buf
    character(len=22), parameter :: tname(n_et) = [character(len=22) :: 'valence band energy', 'core eigenvalue sum', &
      'double counting', 'sphere Hartree', 'exchange-correlation', 'Madelung', 'screened deviations']
    character(len=24), parameter :: tform(n_et) = [character(len=24) :: 'sum_l (q1 + E_nu Q_l)', '', '-int n (v_in + 2Z/r)', &
      '1/2 int n V_H', 'int n eps_xc', '1/2 q V_M', '1/2 q V_F']
    character(len=16) :: unit
    integer :: k, p, pstart, pend, bsize

    unit = 'per site [Ry]'; if (trim(r%method) == 'kspace') unit = 'per atom [Ry]'
    call rep_section(u, 'total energy')
    call wline(u, ' Kohn-Sham form: output density in the input potential; no Fermi smearing entropy')
    call wline(u, ' ' // lj('term', 46) // rj(unit, 22))
    call wline(u, ' ' // repeat('-', 67))
    do k = 1, n_et
      write (buf, '(f22.8)') nz(r%eterm(k), 8)
      call wline(u, ' ' // tname(k) // tform(k) // buf(1:22))
    end do
    call wline(u, ' ' // repeat('-', 67))
    write (buf, '(f22.8)') r%etot
    call wline(u, ' ' // lj('total energy E_tot', 46) // buf(1:22))

    if (.not. (allocated(r%ecomp) .and. r%ncomp > 1)) return
    call wline(u, '')
    if (trim(r%method) == 'kspace') then
      call wline(u, ' per class [Ry per atom of the class]; last column per atom of the cell')
    else
      call wline(u, ' per component [Ry per atom of the component]; last column per site = sum_p w_p E_p')
    end if
    bsize = 3
    pstart = 1
    do while (pstart <= r%ncomp)
      pend = min(r%ncomp, pstart + bsize - 1)
      line = ' ' // lj('term', 22)
      do p = pstart, pend
        line = line(1:len_trim(line)) // rj(r%comp(p)%label, 16)
      end do
      if (pend == r%ncomp) line = line(1:len_trim(line)) // rj('per site', 16)
      call wline(u, line)
      call wline(u, ' ' // repeat('-', len_trim(line) - 1))
      do k = 1, n_et
        line = ' ' // tname(k)
        do p = pstart, pend
          write (buf, '(f16.8)') nz(r%ecomp(k, p), 8); line = line(1:23 + (p - pstart) * 16) // buf(1:16)
        end do
        if (pend == r%ncomp) then
          write (buf, '(f16.8)') nz(r%eterm(k), 8); line = line(1:23 + (pend - pstart + 1) * 16) // buf(1:16)
        end if
        call wline(u, line)
      end do
      line = ' ' // lj('total energy', 22)
      do p = pstart, pend
        write (buf, '(f16.8)') sum(r%ecomp(:, p)); line = line(1:23 + (p - pstart) * 16) // buf(1:16)
      end do
      if (pend == r%ncomp) then
        write (buf, '(f16.8)') r%etot; line = line(1:23 + (pend - pstart + 1) * 16) // buf(1:16)
      end if
      call wline(u, ' ' // repeat('-', len_trim(line) - 1))
      call wline(u, line)
      call wline(u, '')
      pstart = pend + 1
    end do
  end subroutine rep_energy


  !> Section 4.8: Converged results per component
  subroutine rep_components(u, r)
    integer, intent(in) :: u
    type(rep_run_t), intent(in) :: r
    character(len=90) :: line
    character(len=1), parameter :: lname(0:2) = ['s', 'p', 'd']
    integer :: p, s, l
    real(dp) :: qtot, ntot

    call rep_section(u, 'converged results')
    write (line, '(a,f14.8,a)') ' E_F                       : ', r%ef, ' Ry'
    call wline(u, line)
    call wline(u, ' cell charge residual      : ' // trim(esf(r%qcell, 2)) // ' e per site   (sum_p w_p qnet_p; 0 = neutral)')
    if (r%nspin == 2) then
      write (line, '(a,f12.5,a)') ' spin moment per site      : ', nz(r%moment_site, 5), ' mu_B'
      call wline(u, line)
    end if

    do p = 1, r%ncomp
      call wline(u, '')
      write (line, '(a,a,a,f5.1,a,i0,a,f7.4,a,f7.4)') ' component ', trim(r%comp(p)%label), &
        '   Z ', r%comp(p)%z, '   valence ', r%comp(p)%valence, '   c ', r%comp(p)%concentration, &
        '   weight per site ', r%comp(p)%weight_per_site
      call wline(u, line)
      write (line, '(a,f8.5,a,f10.5,a,f9.5,a,f9.5,a)') '   S ', r%comp(p)%sws, ' Bohr   electrons in sphere ', &
        r%comp(p)%z + r%comp(p)%qnet, ' e   qnet ', nz(r%comp(p)%qnet, 5), ' e   moment ', nz(r%comp(p)%moment, 5), ' mu_B'
      call wline(u, line)
      call wline(u, '')

      call wline(u, '   potential parameters, gamma representation [Ry]')
      call wline(u, '   spin  l' // rj('E_nu', 12) // rj('C', 12) // rj('sqrt(Delta)', 12) // rj('p', 12) // rj('gamma', 12))
      call wline(u, '   ' // repeat('-', 67))
      do s = 1, r%nspin
        do l = 0, 2
          write (line, '(3x,i4,2x,a1,5f12.6)') s, lname(l), r%comp(p)%enu(l, s), r%comp(p)%c(l, s), &
            r%comp(p)%sqrdel(l, s), r%comp(p)%p(l, s), r%comp(p)%gamma(l, s)
          call wline(u, line)
        end do
      end do
      call wline(u, '')

      if (r%nspin == 2) then
        call wline(u, '   moments and DOS at E_F, per spin channel (row 1-2: spin 1 minus spin 2)')
      else
        call wline(u, '   moments and DOS at E_F, both spins')
      end if
      call wline(u, '   spin  l' // rj('q0 [e]', 11) // rj('q1 [e Ry]', 11) // rj('q2 [e Ry^2]', 12) // rj('Q_l [e]', 11) // &
                    rj('N_l(E_F) [1/Ry]', 17))
      call wline(u, '   ' // repeat('-', 69))
      do s = 1, r%nspin
        do l = 0, 2
          write (line, '(3x,i4,2x,a1,2f11.5,f12.5,f11.5,f17.4)') s, lname(l), nz(r%comp(p)%q(0, l, s), 5), &
            nz(r%comp(p)%q(1, l, s), 5), nz(r%comp(p)%q(2, l, s), 5), nz(r%comp(p)%qch(l, s), 5), nz(r%comp(p)%dosef(l, s), 4)
          call wline(u, line)
        end do
      end do
      if (r%nspin == 2) then
        do l = 0, 2
          write (line, '(3x,a4,2x,a1,2f11.5,f12.5,f11.5,f17.4)') ' 1-2', lname(l), &
            nz(r%comp(p)%q(0, l, 1) - r%comp(p)%q(0, l, 2), 5), nz(r%comp(p)%q(1, l, 1) - r%comp(p)%q(1, l, 2), 5), &
            nz(r%comp(p)%q(2, l, 1) - r%comp(p)%q(2, l, 2), 5), nz(r%comp(p)%qch(l, 1) - r%comp(p)%qch(l, 2), 5), &
            nz(r%comp(p)%dosef(l, 1) - r%comp(p)%dosef(l, 2), 4)
          call wline(u, line)
        end do
      end if
      qtot = sum(r%comp(p)%qch(:, 1:r%nspin))
      ntot = sum(r%comp(p)%dosef(:, 1:r%nspin))
      call wline(u, '   ' // repeat('-', 69))
      write (line, '(3x,a,f11.5,f17.4)') lj('sum (valence electrons in the sphere)', 41), qtot, ntot
      call wline(u, line)
    end do
  end subroutine rep_components


  !> Section 4.9: Fermi level card
  subroutine rep_fermi(u, r)
    integer, intent(in) :: u
    type(rep_run_t), intent(in) :: r
    character(len=90) :: line
    real(dp) :: nef, nef_ev, ge, chi
    integer :: p

    nef = r%dosef_site
    nef_ev = nef / ry_to_ev
    ge = gamma_e_const * nef
    chi = chi_p_const * nef

    call rep_section(u, 'Fermi level')
    write (line, '(a,f14.8,a)') ' E_F                                       : ', r%ef, ' Ry'
    call wline(u, line)
    write (line, '(a,f10.4,a,f9.5,a)') ' N(E_F) per site, both spins               : ', nef, ' states/Ry   (', nef_ev, ' states/eV)'
    call wline(u, line)
    if (trim(r%method) == 'kspace') then
      write (line, '(a,f10.4,a)') ' smearing of N(E_F)                        : ', r%k_smear, ' Ry  (k-space BZ integration)'
    else
      write (line, '(a,f10.4,a)') ' broadening of N(E_F)                      : ', r%ed_eta, ' Ry  (G at E_F + i eta)'
    end if
    call wline(u, line)
    write (line, '(a,f10.4,a)') ' Sommerfeld coefficient, bare band         : ', ge, ' mJ/(mol K^2)'
    call wline(u, line)
    call wline(u, ' Pauli susceptibility, bare (no Stoner)    :  ' // trim(esf(chi, 3)) // ' emu/mol')
    call wline(u, ' gamma_e = (pi^2/3) k_B^2 N(E_F) = 0.173247 mJ/(mol K^2) per state/Ry;')
    call wline(u, ' chi_P = mu_B^2 N(E_F) = 2.37605E-06 emu/mol per state/Ry;  mol = mole of sites')
    call wline(u, '')

    call wline(u, ' ' // lj('component', 12) // rj('s [1/Ry]', 11) // rj('p [1/Ry]', 11) // rj('d [1/Ry]', 11) // &
                  rj('total [1/Ry]', 15) // rj('weight', 9))
    call wline(u, ' ' // repeat('-', 68))
    do p = 1, r%ncomp
      write (line, '(1x,a12,3f11.4,f15.4,f9.4)') r%comp(p)%label, sum(r%comp(p)%dosef(0, 1:r%nspin)), &
        sum(r%comp(p)%dosef(1, 1:r%nspin)), sum(r%comp(p)%dosef(2, 1:r%nspin)), sum(r%comp(p)%dosef(:, 1:r%nspin)), &
        r%comp(p)%weight_per_site
      call wline(u, line)
    end do
    call wline(u, ' ' // repeat('-', 68))
    write (line, '(1x,a45,f15.4)') lj('per site (sum_p w_p N_p)', 45), nef
    call wline(u, line)
  end subroutine rep_fermi


  !> Section 4.10: Health checks
  subroutine rep_checks(u, r)
    integer, intent(in) :: u
    type(rep_run_t), intent(in) :: r
    character(len=8) :: st

    call rep_section(u, 'health checks')
    call wline(u, ' ' // lj('check', 44) // '  ' // lj('value', 20) // '  ' // lj('limit', 12) // '  status')
    call wline(u, ' ' // repeat('-', 87))

    if (r%converged) then; st = 'OK'; else; st = 'FAIL'; end if
    call wcheck(u, 'SCF converged (dV)', trim(esf(r%hist_dv(max(1, r%iter)), 1)) // ' Ry', trim(esf(r%tol, 1)) // ' Ry', st)

    if (r%nclamp == 0) then; st = 'OK'; else; st = 'FAIL'; end if
    call wcheck(u, 'E_nu held at its window, last / whole run', trim(itoa(r%nclamp)) // ' / ' // trim(itoa(r%nclamp_run)), '0', st)

    if (trim(r%method) == 'kspace') then
      if (abs(r%eband - r%sumev) <= 1.0e-5_dp) then; st = 'OK'; else; st = 'WARN'; end if
      call wcheck(u, 'band energy check |eband - sumev|', trim(esf(abs(r%eband - r%sumev), 1)) // ' Ry', '1.0E-05 Ry', st)
    else
      if (r%nnoncausal == 0) then; st = 'OK'; else; st = 'WARN'; end if
      call wcheck(u, 'non-causal contour values (-Im G <= 0)', trim(itoa(r%nnoncausal)) // ' of ' // trim(itoa(r%nnode)), '0', st)
      if (r%below_ok) then; st = 'OK'; else; st = 'WARN'; end if
      call wcheck(u, 'valence weight below contour, first / final', trim(esf(r%wbelow1, 1)) // ' / ' // trim(esf(r%wbelow, 1)), &
                  trim(esf(r%wbelow_lim, 1)) // ' e', st)
      if (abs(r%qcell) <= 1.0e-3_dp) then; st = 'OK'; else; st = 'WARN'; end if
      call wcheck(u, 'cell charge residual', trim(esf(r%qcell, 1)) // ' e', '1.0E-03 e', st)
      call wcheck(u, 'recursion depth used, mean', trim(ftoa(r%depth_mean, 2)) // ' of ' // trim(itoa(r%nrec)), '-', 'info')
      call wcheck(u, 'bp tail fallbacks, last iteration', trim(itoa(r%nfall)), '-', 'info')
      call wcheck(u, 'Anderson restarts in the run', trim(itoa(r%mix_nrestart)), '-', 'info')
      if (allocated(r%hist_npass)) call wcheck(u, 'E_F passes, last iteration', trim(itoa(r%hist_npass(max(1, r%iter)))), '-', 'info')
    end if
  end subroutine rep_checks


  subroutine wcheck(u, name, val, lim, st)
    integer, intent(in) :: u
    character(len=*), intent(in) :: name, val, lim, st
    call wline(u, ' ' // lj(name, 44) // '  ' // lj(val, 20) // '  ' // lj(lim, 12) // '  ' // trim(st))
  end subroutine wcheck


  !> Section 4.11: Timing table.  The total is the whole run (r%timer%wall_total, cpu_total as set by the adapter); "other" is
  !> the time outside the timed parts (input, setup outside the SCF, writing the files).
  subroutine rep_timing(u, r)
    integer, intent(in) :: u
    type(rep_run_t), intent(in) :: r
    real(dp) :: wtot, ctot, wsum, csum
    integer :: ib, rate

    wtot = max(r%timer%wall_total, sum(r%timer%wall), 1.0e-6_dp)
    ctot = max(r%timer%cpu_total, sum(r%timer%cpu))

    call rep_section(u, 'timing')
    call wline(u, ' ' // lj('part', 36) // rj('calls', 8) // rj('wall [s]', 12) // rj('cpu [s]', 12) // rj('cpu/wall', 10) // &
               rj('share', 8))
    call wline(u, ' ' // repeat('-', 85))
    wsum = 0.0_dp; csum = 0.0_dp
    do ib = 1, n_tb
      if (r%timer%ncall(ib) > 0 .or. r%timer%wall(ib) > 0.0_dp) then
        wsum = wsum + r%timer%wall(ib)
        csum = csum + r%timer%cpu(ib)
        call trow(tb_name(ib), trim(itoa(r%timer%ncall(ib))), r%timer%wall(ib), r%timer%cpu(ib))
      end if
    end do
    call trow('other (outside the parts above)', '', max(0.0_dp, wtot - wsum), max(0.0_dp, ctot - csum))
    call wline(u, ' ' // repeat('-', 85))
    call trow('total', '', wtot, ctot)

    if (r%timer%work(tb_chains) > 0_i8) then
      rate = int(real(r%timer%work(tb_chains), dp) / max(1.0e-4_dp, r%timer%wall(tb_chains)))
      call wline(u, ' chain evaluations, SCF contour: ' // trim(i8toa(r%timer%work(tb_chains))) // '   (' // trim(itoa(rate)) // &
                 ' per wall second)')
    end if
    call wline(u, ' rank 0 of ' // trim(itoa(r%mpi_ranks)) // '; CPU time is the process time of this rank (all threads)')
  contains
    subroutine trow(name, calls, w, c)
      character(len=*), intent(in) :: name, calls
      real(dp), intent(in) :: w, c
      character(len=16) :: bw, bc, bq, bs
      write (bw, '(f12.2)') w
      write (bc, '(f12.2)') c
      if (w >= 0.01_dp) then
        write (bq, '(f10.1)') c / w
      else
        bq = rj('-', 10)
      end if
      write (bs, '(f7.1,a1)') 100.0_dp * w / wtot, '%'
      call wline(u, ' ' // lj(name, 36) // rj(calls, 8) // bw(1:12) // bc(1:12) // bq(1:10) // bs(1:8))
    end subroutine trow
  end subroutine rep_timing


  !> Section 4.12: Warnings (the list of collect_warnings)
  subroutine rep_warnings(u, r)
    integer, intent(in) :: u
    type(rep_run_t), intent(in) :: r
    character(len=128), allocatable :: w(:)
    integer :: i
    call rep_section(u, 'warnings')
    call collect_warnings(r, w)
    if (size(w) == 0) then
      call wline(u, ' none')
    else
      do i = 1, size(w)
        call write_wrapped(u, ' ', trim(w(i)), '   ')
      end do
    end if
  end subroutine rep_warnings

  !> Section 4.12 (files): Files written
  subroutine rep_files(u, r)
    integer, intent(in) :: u
    type(rep_run_t), intent(in) :: r
    character(len=90) :: line
    character(len=28) :: fn
    call rep_section(u, 'files written')
    if (len_trim(r%dos_file) > 0) then
      fn = trim(r%dos_file)
      write (line, '(1x,a28,2x,a)') fn, 'DOS curve'
      call wline(u, line)
    end if
    if (len_trim(r%pot_file) > 0) then
      fn = trim(r%pot_file)
      write (line, '(1x,a28,2x,a)') fn, 'restart potentials (last iteration)'
      call wline(u, line)
    end if
    if (len_trim(r%params_file) > 0) then
      fn = trim(r%params_file)
      write (line, '(1x,a28,2x,a)') fn, 'potential parameters'
      call wline(u, line)
    end if
    if (len_trim(r%json_file) > 0) then
      fn = trim(r%json_file)
      write (line, '(1x,a28,2x,a)') fn, 'machine-readable summary'
      call wline(u, line)
    end if
    if (len_trim(r%out_file) > 0) then
      fn = trim(r%out_file)
      write (line, '(1x,a28,2x,a)') fn, 'text report'
      call wline(u, line)
    end if
  end subroutine rep_files


  !> Section 4.13: End of run
  subroutine rep_close(u, r)
    integer, intent(in) :: u
    type(rep_run_t), intent(in) :: r
    character(len=90) :: line
    character(len=128), allocatable :: w(:)

    call collect_warnings(r, w)
    call rep_section(u, 'end of run')
    write (line, '(a,a,a,f10.2,a)') ' finished           : ', trim(fmt_dt(r%finish_date, r%finish_time)), &
      '   total wall time ', r%wall_total, ' s'
    call wline(u, line)
    if (.not. r%converged) then
      call wline(u, ' status             : NORMAL TERMINATION   SCF NOT CONVERGED')
    else if (size(w) > 0) then
      call wline(u, ' status             : NORMAL TERMINATION   SCF CONVERGED WITH WARNINGS')
    else
      call wline(u, ' status             : NORMAL TERMINATION   SCF CONVERGED')
    end if
    call wline(u, repeat('=', 90))
  end subroutine rep_close

  !> Write the entire text report to a file.
  subroutine rep_write_report(r, filename)
    type(rep_run_t), intent(in) :: r
    character(len=*), intent(in) :: filename
    integer :: u

    open (newunit=u, file=filename, status='replace', action='write')
    call rep_banner(u, r)
    call rep_echo_input(u, r%infile)
    call rep_method(u, r)
    call rep_structure(u, r)
    if (r%sro%has_sro) call rep_sro(u, r)
    call rep_scf_history(u, r)
    call rep_energy(u, r)
    call rep_components(u, r)
    call rep_fermi(u, r)
    call rep_checks(u, r)
    call rep_timing(u, r)
    call rep_warnings(u, r)
    call rep_files(u, r)
    call rep_close(u, r)
    close (u)
  end subroutine rep_write_report

  !> Write machine-readable JSON summary (<out>.json)
  subroutine rep_write_json(r, filename)
    type(rep_run_t), intent(in) :: r
    character(len=*), intent(in) :: filename
    integer :: u, p, it, k
    character(len=32) :: et_keys(n_et)
    character(len=128), allocatable :: wl(:)

    et_keys = ['band           ', 'core           ', 'double_counting', 'hartree        ', &
               'xc             ', 'madelung       ', 'screened       ']

    open (newunit=u, file=filename, status='replace', action='write')
    write (u, '(a)') '{'
    write (u, '(a)') '  "schema": "asr-summary-1",'
    write (u, '(a)') '  "provenance": {'
    write (u, '(a,a,a)') '    "version": "', trim(program_version), '",'
    write (u, '(a,a,a)') '    "git": "', trim(r%git_rev), '",'
    if (r%git_dirty) then
      write (u, '(a)') '    "dirty": true,'
    else
      write (u, '(a)') '    "dirty": false,'
    end if
    write (u, '(a,a,a)') '    "compiler": "', trim(r%compiler_ver), '",'
    write (u, '(a,a,a)') '    "host": "', trim(r%host), '",'
    write (u, '(a,i0,a)') '    "mpi_ranks": ', r%mpi_ranks, ','
    write (u, '(a,i0,a)') '    "omp_threads": ', r%omp_threads, ','
    write (u, '(a,a,a)') '    "input": "', trim(r%infile), '",'
    write (u, '(a,a,a)') '    "input_fnv64": "', trim(r%fingerprint), '"'
    write (u, '(a)') '  },'

    write (u, '(a)') '  "method": {'
    write (u, '(a,a,a)') '    "route": "', trim(r%route), '",'
    write (u, '(a,i0,a)') '    "levels": ', r%nlevels, ','
    write (u, '(a,i0,a)') '    "nrec": ', r%nrec, ','
    write (u, '(a,i0,a)') '    "rank": ', r%maxrank, ','
    write (u, '(a,i0,a)') '    "nspin": ', r%nspin, ','
    write (u, '(a,a,a)') '    "xc": "', trim(r%xc), '"'
    write (u, '(a)') '  },'

    write (u, '(a)') '  "status": {'
    if (r%converged) then
      write (u, '(a)') '    "converged": true,'
    else
      write (u, '(a)') '    "converged": false,'
    end if
    write (u, '(a,i0,a)') '    "iterations": ', r%iter, ','
    call collect_warnings(r, wl)
    write (u, '(a)', advance='no') '    "warnings": ['
    do k = 1, size(wl)
      write (u, '(a,a,a)', advance='no') '"', trim(json_escape(wl(k))), '"'
      if (k < size(wl)) write (u, '(a)', advance='no') ', '
    end do
    write (u, '(a)') ']'
    write (u, '(a)') '  },'

    write (u, '(a,es24.16,a)') '  "fermi_energy_ry": ', r%ef, ','
    write (u, '(a,es24.16,a)') '  "total_energy_ry_per_site": ', r%etot, ','

    write (u, '(a)') '  "energy_terms_ry_per_site": {'
    do k = 1, n_et
      write (u, '(a,a,a,es24.16)', advance='no') '    "', trim(et_keys(k)), '": ', r%eterm(k)
      if (k < n_et) then; write (u, '(a)') ','; else; write (u, '(a)') ''; end if
    end do
    write (u, '(a)') '  },'

    write (u, '(a)') '  "components": ['
    do p = 1, r%ncomp
      write (u, '(a)') '    {'
      write (u, '(a,a,a)') '      "label": "', trim(r%comp(p)%label), '",'
      write (u, '(a,a,a)') '      "species": "', trim(r%comp(p)%species), '",'
      write (u, '(a,a,a)') '      "sublattice": "', trim(r%comp(p)%subname), '",'
      write (u, '(a,f6.1,a)') '      "z": ', r%comp(p)%z, ','
      write (u, '(a,f8.4,a)') '      "concentration": ', r%comp(p)%concentration, ','
      write (u, '(a,f8.4,a)') '      "weight_per_site": ', r%comp(p)%weight_per_site, ','
      write (u, '(a,es24.16,a)') '      "qnet_e": ', r%comp(p)%qnet, ','
      write (u, '(a,es24.16,a)') '      "moment_mub": ', r%comp(p)%moment, ','
      write (u, '(a,es24.16,a,es24.16,a,es24.16,a)') '      "occupation_e": {"s": ', sum(r%comp(p)%q(0, 0, :)), &
        ', "p": ', sum(r%comp(p)%q(0, 1, :)), ', "d": ', sum(r%comp(p)%q(0, 2, :)), '},'
      write (u, '(a,es24.16,a,es24.16,a,es24.16,a)') '      "dos_ef_per_ry": {"s": ', sum(r%comp(p)%dosef(0, :)), &
        ', "p": ', sum(r%comp(p)%dosef(1, :)), ', "d": ', sum(r%comp(p)%dosef(2, :)), '},'
      write (u, '(a)') '      "energy_terms_ry": {'
      do k = 1, n_et
        write (u, '(a,a,a,es24.16)', advance='no') '        "', trim(et_keys(k)), '": ', r%comp(p)%ecomp(k)
        if (k < n_et) then; write (u, '(a)') ','; else; write (u, '(a)') ''; end if
      end do
      write (u, '(a)') '      }'
      if (p < r%ncomp) then; write (u, '(a)') '    },'; else; write (u, '(a)') '    }'; end if
    end do
    write (u, '(a)') '  ],'

    write (u, '(a,es24.16,a)') '  "dos_ef_per_site_per_ry": ', r%dosef_site, ','

    write (u, '(a)') '  "history": {'
    write (u, '(a)', advance='no') '    "etot_ry": ['
    if (allocated(r%hist_etot)) then
      do it = 1, r%iter
        write (u, '(es24.16)', advance='no') r%hist_etot(it)
        if (it < r%iter) write (u, '(a)', advance='no') ', '
      end do
    end if
    write (u, '(a)') '],'

    write (u, '(a)', advance='no') '    "dv_ry": ['
    if (allocated(r%hist_dv)) then
      do it = 1, r%iter
        write (u, '(es24.16)', advance='no') r%hist_dv(it)
        if (it < r%iter) write (u, '(a)', advance='no') ', '
      end do
    end if
    write (u, '(a)') '],'

    write (u, '(a)', advance='no') '    "ef_ry": ['
    if (allocated(r%hist_ef)) then
      do it = 1, r%iter
        write (u, '(es24.16)', advance='no') r%hist_ef(it)
        if (it < r%iter) write (u, '(a)', advance='no') ', '
      end do
    end if
    write (u, '(a)') '],'

    write (u, '(a)', advance='no') '    "npass": ['
    if (allocated(r%hist_npass)) then
      do it = 1, r%iter
        write (u, '(i0)', advance='no') r%hist_npass(it)
        if (it < r%iter) write (u, '(a)', advance='no') ', '
      end do
    end if
    write (u, '(a)') '],'

    write (u, '(a)', advance='no') '    "wall_s": ['
    if (allocated(r%hist_wall)) then
      do it = 1, r%iter
        write (u, '(f10.2)', advance='no') r%hist_wall(it)
        if (it < r%iter) write (u, '(a)', advance='no') ', '
      end do
    end if
    write (u, '(a)') ']'
    write (u, '(a)') '  },'

    write (u, '(a)') '  "checks": {'
    write (u, '(a,i0,a)') '    "nnoncausal": ', r%nnoncausal, ','
    write (u, '(a,i0,a)') '    "nnode": ', r%nnode, ','
    write (u, '(a,es24.16,a)') '    "wbelow_first_e": ', r%wbelow1, ','
    write (u, '(a,es24.16,a)') '    "wbelow_final_e": ', r%wbelow, ','
    write (u, '(a,i0,a)') '    "nclamp_last": ', r%nclamp, ','
    write (u, '(a,i0,a)') '    "nclamp_run": ', r%nclamp_run, ','
    write (u, '(a,es24.16)') '    "qcell_e": ', r%qcell
    write (u, '(a)') '  },'

    write (u, '(a)') '  "timing_s": {'
    write (u, '(a,f12.2,a)') '    "total_wall": ', r%timer%wall_total, ','
    write (u, '(a,f12.2,a)') '    "total_cpu": ', r%timer%cpu_total, ','
    write (u, '(a)') '    "parts": {'
    do k = 1, n_tb
      write (u, '(a,a,a,f12.2)', advance='no') '      "', trim(tb_name(k)), '": ', r%timer%wall(k)
      if (k < n_tb) then; write (u, '(a)') ','; else; write (u, '(a)') ''; end if
    end do
    write (u, '(a)') '    }'
    write (u, '(a)') '  }'
    write (u, '(a)') '}'
    close (u)
  end subroutine rep_write_json

  !> Write self-describing DOS file header
  subroutine rep_write_dos_header(u, r, npoints, emin, emax)
    integer, intent(in) :: u
    type(rep_run_t), intent(in) :: r
    integer, intent(in) :: npoints
    real(dp), intent(in) :: emin, emax
    character(len=1), parameter :: lname(0:2) = ['s', 'p', 'd']
    integer :: col, p, l, s, ia

    write (u, '(a)') '# ASR density of states'
    write (u, '(a,a,a,a,a,a,a)') '# run        : ', trim(r%title), '   (', trim(r%infile), ', git ', trim(r%git_rev), ')'
    if (trim(r%method) == 'kspace') then
      write (u, '(a,i0,a,f8.6,a)') '# route      : kspace, nk ', r%k_nk, ';  smearing ', r%k_smear, ' Ry'
    else
      write (u, '(a,a,a,i0,a,f8.6,a)') '# route      : recursion ', trim(r%route), ', levels ', r%nlevels, &
        ';  broadening eta ', r%eta, ' Ry (G at E + i eta)'
    end if
    write (u, '(a,f12.8,a)') '# E_F        : ', r%ef, ' Ry'
    write (u, '(a,i0,a,f10.6,a,f10.6,a)') '# mesh       : ', npoints, ' energies, ', emin, ' .. ', emax, ' Ry'
    write (u, '(a)') '# units      : E in Ry; DOS in states/Ry per atom of a component, or per site'
    if (r%nspin == 1) then
      write (u, '(a,i0,a)') '# spin       : nspin ', r%nspin, '; DOS columns are sums over both spins'
    else
      write (u, '(a,i0,a)') '# spin       : nspin ', r%nspin, '; spin polarised'
    end if

    col = 1
    write (u, '(a,i2,a)') '# column ', col, '  : E [Ry]'; col = col + 1

    if (trim(r%method) == 'kspace') then
      ! Columns of k-space DOS
      do ia = 1, r%nat
        do s = 1, r%nspin
          do l = 0, 2
            write (u, '(a,i2,a,i0,a,a,a,i0,a,a1)') '# column ', col, '  : atom ', ia, ' (', &
              trim(r%sites(ia)%species), ') spin ', s, ' ', lname(l)
            col = col + 1
          end do
        end do
      end do
    else
      ! Columns of recursion DOS
      write (u, '(a,i2,a)') '# column ', col, '  : total per site, sum_p w_p N_p'; col = col + 1
      do p = 1, r%ncomp
        do l = 0, 2
          write (u, '(a,i2,a,a,a,a1)') '# column ', col, '  : ', trim(r%comp(p)%label), ' ', lname(l)
          col = col + 1
        end do
      end do
      write (u, '(a,i2,a)') '# column ', col, '  : total per site, spin 1'; col = col + 1
      if (r%nspin == 2) then
        write (u, '(a,i2,a)') '# column ', col, '  : total per site, spin 2'; col = col + 1
        do p = 1, r%ncomp
          do s = 1, 2
            do l = 0, 2
              write (u, '(a,i2,a,a,a,i0,a,a1)') '# column ', col, '  : ', trim(r%comp(p)%label), ' spin ', s, ' ', lname(l)
              col = col + 1
            end do
          end do
        end do
      end if
    end if
  end subroutine rep_write_dos_header


  !> All warnings of the run: the ones the driver collected, then the failed or doubtful health checks.  The report's warning
  !> section, its status line and the JSON summary use this one list.
  subroutine collect_warnings(r, w)
    type(rep_run_t), intent(in) :: r
    character(len=128), allocatable, intent(out) :: w(:)
    character(len=128) :: buf(64)
    integer :: n, i
    n = 0
    if (allocated(r%warnings)) then
      do i = 1, size(r%warnings)
        if (len_trim(r%warnings(i)) > 0) call add(trim(r%warnings(i)))
      end do
    end if
    if (.not. r%converged) call add('[SCF] not converged after ' // trim(itoa(r%iter)) // ' iterations')
    if (r%nclamp > 0) call add('[health checks] E_nu held at its band window in the last iteration (' // &
                               trim(itoa(r%nclamp)) // ' channels): not a fixed point')
    if (r%nnoncausal > 0) call add('[health checks] Green function not causal at ' // trim(itoa(r%nnoncausal)) // ' of ' // &
                                   trim(itoa(r%nnode)) // ' contour node values')
    if (.not. r%below_ok) call add('[health checks] valence weight below the contour bottom (the charges miss it)')
    if (abs(r%qcell) > 1.0e-3_dp) call add('[health checks] cell charge residual above 1.0E-03 e per site')
    if (trim(r%method) == 'kspace' .and. abs(r%eband - r%sumev) > 1.0e-5_dp) &
      call add('[health checks] band energy from the moments differs from the eigenvalue sum')
    allocate (w(n))
    w = buf(1:n)
  contains
    subroutine add(t)
      character(len=*), intent(in) :: t
      if (n < size(buf)) then
        n = n + 1; buf(n) = t
      end if
    end subroutine add
  end subroutine collect_warnings

  !> Compiler flags without the per-file parts of compiler_options() (-c, -o file, -I dir, module directories, file names)
  !> and without repeated flags.
  function clean_flags(opts) result(out)
    character(len=*), intent(in) :: opts
    character(len=len(opts)) :: out
    character(len=512) :: tok
    integer :: p, q, n
    logical :: skipnext
    out = ''; skipnext = .false.
    n = len_trim(opts); p = 1
    do while (p <= n)
      do while (p <= n .and. opts(p:p) == ' '); p = p + 1; end do
      if (p > n) exit
      q = p
      do while (q <= n .and. opts(q:q) /= ' '); q = q + 1; end do
      tok = opts(p:q - 1); p = q
      if (skipnext) then
        skipnext = .false.; cycle
      end if
      select case (trim(tok))
      case ('-c')
        cycle
      case ('-o', '-module', '-J', '-I')
        skipnext = .true.; cycle
      end select
      if (tok(1:1) /= '-') cycle
      if (tok(1:2) == '-I' .or. tok(1:2) == '-J' .or. index(tok, '-module=') == 1) cycle
      if (index(' ' // trim(out) // ' ', ' ' // trim(tok) // ' ') > 0) cycle
      out = trim(out) // ' ' // trim(tok)
    end do
    out = adjustl(out)
  end function clean_flags

  !> date_and_time strings (ccyymmdd, hhmmss.sss) as 'ccyy-mm-dd hh:mm:ss'
  function fmt_dt(d, t) result(s)
    character(len=*), intent(in) :: d, t
    character(len=19) :: s
    s = ''
    if (len_trim(d) < 8 .or. len_trim(t) < 6) return
    s = d(1:4) // '-' // d(5:6) // '-' // d(7:8) // ' ' // t(1:2) // ':' // t(3:4) // ':' // t(5:6)
  end function fmt_dt

  !> absolute path of the working directory: workdir if absolute, else $PWD joined with it
  function abs_dir(workdir) result(s)
    character(len=*), intent(in) :: workdir
    character(len=1024) :: s
    character(len=1024) :: pwd
    integer :: st
    if (len_trim(workdir) > 0) then
      if (workdir(1:1) == '/') then
        s = workdir; return
      end if
    end if
    call get_environment_variable('PWD', pwd, status=st)
    if (st /= 0 .or. len_trim(pwd) == 0) pwd = '.'
    if (len_trim(workdir) == 0 .or. trim(workdir) == '.' .or. trim(workdir) == './') then
      s = pwd
    else if (workdir(1:2) == './') then
      s = trim(pwd) // '/' // workdir(3:)
    else
      s = trim(pwd) // '/' // workdir
    end if
  end function abs_dir

  !> s left-justified in a field of width w (cut when longer)
  function lj(s, w) result(r)
    character(len=*), intent(in) :: s
    integer, intent(in) :: w
    character(len=w) :: r
    r = s
  end function lj

  !> s (trimmed) right-justified in a field of width w; a longer string keeps its first w characters
  function rj(s, w) result(r)
    character(len=*), intent(in) :: s
    integer, intent(in) :: w
    character(len=w) :: r
    integer :: n
    n = len_trim(s)
    if (n >= w) then
      r = s(1:w)
    else
      r = repeat(' ', w - n) // s(1:n)
    end if
  end function rj

  !> x with values that print as zero at dec decimals set to +0 (no '-0.00000' in tables)
  elemental function nz(x, dec) result(y)
    real(dp), intent(in) :: x
    integer, intent(in) :: dec
    real(dp) :: y
    y = x
    if (abs(x) < 0.5_dp * 10.0_dp**(-dec)) y = 0.0_dp
  end function nz

  !> x in es format with d decimals, left-justified (es8.1 for d = 1)
  function esf(x, d) result(s)
    real(dp), intent(in) :: x
    integer, intent(in) :: d
    character(len=24) :: s
    character(len=16) :: f
    write (f, '(a,i0,a,i0,a)') '(es', d + 8, '.', d, ')'
    write (s, f) x
    s = adjustl(s)
  end function esf

  !> JSON string body: backslash and double quote escaped
  function json_escape(t) result(s)
    character(len=*), intent(in) :: t
    character(len=2 * len(t)) :: s
    integer :: i, k
    s = ''; k = 0
    do i = 1, len_trim(t)
      if (t(i:i) == '"' .or. t(i:i) == '\') then
        k = k + 1; s(k:k) = '\'
      end if
      k = k + 1; s(k:k) = t(i:i)
    end do
  end function json_escape

  function i8toa(i) result(s)
    integer(i8), intent(in) :: i
    character(len=24) :: s
    write (s, '(i0)') i
    s = adjustl(s)
  end function i8toa

  ! Helper: integer to string
  function itoa(i) result(s)
    integer, intent(in) :: i
    character(len=16) :: s
    write (s, '(i0)') i
    s = adjustl(s)
  end function itoa

  ! Helper: float to string
  function ftoa(x, dec) result(s)
    real(dp), intent(in) :: x
    integer, intent(in) :: dec
    character(len=24) :: s
    character(len=16) :: fmt
    write (fmt, '(a,i0,a)') '(f16.', dec, ')'
    write (s, fmt) x
    s = adjustl(s)
  end function ftoa

end module asr_report
