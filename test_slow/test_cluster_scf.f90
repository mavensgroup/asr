!> Gate C of notes/DESIGN_structures.md: the class SCF on the site cluster of the structure (asr_crystal_graph) reproduces the
!> class SCF of the periodic half-lattice box (asr_lattices) where that box does not wrap around.  At levels 2 the closure reaches
!> 4 bonds (at most 4 alat with the default cluster_radius 1.2), so a box of 10 cubes holds every reachable site once.  The
!> reference numbers below come from the box code (commit 3353591, supercell = 10).  Small settings, two SCF iterations from
!> free atoms; the numbers are not converged.  Cases: bcc FeCo (spin), L2_1 Co2Fe0.4Cr0.6Al (period 2, three sublattices,
!> partial disorder, spin), fcc FeNiCr with first-shell SRO, B2-type CuZn on two sublattices through route ed3.
!> Checked: the augmented-space node count of every sublattice (exact), E_F, the total energy, every component moment and net
!> charge.  Run with REF=print as the first argument to print the numbers instead.
program test_cluster_scf
  use asr_kinds, only: dp
  use asr_driver, only: run_input_t, read_run_input, class_species, run_input
  use asr_class_scf, only: class_scf, class_scf_result_t
  use asr_testing
  implicit none
  character(len=1), parameter :: nl = new_line('a')
  character(len=*), parameter :: common = '[scf]'//nl//'xc = "vbh"'//nl//'mixing = 0.2'//nl//'max_iterations = 2'//nl// &
    'tolerance = 1.0e3'//nl
  character(len=*), parameter :: asr2 = '[asr]'//nl//'levels = 2'//nl//'recursion_steps = 8'//nl//'contour_points = 16'//nl// &
    'fermi_kt = 0.004'//nl
  character(len=*), parameter :: sp_fe = '[[species]]'//nl//'symbol = "Fe"'//nl//'z = 26'//nl//'nval = 8'//nl// &
    'konfig = [4, 4, 3, 4]'//nl//'occupations = [2.0, 6.0]'//nl
  character(len=*), parameter :: sp_co = '[[species]]'//nl//'symbol = "Co"'//nl//'z = 27'//nl//'nval = 9'//nl// &
    'konfig = [4, 4, 3, 4]'//nl//'occupations = [2.0, 7.0]'//nl
  character(len=*), parameter :: sp_cr = '[[species]]'//nl//'symbol = "Cr"'//nl//'z = 24'//nl//'nval = 6'//nl// &
    'konfig = [4, 4, 3, 4]'//nl//'occupations = [1.0, 5.0]'//nl
  character(len=*), parameter :: sp_ni = '[[species]]'//nl//'symbol = "Ni"'//nl//'z = 28'//nl//'nval = 10'//nl// &
    'konfig = [4, 4, 3, 4]'//nl//'occupations = [2.0, 8.0]'//nl
  character(len=*), parameter :: sp_al = '[[species]]'//nl//'symbol = "Al"'//nl//'z = 13'//nl//'nval = 3'//nl// &
    'konfig = [3, 3, 3, 4]'//nl//'occupations = [2.0, 0.0]'//nl
  character(len=*), parameter :: sp_cu = '[[species]]'//nl//'symbol = "Cu"'//nl//'z = 29'//nl//'nval = 11'//nl// &
    'konfig = [4, 4, 3, 4]'//nl//'occupations = [1.0, 10.0]'//nl
  character(len=*), parameter :: sp_zn = '[[species]]'//nl//'symbol = "Zn"'//nl//'z = 30'//nl//'nval = 12'//nl// &
    'konfig = [4, 4, 3, 4]'//nl//'occupations = [2.0, 10.0]'//nl
  character(len=*), parameter :: feco = '[structure]'//nl//'lattice = "bcc"'//nl//'alat = 2.84'//nl//'centering = "I"'//nl// &
    '[[sublattice]]'//nl//'name = "bcc"'//nl//'positions = [[0.0, 0.0, 0.0]]'//nl//'species = ["Fe", "Co"]'//nl// &
    'concentrations = [0.6, 0.4]'//nl//sp_fe//sp_co//common//'spin = 2'//nl//asr2
  character(len=*), parameter :: heusler = '[structure]'//nl//'lattice = "bcc"'//nl//'alat = 2.8635'//nl//'period = 2'//nl// &
    'centering = "F"'//nl// &
    '[[sublattice]]'//nl//'name = "8c"'//nl//'positions = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]'//nl//'species = ["Co"]'//nl// &
    '[[sublattice]]'//nl//'name = "4a"'//nl//'positions = [[0.0, 0.0, 0.0]]'//nl//'species = ["Al"]'//nl// &
    '[[sublattice]]'//nl//'name = "4b"'//nl//'positions = [[0.5, 0.5, 0.5]]'//nl//'species = ["Fe", "Cr"]'//nl// &
    'concentrations = [0.4, 0.6]'//nl//sp_co//sp_fe//sp_cr//sp_al//common//'spin = 2'//nl//asr2
  character(len=*), parameter :: fenicr = '[structure]'//nl//'lattice = "fcc"'//nl//'alat = 3.59'//nl//'centering = "F"'//nl// &
    '[[sublattice]]'//nl//'name = "fcc"'//nl//'positions = [[0.0, 0.0, 0.0]]'//nl//'species = ["Fe", "Ni", "Cr"]'//nl// &
    'concentrations = [0.66, 0.14, 0.20]'//nl// &
    '[[sro]]'//nl//'pairs = [["Fe", "Ni", 0.3], ["Fe", "Cr", 0.2], ["Ni", "Cr", -0.1]]'//nl// &
    sp_fe//sp_ni//sp_cr//common//'spin = 1'//nl//asr2
  character(len=*), parameter :: cuzn = '[structure]'//nl//'lattice = "bcc"'//nl//'alat = 2.95'//nl//'centering = "P"'//nl// &
    '[[sublattice]]'//nl//'name = "a"'//nl//'positions = [[0.0, 0.0, 0.0]]'//nl//'species = ["Cu", "Zn"]'//nl// &
    'concentrations = [0.8, 0.2]'//nl// &
    '[[sublattice]]'//nl//'name = "b"'//nl//'positions = [[0.5, 0.5, 0.5]]'//nl//'species = ["Zn", "Cu"]'//nl// &
    'concentrations = [0.75, 0.25]'//nl//sp_cu//sp_zn//common//'spin = 1'//nl//asr2//'route = "ed3"'//nl
  logical :: prt
  character(len=16) :: arg
  prt = .false.
  if (command_argument_count() >= 1) then; call get_command_argument(1, arg); prt = trim(arg) == 'REF=print'; end if
  print '(a)', 'test_cluster_scf'
  call run_case('bcc FeCo', feco, [652], -3.6974052776241245e-01_dp, -2.6382141589607245e+03_dp, &
                [2.2615112086433653e+00_dp, 1.6990970325126857e+00_dp], [6.5606624146262504e-02_dp, -9.8410304179537889e-02_dp])
  call run_case('L2_1 Co2Fe0.4Cr0.6Al', heusler, [348, 229, 182], -3.7009722269919343e-01_dp, -2.0812638865115482e+03_dp, &
                [1.0423340363627753e+00_dp, -6.0680807970516101e-02_dp, 2.4235386093174869e+00_dp, 1.9663162734337110e+00_dp], &
                [5.2878261797863502e-01_dp, -1.0034240898547910e+00_dp, -3.0830477613642060e-01_dp, 1.1530036484840522e-01_dp])
  call run_case('fcc FeNiCr SRO', fenicr, [2143], -3.7769246264264722e-01_dp, -2.5222479263826299e+03_dp, &
                [0.0_dp, 0.0_dp, 0.0_dp], &
                [9.8330939188631561e-02_dp, 6.7214038138793342e-01_dp, -7.9499105916784174e-01_dp])
  call run_case('B2 CuZn ed3', cuzn, [652, 652], -1.6865670580707279e-01_dp, -3.4397395544299511e+03_dp, &
                [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
                [2.0783556868894948e-01_dp, -1.8163836123302701e-01_dp, -2.2353607061289793e-01_dp, 1.5084347073269555e-01_dp])
  if (.not. prt) call run_nias_pipeline()
  call report()
contains

  subroutine run_nias_pipeline()
    character(len=*), parameter :: fn = 'build/test_nias_scf_tmp.toml'
    character(len=*), parameter :: nias_txt = &
      '[calculation]'//nl//'name = "NiAs pipeline test"'//nl// &
      '[structure]'//nl//'lattice = "nias"'//nl//'alat = 3.90'//nl//'c_over_a = 1.62'//nl// &
      sp_co//sp_cu//sp_zn// &
      '[[sublattice]]'//nl//'name = "cat"'//nl//'species = ["Co"]'//nl// &
      '[[sublattice]]'//nl//'name = "an"'//nl//'species = ["Cu", "Zn"]'//nl//'concentrations = [0.6, 0.4]'//nl// &
      '[[sro]]'//nl//'sublattices = ["an", "an"]'//nl//'pairs = [["Cu", "Zn", -0.15]]'//nl// &
      '[scf]'//nl//'spin = 1'//nl//'xc = "vbh"'//nl//'mixing = 0.2'//nl//'max_iterations = 1'//nl//'tolerance = 1.0e3'//nl// &
      '[asr]'//nl//'levels = 2'//nl//'recursion_steps = 6'//nl//'contour_points = 8'//nl//'fermi_kt = 0.008'//nl// &
      'cluster_radius = 1.05'//nl// &
      '[spectral]'//nl//'kpath = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]'//nl//'labels = ["G", "M"]'//nl// &
      'nk = 2'//nl//'npoints = 5'//nl//'rmax = 1.1'//nl//'file = "test_nias_spec.dat"'//nl// &
      '[epi]'//nl//'shells = 2'//nl//'lambda_points = 1'//nl//'file = "test_nias_epi.toml"'//nl// &
      '[dos]'//nl//'npoints = 3'//nl// &
      '[output]'//nl//'dos_file = "test_nias_dos.dat"'//nl//'out_file = "test_nias.out"'//nl
    type(run_input_t) :: inp
    type(class_scf_result_t) :: r
    character(len=:), allocatable :: msg
    character(len=256) :: line
    integer :: u, ios
    logical :: ok, all_le_90, has_vec_a
    open (newunit=u, file=fn, status='replace'); write (u, '(a)') nias_txt; close (u)
    call read_run_input(fn, inp, msg)
    call check(len(msg) == 0, 'NiAs pipeline: input accepted')
    if (len(msg) > 0) then; print '(2x,a)', msg; return; end if
    call run_input(inp, ok, r)
    call check(r%ncomp == 3 .and. size(r%nodes) == 2, 'NiAs pipeline: 3 components, 2 sublattices')
    call check(r%screen%z == 6, 'NiAs pipeline: intra-sublattice anion SRO uses z = 6 shell')
    call check(allocated(r%spec_a) .and. all(r%spec_a >= 0.0_dp), 'NiAs pipeline: ball-window spectral function non-negative')
    call check(allocated(r%epi) .and. size(r%epi) == 1, 'NiAs pipeline: one anion-anion EPI pair type in first two shells')
    if (allocated(r%epi) .and. size(r%epi) == 1) then
      call check(r%epi(1)%kap == 2 .and. r%epi(1)%kap2 == 2 .and. r%epi(1)%mult == 6 .and. abs(r%epi(1)%uband(2, 2)) > 1.0e-6_dp, &
                 'NiAs pipeline: non-zero anion-anion EPI on shell 2 (6 sites)')
    end if
    open (newunit=u, file='build/test_nias.out', status='old', action='read', iostat=ios)
    all_le_90 = (ios == 0)
    do while (ios == 0)
      read (u, '(a)', iostat=ios) line
      if (ios == 0 .and. len_trim(line) > 90) all_le_90 = .false.
    end do
    if (all_le_90 .or. ios /= 0) close (u)
    call check(all_le_90, 'NiAs pipeline: all lines of report <= 90 columns')
    open (newunit=u, file='build/test_nias_epi.toml', status='old', action='read', iostat=ios)
    has_vec_a = .false.
    do while (ios == 0)
      read (u, '(a)', iostat=ios) line
      if (ios == 0 .and. index(line, 'vector_alat =') > 0) has_vec_a = .true.
    end do
    close (u)
    call check(has_vec_a, 'NiAs pipeline: pair file writes vector_alat for non-cubic pair')
  end subroutine run_nias_pipeline

  subroutine run_case(tag, txt, nodes, ef, etot, mom, qnet)
    character(len=*), intent(in) :: tag, txt
    integer, intent(in) :: nodes(:)
    real(dp), intent(in) :: ef, etot, mom(:), qnet(:)
    type(run_input_t) :: inp
    type(class_scf_result_t) :: r
    character(len=:), allocatable :: msg
    character(len=*), parameter :: f = 'build/test_cluster_scf_tmp.toml'
    integer :: u, p
    open (newunit=u, file=f, status='replace'); write (u, '(a)') txt; close (u)
    call read_run_input(f, inp, msg)
    call check(len(msg) == 0, tag//': input accepted')
    if (len(msg) > 0) then; print '(2x,a)', msg; return; end if
    call class_species(inp)
    inp%ci%ndos = 0
    call class_scf(inp%ci, r)
    if (prt) then
      print '(2x,a)', tag
      print '(4x,a,*(i0,:,", "))', 'nodes [', r%nodes; print '(4x,a)', ']'
      print '(4x,a,es24.16)', 'ef   ', r%ef
      print '(4x,a,es24.16)', 'etot ', r%etot
      print '(4x,a,*(es24.16,:,", "))', 'mom  ', (r%comp(p)%moment, p=1, r%ncomp)
      print '(4x,a,*(es24.16,:,", "))', 'qnet ', (r%comp(p)%qnet, p=1, r%ncomp)
      return
    end if
    print '(2x,a,a,*(i8))', tag, ': nodes per sublattice', r%nodes
    print '(2x,a,a,2es11.3,a,2es11.3)', tag, ': E_F, etot - reference', r%ef - ef, r%etot - etot, &
      ';  max |moment, qnet - reference|', maxval(abs([(r%comp(p)%moment, p=1, r%ncomp)] - mom)), &
      maxval(abs([(r%comp(p)%qnet, p=1, r%ncomp)] - qnet))
    call check(size(r%nodes) == size(nodes), tag//': one operator per sublattice')
    if (size(r%nodes) == size(nodes)) call check(all(r%nodes == nodes), tag//': node counts of the box that does not wrap')
    call check(abs(r%ef - ef) < 1.0e-8_dp, tag//': E_F of the box code (1e-8 Ry)')
    call check(abs(r%etot - etot) < 1.0e-7_dp, tag//': total energy of the box code (1e-7 Ry)')
    call check(r%ncomp == size(mom), tag//': components')
    if (r%ncomp /= size(mom)) return
    call check(maxval(abs([(r%comp(p)%moment, p=1, r%ncomp)] - mom)) < 1.0e-7_dp, tag//': moments of the box code (1e-7)')
    call check(maxval(abs([(r%comp(p)%qnet, p=1, r%ncomp)] - qnet)) < 1.0e-7_dp, tag//': net charges of the box code (1e-7)')
  end subroutine run_case
end program test_cluster_scf
