!> The one input format (asr_driver): method selection, strict keys, the removed kind key, and the k-space cell of an ordered
!> decorated structure (primitive reduction by the centering translations, one class per sublattice, B2 sphere ratio).
program test_run_input
  use asr_kinds, only: dp
  use asr_driver
  use asr_testing
  implicit none
  type(run_input_t) :: inp
  character(len=:), allocatable :: msg
  integer :: u
  character(len=*), parameter :: f = 'build/test_run_input_tmp.toml'
  character(len=1), parameter :: nl = new_line('a')
  character(len=*), parameter :: fe = '[[species]]'//nl//'symbol = "Fe"'//nl//'z = 26'//nl//'nval = 8'//nl// &
                                      'occupations = [2.0, 6.0]'//nl
  character(len=*), parameter :: co = '[[species]]'//nl//'symbol = "Co"'//nl//'z = 27'//nl//'nval = 9'//nl// &
                                      'occupations = [2.0, 7.0]'//nl
  character(len=*), parameter :: al = '[[species]]'//nl//'symbol = "Al"'//nl//'z = 13'//nl//'nval = 3'//nl// &
                                      'konfig = [3, 3, 3, 4]'//nl//'occupations = [2.0, 0.0]'//nl
  character(len=*), parameter :: bccfe = '[structure]'//nl//'lattice = "bcc"'//nl//'alat = 2.8665'//nl//'centering = "I"'//nl// &
                                         fe//'[[sublattice]]'//nl//'name = "a"'//nl//'positions = [[0.0, 0.0, 0.0]]'//nl// &
                                         'species = ["Fe"]'//nl
  character(len=*), parameter :: b2 = '[structure]'//nl//'lattice = "bcc"'//nl//'alat = 2.857'//nl//fe//co// &
                                      '[[sublattice]]'//nl//'name = "a"'//nl//'positions = [[0.0, 0.0, 0.0]]'//nl//'species = ["Co"]'//nl// &
                                      '[[sublattice]]'//nl//'name = "b"'//nl//'positions = [[0.5, 0.5, 0.5]]'//nl//'species = ["Fe"]'//nl
  character(len=*), parameter :: l21 = '[structure]'//nl//'lattice = "bcc"'//nl//'alat = 2.8635'//nl//'period = 2'//nl// &
                                       'centering = "F"'//nl//co//fe//al// &
                                       '[[sublattice]]'//nl//'name = "8c"'//nl//'positions = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]'//nl// &
                                       'species = ["Co"]'//nl// &
                                       '[[sublattice]]'//nl//'name = "4a"'//nl//'positions = [[0.0, 0.0, 0.0]]'//nl//'species = ["Al"]'//nl// &
                                       '[[sublattice]]'//nl//'name = "4b"'//nl//'positions = [[0.5, 0.5, 0.5]]'//nl//'species = ["Fe"]'//nl
  character(len=*), parameter :: kspace = '[calculation]'//nl//'method = "kspace"'//nl

  print '(a)', 'test_run_input'
  call put(bccfe)
  call read_run_input(f, inp, msg)
  call check(len(msg) == 0, 'pure bcc Fe accepted')
  call check(inp%method == 'asr' .and. inp%ci%dec%ncls == 1 .and. inp%ci%dec%nsite == 2, 'method defaults to asr; bcc graph, 2 sites')
  call check(inp%ci%tol == 1.0e-5_dp .and. inp%ci%nlevels == 4 .and. inp%ci%ed_tail == 'bp', 'asr defaults filled')

  call put('[calculation]'//nl//'kind = "lmto"'//nl//bccfe)
  call read_run_input(f, inp, msg)
  call check(index(msg, 'kind is no longer read') > 0, 'old kind key rejected with a pointer to method')
  call put(bccfe//'[asr]'//nl//'levls = 3'//nl)
  call read_run_input(f, inp, msg)
  call check(index(msg, 'unknown key') > 0, 'typo in [asr] rejected')
  call put(kspace//bccfe//'[asr]'//nl//'levels = 3'//nl)
  call read_run_input(f, inp, msg)
  call check(index(msg, '[asr] is read with method = "asr" only') > 0, '[asr] with method = kspace rejected')
  call put(bccfe//'[kspace]'//nl//'nk = 12'//nl)
  call read_run_input(f, inp, msg)
  call check(index(msg, '[kspace] is read') > 0, '[kspace] with method = asr rejected')
  call put('[structure]'//nl//'lattice = "cell"'//nl//'alat = 3.0'//nl//'vectors = [[1,0,0],[0,1,0],[0,0,1]]'//nl//fe// &
           '[[site]]'//nl//'position = [0.0, 0.0, 0.0]'//nl//'species = "Fe"'//nl)
  call read_run_input(f, inp, msg)
  call check(index(msg, 'method = "kspace" only') > 0, 'explicit cell needs method = kspace')
  call put(kspace//'[structure]'//nl//'lattice = "cell"'//nl//'alat = 3.0'//nl//'site_classes = true'//nl// &
           'vectors = [[1,0,0],[0,1,0],[0,0,1]]'//nl//fe//'[[site]]'//nl//'position = [0.0, 0.0, 0.0]'//nl//'species = "Fe"'//nl)
  call read_run_input(f, inp, msg)
  call check(len(msg) == 0 .and. inp%li%lattice == 'cell' .and. inp%li%site_classes .and. size(inp%li%site_sp) == 1, &
             'explicit cell with method = kspace accepted')
  call put(kspace//'[structure]'//nl//'lattice = "bcc"'//nl//'alat = 2.85'//nl//'centering = "I"'//nl//fe//co// &
           '[[sublattice]]'//nl//'positions = [[0.0, 0.0, 0.0]]'//nl//'species = ["Fe", "Co"]'//nl//'concentrations = [0.5, 0.5]'//nl)
  call read_run_input(f, inp, msg)
  call check(index(msg, 'ordered structure') > 0, 'disorder with method = kspace rejected')

  ! k-space cells of ordered decorated structures
  call put(kspace//bccfe)
  call read_run_input(f, inp, msg)
  call check(len(msg) == 0 .and. size(inp%li%site_sp) == 1 .and. abs(inp%li%cell_vec(1, 1) + 0.5_dp) < 1.0e-12_dp, &
             'bcc Fe -> one-atom bcc primitive cell')
  call check(inp%li%tolerance == 1.0e-6_dp .and. inp%li%nk == 12 .and. inp%li%bz == 'tetra', 'kspace defaults filled')
  call put(kspace//'[structure]'//nl//'lattice = "fcc"'//nl//'alat = 4.09'//nl//'centering = "F"'//nl//fe// &
           '[[sublattice]]'//nl//'positions = [[0.0, 0.0, 0.0]]'//nl//'species = ["Fe"]'//nl)
  call read_run_input(f, inp, msg)
  call check(len(msg) == 0 .and. size(inp%li%site_sp) == 1 .and. abs(inp%li%cell_vec(1, 1)) < 1.0e-12_dp, &
             'fcc -> one-atom fcc primitive cell')
  call put(kspace//b2)
  call read_run_input(f, inp, msg)
  call check(len(msg) == 0 .and. inp%li%lattice == 'cell' .and. size(inp%li%site_sp) == 2, 'B2 -> two-atom simple-cubic cell')
  call put(kspace//b2//'[kspace]'//nl//'sphere_ratio = 1.1'//nl)
  call read_run_input(f, inp, msg)
  call check(len(msg) == 0 .and. inp%li%lattice == 'b2' .and. inp%li%ratio == 1.1_dp .and. inp%li%species(1)%symbol == 'Co', &
             'B2 sphere_ratio -> unequal spheres, origin sublattice first')
  call put(kspace//l21)
  call read_run_input(f, inp, msg)
  call check(len(msg) == 0 .and. size(inp%li%site_sp) == 4 .and. count(inp%li%site_sp == 1) == 2 .and. &
             abs(inp%li%cell_vec(1, 1)) < 1.0e-12_dp .and. abs(inp%li%cell_vec(2, 1) - 1.0_dp) < 1.0e-12_dp, &
             'L2_1 -> four-atom fcc primitive cell (period 2), Co twice')
  call put(kspace//l21//'[kspace]'//nl//'sphere_ratio = 1.1'//nl)
  call read_run_input(f, inp, msg)
  call check(index(msg, 'B2') > 0, 'sphere_ratio outside B2 rejected')

  ! every shipped example is a valid input (fpm test runs from the project root)
  call example('ag_fcc.toml', 'kspace'); call example('fe_bcc.toml', 'kspace'); call example('feco_b2_kspace.toml', 'kspace')
  call example('fecu_b2.toml', 'kspace'); call example('feco_b2_asr.toml', 'asr'); call example('cuzn_b2.toml', 'asr')
  call example('cuzn_alloy.toml', 'asr'); call example('co2fe04cr06al_l21.toml', 'asr')
  call example('fe66ni14cr20_fcc.toml', 'asr')
  call report()
contains
  subroutine example(name, method)
    character(len=*), intent(in) :: name, method
    call read_run_input('example/'//name, inp, msg)
    if (len(msg) > 0) print '(a)', '        '//msg
    call check(len(msg) == 0 .and. inp%method == method, 'example/'//name//' reads (method = '//method//')')
  end subroutine example

  subroutine put(text)
    character(len=*), intent(in) :: text
    open (newunit=u, file=f, status='replace'); write (u, '(a)') text; close (u)
  end subroutine put
end program test_run_input
