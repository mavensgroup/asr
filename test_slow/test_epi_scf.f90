!> Effective pair interactions through the production routine (asr_class_scf::epi_run) and the driver.
!> fcc Cu-Ni-Ni' with Ni' a second label of Ni (the same atom): the pair table must not see the label.  Exact statements:
!>   U(c, 1) = U(1, s) = 0 (reference species);  U(2,2) = U(2,3) = U(3,2) = U(3,3);  V(Ni, Ni') = 0;
!>   the same for the screened Coulomb part.
!> A wrong species index, start vector or frozen site in the chain loop breaks them.  The run goes through the driver, which
!> writes the report and the pair file (build/test_epi_scf_epi.toml).  Small settings (two SCF iterations, levels 2, first
!> neighbours only, one path point): the numbers are not converged.
program test_epi_scf
  use asr_kinds, only: dp
  use asr_driver, only: run_input_t, read_run_input, run_input
  use asr_class_scf, only: class_scf_result_t, wbelow_tol
  use asr_epi, only: epi_exchange
  use asr_testing
  implicit none
  type(run_input_t) :: inp
  type(class_scf_result_t) :: r
  character(len=:), allocatable :: msg
  real(dp), allocatable :: v(:, :)
  real(dp) :: u22
  integer :: u, ios
  logical :: ok, ex
  character(len=256) :: line
  character(len=*), parameter :: f = 'build/test_epi_scf_tmp.toml'
  character(len=1), parameter :: nl = new_line('a')
  character(len=*), parameter :: txt = '[calculation]'//nl//'name = "epi test"'//nl// &
    '[structure]'//nl//'lattice = "fcc"'//nl//'alat = 3.57'//nl//'centering = "F"'//nl// &
    '[[species]]'//nl//'symbol = "Cu"'//nl//'z = 29'//nl//'nval = 11'//nl//'konfig = [4, 4, 3, 4]'//nl//'occupations = [1.0, 10.0]'//nl// &
    '[[species]]'//nl//'symbol = "Ni"'//nl//'z = 28'//nl//'nval = 10'//nl//'konfig = [4, 4, 3, 4]'//nl//'occupations = [2.0, 8.0]'//nl// &
    '[[species]]'//nl//'symbol = "Nx"'//nl//'z = 28'//nl//'nval = 10'//nl//'konfig = [4, 4, 3, 4]'//nl//'occupations = [2.0, 8.0]'//nl// &
    '[[sublattice]]'//nl//'name = "fcc"'//nl//'positions = [[0.0, 0.0, 0.0]]'//nl//'species = ["Cu", "Ni", "Nx"]'//nl// &
    'concentrations = [0.5, 0.3, 0.2]'//nl// &
    '[scf]'//nl//'spin = 1'//nl//'xc = "vbh"'//nl//'mixing = 0.3'//nl//'max_iterations = 2'//nl//'tolerance = 1.0e3'//nl// &
    '[asr]'//nl//'levels = 2'//nl//'recursion_steps = 6'//nl//'supercell = 4'//nl//'contour_points = 8'//nl//'fermi_kt = 0.008'//nl// &
    'cluster_radius = 0.8'//nl//'flip_radius = 2'//nl// &
    '[epi]'//nl//'shells = 2'//nl//'lambda_points = 1'//nl//'file = "test_epi_scf_epi.toml"'//nl// &
    '[dos]'//nl//'npoints = 3'//nl//'[output]'//nl//'dos_file = "test_epi_scf_dos.dat"'//nl//'out_file = "test_epi_scf.out"'//nl

  print '(a)', 'test_epi_scf'
  open (newunit=u, file=f, status='replace'); write (u, '(a)') txt; close (u)
  call read_run_input(f, inp, msg)
  call check(len(msg) == 0, 'fcc Cu-Ni-Ni'' input with [epi] accepted')
  if (len(msg) > 0) then; print '(a)', msg; call report(); stop; end if
  call run_input(inp, ok, r)
  call check(allocated(r%epi), 'pair tables returned')
  if (.not. allocated(r%epi)) then; call report(); stop; end if
  call check(size(r%epi) == 2, 'two pair types: first and second shell of fcc')
  if (size(r%epi) /= 2) then; call report(); stop; end if
  call check(all(r%epi%mult == [12, 6]) .and. all(r%epi%nop == [4, 4]) .and. all(r%epi%kap == 1) .and. all(r%epi%kap2 == 1), &
             'fcc shells: 12 and 6 sites, pair groups of order 4')
  call check(r%epi_wbelow < wbelow_tol, 'no weight of a path site below the contour bottom')
  associate (ub => r%epi(1)%uband, ue => r%epi(1)%ues)
    print '(a,4es16.8)', '  shell 1: U band (2,2) (2,3) (3,2) (3,3) (Ry) ', ub(2, 2), ub(2, 3), ub(3, 2), ub(3, 3)
    print '(a,4es16.8)', '  shell 1: U Coulomb                            ', ue(2, 2), ue(2, 3), ue(3, 2), ue(3, 3)
    u22 = ub(2, 2)
    call check(all(ub(1, :) == 0.0_dp) .and. all(ub(:, 1) == 0.0_dp) .and. all(ue(1, :) == 0.0_dp) .and. all(ue(:, 1) == 0.0_dp), &
               'the rows and columns of the reference species vanish')
    call check(abs(u22) > 1.0e-5_dp, 'the Cu-Ni pair interaction is not trivially small')
    call check_close(maxval(abs(ub(2:3, 2:3) - u22)) / abs(u22), 0.0_dp, 1.0e-7_dp, 'band: U(2,2) = U(2,3) = U(3,2) = U(3,3) for the relabelled species')
    call check_close(maxval(abs(ue(2:3, 2:3) - ue(2, 2))) / abs(ue(2, 2)), 0.0_dp, 1.0e-7_dp, 'screened Coulomb: the same')
    v = epi_exchange(ub + ue)
    call check_close(abs(v(2, 3)) / abs(v(1, 2)), 0.0_dp, 1.0e-7_dp, 'V(Ni, Ni'') = 0')
    call check_close(v(1, 3), v(1, 2), 1.0e-7_dp * abs(v(1, 2)), 'V(Cu, Ni) = V(Cu, Ni'')')
    print '(a,3f10.4)', '  shell 1: V(Cu,Ni) band, Coulomb, total (mRy; settings not converged) ', 1.0e3_dp * u22, 1.0e3_dp * ue(2, 2), 1.0e3_dp * v(1, 2)
    call check(ue(2, 2) > 0.0_dp, 'the screened Coulomb interaction favours unlike neighbours')
  end associate
  call check(abs(r%epi(2)%ues(2, 2)) < abs(r%epi(1)%ues(2, 2)), 'second shell: weaker screened Coulomb interaction')

  ! the driver: pair file
  inquire (file='build/test_epi_scf_epi.toml', exist=ex)
  call check(ex, 'the driver wrote the pair file')
  if (ex) then
    open (newunit=u, file='build/test_epi_scf_epi.toml', status='old', action='read')
    ios = 0; ok = .false.
    do while (ios == 0)
      read (u, '(a)', iostat=ios) line
      if (ios == 0 .and. index(line, '[[pair]]') > 0) ok = .true.
    end do
    close (u)
    call check(ok, 'the pair file holds [[pair]] tables')
  end if
  call report()
end program test_epi_scf
