!> Ordered B2 FeCu (nonmagnetic) generated from scratch: two classes, Ewald Madelung, charge transfer.
!> The Stuttgart converged moments are an ORACLE for comparison only.  The only geometry taken from that
!> run is the Fe:Cu sphere-radius ratio (so the comparison is of the same cell); spheres fill the volume.
program test_asa_b2
  use asr_kinds, only: dp, pi
  use asr_asa
  use asr_atom, only: shell_t
  use asr_lmto_cell
  use asr_testing
  implicit none
  type(asa_input_t) :: inp
  type(asa_result_t) :: res
  type(cell_t) :: cell
  real(dp) :: a_bohr, ratio, sfe, scu, occ(0:2, 2)
  integer :: ic, l
  real(dp), parameter :: ofe(0:2) = [0.65065469_dp, 0.77696927_dp, 6.76341845_dp]
  real(dp), parameter :: ocu(0:2) = [0.70640914_dp, 0.67687944_dp, 9.42566900_dp]
  real(dp), parameter :: oef = -0.18713836_dp

  print '(a)', 'test_asa_b2'
  ratio = 2.91398733_dp / 2.65686809_dp            ! oracle sphere-radius ratio (geometry only)
  a_bohr = 3.0_dp / 0.529177_dp
  scu = (a_bohr**3 / (4.0_dp * pi / 3.0_dp * (1.0_dp + ratio**3)))**(1.0_dp / 3.0_dp); sfe = ratio * scu
  allocate (inp%cls(2))
  inp%cls(1)%z = 26.0_dp; inp%cls(1)%nval = 8; inp%cls(1)%sws = sfe
  inp%cls(2)%z = 29.0_dp; inp%cls(2)%nval = 11; inp%cls(2)%sws = scu
  do ic = 1, 2
    inp%cls(ic)%core = core_shells_for([4, 4, 3, 4]); inp%cls(ic)%nvl = [4, 4, 3]
  end do
  inp%cls(1)%atom_shells = [inp%cls(1)%core, shell_t(4, 0, 2.0_dp), shell_t(3, 2, 6.0_dp)]
  inp%cls(2)%atom_shells = [inp%cls(2)%core, shell_t(4, 0, 1.0_dp), shell_t(3, 2, 10.0_dp)]
  inp%alat = a_bohr; inp%nspin = 1; inp%nk = 12; inp%smear = 0.005_dp; inp%beta = 0.3_dp; inp%maxit = 100; inp%tol = 1.0e-6_dp
  cell%nat = 2; cell%avec = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
  allocate (cell%tau(3, 2), cell%isp(2)); cell%tau(:, 1) = 0.0_dp; cell%tau(:, 2) = 0.5_dp; cell%isp = [1, 2]
  call setup_cell(cell)
  call asa_scf(inp, cell, res)
  call check(res%converged, 'B2 FeCu SCF converged')
  do ic = 1, 2
    do l = 0, 2
      occ(l, ic) = res%cls(ic)%q(0, l, 1) + res%cls(ic)%pp(l, 1)%p * res%cls(ic)%q(2, l, 1)
    end do
  end do
  print '(a,3f8.4,a,3f8.4)', '    Fe s,p,d (ours) ', occ(:, 1), '   oracle ', ofe
  print '(a,3f8.4,a,3f8.4)', '    Cu s,p,d (ours) ', occ(:, 2), '   oracle ', ocu
  print '(a,f9.5,a,f9.5,a,f8.4)', '    E_F ours ', res%ef, '  oracle ', oef, '   charge transfer (Fe) ', res%cls(1)%qnet
  call check_close(res%cls(1)%qnet + res%cls(2)%qnet, 0.0_dp, 1.0e-6_dp, 'cell is neutral (sum of sphere charges)')
  call check_close(occ(2, 1), ofe(2), 0.05_dp, 'Fe d occupation vs Stuttgart 6.763')
  call check_close(occ(2, 2), ocu(2), 0.03_dp, 'Cu d occupation vs Stuttgart 9.426')
  call check_close(occ(0, 1), ofe(0), 0.03_dp, 'Fe s occupation vs Stuttgart 0.651')
  call check_close(occ(1, 1), ofe(1), 0.03_dp, 'Fe p occupation vs Stuttgart 0.777')
  call check_close(occ(0, 2), ocu(0), 0.03_dp, 'Cu s occupation vs Stuttgart 0.706')
  call check_close(occ(1, 2), ocu(1), 0.03_dp, 'Cu p occupation vs Stuttgart 0.677')
  call check_close(res%ef, oef, 0.012_dp, 'E_F within 12 mRy of Stuttgart')
  call report()
end program test_asa_b2
