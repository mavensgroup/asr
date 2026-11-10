!> The valence band energy of the k-space ASA total energy must be the occupied eigenvalue sum of the bands it solves.
!> For the generalised problem H z = e O z the identity is, per eigenvector and summed over all sites and orbitals,
!>   e z^T O z = z^T H z = sum_RL [ q1 + E_nu (q0 + p^gamma q2) ],
!> with q0 + p^gamma q2 the sphere charge of the channel.  For the orthogonal second-order Hamiltonian the moments are energy
!> moments and the sum is q1 + E_nu q0.  The identity holds at every iteration and for any occupation weights, so short
!> unconverged runs test it.
program test_band_energy
  use asr_kinds, only: dp, pi
  use asr_asa
  use asr_atom, only: shell_t
  use asr_lmto_cell
  use asr_testing
  implicit none
  type(asa_input_t) :: inp
  type(asa_result_t) :: res
  type(cell_t) :: cell
  real(dp) :: a_bohr, ratio, sfe, scu
  integer :: ic

  print '(a)', 'test_band_energy'
  ratio = 2.91398733_dp / 2.65686809_dp
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
  inp%alat = a_bohr; inp%nk = 6; inp%smear = 0.005_dp; inp%beta = 0.3_dp; inp%maxit = 3; inp%tol = 1.0e-6_dp
  cell%nat = 2; cell%avec = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
  allocate (cell%tau(3, 2), cell%isp(2)); cell%tau(:, 1) = 0.0_dp; cell%tau(:, 2) = 0.5_dp; cell%isp = [1, 2]
  call setup_cell(cell)

  inp%nspin = 1; inp%bz = 'tetra'
  call asa_scf(inp, cell, res)
  print '(a,2f16.10)', '    B2 FeCu, tetrahedra: band energy from moments, eigenvalue sum (Ry/cell) ', res%eband, res%sumev
  call check_close(res%eband, res%sumev, 1.0e-9_dp, 'full ASA, tetrahedra: band energy == occupied eigenvalue sum')

  inp%bz = 'gauss'
  call asa_scf(inp, cell, res)
  call check_close(res%eband, res%sumev, 1.0e-9_dp, 'full ASA, Gaussian weights: band energy == occupied eigenvalue sum')

  inp%nspin = 2; inp%zeta0 = 0.2_dp; inp%bz = 'tetra'
  call asa_scf(inp, cell, res)
  call check_close(res%eband, res%sumev, 1.0e-9_dp, 'full ASA, spin polarised: band energy == occupied eigenvalue sum')

  inp%nspin = 1; inp%zeta0 = 0.0_dp; inp%second_order = .true.
  call asa_scf(inp, cell, res)
  call check_close(res%eband, res%sumev, 1.0e-9_dp, 'second-order H: band energy == occupied eigenvalue sum')
  call report()
end program test_band_energy
