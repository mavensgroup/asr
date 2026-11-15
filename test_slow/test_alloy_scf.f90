!> Alloy SCF from the augmented-space recursion, trivial alloy A = B = Ag (no disorder): must reproduce the ordered
!> ASA SCF (asa_scf).  Exercises orthogonal-representation moments, species-resolved Green functions, contour
!> integration, the Fermi-level search, relaxed core, Hartree/xc/SIM and Anderson mixing.
program test_alloy_scf
  use asr_kinds, only: dp
  use asr_alloy_scf
  use asr_asa
  use asr_atom, only: shell_t
  use asr_lmto_cell
  use asr_testing
  implicit none
  type(alloy_input_t) :: ai
  type(alloy_result_t) :: ar, ae
  type(asa_input_t) :: oi
  type(asa_result_t) :: or
  type(cell_t) :: cell
  integer :: s, l

  print '(a)', 'test_alloy_scf'
  ! ordered reference (first-order ASA)
  allocate (oi%cls(1))
  oi%cls(1)%z = 47.0_dp; oi%cls(1)%nval = 11; oi%cls(1)%core = core_shells_for([5, 5, 4, 4]); oi%cls(1)%nvl = [5, 5, 4]
  oi%cls(1)%atom_shells = [oi%cls(1)%core, shell_t(5, 0, 1.0_dp), shell_t(4, 2, 10.0_dp)]
  oi%alat = 7.72898293_dp; oi%nk = 12; oi%maxit = 80
  cell%nat = 1; cell%avec = reshape([0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp], [3, 3])
  allocate (cell%tau(3, 1), cell%isp(1)); cell%tau = 0.0_dp; cell%isp = 1
  call setup_cell(cell); oi%cls(1)%sws = cell%w * oi%alat
  call asa_scf(oi, cell, or)
  ! alloy with identical species
  ai%lattice = 'fcc'; ai%alat = 7.72898293_dp; ai%x = 0.5_dp; ai%nlevels = 8; ai%lsuper = 30; ai%maxit = 40; ai%tol = 1.0e-5_dp
  do s = 1, 2
    ai%spec(s)%z = 47.0_dp; ai%spec(s)%nval = 11; ai%spec(s)%core = core_shells_for([5, 5, 4, 4]); ai%spec(s)%nvl = [5, 5, 4]
    ai%spec(s)%atom_shells = [ai%spec(s)%core, shell_t(5, 0, 1.0_dp), shell_t(4, 2, 10.0_dp)]
  end do
  call alloy_scf(ai, ar)
  call check(or%converged, 'ordered ASA SCF converged')
  call check(ar%converged, 'alloy (A = B) SCF converged')
  do l = 0, 2
    print '(a,i1,a,3f10.5,a,3f10.5)', '    l=', l, '  alloy sqrdel gamma p: ', abs(ar%cls(1)%pp(l, 1)%sqrdel), ar%cls(1)%pp(l, 1)%gamma, ar%cls(1)%pp(l, 1)%p, &
      '   ordered: ', abs(or%cls(1)%pp(l, 1)%sqrdel), or%cls(1)%pp(l, 1)%gamma, or%cls(1)%pp(l, 1)%p
    call check_close(abs(ar%cls(1)%pp(l, 1)%sqrdel) / abs(or%cls(1)%pp(l, 1)%sqrdel), 1.0_dp, 0.02_dp, 'sqrt(Delta) within 2% of the ordered SCF')
    call check_close(ar%cls(1)%pp(l, 1)%p / or%cls(1)%pp(l, 1)%p, 1.0_dp, 0.02_dp, 'p within 2% of the ordered SCF')
  end do
  call check_close(ar%cls(1)%pp(0, 1)%gamma / or%cls(1)%pp(0, 1)%gamma, 1.0_dp, 0.01_dp, 'gamma_s within 1%')
  call check_close(ar%cls(1)%pp(1, 1)%gamma / or%cls(1)%pp(1, 1)%gamma, 1.0_dp, 0.01_dp, 'gamma_p within 1%')
  call check_close(ar%cls(1)%pp(2, 1)%enu, or%cls(1)%pp(2, 1)%enu, 0.01_dp, 'E_nu,d within 10 mRy of the ordered SCF')
  call check_close(ar%cls(1)%pp(2, 1)%c - ar%cls(1)%pp(2, 1)%enu, or%cls(1)%pp(2, 1)%c - or%cls(1)%pp(2, 1)%enu, 0.005_dp, 'C_d - E_nu,d within 5 mRy')
  call check_close(ar%ef, or%ef, 0.03_dp, 'E_F within 30 mRy (orthogonal vs first-order representation)')
  call check_close(ar%cls(1)%qnet, 0.0_dp, 1.0e-3_dp, 'sphere neutral')
  ! energy-dependent recursion (potential functions P(z) - S), same trivial alloy
  ai%energy_dependent = .true.; ai%lsuper = 12; ai%nlevels = 8
  call alloy_scf(ai, ae, verbose=.true.)
  call check(ae%converged, 'energy-dependent alloy (A = B) SCF converged')
  do l = 0, 2
    print '(a,i1,a,3f10.5,a,3f10.5)', '    ED l=', l, '  sqrdel gamma p: ', abs(ae%cls(1)%pp(l, 1)%sqrdel), ae%cls(1)%pp(l, 1)%gamma, ae%cls(1)%pp(l, 1)%p, &
      '   ordered: ', abs(or%cls(1)%pp(l, 1)%sqrdel), or%cls(1)%pp(l, 1)%gamma, or%cls(1)%pp(l, 1)%p
  end do
  print '(a,2f10.5)', '    ED E_F, ordered E_F: ', ae%ef, or%ef
  call check_close(ae%ef, or%ef, 0.03_dp, 'ED: E_F within 30 mRy of the first-order ordered ASA (8 levels)')
  call check_close(abs(ae%cls(1)%pp(2, 1)%sqrdel) / abs(or%cls(1)%pp(2, 1)%sqrdel), 1.0_dp, 0.05_dp, 'ED: sqrt(Delta_d) within 5%')
  call check_close(ae%cls(1)%qnet, 0.0_dp, 1.0e-3_dp, 'ED: sphere neutral')
  call report()
end program test_alloy_scf
