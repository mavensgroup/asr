!> Self-consistent ASA fcc Ag generated entirely by our own code (radial solver, libxc VBH,
!> relaxed core in the sphere, ASA moments).  The Stuttgart converged parameters appear ONLY as
!> an oracle for comparison; nothing from them is an input.
program test_asa_ag
  use asr_kinds, only: dp
  use asr_asa
  use asr_atom, only: shell_t
  use asr_lmto_cell
  use asr_testing
  implicit none
  type(asa_input_t) :: inp
  type(asa_result_t) :: res
  type(cell_t) :: cell
  integer :: l
  ! oracle (Stuttgart LMASA converged Ag) - comparison only
  real(dp), parameter :: o_sd(0:2) = [0.36457406_dp, 0.36517716_dp, 0.10341022_dp]
  real(dp), parameter :: o_g(0:2) = [0.42477526_dp, 0.11537965_dp, 0.00550351_dp]
  real(dp), parameter :: o_p(0:2) = [0.08230178_dp, 0.03623364_dp, 2.08144602_dp]
  real(dp), parameter :: o_enu(0:2) = [-0.54494317_dp, -0.45471939_dp, -0.51341327_dp]
  real(dp), parameter :: o_ef = -0.18641407_dp, o_vmtz = -0.70030984_dp

  print '(a)', 'test_asa_ag'
  allocate (inp%cls(1))
  inp%cls(1)%z = 47.0_dp; inp%cls(1)%nval = 11
  inp%cls(1)%core = core_shells_for([5, 5, 4, 4])
  inp%cls(1)%nvl = [5, 5, 4]
  inp%cls(1)%atom_shells = [inp%cls(1)%core, shell_t(5, 0, 1.0_dp), shell_t(4, 2, 10.0_dp)]
  inp%alat = 7.72898293_dp
  inp%xc = 'vbh'; inp%sigma = 0.7_dp
  inp%nk = 12; inp%smear = 0.01_dp; inp%beta = 0.3_dp; inp%maxit = 100; inp%tol = 1.0e-6_dp
  cell%nat = 1
  cell%avec = reshape([0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp], [3, 3])
  allocate (cell%tau(3, 1), cell%isp(1)); cell%tau = 0.0_dp; cell%isp = 1
  call setup_cell(cell)
  inp%cls(1)%sws = cell%w * inp%alat           ! equal-sphere cell (Ag experimental lattice constant)
  call asa_scf(inp, cell, res)
  call check(res%converged, 'Ag ASA SCF converged')
  do l = 0, 2
    print '(a,i1,a,3f9.5,a,3f9.5)', '    l=', l, '  sqrdel gamma p (ours): ', abs(res%cls(1)%pp(l, 1)%sqrdel), res%cls(1)%pp(l, 1)%gamma, res%cls(1)%pp(l, 1)%p, &
      '   oracle: ', o_sd(l), o_g(l), o_p(l)
    call check_close(abs(res%cls(1)%pp(l, 1)%sqrdel) / o_sd(l), 1.0_dp, 0.02_dp, 'sqrt(Delta) within 2% of Stuttgart')
    call check_close(res%cls(1)%pp(l, 1)%p / o_p(l), 1.0_dp, 0.02_dp, 'p within 2% of Stuttgart')
  end do
  call check_close(res%cls(1)%pp(0, 1)%gamma / o_g(0), 1.0_dp, 0.01_dp, 'gamma_s within 1%')
  call check_close(res%cls(1)%pp(1, 1)%gamma / o_g(1), 1.0_dp, 0.01_dp, 'gamma_p within 1%')
  call check_close(res%cls(1)%pp(2, 1)%gamma / o_g(2), 1.0_dp, 0.10_dp, 'gamma_d within 10% (tiny number)')
  print '(a,f9.5,a,f9.5,a,f9.5,a,f9.5)', '    E_F-V(S): ours ', res%ef - res%cls(1)%v(res%cls(1)%mesh%n, 1), '  oracle ', o_ef - o_vmtz, '   V(S): ours ', res%cls(1)%v(res%cls(1)%mesh%n, 1), '  oracle ', o_vmtz
  call check_close(res%ef - res%cls(1)%v(res%cls(1)%mesh%n, 1), o_ef - o_vmtz, 0.012_dp, 'E_F - V(S) within 12 mRy of Stuttgart')
  call check_close(res%cls(1)%pp(2, 1)%enu - res%ef, o_enu(2) - o_ef, 0.025_dp, 'E_nu,d - E_F within 25 mRy of Stuttgart')
  call check_close(res%cls(1)%v(res%cls(1)%mesh%n, 1), o_vmtz, 0.025_dp, 'V(S) within 25 mRy of Stuttgart')
  call report()
end program test_asa_ag
