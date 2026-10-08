!> Spin-polarised alloy SCF, trivial alloy A = B = bcc Fe (no disorder): the ferromagnetic moment from the
!> augmented-space recursion route must agree with the ordered ASA SCF (asa_scf).
program test_alloy_scf_spin
  use asr_kinds, only: dp
  use asr_alloy_scf
  use asr_asa
  use asr_atom, only: shell_t
  use asr_lmto_cell
  use asr_testing
  implicit none
  type(alloy_input_t) :: ai
  type(alloy_result_t) :: ar
  type(asa_input_t) :: oi
  type(asa_result_t) :: orr
  type(cell_t) :: cell
  integer :: s
  print '(a)', 'test_alloy_scf_spin'
  allocate (oi%cls(1))
  oi%cls(1)%z = 26.0_dp; oi%cls(1)%nval = 8; oi%cls(1)%core = core_shells_for([4, 4, 3, 4]); oi%cls(1)%nvl = [4, 4, 3]
  oi%cls(1)%atom_shells = [oi%cls(1)%core, shell_t(4, 0, 2.0_dp), shell_t(3, 2, 6.0_dp)]
  oi%alat = 2.8665_dp / 0.529177_dp; oi%nk = 14; oi%maxit = 150; oi%nspin = 2; oi%zeta0 = 0.4_dp; oi%tol = 1.0e-5_dp
  cell%nat = 1; cell%avec = reshape([-0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp], [3, 3])
  allocate (cell%tau(3, 1), cell%isp(1)); cell%tau = 0.0_dp; cell%isp = 1
  call setup_cell(cell); oi%cls(1)%sws = cell%w * oi%alat
  call asa_scf(oi, cell, orr)
  ai%lattice = 'bcc'; ai%alat = oi%alat; ai%x = 0.5_dp; ai%nlevels = 10; ai%lsuper = 36; ai%maxit = 70; ai%tol = 1.0e-5_dp
  ai%nspin = 2; ai%zeta0 = 0.1_dp
  do s = 1, 2
    ai%spec(s)%z = 26.0_dp; ai%spec(s)%nval = 8; ai%spec(s)%core = core_shells_for([4, 4, 3, 4]); ai%spec(s)%nvl = [4, 4, 3]
    ai%spec(s)%atom_shells = [ai%spec(s)%core, shell_t(4, 0, 2.0_dp), shell_t(3, 2, 6.0_dp)]
  end do
  call alloy_scf(ai, ar)
  print '(a,f8.4,a,f8.4,a,i0,a,i0)', '    |moment| ordered ASA = ', abs(orr%cls(1)%moment), '   alloy route (A=B) = ', abs(ar%cls(1)%moment), '   it=', ar%iter, '  nodes=', ar%nodes
  call check(orr%converged, 'ordered Fe SCF converged')
  call check(ar%converged, 'alloy-route Fe SCF converged')
  call check(abs(ar%cls(1)%moment) > 1.9_dp .and. abs(ar%cls(1)%moment) < 2.5_dp, 'alloy-route Fe moment is ferromagnetic ~2.2 muB')
  call check_close(abs(ar%cls(1)%moment), abs(orr%cls(1)%moment), 0.2_dp, 'moment within 0.2 muB of the ordered ASA route')
  call check_close(ar%cls(1)%qnet, 0.0_dp, 1.0e-3_dp, 'sphere neutral')
  call report()
end program test_alloy_scf_spin
