!> probe_fe_alloy <nlevels> <lsuper> [maxit=60] : trivial alloy A=B=bcc Fe, spin-polarised, verbose history.
program probe_fe_alloy
  use asr_kinds, only: dp
  use asr_alloy_scf
  use asr_asa, only: core_shells_for
  use asr_atom, only: shell_t
  implicit none
  type(alloy_input_t) :: ai
  type(alloy_result_t) :: ar
  character(len=16) :: arg
  integer :: s
  call get_command_argument(1, arg); read (arg, *) ai%nlevels
  call get_command_argument(2, arg); read (arg, *) ai%lsuper
  ai%maxit = 60; ai%zeta0 = 0.4_dp
  if (command_argument_count() >= 3) then; call get_command_argument(3, arg); read (arg, *) ai%maxit; end if
  if (command_argument_count() >= 4) then; call get_command_argument(4, arg); read (arg, *) ai%nspin; else; ai%nspin = 2; end if
  if (command_argument_count() >= 5) then; call get_command_argument(5, arg); read (arg, *) ai%beta; end if
  if (command_argument_count() >= 6) then; call get_command_argument(6, arg); read (arg, *) ai%ncont; end if
  if (command_argument_count() >= 7) then; call get_command_argument(7, arg); read (arg, *) ai%zeta0; end if
  ai%lattice = 'bcc'; ai%alat = 2.8665_dp / 0.529177_dp; ai%x = 0.5_dp; ai%tol = 1.0e-5_dp
  do s = 1, 2
    ai%spec(s)%z = 26.0_dp; ai%spec(s)%nval = 8; ai%spec(s)%core = core_shells_for([4, 4, 3, 4]); ai%spec(s)%nvl = [4, 4, 3]
    ai%spec(s)%atom_shells = [ai%spec(s)%core, shell_t(4, 0, 2.0_dp), shell_t(3, 2, 6.0_dp)]
  end do
  call alloy_scf(ai, ar, verbose=.true.)
  print '(a,l1,a,i0,a,f8.4,a,f10.6)', 'converged=', ar%converged, ' iterations=', ar%iter, '  |moment|=', abs(ar%cls(1)%moment), '  E_F=', ar%ef
end program probe_fe_alloy
