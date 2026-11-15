!> probe_alloy_scf [nlevels=8] [lsuper=34] : A = B = Ag (no disorder) must reproduce the ordered SCF.
program probe_alloy_scf
  use asr_kinds, only: dp
  use asr_alloy_scf
  use asr_asa, only: core_shells_for
  use asr_atom, only: shell_t
  implicit none
  type(alloy_input_t) :: inp
  type(alloy_result_t) :: res
  character(len=16) :: arg
  integer :: s, l
  inp%nlevels = 8; inp%lsuper = 34
  if (command_argument_count() >= 1) then; call get_command_argument(1, arg); read (arg, *) inp%nlevels; end if
  if (command_argument_count() >= 2) then; call get_command_argument(2, arg); read (arg, *) inp%lsuper; end if
  inp%lattice = 'fcc'; inp%alat = 7.72898293_dp; inp%x = 0.5_dp
  do s = 1, 2
    inp%spec(s)%z = 47.0_dp; inp%spec(s)%nval = 11
    inp%spec(s)%core = core_shells_for([5, 5, 4, 4]); inp%spec(s)%nvl = [5, 5, 4]
    inp%spec(s)%atom_shells = [inp%spec(s)%core, shell_t(5, 0, 1.0_dp), shell_t(4, 2, 10.0_dp)]
  end do
  inp%maxit = 40; inp%tol = 1.0e-5_dp; inp%ncont = 60
  call alloy_scf(inp, res, verbose=.true.)
  print '(a,l1,a,i0,a,es9.2,a,f10.6)', 'converged=', res%converged, ' it=', res%iter, ' dV=', res%dv, '  E_F=', res%ef
  do l = 0, 2
    print '(a,i1,a,5f11.5)', ' l=', l, '  enu c sqrdel p gamma: ', res%cls(1)%pp(l, 1)%enu, res%cls(1)%pp(l, 1)%c, res%cls(1)%pp(l, 1)%sqrdel, &
      res%cls(1)%pp(l, 1)%p, res%cls(1)%pp(l, 1)%gamma
  end do
  print '(a)', ' ordered first-order ASA (asa_scf):  E_F -0.1737 (tetra),  Ag l=0: enu -0.5340 c -0.4494 sqrdel 0.3649 gamma 0.4249;  l=2: enu -0.4859 c -0.4889 sqrdel -0.1044 gamma 0.0052'
end program probe_alloy_scf
