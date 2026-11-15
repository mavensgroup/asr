!> probe_cuni_scf <nlevels> <lsuper> <rank> <flip_radius> [x_Cu=0.5] [nspin=1] [maxit=40]
!> Disordered fcc Cu_x Ni_(1-x) alloy SCF (a = 3.57 A) from the augmented-space recursion.
program probe_cuni_scf
  use asr_kinds, only: dp
  use asr_alloy_scf
  use asr_asa, only: core_shells_for
  use asr_atom, only: shell_t
  implicit none
  type(alloy_input_t) :: ai
  type(alloy_result_t) :: ar
  character(len=16) :: arg
  integer :: l, s
  call get_command_argument(1, arg); read (arg, *) ai%nlevels
  call get_command_argument(2, arg); read (arg, *) ai%lsuper
  call get_command_argument(3, arg); read (arg, *) ai%maxrank
  call get_command_argument(4, arg); read (arg, *) ai%flip_radius
  ai%x = 0.5_dp; ai%nspin = 1; ai%maxit = 40
  if (command_argument_count() >= 5) then; call get_command_argument(5, arg); read (arg, *) ai%x; end if
  if (command_argument_count() >= 6) then; call get_command_argument(6, arg); read (arg, *) ai%nspin; end if
  if (command_argument_count() >= 7) then; call get_command_argument(7, arg); read (arg, *) ai%maxit; end if
  ai%lattice = 'fcc'; ai%alat = 3.57_dp / 0.529177_dp; ai%tol = 1.0e-5_dp; ai%zeta0 = 0.1_dp
  if (command_argument_count() >= 8) then; call get_command_argument(8, arg); read (arg, *) ai%zeta0; end if
  if (command_argument_count() >= 9) then; call get_command_argument(9, arg); if (trim(arg) == 'ni') ai%alat = 3.524_dp / 0.529177_dp; end if
  ai%spec(1)%z = 29.0_dp; ai%spec(1)%nval = 11
  ai%spec(2)%z = 28.0_dp; ai%spec(2)%nval = 10
  do s = 1, 2
    ai%spec(s)%core = core_shells_for([4, 4, 3, 4]); ai%spec(s)%nvl = [4, 4, 3]
  end do
  ai%spec(1)%atom_shells = [ai%spec(1)%core, shell_t(4, 0, 1.0_dp), shell_t(3, 2, 10.0_dp)]
  if (command_argument_count() >= 9) then
    if (trim(arg) == 'ni') then       ! pure-Ni check: both species Ni
      ai%spec(1)%z = 28.0_dp; ai%spec(1)%nval = 10
      ai%spec(1)%atom_shells = [ai%spec(1)%core, shell_t(4, 0, 2.0_dp), shell_t(3, 2, 8.0_dp)]
    end if
  end if
  ai%spec(2)%atom_shells = [ai%spec(2)%core, shell_t(4, 0, 2.0_dp), shell_t(3, 2, 8.0_dp)]
  call alloy_scf(ai, ar, verbose=.true.)
  print '(a,l1,a,i0,a,i0,a,f10.6)', 'converged=', ar%converged, ' iterations=', ar%iter, '  nodes=', ar%nodes, '  E_F=', ar%ef
  print '(a,2f9.5,a,2f8.4)', 'charge transfer qnet (Cu, Ni) =', ar%cls(1)%qnet, ar%cls(2)%qnet, '   moments (Cu, Ni) =', ar%cls(1)%moment, ar%cls(2)%moment
  do s = 1, 2
    do l = 0, 2
      print '(a,i1,a,i1,a,5f11.5)', ' species ', s, ' l=', l, '  enu c sqrdel p gamma: ', ar%cls(s)%pp(l, 1)%enu, ar%cls(s)%pp(l, 1)%c, ar%cls(s)%pp(l, 1)%sqrdel, ar%cls(s)%pp(l, 1)%p, ar%cls(s)%pp(l, 1)%gamma
    end do
  end do
end program probe_cuni_scf
