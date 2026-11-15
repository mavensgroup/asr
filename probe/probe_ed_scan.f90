!> probe_ed_scan: N(E_F) on a grid for ordered Ag (A = B) with the energy-dependent recursion at iteration 1 (diagnostic).
program probe_ed_scan
  use asr_kinds, only: dp
  use asr_alloy_scf
  use asr_asa, only: core_shells_for
  use asr_atom, only: shell_t
  implicit none
  type(alloy_input_t) :: ai
  type(alloy_result_t) :: ar
  character(len=16) :: arg
  integer :: s
  ai%lattice = 'fcc'; ai%alat = 7.72898293_dp; ai%x = 0.5_dp; ai%nlevels = 8; ai%lsuper = 12; ai%maxit = 40
  if (command_argument_count() >= 1) then; call get_command_argument(1, arg); read (arg, *) ai%nlevels; end if
  if (command_argument_count() >= 2) then; call get_command_argument(2, arg); read (arg, *) ai%ncont; end if
  ai%ed_eta = 0.01_dp
  if (command_argument_count() >= 3) then; call get_command_argument(3, arg); read (arg, *) ai%ed_eta; end if
  ai%energy_dependent = .true.; ai%ed_scan = .true.
  do s = 1, 2
    ai%spec(s)%z = 47.0_dp; ai%spec(s)%nval = 11; ai%spec(s)%core = core_shells_for([5, 5, 4, 4]); ai%spec(s)%nvl = [5, 5, 4]
    ai%spec(s)%atom_shells = [ai%spec(s)%core, shell_t(5, 0, 1.0_dp), shell_t(4, 2, 10.0_dp)]
  end do
  call alloy_scf(ai, ar, verbose=.true.)
end program probe_ed_scan
