!> probe_cuzn <nlevels> <lsuper> <rank> [x_Cu=0.5] [eta=0.02] [outfile] [sim_beta] [ed=0|1] [mix] [ncont] [a_Angstrom]
!> Disordered bcc Cu_x Zn_(1-x) (a = 2.85 A): alloy SCF from the augmented-space recursion, then the species- and
!> l-resolved DOS (states/Ry/atom, both spins) written to a file for comparison with the published curve.
program probe_cuzn
  use asr_kinds, only: dp
  use asr_alloy_scf
  use asr_asa, only: core_shells_for
  use asr_atom, only: shell_t
  implicit none
  type(alloy_input_t) :: ai
  type(alloy_result_t) :: ar
  character(len=256) :: arg, outfile
  integer, parameter :: ne = 901
  real(dp) :: egrid(ne), eta, tot, de
  real(dp), allocatable :: dos(:, :, :, :)
  integer :: s, ie, u, l
  real(dp) :: lat_a
  call get_command_argument(1, arg); read (arg, *) ai%nlevels
  call get_command_argument(2, arg); read (arg, *) ai%lsuper
  call get_command_argument(3, arg); read (arg, *) ai%maxrank
  ai%flip_radius = 99; ai%x = 0.5_dp; eta = 0.02_dp; outfile = 'cuzn_dos.dat'
  if (command_argument_count() >= 4) then; call get_command_argument(4, arg); read (arg, *) ai%x; end if
  if (command_argument_count() >= 5) then; call get_command_argument(5, arg); read (arg, *) eta; end if
  if (command_argument_count() >= 6) then; call get_command_argument(6, outfile); end if
  if (command_argument_count() >= 7) then; call get_command_argument(7, arg); read (arg, *) ai%sim_beta; end if
  if (command_argument_count() >= 8) then; call get_command_argument(8, arg); ai%energy_dependent = (trim(arg) == '1'); end if
  if (command_argument_count() >= 9) then; call get_command_argument(9, arg); read (arg, *) ai%beta; end if
  if (command_argument_count() >= 10) then; call get_command_argument(10, arg); read (arg, *) ai%ncont; end if
  ai%ndos = ne; ai%dos_emin = -0.9_dp; ai%dos_emax = 0.7_dp; ai%dos_eta = eta
  ai%lattice = 'bcc'; ai%alat = 2.85_dp / 0.529177_dp; ai%nspin = 1; ai%maxit = 50; ai%tol = 1.0e-5_dp
  if (command_argument_count() >= 11) then; call get_command_argument(11, arg); read (arg, *) lat_a; ai%alat = lat_a / 0.529177_dp; end if
  ai%spec(1)%z = 29.0_dp; ai%spec(1)%nval = 11
  ai%spec(2)%z = 30.0_dp; ai%spec(2)%nval = 12
  do s = 1, 2
    ai%spec(s)%core = core_shells_for([4, 4, 3, 4]); ai%spec(s)%nvl = [4, 4, 3]
  end do
  ai%spec(1)%atom_shells = [ai%spec(1)%core, shell_t(4, 0, 1.0_dp), shell_t(3, 2, 10.0_dp)]
  ai%spec(2)%atom_shells = [ai%spec(2)%core, shell_t(4, 0, 2.0_dp), shell_t(3, 2, 10.0_dp)]
  call alloy_scf(ai, ar, verbose=.true.)
  print '(a,l1,a,i0,a,i0,a,f10.6)', 'converged=', ar%converged, ' iterations=', ar%iter, '  nodes=', ar%nodes, '  E_F=', ar%ef
  print '(a,2f9.5)', 'charge transfer qnet (Cu, Zn) =', ar%cls(1)%qnet, ar%cls(2)%qnet
  do s = 1, 2
    do l = 0, 2
      print '(a,i1,a,i1,a,5f11.5)', ' species ', s, ' l=', l, '  enu c sqrdel p gamma: ', ar%cls(s)%pp(l, 1)%enu, ar%cls(s)%pp(l, 1)%c, ar%cls(s)%pp(l, 1)%sqrdel, ar%cls(s)%pp(l, 1)%p, ar%cls(s)%pp(l, 1)%gamma
    end do
  end do
  allocate (dos(ne, 0:2, 2, 1))
  de = 1.6_dp / (ne - 1)
  do ie = 1, ne; egrid(ie) = -0.9_dp + (ie - 1) * de; end do
  if (ar%ed) then
    dos = ar%edos
  else
    call alloy_ldos(ar, egrid, eta, dos)
  end if
  open (newunit=u, file=trim(outfile), status='replace')
  write (u, '(a,f8.5)') '# Cu_x Zn_(1-x) bcc a=2.85 A;  E_F = ', ar%ef
  write (u, '(a)') '# E(Ry)  total(states/Ry/atom, both spins)  Cu_s Cu_p Cu_d  Zn_s Zn_p Zn_d  (per spin x2)'
  do ie = 1, ne
    tot = 2.0_dp * (ai%x * sum(dos(ie, :, 1, 1)) + (1.0_dp - ai%x) * sum(dos(ie, :, 2, 1)))
    write (u, '(f10.5,8es13.5)') egrid(ie), tot, (2.0_dp * dos(ie, l, 1, 1), l=0, 2), (2.0_dp * dos(ie, l, 2, 1), l=0, 2)
  end do
  close (u)
  print '(a,a)', 'wrote ', trim(outfile)
end program probe_cuzn
