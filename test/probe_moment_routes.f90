!> probe_moment_routes: spin moment of a pure ferromagnet (bcc Fe | fcc Ni) by three routes at the same LSDA (vBH):
!>   asa  ordered k-space ASA (first-order generalised problem, tetrahedra) -- the reference;
!>   h2k  ordered k-space, orthogonal second-order H = E_nu + h - h o h (same Hamiltonian as 'ei', no recursion);
!>   ei   recursion on the energy-independent second-order H = E_nu + h - h o h (alloy route, A = B, no disorder);
!>   ed   recursion on the energy-dependent P(z) - S, G = lambda + mu g mu (alloy route, A = B);
!>   edo  the same on the ordered real-space bcc graph (B2 with A = B, no augmented space; Fe only).
!> usage: probe_moment_routes fe|ni asa|h2k|ei|ed [levels] [supercell] [contour points] [nk] [avg2|bp] [expand|stuttgart] [rcl (alat)] [sigma] [rinv (alat)] [order 2|3]
program probe_moment_routes
  use asr_kinds, only: dp
  use asr_alloy_scf
  use asr_asa
  use asr_atom, only: shell_t
  use asr_lmto_cell
  implicit none
  type(alloy_input_t) :: ai
  type(alloy_result_t) :: ar
  type(asa_input_t) :: oi
  type(asa_result_t) :: orr
  type(cell_t) :: cell
  character(len=16) :: el, route, arg, term, dens
  real(dp) :: z, alat, q(0:2, 2), mom, ef, rcl, sig, rinv
  integer :: hord
  integer :: s, nval, nlev, lsup, ncont, nk, l, nd, it
  type(shell_t) :: val(2)
  logical :: conv

  if (command_argument_count() < 2) then             ! no arguments (e.g. inside asr_test): print usage, succeed
    print '(a)', 'probe_moment_routes: usage fe|ni asa|h2k|ei|ed|edo [levels] [supercell] [contour] [nk] [avg2|bp] [expand|stuttgart] [rcl] [sigma] [rinv]'
    stop
  end if
  call get_command_argument(1, el); call get_command_argument(2, route)
  nlev = 12; lsup = 24; ncont = 60; nk = 24
  if (command_argument_count() >= 3) then; call get_command_argument(3, arg); read (arg, *) nlev; end if
  if (command_argument_count() >= 4) then; call get_command_argument(4, arg); read (arg, *) lsup; end if
  if (command_argument_count() >= 5) then; call get_command_argument(5, arg); read (arg, *) ncont; end if
  if (command_argument_count() >= 6) then; call get_command_argument(6, arg); read (arg, *) nk; end if
  term = 'avg2'
  if (command_argument_count() >= 7) call get_command_argument(7, term)
  dens = 'expand'
  if (command_argument_count() >= 8) call get_command_argument(8, dens)
  rcl = 1.2_dp; sig = 0.0_dp            ! sigma <= 0: screening 'auto' (standard for ei, hard-core 0.7 for ed/edo)
  if (command_argument_count() >= 9) then; call get_command_argument(9, arg); read (arg, *) rcl; end if
  if (command_argument_count() >= 10) then; call get_command_argument(10, arg); read (arg, *) sig; end if
  rinv = 0.0_dp
  if (command_argument_count() >= 11) then; call get_command_argument(11, arg); read (arg, *) rinv; end if
  hord = 2
  if (command_argument_count() >= 12) then; call get_command_argument(12, arg); read (arg, *) hord; end if
  select case (trim(el))
  case ('fe'); z = 26.0_dp; nval = 8;  alat = 2.8665_dp; val = [shell_t(4, 0, 2.0_dp), shell_t(3, 2, 6.0_dp)]
  case ('ni'); z = 28.0_dp; nval = 10; alat = 3.524_dp;  val = [shell_t(4, 0, 2.0_dp), shell_t(3, 2, 8.0_dp)]
  case default; error stop 'element must be fe or ni'
  end select
  alat = alat / 0.529177210903_dp

  if (trim(route) == 'asa' .or. trim(route) == 'h2k') then
    oi%second_order = (trim(route) == 'h2k')
    allocate (oi%cls(1))
    oi%cls(1)%z = z; oi%cls(1)%nval = nval; oi%cls(1)%core = core_shells_for([4, 4, 3, 4]); oi%cls(1)%nvl = [4, 4, 3]
    oi%cls(1)%atom_shells = [oi%cls(1)%core, val]
    if (sig > 0.0_dp) oi%sigma = sig
    oi%alat = alat; oi%nk = nk; oi%maxit = 200; oi%nspin = 2; oi%zeta0 = 0.4_dp; oi%tol = 1.0e-6_dp
    cell%nat = 1
    if (trim(el) == 'fe') then
      cell%avec = reshape([-0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp], [3, 3])
    else
      cell%avec = reshape([0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp], [3, 3])
    end if
    allocate (cell%tau(3, 1), cell%isp(1)); cell%tau = 0.0_dp; cell%isp = 1
    call setup_cell(cell); oi%cls(1)%sws = cell%w * oi%alat
    call asa_scf(oi, cell, orr)
    q = orr%cls(1)%q(0, :, :); mom = orr%cls(1)%moment; ef = orr%ef; conv = orr%converged; it = orr%iter; nd = nk
  else
    ai%lattice = merge('bcc', 'fcc', trim(el) == 'fe'); ai%alat = alat; ai%x = 0.5_dp
    ai%nlevels = nlev; ai%lsuper = lsup; ai%ncont = ncont; ai%maxit = 150; ai%tol = 1.0e-6_dp
    ai%nspin = 2; ai%zeta0 = 0.1_dp; ai%maxrank = 0
    ai%energy_dependent = (trim(route) == 'ed' .or. trim(route) == 'edo'); ai%terminator = trim(term); ai%density = trim(dens)
    ai%ordered = (trim(route) == 'edo'); ai%rcl = rcl; ai%rinv = rinv; ai%order = hord
    if (sig > 0.0_dp) then; ai%screening = 'hardcore'; ai%sigma = sig; else; ai%screening = 'auto'; end if
    do s = 1, 2
      ai%spec(s)%z = z; ai%spec(s)%nval = nval; ai%spec(s)%core = core_shells_for([4, 4, 3, 4]); ai%spec(s)%nvl = [4, 4, 3]
      ai%spec(s)%atom_shells = [ai%spec(s)%core, val]
    end do
    call alloy_scf(ai, ar, verbose=.true.)
    q = ar%cls(1)%q(0, :, :); mom = ar%cls(1)%moment; ef = ar%ef; conv = ar%converged; it = ar%iter; nd = ar%nodes
  end if
  print '(a,a3,1x,a3,1x,a4,a,i3,a,l1,a,i4,a,f9.5,a,f8.4)', 'RESULT ', el, route, term, ' '//dens(1:4)//' levels', nlev, ' conv ', conv, ' it ', it, &
      '  E_F ', ef, '  moment ', mom
  print '(a,3f6.2,a,i0)', 'RESULT   rcl sigma rinv ', rcl, sig, rinv, '   order ', hord
  print '(a,3f8.4,a,3f8.4)', 'RESULT   occ up  s p d ', (q(l, 1), l=0, 2), '   down ', (q(l, 2), l=0, 2)
end program probe_moment_routes
