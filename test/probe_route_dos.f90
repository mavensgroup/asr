!> probe_route_dos <nlevels> <lsuper> <rank> <outfile> [x_Cu=0.5] [eta=0.02] [nlev_ed]
!> Fixed-potential comparison of the two recursion routes for bcc Cu_x Zn_(1-x): the potentials come from the converged
!> energy-independent alloy SCF; the species/l DOS is then computed (a) from the energy-independent recursion (H = E_nu + h - h o h)
!> and (b) from the energy-dependent recursion (P(z) - S) with the same parameters, graph and rank (depth nlev_ed, default = nlevels).
program probe_route_dos
  use asr_kinds, only: dp, pi
  use asr_alloy_scf
  use asr_asa, only: core_shells_for
  use asr_atom, only: shell_t
  use asr_lmto_ham, only: potpar_t
  use asr_lmto_cell, only: alpha_from_sigma
  use asr_lmto_alloy, only: build_lmto_alloy
  use asr_augspace, only: graph_t
  use asr_augorb, only: orb_alloy_t
  use asr_augsym
  use asr_endep
  use asr_symmetry, only: d4h_ops, orbital_rep
  implicit none
  integer, parameter :: ne = 801, iorb(4) = [1, 3, 5, 7], lorb(4) = [0, 1, 2, 2], lorb9(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]
  real(dp), parameter :: worb(4) = [1.0_dp, 3.0_dp, 3.0_dp, 2.0_dp]
  type(alloy_input_t) :: ai
  type(alloy_result_t) :: ar
  type(potpar_t) :: pp(2)
  type(graph_t) :: g
  type(orb_alloy_t) :: mdl
  type(augsym_op) :: sop
  type(ed_params_t) :: edp(2)
  type(ed_work_t) :: wk
  integer, allocatable :: coord(:, :), perm(:, :), mats(:, :, :)
  real(dp), allocatable :: dmat(:, :, :), dei(:, :, :, :), ded(:, :, :)
  character(len=256) :: arg, outfile
  integer :: s, l, ie, o, u, nlev_ed, i
  real(dp) :: eta, x, egrid(ne), alpha(0:2), alpha9(9), de, te(ne), td(ne), nei, ned, efd
  complex(dp) :: ga, gb

  call get_command_argument(1, arg); read (arg, *) ai%nlevels
  call get_command_argument(2, arg); read (arg, *) ai%lsuper
  call get_command_argument(3, arg); read (arg, *) ai%maxrank
  call get_command_argument(4, outfile)
  x = 0.5_dp; eta = 0.02_dp; nlev_ed = ai%nlevels
  if (command_argument_count() >= 5) then; call get_command_argument(5, arg); read (arg, *) x; end if
  if (command_argument_count() >= 6) then; call get_command_argument(6, arg); read (arg, *) eta; end if
  if (command_argument_count() >= 7) then; call get_command_argument(7, arg); read (arg, *) nlev_ed; end if
  ai%x = x; ai%flip_radius = 99; ai%lattice = 'bcc'; ai%alat = 2.85_dp / 0.529177_dp; ai%nspin = 1; ai%maxit = 60; ai%tol = 1.0e-5_dp
  ai%symmetry = .true.; ai%energy_dependent = .false.
  ai%spec(1)%z = 29.0_dp; ai%spec(1)%nval = 11
  ai%spec(2)%z = 30.0_dp; ai%spec(2)%nval = 12
  do s = 1, 2
    ai%spec(s)%core = core_shells_for([4, 4, 3, 4]); ai%spec(s)%nvl = [4, 4, 3]
  end do
  ai%spec(1)%atom_shells = [ai%spec(1)%core, shell_t(4, 0, 1.0_dp), shell_t(3, 2, 10.0_dp)]
  ai%spec(2)%atom_shells = [ai%spec(2)%core, shell_t(4, 0, 2.0_dp), shell_t(3, 2, 10.0_dp)]
  call alloy_scf(ai, ar, verbose=.false.)
  print '(a,l1,a,i0,a,f10.6)', 'energy-independent SCF converged=', ar%converged, ' its=', ar%iter, '  E_F=', ar%ef
  de = 1.6_dp / (ne - 1)
  do ie = 1, ne; egrid(ie) = -0.9_dp + (ie - 1) * de; end do
  allocate (dei(ne, 0:2, 2, 1), ded(ne, 0:2, 2))
  call alloy_ldos(ar, egrid, eta, dei)

  ! ---- same potentials in the energy-dependent route
  do s = 1, 2
    do l = 0, 2
      pp(s)%enu(l) = ar%cls(s)%pp(l, 1)%enu; pp(s)%c(l) = ar%cls(s)%pp(l, 1)%c; pp(s)%sqrdel(l) = ar%cls(s)%pp(l, 1)%sqrdel
      pp(s)%p(l) = ar%cls(s)%pp(l, 1)%p; pp(s)%gamma(l) = ar%cls(s)%pp(l, 1)%gamma; pp(s)%alpha_rep(l) = ar%cls(s)%pp(l, 1)%gamma
    end do
    do i = 1, 9
      edp(s)%c(i) = pp(s)%c(lorb9(i)); edp(s)%del(i) = pp(s)%sqrdel(lorb9(i))**2; edp(s)%gam(i) = pp(s)%gamma(lorb9(i))
    end do
  end do
  alpha = alloy_alpha(ai)
  do i = 1, 9; alpha9(i) = alpha(lorb9(i)); end do
  call build_lmto_alloy('bcc', ai%lsuper, ai%rcl, ai%sigma, pp(1), pp(2), x, g, coord, mdl, ai%rinv, alpha)
  call d4h_ops(coord, 2 * ai%lsuper, perm, mats)
  call orbital_rep(mats, dmat)
  call build_augsym_operator(g, mdl, perm, dmat, 1, nlev_ed, sop, ai%maxrank)
  call ed_work_init(wk, sop%dim())
  print '(a,i0,a,i0)', 'energy-dependent operator nodes=', sop%nn, '  levels=', nlev_ed
  ded = 0.0_dp
  do o = 1, 4
    call augsym_set_start(sop, iorb(o), dmat)
    do ie = 1, ne
      call ed_species_green(sop, wk, cmplx(egrid(ie), eta, dp), edp(1), edp(2), alpha9, x, iorb(o), nlev_ed, sop%flipr(1), sop%md(1), sop%s, ga, gb)
      ded(ie, lorb(o), 1) = ded(ie, lorb(o), 1) - worb(o) * aimag(ga) / pi
      ded(ie, lorb(o), 2) = ded(ie, lorb(o), 2) - worb(o) * aimag(gb) / pi
    end do
  end do
  do ie = 1, ne
    te(ie) = 2.0_dp * (x * sum(dei(ie, :, 1, 1)) + (1.0_dp - x) * sum(dei(ie, :, 2, 1)))
    td(ie) = 2.0_dp * (x * sum(ded(ie, :, 1)) + (1.0_dp - x) * sum(ded(ie, :, 2)))
  end do
  ! electron counts up to the energy-independent E_F, and each route's own E_F for 11.5 electrons
  nei = 0.0_dp; ned = 0.0_dp
  do ie = 1, ne
    if (egrid(ie) <= ar%ef) then; nei = nei + te(ie) * de; ned = ned + td(ie) * de; end if
  end do
  nei = 0.0_dp; ned = 0.0_dp; efd = 0.0_dp
  do ie = 1, ne
    nei = nei + te(ie) * de; ned = ned + td(ie) * de
  end do
  print '(a,f9.4,a,f9.4,a)', 'integral of total DOS over [-0.9, 0.7]: EI ', nei, '   ED ', ned, '  (electrons; the d-band and part of sp are inside)'
  nei = 0.0_dp; ned = 0.0_dp
  do ie = 1, ne
    if (egrid(ie) <= ar%ef) then; nei = nei + te(ie) * de; ned = ned + td(ie) * de; end if
  end do
  print '(a,f9.4,a,f9.4,a,f9.5)', 'electrons below the EI E_F: EI ', nei, '   ED ', ned, '   (E_F = ', ar%ef
  call peaks('EI', egrid, dei(:, 2, 1, 1), dei(:, 2, 2, 1), ar%ef)
  call peaks('ED', egrid, ded(:, 2, 1), ded(:, 2, 2), ar%ef)
  open (newunit=u, file=trim(outfile), status='replace')
  write (u, '(a,f10.6)') '# fixed-potential route comparison, bcc CuZn;  E_F(EI) = ', ar%ef
  write (u, '(a)') '# E(Ry) tot_EI tot_ED  Cu_d_EI Cu_d_ED  Zn_d_EI Zn_d_ED  Cu_sp_EI Cu_sp_ED  Zn_sp_EI Zn_sp_ED   (states/Ry/atom, both spins)'
  do ie = 1, ne
    write (u, '(f10.5,11es13.5)') egrid(ie), te(ie), td(ie), 2.0_dp * dei(ie, 2, 1, 1), 2.0_dp * ded(ie, 2, 1), 2.0_dp * dei(ie, 2, 2, 1), 2.0_dp * ded(ie, 2, 2), &
      2.0_dp * (dei(ie, 0, 1, 1) + dei(ie, 1, 1, 1)), 2.0_dp * (ded(ie, 0, 1) + ded(ie, 1, 1)), 2.0_dp * (dei(ie, 0, 2, 1) + dei(ie, 1, 2, 1)), 2.0_dp * (ded(ie, 0, 2) + ded(ie, 1, 2))
  end do
  close (u)
  print '(a,a)', 'wrote ', trim(outfile)
contains
  subroutine peaks(tag, e, dcu, dzn, ef)
    character(len=*), intent(in) :: tag
    real(dp), intent(in) :: e(:), dcu(:), dzn(:), ef
    print '(a,a,a,f8.4,a,f8.4,a,f8.4)', tag, ' route:  Cu d peak at ', '', e(maxloc(dcu, 1)) - ef, '   Zn d peak at ', e(maxloc(dzn, 1)) - ef, '  (relative to the EI E_F ', ef
  end subroutine peaks
end program probe_route_dos
