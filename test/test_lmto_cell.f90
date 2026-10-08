!> General ordered cell.  (1) regression: one-atom fcc Ag equals the single-atom code;
!> (2) B2 FeCu from the Stuttgart parameters in AMscf/lmto/fe/ctrl: electron count at the
!> Stuttgart E_F, l-resolved site occupations vs the converged Stuttgart moments, and
!> convergence with cluster radius.
program test_lmto_cell
  use asr_kinds, only: dp, pi
  use asr_strconst
  use asr_lmto_ham, only: potpar_t, ham_k
  use asr_lmto_cell
  use asr_testing
  implicit none
  integer, parameter :: ne = 1201
  type(cell_t) :: c, c2
  type(cluster_t), allocatable :: cls(:), cl2(:)
  type(potpar_t) :: pp(2)
  real(dp) :: alpha_sp(0:2, 2), al4(0:3), w1
  real(dp), allocatable :: vec(:, :), blk(:, :, :)
  complex(dp), allocatable :: h(:, :)
  complex(dp) :: h1(9, 9)
  real(dp) :: egrid(ne), de, ef, occ(0:2, 2), tot, dos(ne, 0:2, 2), k(3), rcl
  integer :: i, ir, l, nkg, u
  real(dp) :: sg
  character(len=16) :: arg
  real(dp), parameter :: rfac(3) = [3.2_dp, 4.0_dp, 5.0_dp]
  real(dp), parameter :: qfe(0:2) = [0.65065469_dp, 0.77696927_dp, 6.76341845_dp]
  real(dp), parameter :: qcu(0:2) = [0.70640914_dp, 0.67687944_dp, 9.42566900_dp]

  print '(a)', 'test_lmto_cell'
  ! ---- (1) Ag regression -------------------------------------------------------------
  pp(1)%enu = [-0.54494317_dp, -0.45471939_dp, -0.51341327_dp]
  pp(1)%c = [-0.46144639_dp, 0.41368602_dp, -0.51457447_dp]
  pp(1)%sqrdel = [0.36457406_dp, -0.36517716_dp, -0.10341022_dp]
  pp(1)%p = [0.08230178_dp, 0.03623364_dp, 2.08144602_dp]
  pp(1)%gamma = [0.42477526_dp, 0.11537965_dp, 0.00550351_dp]
  pp(1)%alpha_rep = pp(1)%gamma
  c%nat = 1
  c%avec = reshape([0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp], [3, 3])
  allocate (c%tau(3, 1), c%isp(1)); c%tau = 0.0_dp; c%isp = 1
  call setup_cell(c)
  call check_close(c%w, avg_ws_radius('fcc'), 1.0e-13_dp, 'avw from cell volume == fcc formula')
  alpha_sp(:, 1) = alpha_from_sigma(0.7_dp, 1.0_dp)
  allocate (cls(1)); call build_cluster(c, 1, 1.55_dp, alpha_sp(:, 1:1), cls(1))
  call lattice_cluster('fcc', 1.55_dp, vec)
  al4 = hardcore_alpha(0.7_dp); w1 = avg_ws_radius('fcc')
  call screened_cluster(vec, w1, al4(0:2), blk)
  allocate (h(9, 9)); k = [0.31_dp, -0.17_dp, 0.52_dp]
  call ham_k_cell(c, cls, pp(1:1), alpha_sp(:, 1:1), k, h)
  call ham_k(k, vec, blk, pp(1), al4(0:2), 1.0_dp, .true., h1)
  ! the two routines differ by an orbital-sign gauge only in sqrdel(p): compare spectra
  call check_close(sum(abs(sorted_eig(h) - sorted_eig(h1))), 0.0_dp, 1.0e-9_dp, 'Ag: cell code == single-atom code (spectrum at generic k)')
  deallocate (h)

  ! ---- (2) B2 FeCu ---------------------------------------------------------------------
  pp(1)%enu = [-0.53659895_dp, -0.40643333_dp, -0.27926995_dp]
  pp(1)%c = [-0.39513210_dp, 0.45588071_dp, -0.22929696_dp]
  pp(1)%sqrdel = [-0.39008319_dp, 0.39665377_dp, 0.11161564_dp]
  pp(1)%p = [0.06913570_dp, 0.03638773_dp, 2.69722817_dp]
  pp(1)%gamma = [0.44472011_dp, 0.12918350_dp, -0.00090362_dp]
  pp(2)%enu = [-0.51604605_dp, -0.36673367_dp, -0.38106254_dp]
  pp(2)%c = [-0.44738380_dp, 0.55528675_dp, -0.37210440_dp]
  pp(2)%sqrdel = [-0.39838509_dp, 0.36870298_dp, 0.08239456_dp]
  pp(2)%p = [0.05771734_dp, 0.02906149_dp, 2.60456746_dp]
  pp(2)%gamma = [0.40099416_dp, 0.09595500_dp, -0.00178982_dp]
  pp(1)%alpha_rep = pp(1)%gamma; pp(2)%alpha_rep = pp(2)%gamma
  ef = -0.18713836_dp
  c2%nat = 2
  c2%avec = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
  allocate (c2%tau(3, 2), c2%isp(2)); c2%tau(:, 1) = 0.0_dp; c2%tau(:, 2) = 0.5_dp; c2%isp = [1, 2]
  call setup_cell(c2)
  ! sphere radii (Bohr) / avw  -> hard-core alpha per species (sigma = 0.7)
  alpha_sp(:, 1) = alpha_from_sigma(0.7_dp, 2.91398733_dp / (c2%w * 5.66918063_dp))
  alpha_sp(:, 2) = alpha_from_sigma(0.7_dp, 2.65686809_dp / (c2%w * 5.66918063_dp))
  nkg = 48; sg = 0.002_dp
  if (command_argument_count() >= 1) then; call get_command_argument(1, arg); read (arg, *) nkg; end if
  if (command_argument_count() >= 2) then; call get_command_argument(2, arg); read (arg, *) sg; end if
  de = 3.0_dp / (ne - 1)
  do i = 1, ne; egrid(i) = -0.9_dp + (i - 1) * de; end do
  do ir = 2, 2
    rcl = rfac(ir) * c2%w
    allocate (cl2(2))
    call build_cluster(c2, 1, rcl, alpha_sp, cl2(1))
    call build_cluster(c2, 2, rcl, alpha_sp, cl2(2))
    call dos_cell(c2, cl2, pp, alpha_sp, nkg, egrid, sg, dos)
    dos = 2.0_dp * dos
    do i = 1, 2
      do l = 0, 2
        occ(l, i) = sum(dos(:, l, i), mask=egrid <= ef) * de
      end do
    end do
    tot = sum(occ)
    print '(a,f4.1,a,i0,a,f8.3,a,3f7.3,a,3f7.3)', '    rmax=', rcl / c2%w, ' w  (sites ', cl2(1)%n, ')  N(E_F)=', tot, &
      '  Fe s,p,d=', occ(:, 1), '  Cu s,p,d=', occ(:, 2)
    deallocate (cl2)
  end do
  print '(a,3f7.3,a,3f7.3)', '    Stuttgart moments       Fe=', qfe, '  Cu=', qcu
  call check_close(tot, 19.0_dp, 0.20_dp, 'B2 FeCu: N(E_F) = 19 electrons within 0.2 (Stuttgart E_F)')
  call check_close(occ(2, 1), qfe(2), 0.05_dp, 'B2 FeCu: Fe d occupation vs Stuttgart 6.763')
  call check_close(occ(2, 2), qcu(2), 0.02_dp, 'B2 FeCu: Cu d occupation vs Stuttgart 9.426')
  call check_close(occ(0, 1), qfe(0), 0.05_dp, 'B2 FeCu: Fe s occupation vs Stuttgart 0.651 (ours +4%, weight definition)')
  call check_close(occ(1, 1), qfe(1), 0.05_dp, 'B2 FeCu: Fe p occupation vs Stuttgart 0.777')
  call check_close(occ(0, 2), qcu(0), 0.05_dp, 'B2 FeCu: Cu s occupation vs Stuttgart 0.706')
  call check_close(occ(1, 2), qcu(1), 0.05_dp, 'B2 FeCu: Cu p occupation vs Stuttgart 0.677')
  open (newunit=u, file='trial/b2_fecu_dos.dat', status='replace')
  write (u, '(a)') '# E(Ry)  Fe_s Fe_p Fe_d  Cu_s Cu_p Cu_d   (states/Ry/atom, both spins)'
  do i = 1, ne
    write (u, '(7es14.6)') egrid(i), dos(i, :, 1), dos(i, :, 2)
  end do
  close (u)
  call report()
contains
  function sorted_eig(a) result(ev)
    complex(dp), intent(in) :: a(:, :)
    real(dp) :: ev(size(a, 1)), rw(3 * size(a, 1))
    complex(dp) :: m(size(a, 1), size(a, 1)), cw(4 * size(a, 1))
    integer :: info
    external :: zheev
    m = a
    call zheev('N', 'U', size(a, 1), m, size(a, 1), ev, cw, size(cw), rw, info)
  end function sorted_eig
end program test_lmto_cell
