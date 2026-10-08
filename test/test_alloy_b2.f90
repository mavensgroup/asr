!> test_alloy_b2: ordered B2 CuZn (a = 2.92055 A).  (1) real-space ordered SCF (alloy_scf, ordered);
!> (2) at the SAME converged parameters, the second-order Hamiltonian H = E + h - h o h in k-space (ham_k_cell, nk=24):
!> Fermi level from the electron count and the Cu d / Zn d centroids, to compare with the recursion.  Separates
!> recursion (depth / terminator) error from the H^(2) vs first-order-ASA representation difference.
!> (3) the third-order recursion Hamiltonian (order = 3) against the self-consistent first-order k-space ASA (asa_scf):
!> H^(2) misses the overlap term h p h and puts E_F ~35 mRy too high; H^(3) must close that gap.
program test_alloy_b2
  use asr_kinds, only: dp
  use asr_alloy_scf
  use asr_asa, only: core_shells_for, asa_input_t, asa_result_t, asa_scf
  use asr_lmto_driver, only: atom_shells
  use asr_input_lmto, only: lmto_species_t
  use asr_lmto_cell
  use asr_lmto_ham, only: potpar_t
  use asr_testing
  implicit none
  type(alloy_input_t) :: inp
  type(alloy_result_t) :: res, res3
  type(asa_input_t) :: ki
  type(asa_result_t) :: kr
  type(cell_t) :: kcell
  type(lmto_species_t) :: sp(2)
  type(cell_t) :: cell
  type(cluster_t), allocatable :: cls(:)
  type(potpar_t) :: pps(2)
  character(len=16) :: arg
  real(dp), allocatable :: eb(:), w(:, :), alpha_sp(:, :)
  complex(dp), allocatable :: h(:, :), cw(:)
  real(dp) :: bv(3, 3), k(3), avw, lo, hi, ef, nn, alat, rwork(100)
  integer :: nk, i, j, l, ix, iy, iz, ik, nd, info, nst, it, lw
  external :: zheev

  inp%nlevels = 16; inp%lsuper = 12; alat = 2.92055_dp
  if (command_argument_count() >= 1) then; call get_command_argument(1, arg); read (arg, *) inp%nlevels; end if
  if (command_argument_count() >= 2) then; call get_command_argument(2, arg); read (arg, *) inp%lsuper; end if
  if (command_argument_count() >= 3) then; call get_command_argument(3, arg); read (arg, *) alat; end if
  sp(1)%symbol = 'Cu'; sp(1)%z = 29; sp(1)%nval = 11; sp(1)%konfig = [4, 4, 3, 4]; sp(1)%occ_s = 1.0_dp; sp(1)%occ_d = 10.0_dp
  sp(2)%symbol = 'Zn'; sp(2)%z = 30; sp(2)%nval = 12; sp(2)%konfig = [4, 4, 3, 4]; sp(2)%occ_s = 2.0_dp; sp(2)%occ_d = 10.0_dp
  inp%lattice = 'bcc'; inp%ordered = .true.; inp%alat = alat / 0.529177210903_dp; inp%x = 0.5_dp
  inp%maxit = 60; inp%tol = 1.0e-6_dp; inp%ncont = 80; inp%beta = 0.3_dp
  do i = 1, 2
    inp%spec(i)%z = sp(i)%z; inp%spec(i)%nval = sp(i)%nval; inp%spec(i)%nvl = sp(i)%konfig(0:2)
    inp%spec(i)%core = core_shells_for(sp(i)%konfig); inp%spec(i)%atom_shells = atom_shells(sp(i))
  end do
  call alloy_scf(inp, res, verbose=.false.)
  print '(a,l1,a,i0,a,f10.6,a,2f9.5)', ' recursion: converged=', res%converged, ' it=', res%iter, '  E_F=', res%ef, '  qnet(Cu,Zn)=', res%cls(1)%qnet, res%cls(2)%qnet

  ! k-space H^(2) at the same parameters
  cell%nat = 2
  cell%avec = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
  allocate (cell%tau(3, 2), cell%isp(2)); cell%tau = 0.0_dp; cell%tau(:, 2) = 0.5_dp; cell%isp = [1, 2]
  call setup_cell(cell)
  avw = cell%w * inp%alat
  allocate (alpha_sp(0:2, 2), cls(2))
  do i = 1, 2
    alpha_sp(:, i) = alloy_alpha(inp)
    do l = 0, 2
      pps(i)%enu(l) = res%cls(i)%pp(l, 1)%enu; pps(i)%c(l) = res%cls(i)%pp(l, 1)%c; pps(i)%sqrdel(l) = res%cls(i)%pp(l, 1)%sqrdel
      pps(i)%p(l) = res%cls(i)%pp(l, 1)%p; pps(i)%gamma(l) = res%cls(i)%pp(l, 1)%gamma; pps(i)%alpha_rep(l) = res%cls(i)%pp(l, 1)%gamma
    end do
  end do
  call build_cluster(cell, 1, 4.0_dp * cell%w, alpha_sp, cls(1)); call build_cluster(cell, 2, 4.0_dp * cell%w, alpha_sp, cls(2))
  nk = 24; nd = 18; nst = nk**3 * nd; bv = recip_vectors(cell%avec); lw = 4 * nd
  allocate (eb(nst), w(nst, 2), h(nd, nd), cw(lw))
  !$omp parallel do collapse(2) private(h, cw, k, ik, info, rwork, ix) schedule(static)
  do iz = 0, nk - 1
    do iy = 0, nk - 1
      do ix = 0, nk - 1
        k = (ix + 0.5_dp) / nk * bv(:, 1) + (iy + 0.5_dp) / nk * bv(:, 2) + (iz + 0.5_dp) / nk * bv(:, 3)
        ik = ix + nk * (iy + nk * iz)
        call ham_k_cell(cell, cls, pps, alpha_sp, k, h)
        h = 0.5_dp * (h + conjg(transpose(h)))
        call zheev('V', 'U', nd, h, nd, eb(nd * ik + 1:nd * ik + nd), cw, lw, rwork, info)
        do j = 1, nd
          w(nd * ik + j, 1) = sum(abs(h(1:9, j))**2); w(nd * ik + j, 2) = sum(abs(h(10:18, j))**2)
        end do
      end do
    end do
  end do
  !$omp end parallel do
  lo = minval(eb) - 1.0_dp; hi = maxval(eb) + 1.0_dp
  do it = 1, 200
    ef = 0.5_dp * (lo + hi); nn = 0.0_dp
    do i = 1, nst; if (eb(i) < ef) nn = nn + 2.0_dp / real(nk**3, dp); end do
    if (nn > 23.0_dp) then; hi = ef; else; lo = ef; end if
  end do
  print '(a,f10.6,a)', 'k-space H^(2), nk=24, 23 e/cell:  E_F=', ef, '   (step-function count; resolution ~ 1/nk^3 per state)'
  print '(a,f10.6)', 'difference recursion - k-space(H2) E_F (Ry): ', res%ef - ef
  call check(res%converged, 'ordered B2 real-space SCF converged')
  call check(abs(res%cls(1)%qnet + res%cls(2)%qnet) < 1.0e-4_dp, 'cell neutral')
  call check(res%cls(1)%qnet > 0.12_dp .and. res%cls(1)%qnet < 0.15_dp, 'Cu donates 0.12-0.15 e (first-order k-space SCF: 0.135)')
  call check(abs(res%ef - ef) < 0.008_dp, 'E_F within 8 mRy of k-space H^(2) at the same parameters (depth 16; 2.6 mRy at depth 24)')

  ! (3) third-order Hamiltonian vs the self-consistent first-order k-space ASA
  allocate (ki%cls(2))
  do i = 1, 2
    ki%cls(i) = inp%spec(i)
  end do
  ki%alat = inp%alat; ki%nk = 16; ki%nspin = 1; ki%maxit = 100; ki%tol = 1.0e-6_dp
  kcell%nat = 2; kcell%avec = cell%avec
  allocate (kcell%tau(3, 2), kcell%isp(2)); kcell%tau = cell%tau; kcell%isp = [1, 2]
  call setup_cell(kcell)
  do i = 1, 2; ki%cls(i)%sws = kcell%w * ki%alat; end do
  call asa_scf(ki, kcell, kr)
  inp%order = 3
  call alloy_scf(inp, res3, verbose=.false.)
  print '(a,f10.6,a,f9.5)', ' k-space ASA SCF (nk 16):          E_F=', kr%ef, '  qnet(Cu)=', kr%cls(1)%qnet
  print '(a,f10.6,a,f9.5)', ' recursion order 2:                E_F=', res%ef, '  qnet(Cu)=', res%cls(1)%qnet
  print '(a,f10.6,a,f9.5)', ' recursion order 3:                E_F=', res3%ef, '  qnet(Cu)=', res3%cls(1)%qnet
  call check(res3%converged, 'order-3 recursion SCF converged')
  call check(abs(res3%ef - kr%ef) < 0.006_dp, 'order 3: E_F within 6 mRy of the first-order k-space ASA (order 2 is ~35 mRy off)')
  call check(abs(res3%cls(1)%qnet - kr%cls(1)%qnet) < 0.03_dp, 'order 3: Cu charge transfer within 0.03 e of the k-space ASA')
  call report()
end program test_alloy_b2
