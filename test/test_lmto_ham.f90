!> Ordered fcc TB-LMTO Hamiltonian from the Stuttgart converged Ag parameters.
!> Independent consistency checks (cubic symmetry, Hermiticity) and the
!> reference-DOS comparison that IS verified (d band); the sp region is NOT
!> asserted (known open difference: combined correction / f downfolding).
program test_lmto_ham
  use asr_kinds, only: dp
  use asr_strconst
  use asr_lmto_ham
  use asr_testing
  implicit none
  integer, parameter :: ne = 801
  type(potpar_t) :: pp
  real(dp), allocatable :: vec(:, :), blk(:, :, :)
  real(dp) :: alpha(0:2), al4(0:3), w, egrid(ne), de, dos(ne, 0:2), ref(ne), bvec(3, 3), ex
  real(dp) :: ev(9), rwork(100), nd, cd, nr, cr
  complex(dp) :: h(9, 9), cw(200)
  integer :: i, u, info
  real(dp) :: occ(0:2), nef
  external :: zheev

  print '(a)', 'test_lmto_ham'
  pp%enu = [-0.54494317_dp, -0.45471939_dp, -0.51341327_dp]
  pp%c = [-0.46144639_dp, 0.41368602_dp, -0.51457447_dp]
  pp%sqrdel = [0.36457406_dp, -0.36517716_dp, -0.10341022_dp]
  pp%p = [0.08230178_dp, 0.03623364_dp, 2.08144602_dp]
  pp%gamma = [0.42477526_dp, 0.11537965_dp, 0.00550351_dp]
  pp%alpha_rep = pp%gamma          ! Stuttgart ctrl stores the gamma representation
  w = avg_ws_radius('fcc'); al4 = hardcore_alpha(0.7_dp); alpha = al4(0:2)
  call lattice_cluster('fcc', 1.55_dp, vec)
  call screened_cluster(vec, w, alpha, blk)

  ! Gamma point: cubic degeneracies p(3), t2g(3), eg(2); Hermiticity
  call ham_k([0.0_dp, 0.0_dp, 0.0_dp], vec, blk, pp, alpha, 1.0_dp, .true., h)
  call check_close(maxval(abs(h - conjg(transpose(h)))), 0.0_dp, 1.0e-12_dp, 'H(Gamma) Hermitian')
  call zheev('V', 'U', 9, h, 9, ev, cw, 200, rwork, info)
  print '(a,9f9.4)', '    Gamma eigenvalues: ', ev
  call check_close(real(count(abs(ev(2:9) - ev(1:8)) < 1.0e-9_dp), dp), 5.0_dp, 0.5_dp, 'Gamma: 3(p)+2(t2g)+1(eg) degeneracy pairs (5 equal neighbours)')
  call ham_k([0.31_dp, -0.17_dp, 0.52_dp], vec, blk, pp, alpha, 1.0_dp, .true., h)
  call check_close(maxval(abs(h - conjg(transpose(h)))), 0.0_dp, 1.0e-12_dp, 'H(k) Hermitian at generic k')

  ! reference DOS: d window
  bvec = reshape([-1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp], [3, 3])
  de = 2.6_dp / (ne - 1)
  do i = 1, ne; egrid(i) = -0.9_dp + (i - 1) * de; end do
  call dos_on_grid(24, bvec, vec, blk, pp, alpha, 1.0_dp, .true., egrid, 0.004_dp, dos)
  dos = 2.0_dp * dos
  open (newunit=u, file='data/reference/ag_fcc_stuttgart/dos_d.dat', status='old')
  do i = 1, ne; read (u, *) ex, ref(i); end do
  close (u)
  ref = ref / 4.0_dp
  nd = sum(dos(:, 2), mask=egrid > -0.85_dp .and. egrid < -0.25_dp) * de
  nr = sum(ref, mask=egrid > -0.85_dp .and. egrid < -0.25_dp) * de
  cd = sum(egrid * dos(:, 2), mask=egrid > -0.85_dp .and. egrid < -0.25_dp) * de / nd
  cr = sum(egrid * ref, mask=egrid > -0.85_dp .and. egrid < -0.25_dp) * de / nr
  print '(a,2f8.3,a,2f9.4)', '    d window weight (ours, ref): ', nd, nr, '   centroid (ours, ref): ', cd, cr
  call check_close(nd / nr, 1.0_dp, 0.02_dp, 'd-band weight within 2% of Stuttgart')
  call check_close(cd, cr, 0.002_dp, 'd-band centroid within 2 mRy of Stuttgart')
  ! l-resolved occupation up to the Stuttgart E_F, and N(E_F)  (E < +0.3 Ry only: above that
  ! the second-order H without combined correction is known to differ)
  do i = 0, 2
    occ(i) = sum(dos(:, i), mask=egrid <= -0.18641407_dp) * de
  end do
  nef = sum(sum(dos, dim=2), mask=abs(egrid + 0.18641407_dp) < 0.015_dp) / real(count(abs(egrid + 0.18641407_dp) < 0.015_dp), dp)
  call check_close(occ(0), 0.6834_dp, 0.015_dp, 's occupation at E_F vs Stuttgart (0.683)')
  call check_close(occ(1), 0.6483_dp, 0.025_dp, 'p occupation at E_F vs Stuttgart (0.648)')
  call check_close(occ(2), 9.5645_dp, 0.10_dp, 'd occupation at E_F vs Stuttgart (9.565)')
  call check_close(nef, 3.177_dp, 0.10_dp, 'DOS(E_F) per atom vs Stuttgart (3.177 /Ry)')
  call report()
end program test_lmto_ham
