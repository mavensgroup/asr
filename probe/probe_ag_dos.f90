!> Compare our ordered fcc LMTO DOS (Stuttgart converged Ag parameters from
!> AMscf/lmto/fcc/ctrl) with the Stuttgart LMASA DOS (Ag-s/p/d).  Diagnostic
!> program: prints; writes trial/ag_dos_<variant>.dat.
program probe_ag_dos
  use asr_kinds, only: dp
  use asr_strconst
  use asr_lmto_ham
  implicit none
  integer, parameter :: ne = 801
  type(potpar_t) :: pp
  real(dp), allocatable :: vec(:, :), blk(:, :, :)
  real(dp) :: alpha(0:2), egrid(ne), dos(ne, 0:2), w, bvec(3, 3), de, ef, occ(0:2), ssg
  real(dp) :: ref(ne, 0:3), ex
  integer :: i, l, u, iv, is
  logical :: ab
  character(len=8) :: names(0:3) = ["s.dat   ", "p.dat   ", "d.dat   ", "f.dat   "]
  character(len=64) :: fn
  real(dp) :: dmax, dd, al4(0:3)

  pp%enu = [-0.54494317_dp, -0.45471939_dp, -0.51341327_dp]
  pp%c = [-0.46144639_dp, 0.41368602_dp, -0.51457447_dp]
  pp%sqrdel = [0.36457406_dp, -0.36517716_dp, -0.10341022_dp]
  pp%p = [0.08230178_dp, 0.03623364_dp, 2.08144602_dp]
  pp%gamma = [0.42477526_dp, 0.11537965_dp, 0.00550351_dp]
  pp%alpha_rep = pp%gamma          ! Stuttgart ctrl stores the gamma representation

  w = avg_ws_radius('fcc')
  al4 = hardcore_alpha(0.7_dp); alpha = al4(0:2)
  call lattice_cluster('fcc', 1.55_dp, vec)
  call screened_cluster(vec, w, alpha, blk)
  bvec = reshape([-1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp], [3, 3])
  de = 2.6_dp / (ne - 1)
  do i = 1, ne; egrid(i) = -0.9_dp + (i - 1) * de; end do
  ef = -0.18641407_dp

  do l = 0, 3
    open (newunit=u, file='/home/rudra/Devel/ASR/asr_new/data/reference/ag_fcc_stuttgart/dos_'//trim(names(l)), status='old')
    do i = 1, ne; read (u, *) ex, ref(i, l); end do
    close (u)
  end do
  ref = ref / 4.0_dp          ! per atom, both spins
  print '(a,3f8.4)', 'reference occupations s,p,d per atom (E<EF): ', &
    (sum(ref(:, l), mask=egrid <= ef) * de, l=0, 2)
  print '(a)', 'variant                         occ_s   occ_p   occ_d   Ntot   max|dDOS|(spd,per atom)'
  do iv = 1, 4
    ssg = merge(1.0_dp, -1.0_dp, iv <= 2); ab = (mod(iv, 2) == 0)
    call dos_on_grid(24, bvec, vec, blk, pp, alpha, ssg, ab, egrid, 0.004_dp, dos)
    dos = 2.0_dp * dos              ! both spins, like the reference
    do l = 0, 2
      occ(l) = sum(dos(:, l), mask=egrid <= ef) * de
    end do
    dmax = 0.0_dp
    do i = 1, ne
      dd = abs(sum(dos(i, :)) - sum(ref(i, 0:2)))
      dmax = max(dmax, dd)
    end do
    write (fn, '(a,i0,a)') '/home/rudra/Devel/ASR/asr_new/trial/ag_dos_v', iv, '.dat'
    open (newunit=u, file=trim(fn), status='replace')
    write (u, '(a)') '# E  ours_s ours_p ours_d  ref_s ref_p ref_d'
    do i = 1, ne
      write (u, '(7es14.6)') egrid(i), dos(i, :), ref(i, 0:2)
    end do
    close (u)
    print '(a,i1,a,f4.0,a,l1,a,4f8.3,f9.2)', 'v', iv, ' sign S=', ssg, ' abs(sqrdel)=', ab, '  ', occ, sum(occ), dmax
  end do
  print '(a,f9.2)', 'ref peak total (per atom): ', maxval(sum(ref(:, 0:2), dim=2))
end program probe_ag_dos
