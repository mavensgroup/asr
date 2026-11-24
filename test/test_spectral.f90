!> Bloch spectral function (asr_spectral, asr_edc::edc_green_proj, the target lookup of asr_augcls).
!>  (A) helpers: the k path; the unfolded Bloch sum is exact for a lattice Green function (simple cubic s band on a periodic
!>      cell: A(k) = -Im 1/(z - eps_k) / pi at every commensurate k); the translation test of the decoration (fcc, B2 yes;
!>      L1_2 no); the dressing weights (frozen centre, one species, equal Delta); the Fourier window (values; the ball
!>      window has a non-negative Fourier series on fcc and bcc, the sharp one does not).
!>  (B) readout on a random bcc binary (rank 1, levels 2): the projections of the chains on the target nodes are the same with
!>      the operator reduced by D4h, by D2h and by the identity alone (the same Krylov process in three bases), for every start
!>      orbital that the group allows; the projection on the origin equals g_11; the cubic rotation (x, y, z) -> (y, z, x) maps
!>      G_{px,px}(d) onto G_{pz,pz}(P d) (orbital order s | y z x | xy yz 3z2-r2 xz x2-y2).
program test_spectral
  use asr_kinds, only: dp, pi
  use asr_spectral
  use asr_sublattice, only: decoration_t
  use asr_augspace, only: graph_t
  use asr_augcls
  use asr_edc
  use asr_lmto_alloy, only: orb_alloy_t
  use asr_chain, only: ed_depth_tol, ed_tail, ed_project
  use asr_lmto_alloy, only: build_lmto_alloy
  use asr_symmetry
  use asr_lmto_ham, only: potpar_t
  use asr_testing
  implicit none
  integer, parameter :: lsup = 4, nlev = 24, nper = 6
  real(dp), allocatable :: kv(:, :), dist(:), disp(:, :), v(:)
  integer, allocatable :: icorner(:)
  type(decoration_t) :: dec
  real(dp) :: kp(3, 3), tau(9, 3), dc(9), del(9, 2), err, k(3), ek, gmax, wmin_b, wmin_s
  complex(dp) :: z, gsum
  integer :: i, n, ix, iy, iz, jx, jy, jz, nt

  print '(a)', 'test_spectral'
  ! ---- (A) k path
  kp = reshape([0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.5_dp, 0.0_dp], [3, 3])
  call spectral_kpath(kp, 4, kv, dist, icorner)
  call check(size(dist) == 9 .and. all(icorner == [1, 5, 9]), 'k path: 4 points per segment, corners at 1, 5, 9')
  call check_close(dist(9), 1.5_dp, 1.0e-12_dp, 'k path: length 1 + 0.5')
  call check_close(maxval(abs(kv(:, 9) - kp(:, 3))), 0.0_dp, 1.0e-14_dp, 'k path: last point is the last corner')

  ! ---- (A) unfolded Bloch sum = exact lattice spectral function (sc, nearest neighbours, periodic cell of nper^3 sites)
  z = cmplx(0.13_dp, 0.07_dp, dp)
  n = nper**3; allocate (disp(3, n), v(n))
  i = 0
  do ix = 0, nper - 1; do iy = 0, nper - 1; do iz = 0, nper - 1
        i = i + 1
        disp(:, i) = real([ix, iy, iz] - nper * nint(real([ix, iy, iz], dp) / nper), dp)   ! minimum image (half sites once)
        gsum = 0.0_dp
        do jx = 0, nper - 1; do jy = 0, nper - 1; do jz = 0, nper - 1
              k = real([jx, jy, jz], dp) / nper
              ek = -0.5_dp * (cos(2 * pi * k(1)) + cos(2 * pi * k(2)) + cos(2 * pi * k(3)))
              gsum = gsum + exp(cmplx(0.0_dp, 2 * pi * dot_product(k, disp(:, i)), dp)) / (z - ek)
            end do; end do; end do
        v(i) = -aimag(gsum / n) / pi
      end do; end do; end do
  err = 0.0_dp
  do jx = 0, nper - 1; do jy = 0, nper - 1; do jz = 0, nper - 1
        k = real([jx, jy, jz], dp) / nper
        ek = -0.5_dp * (cos(2 * pi * k(1)) + cos(2 * pi * k(2)) + cos(2 * pi * k(3)))
        err = max(err, abs(spectral_unfold(disp, v, k) + aimag(1.0_dp / (z - ek)) / pi))
      end do; end do; end do
  call check_close(err, 0.0_dp, 1.0e-12_dp, 'unfolded Bloch sum = -Im 1/(z - eps_k) / pi (sc, every commensurate k)')

  ! ---- (A) translations of the decoration
  call set_dec(dec, reshape([0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1], [3, 4]), [1, 1, 1, 1], 1)
  call check(spectral_translations_ok(dec), 'translations: random fcc (one class, 4 sites per cube)')
  call set_dec(dec, reshape([0, 0, 0, 1, 1, 1], [3, 2]), [1, 2], 2)
  call check(spectral_translations_ok(dec), 'translations: B2 (two classes, one site each)')
  call set_dec(dec, reshape([0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1], [3, 4]), [1, 2, 2, 2], 2)
  call check(.not. spectral_translations_ok(dec), 'translations: L1_2 is refused (face sites related by rotations only)')

  ! ---- (A) the Fourier window
  call check(spectral_window('ball', 0.0_dp, 2.0_dp) == 1.0_dp .and. spectral_window('ball', 2.0_dp, 2.0_dp) == 0.0_dp .and. &
             spectral_window('ball', 2.5_dp, 2.0_dp) == 0.0_dp, 'ball window: w(0) = 1, w(rmax) = 0, zero beyond')
  call check_close(spectral_window('ball', 1.0_dp, 2.0_dp), 0.3125_dp, 1.0e-15_dp, 'ball window: w(rmax/2) = 5/16 (two-ball overlap)')
  call check(spectral_window('sharp', 2.0_dp, 2.0_dp) == 1.0_dp .and. spectral_window('sharp', 2.01_dp, 2.0_dp) == 0.0_dp, &
             'sharp window: 1 up to rmax')
  call window_series(1, 'ball', wmin_b); call window_series(1, 'sharp', wmin_s)
  call check(wmin_b > -1.0e-12_dp .and. wmin_s < -0.01_dp, 'fcc, rmax 1, 1.5, 2, 3: Fourier series of the ball window >= 0 at &
             &random k; the sharp one goes negative')
  call window_series(2, 'ball', wmin_b); call window_series(2, 'sharp', wmin_s)
  call check(wmin_b > -1.0e-12_dp .and. wmin_s < -0.01_dp, 'bcc, the same')

  ! ---- (A) dressing weights
  dc = 0.04_dp; del(:, 1) = 0.04_dp; del(:, 2) = 0.09_dp
  call spectral_dressing(dc, del, [1, 2], [0.6_dp, 0.4_dp], .true., tau(:, 1:2))
  call check(all(tau(:, 1) == 1.0_dp) .and. all(tau(:, 2) == 0.0_dp), 'dressing: frozen centre gives tau = (1, 0)')
  call spectral_dressing(dc, del, [2], [1.0_dp], .false., tau(:, 1:1))
  call check_close(maxval(abs(tau(:, 1) - sqrt(0.04_dp / 0.09_dp))), 0.0_dp, 1.0e-14_dp, 'dressing: one species gives T = (D_c/D_s)^1/2')
  del(:, 2) = 0.04_dp
  call spectral_dressing(dc, del, [1, 2], [0.6_dp, 0.4_dp], .false., tau(:, 1:2))
  call check_close(maxval(abs(tau(:, 1) - 1.0_dp)) + maxval(abs(tau(:, 2))), 0.0_dp, 1.0e-14_dp, &
                   'dressing: equal Delta gives tau = (1, 0) (no off-diagonal disorder)')

  ! ---- (B) readout on a random bcc binary
  call readout_test()
  call report()

contains

  subroutine set_dec(d, pts, pcls, ncls)
    type(decoration_t), intent(inout) :: d
    integer, intent(in) :: pts(:, :), pcls(:), ncls
    integer :: ii
    d%m = 2; d%ncls = ncls; d%nsite = size(pcls)
    if (allocated(d%tab)) deallocate (d%tab, d%pts, d%pcls)
    allocate (d%tab(0:1, 0:1, 0:1)); d%tab = 0
    d%pts = pts; d%pcls = pcls
    do ii = 1, size(pcls); d%tab(pts(1, ii), pts(2, ii), pts(3, ii)) = pcls(ii); end do
  end subroutine set_dec

  subroutine readout_test()
    type(potpar_t) :: pp(2)
    real(dp) :: alpha(0:2), d1(9, 9, 1)
    type(graph_t) :: g
    type(orb_alloy_t) :: m0
    type(class_alloy_t) :: mdl
    type(augcls_op) :: op4, op2, op1
    type(edn_t) :: en
    type(edn_work_t) :: w4, w2, w1
    integer, allocatable :: coord(:, :), cls(:), perm(:, :), mats(:, :, :), perm1(:, :), kop(:), tsite(:)
    integer, allocatable :: tn4(:, :), tg4(:, :), tn2(:, :), tg2(:, :), tn1(:, :), tg1(:, :)
    real(dp), allocatable :: d4(:, :, :), d2(:, :, :), tw(:, :, :), tdisp(:, :)
    logical, allocatable :: fl(:)
    complex(dp), allocatable :: gp4(:, :, :), gp2(:, :, :), gp1(:, :, :), gq(:, :)
    complex(dp) :: gc(3), zz
    real(dp) :: e41, e21, e0, erot, pd(3)
    integer :: l, ii, jj, l0, t, t0, tr, nop2, h(3), ns
    integer, parameter :: nsv(2) = [8, 24]
    logical :: oned
    do l = 0, 2
      pp(1)%enu(l) = 0.1_dp * l; pp(1)%c(l) = 0.3_dp - 0.15_dp * l; pp(1)%sqrdel(l) = 0.30_dp - 0.07_dp * l
      pp(1)%p(l) = 0.9_dp - 0.1_dp * l; pp(1)%gamma(l) = 0.25_dp - 0.1_dp * l
      pp(2)%enu(l) = 0.05_dp * l; pp(2)%c(l) = 0.45_dp - 0.1_dp * l; pp(2)%sqrdel(l) = 0.24_dp - 0.06_dp * l
      pp(2)%p(l) = 0.8_dp - 0.1_dp * l; pp(2)%gamma(l) = 0.2_dp - 0.08_dp * l
    end do
    pp(1)%alpha_rep = pp(1)%gamma; pp(2)%alpha_rep = pp(2)%gamma
    alpha = [0.3485_dp, 0.05303_dp, 0.010714_dp]
    call build_lmto_alloy('bcc', lsup, 1.2_dp, 0.7_dp, pp(1), pp(2), 0.5_dp, g, coord, m0, alpha_in=alpha)
    mdl%nl = 9; mdl%nsp = 2; mdl%ncls = 1
    allocate (mdl%nsc(1), mdl%spc(2, 1), mdl%conc(2, 1))
    mdl%nsc = 2; mdl%spc(:, 1) = [1, 2]; mdl%conc(:, 1) = [0.6_dp, 0.4_dp]
    call augcls_set_blocks(mdl, m0%block, m0%block0)
    call edc_setup(mdl, pp, alpha, 1, en)
    allocate (cls(g%n), fl(g%n)); cls = 1; fl = .true.; fl(1) = .false.
    call class_ops(coord, 2 * lsup, cls, perm, mats)
    call orbital_rep(mats, d4)
    allocate (kop(size(perm, 2))); nop2 = 0
    do jj = 1, size(perm, 2)
      if (count(mats(:, :, jj) /= 0) == 3 .and. all([(mats(ii, ii, jj) /= 0, ii=1, 3)])) then; nop2 = nop2 + 1; kop(nop2) = jj; end if
    end do
    call orbital_rep(mats(:, :, kop(1:nop2)), d2)
    allocate (perm1(g%n, 1)); perm1(:, 1) = [(ii, ii=1, g%n)]
    d1 = 0.0_dp; do l = 1, 9; d1(l, l, 1) = 1.0_dp; end do
    call check(size(perm, 2) == 16 .and. nop2 == 8, 'bcc site group: D4h (16 operations), D2h (8)')
    ! targets within 1 alat (origin, 8 + 6 neighbours)
    nt = 0; allocate (tsite(g%n), tdisp(3, g%n))
    do ii = 1, g%n
      h = modulo(coord(:, ii) + lsup, 2 * lsup) - lsup
      if (sum(h**2) <= 4) then; nt = nt + 1; tsite(nt) = ii; tdisp(:, nt) = 0.5_dp * h; end if
    end do
    tsite = tsite(1:nt); tdisp = tdisp(:, 1:nt)
    call check(nt == 15, 'targets: origin and the first two bcc shells (15 sites)')
    call build_augcls_operator(g, mdl, cls, perm, d4, 1, 2, op4, maxrank=1, flippable=fl, tsite=tsite, tnode=tn4, tgop=tg4)
    call build_augcls_operator(g, mdl, cls, perm(:, kop(1:nop2)), d2, 1, 2, op2, maxrank=1, flippable=fl, tsite=tsite, &
                               tnode=tn2, tgop=tg2)
    call build_augcls_operator(g, mdl, cls, perm1, d1, 1, 2, op1, maxrank=1, flippable=fl, tsite=tsite, tnode=tn1, tgop=tg1)
    print '(a,3(i0,1x))', '  nodes D4h / D2h / identity: ', op4%nn, op2%nn, op1%nn
    t0 = 0; do t = 1, nt; if (tsite(t) == 1) t0 = t; end do
    call check(all(pack(tn1, spread([(t /= t0, t=1, nt)], 1, 2)) > 0) .and. tn1(1, t0) > 0 .and. tn1(2, t0) == 0 .and. &
               all(pack(tn2, spread([(t /= t0, t=1, nt)], 1, 2)) > 0) .and. all(pack(tn4, spread([(t /= t0, t=1, nt)], 1, 2)) > 0), &
               'targets: both local states of every neighbour are in the node set; the frozen origin has the reference state only')
    allocate (tw(9, 2, nt))
    do t = 1, nt
      call spectral_dressing(en%dc, en%del, [1, 2], [0.6_dp, 0.4_dp], tsite(t) == 1, tw(:, :, t))
    end do
    call augcls_set_params(op4, mdl); call augcls_freeze_site(op4, 1, 1)
    call augcls_set_params(op2, mdl); call augcls_freeze_site(op2, 1, 1)
    call augcls_set_params(op1, mdl); call augcls_freeze_site(op1, 1, 1)
    call edn_work_init(w4, op4%dim()); call edn_work_init(w2, op2%dim()); call edn_work_init(w1, op1%dim())
    ed_depth_tol = 0.0_dp; ed_tail = 'none'; ed_project = .true.        ! the tail-free Galerkin value: no eigenvalues of the chain
    allocate (gp4(9, nt, 9), gp2(9, nt, 9), gp1(9, nt, 9), gq(9, nt)); gp4 = 0.0_dp
    zz = cmplx(0.35_dp, 0.08_dp, dp)
    do ns = 1, 2
    e41 = 0.0_dp; e21 = 0.0_dp; e0 = 0.0_dp
    do l0 = 1, 9
      oned = all(abs(abs([(d4(l0, l0, ii), ii=1, size(perm, 2))]) - 1.0_dp) < 1.0e-9_dp)
      call augcls_set_start(op2, l0, d2); call augcls_set_start(op1, l0, d1)
      call edc_green_proj(op2, w2, en, zz, l0, nsv(ns), tn2, tg2, tw, gc(2), gp2(:, :, l0), gq)
      call edc_green_proj(op1, w1, en, zz, l0, nsv(ns), tn1, tg1, tw, gc(1), gp1(:, :, l0), gq)
      e21 = max(e21, maxval(abs(gp2(:, :, l0) - gp1(:, :, l0))), abs(gc(2) - gc(1)))
      e0 = max(e0, abs(gp2(l0, t0, l0) - gc(2)))
      if (oned) then
        call augcls_set_start(op4, l0, d4)
        call edc_green_proj(op4, w4, en, zz, l0, nsv(ns), tn4, tg4, tw, gc(3), gp4(:, :, l0), gq)
        e41 = max(e41, maxval(abs(gp4(:, :, l0) - gp1(:, :, l0))))
      end if
    end do
    gmax = 0.0_dp
    do t = 1, nt
      if (t /= t0) gmax = max(gmax, maxval(abs(gp1(:, t, :))))
    end do
    print '(a,i0,a,es10.2,a,es10.2,a,es10.2,a,es10.2)', '  steps ', nsv(ns), ': max |G_0R|, R /= 0: ', gmax, ';  D2h - identity ', &
      e21, ';  D4h - identity ', e41, ';  origin - g_11 ', e0
    if (ns == 1) then          ! short chain: the three bases agree to rounding (3e-12 at 8 steps, 6e-10 at 24: it grows along the chain)
      call check(gmax > 1.0e-3_dp, 'off-site elements are not trivially zero')
      call check_close(e21, 0.0_dp, 1.0e-11_dp, '8 steps: D2h-reduced projections = unreduced projections (all 9 start orbitals)')
      call check_close(e41, 0.0_dp, 1.0e-11_dp, '8 steps: D4h-reduced = unreduced (s, p_z, d_xy, d_3z2, d_x2-y2)')
      call check_close(e0, 0.0_dp, 1.0e-13_dp, '8 steps: projection on the origin = g_11 of the chain')
    else                       ! 24 steps: rounding grows along the chain (no reorthogonalisation), still far below the chain error
      call check_close(e21 / gmax, 0.0_dp, 1.0e-8_dp, '24 steps: D2h = unreduced to rounding (relative)')
      call check_close(e41 / gmax, 0.0_dp, 1.0e-8_dp, '24 steps: D4h = unreduced to rounding (relative)')
      call check_close(e0, 0.0_dp, 1.0e-9_dp, '24 steps: projection on the origin = g_11 to rounding')
    end if
    end do
    ! cubic rotation P: (d1, d2, d3) -> (d2, d3, d1) maps p_x onto p_z: G_{px,px}(d) = G_{pz,pz}(P d)
    erot = 0.0_dp
    do t = 1, nt
      pd = [tdisp(2, t), tdisp(3, t), tdisp(1, t)]
      tr = 0
      do jj = 1, nt
        if (all(abs(tdisp(:, jj) - pd) < 1.0e-9_dp)) tr = jj
      end do
      erot = max(erot, abs(gp2(4, t, 4) - gp2(3, tr, 3)))
    end do
    call check_close(erot, 0.0_dp, 1.0e-11_dp, 'cubic rotation: G_{px,px}(d) = G_{pz,pz}(P d)')
  end subroutine readout_test

  !> min over random k and rmax of sum_R w(|R|) cos(2 pi k . R) / sum_R w(|R|), R on the fcc (lat 1) or bcc (lat 2) lattice
  !> (units alat).  Ball window: >= 0 (asr_spectral header, Poisson summation of a positive-definite function).
  subroutine window_series(lat, kind, wmin)
    integer, intent(in) :: lat
    character(len=*), intent(in) :: kind
    real(dp), intent(out) :: wmin
    real(dp), parameter :: rm(4) = [1.0_dp, 1.5_dp, 2.0_dp, 3.0_dp]
    real(dp), allocatable :: rr(:, :), ww(:)
    real(dp) :: kk(3), a, s0, h(3)
    integer :: ir, i1, i2, i3, n_, ik_
    wmin = huge(1.0_dp)
    do ir = 1, size(rm)
      allocate (rr(3, 20000), ww(20000)); n_ = 0
      do i1 = -8, 8
        do i2 = -8, 8
          do i3 = -8, 8
            h = 0.5_dp * real([i1, i2, i3], dp)                   ! half-lattice points; keep the lattice ones
            if (lat == 1 .and. mod(i1 + i2 + i3, 2) /= 0) cycle  ! fcc: even sum
            if (lat == 2 .and. (mod(i1 - i2, 2) /= 0 .or. mod(i2 - i3, 2) /= 0)) cycle   ! bcc: all even or all odd
            if (norm2(h) > rm(ir) + 1.0e-9_dp) cycle
            n_ = n_ + 1; rr(:, n_) = h; ww(n_) = spectral_window(kind, norm2(h), rm(ir))
          end do
        end do
      end do
      s0 = sum(ww(1:n_))
      do ik_ = 1, 400
        call random_number(kk); kk = 2.0_dp * kk - 1.0_dp
        a = spectral_unfold(rr(:, 1:n_), ww(1:n_), kk)
        wmin = min(wmin, a / s0)
      end do
      deallocate (rr, ww)
    end do
  end subroutine window_series
end program test_spectral
