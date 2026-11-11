!> Site-class operator asr_augcls + energy-dependent route asr_edc:
!>  (A) ring of 4 sites, three classes (ordered species 1 | binary {2,3} | ternary {1,3,4}, species shared between classes), start
!>      site frozen to each allowed species: recursion == exact configuration average of [E - H(E)]^-1_00 (identity group);
!>  (B) decorated bcc (L2_1 pattern: two rock-salt sublattices + the body-centre sublattice; ordered | ternary | ordered), origin on
!>      each sublattice (O_h and T_d site groups): class-reduced chain == unreduced chain, all four start orbitals;
!>  (C) regressions, one class: binary == asr_augsym/asr_ed (same D4h group, same node count), ternary == unreduced asr_augn/asr_edn;
!>  (D) short-range order (first-shell class of the frozen start site with conditional concentrations, augcls_set_class_conc):
!>      ring of 5 == exact average over the conditional product measure, including species of zero conditional weight; fcc ternary:
!>      reduced == unreduced with the shell class, and shell concentrations = x reproduce the operator without a shell class.
program test_augcls
  use asr_kinds, only: dp
  use asr_augspace, only: graph_t
  use asr_augcls
  use asr_edc
  use asr_augorb, only: orb_alloy_t
  use asr_augsym
  use asr_ed, only: asr_ed_t, asr_ed_setup, asr_ed_green
  use asr_endep, only: ed_work_t, ed_work_init
  use asr_augn, only: multi_alloy_t, augn_op, build_augn_operator, augn_freeze_site
  use asr_edn, only: edn_setup, edn_green
  use asr_lmto_alloy, only: build_lmto_alloy
  use asr_symmetry
  use asr_lmto_ham, only: potpar_t, to_alpha, o_alpha
  use asr_testing
  implicit none
  integer, parameter :: orb_l(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2], l0s(4) = [1, 3, 5, 7]
  type(potpar_t) :: pp(4), qq(4)
  real(dp) :: alpha(0:2), oo(0:2, 4)
  integer :: l, sp
  character(len=16) :: arg
  logical :: ringonly = .false.

  print '(a)', 'test_augcls'
  do l = 0, 2
    pp(1)%enu(l) = 0.1_dp * l; pp(1)%c(l) = 0.3_dp - 0.15_dp * l; pp(1)%sqrdel(l) = 0.30_dp - 0.07_dp * l
    pp(1)%p(l) = 0.9_dp - 0.1_dp * l; pp(1)%gamma(l) = 0.25_dp - 0.1_dp * l
    pp(2)%enu(l) = 0.05_dp * l; pp(2)%c(l) = 0.45_dp - 0.1_dp * l; pp(2)%sqrdel(l) = 0.27_dp - 0.06_dp * l
    pp(2)%p(l) = 0.8_dp - 0.1_dp * l; pp(2)%gamma(l) = 0.2_dp - 0.08_dp * l
    pp(3)%enu(l) = -0.07_dp + 0.08_dp * l; pp(3)%c(l) = 0.2_dp - 0.05_dp * l; pp(3)%sqrdel(l) = 0.22_dp - 0.04_dp * l
    pp(3)%p(l) = 0.7_dp - 0.1_dp * l; pp(3)%gamma(l) = 0.15_dp - 0.06_dp * l
    pp(4)%enu(l) = 0.02_dp + 0.03_dp * l; pp(4)%c(l) = 0.6_dp - 0.2_dp * l; pp(4)%sqrdel(l) = 0.25_dp - 0.05_dp * l
    pp(4)%p(l) = 0.85_dp - 0.1_dp * l; pp(4)%gamma(l) = 0.22_dp - 0.07_dp * l
  end do
  do sp = 1, 4; pp(sp)%alpha_rep = pp(sp)%gamma; end do
  alpha = [0.3485_dp, 0.05303_dp, 0.010714_dp]
  do sp = 1, 4
    qq(sp) = to_alpha(pp(sp), alpha); oo(:, sp) = o_alpha(qq(sp), alpha)
  end do
  ! optional argument "ring": only the two ring checks (seconds instead of minutes; compiler checks on the cluster)
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg)
    ringonly = trim(arg) == 'ring'
  end if
  call ring_brute_force()
  if (.not. ringonly) call decorated_bcc()
  if (.not. ringonly) call regressions()
  call ring_sro()
  if (.not. ringonly) call fcc_sro()
  call report()

contains

  !> three classes: 1 = ordered {1}, 2 = binary {2,3} (0.4, 0.6), 3 = ternary {1,3,4} (0.3, 0.5, 0.2)
  subroutine three_classes(m)
    type(class_alloy_t), intent(inout) :: m
    m%nl = 9; m%nsp = 4; m%ncls = 3
    allocate (m%nsc(3), m%spc(3, 3), m%conc(3, 3)); m%spc = 0; m%conc = 0.0_dp
    m%nsc = [1, 2, 3]
    m%spc(1:1, 1) = [1];       m%conc(1:1, 1) = [1.0_dp]
    m%spc(1:2, 2) = [2, 3];    m%conc(1:2, 2) = [0.4_dp, 0.6_dp]
    m%spc(1:3, 3) = [1, 3, 4]; m%conc(1:3, 3) = [0.3_dp, 0.5_dp, 0.2_dp]
  end subroutine three_classes

  subroutine ring_brute_force()
    integer, parameter :: n = 4, nd = 9 * n, nlev = 400
    integer, parameter :: rcls(n) = [1, 2, 3, 2]
    type(graph_t) :: g
    type(class_alloy_t) :: mdl
    type(augcls_op) :: op
    type(edn_t) :: ed
    type(edn_work_t) :: wk
    integer :: perm(n, 1), i, l2, ic, cen, kz, m, cf(n), q, c, ncfg, j, ib, l0, kc, nseed
    integer, allocatable :: seed(:)
    real(dp) :: dmat(9, 9, 1), u, pc, err, psum, perr
    complex(dp) :: z, gre, gden, a(nlev), b(nlev)
    logical :: fl(n)
    allocate (g%nbr_ptr(n + 1), g%nbr(2 * n)); g%n = n
    do i = 1, n
      g%nbr_ptr(i) = 2 * (i - 1) + 1
      g%nbr(2 * i - 1) = mod(i, n) + 1
      g%nbr(2 * i) = mod(i - 2 + n, n) + 1
    end do
    g%nbr_ptr(n + 1) = 2 * n + 1
    perr = 0.0_dp
    call three_classes(mdl)
    allocate (mdl%block(9, 9, 2), mdl%block0(9, 9))
    call random_seed(size=nseed); allocate (seed(nseed)); seed = 20261111 + 7 * [(i, i=1, nseed)]
    call random_seed(put=seed)                           ! fixed seed: a failure can be reproduced
    do l = 1, 9; do l2 = 1, 9
      call random_number(u); mdl%block(l, l2, 1) = 0.5_dp * (u - 0.5_dp)
    end do; end do
    mdl%block(:, :, 2) = transpose(mdl%block(:, :, 1))
    do l = 1, 9; do l2 = l, 9
      call random_number(u); mdl%block0(l, l2) = 0.6_dp * (u - 0.5_dp); mdl%block0(l2, l) = mdl%block0(l, l2)
    end do; end do
    perm(:, 1) = [(i, i=1, n)]
    dmat = 0.0_dp; do l = 1, 9; dmat(l, l, 1) = 1.0_dp; end do
    err = 0.0_dp
    do ic = 1, 3                                         ! start site = first ring site of class ic (rotate the ring labels)
      ib = findloc(rcls, ic, 1)
      do j = 1, mdl%nsc(ic)
        cen = mdl%spc(j, ic)
        call edc_setup(mdl, pp, alpha, cen, ed)
        fl = .true.; fl(1) = .false.
        call build_augcls_operator(g, mdl, cshift(rcls, ib - 1), perm, dmat, 1, 80, op, flippable=fl)
        call augcls_freeze_site(op, 1, cen)
        call edn_work_init(wk, op%dim())
        do l0 = 1, 9, 4
          do kz = 1, 3
            z = cmplx(-0.3_dp + 0.3_dp * kz, 0.3_dp, dp)
            call edc_green(op, wk, ed, z, l0, nlev, gre, a, b, m)
            gre = cf_complex(a, b, m, z)
            gden = (0.0_dp, 0.0_dp)
            ncfg = 1                                       ! plain loop: the implied-do form gave a wrong count with gfortran 14
            do i = 2, n; ncfg = ncfg * mdl%nsc(rcls(1 + mod(i + ib - 2, n))); end do
            psum = 0.0_dp
            do c = 0, ncfg - 1
              q = c; pc = 1.0_dp; cf(1) = cen
              do i = 2, n
                kc = rcls(1 + mod(i + ib - 2, n))      ! a plain variable: gfortran 14 builds failed this check (err 18.7) with an associate here
                cf(i) = mdl%spc(mod(q, mdl%nsc(kc)) + 1, kc); pc = pc * mdl%conc(mod(q, mdl%nsc(kc)) + 1, kc); q = q / mdl%nsc(kc)
              end do
              gden = gden + pc * dense_g(mdl, g, cf, cen, z, l0)
              psum = psum + pc
            end do
            perr = max(perr, abs(psum - 1.0_dp))
            err = max(err, abs(gre - gden))
            if (abs(gre - gden) > 1.0e-6_dp) print '(a,3i3,i5,a,4es12.4)', '    case (class, centre, orbital), steps ', ic, cen, l0, m, &
              ': recursion, exact ', gre, gden
          end do
        end do
        print '(a,i0,a,i0,a,i0,a)', '  ring: start class ', ic, ', central species ', cen, ': ', op%nn, ' nodes'
      end do
    end do
    print '(a,es10.3)', '  ring: max |recursion - exact average| ', err
    call check_close(perr, 0.0_dp, 1.0e-12_dp, 'three classes: configuration weights of the exact average sum to 1')
    call check_close(err, 0.0_dp, 1.0e-6_dp, 'three classes (ordered/binary/ternary): ED recursion == exact configuration average')
  end subroutine ring_brute_force

  !> [E - H(E)]^-1_{l0,l0} at site 1 for one configuration (species per site), H(E) = E - T (E - E_nu - h) T - T h o h T
  complex(dp) function dense_g(mdl, g, cfg, ce, zz, l0)
    type(class_alloy_t), intent(in) :: mdl
    type(graph_t), intent(in) :: g
    integer, intent(in) :: cfg(:), ce, l0
    complex(dp), intent(in) :: zz
    integer :: n, nd, ii, jj, kk, ll, kq, info
    real(dp), allocatable :: s(:, :), hm(:, :), hoh(:, :), dl(:), cm(:), en(:), om(:), dc(:), tt(:)
    complex(dp), allocatable :: ae(:, :), rhs(:)
    integer, allocatable :: ipiv(:)
    external :: zgesv
    n = g%n; nd = 9 * n
    allocate (s(nd, nd), hm(nd, nd), hoh(nd, nd), dl(nd), cm(nd), en(nd), om(nd), dc(nd), tt(nd), ae(nd, nd), rhs(nd), ipiv(nd))
    s = 0.0_dp
    do ii = 1, n
      s(9 * (ii - 1) + 1:9 * ii, 9 * (ii - 1) + 1:9 * ii) = mdl%block0
      do kk = 1, g%nbr_ptr(ii + 1) - g%nbr_ptr(ii)
        jj = g%nbr(g%nbr_ptr(ii) + kk - 1)
        s(9 * (ii - 1) + 1:9 * ii, 9 * (jj - 1) + 1:9 * jj) = s(9 * (ii - 1) + 1:9 * ii, 9 * (jj - 1) + 1:9 * jj) + mdl%block(:, :, kk)
      end do
    end do
    do ii = 1, n
      do ll = 1, 9
        kq = 9 * (ii - 1) + ll
        dl(kq) = qq(cfg(ii))%sqrdel(orb_l(ll))**2
        cm(kq) = qq(cfg(ii))%c(orb_l(ll)) - qq(cfg(ii))%enu(orb_l(ll))
        en(kq) = qq(cfg(ii))%enu(orb_l(ll)); om(kq) = oo(orb_l(ll), cfg(ii))
        dc(kq) = qq(ce)%sqrdel(orb_l(ll))**2
      end do
    end do
    do ii = 1, nd
      do jj = 1, nd; hm(ii, jj) = sqrt(dl(ii)) * s(ii, jj) * sqrt(dl(jj)); end do
      hm(ii, ii) = hm(ii, ii) + cm(ii)
    end do
    do jj = 1, nd; do ii = 1, nd; hoh(ii, jj) = sum(hm(ii, :) * om * hm(:, jj)); end do; end do
    tt = sqrt(dc / dl)
    do jj = 1, nd; do ii = 1, nd
      ae(ii, jj) = tt(ii) * tt(jj) * (merge(zz - en(ii), (0.0_dp, 0.0_dp), ii == jj) - hm(ii, jj) + hoh(ii, jj))
    end do; end do
    rhs = (0.0_dp, 0.0_dp); rhs(l0) = (1.0_dp, 0.0_dp)
    call zgesv(nd, 1, ae, nd, ipiv, rhs, nd, info)
    if (info /= 0) error stop 'dense_g: zgesv failed'
    dense_g = rhs(l0)
  end function dense_g

  complex(dp) function cf_complex(aa, bb, mm, zz)
    complex(dp), intent(in) :: aa(:), bb(:), zz
    integer, intent(in) :: mm
    complex(dp) :: t
    integer :: j
    t = (0.0_dp, 0.0_dp)
    do j = mm, 1, -1
      t = 1.0_dp / (zz - aa(j) - t)
      if (j > 1) t = bb(j - 1)**2 * t
    end do
    cf_complex = t
  end function cf_complex

  !> L2_1 pattern on the bcc graph (half-lattice units, period 4): all odd -> class 3; all even with (x+y+z)/2 even -> class 1,
  !> odd -> class 2; coordinates shifted by `sh` so that the origin sits on the chosen sublattice.
  subroutine l21_classes(coord, sh, cls)
    integer, intent(in) :: coord(:, :), sh(3)
    integer, allocatable, intent(out) :: cls(:)
    integer :: i, c(3)
    allocate (cls(size(coord, 2)))
    do i = 1, size(coord, 2)
      c = modulo(coord(:, i) + sh, 4)
      if (modulo(c(1), 2) == 1) then
        cls(i) = 3
      else if (modulo(sum(c) / 2, 2) == 0) then
        cls(i) = 1
      else
        cls(i) = 2
      end if
    end do
  end subroutine l21_classes

  !> classes for the decorated bcc test: 1 = ordered {1}, 2 = ternary {2,3,4} (0.3, 0.5, 0.2), 3 = ordered {3}
  subroutine heusler_like(m)
    type(class_alloy_t), intent(inout) :: m
    m%nl = 9; m%nsp = 4; m%ncls = 3
    allocate (m%nsc(3), m%spc(3, 3), m%conc(3, 3)); m%spc = 0; m%conc = 0.0_dp
    m%nsc = [1, 3, 1]
    m%spc(1:1, 1) = [1];       m%conc(1:1, 1) = [1.0_dp]
    m%spc(1:3, 2) = [2, 3, 4]; m%conc(1:3, 2) = [0.3_dp, 0.5_dp, 0.2_dp]
    m%spc(1:1, 3) = [3];       m%conc(1:1, 3) = [1.0_dp]
  end subroutine heusler_like

  subroutine decorated_bcc()
    integer, parameter :: nlev = 3, shifts(3, 3) = reshape([0, 0, 0, 2, 0, 0, 1, 1, 1], [3, 3])
    type(graph_t) :: g
    type(orb_alloy_t) :: m0
    type(class_alloy_t) :: mdl
    type(augcls_op) :: opr, opu
    type(edn_t) :: ed
    type(edn_work_t) :: w1, w2
    integer, allocatable :: coord(:, :), cls(:), perm(:, :), mats(:, :, :), pid(:, :)
    real(dp), allocatable :: d(:, :, :), did(:, :, :)
    logical, allocatable :: fl(:)
    logical :: cubic
    complex(dp) :: z, g1, g2, a1(nlev), b1(nlev), a2(nlev), b2(nlev)
    real(dp) :: err, errg
    integer :: ish, ic, j, cen, io, kz, m1, m2, i
    call build_lmto_alloy('bcc', 4, 1.05_dp, 0.7_dp, pp(1), pp(2), 0.5_dp, g, coord, m0, alpha_in=alpha)
    call heusler_like(mdl)
    mdl%block = m0%block; mdl%block0 = m0%block0
    allocate (fl(g%n)); fl = .true.; fl(1) = .false.
    allocate (pid(g%n, 1), did(9, 9, 1)); pid(:, 1) = [(i, i=1, g%n)]; did = 0.0_dp
    do i = 1, 9; did(i, i, 1) = 1.0_dp; end do
    err = 0.0_dp; errg = 0.0_dp
    do ish = 1, 3
      call l21_classes(coord, shifts(:, ish), cls)
      ic = cls(1)
      call class_ops(coord, 8, cls, perm, mats, cubic)
      call orbital_rep(mats, d)
      print '(a,i0,a,i0,a,l1)', '  decorated bcc: origin class ', ic, ', group order ', size(perm, 2), ', cubic site ', cubic
      call check(cubic .and. size(perm, 2) == merge(8, 16, ic == 3), 'class group: D4h on the rock-salt sublattices, D2d on the T_d one')
      do j = 1, mdl%nsc(ic)
        cen = mdl%spc(j, ic)
        call edc_setup(mdl, pp, alpha, cen, ed)
        call build_augcls_operator(g, mdl, cls, perm, d, 1, nlev, opr, maxrank=1, flippable=fl)
        call build_augcls_operator(g, mdl, cls, pid, did, 1, nlev, opu, maxrank=1, flippable=fl)
        call augcls_freeze_site(opr, 1, cen); call augcls_freeze_site(opu, 1, cen)
        print '(a,i0,a,i0,a,i0)', '    central species ', cen, ': nodes unreduced ', opu%nn, '  reduced ', opr%nn
        call edn_work_init(w1, opu%dim()); call edn_work_init(w2, opr%dim())
        do io = 1, 4
          call augcls_set_start(opr, l0s(io), d)
          do kz = 1, 2
            z = cmplx(-0.2_dp + 0.3_dp * kz, 0.15_dp, dp)
            call edc_green(opu, w1, ed, z, l0s(io), nlev, g1, a1, b1, m1)
            call edc_green(opr, w2, ed, z, l0s(io), nlev, g2, a2, b2, m2)
            err = max(err, maxval(abs(a1(1:m1) - a2(1:m2))), maxval(abs(b1(1:m1 - 1) - b2(1:m2 - 1))))
            errg = max(errg, abs(g1 - g2))
          end do
        end do
      end do
    end do
    print '(a,es9.2,a,es9.2)', '  decorated bcc: max |coefficient difference| ', err, '   max |G difference| ', errg
    call check_close(err, 0.0_dp, 1.0e-9_dp, 'decorated bcc: class-reduced == unreduced chain coefficients')
    call check_close(errg, 0.0_dp, 1.0e-9_dp, 'decorated bcc: class-reduced == unreduced Green function')
  end subroutine decorated_bcc

  subroutine regressions()
    integer, parameter :: nlev = 3
    type(graph_t) :: g
    type(orb_alloy_t) :: m0
    type(class_alloy_t) :: mb, mt
    type(augcls_op) :: opc
    type(augsym_op) :: sop
    type(asr_ed_t) :: edb
    type(edn_t) :: ed
    type(edn_work_t) :: w1
    type(ed_work_t) :: w2
    type(multi_alloy_t) :: mn
    type(augn_op) :: opn
    integer, allocatable :: coord(:, :), cls(:), perm(:, :), mats(:, :, :), p16(:, :), m16(:, :, :)
    real(dp), allocatable :: d(:, :, :)
    logical, allocatable :: fl(:)
    complex(dp) :: z, g1, g2, a1(nlev), b1(nlev), a2(nlev), b2(nlev)
    real(dp) :: errb, errt
    integer :: cen, io, kz, m1, m2
    type(potpar_t) :: p3(3)
    call build_lmto_alloy('fcc', 6, 1.2_dp, 0.7_dp, pp(1), pp(2), 0.4_dp, g, coord, m0, alpha_in=alpha)
    allocate (cls(g%n), fl(g%n)); cls = 1; fl = .true.; fl(1) = .false.
    call class_ops(coord, 12, cls, perm, mats)
    call d4h_ops(coord, 12, p16, m16)
    call check(all(perm == p16), 'one class: class_ops == d4h_ops')
    call orbital_rep(mats, d)
    ! binary: A = species 1 (x = 0.4), B = species 2
    mb%nl = 9; mb%nsp = 2; mb%ncls = 1
    allocate (mb%nsc(1), mb%spc(2, 1), mb%conc(2, 1)); mb%nsc = 2; mb%spc(:, 1) = [1, 2]; mb%conc(:, 1) = [0.4_dp, 0.6_dp]
    mb%block = m0%block; mb%block0 = m0%block0
    errb = 0.0_dp
    do cen = 1, 2
      call edc_setup(mb, pp(1:2), alpha, cen, ed)
      call asr_ed_setup(m0, pp(1), pp(2), alpha, cen == 1, edb)
      call build_augcls_operator(g, mb, cls, perm, d, 1, nlev, opc, maxrank=1, flippable=fl)
      call build_augsym_operator(g, m0, perm, d, 1, nlev, sop, maxrank=1, flippable=fl)
      call augcls_freeze_site(opc, 1, cen); call augsym_freeze_site(sop, 1, edb%mc)
      print '(a,i0,a,i0,a,i0)', '  binary regression, central ', cen, ': nodes augcls ', opc%nn, '  augsym ', sop%nn
      call check(opc%nn == sop%nn, 'binary: same node count as asr_augsym')
      call edn_work_init(w1, opc%dim()); call ed_work_init(w2, sop%dim())
      do io = 1, 4
        call augcls_set_start(opc, l0s(io), d); call augsym_set_start(sop, l0s(io), d)
        do kz = 1, 2
          z = cmplx(-0.2_dp + 0.3_dp * kz, 0.15_dp, dp)
          call edc_green(opc, w1, ed, z, l0s(io), nlev, g1, a1, b1, m1)
          call asr_ed_green(sop, w2, edb, z, l0s(io), nlev, g2, a2, b2, m2)
          errb = max(errb, maxval(abs(a1(1:m1) - a2(1:m2))), maxval(abs(b1(1:m1 - 1) - b2(1:m2 - 1))), abs(g1 - g2))
        end do
      end do
    end do
    print '(a,es9.2)', '  binary regression: max |augcls - augsym| ', errb
    call check_close(errb, 0.0_dp, 1.0e-9_dp, 'one binary class == asr_augsym / asr_ed')
    ! ternary: species 2, 3, 4 with (0.3, 0.5, 0.2), reduced augcls vs unreduced augn
    mt%nl = 9; mt%nsp = 3; mt%ncls = 1
    allocate (mt%nsc(1), mt%spc(3, 1), mt%conc(3, 1)); mt%nsc = 3; mt%spc(:, 1) = [1, 2, 3]; mt%conc(:, 1) = [0.3_dp, 0.5_dp, 0.2_dp]
    mt%block = m0%block; mt%block0 = m0%block0
    mn%nl = 9; mn%ns = 3; mn%conc = mt%conc(:, 1); mn%block = m0%block; mn%block0 = m0%block0
    p3 = pp(2:4)
    errt = 0.0_dp
    do cen = 1, 3
      call edc_setup(mt, p3, alpha, cen, ed)
      call edn_setup(mn, p3, alpha, cen, ed)
      call build_augcls_operator(g, mt, cls, perm, d, 1, nlev, opc, maxrank=1, flippable=fl)
      call build_augn_operator(g, mn, 1, nlev, opn, maxrank=1, flippable=fl)
      call augcls_freeze_site(opc, 1, cen); call augn_freeze_site(opn, 1, cen)
      print '(a,i0,a,i0,a,i0)', '  ternary regression, central ', cen, ': nodes augcls ', opc%nn, '  augn (unreduced) ', opn%nn
      call edn_work_init(w1, opc%dim())
      block
        type(edn_work_t) :: w3
        call edn_work_init(w3, opn%dim())
        do io = 1, 4
          call augcls_set_start(opc, l0s(io), d)
          do kz = 1, 2
            z = cmplx(-0.2_dp + 0.3_dp * kz, 0.15_dp, dp)
            call edc_green(opc, w1, ed, z, l0s(io), nlev, g1, a1, b1, m1)
            call edn_green(opn, w3, ed, z, l0s(io), nlev, g2, a2, b2, m2)
            errt = max(errt, maxval(abs(a1(1:m1) - a2(1:m2))), maxval(abs(b1(1:m1 - 1) - b2(1:m2 - 1))), abs(g1 - g2))
          end do
        end do
      end block
    end do
    print '(a,es9.2)', '  ternary regression: max |augcls - augn| ', errt
    call check_close(errt, 0.0_dp, 1.0e-9_dp, 'one ternary class (reduced) == unreduced asr_augn / asr_edn')
  end subroutine regressions

  !> ring of 5: site 1 (class 3, ternary {1,3,4}) frozen; its neighbours are shell classes 4 (copy of 3, on site 2) and 5 (copy of
  !> class 2 = binary {2,3}, on site 5); site 3 ordered, site 4 binary uncorrelated.  Conditional concentrations per central
  !> species, some zero.
  subroutine ring_sro()
    integer, parameter :: n = 5, nlev = 400, rcls(n) = [3, 4, 1, 2, 5]
    real(dp), parameter :: c4(3, 3) = reshape([0.0_dp, 0.7_dp, 0.3_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.2_dp, 0.6_dp, 0.2_dp], [3, 3])
    real(dp), parameter :: c5(2, 3) = reshape([0.1_dp, 0.9_dp, 1.0_dp, 0.0_dp, 0.55_dp, 0.45_dp], [2, 3])
    type(graph_t) :: g
    type(class_alloy_t) :: mdl, mb
    type(augcls_op) :: op
    type(edn_t) :: ed
    type(edn_work_t) :: wk
    integer :: perm(n, 1), i, l2, j, cen, kz, m, cf(n), q, c, ncfg, l0
    real(dp) :: dmat(9, 9, 1), u, pc, err
    complex(dp) :: z, gre, gden, a(nlev), b(nlev)
    logical :: fl(n)
    allocate (g%nbr_ptr(n + 1), g%nbr(2 * n)); g%n = n
    do i = 1, n
      g%nbr_ptr(i) = 2 * (i - 1) + 1
      g%nbr(2 * i - 1) = mod(i, n) + 1
      g%nbr(2 * i) = mod(i - 2 + n, n) + 1
    end do
    g%nbr_ptr(n + 1) = 2 * n + 1
    call three_classes(mb)
    mdl%nl = 9; mdl%nsp = 4; mdl%ncls = 5
    allocate (mdl%nsc(5), mdl%spc(3, 5), mdl%conc(3, 5))
    mdl%nsc(1:3) = mb%nsc; mdl%spc(:, 1:3) = mb%spc; mdl%conc(:, 1:3) = mb%conc
    mdl%nsc(4) = mb%nsc(3); mdl%spc(:, 4) = mb%spc(:, 3); mdl%conc(:, 4) = mb%conc(:, 3)
    mdl%nsc(5) = mb%nsc(2); mdl%spc(:, 5) = mb%spc(:, 2); mdl%conc(:, 5) = mb%conc(:, 2)
    allocate (mdl%block(9, 9, 2), mdl%block0(9, 9))
    do l = 1, 9; do l2 = 1, 9
      call random_number(u); mdl%block(l, l2, 1) = 0.5_dp * (u - 0.5_dp)
    end do; end do
    mdl%block(:, :, 2) = transpose(mdl%block(:, :, 1))
    do l = 1, 9; do l2 = l, 9
      call random_number(u); mdl%block0(l, l2) = 0.6_dp * (u - 0.5_dp); mdl%block0(l2, l) = mdl%block0(l, l2)
    end do; end do
    perm(:, 1) = [(i, i=1, n)]
    dmat = 0.0_dp; do l = 1, 9; dmat(l, l, 1) = 1.0_dp; end do
    fl = .true.; fl(1) = .false.
    err = 0.0_dp
    do j = 1, 3
      cen = mdl%spc(j, 3)
      call edc_setup(mdl, pp, alpha, cen, ed)
      call build_augcls_operator(g, mdl, rcls, perm, dmat, 1, 80, op, flippable=fl)
      call augcls_freeze_site(op, 1, cen)
      call augcls_set_class_conc(op, 4, c4(:, j)); call augcls_set_class_conc(op, 5, c5(:, j))
      call edn_work_init(wk, op%dim())
      do l0 = 1, 9, 4
        do kz = 1, 3
          z = cmplx(-0.3_dp + 0.3_dp * kz, 0.3_dp, dp)
          call edc_green(op, wk, ed, z, l0, nlev, gre, a, b, m)
          gre = cf_complex(a, b, m, z)
          gden = (0.0_dp, 0.0_dp)
          ncfg = product(mdl%nsc(rcls(2:n)))
          do c = 0, ncfg - 1
            q = c; pc = 1.0_dp; cf(1) = cen
            do i = 2, n
              associate (k => rcls(i))
                cf(i) = mdl%spc(mod(q, mdl%nsc(k)) + 1, k)
                select case (k)
                case (4); pc = pc * c4(mod(q, mdl%nsc(k)) + 1, j)
                case (5); pc = pc * c5(mod(q, mdl%nsc(k)) + 1, j)
                case default; pc = pc * mdl%conc(mod(q, mdl%nsc(k)) + 1, k)
                end select
                q = q / mdl%nsc(k)
              end associate
            end do
            if (pc > 0.0_dp) gden = gden + pc * dense_g(mdl, g, cf, cen, z, l0)
          end do
          err = max(err, abs(gre - gden))
        end do
      end do
    end do
    print '(a,es10.3)', '  ring SRO: max |recursion - exact conditional average| ', err
    call check_close(err, 0.0_dp, 1.0e-6_dp, 'SRO shell classes (zero weights included): ED recursion == exact conditional average')
  end subroutine ring_sro

  !> fcc graph, one ternary class {2,3,4}; class 2 = first shell of the origin (geometric, minimum image).
  subroutine fcc_sro()
    integer, parameter :: nlev = 3
    real(dp), parameter :: cs(3, 3) = reshape([0.0_dp, 0.8_dp, 0.2_dp, 0.4_dp, 0.4_dp, 0.2_dp, 0.3_dp, 0.3_dp, 0.4_dp], [3, 3])
    type(graph_t) :: g
    type(orb_alloy_t) :: m0
    type(class_alloy_t) :: mdl
    type(augcls_op) :: opr, opu, op0
    type(edn_t) :: ed
    type(edn_work_t) :: w1, w2, w0
    integer, allocatable :: coord(:, :), cls(:), cls0(:), perm(:, :), perm0(:, :), mats(:, :, :), pid(:, :)
    real(dp), allocatable :: d(:, :, :), did(:, :, :)
    logical, allocatable :: fl(:)
    complex(dp) :: z, g1, g2, g0, a1(nlev), b1(nlev), a2(nlev), b2(nlev)
    real(dp) :: err, err0
    integer :: i, j, cen, io, kz, m1, m2, r2, r2min, nsh
    call build_lmto_alloy('fcc', 4, 1.05_dp, 0.7_dp, pp(1), pp(2), 0.5_dp, g, coord, m0, alpha_in=alpha)
    mdl%nl = 9; mdl%nsp = 4; mdl%ncls = 2
    allocate (mdl%nsc(2), mdl%spc(3, 2), mdl%conc(3, 2))
    mdl%nsc = 3; mdl%spc(:, 1) = [2, 3, 4]; mdl%spc(:, 2) = [2, 3, 4]
    mdl%conc(:, 1) = [0.3_dp, 0.5_dp, 0.2_dp]; mdl%conc(:, 2) = mdl%conc(:, 1)
    mdl%block = m0%block; mdl%block0 = m0%block0
    allocate (fl(g%n)); fl = .true.; fl(1) = .false.
    allocate (cls0(g%n)); cls0 = 1; cls = cls0
    r2min = huge(1)
    do i = 2, g%n; r2min = min(r2min, sum((modulo(coord(:, i) + 4, 8) - 4)**2)); end do
    nsh = 0
    do i = 2, g%n
      r2 = sum((modulo(coord(:, i) + 4, 8) - 4)**2)
      if (r2 == r2min) then; cls(i) = 2; nsh = nsh + 1; end if
    end do
    call check(nsh == 12, 'fcc: 12 first-shell sites')
    call class_ops(coord, 8, cls0, perm0, mats)
    call class_ops(coord, 8, cls, perm, mats)
    call orbital_rep(mats, d)
    call check(size(perm, 2) == size(perm0, 2), 'fcc: the shell class keeps the site group')
    allocate (pid(g%n, 1), did(9, 9, 1)); pid(:, 1) = [(i, i=1, g%n)]; did = 0.0_dp
    do i = 1, 9; did(i, i, 1) = 1.0_dp; end do
    err = 0.0_dp; err0 = 0.0_dp
    do j = 1, 3
      cen = mdl%spc(j, 1)
      call edc_setup(mdl, pp, alpha, cen, ed)
      call build_augcls_operator(g, mdl, cls, perm, d, 1, nlev, opr, maxrank=1, flippable=fl)
      call build_augcls_operator(g, mdl, cls, pid, did, 1, nlev, opu, maxrank=1, flippable=fl)
      call build_augcls_operator(g, mdl, cls0, perm0, d, 1, nlev, op0, maxrank=1, flippable=fl)
      call augcls_freeze_site(opr, 1, cen); call augcls_freeze_site(opu, 1, cen); call augcls_freeze_site(op0, 1, cen)
      if (j == 1) print '(a,i0,a,i0,a,i0)', '  fcc SRO: nodes unreduced ', opu%nn, '  reduced ', opr%nn, '  no shell class ', op0%nn
      call edn_work_init(w1, opu%dim()); call edn_work_init(w2, opr%dim()); call edn_work_init(w0, op0%dim())
      do io = 1, 4
        call augcls_set_start(opr, l0s(io), d); call augcls_set_start(op0, l0s(io), d)
        do kz = 1, 2
          z = cmplx(-0.2_dp + 0.3_dp * kz, 0.15_dp, dp)
          call augcls_set_class_conc(opr, 2, mdl%conc(:, 1))            ! shell = x: no SRO
          call edc_green(opr, w2, ed, z, l0s(io), nlev, g2, a2, b2, m2)
          call edc_green(op0, w0, ed, z, l0s(io), nlev, g0)
          err0 = max(err0, abs(g2 - g0))
          call augcls_set_class_conc(opr, 2, cs(:, j)); call augcls_set_class_conc(opu, 2, cs(:, j))
          call edc_green(opu, w1, ed, z, l0s(io), nlev, g1, a1, b1, m1)
          call edc_green(opr, w2, ed, z, l0s(io), nlev, g2, a2, b2, m2)
          err = max(err, maxval(abs(a1(1:m1) - a2(1:m2))), maxval(abs(b1(1:m1 - 1) - b2(1:m2 - 1))), abs(g1 - g2))
        end do
      end do
    end do
    print '(a,es9.2,a,es9.2)', '  fcc SRO: max |reduced - unreduced| ', err, '   max |G(shell = x) - G(no shell class)| ', err0
    call check_close(err, 0.0_dp, 1.0e-9_dp, 'fcc SRO: class-reduced == unreduced with conditional shell concentrations')
    call check_close(err0, 0.0_dp, 1.0e-12_dp, 'fcc SRO: shell concentrations = x reproduce the operator without a shell class')
  end subroutine fcc_sro
end program test_augcls
