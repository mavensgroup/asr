!> n-species energy-dependent ASR (asr_edn on asr_augn) vs brute force on a ring of N sites, 9 orbitals, THREE species with
!> concentrations (0.3, 0.5, 0.2): for every central species c the start site is frozen to c and the others are averaged:
!>   (1) recursion == exact configuration average of the dense [E - H(E)]^-1_00, H(E) = E - T (E - E_nu - h) T - T h o h T;
!>   (2) that average == the same average of [E - H^(2)]^-1_00 (E - H(E) = T (E - H^(2)) T, T_0 = 1 at the central site).
program test_asr_ed_n
  use asr_kinds, only: dp
  use asr_augspace, only: graph_t
  use asr_augn
  use asr_lmto_ham, only: potpar_t, to_alpha, o_alpha
  use asr_edn
  use asr_testing
  implicit none
  integer, parameter :: n = 4, nl = 9, nd = n * nl, nlev = 500, ns = 3
  integer, parameter :: orb_l(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]
  type(graph_t) :: g
  type(multi_alloy_t) :: mdl
  type(augn_op) :: op
  type(potpar_t) :: pp(ns), qq(ns)
  type(edn_t) :: ed
  type(edn_work_t) :: wk
  real(dp) :: alpha(0:2), conc(ns), oo(0:2, ns), pc, u
  complex(dp) :: z, gre, gdense, gh2, a(nlev), b(nlev)
  real(dp) :: err1, err2, err3
  integer :: i, l, l2, c, m, l0, kz, cen, sp, cf(n), q
  logical :: fl(n)
  external :: zgesv

  print '(a)', 'test_asr_ed_n'
  conc = [0.3_dp, 0.5_dp, 0.2_dp]
  do l = 0, 2
    pp(1)%enu(l) = 0.1_dp * l; pp(1)%c(l) = 0.3_dp - 0.15_dp * l; pp(1)%sqrdel(l) = 0.30_dp - 0.07_dp * l
    pp(1)%p(l) = 0.9_dp - 0.1_dp * l; pp(1)%gamma(l) = 0.25_dp - 0.1_dp * l
    pp(2)%enu(l) = 0.05_dp * l; pp(2)%c(l) = 0.45_dp - 0.1_dp * l; pp(2)%sqrdel(l) = 0.27_dp - 0.06_dp * l
    pp(2)%p(l) = 0.8_dp - 0.1_dp * l; pp(2)%gamma(l) = 0.2_dp - 0.08_dp * l
    pp(3)%enu(l) = -0.07_dp + 0.08_dp * l; pp(3)%c(l) = 0.2_dp - 0.05_dp * l; pp(3)%sqrdel(l) = 0.22_dp - 0.04_dp * l
    pp(3)%p(l) = 0.7_dp - 0.1_dp * l; pp(3)%gamma(l) = 0.15_dp - 0.06_dp * l
  end do
  do sp = 1, ns; pp(sp)%alpha_rep = pp(sp)%gamma; end do
  alpha = [0.3485_dp, 0.05303_dp, 0.010714_dp]
  do sp = 1, ns
    qq(sp) = to_alpha(pp(sp), alpha); oo(:, sp) = o_alpha(qq(sp), alpha)
  end do
  allocate (g%nbr_ptr(n + 1), g%nbr(2 * n)); g%n = n
  do i = 1, n
    g%nbr_ptr(i) = 2 * (i - 1) + 1
    g%nbr(2 * i - 1) = mod(i, n) + 1
    g%nbr(2 * i) = mod(i - 2 + n, n) + 1
  end do
  g%nbr_ptr(n + 1) = 2 * n + 1
  mdl%nl = nl; mdl%ns = ns
  allocate (mdl%conc(ns)); mdl%conc = conc
  allocate (mdl%block(nl, nl, 2), mdl%block0(nl, nl))
  call random_seed()
  do l = 1, nl; do l2 = 1, nl
      call random_number(u); mdl%block(l, l2, 1) = 0.5_dp * (u - 0.5_dp)
  end do; end do
  mdl%block(:, :, 2) = transpose(mdl%block(:, :, 1))
  do l = 1, nl; do l2 = l, nl
      call random_number(u); mdl%block0(l, l2) = 0.6_dp * (u - 0.5_dp); mdl%block0(l2, l) = mdl%block0(l, l2)
  end do; end do

  err1 = 0.0_dp; err2 = 0.0_dp; err3 = 0.0_dp
  do cen = 1, ns
    call edn_setup(mdl, pp, alpha, cen, ed)
    fl = .true.; fl(1) = .false.
    call build_augn_operator(g, mdl, 1, 80, op, flippable=fl)
    call augn_freeze_site(op, 1, cen)
    print '(a,i0,a,i0,a)', '  central species ', cen, ': ', op%nn, ' augmented nodes'
    call check(op%nn == n * ns**(n - 1), 'node set = N ns^(N-1) (start site frozen)')
    call edn_work_init(wk, op%dim())
    do l0 = 1, 9, 4
      do kz = 1, 3
        z = cmplx(-0.3_dp + 0.3_dp * kz, 0.3_dp, dp)
        call edn_green(op, wk, ed, z, l0, nlev, gre, a, b, m)
        gre = cf_complex(a, b, m, z)
        gdense = (0.0_dp, 0.0_dp); gh2 = (0.0_dp, 0.0_dp)
        do c = 0, ns**(n - 1) - 1
          q = c; pc = 1.0_dp; cf(1) = cen
          do i = 2, n
            cf(i) = mod(q, ns) + 1; q = q / ns; pc = pc * conc(cf(i))
          end do
          call dense_green(cf, cen, z, l0, gdense, gh2, pc)
        end do
        err1 = max(err1, abs(gre - gdense)); err2 = max(err2, abs(gdense - gh2)); err3 = max(err3, abs(gre - gh2))
      end do
    end do
  end do
  print '(a,es10.3,a,es10.3,a,es10.3)', '  |rec - dense H(E)| ', err1, '   |dense H(E) - dense H2| ', err2, '   |rec - dense H2| ', err3
  call check_close(err2, 0.0_dp, 1.0e-9_dp, 'E - H(E) = T (E - H2) T: dense average == H2 average (three species)')
  call check_close(err1, 0.0_dp, 1.0e-6_dp, 'three-species ED recursion == exact configuration average of [E - H(E)]^-1')
  call report()
contains
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

  subroutine dense_green(cfg, ce, zz, l0_, gh, gh2_, w)
    integer, intent(in) :: cfg(n), ce, l0_
    complex(dp), intent(in) :: zz
    complex(dp), intent(inout) :: gh, gh2_
    real(dp), intent(in) :: w
    real(dp) :: dl(nd), cm(nd), en(nd), om(nd), tt(nd), hm(nd, nd), s(nd, nd), dc(nd), hoh(nd, nd), h2(nd, nd)
    complex(dp) :: ae(nd, nd), rhs(nd), a2(nd, nd), r2(nd)
    integer :: ii, jj, kk, ll, ll2, ipiv(nd), info, spc, kq
    s = 0.0_dp
    do ii = 1, n
      do ll = 1, nl; do ll2 = 1, nl
          s(nl * (ii - 1) + ll, nl * (ii - 1) + ll2) = mdl%block0(ll, ll2)
      end do; end do
      do kk = 1, 2
        jj = g%nbr(g%nbr_ptr(ii) + kk - 1)
        do ll = 1, nl; do ll2 = 1, nl
            s(nl * (ii - 1) + ll, nl * (jj - 1) + ll2) = mdl%block(ll, ll2, kk)
        end do; end do
      end do
    end do
    do ii = 1, n
      spc = cfg(ii)
      do ll = 1, nl
        kq = nl * (ii - 1) + ll
        dl(kq) = qq(spc)%sqrdel(orb_l(ll))**2
        cm(kq) = qq(spc)%c(orb_l(ll)) - qq(spc)%enu(orb_l(ll))
        en(kq) = qq(spc)%enu(orb_l(ll)); om(kq) = oo(orb_l(ll), spc)
        dc(kq) = qq(ce)%sqrdel(orb_l(ll))**2
      end do
    end do
    hm = 0.0_dp
    do ii = 1, nd
      do jj = 1, nd; hm(ii, jj) = sqrt(dl(ii)) * s(ii, jj) * sqrt(dl(jj)); end do
      hm(ii, ii) = hm(ii, ii) + cm(ii)
    end do
    hoh = 0.0_dp
    do ii = 1, nd; do jj = 1, nd
        hoh(ii, jj) = sum(hm(ii, :) * om * hm(:, jj))
    end do; end do
    h2 = hm - hoh
    do ii = 1, nd; h2(ii, ii) = h2(ii, ii) + en(ii); end do
    tt = sqrt(dc / dl)
    do ii = 1, nd; do jj = 1, nd                     ! E - H(E) = T (E - E_nu - h) T + T hoh T
        ae(ii, jj) = tt(ii) * tt(jj) * (merge(zz - en(ii), (0.0_dp, 0.0_dp), ii == jj) - hm(ii, jj)) + tt(ii) * tt(jj) * hoh(ii, jj)
    end do; end do
    rhs = (0.0_dp, 0.0_dp); rhs(l0_) = (1.0_dp, 0.0_dp)
    call zgesv(nd, 1, ae, nd, ipiv, rhs, nd, info)
    gh = gh + w * rhs(l0_)
    a2 = -h2; do ii = 1, nd; a2(ii, ii) = a2(ii, ii) + zz; end do
    r2 = (0.0_dp, 0.0_dp); r2(l0_) = (1.0_dp, 0.0_dp)
    call zgesv(nd, 1, a2, nd, ipiv, r2, nd, info)
    gh2_ = gh2_ + w * r2(l0_)
  end subroutine dense_green
end program test_asr_ed_n
