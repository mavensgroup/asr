!> AMscf-form energy-dependent ASR (asr_ed) vs brute force on a ring of N sites (9 orbitals per site), start site frozen to
!> the central species c, the other sites A (concentration x) or B.  Krylov space (nearly) exhausted => the recursion must equal
!>   (1) the exact configuration average of the dense  [E - H(E)]^-1_{00},  H(E) = E - T (E - E_nu - h) T - T h o h T,
!>   (2) the same average of [E - H^(2)]^-1_{00} (E - H(E) = T (E - H^(2)) T and T_0 = 1 at the central site),
!> for both central species, with the nonrandom S' = Delta_c^1/2 S Delta_c^1/2 of the AMscf form.
program test_asr_ed
  use asr_kinds, only: dp
  use asr_augspace, only: graph_t
  use asr_augorb
  use asr_lmto_ham, only: potpar_t, to_alpha, o_alpha
  use asr_endep, only: ed_work_t, ed_work_init
  use asr_ed
  use asr_testing
  implicit none
  integer, parameter :: n = 4, nl = 9, nd = n * nl, nlev = 200
  integer, parameter :: orb_l(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]
  type(graph_t) :: g
  type(orb_alloy_t) :: mdl
  type(augorb_op) :: op
  type(potpar_t) :: pa, pb, qq(2)
  type(asr_ed_t) :: ed
  type(ed_work_t) :: wk
  real(dp) :: alpha(0:2), x, pc, oo(0:2, 2), blk1(nl, nl)
  complex(dp) :: z, gre, gdense, gh2, a(nlev), b(nlev)
  real(dp) :: err1, err2, err3, u
  integer :: i, k, l, l2, c, m, l0, kz, cen
  logical :: isa(n), ctr_a
  external :: zgesv

  print '(a)', 'test_asr_ed'
  do l = 0, 2
    pa%enu(l) = 0.1_dp * l; pa%c(l) = 0.3_dp - 0.15_dp * l; pa%sqrdel(l) = 0.30_dp - 0.07_dp * l
    pa%p(l) = 0.9_dp - 0.1_dp * l; pa%gamma(l) = 0.25_dp - 0.1_dp * l; pa%alpha_rep(l) = pa%gamma(l)
    pb%enu(l) = 0.05_dp * l; pb%c(l) = 0.45_dp - 0.1_dp * l; pb%sqrdel(l) = 0.27_dp - 0.06_dp * l
    pb%p(l) = 0.8_dp - 0.1_dp * l; pb%gamma(l) = 0.2_dp - 0.08_dp * l; pb%alpha_rep(l) = pb%gamma(l)
  end do
  alpha = [0.3485_dp, 0.05303_dp, 0.010714_dp]
  qq(1) = to_alpha(pa, alpha); qq(2) = to_alpha(pb, alpha)
  oo(:, 1) = o_alpha(qq(1), alpha); oo(:, 2) = o_alpha(qq(2), alpha)
  x = 0.35_dp
  ! ring: slot 1 -> R+1, slot 2 -> R-1
  allocate (g%nbr_ptr(n + 1), g%nbr(2 * n)); g%n = n
  do i = 1, n
    g%nbr_ptr(i) = 2 * (i - 1) + 1
    g%nbr(2 * i - 1) = mod(i, n) + 1
    g%nbr(2 * i) = mod(i - 2 + n, n) + 1
  end do
  g%nbr_ptr(n + 1) = 2 * n + 1
  mdl%nl = nl; mdl%x = x
  call random_seed()
  do l = 1, nl; do l2 = 1, nl
      call random_number(u); blk1(l, l2) = 0.5_dp * (u - 0.5_dp)
  end do; end do
  allocate (mdl%block(nl, nl, 2), mdl%block0(nl, nl))
  mdl%block(:, :, 1) = blk1; mdl%block(:, :, 2) = transpose(blk1)
  do l = 1, nl; do l2 = l, nl
      call random_number(u); mdl%block0(l, l2) = 0.6_dp * (u - 0.5_dp); mdl%block0(l2, l) = mdl%block0(l, l2)
  end do; end do

  err1 = 0.0_dp; err2 = 0.0_dp; err3 = 0.0_dp
  do cen = 1, 2
    ctr_a = (cen == 1)
    call asr_ed_setup(mdl, pa, pb, alpha, ctr_a, ed)
    block
      logical :: flippable(n)
      flippable = .true.; flippable(1) = .false.
      call build_augorb_operator(g, mdl, 1, 60, op, flippable=flippable)
    end block
    call augorb_freeze_site(op, 1, ed%mc)
    print '(a,i0,a,i0,a)', '  central species ', cen, ': ', op%nn, ' augmented nodes'
    call check(op%nn == n * 2**(n - 1), 'node set = N 2^(N-1) (start site frozen)')
    call ed_work_init(wk, op%dim())
    do l0 = 1, 9, 4
      do kz = 1, 4
        z = cmplx(-0.3_dp + 0.25_dp * kz, 0.3_dp, dp)
        call asr_ed_green(op, wk, ed, z, l0, nlev, gre, a, b, m)
        gre = cf_complex(a, b, m, z)
        ! brute force
        gdense = (0.0_dp, 0.0_dp); gh2 = (0.0_dp, 0.0_dp)
        do c = 0, 2**(n - 1) - 1
          pc = 1.0_dp
          isa(1) = ctr_a
          do i = 2, n
            isa(i) = btest(c, i - 2)
            if (isa(i)) then; pc = pc * x; else; pc = pc * (1.0_dp - x); end if
          end do
          call dense_green(isa, ctr_a, z, l0, gdense, gh2, pc)
        end do
        err1 = max(err1, abs(gre - gdense))
        err2 = max(err2, abs(gdense - gh2))
        err3 = max(err3, abs(gre - gh2))
      end do
    end do
  end do
  print '(a,es10.3,a,es10.3,a,es10.3)', '  |rec - dense H(E)| ', err1, '   |dense H(E) - dense H2| ', err2, '   |rec - dense H2| ', err3
  call check_close(err2, 0.0_dp, 1.0e-9_dp, 'E - H(E) = T (E - H2) T: dense H(E) average == dense H2 average (central T = 1)')
  call check_close(err1, 0.0_dp, 1.0e-6_dp, 'AMscf-form augmented-space recursion == exact configuration average of [E - H(E)]^-1')
  call report()
contains
  !> continued fraction of the complex-symmetric chain (a, b), no tail: G = 1/(z - a1 - b1^2/(z - a2 - ...))
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

  subroutine dense_green(isa_, ca, zz, l0_, gh, gh2_, w)
    ! accumulates w * [E - H(E)]^-1_{l0,l0} into gh and w * [E - H2]^-1_{l0,l0} into gh2_
    logical, intent(in) :: isa_(n), ca
    complex(dp), intent(in) :: zz
    integer, intent(in) :: l0_
    complex(dp), intent(inout) :: gh, gh2_
    real(dp), intent(in) :: w
    real(dp) :: dl(nd), cm(nd), en(nd), om(nd), tt(nd), hm(nd, nd), s(nd, nd), dc(nd), hoh(nd, nd), h2(nd, nd)
    complex(dp) :: ae(nd, nd), rhs(nd), a2(nd, nd), r2(nd)
    integer :: ii, jj, kk, ll, ll2, ipiv(nd), info, spc
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
      spc = merge(1, 2, isa_(ii))
      do ll = 1, nl
        k = nl * (ii - 1) + ll
        dl(k) = qq(spc)%sqrdel(orb_l(ll))**2                    ! Delta_R
        cm(k) = qq(spc)%c(orb_l(ll)) - qq(spc)%enu(orb_l(ll))
        en(k) = qq(spc)%enu(orb_l(ll)); om(k) = oo(orb_l(ll), spc)
        dc(k) = qq(merge(1, 2, ca))%sqrdel(orb_l(ll))**2          ! Delta_c
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
    ! H(E) = E - T (E - E_nu - h) T - T hoh T
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
end program test_asr_ed
