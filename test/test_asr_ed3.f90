!> Gate 1 of docs/DERIVATION_ed3.md: the full-ASA energy-dependent augmented-space operator against brute force on a ring of
!> N sites (9 orbitals), binary alloy, start site frozen to the central species c, all other sites A (x) or B.
!> Brute force: K = E O - H of bands_core (src/asr_asa.f90), O = 1 + o h + h o + h p h,
!> H = E_nu + h + E_nu o h + h o E_nu + h (o + p E_nu) h, h = c + d S d (alpha representation, p from to_alpha).
!> Augmented space: K' = T K T = et - h' + u h' + h' u - h' Qt h', h' = g + S' (src operator asr_augorb with o = E_nu = 0
!> gives h'; mapply gives the configuration operator M), et = Dc/Ds eps, u = o eps, Qt = Ds/Dc (o - p eps), eps = E - E_nu.
!> Checks (every central species, s/p/d start orbitals, four complex energies):
!>   1. K'_aug is complex symmetric.
!>   2. Dense solve of K'_aug: <u|K'^-1|u> = <G11>, g_c G_uu + <w|K'^-1|u> = <(hG)11>,
!>      g_c^2 G_uu + 2 g_c G_wu + <w|K'^-1|w> = <(hGh)11>, with u = |1 L, no flip> and w = S' u (configuration independent).
!>   3. Two Lanczos chains (from u, projected on w; from w, projected on u) == the dense solve; G_uw == G_wu.
!>   4. Regression: u = 0, p = 0, Qt -> -q gives the asr_edc operator: <u|K'^-1|u> = <[E - H^(2)]^-1_11>.
!>   5. Canonical form D - (V S' + S' V) - S' Qt S', V = (1/2 + c o) - (o + c p) eps, == the factored product.
!>   6. Single-site dressing (h z)_1 = eps_1 z_1 / (1 - eps_1 o_1): exact at p^gamma = 0 (checked); its error at p^gamma /= 0
!>      is printed (measurement, not a gate).
program test_asr_ed3
  use asr_kinds, only: dp
  use asr_augspace, only: graph_t
  use asr_augorb
  use asr_lmto_ham, only: potpar_t, to_alpha, o_alpha
  use asr_testing
  implicit none
  integer, parameter :: n = 4, nl = 9, nd = n * nl
  integer, parameter :: orb_l(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]
  type(graph_t) :: g
  type(orb_alloy_t) :: mdl
  type(augorb_op) :: op
  type(potpar_t) :: pa, pb, qq(2)
  real(dp) :: alpha(0:2), x, oo(0:2, 2), blk1(nl, nl), u, eta
  real(dp) :: dsp(9, 2), csp(9, 2), esp(9, 2), osp(9, 2), psp(9, 2), dcen(9), gsp(9, 2), qsp(9, 2)
  complex(dp) :: et(9, 2), uu(9, 2), qt(9, 2), zero9(9, 2), z
  complex(dp), allocatable :: amat(:, :), afac(:, :), rhs(:, :), uvec(:), wvec(:), ej(:), cvec(:), rv(:)
  complex(dp) :: guu, gwu, gww, guw, faug, maug, b11, bh, bhh, b2, cuu, cwu, cww, cuw
  real(dp) :: esym, eg, ef, em, ech, esw, ereg, ecan, ewu, tmp, eq0, eqmax, eqrel, peps
  integer :: i, l, l2, cen, l0, kz, ndim, j, info, il
  integer, allocatable :: ipiv(:)
  integer :: lstart(3) = [1, 2, 5]
  external :: zgesv

  print '(a)', 'test_asr_ed3'
  do l = 0, 2
    pa%enu(l) = 0.1_dp * l; pa%c(l) = 0.3_dp - 0.15_dp * l; pa%sqrdel(l) = 0.30_dp - 0.07_dp * l
    pa%p(l) = 0.9_dp - 0.1_dp * l; pa%gamma(l) = 0.25_dp - 0.1_dp * l; pa%alpha_rep(l) = pa%gamma(l)
    pb%enu(l) = 0.05_dp * l; pb%c(l) = 0.45_dp - 0.1_dp * l; pb%sqrdel(l) = 0.27_dp - 0.06_dp * l
    pb%p(l) = 0.8_dp - 0.1_dp * l; pb%gamma(l) = 0.2_dp - 0.08_dp * l; pb%alpha_rep(l) = pb%gamma(l)
  end do
  alpha = [0.3485_dp, 0.05303_dp, 0.010714_dp]
  qq(1) = to_alpha(pa, alpha); qq(2) = to_alpha(pb, alpha)
  oo(:, 1) = o_alpha(qq(1), alpha); oo(:, 2) = o_alpha(qq(2), alpha)
  do i = 1, 2
    do l = 1, 9
      dsp(l, i) = qq(i)%sqrdel(orb_l(l))**2; csp(l, i) = qq(i)%c(orb_l(l)) - qq(i)%enu(orb_l(l))
      esp(l, i) = qq(i)%enu(orb_l(l)); osp(l, i) = oo(orb_l(l), i); psp(l, i) = qq(i)%p(orb_l(l))
    end do
  end do
  call check(maxval(abs(psp - reshape([(pa%p(orb_l(l)) + osp(l, 1)**2, l=1, 9), (pb%p(orb_l(l)) + osp(l, 2)**2, l=1, 9)], &
                                      [9, 2]))) < 1.0e-12_dp, 'to_alpha: p^alpha = p^gamma + o^2 (from the gamma representation)')
  x = 0.35_dp
  allocate (g%nbr_ptr(n + 1), g%nbr(2 * n)); g%n = n          ! ring: slot 1 -> R+1, slot 2 -> R-1
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
  allocate (mdl%c0(9), mdl%c1(9), mdl%d0(9), mdl%d1(9), mdl%o0(9), mdl%o1(9), mdl%e0(9), mdl%e1(9))
  zero9 = (0.0_dp, 0.0_dp)

  esym = 0; eg = 0; ef = 0; em = 0; ech = 0; esw = 0; ereg = 0; ecan = 0; ewu = 0
  do cen = 1, 2
    dcen = dsp(:, cen)
    do i = 1, 2; gsp(:, i) = dcen / dsp(:, i) * csp(:, i); qsp(:, i) = osp(:, i) * dsp(:, i) / dcen; end do
    mdl%c0 = gsp(:, 2); mdl%c1 = gsp(:, 1) - gsp(:, 2); mdl%d0 = sqrt(dcen); mdl%d1 = 0.0_dp
    mdl%o0 = 0.0_dp; mdl%o1 = 0.0_dp; mdl%e0 = 0.0_dp; mdl%e1 = 0.0_dp
    block
      logical :: flippable(n)
      flippable = .true.; flippable(1) = .false.
      call build_augorb_operator(g, mdl, 1, 60, op, flippable=flippable)
    end block
    call augorb_freeze_site(op, 1, merge(1.0_dp, 0.0_dp, cen == 1))
    ndim = op%dim()
    print '(a,i0,a,i0,a)', '  central species ', cen, ': ', op%nn, ' augmented nodes'
    call check(op%nn == n * 2**(n - 1), 'node set = N 2^(N-1) (start site frozen)')
    if (allocated(amat)) deallocate (amat, afac, rhs, uvec, wvec, ej, cvec, rv, ipiv)
    allocate (amat(ndim, ndim), afac(ndim, ndim), rhs(ndim, 2), uvec(ndim), wvec(ndim), ej(ndim), cvec(ndim), rv(ndim), &
              ipiv(ndim))
    do il = 1, 3
      l0 = lstart(il)
      do kz = 1, 4
        z = cmplx(-0.3_dp + 0.25_dp * kz, 0.3_dp, dp)
        do i = 1, 2
          et(:, i) = dcen / dsp(:, i) * (z - esp(:, i))
          uu(:, i) = osp(:, i) * (z - esp(:, i))
          qt(:, i) = dsp(:, i) / dcen * (osp(:, i) - psp(:, i) * (z - esp(:, i)))
        end do
        ! (1) dense augmented K'
        do j = 1, ndim
          ej = (0.0_dp, 0.0_dp); ej(j) = (1.0_dp, 0.0_dp)
          call kprime(ej, et, uu, qt, amat(:, j))
        end do
        esym = max(esym, maxval(abs(amat - transpose(amat))))
        ! (2) dense solve against brute force
        uvec = (0.0_dp, 0.0_dp); uvec(l0) = (1.0_dp, 0.0_dp)
        call hprime(uvec, wvec); wvec = wvec - gsp(l0, cen) * uvec
        rhs(:, 1) = uvec; rhs(:, 2) = wvec; afac = amat
        call zgesv(ndim, 2, afac, ndim, ipiv, rhs, ndim, info)
        if (info /= 0) error stop 'zgesv (augmented K'')'
        guu = rhs(l0, 1); gwu = sum(wvec * rhs(:, 1)); gww = sum(wvec * rhs(:, 2)); guw = rhs(l0, 2)
        faug = gsp(l0, cen) * guu + gwu
        maug = gsp(l0, cen)**2 * guu + 2.0_dp * gsp(l0, cen) * gwu + gww
        call brute(cen, z, l0, .false., b11, bh, bhh, b2)
        eg = max(eg, abs(guu - b11)); ef = max(ef, abs(faug - bh)); em = max(em, abs(maug - bhh))
        ewu = max(ewu, abs(gwu - guw))
        ! (3) two Lanczos chains (full re-orthogonalisation in the bilinear form)
        eta = sqrt(real(sum(wvec * wvec), dp))
        call chain(amat, uvec, wvec, cuu, cwu)
        call chain(amat, wvec / eta, uvec, cww, cuw)
        cww = eta**2 * cww; cuw = eta * cuw
        ech = max(ech, abs(cuu - guu), abs(cwu - gwu), abs(cww - gww))
        esw = max(esw, abs(cuw - cwu))
        ! (4) regression: u = 0, p = 0, Qt -> -q  is the asr_edc operator  et - h' + h' q h'
        do j = 1, ndim
          ej = (0.0_dp, 0.0_dp); ej(j) = (1.0_dp, 0.0_dp)
          call kprime(ej, et, zero9, cmplx(-qsp, 0.0_dp, dp), afac(:, j))
        end do
        rhs(:, 1) = uvec
        call zgesv(ndim, 1, afac, ndim, ipiv, rhs(:, 1), ndim, info)
        if (info /= 0) error stop 'zgesv (regression)'
        ereg = max(ereg, abs(rhs(l0, 1) - b2))
        ! (5) canonical form against the factored product, on a random complex vector
        do j = 1, ndim
          call random_number(u); call random_number(tmp); cvec(j) = cmplx(u - 0.5_dp, tmp - 0.5_dp, dp)
        end do
        call kcanon(cvec, z, rv)
        ecan = max(ecan, maxval(abs(rv - matmul(amat, cvec))))
      end do
    end do
  end do
  print '(a,es10.3)', '  max |K''_aug - K''_aug^T|               ', esym
  print '(a,3es10.3)', '  dense K''_aug vs brute G11, (hG)11, (hGh)11', eg, ef, em
  print '(a,es10.3)', '  |G_wu - G_uw| (dense)                   ', ewu
  print '(a,2es10.3)', '  two chains vs dense, |G_uw - G_wu| chains', ech, esw
  print '(a,es10.3)', '  regression (asr_edc operator) vs <[E - H2]^-1>', ereg
  print '(a,es10.3)', '  canonical form vs factored product      ', ecan
  call check_close(esym, 0.0_dp, 1.0e-13_dp, 'K''_aug is complex symmetric')
  call check_close(eg, 0.0_dp, 1.0e-10_dp, 'augmented <u|K''^-1|u> == configuration average of G11 = [E O - H]^-1_11')
  call check_close(ef, 0.0_dp, 1.0e-10_dp, 'augmented g_c G_uu + G_wu == average of (h G)11 (w = S'' u deterministic)')
  call check_close(em, 0.0_dp, 1.0e-10_dp, 'augmented g_c^2 G_uu + 2 g_c G_wu + G_ww == average of (h G h)11')
  call check_close(ewu, 0.0_dp, 1.0e-10_dp, 'G_wu == G_uw (complex symmetry)')
  call check_close(ech, 0.0_dp, 1.0e-8_dp, 'two Lanczos chains with projection == dense solve')
  call check_close(esw, 0.0_dp, 1.0e-8_dp, 'cross term from chain u == cross term from chain w')
  call check_close(ereg, 0.0_dp, 1.0e-10_dp, 'u = 0, p = 0, Qt -> -q: K'' is the asr_edc operator (== <[E - H2]^-1_11>)')
  call check_close(ecan, 0.0_dp, 1.0e-12_dp, 'canonical form (V = 1/2 + c o - (o + c p) eps) == factored product')

  ! (6) single-site dressing: exact at p^gamma = 0, measured at p^gamma /= 0
  do cen = 1, 2
    call dressing(cen, .true., eq0, eqrel, peps)
    call check_close(eq0, 0.0_dp, 1.0e-10_dp, 'single-site dressing exact at p^gamma = 0 (sphere charge of site 1)')
    call dressing(cen, .false., eqmax, eqrel, peps)
    print '(a,i0,a,es10.3,a,es10.3,a,f7.4)', '  central ', cen, ': dressing error at p^gamma /= 0: |dQ| max ', eqmax, &
      '  relative ', eqrel, '  max p^gamma eps^2 (occupied) ', peps
  end do
  call report()

contains

  !> y = h' x (complex x; asr_augorb with o = E_nu = 0 applies h' = g + S' to real vectors)
  subroutine hprime(xv, yv)
    complex(dp), intent(in) :: xv(:)
    complex(dp), intent(out) :: yv(:)
    real(dp) :: r1(size(xv)), r2(size(xv)), y1(size(xv)), y2(size(xv))
    r1 = real(xv); r2 = aimag(xv)
    call op%apply(r1, y1); call op%apply(r2, y2)
    yv = cmplx(y1, y2, dp)
  end subroutine hprime

  !> y = X x for the random site-diagonal X = X_B + (X_A - X_B) M (per orbital, complex coefficients)
  subroutine rdiag(xs, xv, yv)
    complex(dp), intent(in) :: xs(9, 2), xv(:)
    complex(dp), intent(out) :: yv(:)
    real(dp) :: r1(size(xv)), r2(size(xv)), y1(size(xv)), y2(size(xv))
    integer :: ii, ll
    r1 = real(xv); r2 = aimag(xv)
    call op%mapply(r1, y1); call op%mapply(r2, y2)
    do ii = 1, op%nn
      do ll = 1, 9
        yv(9 * (ii - 1) + ll) = xs(ll, 2) * xv(9 * (ii - 1) + ll) + (xs(ll, 1) - xs(ll, 2)) * cmplx(y1(9 * (ii - 1) + ll), &
                                y2(9 * (ii - 1) + ll), dp)
      end do
    end do
  end subroutine rdiag

  !> y = K' x = et x - h'x + u h'x + h'(u x - Qt h'x)
  subroutine kprime(xv, ets, us, qts, yv)
    complex(dp), intent(in) :: xv(:), ets(9, 2), us(9, 2), qts(9, 2)
    complex(dp), intent(out) :: yv(:)
    complex(dp) :: t1(size(xv)), t2(size(xv)), t3(size(xv)), t4(size(xv))
    call hprime(xv, t1)
    call rdiag(qts, t1, t2)
    call rdiag(us, xv, t3)
    t3 = t3 - t2
    call hprime(t3, t4)
    call rdiag(ets, xv, yv)
    call rdiag(us, t1, t2)
    yv = yv - t1 + t2 + t4
  end subroutine kprime

  !> canonical form  y = D x - (V S'x + S'(V x)) - S'(Qt S'x),  S'x = h'x - g x
  subroutine kcanon(xv, zz, yv)
    complex(dp), intent(in) :: xv(:), zz
    complex(dp), intent(out) :: yv(:)
    complex(dp) :: dd(9, 2), vv(9, 2), gg(9, 2), qq_(9, 2), ep(9, 2)
    complex(dp) :: sx(size(xv)), t1(size(xv)), t2(size(xv)), t3(size(xv))
    integer :: ii
    do ii = 1, 2
      ep(:, ii) = zz - esp(:, ii)
      dd(:, ii) = dcen / dsp(:, ii) * (ep(:, ii) - csp(:, ii) + 2.0_dp * osp(:, ii) * csp(:, ii) * ep(:, ii) &
                  - csp(:, ii)**2 * (osp(:, ii) - psp(:, ii) * ep(:, ii)))
      vv(:, ii) = (0.5_dp + csp(:, ii) * osp(:, ii)) - (osp(:, ii) + csp(:, ii) * psp(:, ii)) * ep(:, ii)
      gg(:, ii) = cmplx(gsp(:, ii), 0.0_dp, dp)
      qq_(:, ii) = dsp(:, ii) / dcen * (osp(:, ii) - psp(:, ii) * ep(:, ii))
    end do
    call hprime(xv, sx); call rdiag(gg, xv, t1); sx = sx - t1            ! S' x
    call rdiag(dd, xv, yv)
    call rdiag(vv, sx, t1); yv = yv - t1                                  ! - V S' x
    call rdiag(vv, xv, t2); call hprime(t2, t3); call rdiag(gg, t2, t1); t3 = t3 - t1; yv = yv - t3   ! - S' V x
    call rdiag(qq_, sx, t2); call hprime(t2, t3); call rdiag(gg, t2, t1); t3 = t3 - t1; yv = yv - t3  ! - S' Qt S' x
  end subroutine kcanon

  !> complex-symmetric Lanczos on the dense a from s0 (s0^T s0 = 1), full re-orthogonalisation; gss = [J^-1]_11,
  !> gts = sum_k (t^T v_k) [J^-1]_k1  (A^-1 s0 = V J^-1 e1 once the Krylov space is invariant)
  subroutine chain(a, s0, t, gss, gts)
    complex(dp), intent(in) :: a(:, :), s0(:), t(:)
    complex(dp), intent(out) :: gss, gts
    complex(dp), allocatable :: v(:, :), aa(:), bb(:), jm(:, :), y(:), wv(:)
    integer :: k, m, jj, np, inf
    integer, allocatable :: ip(:)
    np = size(s0)
    allocate (v(np, np + 1), aa(np), bb(np), wv(np))
    v(:, 1) = s0; m = np
    do k = 1, np
      wv = matmul(a, v(:, k))
      aa(k) = sum(v(:, k) * wv)
      wv = wv - aa(k) * v(:, k)
      if (k > 1) wv = wv - bb(k - 1) * v(:, k - 1)
      do jj = 1, k; wv = wv - sum(v(:, jj) * wv) * v(:, jj); end do
      do jj = 1, k; wv = wv - sum(v(:, jj) * wv) * v(:, jj); end do
      bb(k) = sqrt(sum(wv * wv))
      if (abs(bb(k)) < 1.0e-10_dp) then; m = k; exit; end if
      v(:, k + 1) = wv / bb(k)
    end do
    allocate (jm(m, m), y(m), ip(m)); jm = (0.0_dp, 0.0_dp)
    do k = 1, m
      jm(k, k) = aa(k)
      if (k < m) then; jm(k, k + 1) = bb(k); jm(k + 1, k) = bb(k); end if
    end do
    y = (0.0_dp, 0.0_dp); y(1) = (1.0_dp, 0.0_dp)
    call zgesv(m, 1, jm, m, ip, y, m, inf)
    if (inf /= 0) error stop 'zgesv (chain)'
    gss = y(1)
    gts = (0.0_dp, 0.0_dp)
    do k = 1, m; gts = gts + sum(t * v(:, k)) * y(k); end do
  end subroutine chain

  !> dense ring matrices for one configuration: h, O, H (bands_core form); pgz: p = o^2 (p^gamma = 0)
  subroutine dense_ops(isa_, pgz, hm, om, hmat, en, ov, pv)
    logical, intent(in) :: isa_(n), pgz
    real(dp), intent(out) :: hm(nd, nd), om(nd, nd), hmat(nd, nd), en(nd), ov(nd), pv(nd)
    real(dp) :: s(nd, nd), dl(nd), cm(nd)
    integer :: ii, jj, kk, ll, ll2, spc, k
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
        dl(k) = dsp(ll, spc); cm(k) = csp(ll, spc); en(k) = esp(ll, spc); ov(k) = osp(ll, spc)
        pv(k) = psp(ll, spc); if (pgz) pv(k) = osp(ll, spc)**2
      end do
    end do
    do ii = 1, nd
      do jj = 1, nd; hm(ii, jj) = sqrt(dl(ii)) * s(ii, jj) * sqrt(dl(jj)); end do
      hm(ii, ii) = hm(ii, ii) + cm(ii)
    end do
    do jj = 1, nd
      do ii = 1, nd
        om(ii, jj) = (ov(ii) + ov(jj)) * hm(ii, jj) + sum(hm(ii, :) * pv * hm(:, jj))
        hmat(ii, jj) = (1.0_dp + en(ii) * ov(ii) + en(jj) * ov(jj)) * hm(ii, jj) + sum(hm(ii, :) * (ov + pv * en) * hm(:, jj))
      end do
      om(jj, jj) = om(jj, jj) + 1.0_dp
      hmat(jj, jj) = hmat(jj, jj) + en(jj)
    end do
  end subroutine dense_ops

  !> configuration averages of G11, (hG)11, (hGh)11 of K = E O - H and of [E - H^(2)]^-1_11 (orbital l0_ on site 1)
  subroutine brute(cen_, zz, l0_, pgz, g11, gh, ghh, g2)
    integer, intent(in) :: cen_, l0_
    complex(dp), intent(in) :: zz
    logical, intent(in) :: pgz
    complex(dp), intent(out) :: g11, gh, ghh, g2
    real(dp) :: hm(nd, nd), om(nd, nd), hmat(nd, nd), en(nd), ov(nd), pv(nd), h2(nd, nd), w
    complex(dp) :: kk(nd, nd), r(nd, 2), a2(nd, nd), r2(nd)
    integer :: c, ii, jj, ip(nd), inf
    logical :: isa_(n)
    g11 = 0; gh = 0; ghh = 0; g2 = 0
    do c = 0, 2**(n - 1) - 1
      w = 1.0_dp; isa_(1) = (cen_ == 1)
      do ii = 2, n
        isa_(ii) = btest(c, ii - 2)
        w = w * merge(x, 1.0_dp - x, isa_(ii))
      end do
      call dense_ops(isa_, pgz, hm, om, hmat, en, ov, pv)
      kk = zz * om - hmat
      r = (0.0_dp, 0.0_dp); r(l0_, 1) = (1.0_dp, 0.0_dp); r(:, 2) = hm(:, l0_)
      call zgesv(nd, 2, kk, nd, ip, r, nd, inf)
      if (inf /= 0) error stop 'zgesv (brute)'
      g11 = g11 + w * r(l0_, 1)
      gh = gh + w * sum(hm(l0_, :) * r(:, 1))
      ghh = ghh + w * sum(hm(:, l0_) * r(:, 2))
      do jj = 1, nd; do ii = 1, nd; h2(ii, jj) = hm(ii, jj) - sum(hm(ii, :) * ov * hm(:, jj)); end do; end do
      a2 = -h2
      do ii = 1, nd; a2(ii, ii) = a2(ii, ii) + zz - en(ii); end do
      r2 = (0.0_dp, 0.0_dp); r2(l0_) = (1.0_dp, 0.0_dp)
      call zgesv(nd, 1, a2, nd, ip, r2, nd, inf)
      g2 = g2 + w * r2(l0_)
    end do
  end subroutine brute

  !> Sphere charge of site 1 summed over the occupied half of the states, exact (z, h z) against the single-site dressing
  !> (h z)_1 = eps_1 z_1/(1 - eps_1 o_1); configuration averaged per orbital.  dqmax: max |dQ|, rel: max |dQ|/Q.
  subroutine dressing(cen_, pgz, dqmax, rel, pe2)
    integer, intent(in) :: cen_
    logical, intent(in) :: pgz
    real(dp), intent(out) :: dqmax, rel, pe2
    real(dp) :: hm(nd, nd), om(nd, nd), hmat(nd, nd), en(nd), ov(nd), pv(nd), ev(nd), hz(nd, nd), w, ep, f, zl, q(9), qd(9)
    real(dp), allocatable :: work(:)
    integer :: c, ii, k, ll, inf
    logical :: isa_(n)
    external :: dsygv
    allocate (work(64 * nd))
    q = 0; qd = 0; pe2 = 0
    do c = 0, 2**(n - 1) - 1
      w = 1.0_dp; isa_(1) = (cen_ == 1)
      do ii = 2, n
        isa_(ii) = btest(c, ii - 2)
        w = w * merge(x, 1.0_dp - x, isa_(ii))
      end do
      call dense_ops(isa_, pgz, hm, om, hmat, en, ov, pv)
      call dsygv(1, 'V', 'U', nd, hmat, nd, om, nd, ev, work, size(work), inf)
      if (inf /= 0) error stop 'dsygv (dressing)'
      hz = matmul(hm, hmat)
      do k = 1, nd / 2
        do ll = 1, 9
          zl = hmat(ll, k); ep = ev(k) - en(ll); f = ep / (1.0_dp - ep * ov(ll))
          q(ll) = q(ll) + w * (zl**2 + 2.0_dp * ov(ll) * zl * hz(ll, k) + pv(ll) * hz(ll, k)**2)
          qd(ll) = qd(ll) + w * zl**2 * (1.0_dp + 2.0_dp * ov(ll) * f + pv(ll) * f**2)
          pe2 = max(pe2, (pv(ll) - ov(ll)**2) * ep**2)
        end do
      end do
    end do
    dqmax = maxval(abs(qd - q)); rel = maxval(abs(qd - q) / max(abs(q), 1.0e-12_dp))
  end subroutine dressing
end program test_asr_ed3
