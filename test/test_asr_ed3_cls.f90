!> Gate 1, part B, of docs/DERIVATION_ed3.md: the full-ASA energy-dependent operator on the site-class operator asr_augcls
!> (N species, site classes, first-shell SRO, point-group reduction).  K' = et - h' + u h' + h' u - h' Qt h' with
!> h' = g + S' (augcls apply with o = E_nu = 0), random diagonals by augcls_rand_c (complex, normalised coordinates).
!>  (B1) ring of 5, four species: classes ordered {1} | binary {2,3} | ternary {1,3,4}; site 1 (ternary) frozen to each of its
!>       species; its two ring neighbours are first-shell classes with conditional concentrations (zero weights included).
!>       Dense augmented K' and the two Lanczos chains == exact conditional configuration average of G11, (hG)11, (hGh)11 of
!>       K = E O - H (bands_core form); complex symmetry; u = 0, p = 0, Qt -> -q == <[E - H2]^-1_11>; canonical form.
!>  (B2) fcc, ternary {2,3,4}, conditional first-shell class: D4h-reduced chains == unreduced chains for G_uu, G_wu, G_ww, G_uw
!>       (the start vector w = h'u - g_c u is built by the operator itself, so it is in normalised orbit coordinates).
program test_asr_ed3_cls
  use asr_kinds, only: dp
  use asr_augspace, only: graph_t
  use asr_augcls
  use asr_augorb, only: orb_alloy_t
  use asr_lmto_alloy, only: build_lmto_alloy
  use asr_symmetry
  use asr_lmto_ham, only: potpar_t, to_alpha, o_alpha
  use asr_testing
  implicit none
  integer, parameter :: orb_l(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2], l0s(4) = [1, 3, 5, 7]
  type(potpar_t) :: pp(4), qq(4)
  real(dp) :: alpha(0:2), oo(0:2, 4)
  real(dp) :: dsp(9, 4), csp(9, 4), esp(9, 4), osp(9, 4), psp(9, 4)
  integer :: l, sp
  external :: zgesv

  print '(a)', 'test_asr_ed3_cls'
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
    do l = 1, 9
      dsp(l, sp) = qq(sp)%sqrdel(orb_l(l))**2; csp(l, sp) = qq(sp)%c(orb_l(l)) - qq(sp)%enu(orb_l(l))
      esp(l, sp) = qq(sp)%enu(orb_l(l)); osp(l, sp) = oo(orb_l(l), sp); psp(l, sp) = qq(sp)%p(orb_l(l))
    end do
  end do
  call ring_sro()
  call fcc_reduced()
  call report()

contains

  !> operator parameters for central species cen: c = g_s = Dc/Ds (C - E_nu)_s, d = sqrt(Dc), o = E_nu = 0  (apply == h')
  subroutine set_hprime(mdl, cen)
    type(class_alloy_t), intent(inout) :: mdl
    integer, intent(in) :: cen
    integer :: s
    if (allocated(mdl%c)) deallocate (mdl%c, mdl%d, mdl%o, mdl%e)
    allocate (mdl%c(9, 4), mdl%d(9, 4), mdl%o(9, 4), mdl%e(9, 4))
    do s = 1, 4
      mdl%c(:, s) = dsp(:, cen) / dsp(:, s) * csp(:, s); mdl%d(:, s) = sqrt(dsp(:, cen))
    end do
    mdl%o = 0.0_dp; mdl%e = 0.0_dp
  end subroutine set_hprime

  !> energy-dependent species coefficients et, u, Qt (and -q for the regression) for central species cen
  subroutine coefs(cen, zz, et, uu, qt, mq)
    integer, intent(in) :: cen
    complex(dp), intent(in) :: zz
    complex(dp), intent(out) :: et(9, 4), uu(9, 4), qt(9, 4), mq(9, 4)
    integer :: s
    do s = 1, 4
      et(:, s) = dsp(:, cen) / dsp(:, s) * (zz - esp(:, s))
      uu(:, s) = osp(:, s) * (zz - esp(:, s))
      qt(:, s) = dsp(:, s) / dsp(:, cen) * (osp(:, s) - psp(:, s) * (zz - esp(:, s)))
      mq(:, s) = -osp(:, s) * dsp(:, s) / dsp(:, cen)
    end do
  end subroutine coefs

  subroutine hprime(op, xv, yv)
    type(augcls_op), intent(in) :: op
    complex(dp), intent(in) :: xv(:)
    complex(dp), intent(out) :: yv(:)
    real(dp) :: r1(size(xv)), r2(size(xv)), y1(size(xv)), y2(size(xv))
    r1 = real(xv); r2 = aimag(xv)
    call op%apply(r1, y1); call op%apply(r2, y2)
    yv = cmplx(y1, y2, dp)
  end subroutine hprime

  !> y = K' x = et x - h'x + u h'x + h'(u x - Qt h'x)
  subroutine kprime(op, xv, et, uu, qt, yv)
    type(augcls_op), intent(in) :: op
    complex(dp), intent(in) :: xv(:), et(:, :), uu(:, :), qt(:, :)
    complex(dp), intent(out) :: yv(:)
    complex(dp) :: t1(size(xv)), t2(size(xv)), t3(size(xv)), t4(size(xv))
    call hprime(op, xv, t1)
    call augcls_rand_c(op, qt, t1, t2)
    call augcls_rand_c(op, uu, xv, t3)
    t3 = t3 - t2
    call hprime(op, t3, t4)
    call augcls_rand_c(op, et, xv, yv)
    call augcls_rand_c(op, uu, t1, t2)
    yv = yv - t1 + t2 + t4
  end subroutine kprime

  !> canonical form  D x - (V S'x + S'(V x)) - S'(Qt S'x),  S'x = h'x - g x
  subroutine kcanon(op, cen, xv, zz, yv)
    type(augcls_op), intent(in) :: op
    integer, intent(in) :: cen
    complex(dp), intent(in) :: xv(:), zz
    complex(dp), intent(out) :: yv(:)
    complex(dp) :: dd(9, 4), vv(9, 4), gg(9, 4), qt(9, 4), ep(9, 4)
    complex(dp) :: sx(size(xv)), t1(size(xv)), t2(size(xv)), t3(size(xv))
    integer :: s
    do s = 1, 4
      ep(:, s) = zz - esp(:, s)
      dd(:, s) = dsp(:, cen) / dsp(:, s) * (ep(:, s) - csp(:, s) + 2.0_dp * osp(:, s) * csp(:, s) * ep(:, s) &
                 - csp(:, s)**2 * (osp(:, s) - psp(:, s) * ep(:, s)))
      vv(:, s) = (0.5_dp + csp(:, s) * osp(:, s)) - (osp(:, s) + csp(:, s) * psp(:, s)) * ep(:, s)
      gg(:, s) = cmplx(dsp(:, cen) / dsp(:, s) * csp(:, s), 0.0_dp, dp)
      qt(:, s) = dsp(:, s) / dsp(:, cen) * (osp(:, s) - psp(:, s) * ep(:, s))
    end do
    call hprime(op, xv, sx); call augcls_rand_c(op, gg, xv, t1); sx = sx - t1
    call augcls_rand_c(op, dd, xv, yv)
    call augcls_rand_c(op, vv, sx, t1); yv = yv - t1
    call augcls_rand_c(op, vv, xv, t2); call hprime(op, t2, t3); call augcls_rand_c(op, gg, t2, t1); yv = yv - (t3 - t1)
    call augcls_rand_c(op, qt, sx, t2); call hprime(op, t2, t3); call augcls_rand_c(op, gg, t2, t1); yv = yv - (t3 - t1)
  end subroutine kcanon

  !> complex-symmetric Lanczos on K' from s0 (s0^T s0 = 1), nlev steps; reorth: full re-orthogonalisation.
  !> gss = [J^-1]_11, gts = sum_k (t^T v_k) [J^-1]_k1.
  subroutine chain(op, et, uu, qt, s0, t, nlev, reorth, gss, gts, acoef, bcoef, mused)
    type(augcls_op), intent(in) :: op
    complex(dp), intent(in) :: et(:, :), uu(:, :), qt(:, :), s0(:), t(:)
    integer, intent(in) :: nlev
    logical, intent(in) :: reorth
    complex(dp), intent(out) :: gss, gts
    complex(dp), intent(out), optional :: acoef(:), bcoef(:)
    integer, intent(out), optional :: mused
    complex(dp), allocatable :: v(:, :), aa(:), bb(:), jm(:, :), y(:), wv(:), gam(:)
    integer :: k, m, jj, np, inf
    integer, allocatable :: ip(:)
    np = size(s0)
    allocate (v(np, merge(nlev + 1, 2, reorth)), aa(nlev), bb(nlev), wv(np), gam(nlev))
    v(:, 1) = s0; m = nlev
    block
      complex(dp) :: vk(np), vprev(np)
      vk = s0; vprev = (0.0_dp, 0.0_dp)
      do k = 1, nlev
        gam(k) = sum(t * vk)
        call kprime(op, vk, et, uu, qt, wv)
        aa(k) = sum(vk * wv)
        wv = wv - aa(k) * vk
        if (k > 1) wv = wv - bb(k - 1) * vprev
        if (reorth) then
          do jj = 1, k; wv = wv - sum(v(:, jj) * wv) * v(:, jj); end do
          do jj = 1, k; wv = wv - sum(v(:, jj) * wv) * v(:, jj); end do
        end if
        if (k == nlev) exit
        bb(k) = sqrt(sum(wv * wv))
        if (abs(bb(k)) < 1.0e-10_dp) then; m = k; exit; end if
        vprev = vk; vk = wv / bb(k)
        if (reorth) v(:, k + 1) = vk
      end do
    end block
    allocate (jm(m, m), y(m), ip(m)); jm = (0.0_dp, 0.0_dp)
    do k = 1, m
      jm(k, k) = aa(k)
      if (k < m) then; jm(k, k + 1) = bb(k); jm(k + 1, k) = bb(k); end if
    end do
    y = (0.0_dp, 0.0_dp); y(1) = (1.0_dp, 0.0_dp)
    call zgesv(m, 1, jm, m, ip, y, m, inf)
    if (inf /= 0) error stop 'zgesv (chain)'
    gss = y(1); gts = sum(gam(1:m) * y)
    if (present(acoef)) acoef(1:m) = aa(1:m)
    if (present(bcoef)) bcoef(1:m - 1) = bb(1:m - 1)
    if (present(mused)) mused = m
  end subroutine chain

  !> one configuration of the ring (species cf per site): K = E O - H (bands_core form) and E - H2;
  !> returns G11, (hG)11, (hGh)11 of K and [E - H2]^-1_11 for orbital l0 of site 1
  subroutine dense_blocks(mdl, g, cf, zz, l0, g11, gh, ghh, g2)
    type(class_alloy_t), intent(in) :: mdl
    type(graph_t), intent(in) :: g
    integer, intent(in) :: cf(:), l0
    complex(dp), intent(in) :: zz
    complex(dp), intent(out) :: g11, gh, ghh, g2
    integer :: n, nd, ii, jj, kk, ll, kq, inf
    real(dp), allocatable :: s(:, :), hm(:, :), om(:, :), hmat(:, :), dl(:), cm(:), en(:), ov(:), pv(:)
    complex(dp), allocatable :: ka(:, :), r(:, :), a2(:, :), r2(:)
    integer, allocatable :: ip(:)
    n = g%n; nd = 9 * n
    allocate (s(nd, nd), hm(nd, nd), om(nd, nd), hmat(nd, nd), dl(nd), cm(nd), en(nd), ov(nd), pv(nd), ka(nd, nd), &
              r(nd, 2), a2(nd, nd), r2(nd), ip(nd))
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
        dl(kq) = dsp(ll, cf(ii)); cm(kq) = csp(ll, cf(ii)); en(kq) = esp(ll, cf(ii))
        ov(kq) = osp(ll, cf(ii)); pv(kq) = psp(ll, cf(ii))
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
        a2(ii, jj) = -(hm(ii, jj) - sum(hm(ii, :) * ov * hm(:, jj)))
      end do
      om(jj, jj) = om(jj, jj) + 1.0_dp
      hmat(jj, jj) = hmat(jj, jj) + en(jj)
      a2(jj, jj) = a2(jj, jj) + zz - en(jj)
    end do
    ka = zz * om - hmat
    r = (0.0_dp, 0.0_dp); r(l0, 1) = (1.0_dp, 0.0_dp); r(:, 2) = hm(:, l0)
    call zgesv(nd, 2, ka, nd, ip, r, nd, inf)
    if (inf /= 0) error stop 'zgesv (dense K)'
    g11 = r(l0, 1); gh = sum(hm(l0, :) * r(:, 1)); ghh = sum(hm(:, l0) * r(:, 2))
    r2 = (0.0_dp, 0.0_dp); r2(l0) = (1.0_dp, 0.0_dp)
    call zgesv(nd, 1, a2, nd, ip, r2, nd, inf)
    if (inf /= 0) error stop 'zgesv (dense H2)'
    g2 = r2(l0)
  end subroutine dense_blocks

  !> local states of zero-weight species have empty rows and columns (class_mt): they never couple to the start sector.
  !> Put 1 on their diagonal so that the dense matrix is invertible; dec marks them.
  subroutine fill_decoupled(a, dec, nz)
    complex(dp), intent(inout) :: a(:, :)
    logical, intent(out) :: dec(:)
    integer, intent(out) :: nz
    integer :: jj
    nz = 0
    do jj = 1, size(a, 1)
      dec(jj) = all(a(:, jj) == (0.0_dp, 0.0_dp)) .and. all(a(jj, :) == (0.0_dp, 0.0_dp))
      if (dec(jj)) then; a(jj, jj) = (1.0_dp, 0.0_dp); nz = nz + 1; end if
    end do
  end subroutine fill_decoupled

  subroutine ring_sro()
    integer, parameter :: n = 5, rcls(n) = [3, 4, 1, 2, 5]
    real(dp), parameter :: c4(3, 3) = reshape([0.0_dp, 0.7_dp, 0.3_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.2_dp, 0.6_dp, 0.2_dp], [3, 3])
    real(dp), parameter :: c5(2, 3) = reshape([0.1_dp, 0.9_dp, 1.0_dp, 0.0_dp, 0.55_dp, 0.45_dp], [2, 3])
    type(graph_t) :: g
    type(class_alloy_t) :: mdl
    type(augcls_op) :: op
    integer :: perm(n, 1), i, l2, j, cen, kz, cf(n), q, c, ncfg, l0, il, nd, jc, inf
    real(dp) :: dmat(9, 9, 1), u, v, pc, eta, gc
    real(dp) :: esym, eb(3), ech, esw, ereg, ecan
    complex(dp) :: z, et(9, 4), uu(9, 4), qt(9, 4), mq(9, 4), b11, bh, bhh, b2, d11, dh, dhh, d2
    complex(dp) :: guu, gwu, gww, guw, cuu, cwu, cww, cuw
    complex(dp), allocatable :: amat(:, :), afac(:, :), rhs(:, :), uvec(:), wvec(:), ej(:), cvec(:), rv(:)
    integer, allocatable :: ipiv(:)
    logical, allocatable :: dec(:)
    integer :: ndec
    logical :: fl(n)
    allocate (g%nbr_ptr(n + 1), g%nbr(2 * n)); g%n = n
    do i = 1, n
      g%nbr_ptr(i) = 2 * (i - 1) + 1
      g%nbr(2 * i - 1) = mod(i, n) + 1
      g%nbr(2 * i) = mod(i - 2 + n, n) + 1
    end do
    g%nbr_ptr(n + 1) = 2 * n + 1
    mdl%nl = 9; mdl%nsp = 4; mdl%ncls = 5
    allocate (mdl%nsc(5), mdl%spc(3, 5), mdl%conc(3, 5)); mdl%spc = 0; mdl%conc = 0.0_dp
    mdl%nsc = [1, 2, 3, 3, 2]
    mdl%spc(1:1, 1) = [1];       mdl%conc(1:1, 1) = [1.0_dp]
    mdl%spc(1:2, 2) = [2, 3];    mdl%conc(1:2, 2) = [0.4_dp, 0.6_dp]
    mdl%spc(1:3, 3) = [1, 3, 4]; mdl%conc(1:3, 3) = [0.3_dp, 0.5_dp, 0.2_dp]
    mdl%spc(:, 4) = mdl%spc(:, 3); mdl%conc(:, 4) = mdl%conc(:, 3)
    mdl%spc(:, 5) = mdl%spc(:, 2); mdl%conc(:, 5) = mdl%conc(:, 2)
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
    esym = 0; eb = 0; ech = 0; esw = 0; ereg = 0; ecan = 0
    do j = 1, 3
      cen = mdl%spc(j, 3)
      call set_hprime(mdl, cen)
      call build_augcls_operator(g, mdl, rcls, perm, dmat, 1, 80, op, flippable=fl)
      call augcls_freeze_site(op, 1, cen)
      call augcls_set_class_conc(op, 4, c4(:, j)); call augcls_set_class_conc(op, 5, c5(:, j))
      nd = op%dim()
      if (j == 1) print '(a,i0,a)', '  ring SRO: ', op%nn, ' augmented nodes'
      if (allocated(amat)) deallocate (amat, afac, rhs, uvec, wvec, ej, cvec, rv, ipiv)
      allocate (amat(nd, nd), afac(nd, nd), rhs(nd, 2), uvec(nd), wvec(nd), ej(nd), cvec(nd), rv(nd), ipiv(nd))
      if (allocated(dec)) deallocate (dec)
      allocate (dec(nd))
      do il = 1, 2
        l0 = merge(1, 5, il == 1)
        call augcls_set_start(op, l0, dmat)
        gc = dsp(l0, cen) / dsp(l0, cen) * csp(l0, cen)
        do kz = 1, 3
          z = cmplx(-0.3_dp + 0.3_dp * kz, 0.3_dp, dp)
          call coefs(cen, z, et, uu, qt, mq)
          do jc = 1, nd
            ej = (0.0_dp, 0.0_dp); ej(jc) = (1.0_dp, 0.0_dp)
            call kprime(op, ej, et, uu, qt, amat(:, jc))
          end do
          esym = max(esym, maxval(abs(amat - transpose(amat))))
          call fill_decoupled(amat, dec, ndec)
          uvec = (0.0_dp, 0.0_dp); uvec(l0) = (1.0_dp, 0.0_dp)
          call hprime(op, uvec, wvec); wvec = wvec - gc * uvec
          rhs(:, 1) = uvec; rhs(:, 2) = wvec; afac = amat
          call zgesv(nd, 2, afac, nd, ipiv, rhs, nd, inf)
          if (inf /= 0) error stop 'zgesv (ring K'')'
          guu = rhs(l0, 1); gwu = sum(wvec * rhs(:, 1)); gww = sum(wvec * rhs(:, 2)); guw = rhs(l0, 2)
          ! exact conditional average
          b11 = 0; bh = 0; bhh = 0; b2 = 0
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
            if (pc <= 0.0_dp) cycle
            call dense_blocks(mdl, g, cf, z, l0, d11, dh, dhh, d2)
            b11 = b11 + pc * d11; bh = bh + pc * dh; bhh = bhh + pc * dhh; b2 = b2 + pc * d2
          end do
          eb(1) = max(eb(1), abs(guu - b11))
          eb(2) = max(eb(2), abs(gc * guu + gwu - bh))
          eb(3) = max(eb(3), abs(gc**2 * guu + 2.0_dp * gc * gwu + gww - bhh))
          esw = max(esw, abs(gwu - guw))
          ! two chains, full re-orthogonalisation
          eta = sqrt(real(sum(wvec * wvec), dp))
          call chain(op, et, uu, qt, uvec, wvec, nd, .true., cuu, cwu)
          call chain(op, et, uu, qt, wvec / eta, uvec, nd, .true., cww, cuw)
          cww = eta**2 * cww; cuw = eta * cuw
          ech = max(ech, abs(cuu - guu), abs(cwu - gwu), abs(cww - gww), abs(cuw - cwu))
          ! regression u = 0, p = 0, Qt -> -q
          do jc = 1, nd
            ej = (0.0_dp, 0.0_dp); ej(jc) = (1.0_dp, 0.0_dp)
            call kprime(op, ej, et, 0.0_dp * uu, mq, afac(:, jc))
          end do
          call fill_decoupled(afac, dec, ndec)
          rhs(:, 1) = uvec
          call zgesv(nd, 1, afac, nd, ipiv, rhs(:, 1), nd, inf)
          ereg = max(ereg, abs(rhs(l0, 1) - b2))
          ! canonical form
          do jc = 1, nd
            call random_number(u); call random_number(v); cvec(jc) = cmplx(u - 0.5_dp, v - 0.5_dp, dp)
          end do
          call kcanon(op, cen, cvec, z, rv)
          do jc = 1, nd; if (dec(jc)) cvec(jc) = (0.0_dp, 0.0_dp); end do      ! amat carries 1 on decoupled states
          call kcanon(op, cen, cvec, z, rv)
          ecan = max(ecan, maxval(abs(rv - matmul(amat, cvec))))
        end do
      end do
    end do
    print '(a,es10.3)', '  ring SRO: max |K''_aug - K''_aug^T|                 ', esym
    print '(a,3es10.3)', '  ring SRO: dense vs conditional average G11, hG, hGh', eb
    print '(a,es10.3,a,es10.3)', '  ring SRO: chains vs dense ', ech, '   |G_wu - G_uw| ', esw
    print '(a,es10.3,a,es10.3)', '  ring SRO: regression (asr_edc) ', ereg, '   canonical form ', ecan
    call check_close(esym, 0.0_dp, 1.0e-13_dp, 'ring SRO: K''_aug complex symmetric (N species, classes, conditional shell)')
    call check_close(maxval(eb), 0.0_dp, 1.0e-10_dp, 'ring SRO: G11, (hG)11, (hGh)11 == exact conditional configuration average')
    call check_close(esw, 0.0_dp, 1.0e-10_dp, 'ring SRO: G_wu == G_uw')
    call check_close(ech, 0.0_dp, 1.0e-8_dp, 'ring SRO: two Lanczos chains with projection == dense solve')
    call check_close(ereg, 0.0_dp, 1.0e-10_dp, 'ring SRO: u = 0, p = 0, Qt -> -q == <[E - H2]^-1_11> (asr_edc operator)')
    call check_close(ecan, 0.0_dp, 1.0e-12_dp, 'ring SRO: canonical form == factored product')
  end subroutine ring_sro

  !> fcc ternary {2,3,4} with a conditional first-shell class: D4h-reduced == unreduced for all four blocks
  subroutine fcc_reduced()
    integer, parameter :: nlev = 6, maxlev = 3
    real(dp), parameter :: cs(3, 3) = reshape([0.0_dp, 0.8_dp, 0.2_dp, 0.4_dp, 0.4_dp, 0.2_dp, 0.3_dp, 0.3_dp, 0.4_dp], [3, 3])
    type(graph_t) :: g
    type(orb_alloy_t) :: m0
    type(class_alloy_t) :: mdl
    type(augcls_op) :: opr, opu
    integer, allocatable :: coord(:, :), cls(:), cls0(:), perm(:, :), mats(:, :, :), pid(:, :)
    real(dp), allocatable :: d(:, :, :), did(:, :, :)
    logical, allocatable :: fl(:)
    complex(dp), allocatable :: ur(:), wr(:), uun(:), wun(:)
    complex(dp) :: z, et(9, 4), uu(9, 4), qt(9, 4), mq(9, 4), g1(4), g2(4)
    real(dp) :: err, etar, etau, gc, eeta, gmax
    integer :: i, j, cen, io, kz, r2, r2min, nsh, l0
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
    call class_ops(coord, 8, cls, perm, mats)
    call orbital_rep(mats, d)
    allocate (pid(g%n, 1), did(9, 9, 1)); pid(:, 1) = [(i, i=1, g%n)]; did = 0.0_dp
    do i = 1, 9; did(i, i, 1) = 1.0_dp; end do
    err = 0.0_dp; eeta = 0.0_dp; gmax = 0.0_dp
    do j = 1, 3
      cen = mdl%spc(j, 1)
      call set_hprime(mdl, cen)
      call build_augcls_operator(g, mdl, cls, perm, d, 1, maxlev, opr, maxrank=1, flippable=fl)
      call build_augcls_operator(g, mdl, cls, pid, did, 1, maxlev, opu, maxrank=1, flippable=fl)
      call augcls_freeze_site(opr, 1, cen); call augcls_freeze_site(opu, 1, cen)
      call augcls_set_class_conc(opr, 2, cs(:, j)); call augcls_set_class_conc(opu, 2, cs(:, j))
      if (j == 1) print '(a,i0,a,i0)', '  fcc: nodes unreduced ', opu%nn, '  reduced ', opr%nn
      if (allocated(ur)) deallocate (ur, wr, uun, wun)
      allocate (ur(opr%dim()), wr(opr%dim()), uun(opu%dim()), wun(opu%dim()))
      do io = 1, 4
        l0 = l0s(io)
        call augcls_set_start(opr, l0, d); call augcls_set_start(opu, l0, did)
        gc = csp(l0, cen)
        ur = (0.0_dp, 0.0_dp); ur(l0) = (1.0_dp, 0.0_dp); uun = (0.0_dp, 0.0_dp); uun(l0) = (1.0_dp, 0.0_dp)
        call hprime(opr, ur, wr); wr = wr - gc * ur
        call hprime(opu, uun, wun); wun = wun - gc * uun
        etar = sqrt(real(sum(wr * wr), dp)); etau = sqrt(real(sum(wun * wun), dp))
        eeta = max(eeta, abs(etar - etau))
        do kz = 1, 2
          z = cmplx(-0.2_dp + 0.3_dp * kz, 0.15_dp, dp)
          call coefs(cen, z, et, uu, qt, mq)
          call chain(opr, et, uu, qt, ur, wr, nlev, .false., g1(1), g1(2))
          call chain(opr, et, uu, qt, wr / etar, ur, nlev, .false., g1(3), g1(4))
          call chain(opu, et, uu, qt, uun, wun, nlev, .false., g2(1), g2(2))
          call chain(opu, et, uu, qt, wun / etau, uun, nlev, .false., g2(3), g2(4))
          g1(3) = etar**2 * g1(3); g1(4) = etar * g1(4); g2(3) = etau**2 * g2(3); g2(4) = etau * g2(4)
          err = max(err, maxval(abs(g1 - g2) / abs(g2)))
          gmax = max(gmax, maxval(abs(g2)))
        end do
      end do
    end do
    print '(a,es9.2,a,es9.2,a,es9.2)', '  fcc: max relative |reduced - unreduced| (G_uu, G_wu, G_ww, G_uw) ', err, &
      '   max |G| ', gmax, '   |eta_r - eta_u| ', eeta
    call check_close(eeta, 0.0_dp, 1.0e-12_dp, 'fcc: |S''u| equal in reduced and unreduced coordinates')
    call check_close(err, 0.0_dp, 1.0e-9_dp, 'fcc SRO: D4h-reduced == unreduced for G_uu, G_wu, G_ww, G_uw, relative (w in orbit coordinates)')
  end subroutine fcc_reduced
end program test_asr_ed3_cls
