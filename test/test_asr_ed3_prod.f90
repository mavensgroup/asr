!> Production routine of docs/DERIVATION_ed3.md (asr_ed3 on asr_augcls), run before any desktop gate:
!>  (P1) augcls_hp_c == augcls_apply with o = e = 0 (one application of S instead of two).
!>  (P2) ring of 5, four species, classes, conditional first shell: ed3_blocks == exact conditional configuration average of
!>       G11, (hG)11, (hGh)11 of K = E O - H (bands_core form), for the parameters of ed3_setup.
!>  (P3) ordered ring: Fermi-contour integrals of the blocks == occupied sums over the eigenpairs of the dense pencil (H, O), and the
!>       channel charges xm0 + 2 o xm1 + p xm2 summed over the cell == number of occupied states (section 8, electron count).
!>  (P4) terminated chains (hull tail), bcc CuZn parameters with a conditional shell, reduced operator: against a long chain,
!>       the error of the charge Q falls with the steps; -Im of the block matrix (q0, q1; q1, q2) positive semidefinite on the main
!>       contour line and above (Im z >= 0.1).  Closer to the axis the counts are printed: the cross term q1 is first order there.
program test_asr_ed3_prod
  use asr_kinds, only: dp, pi
  use asr_augspace, only: graph_t
  use asr_augcls
  use asr_augorb, only: orb_alloy_t
  use asr_lmto_alloy, only: build_lmto_alloy
  use asr_symmetry
  use asr_lmto_ham, only: potpar_t, to_alpha, o_alpha
  use asr_ed3
  use asr_endep, only: ed_tail, ed_depth_tol
  use asr_contour, only: contour_nodes_fermi, fermi_nodes_count
  use asr_testing
  implicit none
  integer, parameter :: orb_l(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2], iorb(4) = [1, 3, 5, 7]
  type(potpar_t) :: pp(4), qq(4)
  real(dp) :: alpha(0:2), oo(0:2, 4)
  real(dp) :: dsp(9, 4), csp(9, 4), esp(9, 4), osp(9, 4), psp(9, 4)
  integer :: l, sp
  external :: zgesv, dsygv

  print '(a)', 'test_asr_ed3_prod'
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
  ed_depth_tol = 0.0_dp
  ed_tail = 'none'
  call ring_sro()
  call ring_ordered()
  ed_tail = 'hull'
  call terminated()
  call report()

contains

  !> dense real matrices of one ring configuration: h, O, H of bands_core
  subroutine dense_ham(mdl, g, cf, hm, om, hmat)
    type(class_alloy_t), intent(in) :: mdl
    type(graph_t), intent(in) :: g
    integer, intent(in) :: cf(:)
    real(dp), intent(out) :: hm(:, :), om(:, :), hmat(:, :)
    integer :: n, nd, ii, jj, kk, ll, kq
    real(dp), allocatable :: s(:, :), dl(:), cm(:), en(:), ov(:), pv(:)
    n = g%n; nd = 9 * n
    allocate (s(nd, nd), dl(nd), cm(nd), en(nd), ov(nd), pv(nd))
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
      end do
      om(jj, jj) = om(jj, jj) + 1.0_dp
      hmat(jj, jj) = hmat(jj, jj) + en(jj)
    end do
  end subroutine dense_ham

  !> G11, (hG)11, (hGh)11 of K = E O - H for orbital l0 of site 1
  subroutine dense_blocks(mdl, g, cf, zz, l0, g11, gh, ghh)
    type(class_alloy_t), intent(in) :: mdl
    type(graph_t), intent(in) :: g
    integer, intent(in) :: cf(:), l0
    complex(dp), intent(in) :: zz
    complex(dp), intent(out) :: g11, gh, ghh
    integer :: nd, inf
    real(dp), allocatable :: hm(:, :), om(:, :), hmat(:, :)
    complex(dp), allocatable :: ka(:, :), r(:, :)
    integer, allocatable :: ip(:)
    nd = 9 * g%n
    allocate (hm(nd, nd), om(nd, nd), hmat(nd, nd), ka(nd, nd), r(nd, 2), ip(nd))
    call dense_ham(mdl, g, cf, hm, om, hmat)
    ka = zz * om - hmat
    r = (0.0_dp, 0.0_dp); r(l0, 1) = (1.0_dp, 0.0_dp); r(:, 2) = hm(:, l0)
    call zgesv(nd, 2, ka, nd, ip, r, nd, inf)
    if (inf /= 0) error stop 'zgesv (dense K)'
    g11 = r(l0, 1); gh = sum(hm(l0, :) * r(:, 1)); ghh = sum(hm(:, l0) * r(:, 2))
  end subroutine dense_blocks

  !> the three blocks of one site orbital k from the dense matrices h, O, H of an ordered cell
  subroutine dense_site(hm, om, hmat, k, zz, dk)
    real(dp), intent(in) :: hm(:, :), om(:, :), hmat(:, :)
    integer, intent(in) :: k
    complex(dp), intent(in) :: zz
    complex(dp), intent(out) :: dk(0:2)
    complex(dp) :: ka(size(hm, 1), size(hm, 1)), r(size(hm, 1), 2)
    integer :: ip(size(hm, 1)), nd, inf
    nd = size(hm, 1)
    ka = zz * om - hmat
    r = (0.0_dp, 0.0_dp); r(k, 1) = (1.0_dp, 0.0_dp); r(:, 2) = hm(:, k)
    call zgesv(nd, 2, ka, nd, ip, r, nd, inf)
    if (inf /= 0) error stop 'zgesv (dense_site)'
    dk(0) = r(k, 1); dk(1) = sum(hm(k, :) * r(:, 1)); dk(2) = sum(hm(:, k) * r(:, 2))
  end subroutine dense_site

  subroutine ring_graph(n, g)
    integer, intent(in) :: n
    type(graph_t), intent(out) :: g
    integer :: i
    allocate (g%nbr_ptr(n + 1), g%nbr(2 * n)); g%n = n
    do i = 1, n
      g%nbr_ptr(i) = 2 * (i - 1) + 1
      g%nbr(2 * i - 1) = mod(i, n) + 1
      g%nbr(2 * i) = mod(i - 2 + n, n) + 1
    end do
    g%nbr_ptr(n + 1) = 2 * n + 1
  end subroutine ring_graph

  subroutine random_blocks(mdl)
    type(class_alloy_t), intent(inout) :: mdl
    real(dp) :: u
    integer :: i, j
    allocate (mdl%block(9, 9, 2), mdl%block0(9, 9))
    do i = 1, 9; do j = 1, 9
      call random_number(u); mdl%block(i, j, 1) = 0.5_dp * (u - 0.5_dp)
    end do; end do
    mdl%block(:, :, 2) = transpose(mdl%block(:, :, 1))
    do i = 1, 9; do j = i, 9
      call random_number(u); mdl%block0(i, j) = 0.6_dp * (u - 0.5_dp); mdl%block0(j, i) = mdl%block0(i, j)
    end do; end do
  end subroutine random_blocks

  subroutine ring_sro()
    integer, parameter :: n = 5, rcls(n) = [3, 4, 1, 2, 5], nlev = 400
    real(dp), parameter :: c4(3, 3) = reshape([0.0_dp, 0.7_dp, 0.3_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.2_dp, 0.6_dp, 0.2_dp], [3, 3])
    real(dp), parameter :: c5(2, 3) = reshape([0.1_dp, 0.9_dp, 1.0_dp, 0.0_dp, 0.55_dp, 0.45_dp], [2, 3])
    type(graph_t) :: g
    type(class_alloy_t) :: mdl
    type(augcls_op) :: op
    type(ed3_t) :: ed
    type(ed3_work_t) :: wk
    integer :: perm(n, 1), i, j, cen, kz, cf(n), q, c, ncfg, l0, il
    real(dp) :: dmat(9, 9, 1), pc, eb(3), eapply, easym, u
    real(dp), allocatable :: xr(:), yr(:)
    complex(dp), allocatable :: xc(:), yc(:)
    complex(dp) :: z, bk(0:3), b11, bh, bhh, d11, dh, dhh
    logical :: fl(n)
    call ring_graph(n, g)
    mdl%nl = 9; mdl%nsp = 4; mdl%ncls = 5
    allocate (mdl%nsc(5), mdl%spc(3, 5), mdl%conc(3, 5)); mdl%spc = 0; mdl%conc = 0.0_dp
    mdl%nsc = [1, 2, 3, 3, 2]
    mdl%spc(1:1, 1) = [1];       mdl%conc(1:1, 1) = [1.0_dp]
    mdl%spc(1:2, 2) = [2, 3];    mdl%conc(1:2, 2) = [0.4_dp, 0.6_dp]
    mdl%spc(1:3, 3) = [1, 3, 4]; mdl%conc(1:3, 3) = [0.3_dp, 0.5_dp, 0.2_dp]
    mdl%spc(:, 4) = mdl%spc(:, 3); mdl%conc(:, 4) = mdl%conc(:, 3)
    mdl%spc(:, 5) = mdl%spc(:, 2); mdl%conc(:, 5) = mdl%conc(:, 2)
    call random_blocks(mdl)
    perm(:, 1) = [(i, i=1, n)]
    dmat = 0.0_dp; do l = 1, 9; dmat(l, l, 1) = 1.0_dp; end do
    fl = .true.; fl(1) = .false.
    eb = 0.0_dp; eapply = 0.0_dp; easym = 0.0_dp
    do j = 1, 3
      cen = mdl%spc(j, 3)
      call ed3_setup(mdl, pp, alpha, cen, ed)
      call build_augcls_operator(g, mdl, rcls, perm, dmat, 1, 80, op, flippable=fl)
      call augcls_freeze_site(op, 1, cen)
      call augcls_set_class_conc(op, 4, c4(:, j)); call augcls_set_class_conc(op, 5, c5(:, j))
      if (j == 1) print '(a,i0,a)', '  ring SRO: ', op%nn, ' augmented nodes'
      call ed3_work_init(wk, op%dim())
      if (allocated(xr)) deallocate (xr, yr, xc, yc)
      allocate (xr(op%dim()), yr(op%dim()), xc(op%dim()), yc(op%dim()))
      do i = 1, op%dim(); call random_number(u); xr(i) = u - 0.5_dp; end do
      call op%apply(xr, yr)
      xc = cmplx(xr, 0.7_dp * xr, dp)
      call augcls_hp_c(op, xc, yc)
      eapply = max(eapply, maxval(abs(yc - cmplx(yr, 0.7_dp * yr, dp))))
      do il = 1, 2
        l0 = merge(1, 5, il == 1)
        call augcls_set_start(op, l0, dmat)
        call ed3_start(op, wk, ed, l0)
        do kz = 1, 3
          z = cmplx(-0.3_dp + 0.3_dp * kz, 0.3_dp, dp)
          call ed3_blocks(op, wk, ed, z, l0, nlev, bk)
          b11 = 0; bh = 0; bhh = 0
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
            call dense_blocks(mdl, g, cf, z, l0, d11, dh, dhh)
            b11 = b11 + pc * d11; bh = bh + pc * dh; bhh = bhh + pc * dhh
          end do
          eb(1) = max(eb(1), abs(bk(0) - b11)); eb(2) = max(eb(2), abs(bk(1) - bh)); eb(3) = max(eb(3), abs(bk(2) - bhh))
          easym = max(easym, abs(bk(3)))
        end do
      end do
    end do
    print '(a,3es10.3,a,es10.3)', '  ring SRO: production vs conditional average G11, hG, hGh', eb, '   |dq1| ', easym
    call check_close(eapply, 0.0_dp, 1.0e-13_dp, 'augcls_hp_c == augcls_apply with o = e = 0')
    call check_close(maxval(eb), 0.0_dp, 1.0e-9_dp, 'ring SRO: ed3_blocks == exact conditional configuration average')
    call check_close(easym, 0.0_dp, 1.0e-9_dp, 'ring SRO: the two chains agree on q1')
  end subroutine ring_sro

  !> ordered ring ABABAB, Fermi contour against the eigenpairs of the dense pencil, both sublattices
  subroutine ring_ordered()
    integer, parameter :: n = 4, nd = 9 * n, npl = 2, nmain = 64, nstep = 48, nlev = 400
    real(dp), parameter :: kt = 0.02_dp
    type(graph_t) :: g
    type(class_alloy_t) :: mdl
    type(augcls_op) :: op
    type(ed3_t) :: ed
    type(ed3_work_t) :: wk
    integer :: rcls(n), perm(n, 1), i, a, cen, cf(n), nq, j, lw, inf, iv, k
    real(dp) :: dmat(9, 9, 1), hm(nd, nd), om(nd, nd), hmat(nd, nd), ev(nd), work(8 * nd), cc(nd), zv(nd), fo(nd), mu, xs(0:2), xc(0:2)
    real(dp) :: emom, nocc, ncount, qs, xd(0:2), econt, echain, om0(nd, nd), hmat0(nd, nd)
    complex(dp) :: dk(0:2)
    complex(dp), allocatable :: zq(:), wq(:)
    complex(dp) :: bk(0:3)
    logical :: fl(n)
    call ring_graph(n, g)
    mdl%nl = 9; mdl%nsp = 2; mdl%ncls = 2
    allocate (mdl%nsc(2), mdl%spc(1, 2), mdl%conc(1, 2)); mdl%nsc = 1; mdl%spc(1, :) = [1, 2]; mdl%conc = 1.0_dp
    call random_blocks(mdl)
    perm(:, 1) = [(i, i=1, n)]
    dmat = 0.0_dp; do l = 1, 9; dmat(l, l, 1) = 1.0_dp; end do
    fl = .true.; fl(1) = .false.
    cf = [1, 2, 1, 2]
    call dense_ham(mdl, g, cf, hm, om, hmat)
    om0 = om; hmat0 = hmat
    lw = size(work)
    call dsygv(1, 'V', 'U', nd, hmat, nd, om, nd, ev, work, lw, inf)
    if (inf /= 0) error stop 'dsygv (ordered ring)'
    mu = 0.5_dp * (ev(nd / 2) + ev(nd / 2 + 1))
    do i = 1, nd; fo(i) = 1.0_dp / (exp((ev(i) - mu) / kt) + 1.0_dp); end do
    nocc = sum(fo)
    print '(a,2f9.4,a,f9.4,a,f8.3)', '  ordered ring: pencil spectrum ', ev(1), ev(nd), '  mu ', mu, '  N ', nocc
    nq = fermi_nodes_count(npl, nmain, nstep)
    allocate (zq(nq), wq(nq))
    call contour_nodes_fermi(ev(1) - 1.0_dp, mu, kt, npl, nmain, nstep, zq, wq)
    emom = 0.0_dp; ncount = 0.0_dp; econt = 0.0_dp; echain = 0.0_dp
    ed_depth_tol = 1.0e-13_dp
    do a = 1, 2                                         ! origin on sublattice a, central species a
      cen = a
      rcls = cf; if (a == 2) rcls = 3 - cf
      call ed3_setup(mdl, pp(1:2), alpha, cen, ed)
      call build_augcls_operator(g, mdl, rcls, perm, dmat, 1, 80, op, flippable=fl)
      call augcls_freeze_site(op, 1, cen)
      call ed3_work_init(wk, op%dim())
      do iv = 1, 9
        call augcls_set_start(op, iv, dmat)
        call ed3_start(op, wk, ed, iv)
        xc = 0.0_dp; xd = 0.0_dp
        do j = 1, nq
          call ed3_blocks(op, wk, ed, zq(j), iv, nlev, bk)
          do k = 0, 2; xc(k) = xc(k) - aimag(wq(j) * bk(k)) / pi; end do
          call dense_site(hm, om0, hmat0, 9 * (a - 1) + iv, zq(j), dk)
          do k = 0, 2; xd(k) = xd(k) - aimag(wq(j) * dk(k)) / pi; end do
          echain = max(echain, maxval(abs(bk(0:2) - dk)))
        end do
        ! the same three sums from the eigenpairs; orbital iv of a site of sublattice a
        xs = 0.0_dp
        do i = 1, nd
          zv = hmat(:, i)
          cc = matmul(hm, zv)
          k = 9 * (a - 1) + iv
          xs(0) = xs(0) + fo(i) * zv(k)**2
          xs(1) = xs(1) + fo(i) * zv(k) * cc(k)
          xs(2) = xs(2) + fo(i) * cc(k)**2
        end do
        emom = max(emom, maxval(abs(xc - xs)))
        econt = max(econt, maxval(abs(xd - xs)))
        qs = xc(0) + 2.0_dp * ed%o(iv, cen) * xc(1) + ed%p(iv, cen) * xc(2)
        ncount = ncount + 2.0_dp * qs                    ! two cells of two sites
      end do
    end do
    ed_depth_tol = 0.0_dp
    print '(a,es10.3,a,es10.3)', '  ordered ring: dense G on the contour vs eigenpairs ', econt, '   chain vs dense G, max over nodes ', echain
    print '(a,es10.3,a,2f12.7)', '  ordered ring: contour vs eigenpairs, max error ', emom, '   count: charges, occupied states ', ncount, nocc
    call check_close(emom, 0.0_dp, 1.0e-8_dp, 'ordered ring: contour blocks == occupied sums of the dense pencil')
    call check_close(ncount, nocc, 1.0e-7_dp, 'ordered ring: channel charges over the cell == number of occupied states')
  end subroutine ring_ordered

  subroutine terminated()
    integer, parameter :: nre = 5, ni = 4, ndep = 2, deps(ndep) = [24, 36], nref = 240
    real(dp), parameter :: etas(ni) = [0.01_dp, 0.03_dp, 0.1_dp, 0.3_dp], res(nre) = [-0.75_dp, -0.5_dp, -0.3_dp, -0.13_dp, 0.1_dp]
    real(dp), parameter :: cs(2, 2) = reshape([0.35_dp, 0.65_dp, 0.65_dp, 0.35_dp], [2, 2])
    type(graph_t) :: g
    type(orb_alloy_t) :: m0
    type(class_alloy_t) :: mdl
    type(augcls_op) :: op
    type(ed3_t) :: ed
    type(ed3_work_t) :: wk
    type(potpar_t) :: pz(2)
    integer, allocatable :: coord(:, :), cls(:), perm(:, :), mats(:, :, :)
    real(dp), allocatable :: d(:, :, :)
    logical, allocatable :: fl(:)
    complex(dp) :: z, bk(0:3), rf(0:3), qz, qr, q0z
    real(dp) :: oq, pq, eq(ni, ndep), sq(ni), a_, b_, c_
    integer :: i, j, cen, io, ir, ii, nd, r2, r2min, iv, nneg(ni, ndep), npsd(ni, ndep)
    pz(1)%enu = [-0.434080_dp, -0.414749_dp, -0.424591_dp]; pz(1)%c = [-0.489840_dp, 0.517506_dp, -0.374831_dp]
    pz(1)%sqrdel = [0.406986_dp, 0.396185_dp, 0.097083_dp]; pz(1)%p = [0.063280_dp, 0.028936_dp, 2.290907_dp]
    pz(1)%gamma = [0.417587_dp, 0.111352_dp, 0.000411_dp]
    pz(2)%enu = [-0.526805_dp, -0.410600_dp, -0.799547_dp]; pz(2)%c = [-0.635378_dp, 0.361782_dp, -0.717302_dp]
    pz(2)%sqrdel = [0.394422_dp, 0.382562_dp, 0.081569_dp]; pz(2)%p = [0.077548_dp, 0.035406_dp, 2.794140_dp]
    pz(2)%gamma = [0.412419_dp, 0.108093_dp, 0.005836_dp]
    pz(1)%alpha_rep = pz(1)%gamma; pz(2)%alpha_rep = pz(2)%gamma
    call build_lmto_alloy('bcc', 4, 1.05_dp, 0.7_dp, pz(1), pz(2), 0.5_dp, g, coord, m0, alpha_in=alpha)
    mdl%nl = 9; mdl%nsp = 2; mdl%ncls = 2
    allocate (mdl%nsc(2), mdl%spc(2, 2), mdl%conc(2, 2))
    mdl%nsc = 2; mdl%spc(:, 1) = [1, 2]; mdl%spc(:, 2) = [1, 2]; mdl%conc = 0.5_dp
    mdl%block = m0%block; mdl%block0 = m0%block0
    allocate (fl(g%n), cls(g%n)); fl = .true.; fl(1) = .false.; cls = 1
    r2min = huge(1)
    do i = 2, g%n; r2min = min(r2min, sum((modulo(coord(:, i) + 4, 8) - 4)**2)); end do
    do i = 2, g%n
      r2 = sum((modulo(coord(:, i) + 4, 8) - 4)**2)
      if (r2 == r2min) cls(i) = 2
    end do
    call class_ops(coord, 8, cls, perm, mats)
    call orbital_rep(mats, d)
    eq = 0.0_dp; sq = 0.0_dp; nneg = 0; npsd = 0
    do j = 1, 2
      cen = j
      call ed3_setup(mdl, pz, alpha, cen, ed)
      call build_augcls_operator(g, mdl, cls, perm, d, 1, 4, op, maxrank=1, flippable=fl)
      call augcls_freeze_site(op, 1, cen)
      call augcls_set_class_conc(op, 2, cs(:, j))
      call ed3_work_init(wk, op%dim())
      do io = 3, 4                                       ! d channels
        iv = iorb(io)
        call augcls_set_start(op, iv, d)
        call ed3_start(op, wk, ed, iv)
        oq = ed%o(iv, cen); pq = ed%p(iv, cen)
        do ir = 1, nre
          do ii = 2, ni
            z = cmplx(res(ir), etas(ii), dp)
            call ed3_blocks(op, wk, ed, z, iv, nref, rf)
            qr = rf(0) + 2.0_dp * oq * rf(1) + pq * rf(2)
            sq(ii) = sq(ii) + abs(qr)
            do nd = 1, ndep
              call ed3_blocks(op, wk, ed, z, iv, deps(nd), bk)
              qz = bk(0) + 2.0_dp * oq * bk(1) + pq * bk(2)
              eq(ii, nd) = eq(ii, nd) + abs(qz - qr)
            end do
          end do
          z = cmplx(res(ir), etas(1), dp)                ! near the axis: bounds only (Matsubara nodes)
          do ii = 1, ni
            z = cmplx(res(ir), etas(ii), dp)
            do nd = 1, ndep
              call ed3_blocks(op, wk, ed, z, iv, deps(nd), bk)
              q0z = bk(0) + 2.0_dp * oq * bk(1) + oq**2 * bk(2)
              a_ = -aimag(q0z); b_ = -aimag(bk(1) + oq * bk(2)); c_ = -aimag(bk(2))
              if (a_ <= 0.0_dp .or. c_ <= 0.0_dp) then
                nneg(ii, nd) = nneg(ii, nd) + 1
              else if (b_**2 > a_ * c_) then
                npsd(ii, nd) = npsd(ii, nd) + 1
              end if
            end do
          end do
        end do
      end do
    end do
    print '(a)', '  terminated: steps  Im z  mean|dQ|/mean|Q|  q0 or q2 not positive  not PSD  (of 20 energies)'
    do nd = 1, ndep
      do ii = 1, ni
        print '(a,i6,f7.2,es14.3,i12,i12)', '  terminated:', deps(nd), etas(ii), eq(ii, nd) / max(sq(ii), tiny(1.0_dp)), nneg(ii, nd), npsd(ii, nd)
      end do
    end do
    call check(all(eq(2:ni, 2) < eq(2:ni, 1)), 'terminated: the error of Q falls from 24 to 36 steps')
    call check(eq(3, 2) / sq(3) < 1.0e-4_dp .and. eq(4, 2) / sq(4) < 1.0e-8_dp, 'terminated: Q converged on the main line (Im z 0.1, 0.3) at 36 steps')
    call check(all(nneg(3:ni, :) == 0) .and. all(npsd(3:ni, :) == 0), 'terminated: block -Im (q0, q1; q1, q2) positive semidefinite at Im z >= 0.1')
    call check(all(nneg(:, :) == 0), 'terminated: q0 and q2 are diagonal elements, -Im positive at all Im z')
  end subroutine terminated
end program test_asr_ed3_prod
