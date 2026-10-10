!> n-component augmented space (asr_augn) vs brute-force enumeration of all ns^N configurations on a ring of N sites with nl orbitals
!> per site and H = E_nu + h - h o h (random c, d, o, E_nu per species, h = c + d S d):
!>  (1) ns = 3: Krylov space exhausted => recursion == exact configuration average of the resolvent (all sites random);
!>  (2) ns = 3 with the start site frozen to species 2 (central-species recursion): == average over the other sites, site 1 fixed;
!>  (3) ns = 2: chain coefficients identical to the binary operator asr_augorb;
!>  (4) the complex-parameter application augn_rand_c reproduces the real one.
program test_augn
  use asr_kinds, only: dp
  use asr_augspace, only: graph_t
  use asr_augorb, only: orb_alloy_t, augorb_op, build_augorb_operator
  use asr_augn
  use asr_recursion, only: lanczos
  use asr_greens, only: cf_green
  use asr_testing
  implicit none
  integer, parameter :: n = 4, nl = 2, nd = n * nl
  type(graph_t) :: g
  type(multi_alloy_t) :: mdl
  type(augn_op) :: op
  real(dp), allocatable :: a(:), b(:), v0(:), xv(:), yv(:), zv(:)
  complex(dp), allocatable :: cx(:), cw(:)
  complex(dp) :: z, gex, grec
  complex(dp), allocatable :: cpar(:, :)
  real(dp) :: u, err, err_f, err_b, pc, sym
  integer :: i, j, k, l, l2, s, m, kz, c, cfg(n), ncfg, nn
  logical :: fl(n)
  real(dp) :: fpar(nl, 3)

  print '(a)', 'test_augn'
  allocate (g%nbr_ptr(n + 1), g%nbr(2 * n)); g%n = n
  do i = 1, n
    g%nbr_ptr(i) = 2 * (i - 1) + 1
    g%nbr(2 * i - 1) = mod(i, n) + 1
    g%nbr(2 * i) = mod(i - 2 + n, n) + 1
  end do
  g%nbr_ptr(n + 1) = 2 * n + 1
  call random_seed()
  call make_model(3, [0.3_dp, 0.5_dp, 0.2_dp])

  ! ---- (1) all sites random, ns = 3
  call build_augn_operator(g, mdl, 1, 60, op)
  print '(a,i0,a,i0)', '  ns = 3 nodes: ', op%nn, '   (N 3^N = ', n * 3**n, ')'
  call check(op%nn == n * 3**n, 'all N ns^N nodes reached')
  allocate (xv(op%dim()), yv(op%dim()), zv(op%dim()), a(op%dim()), b(op%dim()), v0(op%dim()))
  call random_number(xv); call random_number(yv)
  call op%apply(xv, zv); call op%apply(yv, v0)
  sym = abs(dot_product(yv, zv) - dot_product(xv, v0))
  call check_close(sym, 0.0_dp, 1.0e-10_dp, 'matrix-free operator is symmetric')
  v0 = 0.0_dp; v0(1) = 1.0_dp
  call lanczos(op, v0, op%dim(), a, b, m, reorth=.true., tiny=1.0e-8_dp)
  err = 0.0_dp
  do kz = 1, 6
    z = cmplx(-1.5_dp + 0.5_dp * kz, 0.15_dp, dp)
    gex = (0.0_dp, 0.0_dp)
    ncfg = 3**n
    do c = 0, ncfg - 1
      call decode(c, cfg, pc)
      gex = gex + pc * dense_g(cfg, z, 1)
    end do
    grec = cf_green(a, b, m, z, (0.0_dp, 0.0_dp))
    err = max(err, abs(gex - grec))
  end do
  print '(a,i0,a,es9.2)', '  Krylov dimension ', m, '   max |recursion - exact average| = ', err
  call check_close(err, 0.0_dp, 1.0e-8_dp, 'ns = 3 augmented-space recursion == exact configuration average')

  ! ---- (4) complex parameters: X(e + i f) x = X(e) x + i X(f) x for real x
  allocate (cpar(nl, 3), cx(op%dim()), cw(op%dim()))
  call random_number(fpar)
  cpar = cmplx(mdl%e, fpar, dp)
  cx = cmplx(xv, 0.0_dp, dp)
  call augn_rand_c(op, cpar, cx, cw)
  call augn_rand(op, mdl%e, xv, zv)
  call augn_rand(op, fpar, xv, v0)
  err_b = max(maxval(abs(real(cw) - zv)), maxval(abs(aimag(cw) - v0)))
  call check_close(err_b, 0.0_dp, 1.0e-13_dp, 'complex-parameter application == real application (re/im)')

  ! ---- (2) frozen start site, species 2
  fl = .true.; fl(1) = .false.
  call build_augn_operator(g, mdl, 1, 60, op, flippable=fl)
  call augn_freeze_site(op, 1, 2)
  print '(a,i0,a,i0)', '  frozen start site: nodes ', op%nn, '   (N 3^(N-1) = ', n * 3**(n - 1)
  call check(op%nn == n * 3**(n - 1), 'frozen site: N ns^(N-1) nodes')
  deallocate (a, b, v0); allocate (a(op%dim()), b(op%dim()), v0(op%dim()))
  v0 = 0.0_dp; v0(1) = 1.0_dp
  call lanczos(op, v0, op%dim(), a, b, m, reorth=.true., tiny=1.0e-8_dp)
  err_f = 0.0_dp
  do kz = 1, 6
    z = cmplx(-1.5_dp + 0.5_dp * kz, 0.15_dp, dp)
    gex = (0.0_dp, 0.0_dp)
    do c = 0, 3**(n - 1) - 1
      call decode(c, cfg(2:n), pc)
      cfg(1) = 2
      gex = gex + pc * dense_g(cfg, z, 1)
    end do
    grec = cf_green(a, b, m, z, (0.0_dp, 0.0_dp))
    err_f = max(err_f, abs(gex - grec))
  end do
  print '(a,es9.2)', '  frozen central species: max |recursion - exact| = ', err_f
  call check_close(err_f, 0.0_dp, 1.0e-8_dp, 'frozen start site (species 2) == exact average over the other sites')

  ! ---- (3) ns = 2 vs the binary operator
  call make_model(2, [0.35_dp, 0.65_dp])
  call build_augn_operator(g, mdl, 1, 8, op)
  block
    type(orb_alloy_t) :: bm
    type(augorb_op) :: bop
    real(dp), allocatable :: a2(:), b2(:), v2(:), v3(:)
    integer :: m2, m3
    bm%nl = nl; bm%x = mdl%conc(1)
    bm%c0 = mdl%c(:, 2); bm%c1 = mdl%c(:, 1) - mdl%c(:, 2)
    bm%d0 = mdl%d(:, 2); bm%d1 = mdl%d(:, 1) - mdl%d(:, 2)
    bm%o0 = mdl%o(:, 2); bm%o1 = mdl%o(:, 1) - mdl%o(:, 2)
    bm%e0 = mdl%e(:, 2); bm%e1 = mdl%e(:, 1) - mdl%e(:, 2)
    bm%block = mdl%block; bm%block0 = mdl%block0
    call build_augorb_operator(g, bm, 1, 8, bop)
    if (allocated(a)) deallocate (a, b)
    allocate (a2(8), b2(8), a(8), b(8), v2(bop%dim()), v3(op%dim()))
    v2 = 0.0_dp; v2(1) = 1.0_dp; v3 = 0.0_dp; v3(1) = 1.0_dp
    call lanczos(bop, v2, 8, a2, b2, m2)
    call lanczos(op, v3, 8, a, b, m3)
    err = max(maxval(abs(a(1:8) - a2(1:8))), maxval(abs(b(1:7) - b2(1:7))))
    print '(a,i0,a,i0,a,es9.2)', '  ns = 2: nodes ', op%nn, ' vs binary ', bop%nn, '   max |coefficient difference| = ', err
    call check(op%nn == bop%nn, 'ns = 2: same node count as asr_augorb')
    call check_close(err, 0.0_dp, 1.0e-10_dp, 'ns = 2: chain coefficients == asr_augorb')
  end block
  call report()
contains
  subroutine make_model(nsp, x)
    integer, intent(in) :: nsp
    real(dp), intent(in) :: x(:)
    integer :: ll, ss, l1, l22
    mdl%nl = nl; mdl%ns = nsp
    if (allocated(mdl%conc)) deallocate (mdl%conc, mdl%c, mdl%d, mdl%o, mdl%e)
    if (allocated(mdl%block)) deallocate (mdl%block, mdl%block0)
    allocate (mdl%conc(nsp), mdl%c(nl, nsp), mdl%d(nl, nsp), mdl%o(nl, nsp), mdl%e(nl, nsp))
    mdl%conc = x
    do ss = 1, nsp
      do ll = 1, nl
        call random_number(u); mdl%c(ll, ss) = 0.8_dp * (u - 0.5_dp)
        call random_number(u); mdl%d(ll, ss) = 0.3_dp + 0.4_dp * u
        call random_number(u); mdl%o(ll, ss) = -0.2_dp * u
        call random_number(u); mdl%e(ll, ss) = 0.4_dp * (u - 0.5_dp)
      end do
    end do
    allocate (mdl%block(nl, nl, 2), mdl%block0(nl, nl))
    do l1 = 1, nl; do l22 = 1, nl
        call random_number(u); mdl%block(l1, l22, 1) = 0.8_dp * (u - 0.5_dp)
    end do; end do
    mdl%block(:, :, 2) = transpose(mdl%block(:, :, 1))
    do l1 = 1, nl; do l22 = l1, nl
        call random_number(u); mdl%block0(l1, l22) = 0.7_dp * (u - 0.5_dp); mdl%block0(l22, l1) = mdl%block0(l1, l22)
    end do; end do
  end subroutine make_model

  !> mixed-radix decode of c into species 1..ns per site (sites first_site..n of cfg) and its probability
  subroutine decode(c0, cf, prob)
    integer, intent(in) :: c0
    integer, intent(inout) :: cf(:)
    real(dp), intent(out) :: prob
    integer :: q, ii
    q = c0; prob = 1.0_dp
    do ii = 1, size(cf)
      cf(ii) = mod(q, mdl%ns) + 1; q = q / mdl%ns
      prob = prob * mdl%conc(cf(ii))
    end do
  end subroutine decode

  !> dense G_{(site0, orbital 1), same}(z) for the configuration cfg (species per site): (z - H)^-1, H = E_nu + h - h o h
  complex(dp) function dense_g(cf, zz, site0)
    integer, intent(in) :: cf(n), site0
    complex(dp), intent(in) :: zz
    real(dp) :: sm(nd, nd), hm(nd, nd), om(nd), dm(nd), cm(nd), em(nd), hh(nd, nd), ev(nd), wk2(300), vec(nd)
    complex(dp) :: gg
    integer :: ii, jj, kk, ll, ll2, rr, info, sp
    external :: dsyev
    sm = 0.0_dp
    do ii = 1, n
      do ll = 1, nl; do ll2 = 1, nl
          sm(nl * (ii - 1) + ll, nl * (ii - 1) + ll2) = mdl%block0(ll, ll2)
      end do; end do
      do kk = 1, 2
        jj = g%nbr(g%nbr_ptr(ii) + kk - 1)
        do ll = 1, nl; do ll2 = 1, nl
            sm(nl * (ii - 1) + ll, nl * (jj - 1) + ll2) = mdl%block(ll, ll2, kk)
        end do; end do
      end do
      sp = cf(ii)
      do ll = 1, nl
        rr = nl * (ii - 1) + ll
        dm(rr) = mdl%d(ll, sp); cm(rr) = mdl%c(ll, sp); om(rr) = mdl%o(ll, sp); em(rr) = mdl%e(ll, sp)
      end do
    end do
    hm = 0.0_dp
    do ii = 1, nd
      do jj = 1, nd; hm(ii, jj) = dm(ii) * sm(ii, jj) * dm(jj); end do
      hm(ii, ii) = hm(ii, ii) + cm(ii)
    end do
    do ii = 1, nd; do jj = 1, nd
        hh(ii, jj) = hm(ii, jj) - sum(hm(ii, :) * om * hm(:, jj))
    end do; end do
    do ii = 1, nd; hh(ii, ii) = hh(ii, ii) + em(ii); end do
    hh = 0.5_dp * (hh + transpose(hh))
    call dsyev('V', 'U', nd, hh, nd, ev, wk2, size(wk2), info)
    gg = (0.0_dp, 0.0_dp)
    rr = nl * (site0 - 1) + 1
    do jj = 1, nd
      gg = gg + hh(rr, jj)**2 / (zz - ev(jj))
    end do
    dense_g = gg
  end function dense_g

end program test_augn
