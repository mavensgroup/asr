!> Multi-orbital TB-LMTO-form augmented space (H = E + h - h o h with random c, d, o; and the third-order
!> H + h o h o h - (h p h h + h h p h)/2 with random p) vs brute-force
!> enumeration of all 2^N configurations on a ring of 6 sites with 2 orbitals per site.
!> Krylov space exhausted => recursion must equal the exact configuration average.
program test_augorb
  use asr_kinds, only: dp
  use asr_augspace, only: graph_t
  use asr_augorb
  use asr_recursion, only: lanczos
  use asr_greens, only: cf_green
  use asr_testing
  implicit none
  integer, parameter :: n = 6, nl = 2
  type(graph_t) :: g
  type(orb_alloy_t) :: mdl
  type(augorb_op) :: op
  real(dp) :: hmat(n * nl, n * nl), sm(n * nl, n * nl), dm(n * nl), cm(n * nl), om(n * nl, n * nl), h(n * nl, n * nl)
  real(dp) :: hh(n * nl, n * nl), ev(n * nl), work(200), pc, err, sym, pm(n * nl, n * nl), hp(n * nl, n * nl)
  integer :: ord
  real(dp), allocatable :: a(:), b(:), v0(:), xv(:), yv(:), zv(:), xr(:), yr(:)
  complex(dp) :: z, gex, grec
  integer :: i, j, k, c, info, m, nd, kz, l, l2, rr, s2, r1
  logical :: isa(n)
  external :: dsyev

  print '(a)', 'test_augorb'
  ! ring: slot 1 -> R+1, slot 2 -> R-1
  allocate (g%nbr_ptr(n + 1), g%nbr(2 * n)); g%n = n
  do i = 1, n
    g%nbr_ptr(i) = 2 * (i - 1) + 1
    g%nbr(2 * i - 1) = mod(i, n) + 1
    g%nbr(2 * i) = mod(i - 2 + n, n) + 1
  end do
  g%nbr_ptr(n + 1) = 2 * n + 1
  mdl%nl = nl; mdl%x = 0.35_dp
  mdl%c0 = [-0.3_dp, 0.4_dp]; mdl%c1 = [0.5_dp, -0.25_dp]
  mdl%d0 = [0.6_dp, 0.45_dp]; mdl%d1 = [0.2_dp, 0.3_dp]
  mdl%o0 = [-0.15_dp, -0.1_dp]; mdl%o1 = [0.05_dp, 0.08_dp]
  mdl%e0 = [0.1_dp, -0.2_dp]; mdl%e1 = [0.07_dp, -0.04_dp]
  mdl%p0 = [0.3_dp, 0.12_dp]; mdl%p1 = [-0.1_dp, 0.2_dp]
  allocate (mdl%block(nl, nl, 2))
  mdl%block(:, :, 1) = reshape([0.9_dp, -0.3_dp, 0.4_dp, 0.7_dp], [2, 2])
  mdl%block(:, :, 2) = transpose(mdl%block(:, :, 1))          ! S_{R,R-1} = S_{R-1,R}^T (Hermitian)
  allocate (mdl%block0(nl, nl)); mdl%block0 = reshape([0.8_dp, 0.15_dp, 0.15_dp, -0.5_dp], [2, 2])   ! on-site S_RR

  do ord = 2, 3
  print '(a,i0)', '  Hamiltonian order ', ord
  call build_augorb_operator(g, mdl, 1, 40, op, order=ord)
  nd = op%dim()
  print '(a,i0,a,i0)', '    augmented nodes: ', op%nn, '   (N 2^N = ', n * 2**n, ')'
  call check(op%nn == n * 2**n, 'all N 2^N nodes reached')
  ! matrix-free operator must be symmetric: <y|Hx> = <x|Hy>
  allocate (xv(nd), yv(nd), zv(nd), xr(nd), yr(nd))
  call random_number(xv); call random_number(yv)
  call op%apply(xv, zv); call op%apply(yv, xr)
  sym = abs(dot_product(yv, zv) - dot_product(xv, xr))
  call check_close(sym, 0.0_dp, 1.0e-10_dp, 'matrix-free operator is symmetric')

  allocate (a(nd), b(nd), v0(nd)); v0 = 0.0_dp; v0(1) = 1.0_dp       ! node 1 = (site 1, {}), orbital 1
  call lanczos(op, v0, nd, a, b, m, reorth=.true., tiny=1.0e-8_dp)
  print '(a,i0)', '    Krylov dimension used: ', m

  err = 0.0_dp
  do kz = 1, 10
    z = cmplx(-2.0_dp + 0.4_dp * kz, 0.1_dp, dp)
    gex = (0.0_dp, 0.0_dp)
    do c = 0, 2**n - 1
      pc = 1.0_dp
      do i = 1, n
        isa(i) = btest(c, i - 1)
        if (isa(i)) then; pc = pc * mdl%x; else; pc = pc * (1.0_dp - mdl%x); end if
      end do
      sm = 0.0_dp; dm = 0.0_dp; cm = 0.0_dp; om = 0.0_dp; pm = 0.0_dp
      do i = 1, n
        do l = 1, nl
          rr = nl * (i - 1) + l
          dm(rr) = mdl%d0(l) + merge(mdl%d1(l), 0.0_dp, isa(i))
          cm(rr) = mdl%c0(l) + merge(mdl%c1(l), 0.0_dp, isa(i))
          om(rr, rr) = mdl%o0(l) + merge(mdl%o1(l), 0.0_dp, isa(i))
          pm(rr, rr) = mdl%p0(l) + merge(mdl%p1(l), 0.0_dp, isa(i))
        end do
        do l = 1, nl; do l2 = 1, nl
            sm(nl * (i - 1) + l, nl * (i - 1) + l2) = mdl%block0(l, l2)
        end do; end do
        do k = 1, 2
          j = g%nbr(g%nbr_ptr(i) + k - 1)
          do l = 1, nl; do l2 = 1, nl
              sm(nl * (i - 1) + l, nl * (j - 1) + l2) = mdl%block(l, l2, k)
          end do; end do
        end do
      end do
      h = 0.0_dp
      do i = 1, n * nl
        do j = 1, n * nl
          h(i, j) = dm(i) * sm(i, j) * dm(j)
        end do
        h(i, i) = h(i, i) + cm(i)
      end do
      hh = h - matmul(h, matmul(om, h))
      if (ord == 3) then
        hp = matmul(h, matmul(pm, h))
        hh = hh + matmul(h, matmul(om, matmul(h, matmul(om, h)))) - 0.5_dp * (matmul(hp, h) + matmul(h, hp))
      end if
      do i = 1, n
        do l = 1, nl; hh(nl * (i - 1) + l, nl * (i - 1) + l) = hh(nl * (i - 1) + l, nl * (i - 1) + l) + mdl%e0(l) + merge(mdl%e1(l), 0.0_dp, isa(i)); end do
      end do
      hmat = 0.5_dp * (hh + transpose(hh))
      call dsyev('V', 'U', n * nl, hmat, n * nl, ev, work, size(work), info)
      do j = 1, n * nl
        gex = gex + pc * hmat(1, j)**2 / (z - ev(j))
      end do
    end do
    grec = cf_green(a, b, m, z, (0.0_dp, 0.0_dp))
    err = max(err, abs(gex - grec))
  end do
  call check_close(err, 0.0_dp, 1.0e-9_dp, 'TB-LMTO augmented-space recursion == exact configuration average')
  deallocate (a, b, v0, xv, yv, zv, xr, yr)
  end do
  call report()
end program test_augorb
