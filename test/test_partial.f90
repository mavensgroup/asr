!> Species-resolved Green functions: <<n_R G_RR>> from the recursion (transition element) vs enumeration of all
!> configurations weighted by the occupation of the start site (ring of 6 sites, 2 orbitals, TB-LMTO form).
program test_partial
  use asr_kinds, only: dp
  use asr_augspace, only: graph_t
  use asr_augorb
  use asr_partial
  use asr_greens, only: cf_green
  use asr_testing
  implicit none
  integer, parameter :: n = 6, nl = 2
  type(graph_t) :: g
  type(orb_alloy_t) :: mdl
  type(augorb_op) :: op
  real(dp) :: hmat(n * nl, n * nl), sm(n * nl, n * nl), dm(n * nl), cm(n * nl), om(n * nl, n * nl), h(n * nl, n * nl)
  real(dp) :: hh(n * nl, n * nl), ev(n * nl), work(200), pc, err, errb
  real(dp), allocatable :: a(:), b(:), v0(:), u(:, :), c(:), v(:)
  complex(dp) :: z, gA_ex, gB_ex, gavg_ex, fA, gavg
  integer :: i, j, k, cfg, info, m, nd, kz, l, l2, rr, iorb
  logical :: isa(n)
  external :: dsyev

  print '(a)', 'test_partial'
  allocate (g%nbr_ptr(n + 1), g%nbr(2 * n)); g%n = n
  do i = 1, n
    g%nbr_ptr(i) = 2 * (i - 1) + 1; g%nbr(2 * i - 1) = mod(i, n) + 1; g%nbr(2 * i) = mod(i - 2 + n, n) + 1
  end do
  g%nbr_ptr(n + 1) = 2 * n + 1
  mdl%nl = nl; mdl%x = 0.35_dp
  mdl%c0 = [-0.3_dp, 0.4_dp]; mdl%c1 = [0.5_dp, -0.25_dp]
  mdl%d0 = [0.6_dp, 0.45_dp]; mdl%d1 = [0.2_dp, 0.3_dp]
  mdl%o0 = [-0.15_dp, -0.1_dp]; mdl%o1 = [0.05_dp, 0.08_dp]
  mdl%e0 = [0.1_dp, -0.2_dp]; mdl%e1 = [0.07_dp, -0.04_dp]
  allocate (mdl%block(nl, nl, 2), mdl%block0(nl, nl))
  mdl%block(:, :, 1) = reshape([0.9_dp, -0.3_dp, 0.4_dp, 0.7_dp], [2, 2])
  mdl%block(:, :, 2) = transpose(mdl%block(:, :, 1))
  mdl%block0 = reshape([0.8_dp, 0.15_dp, 0.15_dp, -0.5_dp], [2, 2])
  call build_augorb_operator(g, mdl, 1, 40, op)
  nd = op%dim()
  allocate (a(nd), b(nd), v0(nd), v(nd), c(nd)); v0 = 0.0_dp
  iorb = 2
  v0(iorb) = 1.0_dp                               ! node 1 = (site 1, {}), orbital iorb
  call lanczos_basis(op, v0, nd, a, b, u, m)
  ! v = M_R u1 = x u1 + s |R,{R}>   (orbital iorb)
  v = 0.0_dp
  v(iorb) = op%md(1)
  v(nl * (op%flip(1) - 1) + iorb) = v(nl * (op%flip(1) - 1) + iorb) + op%s
  do k = 1, m; c(k) = dot_product(u(:, k), v); end do

  err = 0.0_dp; errb = 0.0_dp
  do kz = 1, 8
    z = cmplx(-1.5_dp + 0.4_dp * kz, 0.12_dp, dp)
    gA_ex = (0.0_dp, 0.0_dp); gB_ex = (0.0_dp, 0.0_dp)
    do cfg = 0, 2**n - 1
      pc = 1.0_dp
      do i = 1, n
        isa(i) = btest(cfg, i - 1)
        if (isa(i)) then; pc = pc * mdl%x; else; pc = pc * (1.0_dp - mdl%x); end if
      end do
      sm = 0.0_dp; dm = 0.0_dp; cm = 0.0_dp; om = 0.0_dp
      do i = 1, n
        do l = 1, nl
          rr = nl * (i - 1) + l
          dm(rr) = mdl%d0(l) + merge(mdl%d1(l), 0.0_dp, isa(i))
          cm(rr) = mdl%c0(l) + merge(mdl%c1(l), 0.0_dp, isa(i))
          om(rr, rr) = mdl%o0(l) + merge(mdl%o1(l), 0.0_dp, isa(i))
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
        do j = 1, n * nl; h(i, j) = dm(i) * sm(i, j) * dm(j); end do
        h(i, i) = h(i, i) + cm(i)
      end do
      hh = h - matmul(h, matmul(om, h))
      do i = 1, n
        do l = 1, nl
          rr = nl * (i - 1) + l
          hh(rr, rr) = hh(rr, rr) + mdl%e0(l) + merge(mdl%e1(l), 0.0_dp, isa(i))
        end do
      end do
      hmat = 0.5_dp * (hh + transpose(hh))
      call dsyev('V', 'U', n * nl, hmat, n * nl, ev, work, size(work), info)
      do j = 1, n * nl
        if (isa(1)) then
          gA_ex = gA_ex + pc * hmat(iorb, j)**2 / (z - ev(j))
        else
          gB_ex = gB_ex + pc * hmat(iorb, j)**2 / (z - ev(j))
        end if
      end do
    end do
    gavg = cf_green(a, b, m, z, (0.0_dp, 0.0_dp))
    fA = transition_green(a, b, m, c, z, (0.0_dp, 0.0_dp))
    err = max(err, abs(fA - gA_ex))                       ! <<n_R G>> = x G_A
    errb = max(errb, abs((gavg - fA) - gB_ex))            ! <<(1-n_R) G>> = (1-x) G_B
    gavg_ex = gA_ex + gB_ex
    err = max(err, abs(gavg - gavg_ex) * 0.0_dp)
  end do
  call check_close(err, 0.0_dp, 1.0e-9_dp, 'A-resolved <<n_R G_RR>> from the transition element == enumeration')
  call check_close(errb, 0.0_dp, 1.0e-9_dp, 'B-resolved <<(1-n_R) G_RR>> == enumeration')
  call report()
end program test_partial
