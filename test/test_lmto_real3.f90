!> Third-order real-space Hamiltonian (lmto_real_op, order = 3) on a small B2 supercell vs dense matrix algebra:
!>   H3 = E + h - h o h + h o h o h - (h p h h + h h p h)/2,
!> with the dense h obtained from the same operator with o = E_nu = 0 (order 2: y = h x).
program test_lmto_real3
  use asr_kinds, only: dp
  use asr_lmto_ham, only: potpar_t
  use asr_augspace, only: graph_t
  use asr_augorb, only: orb_alloy_t
  use asr_lmto_alloy
  use asr_lmto_realization
  use asr_strconst, only: standard_alpha
  use asr_testing
  implicit none
  type(potpar_t) :: pa, pb
  type(graph_t) :: g
  type(orb_alloy_t) :: mdl
  type(lmto_real_op) :: op, oph
  integer, allocatable :: coord(:, :)
  logical, allocatable :: isa(:)
  real(dp), allocatable :: h(:, :), h3(:, :), x(:), y(:), od(:), pd(:), ed(:), t(:, :), hp(:, :), y2(:)
  real(dp) :: err, err2, scale
  integer :: l, n, i, j, s

  print '(a)', 'test_lmto_real3'
  do l = 0, 2                                     ! gamma-representation parameters (alpha_rep = gamma)
    pa%enu(l) = 0.1_dp * l; pa%c(l) = 0.3_dp - 0.15_dp * l; pa%sqrdel(l) = 0.30_dp - 0.1_dp * l
    pa%p(l) = 0.9_dp - 0.1_dp * l; pa%gamma(l) = 0.25_dp - 0.1_dp * l; pa%alpha_rep(l) = pa%gamma(l)
    pb%enu(l) = 0.05_dp * l; pb%c(l) = 0.45_dp - 0.1_dp * l; pb%sqrdel(l) = 0.28_dp - 0.1_dp * l
    pb%p(l) = 0.8_dp - 0.1_dp * l; pb%gamma(l) = 0.2_dp - 0.08_dp * l; pb%alpha_rep(l) = pb%gamma(l)
  end do
  call build_lmto_alloy('bcc', 3, 1.2_dp, 0.7_dp, pa, pb, 0.5_dp, g, coord, mdl, alpha_in=standard_alpha())
  allocate (isa(g%n)); isa = modulo(coord(1, :), 2) == 0
  call build_lmto_ordered(g, mdl, isa, op); op%order = 3
  oph = op; oph%order = 2; oph%o = 0.0_dp; oph%e = 0.0_dp      ! y = h x
  n = 9 * g%n
  allocate (h(n, n), h3(n, n), x(n), y(n), y2(n), od(n), pd(n), ed(n), t(n, n), hp(n, n))
  do j = 1, n
    x = 0.0_dp; x(j) = 1.0_dp
    call oph%apply(x, y); h(:, j) = y
  end do
  do i = 1, g%n
    do l = 1, 9
      od(9 * (i - 1) + l) = op%o(l, i); pd(9 * (i - 1) + l) = op%p(l, i); ed(9 * (i - 1) + l) = op%e(l, i)
    end do
  end do
  print '(a,i0,a,es10.3,a,es10.3)', '    dim ', n, '   max|p| ', maxval(abs(pd)), '   |h - h^T| ', maxval(abs(h - transpose(h)))
  do j = 1, n; t(:, j) = od * h(:, j); end do                 ! o h
  h3 = h - matmul(h, t)                                       ! h - h o h
  h3 = h3 + matmul(h, matmul(t, t))                           ! + h o h o h
  do j = 1, n; hp(:, j) = pd * h(:, j); end do                ! p h
  hp = matmul(h, hp)                                          ! h p h
  h3 = h3 - 0.5_dp * (matmul(hp, h) + matmul(h, hp))
  do i = 1, n; h3(i, i) = h3(i, i) + ed(i); end do
  err = 0.0_dp; err2 = 0.0_dp; scale = maxval(abs(h3))
  do s = 1, 3
    call random_number(x); x = x - 0.5_dp
    call op%apply(x, y)
    err = max(err, maxval(abs(y - matmul(h3, x))) / scale)
    op%order = 2; call op%apply(x, y2); op%order = 3
    err2 = max(err2, maxval(abs(y - y2)) / scale)
  end do
  print '(a,es10.3,a,es10.3)', '    max rel |H3 x - dense| ', err, '   (size of the third-order terms ', err2, ')'
  call check(err < 1.0e-12_dp, 'matrix-free third-order H equals the dense E + h - hoh + hohoh - (hphh + hhph)/2')
  call check(err2 > 1.0e-4_dp, 'the third-order terms are not trivially zero')
  call report()
end program test_lmto_real3
