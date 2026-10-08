!> The numerical realisation operator and the augmented-space operator describe the same Hamiltonian:
!>  (1) a symmetric matrix-free operator; (2) for a fully ordered realisation (x = 1) the recursion on the
!>  realisation equals the augmented-space recursion with the same parameters.
program test_lmto_real
  use asr_kinds, only: dp
  use asr_augspace, only: graph_t
  use asr_augorb
  use asr_lattices, only: periodic_graph_shells
  use asr_lmto_realization
  use asr_recursion, only: lanczos
  use asr_testing
  implicit none
  type(graph_t) :: g
  type(orb_alloy_t) :: mdl
  type(augorb_op) :: aop
  type(lmto_real_op) :: rop
  integer, allocatable :: coord(:, :), bvec(:, :)
  real(dp), allocatable :: xv(:), yv(:), zv(:), wv(:), v1(:), v2(:)
  real(dp) :: a1(5), b1(5), a2(5), b2(5)
  integer :: ns, k, l, m, nd
  print '(a)', 'test_lmto_real'
  call periodic_graph_shells('bcc', 10, 4, g, coord, bvec)
  ns = size(bvec, 2)
  mdl%nl = 3; mdl%x = 1.0_dp
  mdl%c0 = [-0.3_dp, 0.4_dp, 0.1_dp]; mdl%c1 = [0.5_dp, -0.25_dp, 0.2_dp]
  mdl%d0 = [0.6_dp, 0.45_dp, 0.3_dp]; mdl%d1 = [0.2_dp, 0.3_dp, -0.1_dp]
  mdl%o0 = [-0.15_dp, -0.1_dp, -0.2_dp]; mdl%o1 = [0.05_dp, 0.08_dp, 0.03_dp]
  mdl%e0 = [0.1_dp, -0.2_dp, 0.0_dp]; mdl%e1 = [0.07_dp, -0.04_dp, 0.02_dp]
  allocate (mdl%block(3, 3, ns), mdl%block0(3, 3))
  call random_seed()
  do k = 1, ns
    call random_number(mdl%block(:, :, k)); mdl%block(:, :, k) = mdl%block(:, :, k) - 0.5_dp
  end do
  ! Hermiticity: slot k and its inverse vector -v carry transposed blocks
  do k = 1, ns
    do l = 1, ns
      if (all(bvec(:, l) == -bvec(:, k))) mdl%block(:, :, l) = transpose(mdl%block(:, :, k))
    end do
  end do
  mdl%block0 = reshape([0.8_dp, 0.1_dp, 0.0_dp, 0.1_dp, -0.5_dp, 0.2_dp, 0.0_dp, 0.2_dp, 0.3_dp], [3, 3])
  call build_lmto_real(g, mdl, 7, rop)
  nd = rop%dim()
  allocate (xv(nd), yv(nd), zv(nd), wv(nd))
  call random_number(xv); call random_number(yv)
  call rop%apply(xv, zv); call rop%apply(yv, wv)
  call check_close(abs(dot_product(yv, zv) - dot_product(xv, wv)), 0.0_dp, 1.0e-9_dp, 'realisation operator is symmetric')
  ! ordered (all A): realisation recursion == augmented-space recursion (single node per site)
  call build_augorb_operator(g, mdl, 1, 4, aop)
  allocate (v1(aop%dim()), v2(nd)); v1 = 0.0_dp; v1(1) = 1.0_dp; v2 = 0.0_dp; v2(1) = 1.0_dp
  call lanczos(aop, v1, 5, a1, b1, m, reorth=.true.)
  call lanczos(rop, v2, 5, a2, b2, m, reorth=.true.)
  call check_close(max(maxval(abs(a1 - a2)), maxval(abs(b1(1:4) - b2(1:4)))), 0.0_dp, 1.0e-10_dp, 'ordered limit: realisation == augmented space (5 levels)')
  call report()
end program test_lmto_real
