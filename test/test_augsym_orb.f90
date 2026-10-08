!> D4h-reduced multi-orbital augmented space (asr_augsym) vs the unreduced one (asr_augorb): identical Lanczos
!> coefficients for the s, p_z, t2g(xy) and e_g(3z^2-r^2) start orbitals, with far fewer nodes; for the second-order
!> Hamiltonian (depth 2) and the third-order one (depth 1, three hops per level).
program test_augsym_orb
  use asr_kinds, only: dp
  use asr_lmto_ham, only: potpar_t
  use asr_augspace, only: graph_t
  use asr_augorb
  use asr_augsym
  use asr_lmto_alloy
  use asr_symmetry
  use asr_recursion, only: lanczos
  use asr_testing
  implicit none
  type(potpar_t) :: pa, pb
  type(graph_t) :: g
  type(orb_alloy_t) :: mdl
  type(augorb_op) :: op
  type(augsym_op) :: sop
  integer, allocatable :: coord(:, :), perm(:, :), mats(:, :, :)
  real(dp), allocatable :: d(:, :, :), v0(:), v1(:), a1(:), b1(:), a2(:), b2(:)
  real(dp) :: dmax
  integer :: l, io, m, m2, nd, ns, depth, ord
  integer, parameter :: l0(4) = [1, 3, 5, 7]

  print '(a)', 'test_augsym_orb'
  do l = 0, 2
    pa%enu(l) = 0.1_dp * l; pa%c(l) = 0.3_dp - 0.15_dp * l; pa%sqrdel(l) = 0.12_dp + 0.05_dp * l
    pa%p(l) = 0.9_dp - 0.1_dp * l; pa%gamma(l) = 0.2_dp; pa%alpha_rep(l) = 0.2_dp
    pb%enu(l) = 0.05_dp * l; pb%c(l) = 0.45_dp - 0.1_dp * l; pb%sqrdel(l) = 0.1_dp + 0.04_dp * l
    pb%p(l) = 0.8_dp - 0.1_dp * l; pb%gamma(l) = 0.2_dp; pb%alpha_rep(l) = 0.2_dp
  end do
  call build_lmto_alloy('fcc', 6, 1.2_dp, 0.7_dp, pa, pb, 0.4_dp, g, coord, mdl)
  call d4h_ops(coord, 12, perm, mats)
  call orbital_rep(mats, d)
  do ord = 2, 3
  depth = merge(2, 1, ord == 2)
  print '(a,i0,a,i0)', '  Hamiltonian order ', ord, '   depth ', depth
  call build_augorb_operator(g, mdl, 1, depth, op, maxrank=2, order=ord)
  call build_augsym_operator(g, mdl, perm, d, 1, depth, sop, maxrank=2, order=ord)
  print '(a,i0,a,i0,a,f6.1)', '    nodes unreduced ', op%nn, '   reduced ', sop%nn, '   ratio ', real(op%nn) / real(sop%nn)
  nd = op%dim(); ns = sop%dim()
  allocate (v0(nd), v1(ns), a1(depth + 1), b1(depth + 1), a2(depth + 1), b2(depth + 1))
  do io = 1, 4
    call augsym_set_start(sop, l0(io), d)
    v0 = 0.0_dp; v0(l0(io)) = 1.0_dp
    v1 = 0.0_dp; v1(l0(io)) = 1.0_dp
    call lanczos(op, v0, depth + 1, a1, b1, m, reorth=.true.)
    call lanczos(sop, v1, depth + 1, a2, b2, m2, reorth=.true.)
    dmax = max(maxval(abs(a1(1:depth + 1) - a2(1:depth + 1))), maxval(abs(b1(1:depth) - b2(1:depth))))
    print '(a,i0,a,es9.2)', '    start orbital ', l0(io), ': max |coefficient difference| = ', dmax
    call check_close(dmax, 0.0_dp, 1.0e-9_dp, 'reduced == unreduced coefficients')
  end do
  deallocate (v0, v1, a1, b1, a2, b2)
  end do
  call report()
end program test_augsym_orb
