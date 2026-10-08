!> Point-group reduction of the augmented space: recursion coefficients identical to the unreduced ones.
program test_augsym
  use asr_kinds, only: dp
  use asr_operator, only: csr_op
  use asr_augspace
  use asr_lattices
  use asr_symmetry
  use asr_recursion, only: lanczos
  use asr_testing
  implicit none
  character(len=3), parameter :: kinds(3) = ['sc ', 'bcc', 'fcc']
  integer, parameter :: depth = 5
  type(graph_t) :: g
  type(binary_model_t) :: mdl
  type(csr_op) :: full, red
  integer, allocatable :: coord(:, :), perm(:, :)
  real(dp) :: a1(depth), b1(depth), a2(depth), b2(depth), v1(1), dmax
  real(dp), allocatable :: v0(:)
  integer :: ik, m, l
  print '(a)', 'test_augsym'
  mdl%x = 0.4_dp; mdl%epsa = 0.7_dp; mdl%epsb = -0.7_dp; mdl%tAA = -1.1_dp; mdl%tAB = -1.0_dp; mdl%tBB = -0.8_dp
  l = 7
  do ik = 1, 3
    call periodic_graph(trim(kinds(ik)), l, g, coord)
    call cubic_site_perms(coord, 2 * l, perm)
    ! the origin is site 1 (coordinate (0,0,0))
    call build_augmented_operator(g, mdl, 1, depth, full)
    call build_augmented_operator_sym(g, mdl, perm, 1, depth, red)
    allocate (v0(full%n)); v0 = 0.0_dp; v0(1) = 1.0_dp
    call lanczos(full, v0, depth, a1, b1, m, reorth=.true.)
    deallocate (v0); allocate (v0(red%n)); v0 = 0.0_dp; v0(1) = 1.0_dp
    call lanczos(red, v0, depth, a2, b2, m, reorth=.true.)
    deallocate (v0)
    dmax = max(maxval(abs(a1 - a2)), maxval(abs(b1 - b2)))
    print '(a,a,a,i0,a,i0,a,f6.1,a,es9.2)', '    ', kinds(ik), '  states ', full%n, ' -> ', red%n, '   reduction x', real(full%n, dp) / red%n, '   max |coef diff| = ', dmax
    call check(red%n < full%n / 4, trim(kinds(ik))//': at least 4x fewer states')
    call check_close(dmax, 0.0_dp, 1.0e-10_dp, trim(kinds(ik))//': recursion coefficients identical')
  end do
  call report()
end program test_augsym
