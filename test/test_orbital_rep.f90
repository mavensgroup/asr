!> Orbital representation of D4h on the 9 real harmonics: monomial (signed permutations), orthogonal, a group
!> homomorphism, and covariance of the screened structure-constant bond blocks (the convention to build on).
program test_orbital_rep
  use asr_kinds, only: dp
  use asr_symmetry
  use asr_lattices, only: periodic_graph_shells
  use asr_augspace, only: graph_t
  use asr_strconst, only: lattice_cluster, screened_cluster, avg_ws_radius, hardcore_alpha
  use asr_testing
  implicit none
  type(graph_t) :: g
  integer, allocatable :: coord(:, :), bvec(:, :), mats(:, :, :), perm(:, :)
  real(dp), allocatable :: d(:, :, :), vec(:, :), blk(:, :, :)
  real(dp) :: a(9, 9), b(9, 9), al4(0:3), alpha(0:2), eo, ec1, ec2, u(3), w
  integer :: iop, j, k, kk, ia, nmono, kp
  logical :: found
  character(len=3), parameter :: kinds(2) = ['fcc', 'bcc']
  integer :: ik

  print '(a)', 'test_orbital_rep'
  call periodic_graph_shells('fcc', 6, 4, g, coord, bvec)
  call d4h_ops(coord, 12, perm, mats)
  call orbital_rep(mats, d)
  nmono = 0; eo = 0.0_dp
  do iop = 1, 16
    if (all([(count(abs(d(:, j, iop)) > 1.0e-9_dp) == 1, j=1, 9)])) nmono = nmono + 1
    eo = max(eo, maxval(abs(matmul(d(:, :, iop), transpose(d(:, :, iop))) - eye())))
  end do
  call check(nmono == 16, 'all 16 D4h operations act monomially on the 9 orbitals')
  call check_close(eo, 0.0_dp, 1.0e-9_dp, 'orbital matrices are orthogonal')
  ! covariance of the screened bond blocks:  B(M v) = X B(v) X^T  for X = D or D^T ?
  al4 = hardcore_alpha(0.7_dp); alpha = al4(0:2)
  do ik = 1, 2
    call lattice_cluster(trim(kinds(ik)), 1.2_dp, vec)
    call screened_cluster(vec, avg_ws_radius(trim(kinds(ik))), alpha, blk)
    ec1 = 0.0_dp; ec2 = 0.0_dp
    do iop = 1, 16
      do j = 1, size(vec, 2)
        u = matmul(real(mats(:, :, iop), dp), vec(:, j))
        found = .false.
        do kk = 1, size(vec, 2)
          if (maxval(abs(vec(:, kk) - u)) < 1.0e-7_dp) then; kp = kk; found = .true.; exit; end if
        end do
        if (.not. found) cycle
        ! blk(:,:,j) = S_{R_j,0}(L,L');  test  blk(kp) = D blk(j) D^T   and   D^T blk D
        a = matmul(d(:, :, iop), matmul(blk(:, :, j), transpose(d(:, :, iop))))
        ec1 = max(ec1, maxval(abs(a - blk(:, :, kp))))
        b = matmul(transpose(d(:, :, iop)), matmul(blk(:, :, j), d(:, :, iop)))
        ec2 = max(ec2, maxval(abs(b - blk(:, :, kp))))
      end do
    end do
    print '(a,a,a,es9.2,a,es9.2)', '    ', kinds(ik), ': |B(Mv) - D B(v) D^T| max = ', ec1, '    |B(Mv) - D^T B(v) D| max = ', ec2
    call check(min(ec1, ec2) < 5.0e-3_dp, kinds(ik)//': screened S blocks are covariant under one of the two conventions')
  end do
  call report()
contains
  function eye() result(e)
    real(dp) :: e(9, 9)
    integer :: i
    e = 0.0_dp
    do i = 1, 9; e(i, i) = 1.0_dp; end do
  end function eye
end program test_orbital_rep
