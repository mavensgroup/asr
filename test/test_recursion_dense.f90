!> M1 gate 1: recursion == dense diagonalisation when the Krylov space is
!> exhausted (continued fraction without any tail).
program test_recursion_dense
  use asr_kinds, only: dp
  use asr_operator, only: csr_op, csr_from_triplets
  use asr_recursion, only: lanczos
  use asr_greens, only: cf_green
  use asr_testing
  implicit none
  integer, parameter :: n = 60
  real(dp) :: h(n, n), hh(n, n), ev(n), v0(n), work(8 * n), a(n), b(n), r
  integer :: i, j, nt, info, m, k
  integer, allocatable :: ti(:), tj(:)
  real(dp), allocatable :: tv(:)
  type(csr_op) :: op
  complex(dp) :: z, gd, gr
  real(dp) :: maxerr, mom_dense, mom_rec

  print '(a)', 'test_recursion_dense'
  call random_seed()
  call random_number(h)
  h = 0.5_dp * (h + transpose(h))
  do i = 1, n; do j = 1, n
      if (abs(i - j) > 7 .and. mod(i + j, 3) /= 0) h(i, j) = 0.0_dp
  end do; end do
  h = 0.5_dp * (h + transpose(h))
  nt = count(h /= 0.0_dp)
  allocate (ti(nt), tj(nt), tv(nt)); k = 0
  do i = 1, n; do j = 1, n
      if (h(i, j) /= 0.0_dp) then; k = k + 1; ti(k) = i; tj(k) = j; tv(k) = h(i, j); end if
  end do; end do
  call csr_from_triplets(n, ti, tj, tv, op)

  v0 = 0.0_dp; v0(7) = 1.0_dp
  call lanczos(op, v0, n, a, b, m, reorth=.true.)
  call check(m <= n, 'krylov dimension <= n')

  hh = h
  call dsyev('V', 'U', n, hh, n, ev, work, size(work), info)
  call check(info == 0, 'dsyev converged')

  maxerr = 0.0_dp
  do k = 1, 25
    z = cmplx(-3.0_dp + 0.25_dp * k, 0.1_dp, dp)
    ! dense  <7|(z-H)^-1|7>
    gd = (0.0_dp, 0.0_dp)
    do j = 1, n
      gd = gd + hh(7, j)**2 / (z - ev(j))
    end do
    gr = cf_green(a, b, m, z, (0.0_dp, 0.0_dp))
    maxerr = max(maxerr, abs(gd - gr))
  end do
  call check_close(maxerr, 0.0_dp, 1.0e-9_dp, 'G(z) recursion vs dense')

  ! second moment <7|H^2|7> = a1^2 + b1^2
  mom_dense = dot_product(h(:, 7), h(:, 7))
  mom_rec = a(1)**2 + b(1)**2
  r = abs(mom_dense - mom_rec)
  call check_close(r, 0.0_dp, 1.0e-12_dp, 'second moment')
  call report()
end program test_recursion_dense
