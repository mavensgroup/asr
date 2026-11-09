!> asr_mpi: the serial stubs (one rank, rank 0, sums leave the data alone) and, in an MPI build run on any number of ranks
!> (mpirun -n N), the sums over ranks: exact for small integers in complex, the disjoint-share sum of chain_greens bitwise
!> exact on every rank, 64-bit and default integers.
program test_mpi
  use asr_kinds, only: dp
  use asr_mpi
  use asr_testing
  implicit none
  integer, parameter :: n = 1000
  complex(dp) :: x(n), y(n), yref(n)
  integer(8) :: k8(3)
  integer :: k4(2), i, nr, rk
  logical :: ok

  call par_init()
  nr = par_size(); rk = par_rank()
  print '(a,i0,a)', 'test_mpi (', nr, ' ranks)'
  call check(nr >= 1 .and. rk >= 0 .and. rk < nr, 'rank inside 0 .. size - 1')
  call check(par_root() .eqv. (rk == 0), 'root is rank 0')

  ! sum of (rank + 1, i) over ranks = (N (N + 1) / 2, N i), exact in double precision
  do i = 1, n
    x(i) = cmplx(real(rk + 1, dp), real(i, dp), dp)
  end do
  call par_sum_z(x, n)
  ok = .true.
  do i = 1, n
    ok = ok .and. x(i) == cmplx(real(nr * (nr + 1) / 2, dp), real(nr * i, dp), dp)
  end do
  call check(ok, 'complex sum over ranks exact')

  ! disjoint shares (as in chain_greens): each entry set on one rank, zero on the others; the sum equals the entry bitwise
  do i = 1, n
    yref(i) = cmplx(sin(0.37_dp * i) / 3.0_dp, exp(-1.0e-3_dp * i) / 7.0_dp, dp)
    y(i) = (0.0_dp, 0.0_dp)
    if (mod(i - 1, nr) == rk) y(i) = yref(i)
  end do
  call par_sum_z(y, n)
  call check(all(y == yref), 'disjoint-share sum bitwise equal to the entries')

  k8 = [int(rk, 8) + 1_8, 3_8, 2_8**40]
  call par_sum_i8(k8, 3)
  call check(all(k8 == [int(nr, 8) * (nr + 1) / 2, 3_8 * nr, 2_8**40 * nr]), '64-bit integer sum over ranks')
  k4 = [rk, 1]
  call par_sum_i(k4, 2)
  call check(all(k4 == [nr * (nr - 1) / 2, nr]), 'default integer sum over ranks')

  call par_finalize()
  call check(par_size() == 1 .and. par_rank() == 0, 'after par_finalize: serial values')
  call report()
end program test_mpi
