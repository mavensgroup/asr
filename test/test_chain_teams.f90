!> The thread split of the chain step (asr_chain::chain_teams): one thread per chain while every thread has several chains;
!> teams that share the matrix-vector products once a rank holds fewer chains than threads (many MPI nodes).
program test_chain_teams
  use asr_chain, only: chain_teams
  use asr_testing
  implicit none
  integer :: nout, nin, nt, nth
  logical :: ok

  print '(a)', 'test_chain_teams'
  call chain_teams(256, 48, nout, nin)
  call check(nin == 1 .and. nout == 48, 'one node, 256 chains on 48 threads: one thread per chain (as before)')
  call chain_teams(5000, 28, nout, nin)
  call check(nin == 1 .and. nout == 28, 'desktop, 5000 chains on 28 threads: one thread per chain')
  call chain_teams(64, 48, nout, nin)
  call check(nin == 3 .and. nout == 16, '4 nodes (64 chains per rank, 48 threads): 16 teams of 3')
  call chain_teams(32, 48, nout, nin)
  call check(nin == 3 .and. nout == 16, '8 nodes (32 chains per rank, 48 threads): 16 teams of 3')
  call chain_teams(1, 28, nout, nin)
  call check(nin == 28 .and. nout == 1, 'one chain: all threads on its matrix-vector products')
  call chain_teams(10, 1, nout, nin)
  call check(nin == 1 .and. nout == 1, 'one thread: nothing to split')
  ok = .true.
  do nth = 1, 64
    do nt = 1, 300
      call chain_teams(nt, nth, nout, nin)
      if (nout * nin > nth .or. nout < 1 .or. nin < 1 .or. (nin == 1 .and. nout /= nth)) ok = .false.
    end do
  end do
  call check(ok, 'every split: nout x nin <= threads; nin = 1 keeps the full team')
  call report()
end program test_chain_teams
