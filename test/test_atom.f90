!> Free-atom LDA vs NIST reference values (non-relativistic, VWN):
!>   He: Etot = -2.834836 Ha, eps_1s = -0.570425 Ha;  Ne: Etot = -128.233481 Ha.
program test_atom
  use asr_kinds, only: dp
  use asr_atom
  use asr_testing
  implicit none
  type(atom_result_t) :: r
  type(shell_t) :: he(1), ne(3)
  real(dp), parameter :: ha = 2.0_dp          ! Ry per Ha

  print '(a)', 'test_atom'
  he = [shell_t(1, 0, 2.0_dp)]
  call solve_atom(2.0_dp, he, 'vwn', 40.0_dp, 0.01_dp, .false., r)
  print '(a,l1,a,i0,a,f14.7,a,f12.6)', '    He converged=', r%converged, ' iters=', r%iter, '  Etot(Ha)=', r%etot / ha, '  eps1s(Ha)=', r%eps(1) / ha
  call check(r%converged, 'He SCF converged')
  call check_close(r%etot / ha, -2.834836_dp, 2.0e-5_dp, 'He total energy vs NIST LDA')
  call check_close(r%eps(1) / ha, -0.570425_dp, 2.0e-5_dp, 'He eps_1s vs NIST LDA')

  ne = [shell_t(1, 0, 2.0_dp), shell_t(2, 0, 2.0_dp), shell_t(2, 1, 6.0_dp)]
  call solve_atom(10.0_dp, ne, 'vwn', 30.0_dp, 0.01_dp, .false., r)
  print '(a,l1,a,i0,a,f14.6,a,3f12.6)', '    Ne converged=', r%converged, ' iters=', r%iter, '  Etot(Ha)=', r%etot / ha, '  eps(Ha)=', r%eps / ha
  call check(r%converged, 'Ne SCF converged')
  call check_close(r%etot / ha, -128.233481_dp, 5.0e-4_dp, 'Ne total energy vs NIST LDA')
  call report()
end program test_atom
