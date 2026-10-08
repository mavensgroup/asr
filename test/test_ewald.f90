!> Ewald Madelung sums vs the classic ionic-lattice Madelung constants:
!>   CsCl 1.762675 (per nearest-neighbour distance), NaCl 1.747565; independence of cell choice.
program test_ewald
  use asr_kinds, only: dp
  use asr_ewald
  use asr_testing
  implicit none
  real(dp) :: phi(2, 2), tau(3, 2), a1(3, 3), v, r0, alat
  print '(a)', 'test_ewald'
  alat = 5.0_dp
  ! CsCl: simple cubic, +1 at 0, -1 at body centre
  a1 = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
  tau(:, 1) = 0.0_dp; tau(:, 2) = 0.5_dp
  call madelung_phi(a1, tau, alat, phi)
  v = phi(1, 1) - phi(1, 2)                  ! potential at the + site from q = (+1, -1), per charge
  r0 = alat * sqrt(3.0_dp) / 2.0_dp
  call check_close(-v * r0, 1.762675_dp, 2.0e-6_dp, 'CsCl Madelung constant')
  call check_close(phi(1, 1) - phi(1, 2), phi(2, 2) - phi(2, 1), 1.0e-10_dp, 'sites equivalent (+/-) up to sign')
  ! NaCl: fcc, +1 at 0, -1 at (1/2,0,0)
  a1 = reshape([0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp], [3, 3])
  tau(:, 1) = 0.0_dp; tau(:, 2) = [0.5_dp, 0.0_dp, 0.0_dp]
  call madelung_phi(a1, tau, alat, phi)
  v = phi(1, 1) - phi(1, 2)
  call check_close(-v * alat / 2.0_dp, 1.747565_dp, 2.0e-6_dp, 'NaCl Madelung constant')
  ! same CsCl in a doubled (non-primitive) setting must agree
  a1 = reshape([2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
  block
    real(dp) :: t4(3, 4), p4(4, 4), vv
    t4(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]; t4(:, 2) = [0.5_dp, 0.5_dp, 0.5_dp]
    t4(:, 3) = [1.0_dp, 0.0_dp, 0.0_dp]; t4(:, 4) = [1.5_dp, 0.5_dp, 0.5_dp]
    call madelung_phi(a1, t4, alat, p4)
    vv = p4(1, 1) - p4(1, 2) + p4(1, 3) - p4(1, 4)
    call check_close(-vv * r0, 1.762675_dp, 2.0e-6_dp, 'CsCl in a doubled cell gives the same constant')
  end block
  call report()
end program test_ewald
