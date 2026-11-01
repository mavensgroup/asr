!> asr_screening: lattice Coulomb kernel and the single-site screened electrostatics with first-shell SRO (gate B3).
!>  Ewald M(k): CsCl (bcc H) and rock-salt (sc R) Madelung constants.  screen_lattice: beta round trip, the unscreened limit,
!>  the reference values of the derivation (bcc a = 2.85 A: a = 1.263 Ry gives beta 1.045, lambda 0.212), the lattice maximum
!>  of beta, mesh convergence.  screen_matrix: no SRO = SIM, diag(w) F symmetric, sublattice sum rule, the binary 50:50 closed
!>  form a [1/(G00 + z G1 alpha) - 1], and the unscreened limit = bare first-shell charge.
program test_screening
  use asr_kinds, only: dp, pi
  use asr_screening
  use asr_sro, only: sro_conditional
  use asr_testing
  implicit none
  type(screen_t) :: s, s2
  character(len=:), allocatable :: msg
  real(dp) :: p(3, 3), kv(3, 1), m(1), alat, al, x(3), at(3, 3), ah(3, 3), f(3, 3), w(3), f2(2, 2), ah2(2, 2), ev
  real(dp), allocatable :: cond(:, :)
  integer :: i, j

  print '(a)', 'test_screening'
  alat = 2.85_dp / 0.529177210903_dp
  p = reshape([-1, 1, 1, 1, -1, 1, 1, 1, -1], [3, 3]) * 0.5_dp * alat
  kv(:, 1) = [2.0_dp * pi / alat, 0.0_dp, 0.0_dp]
  call coulomb_k(p, kv, m)
  call check_close(-m(1) * sqrt(3.0_dp) / 2.0_dp * alat / 2.0_dp, 1.7626748_dp, 1.0e-6_dp, 'bcc H point: CsCl Madelung constant')
  p = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3]) * alat
  kv(:, 1) = [1.0_dp, 1.0_dp, 1.0_dp] * pi / alat
  call coulomb_k(p, kv, m)
  call check_close(-m(1) * alat / 2.0_dp, 1.7475646_dp, 1.0e-6_dp, 'sc R point: rock-salt Madelung constant')

  call screen_lattice('bcc', alat, 1.0_dp, s, msg)
  call check(len(msg) == 0 .and. .not. s%unscreened, 'bcc beta = 1 accepted')
  call check_close(-s%w00 * s%r1 / (2.0_dp * s%g00), 1.0_dp, 1.0e-10_dp, 'beta round trip: -W00 R1 / (2 G00) = beta')
  call check_close(s%w00, s%a * (1.0_dp - s%g00), 1.0e-12_dp, 'W00 = a (1 - G00)')
  call check(s%z == 8 .and. abs(s%r1 - sqrt(3.0_dp) / 2.0_dp * alat) < 1.0e-10_dp, 'bcc first shell: z = 8 at sqrt(3)/2 a')
  call check_close(s%a_min * s%r1 / 2.0_dp, 1.7626748_dp, 1.0e-6_dp, 'bcc stability bound: a_min R1 / 2 = CsCl Madelung constant (H)')
  call check(s%beta_max > 1.135_dp .and. s%beta_max < 1.7626748_dp, &
             'bcc beta maximum below a_min R1/2 (G00 stays finite at the instability)')
  call screen_lattice('bcc', alat, 1.045_dp, s2, msg)
  call check_close(s2%a, 1.263_dp, 0.01_dp, 'derivation reference: beta 1.045 <-> a 1.263 Ry')
  call check_close(s2%lambda, 0.212_dp, 0.003_dp, 'derivation reference: lambda 0.212')
  call screen_lattice('bcc', alat, 1.045_dp, s, msg, nk=36)
  call check_close(s%lambda, s2%lambda, 2.0e-3_dp, 'mesh 28 vs 36: lambda converged')
  call screen_lattice('bcc', alat, 1.8_dp, s, msg)
  call check(index(msg, 'maximum') > 0, 'beta above the lattice maximum rejected')
  call screen_lattice('fcc', 3.59_dp / 0.529177210903_dp, 1.0_dp, s, msg)
  call check(len(msg) == 0 .and. s%z == 12, 'fcc beta = 1 accepted, z = 12')
  call check(s%lambda > 0.0_dp .and. s%lambda < 1.0_dp, 'fcc: 0 < lambda < 1')

  ! single-site matrix, one sublattice, ternary with SRO
  call screen_lattice('bcc', alat, 1.0_dp, s, msg)
  x = [0.66_dp, 0.14_dp, 0.20_dp]
  at = 0.0_dp; at(1, 2) = -0.1_dp; at(2, 1) = -0.1_dp; at(1, 3) = 0.2_dp; at(3, 1) = 0.2_dp; at(2, 3) = 0.3_dp; at(3, 2) = 0.3_dp
  call sro_conditional(x, x, at, .true., cond, msg)
  do i = 1, 3; do j = 1, 3; ah(i, j) = s%z * (cond(i, j) - x(j)); end do; end do
  call screen_matrix(s, 0.0_dp * ah, f)
  call check_close(maxval(abs(f - (s%w00 / s%g00) * reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3]))), 0.0_dp, 1.0e-14_dp, &
                   'no SRO: F = (W00/G00) I = -2 beta / R1 (the SIM)')
  call check_close(s%w00 / s%g00, -2.0_dp * s%beta / s%r1, 1.0e-10_dp, 'no SRO: W00/G00 = -2 beta / R1')
  call screen_matrix(s, ah, f)
  w = x
  call check_close(maxval(abs(spread(w, 2, 3) * f - transpose(spread(w, 2, 3) * f))), 0.0_dp, 1.0e-13_dp, &
                   'ternary SRO: diag(w) F symmetric (variational energy)')
  call check_close(maxval(abs(matmul(x, f) - (s%w00 / s%g00) * x)), 0.0_dp, 1.0e-13_dp, &
                   'ternary SRO: x^T F = (W00/G00) x^T (no potential from the sublattice mean)')
  ! binary 50:50: eigenvalue on (1, -1) is a [1/(G00 + z G1 alpha) - 1]
  al = -0.15_dp
  ah2 = reshape([0.5_dp * al, -0.5_dp * al, -0.5_dp * al, 0.5_dp * al], [2, 2]) * s%z
  call screen_matrix(s, ah2, f2)
  ev = f2(1, 1) - f2(1, 2)
  call check_close(ev, s%a * (1.0_dp / (s%g00 + s%z * s%g1 * al) - 1.0_dp), 1.0e-12_dp, &
                   'binary 50:50: F(1,-1) = a [1/(G00 + z G1 alpha) - 1]')
  ! unscreened limit: F = bare first-shell charge
  call screen_lattice('bcc', alat, 0.0_dp, s, msg)
  call check(len(msg) == 0 .and. s%unscreened, 'beta = 0: unscreened limit')
  call screen_matrix(s, ah, f)
  call check_close(maxval(abs(f - (2.0_dp / s%r1) * ah)), 0.0_dp, 1.0e-14_dp, 'beta = 0: F = (2/R1) Ahat (bare shell charge)')
  call report()
end program test_screening
