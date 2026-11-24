!> Ewald Madelung sums vs the classic ionic-lattice Madelung constants:
!>   CsCl 1.762675 (per nearest-neighbour distance), NaCl 1.747565; independence of cell choice.
program test_ewald
  use asr_kinds, only: dp
  use asr_crystal, only: structure_t
  use asr_prototypes, only: prototype_structure
  use asr_sublattice, only: structure_madelung
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
  block
    use asr_crystal, only: structure_init
    type(structure_t) :: st, st2
    real(dp), allocatable :: ps(:, :), ps2(:, :)
    real(dp) :: a2(3, 3), f8(3, 8), p4(4, 4), q4(4), vc, va
    integer :: b
    character(len=:), allocatable :: msg
    call prototype_structure('b2', st, msg)
    call structure_madelung(st, alat, ps)
    call check_close(-0.5_dp * (ps(1, 1) - ps(1, 2)) * r0, 1.762675_dp, 2.0e-6_dp, 'structure_madelung on b2 gives CsCl constant')
    call prototype_structure('nias', st, msg, c_over_a=1.6_dp)
    call structure_madelung(st, alat, ps)
    call madelung_phi(st%avec, st%tau, alat, p4)
    q4 = [1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp]
    vc = 2.0_dp * dot_product(p4(1, :), q4); va = 2.0_dp * dot_product(p4(3, :), q4)
    call check_close(max(abs((ps(1, 1) - ps(1, 2)) - vc), abs((ps(2, 1) - ps(2, 2)) - va)), 0.0_dp, 1.0e-12_dp, &
                     'structure_madelung on NiAs: neutral shifts equal site Ewald potentials')
    a2 = st%avec; a2(:, 1) = 2.0_dp * a2(:, 1)
    do b = 1, 4
      f8(:, b) = [0.5_dp * st%frac(1, b), st%frac(2, b), st%frac(3, b)]
      f8(:, b + 4) = [0.5_dp * st%frac(1, b) + 0.5_dp, st%frac(2, b), st%frac(3, b)]
    end do
    call structure_init(st2, a2, f8, [1, 1, 2, 2, 1, 1, 2, 2], msg)
    call structure_madelung(st2, alat, ps2)
    call check_close(maxval(abs((ps2(:, 1) - ps2(:, 2)) - (ps(:, 1) - ps(:, 2)))), 0.0_dp, 1.0e-11_dp, &
                     'structure_madelung on NiAs: doubled cell gives the same neutral shifts')
    call check_close(ps(1, 2), ps(2, 1), 1.0e-12_dp, 'structure_madelung on NiAs: symmetric matrix for equal sublattice sizes')
    call check(ps(1, 1) - ps(1, 2) < -0.1_dp, 'structure_madelung on NiAs: opposite sublattice charges lower energy')
  end block
  call report()
end program test_ewald
