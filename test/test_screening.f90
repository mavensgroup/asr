!> asr_screening: lattice Coulomb kernel and the single-site screened electrostatics with first-shell SRO.
!>  Ewald M(k): CsCl (bcc H) and rock-salt (sc R) Madelung constants.  screen_lattice: beta round trip, the unscreened limit,
!>  the reference values of the derivation (bcc a = 2.85 A: a = 1.263 Ry gives beta 1.045, lambda 0.212), the lattice maximum
!>  of beta, mesh convergence.  screen_matrix: no SRO = SIM (linear form, bit for bit), diag(w) F symmetric, sublattice sum
!>  rule, the unscreened limit = bare first-shell charge.  Closure (screen_gamma): reference values of the numpy minimal-tail
!>  closure that was tested against exact supercells (tools/sro_closure/sro_closure.py `clip`, same nk = 28 mesh),
!>  Gamma = 0 at full clustering, Gamma = G(H) for perfect B2 order (exact), positive and decreasing on the realizable range,
!>  exactly linear on the realizable first-shell range, continuous at its edges and at the condensates, finite where the
!>  first-shell form has its pole (fcc beta 1.2).
!>  screen_update: the drive iteration converges to V = F dq, and to dq = 0 on a mode with Gamma = 0.
program test_screening
  use asr_kinds, only: dp, pi
  use asr_screening
  use asr_sro, only: sro_conditional
  use asr_testing
  implicit none
  type(screen_t) :: s, s2
  character(len=:), allocatable :: msg
  real(dp) :: p(3, 3), kv(3, 1), m(1), alat, afcc, al, x(3), at(3, 3), ah(3, 3), f(3, 3), w(3), f2(2, 2), ah2(2, 2), ev
  real(dp) :: g3(3, 3), w3(3, 3), gi3(3, 3), g2(2, 2), w2(2, 2), gi2(2, 2), lam2(2), gam2(2), lam3(3), gam3(3)
  real(dp) :: h, gp, gm_, gprev, gnow, lam, dmax, dq3(3), u3(3), uo3(3), v3(3)
  real(dp), allocatable :: cond(:, :)
  integer :: i, j, k
  logical :: ok

  print '(a)', 'test_screening'
  alat = 2.85_dp / 0.529177210903_dp
  afcc = 3.59_dp / 0.529177210903_dp
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
  call screen_lattice('fcc', afcc, 1.0_dp, s, msg)
  call check(len(msg) == 0 .and. s%z == 12, 'fcc beta = 1 accepted, z = 12')
  call check(s%lambda > 0.0_dp .and. s%lambda < 1.0_dp, 'fcc: 0 < lambda < 1')
  call check_close(s%c1_ord, -1.0_dp / 3.0_dp, 1.0e-14_dp, 'fcc ordering point X: c1 = -1/3')

  ! ---- single-site matrix, one sublattice, ternary with SRO
  call screen_lattice('bcc', alat, 1.0_dp, s, msg)
  call check_close(s%c1_ord, -1.0_dp, 1.0e-14_dp, 'bcc ordering point H: c1 = -1')
  x = [0.66_dp, 0.14_dp, 0.20_dp]; w = x
  at = 0.0_dp; at(1, 2) = -0.1_dp; at(2, 1) = -0.1_dp; at(1, 3) = 0.2_dp; at(3, 1) = 0.2_dp; at(2, 3) = 0.3_dp; at(3, 2) = 0.3_dp
  call sro_conditional(x, x, at, .true., cond, msg)
  do i = 1, 3; do j = 1, 3; ah(i, j) = s%z * (cond(i, j) - x(j)); end do; end do
  call screen_matrix(s, 0.0_dp * ah, w, fm=f, gm=g3, wm=w3)
  call check_close(maxval(abs(f - (s%w00 / s%g00) * eye3())), 0.0_dp, 0.0_dp, 'no SRO: F = (W00/G00) I = -2 beta / R1 (the SIM), bit for bit')
  call check_close(maxval(abs(g3 - s%g00 * eye3())) + maxval(abs(w3 - s%w00 * eye3())), 0.0_dp, 0.0_dp, &
                   'no SRO: Gamma = G00 I, W = W00 I exactly')
  call check_close(s%w00 / s%g00, -2.0_dp * s%beta / s%r1, 1.0e-10_dp, 'no SRO: W00/G00 = -2 beta / R1')
  call screen_matrix(s, ah, w, fm=f, gm=g3, wm=w3, lam=lam3, gam=gam3)
  call check_close(maxval(abs(spread(w, 2, 3) * f - transpose(spread(w, 2, 3) * f))), 0.0_dp, 1.0e-13_dp, &
                   'ternary SRO: diag(w) F symmetric (variational energy)')
  call check_close(maxval(abs(matmul(x, f) - (s%w00 / s%g00) * x)), 0.0_dp, 1.0e-13_dp, &
                   'ternary SRO: x^T F = (W00/G00) x^T (no potential from the sublattice mean)')
  call check_close(maxval(abs(matmul(g3, [1.0_dp, 1.0_dp, 1.0_dp]) - s%g00)), 0.0_dp, 1.0e-13_dp, &
                   'ternary SRO: Gamma 1 = G00 1 (the sublattice mean is not correlated)')
  call check_close(maxval(abs(matmul(g3, ah) - matmul(ah, g3))), 0.0_dp, 1.0e-12_dp, 'ternary SRO: Gamma is a function of Ahat (commutes)')
  call check_close(maxval(abs(w3 - s%a * (eye3() - g3))), 0.0_dp, 1.0e-12_dp, 'ternary SRO: W = a (I - Gamma)')
  call check_close(maxval(abs(f - matmul(inv3(g3), w3))), 0.0_dp, 1.0e-12_dp, 'ternary SRO: F = Gamma^-1 W')
  dmax = 0.0_dp
  do i = 1, 3; dmax = max(dmax, abs(gam3(i) - screen_gamma(s, lam3(i)))); end do
  call check_close(dmax, 0.0_dp, 1.0e-14_dp, 'ternary SRO: mode Gammas are the scalar closure at the eigenvalues of Ahat')
  ! binary 50:50: eigenvalue on (1, -1) is a [1/Gamma(z alpha) - 1]
  al = -0.15_dp
  ah2 = reshape([0.5_dp * al, -0.5_dp * al, -0.5_dp * al, 0.5_dp * al], [2, 2]) * s%z
  call screen_matrix(s, ah2, [0.5_dp, 0.5_dp], fm=f2)
  ev = f2(1, 1) - f2(1, 2)
  call check_close(ev, s%a * (1.0_dp / screen_gamma(s, s%z * al) - 1.0_dp), 1.0e-12_dp, &
                   'binary 50:50: F(1,-1) = a [1/Gamma(z alpha) - 1]')

  ! ---- closure: reference values (numpy closure, nk = 28; tested against exact supercells)
  call screen_lattice('fcc', afcc, 1.0_dp, s, msg)
  call check_close(screen_gamma(s, 7.92_dp), 0.455347061157_dp, 2.0e-8_dp, 'closure fcc beta 1.0, lambda 7.92 (clipped) = numpy')
  call check_close(screen_gamma(s, 2.0_dp), s%g00 + 2.0_dp * s%g1, 0.0_dp, 'closure fcc beta 1.0, lambda 2: realizable, linear form exactly')
  call check_close(screen_gamma(s, -2.5_dp), 1.620956472772_dp, 2.0e-8_dp, 'closure fcc beta 1.0, lambda -2.5 (ordering, clipped) = numpy')
  call check_close(screen_gamma(s, 5.0_dp), 0.797934787698_dp, 2.0e-8_dp, 'closure fcc beta 1.0, lambda 5 = numpy')
  call check_close(screen_gamma(s, 11.9_dp), 0.010717143990_dp, 2.0e-8_dp, 'closure fcc beta 1.0, lambda 11.9 (near the condensate) = numpy')
  call check_close(s%g00 + s%g1 * 7.92_dp, 0.486728364529_dp, 1.0e-9_dp, 'same lattice functions as numpy (G00 + G1 7.92)')
  ! slope G1 at lambda = 0 (the linear form is the first order of the closure)
  h = 0.01_dp
  gp = screen_gamma(s, h); gm_ = screen_gamma(s, -h)
  call check_close((gp - gm_) / (2.0_dp * h), s%g1, 1.0e-3_dp * abs(s%g1), 'closure: dGamma/dlambda = G1 at lambda = 0')
  call check_close(screen_gamma(s, 0.0_dp), s%g00, 0.0_dp, 'closure: Gamma(0) = G00')
  ! full clustering: Gamma = 0
  call check_close(screen_gamma(s, real(s%z, dp)), 0.0_dp, 1.0e-15_dp, 'closure: Gamma = 0 at full clustering (lambda = z)')
  ! positive and decreasing on the realizable range z c1(X) .. z, continuous at both limits
  call check_positive_decreasing(s, 'fcc beta 1.0')
  call check_close(screen_gamma(s, 3.0_dp + 1.0e-7_dp), s%g00 + 3.0_dp * s%g1, 1.0e-6_dp, 'closure: continuous at lambda = 3 (end of the linear range, fcc)')
  call check_close(screen_gamma(s, -1.0_dp - 1.0e-7_dp), s%g00 - s%g1, 1.0e-6_dp, 'closure: continuous at lambda = -1 (end of the linear range)')
  call check_close(screen_gamma(s, s%z * s%alpha_cl * (1.0_dp + 1.0e-9_dp)), screen_gamma(s, s%z * s%alpha_cl * (1.0_dp - 1.0e-9_dp)), &
                   1.0e-6_dp, 'closure: continuous at the clustering condensate')
  call check_close(screen_gamma(s, s%z * s%alpha_ord * (1.0_dp + 1.0e-9_dp)), screen_gamma(s, s%z * s%alpha_ord * (1.0_dp - 1.0e-9_dp)), &
                   1.0e-6_dp, 'closure: continuous at the ordering condensate')
  ! fcc beta 1.2: the first-shell form has its pole at lambda = 7.67; the closure is finite and positive
  call screen_lattice('fcc', afcc, 1.2_dp, s, msg)
  call check(s%g00 + s%g1 * 7.92_dp < 0.0_dp, 'fcc beta 1.2: first-shell form negative at lambda 7.92 (past its pole)')
  call check_close(screen_gamma(s, 7.92_dp), 0.353869681688_dp, 2.0e-8_dp, 'closure fcc beta 1.2, lambda 7.92: positive, = numpy')
  call check_positive_decreasing(s, 'fcc beta 1.2')
  ! bcc: ordering side (CuZn), and perfect B2 order = G(H) exactly
  call screen_lattice('bcc', alat, 1.0_dp, s, msg)
  call check_close(screen_gamma(s, -1.184_dp), 1.583643444846_dp, 2.0e-8_dp, 'closure bcc beta 1.0, lambda -1.184 (CuZn alpha -0.148) = numpy')
  call check(abs(screen_gamma(s, -1.184_dp) / (s%g00 - 1.184_dp * s%g1) - 1.0_dp) < 2.0e-4_dp, 'closure: CuZn alpha -0.148 within 0.02 % of the linear form')
  call check_close(screen_gamma(s, -2.13_dp), 1.694157745922_dp, 2.0e-8_dp, 'closure bcc beta 1.0, lambda -2.13 (CuZn alpha -0.266) = numpy')
  call check_close(screen_gamma(s, 2.0_dp), 1.122360183532_dp, 2.0e-8_dp, 'closure bcc beta 1.0, lambda 2 (clustering, clipped) = numpy')
  call check_close(screen_gamma(s, -real(s%z, dp)), s%a / (s%a - 2.0_dp * 1.7626748_dp / s%r1), 1.0e-6_dp, &
                   'closure: perfect B2 order (alpha_1 = -1) gives Gamma = G(H) = a / (a + M(H)) exactly')
  call check_positive_decreasing(s, 'bcc beta 1.0')

  ! ---- full clustering in screen_matrix: binary 50:50, alpha = 1
  call screen_lattice('fcc', afcc, 1.0_dp, s, msg)
  al = 1.0_dp
  ah2 = reshape([0.5_dp * al, -0.5_dp * al, -0.5_dp * al, 0.5_dp * al], [2, 2]) * s%z
  call screen_matrix(s, ah2, [0.5_dp, 0.5_dp], gm=g2, wm=w2, ginv=gi2, lam=lam2, gam=gam2)
  call check_close(g2(1, 1) - g2(1, 2), 0.0_dp, 1.0e-12_dp, 'alpha = 1: Gamma = 0 on the mode (1, -1)')
  call check_close(w2(1, 1) - w2(1, 2), s%a, 1.0e-12_dp, 'alpha = 1: W = a on the mode (1, -1) (V = a db, finite)')
  call check_close(gi2(1, 1) - gi2(1, 2), 0.0_dp, 1.0e-12_dp, 'alpha = 1: pseudo-inverse leaves the Gamma = 0 mode out')
  call check_close(g2(1, 1) + g2(1, 2), s%g00, 1.0e-12_dp, 'alpha = 1: Gamma = G00 on the uniform mode')

  ! ---- unscreened limit: F = bare first-shell charge, unchanged
  call screen_lattice('bcc', alat, 0.0_dp, s, msg)
  call check(len(msg) == 0 .and. s%unscreened, 'beta = 0: unscreened limit')
  call screen_matrix(s, ah, w, fm=f)
  call check_close(maxval(abs(f - (2.0_dp / s%r1) * ah)), 0.0_dp, 0.0_dp, 'beta = 0: F = (2/R1) Ahat (bare shell charge), as before')

  ! ---- screen_update: drive iteration with a model response dq = P (c - chi v)
  ! (a) moderate SRO (bcc ternary above): the fixed point is v = F dq
  call screen_lattice('bcc', alat, 1.0_dp, s, msg)
  call screen_matrix(s, ah, w, fm=f, gm=g3, wm=w3, ginv=gi3)
  call drive_fixed_point(g3, w3, x, [0.10_dp, -0.30_dp, -0.06_dp], 2.0_dp, ok, dmax, f)
  call check(ok, 'drive iteration converges (moderate SRO)')
  call check_close(dmax, 0.0_dp, 1.0e-10_dp, 'drive fixed point: v = F dq (the old form) for moderate SRO')
  ! (b) fcc (1,1,0) triad: Fe has only Fe neighbours, Ahat eigenvalue 12, Gamma = 0 on that mode
  call screen_lattice('fcc', afcc, 1.0_dp, s, msg)
  at = 0.0_dp; at(1, 2) = 1.0_dp; at(2, 1) = 1.0_dp; at(1, 3) = 1.0_dp; at(3, 1) = 1.0_dp
  call sro_conditional(x, x, at, .true., cond, msg)
  call check(abs(cond(1, 1) - 1.0_dp) < 1.0e-12_dp, 'triad (1,1,0): P(Fe | Fe) = 1')
  do i = 1, 3; do j = 1, 3; ah(i, j) = s%z * (cond(i, j) - x(j)); end do; end do
  call screen_matrix(s, ah, w, gm=g3, wm=w3, ginv=gi3, lam=lam3, gam=gam3)
  call check(abs(maxval(lam3) - 12.0_dp) < 1.0e-10_dp .and. minval(gam3) < 1.0e-12_dp, 'triad (1,1,0): mode lambda = 12 has Gamma = 0')
  call drive_fixed_point(g3, w3, x, [0.10_dp, -0.30_dp, -0.06_dp], 2.0_dp, ok, dmax)
  call check(ok, 'drive iteration converges with a Gamma = 0 mode (F infinite)')
  call check_close(dmax, 0.0_dp, 1.0e-10_dp, 'drive fixed point: dq has no component on the Gamma = 0 mode')
  ! (c) the start of the SCF.  A potential without an SRO part has u = 0: then u_out = dq and v = W dq, at most a |dq| per
  ! mode (0 <= Gamma <= 2), and v follows the charges: dq = 0 gives v = 0.  The start u = Gamma^+ dq stands for an input
  ! potential that already holds F dq.  It keeps the first charges: with dq = 0 in the next iteration v = W (I - Gamma) u is
  ! still several Ry for the free-atom charges of the triad (1,1,0) in the first iteration.
  dq3 = [1.5062_dp, 0.5614_dp, -5.3635_dp]; dq3 = dq3 - sum(x * dq3)
  u3 = 0.0_dp
  call screen_update(g3, w3, x, [1, 1, 1], dq3, u3, uo3, v3)
  call check_close(maxval(abs(uo3 - dq3)), 0.0_dp, 1.0e-13_dp, 'drive start u = 0: u_out = dq')
  call check_close(maxval(abs(v3 - matmul(w3, dq3))), 0.0_dp, 1.0e-13_dp, 'drive start u = 0: v = W dq')
  call check(sqrt(sum(w * v3**2)) <= s%a * sqrt(sum(w * dq3**2)) * (1.0_dp + 1.0e-12_dp), 'drive start u = 0: |v| <= a |dq| (w norm)')
  call screen_update(g3, w3, x, [1, 1, 1], [0.0_dp, 0.0_dp, 0.0_dp], u3, uo3, v3)
  call check_close(maxval(abs(v3)), 0.0_dp, 0.0_dp, 'drive start u = 0: no charge deviation, no SRO potential')
  u3 = matmul(gi3, dq3)
  call screen_project(x, [1, 1, 1], u3)
  call screen_update(g3, w3, x, [1, 1, 1], [0.0_dp, 0.0_dp, 0.0_dp], u3, uo3, v3)
  call check(maxval(abs(v3)) > 2.0_dp, 'old start u = Gamma^+ dq: the first charges stay in the potential (> 2 Ry at dq = 0)')

  ! ---- basis-resolved screening on general structures (hcp, NiAs, intra-sublattice SRO)
  block
    use asr_crystal, only: structure_t
    use asr_prototypes, only: prototype_structure
    type(structure_t) :: st
    type(screen_t) :: sb, shcp, sn1, sn2
    complex(dp) :: mb1(1, 1, 1), mba(4, 4, 1), mbb(4, 4, 1), mbpm(4, 4, 2)
    real(dp) :: rvs(3, 2), grs(2), wrs(2), gr4(1), wr4(1), kv2(3, 2)
    p = reshape([-1, 1, 1, 1, -1, 1, 1, 1, -1], [3, 3]) * 0.5_dp * alat
    kv(:, 1) = [0.37_dp, -0.21_dp, 0.53_dp]
    call coulomb_k(p, kv, m)
    call coulomb_k_basis(p, reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1]), kv, mb1)
    call check_close(real(mb1(1, 1, 1), dp), m(1), 1.0e-12_dp, 'coulomb_k_basis matches coulomb_k on a one-site lattice')
    call check_close(aimag(mb1(1, 1, 1)), 0.0_dp, 1.0e-14_dp, 'coulomb_k_basis diagonal is real')

    call prototype_structure('nias', st, msg, c_over_a=1.7_dp)
    kv2(:, 1) = [0.31_dp, -0.19_dp, 0.27_dp]; kv2(:, 2) = -kv2(:, 1)
    call coulomb_k_basis(st%avec * alat, st%tau * alat, kv2(:, 1:1), mba, eta=0.35_dp)
    call coulomb_k_basis(st%avec * alat, st%tau * alat, kv2(:, 1:1), mbb, eta=0.60_dp)
    call check_close(maxval(abs(mba - mbb)) / maxval(abs(mba)), 0.0_dp, 1.0e-12_dp, &
                     'coulomb_k_basis off-diagonal independent of Ewald eta')
    call coulomb_k_basis(st%avec * alat, st%tau * alat, kv2, mbpm)
    call check_close(abs(mbpm(3, 3, 1) - mbpm(4, 4, 1)) / abs(mbpm(3, 3, 1)), 0.0_dp, 1.0e-12_dp, &
                     'NiAs M(k): M_33(k) = M_44(k)')
    call check_close(abs(abs(mbpm(1, 3, 1)) - abs(mbpm(1, 4, 2))) / abs(mbpm(1, 3, 1)), 0.0_dp, 1.0e-12_dp, &
                     'NiAs M(k): |M_13(k)| = |M_14(-k)|')
    call check_close(abs(abs(mbpm(2, 3, 1)) - abs(mbpm(2, 4, 2))) / abs(mbpm(2, 3, 1)), 0.0_dp, 1.0e-12_dp, &
                     'NiAs M(k): |M_23(k)| = |M_24(-k)|')

    call prototype_structure('b2', st, msg)
    call screen_structure(st, alat, 1.0_dp, sb, msg)
    call screen_lattice('bcc', alat, 1.0_dp, s2, msg)
    call check_close(abs(sb%a - s2%a) + abs(sb%g00 - s2%g00) + abs(sb%g1 - s2%g1), 0.0_dp, 1.0e-14_dp, &
                     'screen_structure on B2 delegates to bcc screen_lattice')

    call prototype_structure('hcp', st, msg)
    call screen_structure(st, alat, 1.0_dp, shcp, msg, nk=12)
    call check(len(msg) == 0 .and. shcp%z == 12 .and. abs(shcp%r1 - alat) < 1.0e-10_dp, 'hcp ideal: z = 12, R1 = a')
    call check_close(-shcp%w00 * shcp%r1 / (2.0_dp * shcp%g00), 1.0_dp, 1.0e-10_dp, 'hcp: beta round trip')
    call check_close(shcp%w00, shcp%a * (1.0_dp - shcp%g00), 1.0e-12_dp, 'hcp: W00 = a (1 - G00)')
    call check_close(sum(shcp%c1k**2) / real(size(shcp%c1k), dp), 1.0_dp / 12.0_dp, 1.0e-12_dp, &
                     'hcp branches: <c1_m^2> = 1 / z (no branch mixing)')
    call check_close(shcp%c1_ord, -1.0_dp / 3.0_dp, 1.0e-12_dp, 'hcp ideal ordering branch: c1_ord = -1/3')
    call check_close(screen_gamma(shcp, 0.0_dp), shcp%g00, 1.0e-14_dp, 'hcp closure: random limit Gamma(0) = G00')
    call check_close(screen_gamma(shcp, 12.0_dp), 0.0_dp, 1.0e-14_dp, 'hcp closure: complete clustering Gamma(z) = 0')
    call check_close(screen_gamma(shcp, -4.0_dp), shcp%g_ord, 1.0e-12_dp, 'hcp closure: ordering limit Gamma(z c1_ord) = g_ord')
    call check_positive_decreasing(shcp, 'hcp beta 1.0')

    call prototype_structure('nias', st, msg, c_over_a=1.6_dp)
    call screen_structure(st, alat, 1.0_dp, sn1, msg, nk=10, sub_sro=2)
    call check(len(msg) == 0 .and. sn1%z == 6, 'NiAs c/a 1.6 anion SRO shell: z = 6 out of plane')
    call check_close(sn1%r1, alat * sqrt(1.0_dp / 3.0_dp + 0.4_dp**2), 1.0e-10_dp, 'NiAs: R1_nn is the cation-anion bond')
    call check_close(-sn1%w00 * sn1%r1 / (2.0_dp * sn1%g00), 1.0_dp, 1.0e-10_dp, 'NiAs anion: beta round trip')
    call check_close(sum(sn1%c1k**2) / real(size(sn1%c1k), dp), 1.0_dp / 6.0_dp, 1.0e-12_dp, &
                     'NiAs c/a 1.6 anion branches: <c1_m^2> = 1 / z')
    call check_close(sn1%c1_ord, -1.0_dp, 1.0e-12_dp, 'NiAs c/a 1.6 bipartite anion shell: c1_ord = -1 at k = 0')
    call check_close(screen_gamma(sn1, 0.0_dp), sn1%g00, 1.0e-14_dp, 'NiAs c/a 1.6 closure: Gamma(0) = G00')
    call check_close(screen_gamma(sn1, 6.0_dp), 0.0_dp, 1.0e-14_dp, 'NiAs c/a 1.6 closure: Gamma(z) = 0')
    call check_close(screen_gamma(sn1, -6.0_dp), sn1%g_ord, 1.0e-12_dp, 'NiAs c/a 1.6 closure: Gamma(-z) = g_ord')
    ! At c/a = 1.6 < sqrt(8/3) the out-of-plane shell (0.9866 a) is nearly degenerate with the in-plane shell (1.0000 a).
    ! Alternating out-of-plane layers at k = 0 force like in-plane neighbours, so Gamma turns over below lambda = -2 while
    ! staying strictly positive on [-6, 6) and strictly decreasing on [-1, 6] (notes/DERIVATION_closure_basis.md).
    call check_range_positive_decreasing(sn1, -1.0_dp, 6.0_dp, 'nias c/a 1.6 [-1, z]')
    call check(screen_gamma(sn1, -6.0_dp) > 0.0_dp .and. screen_gamma(sn1, -6.0_dp) < screen_gamma(sn1, -2.0_dp), &
               'NiAs c/a 1.6: positive ordering turnover when in-plane shell is frustrated')
    block
      real(dp) :: rvo(3, 6), gro(6), wro(6)
      rvo(:, 1) = st%tau(:, 4) - st%tau(:, 3) - st%avec(:, 3)
      rvo(:, 2) = rvo(:, 1) + st%avec(:, 2)
      rvo(:, 3) = rvo(:, 1) - st%avec(:, 1)
      rvo(:, 4) = st%tau(:, 4) - st%tau(:, 3)
      rvo(:, 5) = rvo(:, 4) + st%avec(:, 2)
      rvo(:, 6) = rvo(:, 4) - st%avec(:, 1)
      rvo = rvo * alat
      call screen_structure_shells(st, alat, sn1, rvo, gro, wro, b0=3, b2=[4, 4, 4, 4, 4, 4], nk=10)
      call check_close(sum(gro) / 6.0_dp, sn1%g1, 1.0e-12_dp, 'NiAs c/a 1.6 inter-site shell: <G(R1_anion)> = G1')
    end block
    ! k = 0.  The G = 0 term of M(k) is the same for every pair of basis sites, so at k = 0 only the neutral block
    ! Q^T M(0) Q is defined (Q: an orthonormal basis of the vectors whose components sum to zero).  The neutral block does
    ! not depend on eta.  It is the limit of the neutral block of M(k) for k -> 0.  The ordering mode of the c/a 1.6 anion
    ! shell is v = (0, 0, 1, -1) / sqrt(2) at k = 0, so g_ord = a v^T Q (a + Q^T M(0) Q)^-1 Q^T v (solved here directly).
    block
      use asr_lapack, only: dgesv
      complex(dp) :: m0a(4, 4, 1), m0b(4, 4, 1), mk1(4, 4, 1), mk2(4, 4, 1)
      real(dp) :: q(4, 3), n0a(3, 3), n0b(3, 3), kz(3, 1), ks1(3, 1), ks2(3, 1), am(3, 3), qv(3), xs(3), d1, d2
      integer :: ipiv(3), info, jq
      q = 0.0_dp
      do jq = 1, 3
        q(1:jq, jq) = 1.0_dp / sqrt(real(jq * (jq + 1), dp)); q(jq + 1, jq) = -real(jq, dp) / sqrt(real(jq * (jq + 1), dp))
      end do
      kz = 0.0_dp
      call coulomb_k_basis(st%avec * alat, st%tau * alat, kz, m0a, eta=0.35_dp, gamma_neutral=.true.)
      call coulomb_k_basis(st%avec * alat, st%tau * alat, kz, m0b, eta=0.60_dp, gamma_neutral=.true.)
      n0a = matmul(transpose(q), matmul(real(m0a(:, :, 1), dp), q))
      n0b = matmul(transpose(q), matmul(real(m0b(:, :, 1), dp), q))
      call check_close(maxval(abs(n0a - n0b)) / maxval(abs(n0a)), 0.0_dp, 1.0e-12_dp, &
                       'NiAs M(0): neutral block independent of Ewald eta')
      call check_close(maxval(abs(aimag(m0a))), 0.0_dp, 1.0e-12_dp, 'NiAs M(0) is real')
      ks1(:, 1) = 1.0e-3_dp * [0.6_dp, -0.48_dp, 0.64_dp]; ks2 = 0.1_dp * ks1
      call coulomb_k_basis(st%avec * alat, st%tau * alat, ks1, mk1, eta=0.35_dp)
      call coulomb_k_basis(st%avec * alat, st%tau * alat, ks2, mk2, eta=0.35_dp)
      d1 = maxval(abs(matmul(transpose(q), matmul(mk1(:, :, 1), q)) - n0a))
      d2 = maxval(abs(matmul(transpose(q), matmul(mk2(:, :, 1), q)) - n0a))
      print '(2x,a,2es10.2)', 'NiAs neutral block |N(k) - N(0)| at |k| = 1e-3, 1e-4 (1/Bohr):', d1, d2
      ! Pairs of sites not related by inversion give terms linear in k, so the difference falls by 10 when |k| does.
      call check(abs(d1 / d2 - 10.0_dp) < 0.05_dp .and. d1 < 1.0e-2_dp * maxval(abs(n0a)), &
                 'NiAs M(k): neutral block is continuous at k = 0 (difference linear in |k|)')
      qv = matmul(transpose(q), [0.0_dp, 0.0_dp, 1.0_dp, -1.0_dp] / sqrt(2.0_dp))
      call coulomb_k_basis(st%avec * alat, st%tau * alat, kz, m0a, gamma_neutral=.true.)
      am = matmul(transpose(q), matmul(real(m0a(:, :, 1), dp), q))
      do jq = 1, 3; am(jq, jq) = am(jq, jq) + sn1%a; end do
      xs = qv
      call dgesv(3, 1, am, 3, ipiv, xs, 3, info)
      call check(info == 0, 'NiAs k = 0 neutral solve')
      print '(2x,a,es24.16,es10.2)', 'NiAs c/a 1.6 g_ord and |g_ord - direct solve|:', sn1%g_ord, &
        abs(sn1%g_ord - sn1%a * dot_product(qv, xs))
      call check_close(sn1%g_ord, sn1%a * dot_product(qv, xs), 1.0e-12_dp, &
                       'NiAs c/a 1.6: g_ord is the k = 0 neutral mode Green function (direct solve)')
    end block

    call prototype_structure('nias', st, msg, c_over_a=1.7_dp)
    call screen_structure(st, alat, 1.0_dp, sn2, msg, nk=10, sub_sro=2)
    call check(len(msg) == 0 .and. sn2%z == 6, 'NiAs c/a 1.7 anion SRO shell: z = 6 in plane')
    call check_close(sn2%c1_ord, -0.5_dp, 1.0e-12_dp, 'NiAs c/a 1.7 in-plane anion shell: c1_ord = -1/2 at K')
    ! The in-plane shell joins each anion site only to itself, so both branches have c1 = -1/2 at K and H.  g_ord is then the
    ! largest eigenvalue of the 2 x 2 anion block of G = a (a + M)^-1 at K or H, not a diagonal element.
    block
      use asr_lapack, only: zgesv
      complex(dp) :: mkh(4, 4, 2), gk4(4, 4), id4(4, 4)
      real(dp) :: fkh(3, 2), kkh(3, 2), bb(3, 3), pp(3, 3), gmx, gdiag
      integer :: ip4(4), info4, jk
      pp = st%avec * alat
      bb(:, 1) = [pp(2, 2) * pp(3, 3) - pp(3, 2) * pp(2, 3), pp(3, 2) * pp(1, 3) - pp(1, 2) * pp(3, 3), &
                  pp(1, 2) * pp(2, 3) - pp(2, 2) * pp(1, 3)]
      bb(:, 2) = [pp(2, 3) * pp(3, 1) - pp(3, 3) * pp(2, 1), pp(3, 3) * pp(1, 1) - pp(1, 3) * pp(3, 1), &
                  pp(1, 3) * pp(2, 1) - pp(2, 3) * pp(1, 1)]
      bb(:, 3) = [pp(2, 1) * pp(3, 2) - pp(3, 1) * pp(2, 2), pp(3, 1) * pp(1, 2) - pp(1, 1) * pp(3, 2), &
                  pp(1, 1) * pp(2, 2) - pp(2, 1) * pp(1, 2)]
      bb = 2.0_dp * pi * bb / dot_product(pp(:, 1), bb(:, 1))
      fkh(:, 1) = [1.0_dp, 1.0_dp, 0.0_dp] / 3.0_dp; fkh(:, 2) = [1.0_dp / 3.0_dp, 1.0_dp / 3.0_dp, 0.5_dp]
      kkh = matmul(bb, fkh)
      call coulomb_k_basis(pp, st%tau * alat, kkh, mkh)
      gmx = 0.0_dp; gdiag = 0.0_dp
      do jk = 1, 2
        id4 = (0.0_dp, 0.0_dp); do j = 1, 4; id4(j, j) = sn2%a; end do
        gk4 = mkh(:, :, jk); do j = 1, 4; gk4(j, j) = gk4(j, j) + sn2%a; end do
        call zgesv(4, 4, gk4, 4, ip4, id4, 4, info4)
        gmx = max(gmx, 0.5_dp * real(id4(3, 3) + id4(4, 4), dp) + &
                       sqrt(0.25_dp * real(id4(3, 3) - id4(4, 4), dp)**2 + abs(id4(3, 4))**2))
        gdiag = max(gdiag, real(id4(3, 3), dp), real(id4(4, 4), dp))
      end do
      print '(2x,a,3f12.8)', 'NiAs c/a 1.7: g_ord, largest eigenvalue of G_anion(K, H), largest diagonal:', sn2%g_ord, gmx, gdiag
      call check_close(sn2%g_ord, gmx, 1.0e-12_dp, &
                       'NiAs c/a 1.7: degenerate ordering branches use the largest eigenvalue of the projected G')
    end block
    ! Cell-choice independence of g_ord.  Simple hexagonal (c/a 1.7, in-plane shell) in a cell doubled along c has two sites
    ! stacked on top of each other.  Both branches have c1 = -1/2 at (K, kz') and G_12 /= 0 there.  The two branches are the
    ! primitive wavevectors (K, kz) and (K, kz + 1/2) folded together, so g_ord must be the largest a / (a + M(K, kz)) of the
    ! one-site lattice (coulomb_k) at kz = 0, 1/4, 1/2 (the folded special points kz' = 0, 1/2).
    block
      use asr_crystal, only: structure_init
      type(structure_t) :: st2, sth
      type(screen_t) :: sd
      real(dp) :: av2(3, 3), bp(3, 3), pp(3, 3), kp(3, 3), mp(3)
      call prototype_structure('hex', sth, msg, c_over_a=1.7_dp)
      av2 = sth%avec; av2(:, 3) = 2.0_dp * av2(:, 3)
      call structure_init(st2, av2, reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp], [3, 2]), [1, 1], msg)
      call check(len(msg) == 0, 'hex doubled along c: structure accepted')
      call screen_structure(st2, alat, 1.0_dp, sd, msg, nk=10)
      call check(len(msg) == 0 .and. sd%z == 6, 'hex doubled along c: in-plane shell, z = 6')
      pp = sth%avec * alat
      bp(:, 1) = [pp(2, 2) * pp(3, 3) - pp(3, 2) * pp(2, 3), pp(3, 2) * pp(1, 3) - pp(1, 2) * pp(3, 3), &
                  pp(1, 2) * pp(2, 3) - pp(2, 2) * pp(1, 3)]
      bp(:, 2) = [pp(2, 3) * pp(3, 1) - pp(3, 3) * pp(2, 1), pp(3, 3) * pp(1, 1) - pp(1, 3) * pp(3, 1), &
                  pp(1, 3) * pp(2, 1) - pp(2, 3) * pp(1, 1)]
      bp(:, 3) = [pp(2, 1) * pp(3, 2) - pp(3, 1) * pp(2, 2), pp(3, 1) * pp(1, 2) - pp(1, 1) * pp(3, 2), &
                  pp(1, 1) * pp(2, 2) - pp(2, 1) * pp(1, 2)]
      bp = 2.0_dp * pi * bp / dot_product(pp(:, 1), bp(:, 1))
      kp = matmul(bp, reshape([1.0_dp / 3.0_dp, 1.0_dp / 3.0_dp, 0.0_dp, 1.0_dp / 3.0_dp, 1.0_dp / 3.0_dp, 0.25_dp, &
                               1.0_dp / 3.0_dp, 1.0_dp / 3.0_dp, 0.5_dp], [3, 3]))
      call coulomb_k(pp, kp, mp)
      print '(2x,a,f12.8,a,3f12.8)', 'hex doubled: g_ord', sd%g_ord, ';  one-site a / (a + M(K, kz = 0, 1/4, 1/2)):', &
        sd%a / (sd%a + mp)
      call check_close(sd%c1_ord, -0.5_dp, 1.0e-12_dp, 'hex doubled along c: c1_ord = -1/2')
      call check_close(sd%g_ord, maxval(sd%a / (sd%a + mp)), 1.0e-12_dp, &
                       'hex doubled along c: g_ord = largest one-site G at the folded K points (cell-choice independence)')
    end block
    call check_positive_decreasing(sn2, 'nias c/a 1.7 anion beta 1.0')
    rvs(:, 1) = 0.0_dp; rvs(:, 2) = [alat, 0.0_dp, 0.0_dp]
    call screen_structure_shells(st, alat, sn2, rvs, grs, wrs, b0=3, b2=[3, 3], nk=10)
    call screen_structure_shells(st, alat, sn2, rvs(:, 1:1), gr4, wr4, b0=4, b2=[4], nk=10)
    call check_close(grs(1), gr4(1), 1.0e-12_dp, 'NiAs screen_structure_shells: G_33(0) = G_44(0)')
    call check_close(grs(1), sn2%g00, 1.0e-12_dp, 'NiAs screen_structure_shells: G(0) = G00 on anion sublattice')
    call check_close(grs(2), sn2%g1, 1.0e-12_dp, 'NiAs screen_structure_shells: G(R1_anion) = G1')
    call check_close(wrs(2), sn2%w1, 1.0e-12_dp, 'NiAs screen_structure_shells: W(R1_anion) = W1')
  end block
  call report()

contains

  function eye3() result(e)
    real(dp) :: e(3, 3)
    integer :: ii
    e = 0.0_dp; do ii = 1, 3; e(ii, ii) = 1.0_dp; end do
  end function eye3

  function inv3(a) result(b)
    real(dp), intent(in) :: a(3, 3)
    real(dp) :: b(3, 3), dt
    b(1, 1) = a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2); b(1, 2) = a(1, 3) * a(3, 2) - a(1, 2) * a(3, 3)
    b(1, 3) = a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2); b(2, 1) = a(2, 3) * a(3, 1) - a(2, 1) * a(3, 3)
    b(2, 2) = a(1, 1) * a(3, 3) - a(1, 3) * a(3, 1); b(2, 3) = a(1, 3) * a(2, 1) - a(1, 1) * a(2, 3)
    b(3, 1) = a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1); b(3, 2) = a(1, 2) * a(3, 1) - a(1, 1) * a(3, 2)
    b(3, 3) = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)
    dt = a(1, 1) * b(1, 1) + a(1, 2) * b(2, 1) + a(1, 3) * b(3, 1)
    b = b / dt
  end function inv3

  !> Gamma(lambda) > 0 and non-increasing on lambda in [z c1(k_o), z) (the realizable range of one sublattice), 0 at z.
  subroutine check_positive_decreasing(sc, label)
    type(screen_t), intent(in) :: sc
    character(len=*), intent(in) :: label
    call check_range_positive_decreasing(sc, sc%z * sc%c1_ord, real(sc%z, dp), label)
  end subroutine check_positive_decreasing

  subroutine check_range_positive_decreasing(sc, lmin, lmax, label)
    type(screen_t), intent(in) :: sc
    real(dp), intent(in) :: lmin, lmax
    character(len=*), intent(in) :: label
    logical :: pos, dec
    pos = .true.; dec = .true.; gprev = huge(1.0_dp)
    do k = 0, 400
      lam = lmin + (lmax - lmin) * k / 400.0_dp
      gnow = screen_gamma(sc, lam)
      if (lam < real(sc%z, dp) - 1.0e-12_dp .and. gnow <= 0.0_dp) pos = .false.
      if (gnow > gprev + 1.0e-12_dp) dec = .false.
      gprev = gnow
    end do
    call check(pos, 'closure '//label//': Gamma > 0 on the realizable range (no pole)')
    call check(dec, 'closure '//label//': Gamma decreases with lambda')
  end subroutine check_range_positive_decreasing

  !> Damped iteration of u with dq = P (c - chi W u) (the potential implied by u acts on a linear model response), then
  !> u_out = screen_update.  ok: converged.  dmax: |v - F dq| if fref is given, else |Gamma Gamma^+ dq - dq| (dq in the range
  !> of Gamma, i.e. no charge on a Gamma = 0 mode).
  subroutine drive_fixed_point(gm, wm, xc, c0, chi, ok, dmax, fref)
    real(dp), intent(in) :: gm(3, 3), wm(3, 3), xc(3), c0(3), chi
    logical, intent(out) :: ok
    real(dp), intent(out) :: dmax
    real(dp), intent(in), optional :: fref(3, 3)
    real(dp) :: u(3), uo(3), v(3), c(3), dq(3), gi(3, 3)
    integer :: it, grp(3)
    grp = 1; c = c0
    call screen_project(xc, grp, c)
    u = 0.0_dp; ok = .false.
    do it = 1, 2000
      dq = c - chi * matmul(wm, u)
      call screen_project(xc, grp, dq)
      call screen_update(gm, wm, xc, grp, dq, u, uo, v)
      if (maxval(abs(uo - u)) < 1.0e-14_dp) then; ok = .true.; exit; end if
      u = u + 0.3_dp * (uo - u)
    end do
    if (present(fref)) then
      dmax = maxval(abs(v - matmul(fref, dq)))
    else
      call screen_matrix(s, ah, w, ginv=gi)
      dmax = maxval(abs(matmul(gm, matmul(gi, dq)) - dq))
    end if
  end subroutine drive_fixed_point
end program test_screening
