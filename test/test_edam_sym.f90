!> AMscf-form energy-dependent ASR: the D4h-reduced operator (asr_augsym) vs the unreduced one (asr_augorb), both with the
!> start site frozen to the central species: identical chain coefficients and Green function for the s, p_z, t2g, e_g start orbitals
!> (fcc), both central species, complex energies.
program test_edam_sym
  use asr_kinds, only: dp
  use asr_lmto_ham, only: potpar_t
  use asr_augspace, only: graph_t
  use asr_augorb
  use asr_augsym
  use asr_endep, only: ed_work_t, ed_work_init
  use asr_edam
  use asr_lmto_alloy
  use asr_symmetry
  use asr_recursion, only: lanczos
  use asr_testing
  implicit none
  type(potpar_t) :: pa, pb
  type(graph_t) :: g
  type(orb_alloy_t) :: mdl
  type(augorb_op) :: op
  type(augsym_op) :: sop
  type(edam_t) :: ed
  type(ed_work_t) :: w1, w2
  integer, allocatable :: coord(:, :), perm(:, :), mats(:, :, :)
  real(dp), allocatable :: d(:, :, :)
  logical, allocatable :: flippable(:)
  complex(dp) :: z, g1, g2, a1(8), b1(8), a2(8), b2(8)
  real(dp) :: alpha(0:2), err, errg
  integer :: l, io, m1, m2, cen, kz, m3
  real(dp), allocatable :: v0(:), a3(:), b3(:)
  real(dp) :: err0
  type(augsym_op) :: sei, sed
  integer, parameter :: l0(4) = [1, 3, 5, 7], nlev = 3

  print '(a)', 'test_edam_sym'
  do l = 0, 2
    pa%enu(l) = 0.1_dp * l; pa%c(l) = 0.3_dp - 0.15_dp * l; pa%sqrdel(l) = 0.12_dp + 0.05_dp * l
    pa%p(l) = 0.9_dp - 0.1_dp * l; pa%gamma(l) = 0.2_dp; pa%alpha_rep(l) = 0.2_dp
    pb%enu(l) = 0.05_dp * l; pb%c(l) = 0.45_dp - 0.1_dp * l; pb%sqrdel(l) = 0.1_dp + 0.04_dp * l
    pb%p(l) = 0.8_dp - 0.1_dp * l; pb%gamma(l) = 0.2_dp; pb%alpha_rep(l) = 0.2_dp
  end do
  alpha = [0.3485_dp, 0.05303_dp, 0.010714_dp]
  call build_lmto_alloy('fcc', 6, 1.2_dp, 0.7_dp, pa, pb, 0.4_dp, g, coord, mdl, alpha_in=alpha)
  call d4h_ops(coord, 12, perm, mats)
  call orbital_rep(mats, d)
  allocate (flippable(g%n)); flippable = .true.; flippable(1) = .false.
  err = 0.0_dp; errg = 0.0_dp
  do cen = 1, 2
    call edam_setup(mdl, pa, pb, alpha, cen == 1, ed)
    call build_augorb_operator(g, mdl, 1, nlev, op, maxrank=1, flippable=flippable)
    call build_augsym_operator(g, mdl, perm, d, 1, nlev, sop, maxrank=1, flippable=flippable)
    call augorb_freeze_site(op, 1, ed%mc)
    call augsym_freeze_site(sop, 1, ed%mc)
    print '(a,i0,a,i0,a,i0)', '  central species ', cen, ': nodes unreduced ', op%nn, '  reduced ', sop%nn
    call ed_work_init(w1, op%dim()); call ed_work_init(w2, sop%dim())
    do io = 1, 4
      call augsym_set_start(sop, l0(io), d)
      do kz = 1, 3
        z = cmplx(-0.2_dp + 0.3_dp * kz, 0.15_dp, dp)
        call edam_green(op, w1, ed, z, l0(io), nlev, g1, a1, b1, m1)
        call edam_green(sop, w2, ed, z, l0(io), nlev, g2, a2, b2, m2)
        err = max(err, maxval(abs(a1(1:m1) - a2(1:m2))), maxval(abs(b1(1:m1 - 1) - b2(1:m2 - 1))))
        errg = max(errg, abs(g1 - g2))
      end do
    end do
  end do
  print '(a,es9.2,a,es9.2)', '  max |coefficient difference| ', err, '   max |G difference| ', errg
  ! ---- ordered limit A = B: T = 1, H(E) = H^(2) at every E => the ED chain equals the energy-independent chain (depth 8)
  call build_lmto_alloy('fcc', 6, 1.2_dp, 0.7_dp, pa, pa, 0.4_dp, g, coord, mdl, alpha_in=alpha)
  call build_augsym_operator(g, mdl, perm, d, 1, 8, sei, 1, order=2)
  flippable = .true.; flippable(1) = .false.
  call edam_setup(mdl, pa, pa, alpha, .true., ed)
  call build_augsym_operator(g, mdl, perm, d, 1, 8, sed, 1, flippable, order=2)
  call augsym_freeze_site(sed, 1, ed%mc)
  call ed_work_init(w2, sed%dim())
  allocate (v0(sei%dim()), a3(8), b3(8)); err0 = 0.0_dp
  do io = 1, 4
    call augsym_set_start(sei, l0(io), d); call augsym_set_start(sed, l0(io), d)
    v0 = 0.0_dp; v0(l0(io)) = 1.0_dp
    call lanczos(sei, v0, 8, a3, b3, m3, reorth=.true.)
    z = cmplx(0.1_dp, 0.2_dp, dp)
    call edam_green(sed, w2, ed, z, l0(io), 8, g2, a2, b2, m2)
    err0 = max(err0, maxval(abs(a2(1:8) - a3(1:8))), maxval(abs(b2(1:7) - b3(1:7))))
  end do
  print '(a,es9.2)', '  ordered limit: max |ED chain - H^(2) chain| ', err0
  call check_close(err0, 0.0_dp, 1.0e-9_dp, 'ordered limit (A = B): AMscf-form chain == H^(2) chain')
  call check_close(err, 0.0_dp, 1.0e-9_dp, 'reduced == unreduced chain coefficients (AMscf-form ED)')
  call check_close(errg, 0.0_dp, 1.0e-9_dp, 'reduced == unreduced Green function')
  call report()
end program test_edam_sym
