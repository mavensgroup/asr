!> Block kernels asr_augblk against the single-energy kernels of asr_augcls, asr_edc and asr_ed3 (speed item S1/S2):
!>  (K) blk_hp == augcls_hp_c, blk_rand_lin == augcls_rand_c with X(z) = z A + B, on a reduced fcc operator (D4h group) with a
!>      conditional shell class;
!>  (E) edc_green_blk == edc_green and ed3_blocks_blk == ed3_blocks for a block of energies, on the ring of 5 (four species,
!>      classes, a conditional shell with species of zero weight, frozen start site) and on the reduced fcc operator, without and
!>      with the adaptive depth test (members then leave the block at different depths), and a short block padded with idle
!>      members.  The ring stays below the step count at
!>      which its small augmented space is exhausted: beyond it the chains carry rounding noise only, and the depth test then
!>      stops single and block chains at different steps.  The fcc chains (40 steps) amplify the rounding-level kernel
!>      difference to about 1e-10.
!> Timings are printed for information only (the laptop runs one thread; timings belong on the desktop).
program test_augcls_blk
  use asr_kinds, only: dp
  use asr_augspace, only: graph_t
  use asr_augcls
  use asr_augblk
  use asr_edc
  use asr_ed3
  use asr_augorb, only: orb_alloy_t
  use asr_endep, only: ed_depth_tol
  use asr_lmto_alloy, only: build_lmto_alloy
  use asr_symmetry
  use asr_lmto_ham, only: potpar_t
  use asr_testing
  implicit none
  integer, parameter :: l0s(4) = [1, 3, 5, 7], nz = 7
  type(potpar_t) :: pp(4)
  real(dp) :: alpha(0:2)
  complex(dp) :: zs(nz)
  integer :: l, sp, j

  print '(a)', 'test_augcls_blk'
  do l = 0, 2
    pp(1)%enu(l) = 0.1_dp * l; pp(1)%c(l) = 0.3_dp - 0.15_dp * l; pp(1)%sqrdel(l) = 0.30_dp - 0.07_dp * l
    pp(1)%p(l) = 0.9_dp - 0.1_dp * l; pp(1)%gamma(l) = 0.25_dp - 0.1_dp * l
    pp(2)%enu(l) = 0.05_dp * l; pp(2)%c(l) = 0.45_dp - 0.1_dp * l; pp(2)%sqrdel(l) = 0.27_dp - 0.06_dp * l
    pp(2)%p(l) = 0.8_dp - 0.1_dp * l; pp(2)%gamma(l) = 0.2_dp - 0.08_dp * l
    pp(3)%enu(l) = -0.07_dp + 0.08_dp * l; pp(3)%c(l) = 0.2_dp - 0.05_dp * l; pp(3)%sqrdel(l) = 0.22_dp - 0.04_dp * l
    pp(3)%p(l) = 0.7_dp - 0.1_dp * l; pp(3)%gamma(l) = 0.15_dp - 0.06_dp * l
    pp(4)%enu(l) = 0.02_dp + 0.03_dp * l; pp(4)%c(l) = 0.6_dp - 0.2_dp * l; pp(4)%sqrdel(l) = 0.25_dp - 0.05_dp * l
    pp(4)%p(l) = 0.85_dp - 0.1_dp * l; pp(4)%gamma(l) = 0.22_dp - 0.07_dp * l
  end do
  do sp = 1, 4; pp(sp)%alpha_rep = pp(sp)%gamma; end do
  alpha = [0.3485_dp, 0.05303_dp, 0.010714_dp]
  do j = 1, nz                                            ! from far above the axis to close to it (different depths)
    zs(j) = cmplx(-0.4_dp + 0.13_dp * j, 0.8_dp / real(j, dp)**1.5_dp, dp)
  end do
  call ring_case()
  call fcc_case()
  call report()

contains

  subroutine ring_case()
    integer, parameter :: n = 5, rcls(n) = [3, 4, 1, 2, 5], nlev = 24
    real(dp), parameter :: c4(3, 3) = reshape([0.0_dp, 0.7_dp, 0.3_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.2_dp, 0.6_dp, 0.2_dp], [3, 3])
    real(dp), parameter :: c5(2, 3) = reshape([0.1_dp, 0.9_dp, 1.0_dp, 0.0_dp, 0.55_dp, 0.45_dp], [2, 3])
    type(graph_t) :: g
    type(class_alloy_t) :: mdl
    type(augcls_op) :: op
    integer :: perm(n, 1), i, jc
    real(dp) :: dmat(9, 9, 1)
    logical :: fl(n)
    allocate (g%nbr_ptr(n + 1), g%nbr(2 * n)); g%n = n
    do i = 1, n
      g%nbr_ptr(i) = 2 * (i - 1) + 1
      g%nbr(2 * i - 1) = mod(i, n) + 1
      g%nbr(2 * i) = mod(i - 2 + n, n) + 1
    end do
    g%nbr_ptr(n + 1) = 2 * n + 1
    mdl%nl = 9; mdl%nsp = 4; mdl%ncls = 5
    allocate (mdl%nsc(5), mdl%spc(3, 5), mdl%conc(3, 5)); mdl%spc = 0; mdl%conc = 0.0_dp
    mdl%nsc = [1, 2, 3, 3, 2]
    mdl%spc(1:1, 1) = [1];       mdl%conc(1:1, 1) = [1.0_dp]
    mdl%spc(1:2, 2) = [2, 3];    mdl%conc(1:2, 2) = [0.4_dp, 0.6_dp]
    mdl%spc(1:3, 3) = [1, 3, 4]; mdl%conc(1:3, 3) = [0.3_dp, 0.5_dp, 0.2_dp]
    mdl%spc(:, 4) = mdl%spc(:, 3); mdl%conc(:, 4) = mdl%conc(:, 3)
    mdl%spc(:, 5) = mdl%spc(:, 2); mdl%conc(:, 5) = mdl%conc(:, 2)
    call random_blocks(mdl)
    perm(:, 1) = [(i, i=1, n)]
    dmat = 0.0_dp; do l = 1, 9; dmat(l, l, 1) = 1.0_dp; end do
    fl = .true.; fl(1) = .false.
    do jc = 1, 3
      call build_augcls_operator(g, mdl, rcls, perm, dmat, 1, 80, op, flippable=fl)
      call compare('ring', mdl, op, mdl%spc(jc, 3), dmat, [1, 5], nlev, c4(:, jc), c5(:, jc), tol=1.0e-11_dp)
    end do
  end subroutine ring_case

  subroutine fcc_case()
    integer, parameter :: nlev = 40
    real(dp), parameter :: cs(3) = [0.1_dp, 0.7_dp, 0.2_dp]
    type(graph_t) :: g
    type(orb_alloy_t) :: m0
    type(class_alloy_t) :: mdl
    type(augcls_op) :: op
    integer, allocatable :: coord(:, :), cls(:), perm(:, :), mats(:, :, :)
    real(dp), allocatable :: d(:, :, :)
    logical, allocatable :: fl(:)
    integer :: i, r2, r2min
    call build_lmto_alloy('fcc', 4, 1.05_dp, 0.7_dp, pp(1), pp(2), 0.5_dp, g, coord, m0, alpha_in=alpha)
    mdl%nl = 9; mdl%nsp = 4; mdl%ncls = 2
    allocate (mdl%nsc(2), mdl%spc(3, 2), mdl%conc(3, 2))
    mdl%nsc = 3; mdl%spc(:, 1) = [2, 3, 4]; mdl%spc(:, 2) = [2, 3, 4]
    mdl%conc(:, 1) = [0.3_dp, 0.5_dp, 0.2_dp]; mdl%conc(:, 2) = mdl%conc(:, 1)
    mdl%block = m0%block; mdl%block0 = m0%block0
    allocate (fl(g%n), cls(g%n)); fl = .true.; fl(1) = .false.; cls = 1
    r2min = huge(1)
    do i = 2, g%n; r2min = min(r2min, sum((modulo(coord(:, i) + 4, 8) - 4)**2)); end do
    do i = 2, g%n
      r2 = sum((modulo(coord(:, i) + 4, 8) - 4)**2)
      if (r2 == r2min) cls(i) = 2
    end do
    call class_ops(coord, 8, cls, perm, mats)
    call orbital_rep(mats, d)
    call build_augcls_operator(g, mdl, cls, perm, d, 1, 3, op, maxrank=1, flippable=fl)
    print '(a,i0,a,i0)', '  fcc: ', op%nn, ' reduced nodes, group order ', size(perm, 2)
    call kernels(mdl, op, d)
    call compare('fcc', mdl, op, 3, d, l0s, nlev, cs, tol=1.0e-9_dp)
  end subroutine fcc_case

  !> (K) kernels on one operator
  subroutine kernels(mdl, op, d)
    type(class_alloy_t), intent(inout) :: mdl
    type(augcls_op), intent(inout) :: op
    real(dp), intent(in) :: d(:, :, :)
    type(ed3_t) :: ed
    type(blk_work_t) :: bw
    real(dp) :: pa(9, 4), pb(9, 4), u, eh, er
    complex(dp), allocatable :: x(:, :), y(:, :), yr(:)
    complex(dp) :: par(9, 4)
    integer :: i, j
    call ed3_setup(mdl, pp, alpha, 3, ed)
    call augcls_set_params(op, mdl); call augcls_freeze_site(op, 1, 3)
    call augcls_set_class_conc(op, 2, [0.1_dp, 0.7_dp, 0.2_dp])
    call augcls_set_start(op, 5, d)
    allocate (op%wk(op%dim(), 8))
    call blk_work_init(bw, 3, op%dim())
    allocate (x(3, op%dim()), y(3, op%dim()), yr(op%dim()))
    do i = 1, op%dim(); do j = 1, 3
      call random_number(u); x(j, i) = cmplx(u - 0.5_dp, 0.3_dp * j * (u - 0.2_dp), dp)
    end do; end do
    do i = 1, 4; do l = 1, 9
      call random_number(u); pa(l, i) = u - 0.3_dp; call random_number(u); pb(l, i) = 0.5_dp * u
    end do; end do
    call blk_hp(op, bw, x, y)
    eh = 0.0_dp
    do j = 1, 3
      call augcls_hp_c(op, x(j, :), yr); eh = max(eh, maxval(abs(y(j, :) - yr)) / maxval(abs(yr)))
    end do
    call blk_rand_lin(op, bw, pa, pb, zs(1:3), x, y)
    er = 0.0_dp
    do j = 1, 3
      par = zs(j) * pa + pb
      call augcls_rand_c(op, par, x(j, :), yr); er = max(er, maxval(abs(y(j, :) - yr)) / maxval(abs(yr)))
    end do
    print '(a,2es10.2)', '  fcc kernels: blk_hp, blk_rand_lin relative difference ', eh, er
    call check_close(eh, 0.0_dp, 1.0e-13_dp, 'blk_hp == augcls_hp_c (reduced fcc, shell class)')
    call check_close(er, 0.0_dp, 1.0e-13_dp, 'blk_rand_lin == augcls_rand_c with X(z) = z A + B')
    deallocate (op%wk)
  end subroutine kernels

  !> (E) chains: block against single energies, both routes, without and with the depth test
  subroutine compare(tag, mdl, op, cen, d, starts, nlev, cshell, cshell2, tol)
    character(*), intent(in) :: tag
    type(class_alloy_t), intent(inout) :: mdl
    type(augcls_op), intent(inout) :: op
    integer, intent(in) :: cen, starts(:), nlev
    real(dp), intent(in) :: d(:, :, :), cshell(:)
    real(dp), intent(in), optional :: cshell2(:)
    real(dp), intent(in) :: tol
    type(edn_t) :: en
    type(edn_work_t) :: wn
    type(ed3_t) :: e3
    type(ed3_work_t) :: w3
    type(blk_work_t) :: bw
    complex(dp) :: g1, gb(nz), gs(nz), bk(0:3), bkb(0:3, nz)
    integer :: io, j, idt, mb(nz), m1, ksh
    real(dp) :: ec, e3d, ts, tb
    integer(8) :: c0, c1, cr
    ksh = merge(4, 2, present(cshell2))
    do idt = 1, 2
      ed_depth_tol = merge(0.0_dp, 1.0e-3_dp, idt == 1)      ! 1e-3: members leave the block at different depths
      ec = 0.0_dp; e3d = 0.0_dp; ts = 0.0_dp; tb = 0.0_dp
      ! edc
      call edc_setup(mdl, pp, alpha, cen, en)
      call augcls_set_params(op, mdl); call augcls_freeze_site(op, 1, cen)
      call augcls_set_class_conc(op, ksh, cshell)
      if (present(cshell2)) call augcls_set_class_conc(op, 5, cshell2)
      if (associated(op%wk)) deallocate (op%wk)
      allocate (op%wk(op%dim(), 8))
      call edn_work_init(wn, op%dim()); call blk_work_init(bw, nz, op%dim())
      do io = 1, size(starts)
        call augcls_set_start(op, starts(io), d)
        call system_clock(c0, cr)
        call edc_green_blk(op, bw, en, zs, starts(io), nlev, gb, mb)
        call system_clock(c1); tb = tb + real(c1 - c0, dp) / real(cr, dp)
        call edc_green_blk(op, bw, en, zs(1:nz - 2), starts(io), nlev, gs(1:nz - 2))      ! short block (idle padding)
        ec = max(ec, maxval(abs(gs(1:nz - 2) - gb(1:nz - 2))))
        do j = 1, nz
          call system_clock(c0)
          call edc_green(op, wn, en, zs(j), starts(io), nlev, g1, mused=m1)
          call system_clock(c1); ts = ts + real(c1 - c0, dp) / real(cr, dp)
          ec = max(ec, abs(gb(j) - g1) / abs(g1))
          if (mb(j) /= m1) ec = 1.0_dp
        end do
      end do
      ! ed3
      call ed3_setup(mdl, pp, alpha, cen, e3)
      call augcls_set_params(op, mdl); call augcls_freeze_site(op, 1, cen)
      call augcls_set_class_conc(op, ksh, cshell)
      if (present(cshell2)) call augcls_set_class_conc(op, 5, cshell2)
      call ed3_work_init(w3, op%dim())
      do io = 1, size(starts)
        call augcls_set_start(op, starts(io), d)
        call ed3_start(op, w3, e3, starts(io))
        call ed3_blocks_blk(op, w3, bw, e3, zs, starts(io), nlev, bkb)
        do j = 1, nz
          call ed3_blocks(op, w3, e3, zs(j), starts(io), nlev, bk)
          e3d = max(e3d, maxval(abs(bkb(0:2, j) - bk(0:2))) / maxval(abs(bk(0:2))))
        end do
      end do
      print '(a,a,a,es8.1,a,2es10.2,a,2f8.3,a)', '  ', tag, ': depth_tol ', ed_depth_tol, '  edc, ed3 relative difference ', ec, e3d, &
        '   edc time single, block ', ts, tb, ' s'
      call check_close(ec, 0.0_dp, tol, tag // ': edc_green_blk == edc_green (same depths)')
      call check_close(e3d, 0.0_dp, tol, tag // ': ed3_blocks_blk == ed3_blocks')
    end do
    ed_depth_tol = 0.0_dp
  end subroutine compare

  subroutine random_blocks(mdl)
    type(class_alloy_t), intent(inout) :: mdl
    real(dp) :: u
    integer :: i, k
    allocate (mdl%block(9, 9, 2), mdl%block0(9, 9))
    do i = 1, 9; do k = 1, 9
      call random_number(u); mdl%block(i, k, 1) = 0.5_dp * (u - 0.5_dp)
    end do; end do
    mdl%block(:, :, 2) = transpose(mdl%block(:, :, 1))
    do i = 1, 9; do k = i, 9
      call random_number(u); mdl%block0(i, k) = 0.6_dp * (u - 0.5_dp); mdl%block0(k, i) = mdl%block0(i, k)
    end do; end do
  end subroutine random_blocks
end program test_augcls_blk
