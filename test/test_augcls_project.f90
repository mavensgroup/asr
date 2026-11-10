!> Symmetric subspace of the group-reduced class operator (asr_augcls) and the chains on it.
!> A reduced vector is valid only if psi(r) = chi(g) D(g)^T psi(r) for every g in the stabiliser of the node r.  The operator maps
!> this subspace V into itself and is symmetric on it; on the rest of the coordinate space it is another, non-symmetric matrix.
!> Rounding takes a Lanczos vector out of V, and the recurrence amplifies that part (ordered B2 cell: 1e-16 at 12 steps, 1e-11 at
!> 24, 1e-5 at 36, percent at 48).  The chains therefore project every new vector (augcls_project_c, blk_project).
!> Ordered B2 on a bcc supercell (reduced by the site group, 216 coordinates) against a dense solve on the unreduced operator
!> (1152 coordinates, exactly symmetric):
!>  (P) the projector is idempotent and symmetric, the start vectors lie in V, and the operator keeps V;
!>  (C) long chains of edc and ed3 equal the dense solve; without the projection ed3 stalls four or more orders above it;
!>  (B) block chains equal single-energy chains at the same length (the 1e-4 difference of B2 ed3 runs with energy_block).
program test_augcls_project
  use asr_kinds, only: dp
  use asr_augspace, only: graph_t
  use asr_augcls
  use asr_augblk
  use asr_edc
  use asr_ed3
  use asr_augorb, only: orb_alloy_t
  use asr_endep, only: ed_depth_tol, ed_tail, ed_project
  use asr_lmto_alloy, only: build_lmto_alloy
  use asr_symmetry
  use asr_lmto_ham, only: potpar_t
  use asr_testing
  implicit none
  integer, parameter :: lsup = 4, nlev = 144, iv = 5
  type(potpar_t) :: pp(2)
  real(dp) :: alpha(0:2), d1(9, 9, 1), u1, e_on(2), e_off(2), eb(2), oc
  type(graph_t) :: g
  type(orb_alloy_t) :: m0
  type(class_alloy_t) :: mdl
  type(augcls_op) :: op, opu
  type(edn_t) :: en
  type(edn_work_t) :: wn
  type(ed3_t) :: e3
  type(ed3_work_t) :: w3, w3u
  type(blk_work_t) :: bw
  integer, allocatable :: coord(:, :), cls(:), perm(:, :), mats(:, :, :), perm1(:, :), ipiv(:)
  real(dp), allocatable :: d(:, :, :), xr(:), yr(:), ur(:), hr(:)
  logical, allocatable :: fl(:)
  complex(dp), allocatable :: km(:, :), ev(:), col(:), rhs(:, :), x(:), y(:), px(:), py(:), ep(:, :)
  complex(dp) :: z, zs(2), gex, q0ex, q2ex, g1, gb(2), bk(0:3), bkb(0:3, 2), et(9, 2), uu(9, 2), qt(9, 2)
  integer :: l, i, nd, ndu, info, iz
  external :: zgesv

  print '(a)', 'test_augcls_project'
  do l = 0, 2
    pp(1)%enu(l) = 0.1_dp * l; pp(1)%c(l) = 0.3_dp - 0.15_dp * l; pp(1)%sqrdel(l) = 0.30_dp - 0.07_dp * l
    pp(1)%p(l) = 0.9_dp - 0.1_dp * l; pp(1)%gamma(l) = 0.25_dp - 0.1_dp * l
    pp(2)%enu(l) = 0.05_dp * l; pp(2)%c(l) = 0.45_dp - 0.1_dp * l; pp(2)%sqrdel(l) = 0.27_dp - 0.06_dp * l
    pp(2)%p(l) = 0.8_dp - 0.1_dp * l; pp(2)%gamma(l) = 0.2_dp - 0.08_dp * l
  end do
  pp(1)%alpha_rep = pp(1)%gamma; pp(2)%alpha_rep = pp(2)%gamma
  alpha = [0.3485_dp, 0.05303_dp, 0.010714_dp]
  call build_lmto_alloy('bcc', lsup, 1.2_dp, 0.7_dp, pp(1), pp(2), 0.5_dp, g, coord, m0, alpha_in=alpha)
  mdl%nl = 9; mdl%nsp = 2; mdl%ncls = 2
  allocate (mdl%nsc(2), mdl%spc(1, 2), mdl%conc(1, 2))
  mdl%nsc = 1; mdl%spc(1, 1) = 1; mdl%spc(1, 2) = 2; mdl%conc = 1.0_dp
  mdl%block = m0%block; mdl%block0 = m0%block0
  allocate (fl(g%n), cls(g%n)); fl = .false.
  do i = 1, g%n; cls(i) = 1 + modulo(coord(1, i), 2); end do
  call class_ops(coord, 2 * lsup, cls, perm, mats)
  call orbital_rep(mats, d)
  allocate (perm1(g%n, 1)); perm1(:, 1) = [(i, i=1, g%n)]
  d1 = 0.0_dp; do l = 1, 9; d1(l, l, 1) = 1.0_dp; end do
  allocate (mdl%c(9, 2), mdl%d(9, 2), mdl%o(9, 2), mdl%e(9, 2))
  call ed3_setup(mdl, pp, alpha, 1, e3)
  call build_augcls_operator(g, mdl, cls, perm, d, 1, 6, op, maxrank=1, flippable=fl)
  call build_augcls_operator(g, mdl, cls, perm1, d1, 1, 6, opu, maxrank=1, flippable=fl)
  call augcls_set_start(op, iv, d); call augcls_set_start(opu, iv, d1)
  nd = op%dim(); ndu = opu%dim()
  print '(a,i0,a,i0,a,i0,a,i0)', '  B2: sites ', g%n, ', reduced nodes ', op%nn, ' (', count(op%stab_ptr(2:) - op%stab_ptr(:op%nn) > 1), &
    ' with a stabiliser), group order ', size(perm, 2)
  call check(all(opu%stab_ptr(2:) - opu%stab_ptr(:opu%nn) == 1), 'unreduced operator: every stabiliser is the identity')
  call ed3_work_init(w3, nd); call ed3_work_init(w3u, ndu); call edn_work_init(wn, nd); call blk_work_init(bw, 2, nd)
  ed_depth_tol = 0.0_dp; ed_tail = 'none'; ed_project = .true.

  ! ---- (P) projector
  call augcls_set_params(op, mdl); call augcls_freeze_site(op, 1, 1)
  call augcls_set_params(opu, mdl); call augcls_freeze_site(opu, 1, 1)
  allocate (x(nd), y(nd), xr(nd), yr(nd))
  call lcg(xr); call lcg(yr); x = cmplx(xr, 0.0_dp, dp); call lcg(xr); x = x + cmplx(0.0_dp, xr, dp)
  y = cmplx(yr, 0.0_dp, dp); call lcg(yr); y = y + cmplx(0.0_dp, yr, dp)
  px = x; call augcls_project_c(op, px); py = y; call augcls_project_c(op, py)
  call check(augcls_leak_c(op, x) > 0.1_dp, 'a generic vector is not in the symmetric subspace')
  call check_close(augcls_leak_c(op, px), 0.0_dp, 1.0e-14_dp, 'P is idempotent')
  call check_close(abs(sum(y * px) - sum(py * x)) / abs(sum(py * x)), 0.0_dp, 1.0e-13_dp, 'P is symmetric in the bilinear form')
  z = cmplx(0.15_dp, 0.03_dp, dp)
  call ed3_coefs(e3, z, et, uu, qt)
  call ed3_kmul(op, w3, et, uu, qt, px, y)
  call check_close(augcls_leak_c(op, y), 0.0_dp, 1.0e-13_dp, 'K'' keeps the symmetric subspace')
  call ed3_kmul(op, w3, et, uu, qt, py, x)
  call check_close(abs(sum(py * y) - sum(px * x)) / abs(sum(px * x)), 0.0_dp, 1.0e-12_dp, 'K'' is symmetric on the symmetric subspace')
  call ed3_kmul(op, w3, et, uu, qt, cmplx(xr, yr, dp), y); call ed3_kmul(op, w3, et, uu, qt, cmplx(yr, xr, dp), x)
  u1 = abs(sum(cmplx(yr, xr, dp) * y) - sum(cmplx(xr, yr, dp) * x)) / abs(sum(cmplx(xr, yr, dp) * x))
  print '(a,es9.2)', '  K'' outside the subspace: relative asymmetry of a generic matrix element ', u1
  call check(u1 > 1.0e-3_dp, 'K'' is not symmetric outside the symmetric subspace')
  call ed3_start(op, w3, e3, iv); call ed3_start(opu, w3u, e3, iv)
  call check_close(max(augcls_leak_c(op, w3%av), augcls_leak_c(op, w3%vv)), 0.0_dp, 1.0e-14_dp, 'ed3 start vectors lie in the symmetric subspace')

  ! ---- (C) chains against the dense solve on the unreduced operator, (B) block against single
  allocate (km(ndu, ndu), ev(ndu), col(ndu), rhs(ndu, 2), ipiv(ndu), ep(9, 2), ur(ndu), hr(ndu))
  e_on = 0.0_dp; e_off = 0.0_dp; eb = 0.0_dp
  do iz = 1, 2
    z = cmplx(0.15_dp, merge(0.03_dp, 0.003_dp, iz == 1), dp); zs = [z, z + (0.01_dp, 0.0_dp)]
    ! edc: z - H(z)
    call edc_setup(mdl, pp, alpha, 1, en)
    call augcls_set_params(op, mdl); call augcls_freeze_site(op, 1, 1)
    call augcls_set_params(opu, mdl); call augcls_freeze_site(opu, 1, 1)
    do l = 1, 9; do i = 1, 2; ep(l, i) = z - en%dc(l) * (z - en%enu(l, i)) / en%del(l, i); end do; end do
    do i = 1, ndu
      ev = (0.0_dp, 0.0_dp); ev(i) = (1.0_dp, 0.0_dp)
      ur = 0.0_dp; ur(i) = 1.0_dp
      call opu%apply(ur, hr)
      call augcls_rand_c(opu, ep, ev, col)
      km(:, i) = -cmplx(hr, 0.0_dp, dp) - col; km(i, i) = km(i, i) + z
    end do
    rhs = (0.0_dp, 0.0_dp); rhs(iv, 1) = (1.0_dp, 0.0_dp)
    call zgesv(ndu, 1, km, ndu, ipiv, rhs, ndu, info); gex = rhs(iv, 1)
    ed_project = .true.;  call edc_green(op, wn, en, z, iv, nlev, g1);  e_on(1) = max(e_on(1), abs(g1 - gex) / abs(gex))
    call edc_green_blk(op, bw, en, zs, iv, nlev, gb);                    eb(1) = max(eb(1), abs(gb(1) - g1) / abs(g1))
    ed_project = .false.; call edc_green(op, wn, en, z, iv, nlev, g1);  e_off(1) = max(e_off(1), abs(g1 - gex) / abs(gex))
    ! ed3: K'(z)
    call ed3_setup(mdl, pp, alpha, 1, e3)
    call augcls_set_params(op, mdl); call augcls_freeze_site(op, 1, 1)
    call augcls_set_params(opu, mdl); call augcls_freeze_site(opu, 1, 1)
    call ed3_start(opu, w3u, e3, iv)
    call ed3_coefs(e3, z, et, uu, qt)
    do i = 1, ndu
      ev = (0.0_dp, 0.0_dp); ev(i) = (1.0_dp, 0.0_dp)
      call ed3_kmul(opu, w3u, et, uu, qt, ev, col); km(:, i) = col
    end do
    if (iz == 1) call check_close(maxval(abs(km - transpose(km))) / maxval(abs(km)), 0.0_dp, 1.0e-14_dp, 'unreduced K'' is complex symmetric')
    rhs(:, 1) = w3u%av; rhs(:, 2) = w3u%vv
    call zgesv(ndu, 2, km, ndu, ipiv, rhs, ndu, info)
    q0ex = sum(w3u%av * rhs(:, 1)); q2ex = sum(w3u%vv * rhs(:, 2)); oc = e3%o(iv, 1)
    ed_project = .true.; call ed3_start(op, w3, e3, iv); call ed3_blocks(op, w3, e3, z, iv, nlev, bk)
    e_on(2) = max(e_on(2), abs(bk(0) + 2.0_dp * oc * bk(1) + oc**2 * bk(2) - q0ex) / abs(q0ex), abs(bk(2) - q2ex) / abs(q2ex))
    call ed3_blocks_blk(op, w3, bw, e3, zs, iv, nlev, bkb)
    eb(2) = max(eb(2), maxval(abs(bkb(0:2, 1) - bk(0:2))) / maxval(abs(bk(0:2))))
    ed_project = .false.; call ed3_start(op, w3, e3, iv); call ed3_blocks(op, w3, e3, z, iv, nlev, bk)
    e_off(2) = max(e_off(2), abs(bk(0) + 2.0_dp * oc * bk(1) + oc**2 * bk(2) - q0ex) / abs(q0ex), abs(bk(2) - q2ex) / abs(q2ex))
    ed_project = .true.
  end do
  print '(a,i0,a,2es10.2,a,2es10.2)', '  ', nlev, ' steps, relative error against the dense solve (edc, ed3): projected ', e_on, ', not projected ', e_off
  print '(a,2es10.2)', '  block against single energy (edc, ed3): ', eb
  call check_close(e_on(1), 0.0_dp, 1.0e-11_dp, 'edc chain on the reduced operator == dense solve on the unreduced operator')
  call check_close(e_on(2), 0.0_dp, 1.0e-11_dp, 'ed3 chains (q0, q2) on the reduced operator == dense solve on the unreduced operator')
  call check(e_off(2) > 1.0e4_dp * e_on(2), 'without the projection the ed3 chain stalls far above the projected one')
  call check_close(eb(1), 0.0_dp, 1.0e-11_dp, 'edc_green_blk == edc_green at 144 steps (ordered cell)')
  call check_close(eb(2), 0.0_dp, 1.0e-11_dp, 'ed3_blocks_blk == ed3_blocks at 144 steps (ordered cell)')
  call report()

contains

  !> deterministic pseudo-random numbers in (-0.5, 0.5)
  subroutine lcg(v)
    real(dp), intent(out) :: v(:)
    integer(8), save :: s = 20261110_8
    integer :: k
    do k = 1, size(v)
      s = mod(s * 6364136223846793005_8 + 1442695040888963407_8, 9223372036854775807_8)
      v(k) = real(mod(abs(s), 1000003_8), dp) / 1000003.0_dp - 0.5_dp
    end do
  end subroutine lcg
end program test_augcls_project
