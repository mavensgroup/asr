!> Effective pair interactions (asr_epi) against exact references on a 16-site bcc cell (dense matrices, no recursion).
!>  (G) pair types, pair groups and the symmetry-adapted start vectors;
!>  (D) disorder at one energy: the frozen-pair chains with the trace formula of asr_epi equal the exact average of
!>      Tr[G dH/dlambda] over all configurations of the random sites (ternary; dH by finite differences of the dense unscaled
!>      H^(2)).  This checks the trace formula, the scaling, the conditional average and the group reduction;
!>  (B) energy blocks: epi_blocks_blk == epi_blocks on the same operator for a block of energies, a short block padded with
!>      idle members, and chains that close inside the block;
!>  (T) thermodynamic integration: U(c, s) from the Fermi contour and Gauss-Legendre in lambda equals the mixed second
!>      difference of the exact grand potentials (ordered medium);
!>  (P) path points: the error of 1 to 4 Gauss-Legendre nodes for species that are close to each other (the case of neighbours
!>      in a transition series) and for the far species of (T).  Result: 2 nodes are enough for close species (0.01 %); 1 node
!>      is not safe (2 % of the largest entry, a factor 2 on a small one); far species need 4 (0.1 %);
!>  (S) the shell values of the screening kernel.
program test_epi
  use asr_kinds, only: dp, pi
  use asr_augspace, only: graph_t
  use asr_augcls
  use asr_augblk, only: blk_work_t, blk_work_init
  use asr_augorb, only: orb_alloy_t
  use asr_edn, only: edn_t
  use asr_epi
  use asr_endep, only: ed_depth_tol, ed_tail, ed_project
  use asr_lmto_alloy, only: build_lmto_alloy
  use asr_lattices, only: periodic_graph
  use asr_symmetry
  use asr_lmto_ham, only: potpar_t
  use asr_contour
  use asr_screening
  use asr_testing
  use asr_lapack, only: zgesv, dsyev
  !$ use omp_lib, only: omp_set_num_threads
  implicit none
  integer, parameter :: lsup = 2, mg = 2 * lsup
  type(potpar_t) :: pp(3)
  type(epi_par_t) :: par(3)
  real(dp) :: alpha(0:2), d1(9, 9, 1)
  type(graph_t) :: g, gf
  type(orb_alloy_t) :: m0
  integer, allocatable :: coord(:, :), coordf(:, :), perm1(:, :), cls(:), perm(:, :), mats(:, :, :)
  real(dp), allocatable :: d(:, :, :)
  type(epi_pair_t), allocatable :: pairs(:)
  character(len=:), allocatable :: msg
  integer :: l, i, jp

  print '(a)', 'test_epi'
  !$ call omp_set_num_threads(1)                     ! 16-site operators: threads only add overhead
  do l = 0, 2
    pp(1)%enu(l) = 0.1_dp * l; pp(1)%c(l) = 0.3_dp - 0.15_dp * l; pp(1)%sqrdel(l) = 0.30_dp - 0.07_dp * l
    pp(1)%p(l) = 0.9_dp - 0.1_dp * l; pp(1)%gamma(l) = 0.25_dp - 0.1_dp * l
    pp(2)%enu(l) = 0.05_dp * l; pp(2)%c(l) = 0.45_dp - 0.1_dp * l; pp(2)%sqrdel(l) = 0.27_dp - 0.06_dp * l
    pp(2)%p(l) = 0.8_dp - 0.1_dp * l; pp(2)%gamma(l) = 0.2_dp - 0.08_dp * l
    pp(3)%enu(l) = 0.02_dp + 0.08_dp * l; pp(3)%c(l) = 0.38_dp - 0.18_dp * l; pp(3)%sqrdel(l) = 0.33_dp - 0.08_dp * l
    pp(3)%p(l) = 0.85_dp - 0.1_dp * l; pp(3)%gamma(l) = 0.22_dp - 0.09_dp * l
  end do
  do i = 1, 3; pp(i)%alpha_rep = pp(i)%gamma; end do
  alpha = [0.3485_dp, 0.05303_dp, 0.010714_dp]
  do i = 1, 3; par(i) = epi_species(pp(i), alpha); end do
  call build_lmto_alloy('bcc', lsup, 1.2_dp, 0.7_dp, pp(1), pp(2), 0.5_dp, g, coord, m0, alpha_in=alpha)
  allocate (perm1(g%n, 1)); perm1(:, 1) = [(i, i=1, g%n)]
  d1 = 0.0_dp; do l = 1, 9; d1(l, l, 1) = 1.0_dp; end do
  ed_depth_tol = 0.0_dp; ed_tail = 'none'; ed_project = .true.
  jp = site_at([1, 1, 1])

  call part_groups()
  call part_disorder()
  call part_ti()
  call part_path()
  call part_screen()
  call report()

contains

  integer function site_at(c)
    integer, intent(in) :: c(3)
    integer :: k
    site_at = 0
    do k = 1, size(coord, 2)
      if (all(modulo(coord(:, k) - c, mg) == 0)) site_at = k
    end do
    if (site_at == 0) error stop 'site_at: not a site'
  end function site_at

  ! ---------------------------------------------------------------- (G)
  subroutine part_groups()
    real(dp) :: vec(9, 9), chi(16, 9), e9(9, 9)
    integer :: lo(9), a, gg, k, neven
    integer, allocatable :: c3(:, :), cl3(:)
    type(graph_t) :: g3
    allocate (cls(g%n)); cls = 1
    call epi_pair_list(coord, mg, cls, 1, pairs, msg)
    call check(len(msg) == 0 .and. size(pairs) == 1, 'bcc: one pair type in the first shell')
    call check(pairs(1)%mult == 8 .and. pairs(1)%r2 == 3 .and. pairs(1)%nop == 2, 'bcc first shell: 8 sites, pair group of order 2')
    call epi_pair_list(coord, mg, cls, 2, pairs, msg)
    call check(len(msg) > 0, 'a shell that does not fit into the supercell is refused')
    cls(jp) = 2
    call epi_pair_ops(coord, mg, cls, perm, mats)
    call orbital_rep(mats, d)
    call epi_adapted(d, vec, chi, lo)
    e9 = matmul(transpose(vec), vec)
    do a = 1, 9; e9(a, a) = e9(a, a) - 1.0_dp; end do
    call check_close(maxval(abs(e9)), 0.0_dp, 1.0e-13_dp, 'adapted start vectors are orthonormal')
    neven = 0
    do a = 1, 9
      if (all(chi(1:size(perm, 2), a) > 0.0_dp)) neven = neven + 1
      do gg = 1, size(perm, 2)
        e9(:, a) = matmul(d(:, :, gg), vec(:, a)) - chi(gg, a) * vec(:, a)
      end do
    end do
    call check_close(maxval(abs(e9)), 0.0_dp, 1.0e-12_dp, 'adapted start vectors are common eigenvectors of the pair group')
    call check(neven == 6, 'bcc (111) pair: 6 even and 3 odd start vectors under the mirror x <-> y')
    call check(all(lo == [1, 2, 2, 2, 5, 5, 5, 5, 5]), 'each start vector lies in one l')
    deallocate (cls)
    ! fcc first shell, and three bcc shells on a larger cell (plain and B2 decorated)
    call periodic_graph('fcc', 2, gf, coordf)
    allocate (cls(gf%n)); cls = 1
    call epi_pair_list(coordf, 4, cls, 1, pairs, msg)
    call check(size(pairs) == 1 .and. pairs(1)%mult == 12 .and. pairs(1)%nop == 4, 'fcc first shell: 12 sites, pair group of order 4')
    deallocate (cls)
    call periodic_graph('bcc', 3, g3, c3)
    allocate (cl3(g3%n)); cl3 = 1
    call epi_pair_list(c3, 6, cl3, 3, pairs, msg)
    call check(size(pairs) == 3, 'bcc: three shells, one pair type each')
    if (size(pairs) == 3) then
      call check(all(pairs%mult == [8, 6, 12]) .and. all(pairs%r2 == [3, 4, 8]), 'bcc shells: 8, 6, 12 sites')
      call check(all(pairs%nop == [2, 4, 4]), 'bcc pair groups: 2 (111), 4 (001, abelian part of C4v), 4 (110)')
    end if
    do k = 1, g3%n; cl3(k) = 1 + modulo(c3(1, k), 2); end do
    call epi_pair_list(c3, 6, cl3, 2, pairs, msg)
    call check(size(pairs) == 2, 'B2: two pair types in two shells')
    if (size(pairs) == 2) call check(all(pairs%cls2 == [2, 1]), 'B2: first shell on the other sublattice, second on the own')
  end subroutine part_groups

  ! ---------------------------------------------------------------- dense unscaled H^(2)
  !> H of the 16-site cell with species spsite(site) (1..3) on every site, site 1 on the path point pt and site jp on species sj
  subroutine dense_h(spsite, pt, sj, h)
    integer, intent(in) :: spsite(:), sj
    type(epi_path_t), intent(in) :: pt
    real(dp), allocatable, intent(out) :: h(:, :)
    type(class_alloy_t) :: md
    type(augcls_op) :: o1
    logical :: fl(g%n)
    real(dp), allocatable :: x(:), y(:)
    integer :: s, k, n
    md%nl = 9; md%nsp = 4; md%ncls = 3
    allocate (md%nsc(3), md%spc(1, 3), md%conc(1, 3)); md%nsc = 1; md%spc(1, :) = [1, 2, 3]; md%conc = 1.0_dp
    md%block = m0%block; md%block0 = m0%block0
    allocate (md%c(9, 4), md%d(9, 4), md%o(9, 4), md%e(9, 4))
    do s = 1, 3
      md%c(:, s) = par(s)%c; md%d(:, s) = par(s)%d; md%o(:, s) = par(s)%o; md%e(:, s) = par(s)%e
    end do
    md%c(:, 4) = pt%c; md%d(:, 4) = pt%d; md%o(:, 4) = pt%o; md%e(:, 4) = pt%e
    fl = .false.
    call build_augcls_operator(g, md, spsite, perm1, d1, 1, 6, o1, maxrank=1, flippable=fl)
    if (o1%nn /= g%n) error stop 'dense_h: node set /= sites'
    call augcls_freeze_any(o1, 1, 4); call augcls_freeze_any(o1, jp, sj)
    n = 9 * o1%nn
    allocate (h(n, n), x(n), y(n))
    do k = 1, n
      x = 0.0_dp; x(k) = 1.0_dp
      call o1%apply(x, y); h(:, k) = y
    end do
    deallocate (o1%wk)
  end subroutine dense_h

  !> Tr[(z - H)^-1 dH/dlambda] + 2 sum_L kappa_L at lambda on the path a -> b (the form of epi_integrand summed over site 0)
  function ref_trace(spsite, a, b, lam, sj, z) result(tr)
    integer, intent(in) :: spsite(:), a, b, sj
    real(dp), intent(in) :: lam
    complex(dp), intent(in) :: z
    complex(dp) :: tr
    real(dp), parameter :: hs = 1.0e-2_dp
    real(dp), allocatable :: h0(:, :), hp(:, :), hm(:, :), hp2(:, :), hm2(:, :)
    complex(dp), allocatable :: km(:, :), rhs(:, :)
    integer, allocatable :: ipiv(:)
    type(epi_path_t) :: pt
    integer :: n, k, info
    call dense_h(spsite, epi_path(par(a), par(b), lam), sj, h0)
    call dense_h(spsite, epi_path(par(a), par(b), lam + hs), sj, hp)
    call dense_h(spsite, epi_path(par(a), par(b), lam - hs), sj, hm)
    call dense_h(spsite, epi_path(par(a), par(b), lam + 2.0_dp * hs), sj, hp2)
    call dense_h(spsite, epi_path(par(a), par(b), lam - 2.0_dp * hs), sj, hm2)
    n = size(h0, 1)
    allocate (km(n, n), rhs(n, n), ipiv(n))
    km = -cmplx(h0, 0.0_dp, dp)
    do k = 1, n; km(k, k) = km(k, k) + z; end do
    rhs = cmplx((-hp2 + 8.0_dp * hp - 8.0_dp * hm + hm2) / (12.0_dp * hs), 0.0_dp, dp)
    call zgesv(n, n, km, n, ipiv, rhs, n, info)
    if (info /= 0) error stop 'ref_trace: zgesv failed'
    tr = (0.0_dp, 0.0_dp)
    do k = 1, n; tr = tr + rhs(k, k); end do
    pt = epi_path(par(a), par(b), lam)
    tr = tr + 2.0_dp * sum(pt%dk)
  end function ref_trace

  !> sum over the nine start vectors of the trace formula from the frozen-pair chains: site 1 on the path a -> b at lam, site jp
  !> on species sj; op is built for mdl (4 species)
  function chain_trace(op, mdl, dd, a, b, lam, sj, z, nlev, tailmax) result(tr)
    type(augcls_op), intent(inout) :: op
    type(class_alloy_t), intent(in) :: mdl
    real(dp), intent(in) :: dd(:, :, :), lam
    integer, intent(in) :: a, b, sj, nlev
    complex(dp), intent(in) :: z
    real(dp), intent(out), optional :: tailmax
    complex(dp) :: tr, x(0:3)
    type(class_alloy_t) :: ml
    type(edn_t) :: ed
    type(epi_work_t) :: wk
    type(epi_path_t) :: pt
    real(dp) :: vec(9, 9), chi(16, 9)
    integer :: lo(9), ia
    ml = mdl
    pt = epi_path(par(a), par(b), lam)
    call epi_setup(ml, par, pt, ed)
    call augcls_set_params(op, ml)
    call augcls_freeze_any(op, 1, 4); call augcls_freeze_any(op, jp, sj)
    call epi_adapted(dd, vec, chi, lo)
    call epi_work_init(wk, op%dim())
    tr = (0.0_dp, 0.0_dp)
    if (present(tailmax)) tailmax = 0.0_dp
    do ia = 1, 9
      call augcls_set_chi(op, chi(1:op%nop, ia))
      call epi_start(op, wk, vec(:, ia))
      call epi_blocks(op, wk, ed, z, nlev, x)
      tr = tr + epi_integrand(pt, lo(ia), z, x)
      if (present(tailmax)) tailmax = max(tailmax, abs(x(3)) / abs(x(1)))
    end do
  end function chain_trace

  ! ---------------------------------------------------------------- (D)
  subroutine part_disorder()
    type(class_alloy_t) :: mdl
    type(augcls_op) :: op, opu
    logical :: fl(g%n)
    integer :: spsite(g%n), r1, r2, s1, s2, sj
    real(dp), parameter :: xc(3) = [0.5_dp, 0.3_dp, 0.2_dp], lam = 0.35_dp
    complex(dp) :: z, ref, tr, tru
    real(dp) :: tmx
    r1 = site_at([2, 0, 0]); r2 = site_at([0, 2, 0])
    allocate (cls(g%n)); cls = 1; cls(r1) = 2; cls(r2) = 2; cls(jp) = 3
    fl = .false.; fl(r1) = .true.; fl(r2) = .true.
    mdl%nl = 9; mdl%nsp = 4; mdl%ncls = 3
    allocate (mdl%nsc(3), mdl%spc(3, 3), mdl%conc(3, 3)); mdl%spc = 0; mdl%conc = 0.0_dp
    mdl%nsc = [1, 3, 1]; mdl%spc(1, 1) = 1; mdl%spc(:, 2) = [1, 2, 3]; mdl%spc(1, 3) = 1
    mdl%conc(1, 1) = 1.0_dp; mdl%conc(:, 2) = xc; mdl%conc(1, 3) = 1.0_dp
    mdl%block = m0%block; mdl%block0 = m0%block0
    call epi_pair_ops(coord, mg, cls, perm, mats)
    call orbital_rep(mats, d)
    call check(size(perm, 2) == 2, 'pair group with two random sites related by the mirror: order 2')
    call build_augcls_operator(g, mdl, cls, perm, d, 1, 8, op, maxrank=2, flippable=fl)
    call build_augcls_operator(g, mdl, cls, perm1, d1, 1, 8, opu, maxrank=2, flippable=fl)
    call check(opu%nn == 9 * g%n, 'unreduced augmented space is complete: 16 sites x 9 configurations')
    print '(a,i0,a,i0)', '  nodes: reduced ', op%nn, ', unreduced ', opu%nn
    z = cmplx(0.12_dp, 0.3_dp, dp)
    do sj = 2, 3
      ref = (0.0_dp, 0.0_dp)
      do s1 = 1, 3
        do s2 = 1, 3
          spsite = 1; spsite(r1) = s1; spsite(r2) = s2
          ref = ref + xc(s1) * xc(s2) * ref_trace(spsite, 1, 2, lam, sj, z)
        end do
      end do
      tr = chain_trace(op, mdl, d, 1, 2, lam, sj, z, 140, tmx)
      tru = chain_trace(opu, mdl, d1, 1, 2, lam, sj, z, 140)
      print '(a,i0,a,2f14.9,a,2es9.2,a,es9.2)', '  partner species ', sj, ': exact average ', ref, '  chains - exact (reduced, unreduced) ', &
        abs(tr - ref), abs(tru - ref), '  X1 mismatch of the two chains ', tmx
      call check_close(abs(tr - ref) / abs(ref), 0.0_dp, 1.0e-8_dp, 'frozen-pair chains (pair group) == exact configuration average of Tr[G dH/dlambda]')
      call check_close(abs(tru - ref) / abs(ref), 0.0_dp, 1.0e-8_dp, 'frozen-pair chains (no symmetry) == exact configuration average')
    end do
    call part_block(op, mdl, d)
    deallocate (cls, op%wk, opu%wk)
  end subroutine part_disorder

  ! ---------------------------------------------------------------- (B)
  !> Block chains against single chains on the operator of (D): the four blocks X0, X1, X2 and the X1 mismatch for every start
  !> vector at five energies.  12 steps: the chains are open and the two agree to rounding (1e-15 to 1e-12).  A work space of
  !> 8 members holds the block of 5 (three idle members).  400 steps: every chain closes, so every member leaves the block
  !> and takes the fraction without a tail (1e-13).  Between the two depths the 16-site cell is close to exhaustion: the
  !> chains then amplify the rounding-level kernel difference (1e-7 at 30 steps), as on the ring of test_augcls_blk.
  subroutine part_block(op, mdl, dd)
    type(augcls_op), intent(inout) :: op
    type(class_alloy_t), intent(in) :: mdl
    real(dp), intent(in) :: dd(:, :, :)
    integer, parameter :: nz = 5
    type(class_alloy_t) :: ml
    type(edn_t) :: ed
    type(epi_work_t) :: wk
    type(blk_work_t) :: bw5, bw8
    type(epi_path_t) :: pt
    complex(dp) :: zs(nz), x(0:3), xb(0:3, nz), xp(0:3, nz)
    real(dp) :: vec(9, 9), chi(16, 9), e5, e8, ec, sc
    integer :: lo(9), ia, j
    do j = 1, nz
      zs(j) = cmplx(-0.3_dp + 0.17_dp * j, 0.6_dp / real(j, dp), dp)
    end do
    ml = mdl
    pt = epi_path(par(1), par(2), 0.35_dp)
    call epi_setup(ml, par, pt, ed)
    call augcls_set_params(op, ml)
    call augcls_freeze_any(op, 1, 4); call augcls_freeze_any(op, jp, 3)
    call epi_adapted(dd, vec, chi, lo)
    call epi_work_init(wk, op%dim())
    call blk_work_init(bw5, nz, op%dim()); call blk_work_init(bw8, 8, op%dim())
    e5 = 0.0_dp; e8 = 0.0_dp; ec = 0.0_dp
    do ia = 1, 9
      call augcls_set_chi(op, chi(1:op%nop, ia))
      call epi_start(op, wk, vec(:, ia))
      call epi_blocks_blk(op, wk, bw5, ed, zs, 12, xb)
      call epi_blocks_blk(op, wk, bw8, ed, zs, 12, xp)
      do j = 1, nz
        call epi_blocks(op, wk, ed, zs(j), 12, x)
        sc = maxval(abs(x(0:2)))
        e5 = max(e5, maxval(abs(xb(:, j) - x)) / sc); e8 = max(e8, maxval(abs(xp(:, j) - x)) / sc)
      end do
      call epi_blocks_blk(op, wk, bw8, ed, zs(1:3), 400, xp(:, 1:3))
      do j = 1, 3
        call epi_blocks(op, wk, ed, zs(j), 400, x)
        ec = max(ec, maxval(abs(xp(0:2, j) - x(0:2))) / maxval(abs(x(0:2))))
      end do
    end do
    call check_close(e5, 0.0_dp, 1.0e-11_dp, 'energy block: epi_blocks_blk == epi_blocks (X0, X1, X2, X1 mismatch), 5 energies, 12 steps')
    call check_close(e8, 0.0_dp, 1.0e-11_dp, 'energy block: the same with a short block (3 idle members)')
    call check_close(ec, 0.0_dp, 1.0e-10_dp, 'energy block: closed chains (400 steps) == single chains')
  end subroutine part_block

  ! ---------------------------------------------------------------- (T)
  !> exact dOmega/dlambda = sum_i f(e_i) <i|dH/dlambda|i> of the dense H on the path a -> b at lam, partner species sj
  function dense_d(spsite, a, b, lam, sj, mu, kt, emin) result(dd)
    integer, intent(in) :: spsite(:), a, b, sj
    real(dp), intent(in) :: lam, mu, kt
    real(dp), intent(out), optional :: emin
    real(dp) :: dd
    real(dp), parameter :: hs = 1.0e-2_dp
    real(dp), allocatable :: h0(:, :), hp(:, :), hm(:, :), hp2(:, :), hm2(:, :), ev(:), work(:)
    integer :: n, k, info
    call dense_h(spsite, epi_path(par(a), par(b), lam), sj, h0)
    call dense_h(spsite, epi_path(par(a), par(b), lam + hs), sj, hp)
    call dense_h(spsite, epi_path(par(a), par(b), lam - hs), sj, hm)
    call dense_h(spsite, epi_path(par(a), par(b), lam + 2.0_dp * hs), sj, hp2)
    call dense_h(spsite, epi_path(par(a), par(b), lam - 2.0_dp * hs), sj, hm2)
    n = size(h0, 1)
    hp = (-hp2 + 8.0_dp * hp - 8.0_dp * hm + hm2) / (12.0_dp * hs)
    allocate (ev(n), work(8 * n))
    call dsyev('V', 'U', n, h0, n, ev, work, size(work), info)
    if (info /= 0) error stop 'dense_d: dsyev failed'
    dd = 0.0_dp
    do k = 1, n
      dd = dd + dot_product(h0(:, k), matmul(hp, h0(:, k))) / (1.0_dp + exp(min(700.0_dp, (ev(k) - mu) / kt)))
    end do
    if (present(emin)) emin = ev(1)
  end function dense_d

  subroutine part_ti()
    integer, parameter :: npl = 2, nmain = 32, nstep = 16
    real(dp), parameter :: kt = 0.04_dp                ! above the level spacing of the 144 states: smooth in lambda
    type(class_alloy_t) :: mdl
    type(augcls_op) :: op
    logical :: fl(g%n)
    integer :: spsite(g%n), c, s, nq, j, il, nlam, n, info
    real(dp) :: om(3, 3), uex(3, 3), uti(3, 3), mu, ebot, xl(8), wl(8), dl(3), acc, dref, err(5), asy(5)
    real(dp), allocatable :: h(:, :), ev(:), work(:)
    complex(dp), allocatable :: zq(:), wq(:)
    type(epi_path_t) :: pt
    allocate (cls(g%n)); cls = 1; cls(jp) = 2
    fl = .false.; spsite = 1
    mdl%nl = 9; mdl%nsp = 4; mdl%ncls = 2
    allocate (mdl%nsc(2), mdl%spc(1, 2), mdl%conc(1, 2)); mdl%nsc = 1; mdl%spc = 1; mdl%conc = 1.0_dp
    mdl%block = m0%block; mdl%block0 = m0%block0
    call epi_pair_ops(coord, mg, cls, perm, mats)
    call orbital_rep(mats, d)
    call build_augcls_operator(g, mdl, cls, perm, d, 1, 6, op, maxrank=1, flippable=fl)
    ! exact grand potentials Omega(c, s) = -kT sum_i ln(1 + exp(-(e_i - mu)/kT)); mu inside the spectrum of the pure cell
    mu = 0.0_dp; ebot = 0.0_dp
    do c = 1, 3
      do s = 1, 3
        pt = epi_path(par(c), par(c), 0.0_dp)
        call dense_h(spsite, pt, s, h)
        n = size(h, 1)
        if (.not. allocated(ev)) allocate (ev(n), work(8 * n))
        call dsyev('N', 'U', n, h, n, ev, work, size(work), info)
        if (c == 1 .and. s == 1) then
          mu = 0.5_dp * (ev(n / 3) + ev(n / 3 + 1)) + 0.013_dp; ebot = ev(1) - 0.4_dp
        end if
        om(c, s) = 0.0_dp
        do j = 1, n
          if (ev(j) < mu) then
            om(c, s) = om(c, s) + (ev(j) - mu) - kt * log(1.0_dp + exp((ev(j) - mu) / kt))
          else
            om(c, s) = om(c, s) - kt * log(1.0_dp + exp(-(ev(j) - mu) / kt))
          end if
        end do
        deallocate (h)
      end do
    end do
    do c = 1, 3; do s = 1, 3
      uex(c, s) = om(c, s) - om(1, s) - om(c, 1) + om(1, 1)
    end do; end do
    ! lambda integration of the exact derivative
    do il = 1, 5
      nlam = merge(8, il, il == 5)
      call gauss_legendre(nlam, xl(1:nlam), wl(1:nlam))
      xl(1:nlam) = 0.5_dp * (1.0_dp + xl(1:nlam)); wl(1:nlam) = 0.5_dp * wl(1:nlam)
      uti = 0.0_dp
      do c = 2, 3
        do j = 1, nlam
          do s = 1, 3; dl(s) = dense_d(spsite, 1, c, xl(j), s, mu, kt); end do
          do s = 2, 3; uti(c, s) = uti(c, s) + wl(j) * (dl(s) - dl(1)); end do
        end do
      end do
      err(il) = maxval(abs(uti - uex)); asy(il) = abs(uti(2, 3) - uti(3, 2))
      print '(a,i0,a,4f11.6,a,es9.2,a,es9.2)', '  ', nlam, ' lambda nodes: U(2,2) U(2,3) U(3,2) U(3,3) (mRy) ', 1.0e3_dp * uti(2:3, 2:3), &
        '   max error (Ry) ', err(il), '  |U(2,3) - U(3,2)| ', asy(il)
    end do
    call check_close(err(4), 0.0_dp, 4.0e-5_dp, 'U(c, s) from 4 lambda nodes == exact grand-potential differences (0.1 %)')
    call check_close(err(5), 0.0_dp, 2.0e-7_dp, 'U(c, s) from 8 lambda nodes == exact grand-potential differences')
    call check(err(4) < err(3) .and. err(3) < err(2), 'the lambda quadrature converges')
    call check(all(asy <= 2.0_dp * err + 1.0e-12_dp) .and. asy(2) > 0.2_dp * err(2), &
               'the asymmetry U(c, s) - U(s, c) measures the quadrature error')
    print '(a,4f11.6)', '  exact:          U(2,2) U(2,3) U(3,2) U(3,3) (mRy) ', 1.0e3_dp * uex(2:3, 2:3)
    call check(maxval(abs(uex(2:3, 2:3))) > 1.0e-4_dp, 'the test pair interactions are not trivially small')
    call check_close(abs(uti(2, 3) - uti(3, 2)), 0.0_dp, 2.0e-7_dp, 'U(c, s) = U(s, c): the path site and the frozen site are interchangeable')
    print '(a,3f11.6)', '  exchange interactions V(1,2) V(1,3) V(2,3) (mRy) ', 1.0e3_dp * [uti(2, 2), uti(3, 3), uti(2, 2) + uti(3, 3) - uti(2, 3) - uti(3, 2)]
    uex = epi_exchange(uti)
    call check_close(uex(1, 2), uti(2, 2), 1.0e-14_dp, 'epi_exchange: V(1, q) = U(q, q) for the reference species 1')
    call check(epi_skip(.true., .true., 3, 2) .and. .not. epi_skip(.true., .true., 3, 1) .and. .not. epi_skip(.true., .true., 2, 3) &
               .and. .not. epi_skip(.true., .true., 3, 3) .and. .not. epi_skip(.false., .true., 3, 2) .and. .not. epi_skip(.true., .false., 3, 2), &
               'epi_skip: only the entries 1 < s < c of one sublattice, only with the option')
    uex = uti; uex(3, 2) = 99.0_dp
    call epi_mirror(uex)
    call check(uex(3, 2) == uti(2, 3) .and. all(uex(2, :) == uti(2, :)) .and. uex(3, 3) == uti(3, 3) .and. all(uex(:, 1) == uti(:, 1)), &
               'epi_mirror: U(3, 2) = U(2, 3), nothing else changes')
    ! Fermi contour of the chains at one path point against the exact derivative
    nq = fermi_nodes_count(npl, nmain, nstep)
    allocate (zq(nq), wq(nq))
    call contour_nodes_fermi(ebot, mu, kt, npl, nmain, nstep, zq, wq)
    call check_close(aimag(sum(wq)), 0.0_dp, 1.0e-9_dp, 'Im sum of the contour weights vanishes (the constant -2 kappa drops out)')
    dref = dense_d(spsite, 1, 3, 0.6_dp, 2, mu, kt, ebot)
    ebot = ebot - 0.4_dp                               ! below the lowest state of this path point
    call contour_nodes_fermi(ebot, mu, kt, npl, nmain, nstep, zq, wq)
    acc = 0.0_dp
    do j = 1, nq
      acc = acc - aimag(wq(j) * chain_trace(op, mdl, d, 1, 3, 0.6_dp, 2, zq(j), 100)) / pi
    end do
    print '(a,f8.4,a,i0,a,2f14.9)', '  mu = ', mu, '  contour nodes ', nq, ':  dOmega/dlambda from the chains, exact ', acc, dref
    call check_close(acc, dref, 1.0e-6_dp, 'dOmega/dlambda: Fermi contour of the chains == sum_i f(e_i) <i|dH/dlambda|i>')
    deallocate (cls, op%wk)
  end subroutine part_ti

  ! ---------------------------------------------------------------- (P)
  !> Path-point convergence for close species.  Species 2 and 3 are species 1 with band centres shifted by +0.03 and -0.04 Ry
  !> and band widths changed by +10 and -10 %.  The midpoint rule (one node) has an error of third order in the species
  !> difference, but the table entries are themselves small differences.  Exact reference: grand-potential differences of
  !> the dense matrices, as in (T).
  subroutine part_path()
    real(dp), parameter :: kt = 0.04_dp
    type(epi_par_t) :: keep(3)
    type(potpar_t) :: q
    integer :: spsite(g%n), c, s, j, il, n, info
    real(dp) :: om(3, 3), uex(3, 3), uti(3, 3), mu, xl(4), wl(4), dl(3), rel(4), sc
    real(dp), allocatable :: h(:, :), ev(:), work(:)
    keep = par
    q = pp(1); q%c = q%c + 0.03_dp; q%enu = q%enu + 0.01_dp; q%sqrdel = 1.05_dp * q%sqrdel; q%gamma = 1.04_dp * q%gamma
    q%alpha_rep = q%gamma; par(2) = epi_species(q, alpha)
    q = pp(1); q%c = q%c - 0.04_dp; q%enu = q%enu - 0.01_dp; q%sqrdel = 0.95_dp * q%sqrdel; q%gamma = 0.97_dp * q%gamma
    q%alpha_rep = q%gamma; par(3) = epi_species(q, alpha)
    spsite = 1; mu = 0.0_dp
    do c = 1, 3
      do s = 1, 3
        call dense_h(spsite, epi_path(par(c), par(c), 0.0_dp), s, h)
        n = size(h, 1)
        if (.not. allocated(ev)) allocate (ev(n), work(8 * n))
        call dsyev('N', 'U', n, h, n, ev, work, size(work), info)
        if (c == 1 .and. s == 1) mu = 0.5_dp * (ev(n / 3) + ev(n / 3 + 1)) + 0.013_dp
        om(c, s) = 0.0_dp
        do j = 1, n
          if (ev(j) < mu) then
            om(c, s) = om(c, s) + (ev(j) - mu) - kt * log(1.0_dp + exp((ev(j) - mu) / kt))
          else
            om(c, s) = om(c, s) - kt * log(1.0_dp + exp(-(ev(j) - mu) / kt))
          end if
        end do
        deallocate (h)
      end do
    end do
    do c = 1, 3; do s = 1, 3
      uex(c, s) = om(c, s) - om(1, s) - om(c, 1) + om(1, 1)
    end do; end do
    sc = maxval(abs(uex(2:3, 2:3)))
    print '(a,4f11.6)', '  close species, exact:    U(2,2) U(2,3) U(3,2) U(3,3) (mRy) ', 1.0e3_dp * uex(2:3, 2:3)
    do il = 1, 4
      call gauss_legendre(il, xl(1:il), wl(1:il))
      xl(1:il) = 0.5_dp * (1.0_dp + xl(1:il)); wl(1:il) = 0.5_dp * wl(1:il)
      uti = 0.0_dp
      do c = 2, 3
        do j = 1, il
          do s = 1, 3; dl(s) = dense_d(spsite, 1, c, xl(j), s, mu, kt); end do
          do s = 2, 3; uti(c, s) = uti(c, s) + wl(j) * (dl(s) - dl(1)); end do
        end do
      end do
      rel(il) = maxval(abs(uti - uex)) / sc
      print '(a,i0,a,4f11.6,a,es9.2)', '  close species, ', il, ' node(s): U(2,2) U(2,3) U(3,2) U(3,3) (mRy) ', 1.0e3_dp * uti(2:3, 2:3), &
        '   max error / max |U| ', rel(il)
    end do
    call check(sc > 1.0e-4_dp, 'close species: the pair interactions are not trivially small')
    call check_close(rel(1), 0.0_dp, 5.0e-2_dp, 'close species: the midpoint rule (1 path point) is within 5 % of the largest U')
    call check(abs(uex(2, 3)) < 0.05_dp * sc .and. rel(1) * sc > 0.5_dp * abs(uex(2, 3)), &
               'close species: the midpoint error is as large as a small table entry (1 path point is not safe)')
    call check_close(rel(2), 0.0_dp, 2.0e-4_dp, 'close species: 2 path points are within 0.02 %')
    call check(rel(4) < rel(2) .and. rel(2) < 1.0e-2_dp * rel(1), 'close species: 2 path points are 100 times better than 1; the quadrature converges')
    par = keep
  end subroutine part_path

  ! ---------------------------------------------------------------- (S)
  subroutine part_screen()
    type(screen_t) :: s
    real(dp) :: rv(3, 3), gr(3), wr(3), al
    al = 5.4_dp
    call screen_lattice('bcc', al, 1.0_dp, s, msg)
    rv(:, 1) = 0.0_dp; rv(:, 2) = 0.5_dp * al * [1.0_dp, 1.0_dp, 1.0_dp]; rv(:, 3) = al * [0.0_dp, 0.0_dp, 1.0_dp]
    call screen_shells('bcc', al, s, rv, gr, wr)
    call check_close(gr(1), s%g00, 1.0e-12_dp, 'screen_shells: G(0) = G00')
    call check_close(gr(2), s%g1, 1.0e-12_dp, 'screen_shells: G(R1) = G1')
    call check_close(wr(1), s%w00, 1.0e-12_dp, 'screen_shells: W(0) = W00')
    call check_close(wr(2), s%w1, 1.0e-12_dp, 'screen_shells: W(R1) = W1')
    call check(wr(2) > 0.0_dp .and. abs(wr(3)) < wr(2), 'the screened kernel is positive on the first shell and smaller on the second')
    print '(a,3f10.5)', '  bcc, beta 1: W(0) W(R1) W(R2) (Ry/e^2) ', wr
    call screen_lattice('bcc', al, 0.0_dp, s, msg)
    call screen_shells('bcc', al, s, rv, gr, wr)
    call check_close(wr(3), 2.0_dp / al, 1.0e-14_dp, 'screen_shells: unscreened kernel is 2 / R')
  end subroutine part_screen
end program test_epi
