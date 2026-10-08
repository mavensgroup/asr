!> Energy-dependent augmented-space recursion (asr_endep) vs exact configuration averages.
!> (1) ordered limit A = B: recursion == exact k-space inverse of (P(z) - S(k)) with the same screened blocks;
!> (2) disorder: species-resolved G from the recursion vs the ensemble average of exact dense inverses of
!>     (P(z) - S) on random periodic supercells (sites of each species and orbitals of each class averaged).
program test_endep
  use asr_kinds, only: dp, pi
  use asr_lmto_ham, only: potpar_t
  use asr_augspace, only: graph_t
  use asr_augorb, only: orb_alloy_t
  use asr_augsym
  use asr_endep
  use asr_lmto_alloy
  use asr_symmetry
  use asr_strconst, only: lattice_cluster, screened_cluster, avg_ws_radius, hardcore_alpha
  use asr_lattices, only: periodic_graph_shells
  use asr_testing
  implicit none
  type(potpar_t) :: pa, pb
  type(graph_t) :: g, gs
  type(orb_alloy_t) :: mdl, mdls
  type(augsym_op) :: sop
  type(ed_params_t) :: spa, spb
  type(ed_work_t) :: wk
  integer, allocatable :: coord(:, :), perm(:, :), mats(:, :, :), cs(:, :), bvs(:, :)
  real(dp), allocatable :: d(:, :, :), vec(:, :), blk(:, :, :)
  complex(dp), allocatable :: amat(:, :), ainv(:, :), gk(:, :)
  real(dp) :: alpha(0:2), al4(0:3), alpha9(9), x, r, err1, err2, kk(3), bv(3, 3), av(3, 3), ph, avg(2, 4), nrm(2, 4)
  complex(dp) :: gavg, ga, gb, z, s9(9, 9), p9a(9), p9b(9), m1(9), m2(9), l1(9), l2(9), gex(2, 4), gref
  integer, parameter :: iorb(4) = [1, 3, 5, 7], lsup = 4, nreal = 2, orb_l(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]
  integer :: maxrank, l, ic, ie, nsite, i, j, k, ir, ik, nk, n1, n2, n3, nd, io, spc, info, kc, cls(9), nlev, ia
  logical, allocatable :: isa(:)
  logical :: ordered
  integer, allocatable :: ipiv(:)
  complex(dp), allocatable :: work(:)
  real(dp), parameter :: za(4) = [-0.1_dp, 0.0_dp, 0.1_dp, 0.25_dp]
  character(len=8) :: dummy

  print '(a)', 'test_endep'
  do l = 0, 2
    pa%enu(l) = 0.1_dp * l; pa%c(l) = 0.3_dp - 0.15_dp * l; pa%sqrdel(l) = 0.30_dp - 0.1_dp * l
    pa%p(l) = 0.9_dp - 0.1_dp * l; pa%gamma(l) = 0.25_dp - 0.1_dp * l; pa%alpha_rep(l) = pa%gamma(l)
    pb%enu(l) = 0.05_dp * l; pb%c(l) = 0.45_dp - 0.1_dp * l; pb%sqrdel(l) = 0.28_dp - 0.1_dp * l
    pb%p(l) = 0.8_dp - 0.1_dp * l; pb%gamma(l) = 0.2_dp - 0.08_dp * l; pb%alpha_rep(l) = pb%gamma(l)
  end do
  ordered = .true.                                    ! default (fast unit test): ordered limit; argument 'full' = disorder run (slow, desktop)
  if (command_argument_count() > 0) then; call get_command_argument(1, dummy); ordered = (trim(dummy) /= 'full'); end if
  if (ordered) pb = pa                                ! ordered limit (A = B)
  al4 = hardcore_alpha(0.7_dp); alpha = al4(0:2)
  do i = 1, 9; alpha9(i) = alpha(orb_l(i)); end do
  cls = [1, 2, 2, 2, 3, 3, 3, 3, 3]
  call load_params(spa, pa); call load_params(spb, pb)
  nlev = 8; maxrank = 4
  if (.not. ordered) nlev = 3
  if (command_argument_count() > 2) then; call get_command_argument(3, dummy); read (dummy, *) maxrank; end if
  if (command_argument_count() > 1) then; call get_command_argument(2, dummy); read (dummy, *) nlev; end if

  ! ---------------- (2) disorder vs ensemble of dense inverses
  x = 0.4_dp
  call build_lmto_alloy('fcc', 6, 1.2_dp, 0.7_dp, pa, pb, x, g, coord, mdl)
  call d4h_ops(coord, 12, perm, mats)
  call orbital_rep(mats, d)
  call build_augsym_operator(g, mdl, perm, d, 1, nlev, sop, maxrank)
  print '(a,i0)', '    reduced nodes: ', sop%nn
  call ed_work_init(wk, sop%dim())
  ! supercell
  call build_lmto_alloy('fcc', lsup, 1.2_dp, 0.7_dp, pa, pb, x, gs, cs, mdls)
  nsite = gs%n; nd = 9 * nsite
  allocate (amat(nd, nd), ainv(nd, nd), isa(nsite), ipiv(nd), work(64 * nd))
  err2 = 0.0_dp
  do ie = 1, size(za)
    z = cmplx(za(ie), 0.4_dp, dp)
    call ed_functions(z, spa, alpha9, p9a, m1, l1)
    call ed_functions(z, spb, alpha9, p9b, m2, l2)
    gex = 0.0_dp; nrm = 0.0_dp
    do ir = 1, nreal
      call random_species(isa, x)
      amat = 0.0_dp
      do i = 1, nsite
        do l = 1, 9
          amat(9 * (i - 1) + l, 9 * (i - 1) + l) = merge(p9a(l), p9b(l), isa(i))
        end do
        do l = 1, 9; do j = 1, 9
          amat(9 * (i - 1) + l, 9 * (i - 1) + j) = amat(9 * (i - 1) + l, 9 * (i - 1) + j) - mdls%block0(l, j)
        end do; end do
        do k = 1, size(mdls%block, 3)
          kc = gs%nbr(gs%nbr_ptr(i) + k - 1)
          do l = 1, 9; do j = 1, 9
            amat(9 * (i - 1) + l, 9 * (kc - 1) + j) = amat(9 * (i - 1) + l, 9 * (kc - 1) + j) - mdls%block(l, j, k)
          end do; end do
        end do
      end do
      ainv = amat
      call zgetrf(nd, nd, ainv, nd, ipiv, info)
      call zgetri(nd, ainv, nd, ipiv, work, size(work), info)
      do i = 1, nsite
        spc = merge(1, 2, isa(i))
        do l = 1, 9
          io = cls(l)
          if (spc == 1) then; gref = l1(l) + m1(l) * ainv(9 * (i - 1) + l, 9 * (i - 1) + l)
          else; gref = l2(l) + m2(l) * ainv(9 * (i - 1) + l, 9 * (i - 1) + l); end if
          if (io == 1) then; gex(spc, 1) = gex(spc, 1) + gref; nrm(spc, 1) = nrm(spc, 1) + 1.0_dp; end if
          if (io == 2) then; gex(spc, 2) = gex(spc, 2) + gref; nrm(spc, 2) = nrm(spc, 2) + 1.0_dp; end if
          if (io == 3 .and. (l == 5 .or. l == 6 .or. l == 8)) then; gex(spc, 3) = gex(spc, 3) + gref; nrm(spc, 3) = nrm(spc, 3) + 1.0_dp; end if
          if (io == 3 .and. (l == 7 .or. l == 9)) then; gex(spc, 4) = gex(spc, 4) + gref; nrm(spc, 4) = nrm(spc, 4) + 1.0_dp; end if
        end do
      end do
    end do
    do io = 1, 4; do spc = 1, 2; gex(spc, io) = gex(spc, io) / nrm(spc, io); end do; end do
    do io = 1, 4
      call augsym_set_start(sop, iorb(io), d)
      call ed_species_green(sop, wk, z, spa, spb, alpha9, x, iorb(io), nlev, sop%flipr(1), sop%md(1), sop%s, ga, gb, gavg)
      r = max(abs(ga - gex(1, io)) / abs(gex(1, io)), abs(gb - gex(2, io)) / abs(gex(2, io)))
      err2 = max(err2, r)
      if (ordered) print '(a,2f10.5,a,2f10.5)', '         g_alpha asr ', gavg, '  exact(first site/orb) ', ainv(iorb(io), iorb(io))
      print '(a,f6.2,a,i1,a,2f10.5,a,2f10.5,a,es9.2)', '    E=', za(ie), ' orb', io, '  A: asr ', ga, '  exact ', gex(1, io), '   rel dev ', r
    end do
  end do
  if (ordered) then
    call check(err2 < 0.02_dp, 'ordered limit: G (8 levels) within 2% of exact inverse of P(z) - S')
  else
    call check(err2 < 0.30_dp, 'disorder: G (3 levels, rank 4) within 30% of ensemble-averaged exact inverses (depth-limited)')
  end if
  call report()
contains
  subroutine load_params(sp, pp)
    type(ed_params_t), intent(out) :: sp
    type(potpar_t), intent(in) :: pp
    integer :: ii
    do ii = 1, 9
      sp%c(ii) = pp%c(orb_l(ii)); sp%del(ii) = pp%sqrdel(orb_l(ii))**2; sp%gam(ii) = pp%gamma(orb_l(ii))
    end do
  end subroutine load_params
  subroutine random_species(flag, conc)
    logical, intent(out) :: flag(:)
    real(dp), intent(in) :: conc
    real(dp) :: u(size(flag))
    call random_number(u)
    flag = u < conc
  end subroutine random_species
end program test_endep
