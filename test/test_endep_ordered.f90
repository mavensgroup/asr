!> Energy-dependent recursion for an ORDERED B2 crystal (no augmented space, asr_endep::ed_ordered_green) vs the exact dense
!> inverse of A(z) = P(z) - S on the same periodic supercell; every start orbital on both sublattices.
program test_endep_ordered
  use asr_kinds, only: dp
  use asr_lmto_ham, only: potpar_t
  use asr_augspace, only: graph_t
  use asr_augorb, only: orb_alloy_t
  use asr_endep
  use asr_lmto_alloy
  use asr_lmto_realization
  use asr_strconst, only: hardcore_alpha
  use asr_testing
  implicit none
  type(potpar_t) :: pa, pb
  type(graph_t) :: g
  type(orb_alloy_t) :: mdl
  type(lmto_real_op) :: rop
  type(ed_params_t) :: spa, spb
  type(ed_work_t) :: wk
  integer, allocatable :: coord(:, :), ipiv(:)
  logical, allocatable :: isa(:)
  complex(dp), allocatable :: amat(:, :), ainv(:, :), work(:)
  real(dp) :: alpha(0:2), al4(0:3), alpha9(9), err
  complex(dp) :: z, p9a(9), p9b(9), m1(9), m2(9), l1(9), l2(9), gr, gex
  integer, parameter :: orb_l(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2], iorb(4) = [1, 3, 5, 7]
  integer :: l, i, j, k, kc, nd, info, nlev, ie, io, ip, site, ls
  real(dp), parameter :: za(3) = [-0.1_dp, 0.1_dp, 0.25_dp]
  character(len=8) :: arg

  print '(a)', 'test_endep_ordered'
  ! third-order potential function: mu^2 = dP/dz and lambda = -P''/(2P') by complex finite differences
  block
    type(ed_params_t) :: sq
    complex(dp) :: z0, h, pp(9), pm(9), pc(9), m2(9), lm(9), dum(9), dd(9)
    real(dp) :: al9(9), e1, e2
    integer :: ii
    do ii = 1, 9
      sq%c(ii) = -0.25_dp + 0.03_dp * ii; sq%del(ii) = 0.02_dp + 0.01_dp * ii; sq%gam(ii) = 0.05_dp - 0.01_dp * ii
      sq%p(ii) = 0.3_dp * ii; sq%enu(ii) = -0.2_dp + 0.01_dp * ii; al9(ii) = 0.02_dp * mod(ii, 3)
    end do
    z0 = cmplx(-0.1_dp, 0.07_dp, dp); h = (1.0e-4_dp, 0.0_dp)
    call ed_functions(z0, sq, al9, pc, m2, lm)
    call ed_functions(z0 + h, sq, al9, pp, dum, dd)
    call ed_functions(z0 - h, sq, al9, pm, dum, dd)
    e1 = maxval(abs((pp - pm) / (2.0_dp * h) - m2) / abs(m2))
    e2 = maxval(abs(-0.5_dp * ((pp - 2.0_dp * pc + pm) / h**2) / ((pp - pm) / (2.0_dp * h)) - lm) / abs(lm))
    print '(a,es10.2,a,es10.2)', '    third-order P(z): rel err mu^2 vs dP/dz ', e1, '   lambda vs -P''''/2P'' ', e2
    call check(e1 < 1.0e-6_dp .and. e2 < 1.0e-4_dp, 'third-order potential function: mu^2 = P'' and lambda = -P''''/(2P'')')
  end block
  do l = 0, 2
    pa%enu(l) = 0.1_dp * l; pa%c(l) = 0.3_dp - 0.15_dp * l; pa%sqrdel(l) = 0.30_dp - 0.1_dp * l
    pa%p(l) = 0.9_dp - 0.1_dp * l; pa%gamma(l) = 0.25_dp - 0.1_dp * l; pa%alpha_rep(l) = pa%gamma(l)
    pb%enu(l) = 0.05_dp * l; pb%c(l) = 0.45_dp - 0.1_dp * l; pb%sqrdel(l) = 0.28_dp - 0.1_dp * l
    pb%p(l) = 0.8_dp - 0.1_dp * l; pb%gamma(l) = 0.2_dp - 0.08_dp * l; pb%alpha_rep(l) = pb%gamma(l)
  end do
  al4 = hardcore_alpha(0.7_dp); alpha = al4(0:2)
  do i = 1, 9; alpha9(i) = alpha(orb_l(i)); end do
  do i = 1, 9
    spa%c(i) = pa%c(orb_l(i)); spa%del(i) = pa%sqrdel(orb_l(i))**2; spa%gam(i) = pa%gamma(orb_l(i))
    spb%c(i) = pb%c(orb_l(i)); spb%del(i) = pb%sqrdel(orb_l(i))**2; spb%gam(i) = pb%gamma(orb_l(i))
  end do
  nlev = 10; ls = 4
  if (command_argument_count() > 0) then; call get_command_argument(1, arg); read (arg, *) nlev; end if
  if (command_argument_count() > 1) then; call get_command_argument(2, arg); read (arg, *) ls; end if
  call build_lmto_alloy('bcc', ls, 1.2_dp, 0.7_dp, pa, pb, 0.5_dp, g, coord, mdl)
  allocate (isa(g%n)); isa = modulo(coord(1, :), 2) == 0
  call build_lmto_ordered(g, mdl, isa, rop)
  nd = 9 * g%n
  print '(a,i0,a,i0)', '    sites ', g%n, '  levels ', nlev
  call ed_work_init(wk, nd)
  allocate (amat(nd, nd), ainv(nd, nd), ipiv(nd), work(64 * nd))
  err = 0.0_dp
  do ie = 1, size(za)
    z = cmplx(za(ie), 0.3_dp, dp)
    call ed_functions(z, spa, alpha9, p9a, m1, l1)
    call ed_functions(z, spb, alpha9, p9b, m2, l2)
    amat = 0.0_dp
    do i = 1, g%n
      do l = 1, 9
        amat(9 * (i - 1) + l, 9 * (i - 1) + l) = merge(p9a(l), p9b(l), isa(i))
      end do
      do l = 1, 9; do j = 1, 9
        amat(9 * (i - 1) + l, 9 * (i - 1) + j) = amat(9 * (i - 1) + l, 9 * (i - 1) + j) - mdl%block0(l, j)
      end do; end do
      do k = 1, size(mdl%block, 3)
        kc = g%nbr(g%nbr_ptr(i) + k - 1)
        do l = 1, 9; do j = 1, 9
          amat(9 * (i - 1) + l, 9 * (kc - 1) + j) = amat(9 * (i - 1) + l, 9 * (kc - 1) + j) - mdl%block(l, j, k)
        end do; end do
      end do
    end do
    ainv = amat
    call zgetrf(nd, nd, ainv, nd, ipiv, info)
    call zgetri(nd, ainv, nd, ipiv, work, size(work), info)
    do ip = 1, 2
      site = 1
      if (ip == 2) then
        do i = 1, g%n; if (.not. isa(i)) then; site = i; exit; end if; end do
      end if
      do io = 1, 4
        l = iorb(io)
        if (ip == 1) then; gex = l1(l) + m1(l) * ainv(9 * (site - 1) + l, 9 * (site - 1) + l)
        else; gex = l2(l) + m2(l) * ainv(9 * (site - 1) + l, 9 * (site - 1) + l); end if
        call ed_ordered_green(rop, wk, z, spa, spb, alpha9, site, l, nlev, gr)
        err = max(err, abs(gr - gex) / abs(gex))
        print '(a,f6.2,a,i1,a,i1,a,2f10.5,a,2f10.5)', '    E=', za(ie), ' sub', ip, ' orb', io, '  asr ', gr, '  exact ', gex
      end do
    end do
  end do
  print '(a,es10.2)', '    max relative deviation: ', err
  call check(err < 0.02_dp, 'ordered B2 ED recursion within 2% of the exact inverse of P(z) - S')
  call report()
end program test_endep_ordered
