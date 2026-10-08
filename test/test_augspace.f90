!> Augmented space (binary alloy, single orbital) vs brute-force enumeration of ALL 2^N configurations
!> on a small cluster.  With the Krylov space exhausted the recursion on the augmented operator must equal
!> the exact configuration average  <<G_00(z)>> = sum_C p(C) G_00[C](z)  to machine precision.
!> Diagonal AND off-diagonal disorder, asymmetric concentration.
program test_augspace
  use asr_kinds, only: dp
  use asr_operator, only: csr_op
  use asr_augspace
  use asr_recursion, only: lanczos
  use asr_greens, only: cf_green
  use asr_testing
  implicit none
  integer, parameter :: nside = 2, n = nside**3
  type(graph_t) :: g
  type(binary_model_t) :: mdl
  type(csr_op) :: op
  real(dp) :: hmat(n, n), ev(n), work(8 * n), w, pc, exact_err, sym
  real(dp), allocatable :: a(:), b(:), dense(:, :), v0(:), e1(:)
  complex(dp) :: z, gex, grec
  integer :: i, j, k, kz, c, info, m, ix, iy, iz, jx, jy, jz, nst
  logical :: isa(n)
  external :: dsyev

  print '(a)', 'test_augspace'
  ! 2x2x2 open cubic cluster
  allocate (g%nbr_ptr(n + 1), g%nbr(3 * n)); g%n = n
  k = 0; g%nbr_ptr(1) = 1
  do i = 1, n
    ix = mod(i - 1, nside); iy = mod((i - 1) / nside, nside); iz = (i - 1) / nside**2
    do j = 1, n
      jx = mod(j - 1, nside); jy = mod((j - 1) / nside, nside); jz = (j - 1) / nside**2
      if (abs(ix - jx) + abs(iy - jy) + abs(iz - jz) == 1) then; k = k + 1; g%nbr(k) = j; end if
    end do
    g%nbr_ptr(i + 1) = k + 1
  end do
  mdl%x = 0.3_dp; mdl%epsa = 0.6_dp; mdl%epsb = -0.4_dp
  mdl%tAA = -1.0_dp; mdl%tAB = -0.7_dp; mdl%tBB = -0.45_dp

  call build_augmented_operator(g, mdl, 1, 60, op)
  nst = op%n
  print '(a,i0,a,i0)', '    augmented states: ', nst, '   (N 2^N = ', n * 2**n, ')'
  call check(nst == n * 2**n, 'all N 2^N augmented states reached')
  ! operator symmetric
  allocate (dense(nst, nst)); dense = 0.0_dp
  do i = 1, nst
    do k = op%ia(i), op%ia(i + 1) - 1
      dense(i, op%ja(k)) = dense(i, op%ja(k)) + op%a(k)
    end do
  end do
  sym = maxval(abs(dense - transpose(dense)))
  call check_close(sym, 0.0_dp, 1.0e-13_dp, 'augmented operator is symmetric')

  allocate (a(nst), b(nst), v0(nst)); v0 = 0.0_dp; v0(1) = 1.0_dp
  call lanczos(op, v0, nst, a, b, m, reorth=.true., tiny=1.0e-9_dp)
  print '(a,i0)', '    Krylov dimension used: ', m

  exact_err = 0.0_dp
  do kz = 1, 12
    z = cmplx(-3.0_dp + 0.5_dp * kz, 0.15_dp, dp)
    ! brute force over all 2^N configurations: bit set = site is A
    gex = (0.0_dp, 0.0_dp)
    do c = 0, 2**n - 1
      pc = 1.0_dp
      do i = 1, n
        isa(i) = btest(c, i - 1)
        if (isa(i)) then; pc = pc * mdl%x; else; pc = pc * (1.0_dp - mdl%x); end if
      end do
      hmat = 0.0_dp
      do i = 1, n
        if (isa(i)) then; hmat(i, i) = mdl%epsa; else; hmat(i, i) = mdl%epsb; end if
        do k = g%nbr_ptr(i), g%nbr_ptr(i + 1) - 1
          j = g%nbr(k)
          if (isa(i) .and. isa(j)) then; hmat(i, j) = mdl%tAA
          else if (isa(i) .neqv. isa(j)) then; hmat(i, j) = mdl%tAB
          else; hmat(i, j) = mdl%tBB; end if
        end do
      end do
      call dsyev('V', 'U', n, hmat, n, ev, work, size(work), info)
      w = 0.0_dp
      do j = 1, n
        gex = gex + pc * hmat(1, j)**2 / (z - ev(j))
      end do
    end do
    grec = cf_green(a, b, m, z, (0.0_dp, 0.0_dp))
    exact_err = max(exact_err, abs(gex - grec))
    e1 = [0.0_dp]
  end do
  call check_close(exact_err, 0.0_dp, 1.0e-9_dp, 'augmented-space recursion == exact configuration average (all 256 configs)')
  call report()
end program test_augspace
