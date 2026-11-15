!> probe_asr_depth [lattice=sc] [x=0.5] [delta=3.0] [offdiag=0] [maxdepth=6] [eta=0.15]
!> Depth-limited augmented-space recursion on a disordered binary alloy vs KPM truth (large supercell),
!> for every terminator.  Error = max|rho_ASR - rho_KPM| / max rho_KPM on the energy window.
program probe_asr_depth
  use asr_kinds, only: dp, pi
  use asr_operator, only: csr_op
  use asr_augspace
  use asr_lattices
  use asr_symmetry
  use asr_recursion, only: lanczos
  use asr_terminator
  use asr_kpm
  implicit none
  type(graph_t) :: ga, gk
  integer, allocatable :: coord(:, :), perm(:, :)
  type(binary_model_t) :: mdl
  type(csr_op) :: op, hk
  character(len=16) :: kind, arg
  real(dp) :: x, delta, eta, offd, emax, c0, s0, err, scale
  integer :: maxdepth, depth, ie, im, ij, nst, m, nmom, nvec, nreal, ir, i
  integer, parameter :: ne = 1001
  real(dp) :: egrid(ne), rhok(ne), rhokb(ne), rhosum(ne), mu(0:799), rho_asr(ne), rhok2(ne), rhokb2(ne), rhosum2(ne)
  real(dp), allocatable :: a(:), b(:), v0(:), mus(:)
  complex(dp) :: z, g
  type(terminator_t) :: t
  character(len=13), parameter :: methods(4) = ['sqrt-last    ', 'sqrt-avg2    ', 'beer-pettifor', 'periodic2    ']
  character(len=11), parameter :: joins(2) = ['butt       ', 'luchini-nex']
  integer :: iarg, maxrank
  kind = 'sc'; x = 0.5_dp; delta = 3.0_dp; offd = 0.0_dp; maxdepth = 6; eta = 0.15_dp
  iarg = command_argument_count()
  if (iarg >= 1) call get_command_argument(1, kind)
  if (iarg >= 2) then; call get_command_argument(2, arg); read (arg, *) x; end if
  if (iarg >= 3) then; call get_command_argument(3, arg); read (arg, *) delta; end if
  if (iarg >= 4) then; call get_command_argument(4, arg); read (arg, *) offd; end if
  if (iarg >= 5) then; call get_command_argument(5, arg); read (arg, *) maxdepth; end if
  if (iarg >= 6) then; call get_command_argument(6, arg); read (arg, *) eta; end if
  maxrank = huge(1)
  if (iarg >= 7) then; call get_command_argument(7, arg); read (arg, *) maxrank; end if
  mdl%x = x; mdl%epsa = 0.5_dp * delta; mdl%epsb = -0.5_dp * delta
  mdl%tAA = -1.0_dp - offd; mdl%tAB = -1.0_dp; mdl%tBB = -1.0_dp + offd
  emax = 14.0_dp
  do ie = 1, ne; egrid(ie) = -emax + 2 * emax * (ie - 1) / real(ne - 1, dp); end do

  ! ---- truth: KPM on L^3 supercells, averaged over realisations
  call periodic_graph(trim(kind), 48, gk)
  c0 = 0.0_dp; s0 = 1.05_dp * (12.0_dp + 0.5_dp * delta + abs(offd) * 12.0_dp)
  nreal = 6; nvec = 60; nmom = 800
  rhosum = 0.0_dp; rhosum2 = 0.0_dp
  do ir = 1, nreal
    call random_alloy_op(gk, mdl, 100 + ir, hk)
    call kpm_moments(hk, c0, s0, nmom, nvec, 7 + ir, mu)
    call kpm_dos(mu, c0, s0, egrid, rhok)
    rhosum = rhosum + rhok
    call random_alloy_op(gk, mdl, 200 + ir, hk)                 ! independent second truth set (noise estimate)
    call kpm_moments(hk, c0, s0, nmom, nvec, 70 + ir, mu)
    call kpm_dos(mu, c0, s0, egrid, rhok2)
    rhosum2 = rhosum2 + rhok2
  end do
  rhok = rhosum / nreal; rhok2 = rhosum2 / nreal
  call lorentz_broaden(egrid, rhok, eta, rhokb)
  call lorentz_broaden(egrid, rhok2, eta, rhokb2)
  scale = maxval(rhokb)
  print '(a,es10.3)', 'KPM truth noise (two independent sets): max|diff|/peak = ', maxval(abs(rhokb - rhokb2)) / scale
  rhokb = 0.5_dp * (rhokb + rhokb2)
  print '(a,a,a,f5.2,a,f6.2,a,f5.2,a,f6.3,a,i0,a)', 'lattice=', trim(kind), '  x=', x, '  delta=', delta, '  offdiag=', offd, '  eta=', eta, '   (KPM truth: ', gk%n * nreal, ' sites)   [point-group reduced states]'
  ! ---- ASR
  call periodic_graph(trim(kind), 2 * maxdepth + 4, ga, coord)
  call cubic_site_perms(coord, 2 * (2 * maxdepth + 4), perm)
  print '(a)', 'depth  states   | max|dDOS|/peak  for:  sqrt-last / sqrt-avg2 / BP / BP+LN / periodic2'
  do depth = 2, maxdepth
    call build_augmented_operator_sym(ga, mdl, perm, 1, depth, op, maxrank)
    nst = op%n
    allocate (a(depth), b(depth), v0(nst)); v0 = 0.0_dp; v0(1) = 1.0_dp
    call lanczos(op, v0, depth, a, b, m, reorth=.true.)
    write (*, '(i4,i9,3x)', advance='no') depth, nst
    do im = 1, 5
      select case (im)
      case (1); t = make_terminator('sqrt-last', 'butt', a, b, depth, max(4, depth / 2))
      case (2); t = make_terminator('sqrt-avg2', 'butt', a, b, depth, max(4, depth / 2))
      case (3); t = make_terminator('beer-pettifor', 'butt', a, b, depth, max(4, depth / 2))
      case (4); t = make_terminator('beer-pettifor', 'luchini-nex', a, b, depth, max(4, depth / 2))
      case (5)
        if (mod(depth, 2) == 0) then
          t = make_terminator('periodic2', 'butt', a, b, depth, 0)
        else
          write (*, '(a11)', advance='no') '   n/a'; cycle
        end if
      end select
      err = 0.0_dp
      do ie = 1, ne
        z = cmplx(egrid(ie), eta, dp)
        g = terminated_green(t, z)
        err = max(err, abs(-aimag(g) / pi - rhokb(ie)))
      end do
      write (*, '(es11.3)', advance='no') err / scale
    end do
    write (*, *)
    deallocate (a, b, v0)
  end do
end program probe_asr_depth
