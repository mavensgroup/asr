!> probe_flip  <x> <delta> <d1> <maxdepth> [eta] [L_truth]   (scans a grid of rank cutoffs and flip radii)
!> Single-orbital H = E + c + d S d  (rank-1 random hopping) on the sc lattice in the augorb machinery:
!> ASR with rank cutoff, flip radius rf (graph distance from the start; beyond it the medium is the
!> virtual average), depth n  vs KPM truth on random realisations.  Error = max|dLDOS|/peak.
program probe_flip
  use asr_kinds, only: dp, pi
  use asr_augspace, only: graph_t, graph_distances
  use asr_augorb
  use asr_lattices, only: periodic_graph_shells
  use asr_lmto_realization
  use asr_recursion, only: lanczos
  use asr_terminator
  use asr_kpm
  implicit none
  type(graph_t) :: g, gt
  type(orb_alloy_t) :: mdl
  type(augorb_op) :: op
  type(lmto_real_op) :: rop
  integer, allocatable :: coord(:, :), bvec(:, :), dist(:)
  logical, allocatable :: flippable(:)
  real(dp), allocatable :: v0(:), a(:), b(:), mu(:)
  real(dp) :: x, delta, d1, eta, scale, err, c0, s0
  integer :: rank, rf, maxdepth, lt, depth, irk, irf, ie, im, m, nd, ir, nreal, nvec, nmom, ns, k
  integer, parameter :: ranks(3) = [2, 3, 4], rfs(4) = [1, 2, 3, 4]
  integer, parameter :: ne = 801
  real(dp) :: egrid(ne), rho(ne), rhosum(ne), rhob(ne), rhob2(ne), rhosum2(ne)
  complex(dp) :: z, gg
  type(terminator_t) :: t
  character(len=16) :: arg
  integer :: iarg
  iarg = command_argument_count()
  call get_command_argument(1, arg); read (arg, *) x
  call get_command_argument(2, arg); read (arg, *) delta
  call get_command_argument(3, arg); read (arg, *) d1
  call get_command_argument(4, arg); read (arg, *) maxdepth
  eta = 0.15_dp; lt = 36
  if (iarg >= 5) then; call get_command_argument(5, arg); read (arg, *) eta; end if
  if (iarg >= 6) then; call get_command_argument(6, arg); read (arg, *) lt; end if

  ! model: eps = c0(n) with delta; hopping t = -d_R d_R'  (d = 1 + d1 n)
  mdl%nl = 1; mdl%x = x
  mdl%c0 = [-0.5_dp * delta]; mdl%c1 = [delta]
  mdl%d0 = [1.0_dp]; mdl%d1 = [d1]
  mdl%o0 = [0.0_dp]; mdl%o1 = [0.0_dp]; mdl%e0 = [0.0_dp]; mdl%e1 = [0.0_dp]
  call periodic_graph_shells('sc', lt, 4, gt, coord, bvec)
  ns = size(bvec, 2)
  allocate (mdl%block(1, 1, ns)); mdl%block = -1.0_dp
  allocate (mdl%block0(1, 1)); mdl%block0 = 0.0_dp
  ! ---- truth
  c0 = 0.0_dp; s0 = 1.05_dp * (6.0_dp * (1.0_dp + abs(d1))**2 + delta); nmom = 800; nvec = 40; nreal = 4
  do ie = 1, ne; egrid(ie) = -s0 + 2 * s0 * (ie - 1) / real(ne - 1, dp); end do
  allocate (mu(0:nmom - 1)); rhosum = 0.0_dp; rhosum2 = 0.0_dp
  do ir = 1, nreal
    call build_lmto_real(gt, mdl, 500 + ir, rop)
    call kpm_moments(rop, c0, s0, nmom, nvec, 21 + ir, mu); call kpm_dos(mu, c0, s0, egrid, rho); rhosum = rhosum + rho
    call build_lmto_real(gt, mdl, 900 + ir, rop)
    call kpm_moments(rop, c0, s0, nmom, nvec, 51 + ir, mu); call kpm_dos(mu, c0, s0, egrid, rho); rhosum2 = rhosum2 + rho
  end do
  call lorentz_broaden(egrid, rhosum / nreal, eta, rhob); call lorentz_broaden(egrid, rhosum2 / nreal, eta, rhob2)
  scale = maxval(rhob)
  print '(a,es10.3,a,i0,a)', '# truth noise = ', maxval(abs(rhob - rhob2)) / scale, '   (', gt%n * nreal * 2, ' sites)'
  rhob = 0.5_dp * (rhob + rhob2)
  print '(a,f5.2,a,f6.2,a,f5.2,a,f6.3)', '# x=', x, ' delta=', delta, ' d1=', d1, ' eta=', eta
  call periodic_graph_shells('sc', maxdepth + 3, 4, g, coord, bvec)
  call graph_distances(g, 1, dist)
  allocate (flippable(g%n))
  do irk = 1, size(ranks)
    do irf = 1, size(rfs)
      rank = ranks(irk); rf = rfs(irf)
      flippable = dist <= rf
      print '(a,i0,a,i0)', '# --- rank <= ', rank, '   flip radius = ', rf
      print '(a)', '# depth   nodes  | max|dDOS|/peak : sqrt-last / sqrt-avg2 / BP+LN'
      do depth = 2, maxdepth
        call build_augorb_operator(g, mdl, 1, depth, op, rank, flippable)
        nd = op%dim()
        allocate (a(depth), b(depth), v0(nd)); v0 = 0.0_dp; v0(1) = 1.0_dp
        call lanczos(op, v0, depth, a, b, m, reorth=.true.)
        write (*, '(i6,i10,3x)', advance='no') depth, op%nn
        do im = 1, 3
          select case (im)
          case (1); t = make_terminator('sqrt-last', 'butt', a, b, depth, max(4, depth / 2))
          case (2); t = make_terminator('sqrt-avg2', 'butt', a, b, depth, max(4, depth / 2))
          case (3); t = make_terminator('beer-pettifor', 'luchini-nex', a, b, depth, max(4, depth / 2))
          end select
          err = 0.0_dp
          do ie = 1, ne
            z = cmplx(egrid(ie), eta, dp); gg = terminated_green(t, z)
            err = max(err, abs(-aimag(gg) / pi - rhob(ie)))
          end do
          write (*, '(es11.3)', advance='no') err / scale
        end do
        write (*, *)
        deallocate (a, b, v0)
        if (op%nn > 4000000) exit
      end do
    end do
  end do
end program probe_flip
