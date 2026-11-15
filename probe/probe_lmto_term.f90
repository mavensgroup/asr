!> Terminator benchmark on the real Ag d/sp bands (ordered fcc, Stuttgart parameters):
!> Lanczos from |0,L> on the k-grid operator, every (asymptote, join, depth) vs the exact
!> orbital Green function on the same grid.  Error metric: max |dLDOS| / max LDOS over the window.
program probe_lmto_term
  use asr_kinds, only: dp, pi
  use asr_strconst
  use asr_lmto_ham
  use asr_recursion, only: lanczos
  use asr_terminator
  implicit none
  integer, parameter :: ne = 301, nmax = 60
  type(potpar_t) :: pp
  type(kgrid_op) :: op
  real(dp), allocatable :: vec(:, :), blk(:, :, :), v(:)
  real(dp) :: alpha(0:2), al4(0:3), w, bvec(3, 3), egrid(ne), a(nmax), b(nmax), eta, err, scale
  complex(dp) :: gex(ne), g, z(ne)
  integer :: nk, orbs(4), io, im, ij, id, ie, m, nd, ni
  character(len=5) :: lab(4) = ['s    ', 'p_x  ', 't2g  ', 'eg   ']
  integer, parameter :: depths(4) = [10, 20, 40, 60]
  character(len=13), parameter :: methods(4) = ['sqrt-last    ', 'sqrt-avg2    ', 'beer-pettifor', 'periodic2    ']
  character(len=11), parameter :: joins(2) = ['butt       ', 'luchini-nex']
  type(terminator_t) :: t
  character(len=16) :: arg

  nk = 24; eta = 0.02_dp
  if (command_argument_count() >= 1) then; call get_command_argument(1, arg); read (arg, *) nk; end if
  if (command_argument_count() >= 2) then; call get_command_argument(2, arg); read (arg, *) eta; end if
  pp%enu = [-0.54494317_dp, -0.45471939_dp, -0.51341327_dp]
  pp%c = [-0.46144639_dp, 0.41368602_dp, -0.51457447_dp]
  pp%sqrdel = [0.36457406_dp, -0.36517716_dp, -0.10341022_dp]
  pp%p = [0.08230178_dp, 0.03623364_dp, 2.08144602_dp]
  pp%gamma = [0.42477526_dp, 0.11537965_dp, 0.00550351_dp]
  pp%alpha_rep = pp%gamma
  w = avg_ws_radius('fcc'); al4 = hardcore_alpha(0.7_dp); alpha = al4(0:2)
  call lattice_cluster('fcc', 1.55_dp, vec)
  call screened_cluster(vec, w, alpha, blk)
  bvec = reshape([-1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp], [3, 3])
  call hk_grid(nk, bvec, vec, blk, pp, alpha, op)
  do ie = 1, ne; egrid(ie) = -0.9_dp + 1.3_dp * (ie - 1) / (ne - 1); z(ie) = cmplx(egrid(ie), eta, dp); end do
  orbs = [1, 2, 5, 7]
  print '(a,i0,a,f6.3,a)', '# Ag fcc LMTO, nk=', nk, ' eta=', eta, ' Ry;  window -0.9..0.4 Ry (occupied region)'
  print '(a)', '# err = max|dLDOS| / max LDOS ; columns depth 10 20 40 60'
  allocate (v(op%dim()))
  do io = 1, 4
    gex = green_exact_orbital(op, orbs(io), z)
    scale = maxval(-aimag(gex) / pi)
    v = 0.0_dp
    v(orbs(io):9 * op%nk3:9) = 1.0_dp            ! real part, orbital `orbs(io)` at every k
    call lanczos(op, v, nmax, a, b, m, reorth=.false.)
    print '(a,a,a,f8.3)', 'orbital ', lab(io), '  peak LDOS =', scale
    do im = 1, 4
      do ij = 1, 2
        write (*, '(a13,1x,a11)', advance='no') methods(im), joins(ij)
        do id = 1, 4
          nd = depths(id); ni = max(4, nd / 2)
          if (mod(nd, 2) /= 0 .and. methods(im) == 'periodic2    ') then; write (*, '(a11)', advance='no') '    n/a'; cycle; end if
          t = make_terminator(trim(methods(im)), trim(joins(ij)), a, b, nd, ni)
          err = 0.0_dp
          do ie = 1, ne
            g = terminated_green(t, z(ie))
            err = max(err, abs(aimag(g) - aimag(gex(ie))) / pi)
          end do
          write (*, '(es11.3)', advance='no') err / scale
        end do
        write (*, *)
      end do
    end do
  end do
end program probe_lmto_term
