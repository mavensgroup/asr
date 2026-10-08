!> term_bench <sc|bcc|fcc|b2> [nk=128] [eta=0.1] [delta=3]
!> b2 = ordered CsCl compound on bcc, on-site +-delta (gapped, period-2 coefficients).
!> Terminator benchmark on an ordered NN band: for every (asymptote, join,
!> depth) report max |LDOS_rec - LDOS_kspace| over the band, k-space being the
!> independent truth.  Output goes to stdout; redirect into trial/ to keep it.
program term_bench
  use asr_kinds, only: dp, pi
  use asr_operator, only: csr_op
  use asr_lattice_toy, only: cubic_cluster
  use asr_recursion, only: lanczos
  use asr_kspace, only: kspace_green, kspace_green_b2, coordination
  use asr_terminator, only: terminator_t, make_terminator, terminated_green
  implicit none
  character(len=16) :: kind, arg
  integer, parameter :: ne = 121, nmax = 60
  integer :: nk, centre, m, ie, id, im, ij, zc, nd, nint
  integer, parameter :: depths(4) = [8, 16, 32, 60]
  character(len=13), parameter :: methods(4) = ['sqrt-last    ', 'sqrt-avg2    ', 'beer-pettifor', 'periodic2    ']
  character(len=11), parameter :: joins(2) = ['butt       ', 'luchini-nex']
  real(dp) :: eta, delta, egrid(ne), a(nmax), b(nmax), err, emax
  real(dp), allocatable :: v(:)
  complex(dp) :: gref(ne), g, z
  type(csr_op) :: op
  type(terminator_t) :: t

  call get_command_argument(1, kind)
  nk = 128; eta = 0.1_dp
  if (command_argument_count() >= 2) then; call get_command_argument(2, arg); read (arg, *) nk; end if
  if (command_argument_count() >= 3) then; call get_command_argument(3, arg); read (arg, *) eta; end if
  delta = 3.0_dp
  if (command_argument_count() >= 4) then; call get_command_argument(4, arg); read (arg, *) delta; end if
  if (trim(kind) == 'b2') then; zc = 8; else; zc = coordination(trim(kind)); end if
  emax = 1.2_dp * zc
  if (trim(kind) == 'b2') emax = 1.1_dp * sqrt(delta**2 + 64.0_dp)
  do ie = 1, ne
    egrid(ie) = -emax + 2 * emax * (ie - 1) / real(ne - 1, dp)
  end do
  if (trim(kind) == 'b2') then
    call cubic_cluster('bcc', 2 * nmax + 6, 1.0_dp, op, centre, b2_delta=delta)
  else
    call cubic_cluster(trim(kind), 2 * nmax + 6, 1.0_dp, op, centre)
  end if
  allocate (v(op%n)); v = 0.0_dp; v(centre) = 1.0_dp
  call lanczos(op, v, nmax, a, b, m)
  if (trim(kind) == 'b2') then
    gref = kspace_green_b2(1.0_dp, delta, egrid, eta, nk)
  else
    gref = kspace_green(trim(kind), 1.0_dp, egrid, eta, nk)
  end if
  print '(a,a,a,i0,a,f6.3)', '# lattice=', trim(kind), ' nk=', nk, ' eta=', eta
  print '(a)', '# max |dLDOS| vs k-space;  columns = depth 8 16 32 60 (even depths: period-2 safe)'
  do im = 1, 4
    do ij = 1, 2
      write (*, '(a13,1x,a11)', advance='no') methods(im), joins(ij)
      do id = 1, 4
        nd = depths(id); nint = max(4, nd / 2)
        t = make_terminator(trim(methods(im)), trim(joins(ij)), a, b, nd, nint)
        err = 0.0_dp
        do ie = 1, ne
          z = cmplx(egrid(ie), eta, dp)
          g = terminated_green(t, z)
          err = max(err, abs(aimag(g) - aimag(gref(ie))) / pi)
        end do
        write (*, '(es11.3)', advance='no') err
      end do
      write (*, *)
    end do
  end do
end program term_bench
