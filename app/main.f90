!> asr <input.toml>
!> Reads one validated TOML record (any directory) and runs the calculation.
!> Output files are written next to the input file.
program asr_main
  use asr_kinds, only: dp
  use asr_input, only: input_t, read_input
  use asr_operator, only: csr_op
  use asr_lattice_toy, only: cubic_cluster
  use asr_recursion, only: lanczos
  use asr_greens, only: cf_green, sqrt_tail, tail_params_avg2, tail_params_last, ldos_from_g
  use asr_input_lmto, only: lmto_input_t, read_lmto_input, peek_kind
  use asr_lmto_driver, only: run_lmto
  use asr_alloy_driver, only: alloy_run_input_t, read_alloy_input, run_alloy
  implicit none
  type(input_t) :: inp
  type(lmto_input_t) :: linp
  type(alloy_run_input_t) :: ainp
  logical :: lok
  character(len=:), allocatable :: msg, infile, outfile
  character(len=512) :: arg
  integer :: centre, m, ie, u
  real(dp) :: e, ainf, binf
  real(dp), allocatable :: a(:), b(:), v(:)
  type(csr_op) :: op
  complex(dp) :: z, g

  if (command_argument_count() /= 1) then
    print '(a)', 'usage: asr <input.toml>'; stop 2
  end if
  call get_command_argument(1, arg); infile = trim(arg)
  if (peek_kind(infile) == 'lmto') then
    call read_lmto_input(infile, linp, msg)
    if (len(msg) > 0) then; print '(a)', 'asr: input error: '//msg; stop 1; end if
    call run_lmto(linp, lok)
    if (.not. lok) stop 3
    stop 0
  end if
  if (peek_kind(infile) == 'alloy') then
    call read_alloy_input(infile, ainp, msg)
    if (len(msg) > 0) then; print '(a)', 'asr: input error: '//msg; stop 1; end if
    call run_alloy(ainp, lok)
    if (.not. lok) stop 3
    stop 0
  end if
  call read_input(infile, inp, msg)
  if (len(msg) > 0) then
    print '(a)', 'asr: input error: '//msg; stop 1
  end if
  print '(a)', 'asr 0.1.0  '//inp%name
  print '(a)', 'input: '//infile

  call cubic_cluster(inp%lattice, 2 * inp%nlevels + 6, inp%hopping, op, centre)
  allocate (a(inp%nlevels), b(inp%nlevels), v(op%n)); v = 0.0_dp; v(centre) = 1.0_dp
  call lanczos(op, v, inp%nlevels, a, b, m)
  if (inp%terminator == 'sqrt-last') then
    call tail_params_last(a, b, inp%nlevels, ainf, binf)
  else
    call tail_params_avg2(a, b, inp%nlevels, ainf, binf)
  end if

  outfile = inp%workdir//inp%dos_file
  open (newunit=u, file=outfile, status='replace', action='write')
  write (u, '(a)') '# '//inp%name
  write (u, '(a,i0,a,a,a,f7.4)') '# nlevels=', inp%nlevels, ' terminator=', inp%terminator, ' eta=', inp%eta
  do ie = 0, inp%npoints - 1
    e = inp%emin + (inp%emax - inp%emin) * ie / (inp%npoints - 1)
    z = cmplx(e, inp%eta, dp)
    g = cf_green(a, b, inp%nlevels, z, sqrt_tail(z, ainf, binf))
    write (u, '(2es16.8)') e, ldos_from_g(g)
  end do
  close (u)
  print '(a)', 'wrote '//outfile
end program asr_main
