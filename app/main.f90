!> asr <input.toml>
!> Reads one validated TOML record (any directory) and runs the calculation it describes (asr_driver: method = "asr" or
!> "kspace").  Output files are written next to the input file.  An MPI build (-DUSE_MPI, asr_mpi) runs the same program on
!> every rank; only rank 0 prints and writes files.
program asr_main
  use asr_driver, only: run_input_t, read_run_input, run_input
  use asr_mpi, only: par_init, par_finalize, par_size
  implicit none
  type(run_input_t) :: inp
  logical :: ok
  character(len=:), allocatable :: msg
  character(len=512) :: arg
  integer :: code

  call par_init()
  code = 0
  if (command_argument_count() /= 1) then
    print '(a)', 'usage: asr <input.toml>'; code = 2
  else
    if (par_size() > 1) print '(a,i0)', 'asr: MPI ranks ', par_size()
    call get_command_argument(1, arg)
    call read_run_input(trim(arg), inp, msg)
    if (len(msg) > 0) then
      print '(a)', 'asr: input error: '//msg; code = 1
    else
      call run_input(inp, ok)
      if (.not. ok) code = 3
    end if
  end if
  call par_finalize()
  if (code /= 0) stop code
end program asr_main
