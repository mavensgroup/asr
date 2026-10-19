!> asr <input.toml>
!> Reads one validated TOML record (any directory) and runs the calculation it describes (asr_driver: method = "asr" or
!> "kspace").  Output files are written next to the input file.
program asr_main
  use asr_driver, only: run_input_t, read_run_input, run_input
  implicit none
  type(run_input_t) :: inp
  logical :: ok
  character(len=:), allocatable :: msg
  character(len=512) :: arg

  if (command_argument_count() /= 1) then
    print '(a)', 'usage: asr <input.toml>'; stop 2
  end if
  call get_command_argument(1, arg)
  call read_run_input(trim(arg), inp, msg)
  if (len(msg) > 0) then; print '(a)', 'asr: input error: '//msg; stop 1; end if
  call run_input(inp, ok)
  if (.not. ok) stop 3
end program asr_main
