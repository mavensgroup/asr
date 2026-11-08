!> Distributed-memory layer.  Built with -DUSE_MPI, it wraps the few MPI calls the code needs.  Built without it, every routine
!> is a serial stub: one rank, rank 0, sums do nothing.  The rest of the code calls these routines without #ifdef, so one
!> source tree gives both programs.
!>
!> The MPI build includes mpif.h and declares its own explicit interfaces for the routines it calls.  It needs no mpi.mod or
!> mpi_f08.mod, so the MPI library may come from another compiler version (a Fortran module file is tied to the compiler that
!> wrote it).  Link with the C-backed mpif.h library and libmpi (OpenMPI: -lmpi_mpifh -lmpi; MPICH and Intel MPI: -lmpifort
!> -lmpi).
!>
!> Use: par_init at program start, par_finalize at the end.  Ranks other than 0 send standard output to /dev/null (par_init),
!> and only rank 0 writes files (callers test par_root).  par_size and par_rank are 1 and 0 until par_init has run, so library
!> code and tests that never call par_init behave as the serial program.
module asr_mpi
  use asr_kinds, only: dp
#ifdef USE_MPI
  use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
#endif
  implicit none
  private
  public :: par_init, par_finalize, par_rank, par_size, par_root, par_sum_z, par_sum_i8, par_sum_i

  integer :: rnk = 0, nrk = 1
  logical :: started = .false.

#ifdef USE_MPI
  include 'mpif.h'

  interface
    subroutine MPI_Init_thread(required, provided, ierror)
      integer, intent(in) :: required
      integer, intent(out) :: provided, ierror
    end subroutine MPI_Init_thread
    subroutine MPI_Initialized(flag, ierror)
      logical, intent(out) :: flag
      integer, intent(out) :: ierror
    end subroutine MPI_Initialized
    subroutine MPI_Comm_rank(comm, rank, ierror)
      integer, intent(in) :: comm
      integer, intent(out) :: rank, ierror
    end subroutine MPI_Comm_rank
    subroutine MPI_Comm_size(comm, size, ierror)
      integer, intent(in) :: comm
      integer, intent(out) :: size, ierror
    end subroutine MPI_Comm_size
    subroutine MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm, ierror)
      type(*) :: sendbuf                                 !! MPI_IN_PLACE here
      type(*), dimension(*) :: recvbuf
      integer, intent(in) :: count, datatype, op, comm
      integer, intent(out) :: ierror
    end subroutine MPI_Allreduce
    subroutine MPI_Finalize(ierror)
      integer, intent(out) :: ierror
    end subroutine MPI_Finalize
    subroutine MPI_Abort(comm, errorcode, ierror)
      integer, intent(in) :: comm, errorcode
      integer, intent(out) :: ierror
    end subroutine MPI_Abort
  end interface

  interface                                              ! POSIX, to send a rank's standard output to /dev/null
    function c_creat(path, mode) result(fd) bind(C, name='creat')
      import :: c_char, c_int
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), value :: mode
      integer(c_int) :: fd
    end function c_creat
    function c_dup2(oldfd, newfd) result(fd) bind(C, name='dup2')
      import :: c_int
      integer(c_int), value :: oldfd, newfd
      integer(c_int) :: fd
    end function c_dup2
    function c_close(fd) result(ret) bind(C, name='close')
      import :: c_int
      integer(c_int), value :: fd
      integer(c_int) :: ret
    end function c_close
  end interface
#endif

contains

  !> start MPI (threads funnelled: MPI is called only outside OpenMP parallel regions); ranks above 0 lose standard output
  subroutine par_init()
#ifdef USE_MPI
    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
    integer :: ierr, prov
    integer(c_int) :: fd, rc
    logical :: flag
    call MPI_Initialized(flag, ierr)
    if (.not. flag) then
      call MPI_Init_thread(MPI_THREAD_FUNNELED, prov, ierr)
      if (ierr /= MPI_SUCCESS) error stop 'par_init: MPI_Init_thread failed'
      started = .true.
      if (prov < MPI_THREAD_FUNNELED) write (error_unit, '(a)') 'par_init: warning: MPI gives no thread support (FUNNELED asked)'
    end if
    call MPI_Comm_rank(MPI_COMM_WORLD, rnk, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nrk, ierr)
    if (rnk > 0) then
      flush (output_unit)
      fd = c_creat('/dev/null'//c_null_char, 0_c_int)
      if (fd >= 0) then
        rc = c_dup2(fd, 1_c_int); rc = c_close(fd)
      end if
    end if
#endif
  end subroutine par_init

  !> end MPI if par_init started it
  subroutine par_finalize()
#ifdef USE_MPI
    integer :: ierr
    if (started) call MPI_Finalize(ierr)
    started = .false.
#endif
    rnk = 0; nrk = 1
  end subroutine par_finalize

  integer function par_rank()
    par_rank = rnk
  end function par_rank

  integer function par_size()
    par_size = nrk
  end function par_size

  logical function par_root()
    par_root = rnk == 0
  end function par_root

  !> x(1:n) = sum over ranks of x(1:n), the same on every rank (complex double)
  subroutine par_sum_z(x, n)
    complex(dp), intent(inout) :: x(*)
    integer, intent(in) :: n
#ifdef USE_MPI
    integer :: ierr
    if (nrk == 1 .or. n <= 0) return
    call MPI_Allreduce(MPI_IN_PLACE, x, n, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) call par_fail('par_sum_z')
#endif
  end subroutine par_sum_z

  !> x(1:n) = sum over ranks of x(1:n), the same on every rank (64-bit integer)
  subroutine par_sum_i8(x, n)
    integer(8), intent(inout) :: x(*)
    integer, intent(in) :: n
#ifdef USE_MPI
    integer :: ierr
    if (nrk == 1 .or. n <= 0) return
    call MPI_Allreduce(MPI_IN_PLACE, x, n, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) call par_fail('par_sum_i8')
#endif
  end subroutine par_sum_i8

  !> x(1:n) = sum over ranks of x(1:n), the same on every rank (default integer)
  subroutine par_sum_i(x, n)
    integer, intent(inout) :: x(*)
    integer, intent(in) :: n
#ifdef USE_MPI
    integer :: ierr
    if (nrk == 1 .or. n <= 0) return
    call MPI_Allreduce(MPI_IN_PLACE, x, n, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) call par_fail('par_sum_i')
#endif
  end subroutine par_sum_i

#ifdef USE_MPI
  subroutine par_fail(where)
    use, intrinsic :: iso_fortran_env, only: error_unit
    character(len=*), intent(in) :: where
    integer :: ierr
    write (error_unit, '(a,i0)') where//': MPI call failed on rank ', rnk
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end subroutine par_fail
#endif
end module asr_mpi
