module asr_kinds
  use, intrinsic :: iso_fortran_env, only: real64, int32, int64
  implicit none
  private
  integer, parameter, public :: dp = real64
  integer, parameter, public :: i4 = int32
  integer, parameter, public :: i8 = int64
  real(dp), parameter, public :: pi = 3.14159265358979323846264338327950288_dp
  complex(dp), parameter, public :: ci = (0.0_dp, 1.0_dp)
end module asr_kinds
