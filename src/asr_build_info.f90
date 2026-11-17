!> Build provenance generated at compile time by tools/gen_build_info.sh
module asr_build_info
  implicit none
  character(len=*), parameter :: git_rev = 'unknown'
  logical, parameter :: git_dirty = .false.
end module asr_build_info
