!> Small helpers shared by the TOML input readers: sub-tables, strict key checking.
module asr_tomlutil
  use tomlf, only: toml_table, toml_array, toml_key, get_value, len
  implicit none
  private
  public :: section, check_keys

contains

  !> Get a sub-table; an absent section is an empty one (all defaults).
  subroutine section(table, name, child)
    type(toml_table), intent(inout) :: table
    character(len=*), intent(in) :: name
    type(toml_table), pointer, intent(out) :: child
    nullify (child)
    call get_value(table, name, child, requested=.false.)
    if (.not. associated(child)) call get_value(table, name, child)
  end subroutine section

  !> Unknown keys are an error (a typo must never silently fall back to a default).
  subroutine check_keys(table, allowed, where, msg)
    type(toml_table), intent(inout) :: table
    character(len=*), intent(in) :: allowed(:), where
    character(len=:), allocatable, intent(inout) :: msg
    type(toml_key), allocatable :: keys(:)
    integer :: i, j
    logical :: known
    call table%get_keys(keys)
    do i = 1, size(keys)
      known = .false.
      do j = 1, size(allowed)
        if (trim(keys(i)%key) == trim(allowed(j))) known = .true.
      end do
      if (.not. known) then
        msg = 'unknown key "'//trim(keys(i)%key)//'" in '//where; return
      end if
    end do
  end subroutine check_keys
end module asr_tomlutil
