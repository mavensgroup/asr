!> The resolved input of a run: every key the reader used, with its value, in TOML form.  The reader (asr_driver) fills every
!> key it reads with its default into the parsed table (toml-f get_value with a default adds the key), so this table holds the
!> complete input after reading.  Keys that are not in the input file are marked '# default' (found by a second parse of the
!> file).  Sections that the chosen method does not read are left out.  Reals are written in the shortest form that reads back
!> to the same double; values that mean "automatic" (huge) are written as inf.
module asr_inputdump
  use asr_kinds, only: dp
  use tomlf, only: toml_table, toml_array, toml_keyval, toml_value, toml_key, toml_load, toml_error, &
                   len, is_array_of_tables, get_value
  use tomlf_constants, only: toml_type
  implicit none
  private
  public :: resolved_input

  integer, parameter :: mark_col = 46      !! column of the '# default' mark

contains

  !> Lines of the resolved input (separated by new_line('a')).  method: 'asr' or 'kspace'.
  subroutine resolved_input(filled, filename, method, text)
    type(toml_table), intent(inout) :: filled
    character(len=*), intent(in) :: filename, method
    character(len=:), allocatable, intent(out) :: text
    type(toml_table), allocatable :: orig
    type(toml_error), allocatable :: err
    type(toml_key), allocatable :: keys(:)
    class(toml_value), pointer :: v, vo
    type(toml_table), pointer :: to
    integer :: i, j
    text = ''
    call toml_load(orig, filename, error=err)
    if (allocated(err)) allocate (orig)
    call filled%get_keys(keys)
    do i = 1, size(keys)
      if (skip_section(keys(i)%key, method)) cycle
      call filled%get(keys(i)%key, v)
      nullify (vo); if (orig%has_key(keys(i)%key)) call orig%get(keys(i)%key, vo)
      select type (v)
      type is (toml_table)
        nullify (to)
        if (associated(vo)) then
          select type (vo); type is (toml_table); to => vo; end select
        end if
        if (len(text) > 0) call add('')
        call add('[' // keys(i)%key // ']')
        call table_body(v, to, keys(i)%key, method)
      type is (toml_array)
        if (is_array_of_tables(v)) then
          do j = 1, len(v)
            if (len(text) > 0) call add('')
            call add('[[' // keys(i)%key // ']]')
            call elem_table(v, vo, j, keys(i)%key)
          end do
        else
          call add(entry(keys(i)%key, value_string(v), .not. associated(vo)))
        end if
      type is (toml_keyval)
        call add(entry(keys(i)%key, value_string(v), .not. associated(vo)))
      end select
    end do
  contains
    subroutine add(line)
      character(len=*), intent(in) :: line
      if (len(text) > 0) then
        text = text // new_line('a') // line
      else
        text = line
      end if
    end subroutine add

    !> one element of an array of tables (a [[species]], [[sublattice]] or [[sro]] entry)
    subroutine elem_table(arr, arro, j, name)
      type(toml_array), intent(inout) :: arr
      class(toml_value), pointer, intent(in) :: arro
      integer, intent(in) :: j
      character(len=*), intent(in) :: name
      class(toml_value), pointer :: e, eo
      type(toml_table), pointer :: to_
      nullify (to_, eo)
      call arr%get(j, e)
      if (associated(arro)) then
        select type (arro)
        type is (toml_array)
          if (j <= len(arro)) call arro%get(j, eo)
        end select
      end if
      if (associated(eo)) then
        select type (eo); type is (toml_table); to_ => eo; end select
      end if
      select type (e)
      type is (toml_table)
        call table_body(e, to_, name, method)
      end select
    end subroutine elem_table

    !> key = value lines of one table; t0: the same table in the input file (null when the whole table is a default)
    recursive subroutine table_body(t, t0, name, meth)
      type(toml_table), intent(inout) :: t
      type(toml_table), pointer, intent(in) :: t0
      character(len=*), intent(in) :: name, meth
      type(toml_key), allocatable :: ks(:)
      class(toml_value), pointer :: w
      integer :: k
      logical :: given
      call t%get_keys(ks)
      do k = 1, size(ks)
        if (skip_key(name, ks(k)%key, meth)) cycle
        given = .false.
        if (associated(t0)) given = t0%has_key(ks(k)%key)
        call t%get(ks(k)%key, w)
        select type (w)
        type is (toml_keyval)
          call add(entry(ks(k)%key, value_string(w), .not. given))
        type is (toml_array)
          call add(entry(ks(k)%key, value_string(w), .not. given))
        type is (toml_table)
          call add(entry(ks(k)%key, '{ ... }', .not. given))     ! no nested tables in the ASR input
        end select
      end do
    end subroutine table_body
  end subroutine resolved_input

  !> sections the chosen method does not read
  logical function skip_section(name, method)
    character(len=*), intent(in) :: name, method
    skip_section = (method == 'asr' .and. name == 'kspace') .or. (method == 'kspace' .and. name == 'asr')
  end function skip_section

  !> keys the chosen method does not read
  logical function skip_key(section, key, method)
    character(len=*), intent(in) :: section, key, method
    skip_key = method == 'asr' .and. section == 'output' .and. key == 'params_file'
  end function skip_key

  function entry(key, val, dflt) result(line)
    character(len=*), intent(in) :: key, val
    logical, intent(in) :: dflt
    character(len=:), allocatable :: line
    line = key // ' = ' // val
    if (dflt) line = line // repeat(' ', max(2, mark_col - len(line))) // '# default'
  end function entry

  !> TOML text of a value (scalar or array, arrays recursively)
  recursive function value_string(v) result(s)
    class(toml_value), intent(inout) :: v
    character(len=:), allocatable :: s
    class(toml_value), pointer :: e
    character(len=:), allocatable :: sv
    integer :: i
    integer(8) :: iv
    real(dp) :: rv
    logical :: lv
    s = ''
    select type (v)
    type is (toml_keyval)
      select case (v%get_type())
      case (toml_type%string)
        call get_value(v, sv); s = '"' // sv // '"'
      case (toml_type%int)
        call get_value(v, iv); s = itoa8(iv)
      case (toml_type%float)
        call get_value(v, rv); s = real_string(rv)
      case (toml_type%boolean)
        call get_value(v, lv); s = merge('true ', 'false', lv); s = trim(s)
      case default
        s = '?'
      end select
    type is (toml_array)
      s = '['
      do i = 1, len(v)
        call v%get(i, e)
        if (i > 1) s = s // ', '
        s = s // value_string(e)
      end do
      s = s // ']'
    type is (toml_table)
      s = '{ ... }'
    end select
  end function value_string

  !> shortest decimal form that reads back to x (fixed for 1e-4 <= |x| < 1e6, else exponent); |x| >= 1e300: inf
  function real_string(x) result(s)
    real(dp), intent(in) :: x
    character(len=:), allocatable :: s
    character(len=40) :: b
    character(len=16) :: f
    real(dp) :: y
    integer :: d, ios
    if (abs(x) >= 1.0e300_dp) then
      s = merge('inf ', '-inf', x > 0.0_dp); s = trim(s); return
    end if
    if (x == 0.0_dp) then
      s = '0.0'; return
    end if
    if (abs(x) >= 1.0e-4_dp .and. abs(x) < 1.0e6_dp) then
      do d = 1, 17
        write (f, '(a,i0,a)') '(f40.', d, ')'
        write (b, f) x
        read (b, *, iostat=ios) y
        if (ios == 0 .and. y == x) exit
      end do
    else
      do d = 1, 17
        if (abs(log10(abs(x))) < 99.0_dp) then
          write (f, '(a,i0,a,i0,a)') '(es', d + 8, '.', d, 'e2)'
        else
          write (f, '(a,i0,a,i0,a)') '(es', d + 9, '.', d, 'e3)'
        end if
        write (b, f) x
        read (b, *, iostat=ios) y
        if (ios == 0 .and. y == x) exit
      end do
    end if
    s = trim(adjustl(b))
  end function real_string

  function itoa8(i) result(s)
    integer(8), intent(in) :: i
    character(len=:), allocatable :: s
    character(len=24) :: b
    write (b, '(i0)') i
    s = trim(b)
  end function itoa8

end module asr_inputdump
