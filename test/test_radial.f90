!> Radial solver: hydrogenic limits, first-order relativistic shift, energy derivatives.
program test_radial
  use asr_kinds, only: dp
  use asr_radial
  use asr_testing
  implicit none
  type(radmesh_t) :: m
  real(dp), allocatable :: v(:), y(:, :, :), yp(:, :, :), ym(:, :, :), yb(:, :, :)
  real(dp) :: e, z, ex, alpha2, shift_ex, shift, de
  integer :: l, n, nod, isph, i
  integer, parameter :: nz = 2
  real(dp), parameter :: zs(nz) = [1.0_dp, 10.0_dp]

  print '(a)', 'test_radial'
  ! ---- non-relativistic hydrogenic levels E = -Z^2/n^2 (Ry)
  do i = 1, nz
    z = zs(i)
    call make_radmesh(z, 60.0_dp / z, 0.01_dp, m, c_nonrel)
    isph = m%n; if (mod(isph, 2) == 0) isph = isph - 1
    allocate (v(m%n)); v(2:) = -2.0_dp * z / m%r(2:); v(1) = v(2)
    if (allocated(yb)) deallocate (yb)
    allocate (yb(2, 0:0, isph))
    do n = 1, 3
      do l = 0, n - 1
        nod = n - l - 1
        call bound_state(m, v, l, nod, z, c_nonrel, isph, -1.5_dp * z * z, -0.01_dp * z * z, e, yb)
        call check_close(e, -z * z / real(n * n, dp), 2.0e-7_dp * z * z, name('nonrel Z=', z, n, l))
      end do
    end do
    deallocate (v)
  end do

  ! ---- first-order scalar-relativistic shift (mass-velocity + Darwin), Z = 3
  z = 3.0_dp
  alpha2 = (2.0_dp / c_ry)**2                      ! alpha^2 in these units
  call make_radmesh(z, 60.0_dp / z, 0.01_dp, m)
  isph = m%n; if (mod(isph, 2) == 0) isph = isph - 1
  allocate (v(m%n)); v(2:) = -2.0_dp * z / m%r(2:); v(1) = v(2)
  deallocate (yb); allocate (yb(2, 0:0, isph))
  do n = 1, 2
    do l = 0, n - 1
      nod = n - l - 1
      call bound_state(m, v, l, nod, z, c_ry, isph, -1.5_dp * z * z, -0.01_dp * z * z, e, yb)
      ex = -z * z / real(n * n, dp)
      if (l == 0) then
        shift_ex = ex * (z * z * alpha2 / real(n * n, dp)) * (real(n, dp) - 0.75_dp)
      else
        shift_ex = ex * (z * z * alpha2 / real(n * n, dp)) * (real(n, dp) / (l + 0.5_dp) - 0.75_dp)
      end if
      shift = e - ex
      call check_close(shift / shift_ex, 1.0_dp, 0.02_dp, name('rel shift ratio n=', real(n, dp), n, l))
    end do
  end do

  ! ---- energy derivatives vs finite differences (l=1, Z=10, E=-7 Ry, sphere r=3)
  z = 10.0_dp
  call make_radmesh(z, 3.0_dp, 0.01_dp, m)
  isph = m%n; if (mod(isph, 2) == 0) isph = isph - 1
  deallocate (v); allocate (v(m%n)); v(2:) = -2.0_dp * z / m%r(2:) + 0.3_dp; v(1) = v(2)
  allocate (y(2, 0:2, isph), yp(2, 0:2, isph), ym(2, 0:2, isph))
  e = -7.0_dp; de = 1.0e-3_dp
  call integrate_regular(m, v, 1, e, z, c_ry, 2, isph, y)
  call integrate_regular(m, v, 1, e + de, z, c_ry, 2, isph, yp)
  call integrate_regular(m, v, 1, e - de, z, c_ry, 2, isph, ym)
  call check_close(maxval(abs((yp(:, 0, isph) - ym(:, 0, isph)) / (2 * de) - y(:, 1, isph)) / abs(y(:, 1, isph))), 0.0_dp, 1.0e-5_dp, &
                   'dP/dE (variational) vs central difference')
  call check_close(maxval(abs((yp(:, 0, isph) - 2 * y(:, 0, isph) + ym(:, 0, isph)) / de**2 - y(:, 2, isph)) / abs(y(:, 2, isph))), &
                   0.0_dp, 1.0e-4_dp, 'd2P/dE2 (variational) vs second difference')
  call check_close(maxval(abs((yp(:, 1, isph) - ym(:, 1, isph)) / (2 * de) - y(:, 2, isph)) / abs(y(:, 2, isph))), &
                   0.0_dp, 1.0e-5_dp, 'dPdot/dE vs d2P/dE2 consistency')
  call report()
contains
  function name(pre, zz, nn, ll) result(s)
    character(*), intent(in) :: pre
    real(dp), intent(in) :: zz
    integer, intent(in) :: nn, ll
    character(len=48) :: s
    write (s, '(a,f5.1,a,i0,a,i0)') pre, zz, ' n=', nn, ' l=', ll
  end function name
end program test_radial
