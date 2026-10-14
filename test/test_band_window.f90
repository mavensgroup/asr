!> band_window: the window [D_l = l, phi_l(S) = 0] of the valence branch contains E_nu-like energies and the
!> potential parameter C, and the free-atom 3d level lies inside the d window.
program test_band_window
  use asr_kinds, only: dp
  use asr_testing
  use asr_radial
  use asr_atom
  use asr_potpar
  implicit none
  type(atom_result_t) :: at
  type(radmesh_t) :: m
  type(potpar_l_t) :: pp
  real(dp), allocatable :: v(:)
  real(dp) :: eb, et, s, e
  integer :: k, l, nodes(0:2) = [3, 2, 0]
  print '(a)', 'test_band_window'
  ! free Fe atom potential (scalar relativistic), cut at a WS sphere S = 2.667 Bohr
  call solve_atom(26.0_dp, fe_shells(), 'vbh', 30.0_dp, 0.01_dp, .true., at)
  s = 2.6671_dp
  call make_radmesh(26.0_dp, s, 0.012_dp, m, c_ry)
  allocate (v(m%n))
  do k = 2, m%n
    v(k) = lin(at%mesh%r, at%v * at%mesh%r, m%r(k)) / m%r(k)
  end do
  v(1) = v(2)
  do l = 0, 2
    call band_window(m, v, l, nodes(l), 26.0_dp, c_ry, eb, et)
    e = 0.5_dp * (eb + et)
    call potpar_at_enu(m, v, l, e, 26.0_dp, c_ry, s, pp)
    print '(a,i2,a,3f10.5)', '    l', l, '  bottom, C(E_nu = mid), top', eb, pp%c, et
    call check(eb < et, 'window ordered')
    call check(pp%c > eb .and. pp%c < et, 'C inside the window')
  end do
  call band_window(m, v, 2, 0, 26.0_dp, c_ry, eb, et)
  print '(a,f10.5)', '    free-atom 3d level', at%eps(6)
  call check(at%eps(6) > eb .and. at%eps(6) < et, 'free-atom 3d level inside the d window')
  e = 5.0_dp
  call check(clamp_enu(m, v, 2, 0, 26.0_dp, c_ry, e) .and. e < et, 'clamp_enu pulls an E_nu above the d top back into the window')
  call report()
contains
  real(dp) function lin(r, f, r0)
    real(dp), intent(in) :: r(:), f(:), r0
    integer :: j
    j = max(1, min(size(r) - 1, count(r <= r0)))
    lin = f(j) + (f(j + 1) - f(j)) * (r0 - r(j)) / (r(j + 1) - r(j))
  end function lin

  function fe_shells() result(sh)
    type(shell_t), allocatable :: sh(:)
    sh = [shell_t(1, 0, 2.0_dp), shell_t(2, 0, 2.0_dp), shell_t(2, 1, 6.0_dp), shell_t(3, 0, 2.0_dp), shell_t(3, 1, 6.0_dp), &
          shell_t(3, 2, 6.0_dp), shell_t(4, 0, 2.0_dp)]
  end function fe_shells
end program test_band_window
