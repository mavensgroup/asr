!> Exact ASA potential function at complex energy (KKR-ASA, no linearisation).
!>
!> The regular scalar-relativistic radial solution (P, Q) of the sphere is integrated at complex energy z.  With
!> p = P(S), u = M Q(S) S (M = 1 + (z - V(S))/c^2) the logarithmic derivative is D = u/p and
!>      P^0(z) = K (D + l + 1)/(D - l),   K = 2(2l+1) (w/S)^(2l+1),        (same convention as asr_potpar)
!>      P^alpha = P^0 / (1 - alpha P^0) = K A / C,   A = u + (l+1) p,  C = (u - l p) - alpha K A.
!> A and C are entire functions of z; their first two z-derivatives come from the variational equations (exact), and
!>      mu^2 = dP^alpha/dz,   lambda = d ln C / dz = C'/C
!> (lambda cancels the poles of mu^2 g at the real zeros of C and is analytic elsewhere; the textbook -P''/(2P') agrees for
!> the two-parameter P but has spurious poles in the upper half plane at the complex zeros of P' of an exact, meromorphic P).  The physical (sphere-normalised)
!> site Green function is G = lambda + mu^2 g^alpha, and the valence density of channel l is
!>      rho_l(r) = -(1/pi) Im int G_l(z) phi_l(z, r)^2 dz,   int_0^S phi^2 dr = 1  (complex, no conjugation).
module asr_pfun
  use asr_kinds, only: dp
  use asr_radial, only: radmesh_t
  implicit none
  private
  public :: pfun_exact, regular_c


contains

  !> Regular solution y(1:2, k, r) = (P_k, Q_k), k = 0..2 the energy derivatives (variational equations, exact), at complex
  !> energy e, integrated outward to the last mesh point (RK4 in the mesh index, same start and midpoint interpolation as
  !> asr_radial::integrate_regular; the derivative components start at zero).
  subroutine regular_c(m, v, l, e, z, c, y)
    type(radmesh_t), intent(in) :: m
    real(dp), intent(in) :: v(:), z, c
    integer, intent(in) :: l
    complex(dp), intent(in) :: e
    complex(dp), intent(out) :: y(2, 0:2, m%n)
    real(dp), allocatable :: rv(:)
    real(dp) :: r0, rm, um, c1, sx
    complex(dp) :: m0, k1(2, 0:2), k2(2, 0:2), k3(2, 0:2), k4(2, 0:2), yy(2, 0:2)
    integer :: i, n
    n = m%n
    allocate (rv(n))
    rv(2:n) = m%r(2:n) * v(2:n)
    rv(1) = 3.0_dp * rv(2) - 3.0_dp * rv(3) + rv(4)
    y = (0.0_dp, 0.0_dp)
    r0 = m%r(2)
    m0 = 1.0_dp + (e - v(2)) / (c * c)
    if (real(m0) > 3.0_dp) then
      sx = sqrt(real(l * (l + 1), dp) + 1.0_dp - (2.0_dp * z / c)**2)
      y(1, 0, 2) = r0**sx
      y(2, 0, 2) = (sx - 1.0_dp) * c * c / (2.0_dp * z) * r0**sx
    else
      c1 = -z / real(l + 1, dp)
      y(1, 0, 2) = r0**(l + 1) * (1.0_dp + c1 * r0)
      y(2, 0, 2) = r0**l * (real(l, dp) + c1 * real(l + 1, dp) * r0) / m0
    end if
    do i = 2, n - 1
      yy = y(:, :, i)
      call rhs_c(l, e, c, m%r(i), rv(i), yy, k1); k1 = k1 * m%drdt(i)
      rm = m%b * (exp(m%a * (i - 0.5_dp)) - 1.0_dp)
      um = mid(rv, i, n)
      call rhs_c(l, e, c, rm, um, yy + 0.5_dp * k1, k2); k2 = k2 * (m%a * (rm + m%b))
      call rhs_c(l, e, c, rm, um, yy + 0.5_dp * k2, k3); k3 = k3 * (m%a * (rm + m%b))
      call rhs_c(l, e, c, m%r(i + 1), rv(i + 1), yy + k3, k4); k4 = k4 * m%drdt(i + 1)
      y(:, :, i + 1) = yy + (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4) / 6.0_dp
      if (abs(y(1, 0, i + 1)) > 1.0e100_dp) error stop 'regular_c: overflow (energy far outside the valence window)'
    end do
  end subroutine regular_c

  !> complex version of asr_radial::rhs with the energy-derivative (variational) equations up to second order
  pure subroutine rhs_c(l, e, c, r, rvv, y, dy)
    integer, intent(in) :: l
    complex(dp), intent(in) :: e, y(2, 0:2)
    real(dp), intent(in) :: c, r, rvv
    complex(dp), intent(out) :: dy(2, 0:2)
    complex(dp) :: mm, a21, a21d, a21dd
    real(dp) :: v, c2, ll
    v = rvv / r; c2 = c * c; ll = real(l * (l + 1), dp)
    mm = 1.0_dp + (e - v) / c2
    a21 = ll / (mm * r * r) + v - e
    a21d = -ll / (mm * mm * r * r * c2) - 1.0_dp
    a21dd = 2.0_dp * ll / (mm**3 * r * r * c2 * c2)
    dy(1, 0) = mm * y(2, 0) + y(1, 0) / r
    dy(2, 0) = -y(2, 0) / r + a21 * y(1, 0)
    dy(1, 1) = mm * y(2, 1) + y(1, 1) / r + y(2, 0) / c2
    dy(2, 1) = -y(2, 1) / r + a21 * y(1, 1) + a21d * y(1, 0)
    dy(1, 2) = mm * y(2, 2) + y(1, 2) / r + 2.0_dp * y(2, 1) / c2
    dy(2, 2) = -y(2, 2) / r + a21 * y(1, 2) + 2.0_dp * a21d * y(1, 1) + a21dd * y(1, 0)
  end subroutine rhs_c

  pure function mid(rv, i, n) result(u)
    real(dp), intent(in) :: rv(:)
    integer, intent(in) :: i, n
    real(dp) :: u
    if (i >= 2 .and. i + 2 <= n) then
      u = (-rv(i - 1) + 9.0_dp * rv(i) + 9.0_dp * rv(i + 1) - rv(i + 2)) / 16.0_dp
    else if (i + 2 <= n) then
      u = (5.0_dp * rv(i) + 15.0_dp * rv(i + 1) - 5.0_dp * rv(i + 2) + rv(i + 3)) / 16.0_dp
    else
      u = (rv(i - 2) - 5.0_dp * rv(i - 1) + 15.0_dp * rv(i) + 5.0_dp * rv(i + 1)) / 16.0_dp
    end if
  end function mid

  !> P^alpha(z), mu^2 = dP^alpha/dz, lambda = -P''/(2P') for channel l at complex energy z; optionally the sphere-normalised
  !> radial density function phi2(r) = (P^2 + Q^2/c^2)/N (4 pi r^2 |phi|^2 convention of asr_radial, complex).
  !> tsing (optional): the single-site term of the site Green function, so that the channel-l valence density is
  !>   rho_l(r) = -1/pi Im int [ sum_m G_m(z) phi2(r) + (2l+1) tsing(r) ] dz,
  !> tsing = f - (int f) phi2, f = (P_phi P_eta + Q_phi Q_eta/c^2)/W (G = (z - H)^-1) with the irregular solution eta of log-derivative
  !> D* = (l + alpha K (l+1))/(1 - alpha K) at S (where P^alpha has its poles): tsing carries no charge, is real on the real
  !> axis and has no poles (eta ~ phi wherever W = 0).
  subroutine pfun_exact(m, v, l, z, zat, c, w, alpha, pf, mu2, lam, phi2, tsing)
    type(radmesh_t), intent(in) :: m
    real(dp), intent(in) :: v(:), zat, c, w, alpha
    integer, intent(in) :: l
    complex(dp), intent(in) :: z
    complex(dp), intent(out) :: pf, mu2, lam
    complex(dp), intent(out), optional :: phi2(m%n), tsing(m%n)
    complex(dp), allocatable :: y(:, :, :)
    complex(dp) :: mm, u0, u1, u2, a0, a1, a2, c0, c1, c2, f, pd, nrm, wr
    complex(dp), allocatable :: ph(:), eta(:, :), ff(:)
    real(dp) :: s, kk
    integer :: n, i
    n = m%n; s = m%r(n)
    kk = 2.0_dp * (2 * l + 1) * (w / s)**(2 * l + 1)
    allocate (y(2, 0:2, n))
    call regular_c(m, v, l, z, zat, c, y)
    ! u = M Q(S) S with M = 1 + (z - V(S))/c^2 (dM/dz = 1/c^2), and its z-derivatives
    mm = 1.0_dp + (z - v(n)) / (c * c)
    u0 = mm * y(2, 0, n) * s
    u1 = (y(2, 0, n) / (c * c) + mm * y(2, 1, n)) * s
    u2 = (2.0_dp * y(2, 1, n) / (c * c) + mm * y(2, 2, n)) * s
    a0 = u0 + (l + 1) * y(1, 0, n); a1 = u1 + (l + 1) * y(1, 1, n); a2 = u2 + (l + 1) * y(1, 2, n)
    c0 = u0 - l * y(1, 0, n) - alpha * kk * a0
    c1 = u1 - l * y(1, 1, n) - alpha * kk * a1
    c2 = u2 - l * y(1, 2, n) - alpha * kk * a2
    pf = kk * a0 / c0
    f = a1 * c0 - a0 * c1                                        ! (A'C - AC')
    pd = kk * f / c0**2
    mu2 = pd
    lam = c1 / c0                                               ! d ln C/dz: residue +1 at the poles of P^alpha only
    if (present(phi2) .or. present(tsing)) then
      allocate (ph(n), eta(2, n), ff(n))
      do i = 1, n; ph(i) = y(1, 0, i)**2 + y(2, 0, i)**2 / (c * c); end do
      nrm = simpson_c(ph * m%drdt, n)
      ph = ph / nrm
      if (present(phi2)) phi2 = ph
      if (present(tsing)) then
        call irregular_c(m, v, l, z, c, (l + alpha * kk * (l + 1)) / (1.0_dp - alpha * kk), eta)
        wr = y(1, 0, n) * eta(2, n) - y(2, 0, n) * eta(1, n)
        ff(1) = 0.0_dp
        do i = 2, n; ff(i) = (y(1, 0, i) * eta(1, i) + y(2, 0, i) * eta(2, i) / (c * c)) / wr; end do     ! (z - H)^-1 convention
        tsing = ff - simpson_c(ff * m%drdt, n) * ph
      end if
    end if
  end subroutine pfun_exact

  !> Irregular solution (P, Q) at complex energy e with logarithmic derivative dstar at S (P(S) = 1), integrated inward
  !> (stable: it grows towards the origin); y(:, 1) is left at zero (the origin is never used).
  subroutine irregular_c(m, v, l, e, c, dstar, y)
    type(radmesh_t), intent(in) :: m
    real(dp), intent(in) :: v(:), c, dstar
    integer, intent(in) :: l
    complex(dp), intent(in) :: e
    complex(dp), intent(out) :: y(2, m%n)
    real(dp), allocatable :: rv(:)
    real(dp) :: rm, um, s
    complex(dp) :: k1(2), k2(2), k3(2), k4(2), yy(2)
    integer :: i, n
    n = m%n; s = m%r(n)
    allocate (rv(n))
    rv(2:n) = m%r(2:n) * v(2:n)
    rv(1) = 3.0_dp * rv(2) - 3.0_dp * rv(3) + rv(4)
    y = (0.0_dp, 0.0_dp)
    y(1, n) = 1.0_dp
    y(2, n) = dstar / ((1.0_dp + (e - v(n)) / (c * c)) * s)       ! D = M Q S / P
    do i = n, 3, -1
      yy = y(:, i)
      call rhs0(l, e, c, m%r(i), rv(i), yy, k1); k1 = -k1 * m%drdt(i)
      rm = m%b * (exp(m%a * (i - 1.5_dp)) - 1.0_dp)
      um = mid(rv, i - 1, n)
      call rhs0(l, e, c, rm, um, yy + 0.5_dp * k1, k2); k2 = -k2 * (m%a * (rm + m%b))
      call rhs0(l, e, c, rm, um, yy + 0.5_dp * k2, k3); k3 = -k3 * (m%a * (rm + m%b))
      call rhs0(l, e, c, m%r(i - 1), rv(i - 1), yy + k3, k4); k4 = -k4 * m%drdt(i - 1)
      y(:, i - 1) = yy + (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4) / 6.0_dp
    end do
  end subroutine irregular_c

  pure subroutine rhs0(l, e, c, r, rvv, y, dy)
    integer, intent(in) :: l
    complex(dp), intent(in) :: e, y(2)
    real(dp), intent(in) :: c, r, rvv
    complex(dp), intent(out) :: dy(2)
    complex(dp) :: mm
    real(dp) :: v
    v = rvv / r
    mm = 1.0_dp + (e - v) / (c * c)
    dy(1) = mm * y(2) + y(1) / r
    dy(2) = -y(2) / r + (real(l * (l + 1), dp) / (mm * r * r) + v - e) * y(1)
  end subroutine rhs0


  pure function simpson_c(f, n) result(sm)
    complex(dp), intent(in) :: f(:)
    integer, intent(in) :: n
    complex(dp) :: sm
    integer :: i
    sm = f(1) + f(n)
    do i = 2, n - 1
      if (mod(i, 2) == 0) then; sm = sm + 4.0_dp * f(i); else; sm = sm + 2.0_dp * f(i); end if
    end do
    sm = sm / 3.0_dp
  end function simpson_c
end module asr_pfun
