!> Ewald lattice sums for point charges with a neutralising background.
!>   Phi(i,j) = electrostatic potential (per unit charge, units 1/Bohr) at site i due to a unit charge
!>   at site j and all its lattice images, excluding the self term.
!> For a neutral set of charges q_j the potential at i is  V_i = sum_j Phi(i,j) q_j  (the background
!> constant cancels), in Hartree atomic units (multiply by 2 for Rydberg).
module asr_ewald
  use asr_kinds, only: dp, pi
  implicit none
  private
  public :: madelung_phi

contains

  !> avec (columns) and tau in units of alat; alat in Bohr.
  subroutine madelung_phi(avec, tau, alat, phi)
    real(dp), intent(in) :: avec(3, 3), tau(:, :), alat
    real(dp), intent(out) :: phi(size(tau, 2), size(tau, 2))
    real(dp) :: a(3, 3), b(3, 3), omega, eta, rmax, gmax, d(3), r, g2, gv(3), sdir, srec
    integer :: nat, i, j, n1, n2, n3, m1, m2, m3, nr(3), ng(3)
    nat = size(tau, 2)
    a = alat * avec
    omega = abs(a(1, 1) * (a(2, 2) * a(3, 3) - a(3, 2) * a(2, 3)) - a(1, 2) * (a(2, 1) * a(3, 3) - a(3, 1) * a(2, 3)) &
                + a(1, 3) * (a(2, 1) * a(3, 2) - a(3, 1) * a(2, 2)))
    b(:, 1) = cross(a(:, 2), a(:, 3)); b(:, 2) = cross(a(:, 3), a(:, 1)); b(:, 3) = cross(a(:, 1), a(:, 2))
    b = 2.0_dp * pi * b / omega
    eta = sqrt(pi) / omega**(1.0_dp / 3.0_dp)
    rmax = 7.0_dp / eta; gmax = 14.0_dp * eta
    do i = 1, 3
      nr(i) = ceiling(rmax / minval([abs(dot_product(a(:, i), cross(a(:, mod(i, 3) + 1), a(:, mod(i + 1, 3) + 1)))) &
                                     / norm2(cross(a(:, mod(i, 3) + 1), a(:, mod(i + 1, 3) + 1)))])) + 1
      ng(i) = ceiling(gmax * norm2(a(:, i)) / (2.0_dp * pi)) + 1
    end do
    do i = 1, nat
      do j = 1, nat
        d = alat * (tau(:, i) - tau(:, j))
        sdir = 0.0_dp
        do n1 = -nr(1), nr(1); do n2 = -nr(2), nr(2); do n3 = -nr(3), nr(3)
              r = norm2(d + n1 * a(:, 1) + n2 * a(:, 2) + n3 * a(:, 3))
              if (r < 1.0e-10_dp) cycle
              if (r > rmax) cycle
              sdir = sdir + erfc(eta * r) / r
        end do; end do; end do
        srec = 0.0_dp
        do m1 = -ng(1), ng(1); do m2 = -ng(2), ng(2); do m3 = -ng(3), ng(3)
              if (m1 == 0 .and. m2 == 0 .and. m3 == 0) cycle
              gv = m1 * b(:, 1) + m2 * b(:, 2) + m3 * b(:, 3)
              g2 = dot_product(gv, gv)
              if (g2 > gmax**2) cycle
              srec = srec + exp(-g2 / (4.0_dp * eta * eta)) / g2 * cos(dot_product(gv, d))
        end do; end do; end do
        phi(i, j) = sdir + 4.0_dp * pi / omega * srec - pi / (omega * eta * eta)
        if (i == j) phi(i, j) = phi(i, j) - 2.0_dp * eta / sqrt(pi)
      end do
    end do
  end subroutine madelung_phi

  pure function cross(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
  end function cross
end module asr_ewald
