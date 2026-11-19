!> Exact reference for ordered cubic NN s-bands: local Green function from a
!> uniform k-grid, G(z) = <1/(z - eps(k))>.  Independent of the recursion.
module asr_kspace
  use asr_kinds, only: dp, pi
  implicit none
  private
  public :: kspace_green, kspace_green_b2, coordination

contains

  pure integer function coordination(kind)
    character(*), intent(in) :: kind
    select case (kind)
    case ('sc'); coordination = 6
    case ('bcc'); coordination = 8
    case default; coordination = 12
    end select
  end function coordination

  !> G(E_i + i eta) for hopping -t, nk^3 midpoint grid, OpenMP.
  function kspace_green(kind, t, egrid, eta, nk) result(g)
    character(*), intent(in) :: kind
    real(dp), intent(in) :: t, egrid(:), eta
    integer, intent(in) :: nk
    complex(dp) :: g(size(egrid))
    integer :: ix, iy, iz, ie
    real(dp) :: kx, ky, kz, eps
    complex(dp) :: z
    g = (0.0_dp, 0.0_dp)
    !$omp parallel do collapse(2) private(ix, iy, iz, kx, ky, kz, eps, ie, z) reduction(+:g)
    do iz = 0, nk - 1
      do iy = 0, nk - 1
        do ix = 0, nk - 1
          kx = 2 * pi * (ix + 0.5_dp) / nk; ky = 2 * pi * (iy + 0.5_dp) / nk
          kz = 2 * pi * (iz + 0.5_dp) / nk
          select case (kind)
          case ('sc'); eps = -2 * t * (cos(kx) + cos(ky) + cos(kz))
          case ('bcc'); eps = -8 * t * cos(kx / 2) * cos(ky / 2) * cos(kz / 2)
          case default
            eps = -4 * t * (cos(kx / 2) * cos(ky / 2) + cos(ky / 2) * cos(kz / 2) + cos(kz / 2) * cos(kx / 2))
          end select
          do ie = 1, size(egrid)
            z = cmplx(egrid(ie), eta, dp)
            g(ie) = g(ie) + 1.0_dp / (z - eps)
          end do
        end do
      end do
    end do
    !$omp end parallel do
    g = g / real(nk, dp)**3
  end function kspace_green

  !> B2 (CsCl) ordered compound on bcc: on-site +delta on the A site, -delta on B,
  !> NN hopping -t between sublattices.  Local G on the A site:
  !>   G_AA = < (z + delta) / (z^2 - delta^2 - eps_k^2) >,  eps_k = -8 t cx cy cz.
  function kspace_green_b2(t, delta, egrid, eta, nk) result(g)
    real(dp), intent(in) :: t, delta, egrid(:), eta
    integer, intent(in) :: nk
    complex(dp) :: g(size(egrid))
    integer :: ix, iy, iz, ie
    real(dp) :: kx, ky, kz, eps
    complex(dp) :: z
    g = (0.0_dp, 0.0_dp)
    !$omp parallel do collapse(2) private(ix, iy, iz, kx, ky, kz, eps, ie, z) reduction(+:g)
    do iz = 0, nk - 1
      do iy = 0, nk - 1
        do ix = 0, nk - 1
          kx = 2 * pi * (ix + 0.5_dp) / nk; ky = 2 * pi * (iy + 0.5_dp) / nk
          kz = 2 * pi * (iz + 0.5_dp) / nk
          eps = -8 * t * cos(kx / 2) * cos(ky / 2) * cos(kz / 2)
          do ie = 1, size(egrid)
            z = cmplx(egrid(ie), eta, dp)
            g(ie) = g(ie) + (z + delta) / (z * z - delta * delta - eps * eps)
          end do
        end do
      end do
    end do
    !$omp end parallel do
    g = g / real(nk, dp)**3
  end function kspace_green_b2
end module asr_kspace
