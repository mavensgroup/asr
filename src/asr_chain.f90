!> The chain machinery of the energy-dependent recursion (asr_edc, asr_ed3): the terminated continued fraction of a
!> complex-symmetric Lanczos chain (cf_species, with its square-root tails), the adaptive depth, the projection switch of the
!> group-reduced chains, and the chain statistics.  Split from asr_endep, which keeps the binary reference routes.
module asr_chain
  use asr_lapack, only: zgeev
  use asr_kinds, only: dp
  implicit none
  private
  public :: ed_depth_tol, ed_nchain, ed_nlev, ed_nfall, ed_tail, cf_species, ed_project

  !> adaptive depth: from level 8, every 4 levels, stop when the terminated Green function changes by less than ed_depth_tol
  !> (relative); 0 = always use the full depth.  ed_nchain / ed_nlev count chains and levels used (statistics).
  real(dp) :: ed_depth_tol = 0.0_dp
  integer(8) :: ed_nchain = 0, ed_nlev = 0
  !> 'bp' tails whose fixed point failed (zgeev error: last iterate kept; 200 iterations: hull values): statistics
  integer(8) :: ed_nfall = 0
  !> square-root tail of cf_species: 'hull' (band edges = extreme eigenvalues of the finite chain), 'bp' (Beer-Pettifor fixed
  !> point: edges where the terminated fraction has its poles, outside the hull), 'avg2' (last two coefficients).
  !> Warning: on the complex chains of the energy-dependent route 'bp' can lock onto a state split off by the end-site shift; the
  !> terminated fraction then has a pole next to the contour and N(E_F) varies steeply (levels-4 Fe66Ni14Cr20 SCF stall, Cr e_g).
  !> 'hull' is the default of the class SCF; test_terminators checks the Nevanlinna bound |dG/dz| <= -Im G / Im z.
  character(len=8) :: ed_tail = 'hull'
  !> The chains on a group-reduced operator (asr_edc, asr_ed3) project every new Lanczos vector onto the symmetric subspace
  !> (asr_augcls::augcls_project_c).  .false. switches this off: for tests and diagnostics only.
  logical :: ed_project = .true.

contains

  !> g_11 and sum_k c_k [g]_{k1} of the complex-symmetric tridiagonal chain (a, b) with spectral variable z, terminated by a
  !> square-root tail whose band edges are the extreme eigenvalues (real parts) of the finite chain (Beer-Pettifor type): the
  !> tail band then contains the spectral variable wherever the energy lies inside the band, giving a continuum instead of
  !> m isolated poles.  tail_in overrides ed_tail (thread-safe tail scans); edges returns the tail's (centre, quarter width).
  subroutine cf_species(a, b, c, m, z, g11, f, tail_in, edges)
    complex(dp), intent(in) :: a(:), b(:), c(:), z
    integer, intent(in) :: m
    complex(dp), intent(out) :: g11, f
    character(len=*), intent(in), optional :: tail_in
    complex(dp), intent(out), optional :: edges(2)
    complex(dp) :: gt(m + 1), prod, tail, bb, d, s, g1, g2, aa
    complex(dp) :: tm(m, m), ev(m), vl(1, 1), vr(1, 1), wrk(8 * m + 8)
    real(dp) :: rw(2 * m)
    integer :: j, k, info, ilo, ihi
    character(len=8) :: tl
    tl = ed_tail; if (present(tail_in)) tl = tail_in
    if (m >= 2 .and. trim(tl) /= 'none') then
      tm = (0.0_dp, 0.0_dp)
      do j = 1, m
        tm(j, j) = a(j)
        if (j < m) then; tm(j, j + 1) = b(j); tm(j + 1, j) = b(j); end if
      end do
      call zgeev('N', 'N', m, tm, m, ev, vl, 1, vr, 1, wrk, size(wrk), rw, info)
    else
      info = 1
    end if
    if (info == 0) then
      ilo = minloc(real(ev), 1); ihi = maxloc(real(ev), 1)
      aa = 0.5_dp * (ev(ilo) + ev(ihi)); bb = 0.25_dp * (ev(ihi) - ev(ilo))
      if (trim(tl) == 'bp') call bp_fixed_point()
    else
      aa = a(m); bb = (0.0_dp, 0.0_dp); if (m > 1) bb = b(m - 1)
    end if
    if (trim(tl) == 'avg2' .and. m >= 3) then
      aa = 0.5_dp * (a(m) + a(m - 1)); bb = sqrt(0.5_dp * (b(m - 1)**2 + b(m - 2)**2))
    end if
    if (trim(tl) == 'none') then; aa = (0.0_dp, 0.0_dp); bb = (0.0_dp, 0.0_dp); end if   ! exact finite spectra (unit tests)
    if (present(edges)) edges = [aa, bb]
    d = z - aa
    if (abs(bb) < 1.0e-14_dp) then
      tail = (0.0_dp, 0.0_dp)
    else
      s = sqrt(d * d - 4.0_dp * bb * bb)
      g1 = (d - s) / (2.0_dp * bb * bb); g2 = (d + s) / (2.0_dp * bb * bb)
      if (abs(bb * g1) <= abs(bb * g2)) then; tail = bb * bb * g1; else; tail = bb * bb * g2; end if
    end if
    gt(m + 1) = (0.0_dp, 0.0_dp)
    gt(m) = 1.0_dp / (z - a(m) - tail)
    do j = m - 1, 1, -1
      gt(j) = 1.0_dp / (z - a(j) - b(j)**2 * gt(j + 1))
    end do
    g11 = gt(1)
    f = (0.0_dp, 0.0_dp); prod = gt(1)
    do k = 1, m
      f = f + c(k) * prod
      if (k < m) prod = prod * b(k) * gt(k + 1)
    end do
  contains
    !> Beer-Pettifor: A + 2B (A - 2B) = extreme eigenvalue of the chain with a_m -> a_m + B (a_m - B), i.e. the terminated
    !> fraction has its poles exactly at the tail band edges (complex-symmetric version of asr_terminator%bp_asymptote);
    !> keeps the hull values if the fixed point does not converge
    subroutine bp_fixed_point()
      complex(dp) :: ep, em, an, bn
      integer :: it
      do it = 1, 200
        call chain_matrix(a(m) + bb); call zgeev('N', 'N', m, tm, m, ev, vl, 1, vr, 1, wrk, size(wrk), rw, info)
        if (info /= 0) then
          !$omp atomic
          ed_nfall = ed_nfall + 1
          return
        end if
        ep = ev(maxloc(real(ev), 1))
        call chain_matrix(a(m) - bb); call zgeev('N', 'N', m, tm, m, ev, vl, 1, vr, 1, wrk, size(wrk), rw, info)
        if (info /= 0) then
          !$omp atomic
          ed_nfall = ed_nfall + 1
          return
        end if
        em = ev(minloc(real(ev), 1))
        an = 0.5_dp * (ep + em); bn = 0.25_dp * (ep - em)
        if (abs(an - aa) + abs(bn - bb) < 1.0e-12_dp * (1.0_dp + abs(bn))) then
          aa = an; bb = bn; return
        end if
        aa = an; bb = bn
      end do
      ! not converged: back to the hull
      !$omp atomic
      ed_nfall = ed_nfall + 1
      call chain_matrix(a(m)); call zgeev('N', 'N', m, tm, m, ev, vl, 1, vr, 1, wrk, size(wrk), rw, info)
      ilo = minloc(real(ev), 1); ihi = maxloc(real(ev), 1)
      aa = 0.5_dp * (ev(ilo) + ev(ihi)); bb = 0.25_dp * (ev(ihi) - ev(ilo))
    end subroutine bp_fixed_point

    subroutine chain_matrix(am)
      complex(dp), intent(in) :: am
      tm = (0.0_dp, 0.0_dp)
      do j = 1, m
        tm(j, j) = a(j)
        if (j < m) then; tm(j, j + 1) = b(j); tm(j + 1, j) = b(j); end if
      end do
      tm(m, m) = am
    end subroutine chain_matrix
  end subroutine cf_species

end module asr_chain
