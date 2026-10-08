!> Linear tetrahedron method (Bloechl, Jepsen, Andersen 1994) on a regular k-grid in fractional
!> coordinates (works for any lattice).  Each grid cube is split into 6 tetrahedra along the main
!> diagonal; the grid is periodic.  Weights are per (band, k) and include the 1/Nk volume, so
!> sum(w) = number of occupied states per spin per cell (fraction of occupied band-volume).
module asr_tetra
  use asr_kinds, only: dp
  implicit none
  private
  public :: tetra_weights, tetra_fermi

  integer, parameter :: ctet(4, 6) = reshape([0, 1, 3, 7,  0, 1, 5, 7,  0, 2, 3, 7,  0, 2, 6, 7,  0, 4, 5, 7,  0, 4, 6, 7], [4, 6])

contains

  !> e(nb, nk**3): band energies, k index = 1 + ix + nk*(iy + nk*iz).
  !> blochl_corr: add the 1/40 curvature correction (improves smooth bands).
  subroutine tetra_weights(nk, e, ef, w, blochl_corr)
    integer, intent(in) :: nk
    real(dp), intent(in) :: e(:, :), ef
    real(dp), intent(out) :: w(size(e, 1), size(e, 2))
    logical, intent(in), optional :: blochl_corr
    integer :: nb, ix, iy, iz, it, ib, c, kidx(4), co(3)
    real(dp) :: vt, es(4), ww(4)
    integer :: ord(4)
    logical :: corr
    nb = size(e, 1)
    corr = .false.; if (present(blochl_corr)) corr = blochl_corr
    w = 0.0_dp
    vt = 1.0_dp / (6.0_dp * real(nk, dp)**3)
    do iz = 0, nk - 1
      do iy = 0, nk - 1
        do ix = 0, nk - 1
          do it = 1, 6
            do c = 1, 4
              co = [iand(ctet(c, it), 1), iand(ishft(ctet(c, it), -1), 1), iand(ishft(ctet(c, it), -2), 1)]
              kidx(c) = 1 + modulo(ix + co(1), nk) + nk * (modulo(iy + co(2), nk) + nk * modulo(iz + co(3), nk))
            end do
            do ib = 1, nb
              es = [e(ib, kidx(1)), e(ib, kidx(2)), e(ib, kidx(3)), e(ib, kidx(4))]
              call sort4(es, ord)
              call one_tet(es(ord(1)), es(ord(2)), es(ord(3)), es(ord(4)), ef, vt, ww, corr)
              do c = 1, 4
                w(ib, kidx(ord(c))) = w(ib, kidx(ord(c))) + ww(c)
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine tetra_weights

  !> Fermi level so that  fac * (occupied volume) = nelec.  The tetrahedra are sorted once and the
  !> occupied volume at trial energies uses its closed form (no per-corner weights needed).
  subroutine tetra_fermi(nk, e, nelec, fac, ef, blochl_corr)
    integer, intent(in) :: nk
    real(dp), intent(in) :: e(:, :), nelec, fac
    real(dp), intent(out) :: ef
    logical, intent(in), optional :: blochl_corr
    real(dp), allocatable :: es(:, :)
    real(dp) :: lo, hi, nn, vt, t(4)
    integer :: nb, ix, iy, iz, it, ib, c, kidx(4), co(3), n, ntot, itr, ord(4)
    nb = size(e, 1)
    ntot = nk**3 * 6 * nb
    allocate (es(4, ntot))
    n = 0
    do iz = 0, nk - 1
      do iy = 0, nk - 1
        do ix = 0, nk - 1
          do it = 1, 6
            do c = 1, 4
              co = [iand(ctet(c, it), 1), iand(ishft(ctet(c, it), -1), 1), iand(ishft(ctet(c, it), -2), 1)]
              kidx(c) = 1 + modulo(ix + co(1), nk) + nk * (modulo(iy + co(2), nk) + nk * modulo(iz + co(3), nk))
            end do
            do ib = 1, nb
              t = [e(ib, kidx(1)), e(ib, kidx(2)), e(ib, kidx(3)), e(ib, kidx(4))]
              call sort4(t, ord)
              n = n + 1
              es(:, n) = [t(ord(1)), t(ord(2)), t(ord(3)), t(ord(4))]
            end do
          end do
        end do
      end do
    end do
    vt = 1.0_dp / (6.0_dp * real(nk, dp)**3)
    lo = minval(e) - 1.0e-6_dp; hi = maxval(e) + 1.0e-6_dp
    do itr = 1, 100
      ef = 0.5_dp * (lo + hi)
      nn = 0.0_dp
      !$omp parallel do reduction(+:nn)
      do it = 1, ntot
        nn = nn + occupied_fraction(es(1, it), es(2, it), es(3, it), es(4, it), ef)
      end do
      !$omp end parallel do
      nn = fac * vt * nn
      if (nn > nelec) then; hi = ef; else; lo = ef; end if
      if (hi - lo < 1.0e-12_dp) exit
    end do
    ef = 0.5_dp * (lo + hi)
    if (present(blochl_corr)) then; end if
  end subroutine tetra_fermi

  !> Fraction of a tetrahedron (corners sorted e1<=..<=e4) with energy below e.
  pure function occupied_fraction(e1, e2, e3, e4, e) result(f)
    real(dp), intent(in) :: e1, e2, e3, e4, e
    real(dp) :: f
    if (e <= e1) then
      f = 0.0_dp
    else if (e < e2) then
      f = (e - e1)**3 / ((e2 - e1) * (e3 - e1) * (e4 - e1))
    else if (e < e3) then
      f = ((e2 - e1)**2 + 3.0_dp * (e2 - e1) * (e - e2) + 3.0_dp * (e - e2)**2 &
           - (e3 - e1 + e4 - e2) / ((e3 - e2) * (e4 - e2)) * (e - e2)**3) / ((e3 - e1) * (e4 - e1))
    else if (e < e4) then
      f = 1.0_dp - (e4 - e)**3 / ((e4 - e1) * (e4 - e2) * (e4 - e3))
    else
      f = 1.0_dp
    end if
  end function occupied_fraction

  pure subroutine sort4(a, ord)
    real(dp), intent(in) :: a(4)
    integer, intent(out) :: ord(4)
    integer :: i, j, t
    ord = [1, 2, 3, 4]
    do i = 2, 4
      t = ord(i); j = i - 1
      do while (j >= 1)
        if (a(ord(j)) <= a(t)) exit
        ord(j + 1) = ord(j); j = j - 1
      end do
      ord(j + 1) = t
    end do
  end subroutine sort4

  !> Weights of the 4 sorted corners (e1<=e2<=e3<=e4) of one tetrahedron of volume fraction vt.
  pure subroutine one_tet(e1, e2, e3, e4, e, vt, w, corr)
    real(dp), intent(in) :: e1, e2, e3, e4, e, vt
    real(dp), intent(out) :: w(4)
    logical, intent(in) :: corr
    real(dp) :: c, c1, c2, c3, dos, d21, d31, d41, d32, d42, d43
    w = 0.0_dp
    dos = 0.0_dp
    if (e <= e1) then
      return
    else if (e < e2) then
      d21 = e2 - e1; d31 = e3 - e1; d41 = e4 - e1
      c = vt * (e - e1)**3 / (d21 * d31 * d41)
      w(1) = c * (4.0_dp - (e - e1) * (1.0_dp / d21 + 1.0_dp / d31 + 1.0_dp / d41)) / 4.0_dp
      w(2) = c * (e - e1) / d21 / 4.0_dp
      w(3) = c * (e - e1) / d31 / 4.0_dp
      w(4) = c * (e - e1) / d41 / 4.0_dp
      dos = 3.0_dp * vt * (e - e1)**2 / (d21 * d31 * d41)
    else if (e < e3) then
      d21 = e2 - e1; d31 = e3 - e1; d41 = e4 - e1; d32 = e3 - e2; d42 = e4 - e2
      c1 = vt * (e - e1)**2 / (d41 * d31)
      c2 = vt * (e - e1) * (e - e2) * (e3 - e) / (d41 * d32 * d31)
      c3 = vt * (e - e2)**2 * (e4 - e) / (d42 * d32 * d41)
      w(1) = (c1 + (c1 + c2) * (e3 - e) / d31 + (c1 + c2 + c3) * (e4 - e) / d41) / 4.0_dp
      w(2) = (c1 + c2 + c3 + (c2 + c3) * (e3 - e) / d32 + c3 * (e4 - e) / d42) / 4.0_dp
      w(3) = ((c1 + c2) * (e - e1) / d31 + (c2 + c3) * (e - e2) / d32) / 4.0_dp
      w(4) = ((c1 + c2 + c3) * (e - e1) / d41 + c3 * (e - e2) / d42) / 4.0_dp
      dos = vt / (d31 * d41) * (3.0_dp * d21 + 6.0_dp * (e - e2) - 3.0_dp * (d31 + d42) * (e - e2)**2 / (d32 * d42))
    else if (e < e4) then
      d41 = e4 - e1; d42 = e4 - e2; d43 = e4 - e3
      c = vt * (e4 - e)**3 / (d41 * d42 * d43)
      w(1) = vt / 4.0_dp - c * (e4 - e) / d41 / 4.0_dp
      w(2) = vt / 4.0_dp - c * (e4 - e) / d42 / 4.0_dp
      w(3) = vt / 4.0_dp - c * (e4 - e) / d43 / 4.0_dp
      w(4) = vt / 4.0_dp - c * (4.0_dp - (e4 - e) * (1.0_dp / d41 + 1.0_dp / d42 + 1.0_dp / d43)) / 4.0_dp
      dos = 3.0_dp * vt * (e4 - e)**2 / (d41 * d42 * d43)
    else
      w = vt / 4.0_dp
    end if
    if (corr .and. dos /= 0.0_dp) then
      w(1) = w(1) + dos / 40.0_dp * ((e2 - e1) + (e3 - e1) + (e4 - e1))
      w(2) = w(2) + dos / 40.0_dp * ((e1 - e2) + (e3 - e2) + (e4 - e2))
      w(3) = w(3) + dos / 40.0_dp * ((e1 - e3) + (e2 - e3) + (e4 - e3))
      w(4) = w(4) + dos / 40.0_dp * ((e1 - e4) + (e2 - e4) + (e3 - e4))
    end if
  end subroutine one_tet
end module asr_tetra
