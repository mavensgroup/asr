!> Second-order TB-LMTO-ASA Hamiltonian for an ordered one-atom cell, spd.
!>   H = E_nu + h - h o h,   h = (C - E_nu) + sqrt(Delta) S^alpha sqrt(Delta)
!> in the alpha (screened) representation.  o^alpha = (gamma-alpha)/((gamma-alpha)(E_nu-C)-Delta)
!> (as in the Stuttgart makpph).
module asr_lmto_ham
  use asr_kinds, only: dp, pi, ci
  use asr_operator, only: linop
  implicit none
  private
  public :: potpar_t, orbital_l, ham_k, dos_on_grid, o_alpha, to_alpha, hk_grid, kgrid_op, green_exact_orbital

  !> Potential parameters for l = 0,1,2 *in the representation alpha_rep*.
  !> Stuttgart ctrl files store them in the gamma representation (alpha = gamma);
  !> to_alpha converts to the screening alpha of the structure constants.
  type :: potpar_t
    real(dp) :: enu(0:2) = 0.0_dp, c(0:2) = 0.0_dp, sqrdel(0:2) = 0.0_dp
    real(dp) :: p(0:2) = 0.0_dp, gamma(0:2) = 0.0_dp
    real(dp) :: alpha_rep(0:2) = 0.0_dp
  end type potpar_t
  integer, parameter :: orbital_l(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]

  !> H(k) on an nk^3 uniform grid; one complex 9x9 per k, as a linear operator on
  !> real vectors of length 2*9*Nk (re/im doubled).  Lanczos from |0,L> on this
  !> operator gives the recursion coefficients of the crystal (on the grid measure).
  type, extends(linop) :: kgrid_op
    integer :: nk3 = 0
    complex(dp), allocatable :: hk(:, :, :)
  contains
    procedure :: apply => kgrid_apply
    procedure :: dim => kgrid_dim
  end type kgrid_op

contains

  subroutine kgrid_apply(self, x, y)
    class(kgrid_op), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    integer :: ik, i, n
    complex(dp) :: xc(9), yc(9)
    n = self%nk3
    !$omp parallel do private(ik, i, xc, yc) schedule(static)
    do ik = 1, n
      do i = 1, 9
        xc(i) = cmplx(x(9 * (ik - 1) + i), x(9 * n + 9 * (ik - 1) + i), dp)
      end do
      yc = matmul(self%hk(:, :, ik), xc)
      do i = 1, 9
        y(9 * (ik - 1) + i) = real(yc(i), dp)
        y(9 * n + 9 * (ik - 1) + i) = aimag(yc(i))
      end do
    end do
    !$omp end parallel do
  end subroutine kgrid_apply

  pure integer function kgrid_dim(self)
    class(kgrid_op), intent(in) :: self
    kgrid_dim = 18 * self%nk3
  end function kgrid_dim

  subroutine hk_grid(nk, bvec, vec, blk, pp, alpha, op)
    integer, intent(in) :: nk
    real(dp), intent(in) :: bvec(3, 3), vec(:, :), blk(:, :, :), alpha(0:2)
    type(potpar_t), intent(in) :: pp
    type(kgrid_op), intent(out) :: op
    integer :: ix, iy, iz, ik
    real(dp) :: k(3)
    op%nk3 = nk**3
    allocate (op%hk(9, 9, op%nk3))
    !$omp parallel do collapse(2) private(ix, iy, iz, k, ik)
    do iz = 0, nk - 1
      do iy = 0, nk - 1
        do ix = 0, nk - 1
          ik = 1 + ix + nk * (iy + nk * iz)
          k = (ix + 0.5_dp) / nk * bvec(:, 1) + (iy + 0.5_dp) / nk * bvec(:, 2) + (iz + 0.5_dp) / nk * bvec(:, 3)
          call ham_k(k, vec, blk, pp, alpha, 1.0_dp, .true., op%hk(:, :, ik))
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine hk_grid

  !> Exact G_LL(z) = (1/Nk) sum_k sum_b |U_Lb(k)|^2 / (z - e_b(k)) on the same grid.
  function green_exact_orbital(op, orb, z) result(g)
    type(kgrid_op), intent(in) :: op
    integer, intent(in) :: orb
    complex(dp), intent(in) :: z(:)
    complex(dp) :: g(size(z))
    integer :: ik, ib, ie, info
    real(dp) :: ev(9), rwork(100)
    complex(dp) :: h(9, 9), cw(200)
    external :: zheev
    g = (0.0_dp, 0.0_dp)
    !$omp parallel do private(ik, ib, ie, ev, h, cw, rwork, info) reduction(+:g)
    do ik = 1, op%nk3
      h = op%hk(:, :, ik)
      call zheev('V', 'U', 9, h, 9, ev, cw, 200, rwork, info)
      do ib = 1, 9
        do ie = 1, size(z)
          g(ie) = g(ie) + abs(h(orb, ib))**2 / (z(ie) - ev(ib))
        end do
      end do
    end do
    !$omp end parallel do
    g = g / real(op%nk3, dp)
  end function green_exact_orbital

  !> Representation change (Stuttgart pptrns, Varenna p88):
  !>   x = 1 + (C - Enu)(alpha_new - alpha_old)/Delta_old
  !>   C' = Enu + (C-Enu) x,  sqrt(Delta)' = sqrt(Delta) x,  p' = p - o_old^2 + o_new^2
  pure function to_alpha(pp, alpha_new) result(q)
    type(potpar_t), intent(in) :: pp
    real(dp), intent(in) :: alpha_new(0:2)
    type(potpar_t) :: q
    real(dp) :: del, x, oold, onew
    integer :: l
    q = pp
    do l = 0, 2
      del = pp%sqrdel(l)**2
      oold = (pp%gamma(l) - pp%alpha_rep(l)) / ((pp%gamma(l) - pp%alpha_rep(l)) * (pp%enu(l) - pp%c(l)) - del)
      x = 1.0_dp + (pp%c(l) - pp%enu(l)) * (alpha_new(l) - pp%alpha_rep(l)) / del
      q%c(l) = pp%enu(l) + (pp%c(l) - pp%enu(l)) * x
      q%sqrdel(l) = pp%sqrdel(l) * x
      del = q%sqrdel(l)**2
      onew = (pp%gamma(l) - alpha_new(l)) / ((pp%gamma(l) - alpha_new(l)) * (pp%enu(l) - q%c(l)) - del)
      q%p(l) = pp%p(l) - oold * oold + onew * onew
      q%alpha_rep(l) = alpha_new(l)
    end do
  end function to_alpha

  !> o^alpha for parameters already in the representation `alpha`.
  pure function o_alpha(pp, alpha) result(o)
    type(potpar_t), intent(in) :: pp
    real(dp), intent(in) :: alpha(0:2)
    real(dp) :: o(0:2)
    integer :: l
    do l = 0, 2
      o(l) = (pp%gamma(l) - alpha(l)) / ((pp%gamma(l) - alpha(l)) * (pp%enu(l) - pp%c(l)) - pp%sqrdel(l)**2)
    end do
  end function o_alpha

  !> H(k), k in Cartesian units of 2pi/alat; blk(:,:,iR) = S^alpha_{R0} blocks for sites vec(:,iR) (alat).
  !> ssign multiplies S (convention test); use_abs takes |sqrdel|.
  subroutine ham_k(k, vec, blk, pp, alpha, ssign, use_abs, h)
    real(dp), intent(in) :: k(3), vec(:, :), blk(:, :, :), alpha(0:2), ssign
    type(potpar_t), intent(in) :: pp
    logical, intent(in) :: use_abs
    complex(dp), intent(out) :: h(9, 9)
    complex(dp) :: s(9, 9), hh(9, 9)
    type(potpar_t) :: q
    real(dp) :: o(0:2), sd(9), cn(9), en(9), ov(9)
    integer :: i, j, ir, l
    s = (0.0_dp, 0.0_dp)
    do ir = 1, size(vec, 2)
      s = s + transpose(blk(:, :, ir)) * exp(ci * 2.0_dp * pi * dot_product(k, vec(:, ir)))
    end do
    s = 0.5_dp * ssign * (s + conjg(transpose(s)))   ! finite-cluster rows are only Hermitian to truncation error
    q = to_alpha(pp, alpha)
    o = o_alpha(q, alpha)
    do i = 1, 9
      l = orbital_l(i)
      sd(i) = q%sqrdel(l); if (use_abs) sd(i) = abs(sd(i))
      cn(i) = q%c(l) - q%enu(l); en(i) = q%enu(l); ov(i) = o(l)
    end do
    do j = 1, 9; do i = 1, 9
        hh(i, j) = sd(i) * s(i, j) * sd(j)
    end do; end do
    do i = 1, 9; hh(i, i) = hh(i, i) + cn(i); end do
    ! H = E + h - h o h
    do j = 1, 9; do i = 1, 9
        h(i, j) = hh(i, j) - sum(hh(i, :) * ov * hh(:, j))
    end do; end do
    do i = 1, 9; h(i, i) = h(i, i) + en(i); end do
  end subroutine ham_k

  !> Gaussian-smeared l-resolved DOS (states/Ry/atom/spin) on `egrid` from an
  !> nk^3 uniform grid over the reciprocal primitive cell (bvec columns, units 2pi/alat).
  subroutine dos_on_grid(nk, bvec, vec, blk, pp, alpha, ssign, use_abs, egrid, sigma, dos_l)
    integer, intent(in) :: nk
    real(dp), intent(in) :: bvec(3, 3), vec(:, :), blk(:, :, :), alpha(0:2), ssign, egrid(:), sigma
    type(potpar_t), intent(in) :: pp
    logical, intent(in) :: use_abs
    real(dp), intent(out) :: dos_l(size(egrid), 0:2)
    integer :: ix, iy, iz, ie, ib, l, n, info
    real(dp) :: k(3), ev(9), w(9, 0:2), rwork(100), x, g
    complex(dp) :: h(9, 9), cw(200)
    external :: zheev
    dos_l = 0.0_dp
    !$omp parallel do collapse(2) private(ix,iy,iz,k,h,ev,cw,rwork,info,ib,l,w,ie,x,g) reduction(+:dos_l)
    do iz = 0, nk - 1
      do iy = 0, nk - 1
        do ix = 0, nk - 1
          k = (ix + 0.5_dp) / nk * bvec(:, 1) + (iy + 0.5_dp) / nk * bvec(:, 2) + (iz + 0.5_dp) / nk * bvec(:, 3)
          call ham_k(k, vec, blk, pp, alpha, ssign, use_abs, h)
          call zheev('V', 'U', 9, h, 9, ev, cw, 200, rwork, info)
          do ib = 1, 9
            do l = 0, 2
              w(ib, l) = 0.0_dp
            end do
            do n = 1, 9
              w(ib, orbital_l(n)) = w(ib, orbital_l(n)) + abs(h(n, ib))**2
            end do
            do ie = 1, size(egrid)
              x = (egrid(ie) - ev(ib)) / sigma
              if (abs(x) > 6.0_dp) cycle
              g = exp(-0.5_dp * x * x) / (sigma * sqrt(2.0_dp * pi))
              do l = 0, 2
                dos_l(ie, l) = dos_l(ie, l) + g * w(ib, l)
              end do
            end do
          end do
        end do
      end do
    end do
    !$omp end parallel do
    dos_l = dos_l / real(nk, dp)**3
  end subroutine dos_on_grid
end module asr_lmto_ham
