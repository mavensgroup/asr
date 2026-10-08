!> General ordered cell (any lattice, any basis, per-species parameters and
!> screening) for the second-order TB-LMTO-ASA Hamiltonian, spd.
!> Structure constants: bare canonical S0 on a finite cluster around each basis
!> atom, screened with per-site hard-core alpha:  (1 - S0 alpha) S^alpha = S0.
module asr_lmto_cell
  use asr_kinds, only: dp, pi, ci
  use asr_strconst, only: bare_s0_block
  use asr_lmto_ham, only: potpar_t, to_alpha, o_alpha
  implicit none
  private
  public :: cell_t, cluster_t, build_cluster, ham_k_cell, dos_cell, recip_vectors, alpha_from_sigma, setup_cell, hsmall_k

  integer, parameter :: orb_l(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]

  type :: cell_t
    real(dp) :: avec(3, 3) = 0.0_dp          !! primitive vectors (columns), units of alat
    integer :: nat = 0
    real(dp), allocatable :: tau(:, :)       !! (3,nat), units of alat
    integer, allocatable :: isp(:)           !! species of each atom
    real(dp) :: w = 0.0_dp                   !! average WS radius, alat units
  end type cell_t

  type :: cluster_t
    integer :: n = 0
    real(dp), allocatable :: d(:, :)         !! position - centre (3,n), alat
    integer, allocatable :: iat(:)           !! basis atom of each site
    real(dp), allocatable :: s(:, :, :)      !! S^alpha_{site, centre}, (9,9,n)
  end type cluster_t

contains

  !> alpha_l = (sigma wsr)^(2l+1) / (2(2l+1)), wsr = sphere radius / avw.
  pure function alpha_from_sigma(sigma, wsr) result(al)
    real(dp), intent(in) :: sigma, wsr
    real(dp) :: al(0:2)
    integer :: l
    do l = 0, 2
      al(l) = (sigma * wsr)**(2 * l + 1) / (2.0_dp * (2 * l + 1))
    end do
  end function alpha_from_sigma

  !> Fill avw from the cell volume and atom count (equal-volume WS spheres).
  subroutine setup_cell(c)
    type(cell_t), intent(inout) :: c
    real(dp) :: v
    v = abs(c%avec(1, 1) * (c%avec(2, 2) * c%avec(3, 3) - c%avec(3, 2) * c%avec(2, 3)) &
            - c%avec(1, 2) * (c%avec(2, 1) * c%avec(3, 3) - c%avec(3, 1) * c%avec(2, 3)) &
            + c%avec(1, 3) * (c%avec(2, 1) * c%avec(3, 2) - c%avec(3, 1) * c%avec(2, 2)))
    c%w = (3.0_dp * v / (4.0_dp * pi * c%nat))**(1.0_dp / 3.0_dp)
  end subroutine setup_cell

  !> Reciprocal primitive vectors (columns) in units of 2pi/alat.
  function recip_vectors(avec) result(b)
    real(dp), intent(in) :: avec(3, 3)
    real(dp) :: b(3, 3), v
    b(:, 1) = cross(avec(:, 2), avec(:, 3))
    b(:, 2) = cross(avec(:, 3), avec(:, 1))
    b(:, 3) = cross(avec(:, 1), avec(:, 2))
    v = dot_product(avec(:, 1), b(:, 1))
    b = b / v
  end function recip_vectors

  pure function cross(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
  end function cross

  !> Cluster of all sites within rmax (alat) of atom b, screened with alpha_sp(:, species).
  subroutine build_cluster(c, b, rmax, alpha_sp, cl)
    type(cell_t), intent(in) :: c
    integer, intent(in) :: b
    real(dp), intent(in) :: rmax, alpha_sp(0:, :)
    type(cluster_t), intent(out) :: cl
    integer :: n1, n2, n3, a, nmax, n, i, j, nn, info, q, l
    real(dp) :: p(3), dist
    real(dp), allocatable :: tmp(:, :), s0(:, :), m(:, :), x(:, :)
    integer, allocatable :: tia(:), ipiv(:)
    real(dp) :: b9(9, 9), ad(9)
    external :: dgesv
    nmax = ceiling(rmax / minval([(norm2(c%avec(:, i)), i=1, 3)])) + 2
    allocate (tmp(3, c%nat * (2 * nmax + 1)**3), tia(c%nat * (2 * nmax + 1)**3)); n = 0
    do n3 = -nmax, nmax; do n2 = -nmax, nmax; do n1 = -nmax, nmax
          do a = 1, c%nat
            p = n1 * c%avec(:, 1) + n2 * c%avec(:, 2) + n3 * c%avec(:, 3) + c%tau(:, a) - c%tau(:, b)
            dist = norm2(p)
            if (dist > rmax + 1.0e-9_dp) cycle
            n = n + 1; tmp(:, n) = p; tia(n) = a
          end do
    end do; end do; end do
    ! centre first
    cl%n = n
    allocate (cl%d(3, n), cl%iat(n))
    j = 0
    do i = 1, n
      if (norm2(tmp(:, i)) < 1.0e-10_dp .and. tia(i) == b) then
        j = j + 1; cl%d(:, 1) = tmp(:, i); cl%iat(1) = tia(i)
      end if
    end do
    if (j /= 1) error stop 'build_cluster: centre not unique'
    q = 1
    do i = 1, n
      if (norm2(tmp(:, i)) < 1.0e-10_dp .and. tia(i) == b) cycle
      q = q + 1; cl%d(:, q) = tmp(:, i); cl%iat(q) = tia(i)
    end do
    nn = 9 * n
    allocate (s0(nn, nn), m(nn, nn), x(nn, 9), ipiv(nn), cl%s(9, 9, n))
    s0 = 0.0_dp
    do j = 1, n
      do i = 1, n
        if (i == j) cycle
        call bare_s0_block(cl%d(:, i) - cl%d(:, j), c%w, b9)
        s0(9 * (i - 1) + 1:9 * i, 9 * (j - 1) + 1:9 * j) = b9
      end do
    end do
    m = -s0
    do j = 1, n
      do l = 1, 9
        ad(l) = alpha_sp(orb_l(l), c%isp(cl%iat(j)))
        m(:, 9 * (j - 1) + l) = m(:, 9 * (j - 1) + l) * ad(l)
      end do
    end do
    do i = 1, nn; m(i, i) = m(i, i) + 1.0_dp; end do
    x = s0(:, 1:9)
    call dgesv(nn, 9, m, nn, ipiv, x, nn, info)
    if (info /= 0) error stop 'build_cluster: singular'
    do i = 1, n
      cl%s(:, :, i) = x(9 * (i - 1) + 1:9 * i, :)
    end do
  end subroutine build_cluster

  !> H(k) (9nat x 9nat), k Cartesian in 2pi/alat.  pps(species) in the representation
  !> pps%alpha_rep; transformed here to alpha_sp(:,species).
  subroutine ham_k_cell(c, cls, pps, alpha_sp, k, h)
    type(cell_t), intent(in) :: c
    type(cluster_t), intent(in) :: cls(:)
    type(potpar_t), intent(in) :: pps(:)
    real(dp), intent(in) :: alpha_sp(0:, :), k(3)
    complex(dp), intent(out) :: h(:, :)
    complex(dp), allocatable :: s(:, :), hh(:, :), t(:, :)
    type(potpar_t) :: q
    real(dp), allocatable :: sd(:), cn(:), en(:), ov(:)
    real(dp) :: ol(0:2)
    integer :: ib, i, a, nd, ia, l, r0, c0, isp
    nd = 9 * c%nat
    allocate (s(nd, nd), hh(nd, nd), t(nd, nd), sd(nd), cn(nd), en(nd), ov(nd))
    s = (0.0_dp, 0.0_dp)
    do ib = 1, c%nat
      do i = 1, cls(ib)%n
        a = cls(ib)%iat(i)
        r0 = 9 * (ib - 1); c0 = 9 * (a - 1)
        s(r0 + 1:r0 + 9, c0 + 1:c0 + 9) = s(r0 + 1:r0 + 9, c0 + 1:c0 + 9) &
          + transpose(cls(ib)%s(:, :, i)) * exp(ci * 2.0_dp * pi * dot_product(k, cls(ib)%d(:, i)))
      end do
    end do
    s = 0.5_dp * (s + conjg(transpose(s)))
    do ia = 1, c%nat
      isp = c%isp(ia)
      q = to_alpha(pps(isp), alpha_sp(:, isp))
      ol = o_alpha(q, alpha_sp(:, isp))
      do l = 1, 9
        i = 9 * (ia - 1) + l
        sd(i) = q%sqrdel(orb_l(l)); cn(i) = q%c(orb_l(l)) - q%enu(orb_l(l))
        en(i) = q%enu(orb_l(l)); ov(i) = ol(orb_l(l))
      end do
    end do
    do a = 1, nd
      do i = 1, nd
        hh(i, a) = sd(i) * s(i, a) * sd(a)
      end do
      hh(a, a) = hh(a, a) + cn(a)
    end do
    ! H = E + h - h o h
    do a = 1, nd
      t(:, a) = ov * hh(:, a)
    end do
    h = hh - matmul(hh, t)
    do i = 1, nd; h(i, i) = h(i, i) + en(i); end do
  end subroutine ham_k_cell

  !> Small h matrix (alpha representation) and the diagonal o, p, E_nu per orbital, for the
  !> first-order ASA generalised problem  (h + h o h) z = (e - E_nu)(1 + o h + h o + h p h) z.
  !> pps(species) are in the representation pps%alpha_rep and are transformed to alpha_sp here.
  subroutine hsmall_k(c, cls, pps, alpha_sp, k, hh, ov, pv, en)
    type(cell_t), intent(in) :: c
    type(cluster_t), intent(in) :: cls(:)
    type(potpar_t), intent(in) :: pps(:)
    real(dp), intent(in) :: alpha_sp(0:, :), k(3)
    complex(dp), intent(out) :: hh(:, :)
    real(dp), intent(out) :: ov(:), pv(:), en(:)
    complex(dp), allocatable :: s(:, :)
    type(potpar_t) :: q
    real(dp) :: ol(0:2), sd(9 * c%nat), cn(9 * c%nat)
    integer :: ib, i, a, nd, ia, l, r0, c0, isp
    nd = 9 * c%nat
    allocate (s(nd, nd)); s = (0.0_dp, 0.0_dp)
    do ib = 1, c%nat
      do i = 1, cls(ib)%n
        a = cls(ib)%iat(i)
        r0 = 9 * (ib - 1); c0 = 9 * (a - 1)
        s(r0 + 1:r0 + 9, c0 + 1:c0 + 9) = s(r0 + 1:r0 + 9, c0 + 1:c0 + 9) &
          + transpose(cls(ib)%s(:, :, i)) * exp(ci * 2.0_dp * pi * dot_product(k, cls(ib)%d(:, i)))
      end do
    end do
    s = 0.5_dp * (s + conjg(transpose(s)))
    do ia = 1, c%nat
      isp = c%isp(ia)
      q = to_alpha(pps(isp), alpha_sp(:, isp))
      ol = o_alpha(q, alpha_sp(:, isp))
      do l = 1, 9
        i = 9 * (ia - 1) + l
        sd(i) = q%sqrdel(orb_l(l)); cn(i) = q%c(orb_l(l)) - q%enu(orb_l(l))
        en(i) = q%enu(orb_l(l)); ov(i) = ol(orb_l(l)); pv(i) = q%p(orb_l(l))
      end do
    end do
    do a = 1, nd
      do i = 1, nd
        hh(i, a) = sd(i) * s(i, a) * sd(a)
      end do
      hh(a, a) = hh(a, a) + cn(a)
    end do
  end subroutine hsmall_k

  !> Gaussian-smeared (atom, l)-resolved DOS, states/Ry/spin; dos(ne, 0:2, nat).
  subroutine dos_cell(c, cls, pps, alpha_sp, nk, egrid, sigma, dos, nelec_at_e)
    type(cell_t), intent(in) :: c
    type(cluster_t), intent(in) :: cls(:)
    type(potpar_t), intent(in) :: pps(:)
    real(dp), intent(in) :: alpha_sp(0:, :), egrid(:), sigma
    integer, intent(in) :: nk
    real(dp), intent(out) :: dos(size(egrid), 0:2, c%nat)
    real(dp), intent(in), optional :: nelec_at_e
    real(dp) :: bv(3, 3), k(3), x, g, wgt(0:2, c%nat)
    real(dp), allocatable :: ev(:), rwork(:)
    complex(dp), allocatable :: h(:, :), cw(:)
    integer :: ix, iy, iz, ib, ie, nd, ia, l, info
    external :: zheev
    bv = recip_vectors(c%avec); nd = 9 * c%nat
    dos = 0.0_dp
    !$omp parallel private(ev, rwork, h, cw, k, x, g, wgt, ib, ie, ia, l, info, iy, ix)
    allocate (ev(nd), rwork(3 * nd), h(nd, nd), cw(4 * nd))
    !$omp do collapse(2) reduction(+:dos)
    do iz = 0, nk - 1
      do iy = 0, nk - 1
        do ix = 0, nk - 1
          k = (ix + 0.5_dp) / nk * bv(:, 1) + (iy + 0.5_dp) / nk * bv(:, 2) + (iz + 0.5_dp) / nk * bv(:, 3)
          call ham_k_cell(c, cls, pps, alpha_sp, k, h)
          call zheev('V', 'U', nd, h, nd, ev, cw, 4 * nd, rwork, info)
          do ib = 1, nd
            wgt = 0.0_dp
            do ia = 1, c%nat
              do l = 1, 9
                wgt(orb_l(l), ia) = wgt(orb_l(l), ia) + abs(h(9 * (ia - 1) + l, ib))**2
              end do
            end do
            do ie = 1, size(egrid)
              x = (egrid(ie) - ev(ib)) / sigma
              if (abs(x) > 6.0_dp) cycle
              g = exp(-0.5_dp * x * x) / (sigma * sqrt(2.0_dp * pi))
              dos(ie, :, :) = dos(ie, :, :) + g * wgt
            end do
          end do
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
    dos = dos / real(nk, dp)**3
    if (present(nelec_at_e)) then; end if
  end subroutine dos_cell
end module asr_lmto_cell
