!> Cubic point-group (O_h, 48 operations) site permutations about the origin site of a periodic
!> sc/bcc/fcc supercell built by asr_lattices::periodic_graph (integer half-lattice coordinates).
module asr_symmetry
  implicit none
  private
  public :: cubic_site_perms, d4h_ops, orbital_rep

contains

  !> coord(3, n): coordinates in 0..m-1 (periodic).  perm(site, op) = image of `site` under op;
  !> op = signed permutation of the coordinate axes about the origin (0,0,0) which must be a site.
  subroutine cubic_site_perms(coord, m, perm, mats)
    integer, intent(in) :: coord(:, :), m
    integer, allocatable, intent(out) :: perm(:, :)
    integer, allocatable, intent(out), optional :: mats(:, :, :)    !! 3x3 signed permutation matrices, mats(:,:,iop)
    integer, allocatable :: idx(:, :, :)
    integer :: n, i, iop, ip, is, sg(3), pp(3), c(3), q(3)
    integer, parameter :: ax(3, 6) = reshape([1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, 3, 2, 1], [3, 6])
    n = size(coord, 2)
    allocate (idx(0:m - 1, 0:m - 1, 0:m - 1)); idx = 0
    do i = 1, n; idx(coord(1, i), coord(2, i), coord(3, i)) = i; end do
    allocate (perm(n, 48))
    if (present(mats)) then; allocate (mats(3, 3, 48)); mats = 0; end if
    iop = 0
    do ip = 1, 6
      pp = ax(:, ip)
      do is = 0, 7
        sg = [1 - 2 * iand(is, 1), 1 - 2 * iand(ishft(is, -1), 1), 1 - 2 * iand(ishft(is, -2), 1)]
        iop = iop + 1
        if (present(mats)) then
          do i = 1, 3; mats(i, pp(i), iop) = sg(i); end do
        end if
        do i = 1, n
          c = coord(:, i)
          q = [modulo(sg(1) * c(pp(1)), m), modulo(sg(2) * c(pp(2)), m), modulo(sg(3) * c(pp(3)), m)]
          if (idx(q(1), q(2), q(3)) == 0) error stop 'cubic_site_perms: image is not a lattice site'
          perm(i, iop) = idx(q(1), q(2), q(3))
        end do
      end do
    end do
  end subroutine cubic_site_perms

  !> The 16 operations of D4h about the origin (signed permutations that keep the z axis): the common
  !> little group of the s, p_z, t2g(xy) and e_g(3z^2-r^2) start orbitals (s uses the monomial subgroup).
  !> perm(site, k), mats(:,:,k) for k = 1..16.
  subroutine d4h_ops(coord, m, perm, mats)
    integer, intent(in) :: coord(:, :), m
    integer, allocatable, intent(out) :: perm(:, :), mats(:, :, :)
    integer, allocatable :: p48(:, :), m48(:, :, :)
    integer :: k, kk
    call cubic_site_perms(coord, m, p48, m48)
    allocate (perm(size(coord, 2), 16), mats(3, 3, 16))
    kk = 0
    do k = 1, 48
      if (abs(m48(3, 3, k)) == 1) then
        kk = kk + 1; perm(:, kk) = p48(:, k); mats(:, :, kk) = m48(:, :, k)
      end if
    end do
    if (kk /= 16) error stop 'd4h_ops: expected 16 operations'
  end subroutine d4h_ops

  !> Orbital representation of each operation on the 9 real harmonics (s | y z x | xy yz 3z2-r2 xz x2-y2):
  !>   Y_k(M p) = sum_j d(k,j,iop) Y_j(p).   Fitted by least squares on generic points; entries are 0, +-1
  !> for the monomial (D4h) operations.
  subroutine orbital_rep(mats, d)
    use asr_kinds, only: dp
    use asr_harmonics, only: ylm_real
    integer, intent(in) :: mats(:, :, :)
    real(dp), allocatable, intent(out) :: d(:, :, :)
    integer, parameter :: np = 24
    real(dp) :: pts(3, np), a(np, 9), b(np, 9), yv(25), mp(3), ata(9, 9), atb(9, 9), x(9, 9)
    integer :: ip, iop, k, ipiv(9), info
    external :: dgesv
    call random_seed()
    do ip = 1, np
      call random_number(pts(:, ip)); pts(:, ip) = pts(:, ip) - 0.5_dp
      call ylm_real(2, pts(1, ip), pts(2, ip), pts(3, ip), yv(1:9)); a(ip, :) = yv(1:9)
    end do
    allocate (d(9, 9, size(mats, 3)))
    ata = matmul(transpose(a), a)
    do iop = 1, size(mats, 3)
      do ip = 1, np
        mp = matmul(real(mats(:, :, iop), dp), pts(:, ip))
        call ylm_real(2, mp(1), mp(2), mp(3), yv(1:9)); b(ip, :) = yv(1:9)
      end do
      atb = matmul(transpose(a), b)
      x = atb
      call dgesv(9, 9, ata, 9, ipiv, x, 9, info)
      ata = matmul(transpose(a), a)
      d(:, :, iop) = transpose(x)
      do k = 1, 9
        where (abs(d(:, k, iop)) < 1.0e-9_dp) d(:, k, iop) = 0.0_dp
      end do
    end do
  end subroutine orbital_rep
end module asr_symmetry
