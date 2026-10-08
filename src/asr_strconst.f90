!> Bare canonical (kappa=0) and screened tight-binding (TB-LMTO) structure
!> constants for s,p,d.  Written from the mstrx2 one-centre expansion
!> (Stuttgart TBSTR) and the Andersen-Jepsen screening  S^a = S0 (1 - a S0)^-1.
!>
!> Convention (mstrx2 / INFO_<lat>):  S0_ss(d) = -2 w/d.  Lengths in units of
!> the lattice constant `alat` unless stated; w = average WS radius.
module asr_strconst
  use asr_kinds, only: dp, pi
  use asr_harmonics, only: ylm_real, gaunt_table, lm_idx
  implicit none
  private
  public :: avg_ws_radius, hardcore_alpha, standard_alpha, lattice_cluster, bare_s0_block, screened_cluster, nl_spd

  integer, parameter :: nl_spd = 9
  integer, parameter :: orb_l(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]
  real(dp), allocatable, save :: gt(:, :, :)

contains

  !> Average WS radius in units of alat for the 1-, 2-, 4-atom cubic cells.
  function avg_ws_radius(kind) result(w)
    character(*), intent(in) :: kind
    real(dp) :: w
    select case (kind)
    case ('sc');  w = (3.0_dp / (4.0_dp * pi))**(1.0_dp / 3.0_dp)
    case ('bcc'); w = (3.0_dp / (8.0_dp * pi))**(1.0_dp / 3.0_dp)
    case ('fcc'); w = (3.0_dp / (16.0_dp * pi))**(1.0_dp / 3.0_dp)
    case default; error stop 'avg_ws_radius: unknown lattice'
    end select
  end function avg_ws_radius

  !> Standard TB-LMTO screening constants (s, p, d) in units of w^(2l+1), Andersen & Jepsen, PRL 53, 2571 (1984);
  !> the representation used by the TB-LMTO-ASR references (Chakrabarti & Mookerjee 2005; thesis eq. 2.60 prints
  !> alpha_d as 0.00107, a typo for 0.010714).  Close to hard spheres with sigma_l = 0.697, 0.683, 0.640.
  pure function standard_alpha() result(al)
    real(dp) :: al(0:2)
    al = [0.3485_dp, 0.05303_dp, 0.010714_dp]
  end function standard_alpha

  !> Hard-sphere screening constants for sphere radius sigma*w:
  !>   alpha_l = sigma^(2l+1) / (2 (2l+1))   (s, p, d, f)
  pure function hardcore_alpha(sigma) result(al)
    real(dp), intent(in) :: sigma
    real(dp) :: al(0:3)
    integer :: l
    do l = 0, 3
      al(l) = sigma**(2 * l + 1) / (2.0_dp * (2 * l + 1))
    end do
  end function hardcore_alpha

  !> Lattice vectors (units of alat) of all sites within `rmax`, origin first.
  subroutine lattice_cluster(kind, rmax, vec)
    character(*), intent(in) :: kind
    real(dp), intent(in) :: rmax
    real(dp), allocatable, intent(out) :: vec(:, :)
    real(dp), allocatable :: tmp(:, :)
    integer :: n, i, j, k, nmax, ip(3)
    real(dp) :: p(3), d
    nmax = ceiling(2.0_dp * rmax) + 1
    allocate (tmp(3, (2 * nmax + 1)**3)); n = 0
    do k = -nmax, nmax; do j = -nmax, nmax; do i = -nmax, nmax
          ip = [i, j, k]
          select case (kind)         ! integer coordinates in half alat
          case ('sc');  if (any(modulo(ip, 2) /= 0)) cycle
          case ('bcc'); if (any(modulo(ip, 2) /= modulo(i, 2))) cycle
          case ('fcc'); if (modulo(sum(ip), 2) /= 0) cycle
          case default; error stop 'lattice_cluster: unknown lattice'
          end select
          p = 0.5_dp * real(ip, dp); d = sqrt(sum(p * p))
          if (d > rmax + 1.0e-9_dp) cycle
          n = n + 1; tmp(:, n) = p
    end do; end do; end do
    allocate (vec(3, n))
    ! origin first, then by distance
    vec = tmp(:, 1:n)
    call sort_by_distance(vec)
  end subroutine lattice_cluster

  subroutine sort_by_distance(v)
    real(dp), intent(inout) :: v(:, :)
    integer :: i, j, n
    real(dp) :: t(3), di, dj
    n = size(v, 2)
    do i = 2, n
      t = v(:, i); di = sum(t * t); j = i - 1
      do while (j >= 1)
        dj = sum(v(:, j) * v(:, j))
        if (dj <= di + 1.0e-12_dp) exit
        v(:, j + 1) = v(:, j); j = j - 1
      end do
      v(:, j + 1) = t
    end do
  end subroutine sort_by_distance

  !> Bare 9x9 block S0_{LL'}(d) for bond vector d (units of alat), w in the same units:
  !>  S0_{LL'} = -8 pi (-1)^{l'} (2l''-1)!! / [(2l-1)!! (2l'-1)!!] (w/d)^{l''+1}
  !>             sum_{m''} G(L,L',l''m'') Y_{l''m''}(dhat),   l'' = l + l'
  subroutine bare_s0_block(d, w, blk)
    real(dp), intent(in) :: d(3), w
    real(dp), intent(out) :: blk(9, 9)
    real(dp) :: yv(25), r, df(0:4), fac, s
    integer :: i, j, li, lj, ll, mm
    if (.not. allocated(gt)) call gaunt_table(2, gt)
    r = sqrt(sum(d * d))
    call ylm_real(4, d(1), d(2), d(3), yv)
    df(0) = 1.0_dp                      ! (2l-1)!!, with (-1)!! = 1
    do i = 1, 4; df(i) = df(i - 1) * (2 * i - 1); end do
    do j = 1, 9
      lj = orb_l(j)
      do i = 1, 9
        li = orb_l(i); ll = li + lj
        s = 0.0_dp
        do mm = -ll, ll
          s = s + gt(i, j, lm_idx(ll, mm)) * yv(lm_idx(ll, mm))
        end do
        fac = -8.0_dp * pi * (-1.0_dp)**lj * df(ll) / (df(li) * df(lj)) * (w / r)**(ll + 1)
        blk(i, j) = fac * s
      end do
    end do
  end subroutine bare_s0_block

  !> Screened S^alpha_{0R} for every cluster site R (blocks 9x9, blk(:,:,iR)),
  !> by solving (1 - S0 alpha) S^alpha = S0 on the finite cluster.
  !> alpha(0:2) for s,p,d.  The row of the origin site is returned.
  subroutine screened_cluster(vec, w, alpha, blk)
    real(dp), intent(in) :: vec(:, :), w, alpha(0:2)
    real(dp), allocatable, intent(out) :: blk(:, :, :)
    integer :: n, nn, i, j, a, b, info
    real(dp), allocatable :: s0(:, :), m(:, :), x(:, :)
    integer, allocatable :: ipiv(:)
    real(dp) :: b9(9, 9), ad(9)
    external :: dgesv
    n = size(vec, 2); nn = 9 * n
    allocate (s0(nn, nn), m(nn, nn), x(nn, 9), ipiv(nn), blk(9, 9, n))
    s0 = 0.0_dp
    do j = 1, n
      do i = 1, n
        if (i == j) cycle
        call bare_s0_block(vec(:, i) - vec(:, j), w, b9)
        s0(9 * (i - 1) + 1:9 * i, 9 * (j - 1) + 1:9 * j) = b9
      end do
    end do
    do a = 1, 9; ad(a) = alpha(orb_l(a)); end do
    ! M = 1 - S0 alpha  (alpha diagonal => scale columns)
    m = -s0
    do j = 1, n; do b = 1, 9
        m(:, 9 * (j - 1) + b) = m(:, 9 * (j - 1) + b) * ad(b)
    end do; end do
    do i = 1, nn; m(i, i) = m(i, i) + 1.0_dp; end do
    ! rhs: columns of S0 belonging to the origin site (site 1)
    x = s0(:, 1:9)
    call dgesv(nn, 9, m, nn, ipiv, x, nn, info)
    if (info /= 0) error stop 'screened_cluster: singular'
    do i = 1, n
      blk(:, :, i) = x(9 * (i - 1) + 1:9 * i, :)
    end do
  end subroutine screened_cluster
end module asr_strconst
