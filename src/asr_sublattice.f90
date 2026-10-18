!> Decorated lattices for the site-class augmented space: the sites of an sc / bcc / fcc graph (asr_lattices, integer
!> half-lattice coordinates) are divided into sublattices (classes) by a decoration cell of `period` conventional cubes.
!> A sublattice is given by its Wyckoff positions (fractional coordinates of the decoration cell) and the centering
!> translations of the cell ('P', 'I', 'F').  Examples on the bcc graph with period 2: L2_1 X2YZ (F: Z (0,0,0), Y (1/2,1/2,1/2),
!> X (1/4,1/4,1/4) + (3/4,3/4,3/4)), D0_3, B2-type partial order; on fcc with period 1: L1_2 (P: (0,0,0) | the three face centres).
!> Every graph site of the decoration cell must belong to exactly one sublattice.
module asr_sublattice
  use asr_kinds, only: dp
  use asr_ewald, only: madelung_phi
  implicit none
  private
  public :: decoration_t, decoration_build, decoration_classes, decoration_madelung

  type :: decoration_t
    character(len=8) :: lattice = 'bcc'
    integer :: period = 1, m = 2, ncls = 0, nsite = 0      !! m = 2 period: cell edge in half-lattice units
    integer, allocatable :: tab(:, :, :)                   !! class of every graph point of the cell (0: not a site)
    integer, allocatable :: rep(:, :)                      !! rep(:, class): a representative point (half-lattice units)
    real(dp), allocatable :: frac(:)                       !! fraction of the sites in each class
    integer, allocatable :: pts(:, :), pcls(:)             !! all sites of the cell and their classes
  end type decoration_t

contains

  pure logical function on_graph(kind, h)
    character(len=*), intent(in) :: kind
    integer, intent(in) :: h(3)
    select case (kind)
    case ('sc'); on_graph = all(modulo(h, 2) == 0)
    case ('bcc'); on_graph = modulo(h(1), 2) == modulo(h(2), 2) .and. modulo(h(2), 2) == modulo(h(3), 2)
    case ('fcc'); on_graph = modulo(sum(h), 2) == 0
    case default; on_graph = .false.
    end select
  end function on_graph

  !> pos(3, k), owner(k): Wyckoff positions (fractional, decoration cell) and their sublattice 1..ncls; centering 'P'|'I'|'F'.
  !> msg is empty on success.
  subroutine decoration_build(lattice, period, centering, pos, owner, ncls, dec, msg)
    character(len=*), intent(in) :: lattice, centering
    integer, intent(in) :: period, owner(:), ncls
    real(dp), intent(in) :: pos(:, :)
    type(decoration_t), intent(out) :: dec
    character(len=:), allocatable, intent(out) :: msg
    real(dp), allocatable :: tr(:, :)
    real(dp) :: f(3)
    integer :: k, it, h(3), i, j, l, n
    msg = ''
    dec%lattice = lattice; dec%period = period; dec%m = 2 * period; dec%ncls = ncls
    select case (centering)
    case ('P'); tr = reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1])
    case ('I'); tr = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp], [3, 2])
    case ('F'); tr = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp], [3, 4])
    case default; msg = "centering must be 'P', 'I' or 'F'"; return
    end select
    allocate (dec%tab(0:dec%m - 1, 0:dec%m - 1, 0:dec%m - 1)); dec%tab = 0
    do k = 1, size(pos, 2)
      do it = 1, size(tr, 2)
        f = modulo(pos(:, k) + tr(:, it), 1.0_dp) * dec%m
        h = modulo(nint(f), dec%m)
        if (maxval(abs(f - nint(f))) > 1.0e-6_dp) then; msg = 'a sublattice position is not a site of the graph lattice'; return; end if
        if (.not. on_graph(lattice, h)) then; msg = 'a sublattice position is not a site of the graph lattice'; return; end if
        if (dec%tab(h(1), h(2), h(3)) /= 0 .and. dec%tab(h(1), h(2), h(3)) /= owner(k)) then
          msg = 'two sublattices share a site'; return
        end if
        dec%tab(h(1), h(2), h(3)) = owner(k)
      end do
    end do
    n = 0
    do l = 0, dec%m - 1; do j = 0, dec%m - 1; do i = 0, dec%m - 1
      if (.not. on_graph(lattice, [i, j, l])) cycle
      n = n + 1
      if (dec%tab(i, j, l) == 0) then; msg = 'a site of the decoration cell belongs to no sublattice'; return; end if
    end do; end do; end do
    dec%nsite = n
    allocate (dec%pts(3, n), dec%pcls(n), dec%rep(3, ncls), dec%frac(ncls)); dec%rep = -1; dec%frac = 0.0_dp
    n = 0
    do l = 0, dec%m - 1; do j = 0, dec%m - 1; do i = 0, dec%m - 1
      if (.not. on_graph(lattice, [i, j, l])) cycle
      n = n + 1; dec%pts(:, n) = [i, j, l]; dec%pcls(n) = dec%tab(i, j, l)
      if (dec%rep(1, dec%pcls(n)) < 0) dec%rep(:, dec%pcls(n)) = [i, j, l]
      dec%frac(dec%pcls(n)) = dec%frac(dec%pcls(n)) + 1.0_dp
    end do; end do; end do
    if (any(dec%frac == 0.0_dp)) then; msg = 'a sublattice has no site'; return; end if
    dec%frac = dec%frac / real(n, dp)
  end subroutine decoration_build

  !> classes of the graph sites (coord: half-lattice units of a periodic supercell of edge mg) with the representative site of
  !> class `origin` moved to the graph origin.  mg must be a multiple of the decoration cell edge.
  subroutine decoration_classes(dec, coord, mg, origin, cls)
    type(decoration_t), intent(in) :: dec
    integer, intent(in) :: coord(:, :), mg, origin
    integer, allocatable, intent(out) :: cls(:)
    integer :: i, h(3)
    if (modulo(mg, dec%m) /= 0) error stop 'decoration_classes: the supercell is not a multiple of the decoration cell'
    allocate (cls(size(coord, 2)))
    do i = 1, size(coord, 2)
      h = modulo(coord(:, i) + dec%rep(:, origin), dec%m)
      cls(i) = dec%tab(h(1), h(2), h(3))
      if (cls(i) == 0) error stop 'decoration_classes: graph site outside the decoration'
    end do
  end subroutine decoration_classes

  !> Madelung shifts per class (Ry) of the class-average net charges qbar(class):  V_k = sum_k' phi(k, k') qbar(k'),
  !> phi(k, k') = 2 <sum_{j in k'} Phi(i, j)>_{i in k}  over the sites of the decoration cell (Ewald, asr_ewald; alat = graph
  !> lattice constant, Bohr).  The average over i is exact when a class is one Wyckoff orbit.  Charges must be neutral per cell.
  subroutine decoration_madelung(dec, alat, phi)
    type(decoration_t), intent(in) :: dec
    real(dp), intent(in) :: alat
    real(dp), allocatable, intent(out) :: phi(:, :)
    real(dp), allocatable :: tau(:, :), ph(:, :)
    real(dp) :: avec(3, 3)
    integer :: k, i, j
    allocate (tau(3, dec%nsite), ph(dec%nsite, dec%nsite), phi(dec%ncls, dec%ncls)); phi = 0.0_dp
    tau = 0.5_dp * real(dec%pts, dp)
    avec = 0.0_dp; do k = 1, 3; avec(k, k) = real(dec%period, dp); end do
    call madelung_phi(avec, tau, alat, ph)
    do i = 1, dec%nsite
      do j = 1, dec%nsite
        phi(dec%pcls(i), dec%pcls(j)) = phi(dec%pcls(i), dec%pcls(j)) + 2.0_dp * ph(i, j)
      end do
    end do
    do k = 1, dec%ncls
      phi(k, :) = phi(k, :) / (dec%frac(k) * dec%nsite)
    end do
  end subroutine decoration_madelung
end module asr_sublattice
