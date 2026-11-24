!> Periodic crystal structures: primitive vectors and basis sites, each with a sublattice label.  This is the only description of
!> a lattice inside the program.  Lengths are Cartesian, in units of alat.  Named lattices are data in asr_prototypes; every
!> geometric and symmetric object (graph, neighbour slots, site types, space group, start orbitals) is derived from structure_t
!> by general code, so a new structure never changes a physics module.  A site of the infinite crystal is (n1, n2, n3, b): integer
!> cell indices and a basis index.  Floating point decides distances only, never whether two sites are the same.
module asr_crystal
  use asr_kinds, only: dp, pi
  implicit none
  private
  public :: structure_t, structure_init, to_frac, to_cart, same_frac, cell_volume, ws_radius, det3, inv3

  !> Tolerance for geometric equality (fractional coordinates and alat units).
  real(dp), parameter, public :: geo_tol = 1.0e-6_dp

  type :: structure_t
    real(dp) :: avec(3, 3) = 0.0_dp       !! primitive vectors (columns), alat units
    real(dp) :: bvec(3, 3) = 0.0_dp       !! dual vectors (columns): b_i . a_j = delta_ij (no factor 2 pi)
    integer :: nb = 0                     !! basis sites per primitive cell
    integer :: nlabel = 0                 !! number of sublattice labels (max of label)
    real(dp), allocatable :: frac(:, :)   !! (3, nb) fractional coordinates, wrapped into [0, 1)
    real(dp), allocatable :: tau(:, :)    !! (3, nb) Cartesian positions, alat units (avec frac)
    integer, allocatable :: label(:)      !! sublattice label of each basis site; a symmetry maps a site onto one with its label
  end type structure_t

contains

  !> Structure from primitive vectors (columns, alat units) and fractional basis positions with labels.  The positions are
  !> wrapped into the cell.  msg is empty on success.
  subroutine structure_init(s, avec, frac, label, msg)
    type(structure_t), intent(out) :: s
    real(dp), intent(in) :: avec(3, 3), frac(:, :)
    integer, intent(in) :: label(:)
    character(len=:), allocatable, intent(out) :: msg
    integer :: b, c
    msg = ''
    if (size(frac, 1) /= 3 .or. size(frac, 2) < 1 .or. size(label) /= size(frac, 2)) then
      msg = 'structure_init: frac must be (3, nb) and label (nb) with nb >= 1'; return
    end if
    if (abs(det3(avec)) < geo_tol) then; msg = 'structure_init: the primitive vectors are linearly dependent'; return; end if
    if (any(label < 1)) then; msg = 'structure_init: labels must be positive'; return; end if
    s%avec = avec
    s%bvec = transpose(inv3(avec))
    s%nb = size(frac, 2)
    s%nlabel = maxval(label)
    allocate (s%frac(3, s%nb), s%tau(3, s%nb), s%label(s%nb))
    s%frac = frac - floor(frac + geo_tol)
    s%tau = matmul(avec, s%frac)
    s%label = label
    do b = 1, s%nb
      do c = b + 1, s%nb
        if (same_frac(s%frac(:, b), s%frac(:, c))) then; msg = 'structure_init: two basis sites coincide'; return; end if
      end do
    end do
  end subroutine structure_init

  !> Fractional coordinates of a Cartesian vector.
  pure function to_frac(s, r) result(f)
    type(structure_t), intent(in) :: s
    real(dp), intent(in) :: r(3)
    real(dp) :: f(3)
    f = matmul(transpose(s%bvec), r)
  end function to_frac

  !> Cartesian vector of fractional coordinates.
  pure function to_cart(s, f) result(r)
    type(structure_t), intent(in) :: s
    real(dp), intent(in) :: f(3)
    real(dp) :: r(3)
    r = matmul(s%avec, f)
  end function to_cart

  !> Do two fractional positions differ by a lattice vector?
  pure logical function same_frac(f1, f2)
    real(dp), intent(in) :: f1(3), f2(3)
    real(dp) :: d(3)
    d = f1 - f2
    same_frac = all(abs(d - anint(d)) < geo_tol)
  end function same_frac

  !> Volume of the primitive cell, alat^3.
  pure real(dp) function cell_volume(s)
    type(structure_t), intent(in) :: s
    cell_volume = abs(det3(s%avec))
  end function cell_volume

  !> Average Wigner-Seitz radius (alat units): the sphere volume equals the cell volume per basis site.
  pure real(dp) function ws_radius(s)
    type(structure_t), intent(in) :: s
    ws_radius = (3.0_dp * cell_volume(s) / (4.0_dp * pi * real(s%nb, dp)))**(1.0_dp / 3.0_dp)
  end function ws_radius

  pure real(dp) function det3(a)
    real(dp), intent(in) :: a(3, 3)
    det3 = a(1, 1) * (a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2)) - a(1, 2) * (a(2, 1) * a(3, 3) - a(2, 3) * a(3, 1)) &
           + a(1, 3) * (a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1))
  end function det3

  !> Inverse of a 3 x 3 matrix (the caller ensures det /= 0).
  pure function inv3(a) result(b)
    real(dp), intent(in) :: a(3, 3)
    real(dp) :: b(3, 3), d
    d = det3(a)
    b(1, 1) = a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2); b(1, 2) = a(1, 3) * a(3, 2) - a(1, 2) * a(3, 3)
    b(1, 3) = a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2); b(2, 1) = a(2, 3) * a(3, 1) - a(2, 1) * a(3, 3)
    b(2, 2) = a(1, 1) * a(3, 3) - a(1, 3) * a(3, 1); b(2, 3) = a(1, 3) * a(2, 1) - a(1, 1) * a(2, 3)
    b(3, 1) = a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1); b(3, 2) = a(1, 2) * a(3, 1) - a(1, 1) * a(3, 2)
    b(3, 3) = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)
    b = b / d
  end function inv3
end module asr_crystal
