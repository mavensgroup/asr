!> Named crystal structures as data.  This is the only module that knows lattice names: a new structure is one more entry in
!> prototype_structure, or an explicit cell (vectors and basis) in the input, and no other module changes.
!> Cubic prototypes: lengths in units of the cubic lattice constant, axes along the cube edges.  Hexagonal prototypes: lengths
!> in units of a, a1 along x, a2 at 120 degrees in the xy plane, c along z.
module asr_prototypes
  use asr_kinds, only: dp
  use asr_crystal, only: structure_t, structure_init, cell_volume
  implicit none
  private
  public :: prototype_structure, prototype_bravais_data, match_cubic_bravais

  !> The names prototype_structure accepts.
  character(len=*), parameter, public :: prototype_names = 'sc bcc fcc b2 hex hcp nias'

contains

  !> Structure of a named prototype.  c_over_a is required for 'hex' and 'nias'; 'hcp' defaults to the ideal sqrt(8/3).
  !> msg is empty on success.
  subroutine prototype_structure(name, s, msg, c_over_a)
    character(len=*), intent(in) :: name
    type(structure_t), intent(out) :: s
    character(len=:), allocatable, intent(out) :: msg
    real(dp), intent(in), optional :: c_over_a
    real(dp), parameter :: h = 0.5_dp, t1 = 1.0_dp / 3.0_dp, t2 = 2.0_dp / 3.0_dp, q1 = 0.25_dp, q3 = 0.75_dp
    real(dp) :: ca
    msg = ''
    select case (trim(name))
    case ('sc')
      call structure_init(s, cubic(), reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1]), [1], msg)
    case ('bcc')
      call structure_init(s, reshape([-h, h, h, h, -h, h, h, h, -h], [3, 3]), reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1]), [1], msg)
    case ('fcc')
      call structure_init(s, reshape([0.0_dp, h, h, h, 0.0_dp, h, h, h, 0.0_dp], [3, 3]), &
                          reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1]), [1], msg)
    case ('b2')               ! CsCl: two simple cubic sublattices
      call structure_init(s, cubic(), reshape([0.0_dp, 0.0_dp, 0.0_dp, h, h, h], [3, 2]), [1, 2], msg)
    case ('hex')              ! simple hexagonal, one site
      if (.not. present(c_over_a)) then; msg = 'prototype hex: c_over_a is required'; return; end if
      call structure_init(s, hexagonal(c_over_a), reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1]), [1], msg)
    case ('hcp')              ! two equivalent sites
      ca = sqrt(8.0_dp / 3.0_dp); if (present(c_over_a)) ca = c_over_a
      call structure_init(s, hexagonal(ca), reshape([0.0_dp, 0.0_dp, 0.0_dp, t1, t2, h], [3, 2]), [1, 1], msg)
    case ('nias')             ! NiAs (B8_1): cation sites 2a (label 1), anion sites 2c (label 2)
      if (.not. present(c_over_a)) then; msg = 'prototype nias: c_over_a is required'; return; end if
      call structure_init(s, hexagonal(c_over_a), &
                          reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, h, t1, t2, q1, t2, t1, q3], [3, 4]), [1, 1, 2, 2], msg)
    case default
      msg = 'unknown prototype '//trim(name)//' (known: '//prototype_names//')'
    end select
  end subroutine prototype_structure

  pure function cubic() result(a)
    real(dp) :: a(3, 3)
    a = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
  end function cubic

  pure function hexagonal(ca) result(a)
    real(dp), intent(in) :: ca
    real(dp) :: a(3, 3)
    a = reshape([1.0_dp, 0.0_dp, 0.0_dp, -0.5_dp, sqrt(3.0_dp) / 2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, ca], [3, 3])
  end function hexagonal

  !> Primitive vectors p_prim (columns, alat units) and zone-boundary ordering wavevectors k_special (columns, Cartesian units of
  !> 2 pi / alat) of a named cubic Bravais prototype ('sc', 'bcc', or 'fcc').  msg is empty on success.
  subroutine prototype_bravais_data(name, p_prim, k_special, msg)
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: p_prim(3, 3), k_special(3, 3)
    character(len=:), allocatable, intent(out) :: msg
    real(dp), parameter :: h = 0.5_dp
    msg = ''
    p_prim = 0.0_dp; k_special = 0.0_dp
    select case (trim(name))
    case ('sc')
      p_prim = cubic()
      k_special = reshape([h, h, h, h, h, 0.0_dp, h, 0.0_dp, 0.0_dp], [3, 3])
    case ('bcc')
      p_prim = reshape([-h, h, h, h, -h, h, h, h, -h], [3, 3])
      k_special = reshape([1.0_dp, 0.0_dp, 0.0_dp, h, h, 0.0_dp, h, h, h], [3, 3])
    case ('fcc')
      p_prim = reshape([0.0_dp, h, h, h, 0.0_dp, h, h, h, 0.0_dp], [3, 3])
      k_special = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, h, 0.0_dp, h, h, h], [3, 3])
    case default
      msg = 'unknown prototype '//trim(name)//' (known cubic Bravais: sc bcc fcc)'
    end select
  end subroutine prototype_bravais_data

  !> True if the basis sites of st (labels ignored) form a single sc, bcc or fcc Bravais lattice aligned with the Cartesian cube
  !> axes of cubic lattice constant scale (in units of st%alat).  Returns the exact prototype primitive vectors p_prim (units of
  !> alat * scale) and high-symmetry points k_special (Cartesian units of 2 pi / (alat * scale)) so that the offset k mesh is
  !> identical to the named cubic prototype.
  logical function match_cubic_bravais(st, p_prim, scale, k_special) result(ok)
    type(structure_t), intent(in) :: st
    real(dp), intent(out) :: p_prim(3, 3), scale, k_special(3, 3)
    character(len=:), allocatable :: msg
    real(dp), allocatable :: sh(:, :)
    real(dp) :: vsite, r1, acub, v(3), r
    integer :: b1, b2, i1, i2, i3, nsh, z
    p_prim = 0.0_dp; scale = 1.0_dp; k_special = 0.0_dp; ok = .false.
    r1 = huge(1.0_dp)
    do b1 = 1, st%nb
      do b2 = 1, st%nb
        do i1 = -2, 2; do i2 = -2, 2; do i3 = -2, 2
          if (b1 == b2 .and. i1 == 0 .and. i2 == 0 .and. i3 == 0) cycle
          v = i1 * st%avec(:, 1) + i2 * st%avec(:, 2) + i3 * st%avec(:, 3) + st%tau(:, b2) - st%tau(:, b1)
          r = norm2(v)
          if (r > 1.0e-10_dp) r1 = min(r1, r)
        end do; end do; end do
      end do
    end do
    allocate (sh(3, 128 * st%nb)); nsh = 0
    do b1 = 1, st%nb
      do b2 = 1, st%nb
        do i1 = -2, 2; do i2 = -2, 2; do i3 = -2, 2
          v = i1 * st%avec(:, 1) + i2 * st%avec(:, 2) + i3 * st%avec(:, 3) + st%tau(:, b2) - st%tau(:, b1)
          if (abs(norm2(v) - r1) < 1.0e-7_dp * r1) then
            nsh = nsh + 1; sh(:, nsh) = v
          end if
        end do; end do; end do
      end do
    end do
    if (mod(nsh, st%nb) /= 0) return
    z = nsh / st%nb
    vsite = cell_volume(st) / real(st%nb, dp)
    if (z == 6) then
      acub = r1
      if (abs(vsite - acub**3) < 1.0e-6_dp * acub**3 .and. cubic_shell(sh(:, 1:nsh), acub, 1.0_dp)) then
        call prototype_bravais_data('sc', p_prim, k_special, msg)
        scale = acub; ok = .true.
      end if
    else if (z == 8) then
      acub = 2.0_dp * r1 / sqrt(3.0_dp)
      if (abs(vsite - 0.5_dp * acub**3) < 1.0e-6_dp * acub**3 .and. cubic_shell(sh(:, 1:nsh), acub, 0.5_dp * sqrt(3.0_dp))) then
        call prototype_bravais_data('bcc', p_prim, k_special, msg)
        scale = acub; ok = .true.
      end if
    else if (z == 12) then
      acub = sqrt(2.0_dp) * r1
      if (abs(vsite - 0.25_dp * acub**3) < 1.0e-6_dp * acub**3 .and. cubic_shell(sh(:, 1:nsh), acub, 0.5_dp * sqrt(2.0_dp))) then
        call prototype_bravais_data('fcc', p_prim, k_special, msg)
        scale = acub; ok = .true.
      end if
    end if
  contains
    pure logical function cubic_shell(sv, a, r0)
      real(dp), intent(in) :: sv(:, :), a, r0
      real(dp) :: h(3)
      integer :: j
      cubic_shell = .true.
      do j = 1, size(sv, 2)
        if (abs(norm2(sv(:, j)) - a * r0) > 1.0e-6_dp * a) then; cubic_shell = .false.; return; end if
        h = 2.0_dp * sv(:, j) / a
        if (maxval(abs(h - anint(h))) > 1.0e-6_dp) then; cubic_shell = .false.; return; end if
      end do
    end function cubic_shell
  end function match_cubic_bravais
end module asr_prototypes
