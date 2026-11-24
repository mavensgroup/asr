!> Named crystal structures as data.  This is the only module that knows lattice names: a new structure is one more entry in
!> prototype_structure, or an explicit cell (vectors and basis) in the input, and no other module changes.
!> Cubic prototypes: lengths in units of the cubic lattice constant, axes along the cube edges.  Hexagonal prototypes: lengths
!> in units of a, a1 along x, a2 at 120 degrees in the xy plane, c along z.
module asr_prototypes
  use asr_kinds, only: dp
  use asr_crystal, only: structure_t, structure_init
  implicit none
  private
  public :: prototype_structure

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
end module asr_prototypes
