!> Species and k-space (TB-LMTO-ASA) run parameters.  Read by asr_driver ([calculation] method = "kspace") and filled from the
!> structure there; run by asr_lmto_driver.
module asr_input_lmto
  use asr_kinds, only: dp
  implicit none
  private
  public :: lmto_species_t, lmto_input_t

  type :: lmto_species_t
    character(len=:), allocatable :: symbol
    real(dp) :: z = 0.0_dp
    integer :: nval = 0
    integer :: konfig(0:3) = [4, 4, 3, 4]
    real(dp) :: occ_s = 2.0_dp, occ_d = 0.0_dp
    real(dp) :: spin_sign = 1.0_dp              !! sign of the starting spin polarisation (zeta0) of this species
  end type lmto_species_t

  type :: lmto_input_t
    character(len=:), allocatable :: name, lattice, xc, bz, dos_file, params_file, out_file, workdir, infile
    real(dp) :: alat_angstrom = 0.0_dp, ratio = 1.0_dp
    type(lmto_species_t), allocatable :: species(:)
    integer :: spin = 1, nk = 12, maxit = 100
    real(dp) :: smear = 0.01_dp, mixing = 0.3_dp, tolerance = 1.0e-6_dp, zeta0 = 0.4_dp, sigma = 0.7_dp
    logical :: relativistic = .true.
    logical :: trace = .false.
    real(dp) :: emin = -1.0_dp, emax = 1.0_dp, dos_smear = 0.01_dp
    integer :: npoints = 801, dos_nk = 24
    real(dp) :: cell_vec(3, 3) = 0.0_dp                !! lattice = "cell": columns a1, a2, a3 (alat units)
    real(dp), allocatable :: site_pos(:, :)            !! (3, nsite), Cartesian, alat units
    integer, allocatable :: site_sp(:)                 !! species index of each site
    logical :: site_classes = .false.
  end type lmto_input_t
end module asr_input_lmto
