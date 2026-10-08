!> LDA exchange-correlation through libxc's C ABI (no libxc Fortran module needed).
!> Input density in electrons/bohr^3; outputs in Rydberg (libxc returns Hartree).
module asr_xc
  use, intrinsic :: iso_c_binding
  use asr_kinds, only: dp
  implicit none
  private
  public :: xc_lda_t, xc_init, xc_eval, xc_eval_pol, xc_finalize

  integer(c_int), parameter :: XC_LDA_X = 1, XC_LDA_C_VWN = 7, XC_LDA_C_PW = 12, XC_LDA_C_PW_MOD = 13, XC_LDA_C_VBH = 17

  interface
    function c_xc_func_alloc() bind(C, name='xc_func_alloc') result(p)
      import :: c_ptr
      type(c_ptr) :: p
    end function
    function c_xc_func_init(p, functional, nspin) bind(C, name='xc_func_init') result(ierr)
      import :: c_ptr, c_int
      type(c_ptr), value :: p
      integer(c_int), value :: functional, nspin
      integer(c_int) :: ierr
    end function
    subroutine c_xc_func_end(p) bind(C, name='xc_func_end')
      import :: c_ptr
      type(c_ptr), value :: p
    end subroutine
    subroutine c_xc_func_free(p) bind(C, name='xc_func_free')
      import :: c_ptr
      type(c_ptr), value :: p
    end subroutine
    subroutine c_xc_lda_exc_vxc(p, np, rho, zk, vrho) bind(C, name='xc_lda_exc_vxc')
      import :: c_ptr, c_size_t, c_double
      type(c_ptr), value :: p
      integer(c_size_t), value :: np
      real(c_double), intent(in) :: rho(*)
      real(c_double), intent(out) :: zk(*), vrho(*)
    end subroutine
  end interface

  type :: xc_lda_t
    type(c_ptr) :: px = c_null_ptr, pc = c_null_ptr          !! unpolarised
    type(c_ptr) :: px2 = c_null_ptr, pc2 = c_null_ptr        !! spin-polarised
    character(len=8) :: name = ''
  end type xc_lda_t

contains

  !> name: 'vwn' | 'pw' | 'pwmod' | 'vbh'  (exchange is always Slater)
  subroutine xc_init(f, name)
    type(xc_lda_t), intent(out) :: f
    character(len=*), intent(in) :: name
    integer(c_int) :: idc, ierr
    select case (trim(name))
    case ('vwn'); idc = XC_LDA_C_VWN
    case ('pw'); idc = XC_LDA_C_PW
    case ('pwmod'); idc = XC_LDA_C_PW_MOD
    case ('vbh'); idc = XC_LDA_C_VBH
    case default; error stop 'xc_init: unknown functional'
    end select
    f%name = name
    f%px = c_xc_func_alloc(); ierr = c_xc_func_init(f%px, XC_LDA_X, 1_c_int)
    if (ierr /= 0) error stop 'xc_init: exchange init failed'
    f%pc = c_xc_func_alloc(); ierr = c_xc_func_init(f%pc, idc, 1_c_int)
    if (ierr /= 0) error stop 'xc_init: correlation init failed'
    f%px2 = c_xc_func_alloc(); ierr = c_xc_func_init(f%px2, XC_LDA_X, 2_c_int)
    if (ierr /= 0) error stop 'xc_init: polarised exchange init failed'
    f%pc2 = c_xc_func_alloc(); ierr = c_xc_func_init(f%pc2, idc, 2_c_int)
    if (ierr /= 0) error stop 'xc_init: polarised correlation init failed'
  end subroutine xc_init

  !> exc: energy per electron (Ry), vxc: potential (Ry), nspin = 1.
  subroutine xc_eval(f, rho, exc, vxc)
    type(xc_lda_t), intent(in) :: f
    real(dp), intent(in) :: rho(:)
    real(dp), intent(out) :: exc(size(rho)), vxc(size(rho))
    real(c_double) :: r(size(rho)), zx(size(rho)), vx(size(rho)), zc(size(rho)), vc(size(rho))
    r = max(rho, 1.0e-30_dp)
    call c_xc_lda_exc_vxc(f%px, int(size(rho), c_size_t), r, zx, vx)
    call c_xc_lda_exc_vxc(f%pc, int(size(rho), c_size_t), r, zc, vc)
    exc = 2.0_dp * (zx + zc); vxc = 2.0_dp * (vx + vc)
  end subroutine xc_eval

  !> Spin-polarised: rho_up, rho_dn (e/bohr^3) -> exc per electron (Ry), vxc_up, vxc_dn (Ry).
  subroutine xc_eval_pol(f, rup, rdn, exc, vup, vdn)
    type(xc_lda_t), intent(in) :: f
    real(dp), intent(in) :: rup(:), rdn(:)
    real(dp), intent(out) :: exc(size(rup)), vup(size(rup)), vdn(size(rup))
    real(c_double) :: r(2 * size(rup)), zx(size(rup)), vx(2 * size(rup)), zc(size(rup)), vc(2 * size(rup))
    integer :: i
    do i = 1, size(rup)
      r(2 * i - 1) = max(rup(i), 1.0e-30_dp); r(2 * i) = max(rdn(i), 1.0e-30_dp)
    end do
    call c_xc_lda_exc_vxc(f%px2, int(size(rup), c_size_t), r, zx, vx)
    call c_xc_lda_exc_vxc(f%pc2, int(size(rup), c_size_t), r, zc, vc)
    exc = 2.0_dp * (zx + zc)
    do i = 1, size(rup)
      vup(i) = 2.0_dp * (vx(2 * i - 1) + vc(2 * i - 1)); vdn(i) = 2.0_dp * (vx(2 * i) + vc(2 * i))
    end do
  end subroutine xc_eval_pol

  subroutine xc_finalize(f)
    type(xc_lda_t), intent(inout) :: f
    if (c_associated(f%px)) then; call c_xc_func_end(f%px); call c_xc_func_free(f%px); f%px = c_null_ptr; end if
    if (c_associated(f%pc)) then; call c_xc_func_end(f%pc); call c_xc_func_free(f%pc); f%pc = c_null_ptr; end if
    if (c_associated(f%px2)) then; call c_xc_func_end(f%px2); call c_xc_func_free(f%px2); f%px2 = c_null_ptr; end if
    if (c_associated(f%pc2)) then; call c_xc_func_end(f%pc2); call c_xc_func_free(f%pc2); f%pc2 = c_null_ptr; end if
  end subroutine xc_finalize
end module asr_xc
