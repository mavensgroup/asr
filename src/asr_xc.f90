!> LDA and GGA (PBE) exchange-correlation through libxc's C ABI (no libxc Fortran module needed).
!> Input density in electrons/bohr^3; outputs in Rydberg (libxc returns Hartree).
!> GGA on a spherical density n(r): sigma = (dn/dr)^2 and v = de/dn - (1/r^2) d/dr [r^2 2 (de/dsigma) dn/dr]
!> (spin-polarised: v_up = vrho_up - div(2 vsigma_uu grad n_up + vsigma_ud grad n_dn), likewise down); the callers pass
!> the radial mesh r(:) (r(1) = 0 is allowed: the point takes the value of its neighbour).
module asr_xc
  use, intrinsic :: iso_c_binding
  use asr_kinds, only: dp
  implicit none
  private
  public :: xc_lda_t, xc_init, xc_eval, xc_eval_pol, xc_finalize, xc_gga_rcut

  !> GGA: the divergence term of the potential is switched off smoothly inside r ~ xc_gga_rcut (factor 1 - exp(-(r/rc)^2)), where
  !> the nuclear cusp makes it ~ 1/r.  With the full term the crystal (ASA) SCF is unstable (bcc Fe: E_F drifts to +2 Ry; the
  !> mixed smooth potential V + 2Z/r is dominated by the c/r term on the dense inner mesh); 0.01 Bohr (inside the 1s shell of 3d
  !> atoms) converges, and the Fe moment changes by < 0.004 muB between 0.01 and 0.1 Bohr.  0 = full GGA potential everywhere.
  real(dp) :: xc_gga_rcut = 0.01_dp

  integer(c_int), parameter :: XC_LDA_X = 1, XC_LDA_C_VWN = 7, XC_LDA_C_PW = 12, XC_LDA_C_PW_MOD = 13, XC_LDA_C_VBH = 17, &
                                XC_GGA_X_PBE = 101, XC_GGA_C_PBE = 130

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
    subroutine c_xc_gga_exc_vxc(p, np, rho, sigma, zk, vrho, vsigma) bind(C, name='xc_gga_exc_vxc')
      import :: c_ptr, c_size_t, c_double
      type(c_ptr), value :: p
      integer(c_size_t), value :: np
      real(c_double), intent(in) :: rho(*), sigma(*)
      real(c_double), intent(out) :: zk(*), vrho(*), vsigma(*)
    end subroutine
  end interface

  type :: xc_lda_t
    type(c_ptr) :: px = c_null_ptr, pc = c_null_ptr          !! unpolarised
    type(c_ptr) :: px2 = c_null_ptr, pc2 = c_null_ptr        !! spin-polarised
    character(len=8) :: name = ''
    logical :: gga = .false.
  end type xc_lda_t

contains

  !> name: 'vwn' | 'pw' | 'pwmod' | 'vbh' (LDA, Slater exchange) | 'pbe' (GGA, PBE exchange and correlation)
  subroutine xc_init(f, name)
    type(xc_lda_t), intent(out) :: f
    character(len=*), intent(in) :: name
    integer(c_int) :: idc, idx, ierr
    idx = XC_LDA_X
    select case (trim(name))
    case ('vwn'); idc = XC_LDA_C_VWN
    case ('pw'); idc = XC_LDA_C_PW
    case ('pwmod'); idc = XC_LDA_C_PW_MOD
    case ('vbh'); idc = XC_LDA_C_VBH
    case ('pbe'); idc = XC_GGA_C_PBE; idx = XC_GGA_X_PBE; f%gga = .true.
    case default; error stop 'xc_init: unknown functional'
    end select
    f%name = name
    f%px = c_xc_func_alloc(); ierr = c_xc_func_init(f%px, idx, 1_c_int)
    if (ierr /= 0) error stop 'xc_init: exchange init failed'
    f%pc = c_xc_func_alloc(); ierr = c_xc_func_init(f%pc, idc, 1_c_int)
    if (ierr /= 0) error stop 'xc_init: correlation init failed'
    f%px2 = c_xc_func_alloc(); ierr = c_xc_func_init(f%px2, idx, 2_c_int)
    if (ierr /= 0) error stop 'xc_init: polarised exchange init failed'
    f%pc2 = c_xc_func_alloc(); ierr = c_xc_func_init(f%pc2, idc, 2_c_int)
    if (ierr /= 0) error stop 'xc_init: polarised correlation init failed'
  end subroutine xc_init

  !> exc: energy per electron (Ry), vxc: potential (Ry), nspin = 1.  r: radial mesh (required for GGA).
  subroutine xc_eval(f, rho, exc, vxc, r)
    type(xc_lda_t), intent(in) :: f
    real(dp), intent(in) :: rho(:)
    real(dp), intent(out) :: exc(size(rho)), vxc(size(rho))
    real(dp), intent(in), optional :: r(:)
    real(c_double) :: rr(size(rho)), zx(size(rho)), vx(size(rho)), zc(size(rho)), vc(size(rho))
    real(c_double) :: sg(size(rho)), sx(size(rho)), sc(size(rho))
    real(dp) :: g(size(rho)), fl(size(rho)), dfl(size(rho))
    integer :: n
    n = size(rho)
    rr = max(rho, 1.0e-30_dp)
    if (.not. f%gga) then
      call c_xc_lda_exc_vxc(f%px, int(n, c_size_t), rr, zx, vx)
      call c_xc_lda_exc_vxc(f%pc, int(n, c_size_t), rr, zc, vc)
      exc = 2.0_dp * (zx + zc); vxc = 2.0_dp * (vx + vc)
      return
    end if
    if (.not. present(r)) error stop 'xc_eval: GGA needs the radial mesh'
    call radial_deriv(r, rho, g)
    sg = g * g
    call c_xc_gga_exc_vxc(f%px, int(n, c_size_t), rr, sg, zx, vx, sx)
    call c_xc_gga_exc_vxc(f%pc, int(n, c_size_t), rr, sg, zc, vc, sc)
    fl = r**2 * 2.0_dp * (sx + sc) * g                       ! r^2 * 2 vsigma * dn/dr
    call radial_deriv(r, fl, dfl)
    exc = 2.0_dp * (zx + zc)
    call core_switch(r, dfl)
    vxc(2:n) = 2.0_dp * (vx(2:n) + vc(2:n) - dfl(2:n) / r(2:n)**2)
    vxc(1) = vxc(2)
  end subroutine xc_eval

  !> Spin-polarised: rho_up, rho_dn (e/bohr^3) -> exc per electron (Ry), vxc_up, vxc_dn (Ry).  r: radial mesh (GGA).
  subroutine xc_eval_pol(f, rup, rdn, exc, vup, vdn, r)
    type(xc_lda_t), intent(in) :: f
    real(dp), intent(in) :: rup(:), rdn(:)
    real(dp), intent(out) :: exc(size(rup)), vup(size(rup)), vdn(size(rup))
    real(dp), intent(in), optional :: r(:)
    real(c_double) :: rr(2 * size(rup)), zx(size(rup)), vx(2 * size(rup)), zc(size(rup)), vc(2 * size(rup))
    real(c_double) :: sg(3 * size(rup)), sx(3 * size(rup)), sc(3 * size(rup))
    real(dp) :: gu(size(rup)), gd(size(rup)), fl(size(rup)), dfl(size(rup)), vs(3, size(rup))
    integer :: i, n
    n = size(rup)
    do i = 1, n
      rr(2 * i - 1) = max(rup(i), 1.0e-30_dp); rr(2 * i) = max(rdn(i), 1.0e-30_dp)
    end do
    if (.not. f%gga) then
      call c_xc_lda_exc_vxc(f%px2, int(n, c_size_t), rr, zx, vx)
      call c_xc_lda_exc_vxc(f%pc2, int(n, c_size_t), rr, zc, vc)
      exc = 2.0_dp * (zx + zc)
      do i = 1, n
        vup(i) = 2.0_dp * (vx(2 * i - 1) + vc(2 * i - 1)); vdn(i) = 2.0_dp * (vx(2 * i) + vc(2 * i))
      end do
      return
    end if
    if (.not. present(r)) error stop 'xc_eval_pol: GGA needs the radial mesh'
    call radial_deriv(r, rup, gu); call radial_deriv(r, rdn, gd)
    do i = 1, n
      sg(3 * i - 2) = gu(i) * gu(i); sg(3 * i - 1) = gu(i) * gd(i); sg(3 * i) = gd(i) * gd(i)
    end do
    call c_xc_gga_exc_vxc(f%px2, int(n, c_size_t), rr, sg, zx, vx, sx)
    call c_xc_gga_exc_vxc(f%pc2, int(n, c_size_t), rr, sg, zc, vc, sc)
    exc = 2.0_dp * (zx + zc)
    do i = 1, n
      vs(:, i) = sx(3 * i - 2:3 * i) + sc(3 * i - 2:3 * i)
    end do
    fl = r**2 * (2.0_dp * vs(1, :) * gu + vs(2, :) * gd)
    call radial_deriv(r, fl, dfl); call core_switch(r, dfl)
    do i = 2, n
      vup(i) = 2.0_dp * (vx(2 * i - 1) + vc(2 * i - 1) - dfl(i) / r(i)**2)
    end do
    fl = r**2 * (2.0_dp * vs(3, :) * gd + vs(2, :) * gu)
    call radial_deriv(r, fl, dfl); call core_switch(r, dfl)
    do i = 2, n
      vdn(i) = 2.0_dp * (vx(2 * i) + vc(2 * i) - dfl(i) / r(i)**2)
    end do
    vup(1) = vup(2); vdn(1) = vdn(2)
  end subroutine xc_eval_pol

  pure subroutine core_switch(r, d)
    real(dp), intent(in) :: r(:)
    real(dp), intent(inout) :: d(:)
    if (xc_gga_rcut > 0.0_dp) d = d * (1.0_dp - exp(-(r / xc_gga_rcut)**2))
  end subroutine core_switch

  !> df/dr on a non-uniform mesh: three-point formula (exact for quadratics), one-sided at the ends
  pure subroutine radial_deriv(r, f, df)
    real(dp), intent(in) :: r(:), f(:)
    real(dp), intent(out) :: df(size(f))
    real(dp) :: h1, h2
    integer :: i, n
    n = size(f)
    do i = 2, n - 1
      h1 = r(i) - r(i - 1); h2 = r(i + 1) - r(i)
      df(i) = (-h2 / (h1 * (h1 + h2))) * f(i - 1) + ((h2 - h1) / (h1 * h2)) * f(i) + (h1 / (h2 * (h1 + h2))) * f(i + 1)
    end do
    h1 = r(2) - r(1); h2 = r(3) - r(2)
    df(1) = (-(2.0_dp * h1 + h2) / (h1 * (h1 + h2))) * f(1) + ((h1 + h2) / (h1 * h2)) * f(2) - (h1 / (h2 * (h1 + h2))) * f(3)
    h1 = r(n - 1) - r(n - 2); h2 = r(n) - r(n - 1)
    df(n) = (h2 / (h1 * (h1 + h2))) * f(n - 2) - ((h1 + h2) / (h1 * h2)) * f(n - 1) + ((h1 + 2.0_dp * h2) / (h2 * (h1 + h2))) * f(n)
  end subroutine radial_deriv

  subroutine xc_finalize(f)
    type(xc_lda_t), intent(inout) :: f
    if (c_associated(f%px)) then; call c_xc_func_end(f%px); call c_xc_func_free(f%px); f%px = c_null_ptr; end if
    if (c_associated(f%pc)) then; call c_xc_func_end(f%pc); call c_xc_func_free(f%pc); f%pc = c_null_ptr; end if
    if (c_associated(f%px2)) then; call c_xc_func_end(f%px2); call c_xc_func_free(f%px2); f%px2 = c_null_ptr; end if
    if (c_associated(f%pc2)) then; call c_xc_func_end(f%pc2); call c_xc_func_free(f%pc2); f%pc2 = c_null_ptr; end if
  end subroutine xc_finalize
end module asr_xc
