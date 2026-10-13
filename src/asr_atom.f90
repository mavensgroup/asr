!> Free-atom self-consistent LDA (all electrons, scalar-relativistic or non-relativistic).
!> Spherical, non-spin-polarised; Rydberg units.  Used to validate the radial/Poisson/XC stack
!> and to seed crystal calculations.
module asr_atom
  use asr_kinds, only: dp, pi
  use asr_radial
  use asr_xc
  implicit none
  private
  public :: shell_t, atom_result_t, solve_atom, cumint, poisson_hartree

  type :: shell_t
    integer :: n, l
    real(dp) :: occ
  end type shell_t

  type :: atom_result_t
    type(radmesh_t) :: mesh
    real(dp), allocatable :: v(:), rho(:), nr(:)     !! V_true(r), rho(r) (e/bohr^3), nr = 4 pi r^2 rho
    real(dp), allocatable :: eps(:)
    real(dp) :: etot = 0.0_dp                        !! Ry
    integer :: iter = 0
    logical :: converged = .false.
  end type atom_result_t

contains

  !> F(i) = int_{r(1)}^{r(i)} f dr, f already given as f(r) (not multiplied by dr/ds), 4th order.
  pure subroutine cumint(m, f, cum)
    type(radmesh_t), intent(in) :: m
    real(dp), intent(in) :: f(:)
    real(dp), intent(out) :: cum(size(f))
    real(dp) :: g(size(f))
    integer :: i, n
    n = size(f)
    g = f * m%drdt
    cum(1) = 0.0_dp
    do i = 1, n - 1
      if (i >= 2 .and. i + 2 <= n) then
        cum(i + 1) = cum(i) + (-g(i - 1) + 13.0_dp * g(i) + 13.0_dp * g(i + 1) - g(i + 2)) / 24.0_dp
      else if (i + 3 <= n) then
        cum(i + 1) = cum(i) + (9.0_dp * g(i) + 19.0_dp * g(i + 1) - 5.0_dp * g(i + 2) + g(i + 3)) / 24.0_dp
      else
        cum(i + 1) = cum(i) + (g(i - 2) - 5.0_dp * g(i - 1) + 19.0_dp * g(i) + 9.0_dp * g(i + 1)) / 24.0_dp
      end if
    end do
  end subroutine cumint

  !> Hartree potential (Ry) of a spherical density with radial charge n(r) = 4 pi r^2 rho:
  !> V_H(r) = 2 [ Q(r)/r + int_r^inf n(r')/r' dr' ]
  subroutine poisson_hartree(m, nr, vh, qtot)
    type(radmesh_t), intent(in) :: m
    real(dp), intent(in) :: nr(:)
    real(dp), intent(out) :: vh(size(nr))
    real(dp), intent(out), optional :: qtot
    real(dp) :: q(size(nr)), w(size(nr)), wc(size(nr))
    integer :: i, n
    n = size(nr)
    call cumint(m, nr, q)
    w = 0.0_dp
    do i = 2, n; w(i) = nr(i) / m%r(i); end do
    w(1) = w(2)
    call cumint(m, w, wc)
    vh(1) = 2.0_dp * wc(n)
    do i = 2, n
      vh(i) = 2.0_dp * (q(i) / m%r(i) + (wc(n) - wc(i)))
    end do
    if (present(qtot)) qtot = q(n)
  end subroutine poisson_hartree

  subroutine solve_atom(z, shells, xcname, rmax, a, relativistic, res, maxit, tol, mix)
    real(dp), intent(in) :: z, rmax, a
    type(shell_t), intent(in) :: shells(:)
    character(len=*), intent(in) :: xcname
    logical, intent(in) :: relativistic
    type(atom_result_t), intent(out) :: res
    integer, intent(in), optional :: maxit
    real(dp), intent(in), optional :: tol, mix
    type(xc_lda_t) :: xf
    real(dp), allocatable :: v(:), vnew(:), vh(:), exc(:), vxc(:), rho(:), nrr(:), y(:, :, :), tmp(:), cum(:)
    real(dp) :: c, beta, dv, eps, e1, e2, etot, ehar, exc_int, vxc_int, sumeps
    integer :: it, nit, is, n, isph, i, nod
    real(dp) :: tl

    nit = 200; if (present(maxit)) nit = maxit
    tl = 1.0e-9_dp; if (present(tol)) tl = tol
    beta = 0.3_dp; if (present(mix)) beta = mix
    c = c_nonrel; if (relativistic) c = c_ry
    call make_radmesh(z, rmax, a, res%mesh, c)
    n = res%mesh%n
    isph = n; if (mod(isph, 2) == 0) isph = isph - 1
    allocate (v(n), vnew(n), vh(n), exc(n), vxc(n), rho(n), nrr(n), tmp(n), cum(n), res%eps(size(shells)))
    allocate (y(2, 0:0, isph))
    call xc_init(xf, xcname)
    ! start: Thomas-Fermi-like screened Coulomb
    v(2:) = -2.0_dp * z / res%mesh%r(2:) * exp(-res%mesh%r(2:) / (0.6_dp * z**(-1.0_dp / 3.0_dp) * 1.0_dp * 3.0_dp))
    v(1) = v(2)
    do i = 2, n
      v(i) = -2.0_dp * z / res%mesh%r(i) * (1.0_dp + 0.0_dp)
      v(i) = v(i) * exp(-res%mesh%r(i) * z**(1.0_dp / 3.0_dp) / 4.0_dp) - 1.0e-3_dp
    end do
    v(1) = v(2)
    res%converged = .false.
    do it = 1, nit
      nrr = 0.0_dp; sumeps = 0.0_dp
      do is = 1, size(shells)
        nod = shells(is)%n - shells(is)%l - 1
        call bound_state(res%mesh, v, shells(is)%l, nod, z, c, isph, -1.6_dp * z * z - 20.0_dp, -1.0e-4_dp, eps, y)
        res%eps(is) = eps
        sumeps = sumeps + shells(is)%occ * eps
        nrr(1:isph) = nrr(1:isph) + shells(is)%occ * (y(1, 0, 1:isph)**2 + y(2, 0, 1:isph)**2 / c**2)
      end do
      ! charge density and potentials from the output density
      rho = 0.0_dp
      do i = 2, n; rho(i) = nrr(i) / (4.0_dp * pi * res%mesh%r(i)**2); end do
      rho(1) = rho(2)
      call poisson_hartree(res%mesh, nrr, vh)
      call xc_eval(xf, rho, exc, vxc, res%mesh%r)
      vnew(2:) = -2.0_dp * z / res%mesh%r(2:) + vh(2:) + vxc(2:)
      vnew(1) = vnew(2)
      dv = maxval(abs(vnew(2:n) - v(2:n)) * res%mesh%r(2:n) / (1.0_dp + res%mesh%r(2:n)))
      ! total energy (Ry):  E = sum eps - 1/2 int n V_H - int n vxc + int n exc
      tmp = nrr * vh;  call cumint(res%mesh, tmp, cum); ehar = 0.5_dp * cum(n)
      tmp = nrr * vxc; call cumint(res%mesh, tmp, cum); vxc_int = cum(n)
      tmp = nrr * exc; call cumint(res%mesh, tmp, cum); exc_int = cum(n)
      etot = sumeps - ehar - vxc_int + exc_int
      res%iter = it
      if (dv < tl) then
        res%converged = .true.
        v = vnew
        exit
      end if
      v = (1.0_dp - beta) * v + beta * vnew
    end do
    call xc_finalize(xf)
    res%etot = etot
    allocate (res%v(n), res%rho(n), res%nr(n))
    res%v = v; res%rho = rho; res%nr = nrr
    e1 = 0.0_dp; e2 = 0.0_dp
  end subroutine solve_atom
end module asr_atom
