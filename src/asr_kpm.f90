!> Kernel polynomial method (Chebyshev) density of states with stochastic trace and the Jackson kernel.
!> Independent of the continued-fraction machinery: used as ground truth for large disordered systems.
module asr_kpm
  use asr_kinds, only: dp, pi
  use asr_operator, only: linop
  implicit none
  private
  public :: kpm_moments, kpm_dos, lorentz_broaden, spectral_bounds

contains

  !> Extremal eigenvalues of op (estimated by nstep Lanczos steps from a random vector), widened by `margin`.
  subroutine spectral_bounds(op, nstep, margin, emin, emax)
    class(linop), intent(in) :: op
    integer, intent(in) :: nstep
    real(dp), intent(in) :: margin
    real(dp), intent(out) :: emin, emax
    real(dp), allocatable :: v(:), u(:), up(:), w(:)
    real(dp) :: a(nstep), b(nstep), d(nstep), e(nstep), bn
    integer :: i, dim, info
    external :: dsterf
    dim = op%dim()
    allocate (v(dim), u(dim), up(dim), w(dim))
    call random_number(v); v = v - 0.5_dp
    u = v / sqrt(dot_product(v, v)); up = 0.0_dp; bn = 0.0_dp
    do i = 1, nstep
      call op%apply(u, w)
      a(i) = dot_product(u, w)
      w = w - a(i) * u - bn * up
      bn = sqrt(dot_product(w, w)); b(i) = bn
      up = u; u = w / max(bn, 1.0e-300_dp)
    end do
    d = a; e = 0.0_dp; e(1:nstep - 1) = b(1:nstep - 1)
    call dsterf(nstep, d, e, info)
    emin = d(1) - margin; emax = d(nstep) + margin
  end subroutine spectral_bounds

  !> mu_n = (1/R) sum_r <r|T_n(Hs)|r>, Hs = (H - c)/s with spectrum inside (-1,1);  n = 0..nmom-1.
  subroutine kpm_moments(op, c, s, nmom, nvec, seed, mu)
    class(linop), intent(in) :: op
    real(dp), intent(in) :: c, s
    integer, intent(in) :: nmom, nvec, seed
    real(dp), intent(out) :: mu(0:nmom - 1)
    real(dp), allocatable :: r(:), t0(:), t1(:), t2(:), w(:)
    integer :: n, iv, k, dim
    integer, allocatable :: sd(:)
    dim = op%dim()
    allocate (r(dim), t0(dim), t1(dim), t2(dim), w(dim))
    call random_seed(size=k); allocate (sd(k)); sd = seed + 101 * [(iv, iv=1, k)]; call random_seed(put=sd)
    mu = 0.0_dp
    do iv = 1, nvec
      call random_number(r); r = merge(1.0_dp, -1.0_dp, r > 0.5_dp)         ! Rademacher random vector
      t0 = r
      call op%apply(r, w); t1 = (w - c * r) / s
      mu(0) = mu(0) + dot_product(r, t0) / dim
      if (nmom > 1) mu(1) = mu(1) + dot_product(r, t1) / dim
      do n = 2, nmom - 1
        call op%apply(t1, w)
        t2 = 2.0_dp * (w - c * t1) / s - t0
        mu(n) = mu(n) + dot_product(r, t2) / dim
        t0 = t1; t1 = t2
      end do
    end do
    mu = mu / real(nvec, dp)
  end subroutine kpm_moments

  !> DOS(E) (states per site per unit energy) on the grid e, Jackson kernel.
  subroutine kpm_dos(mu, c, s, e, rho)
    real(dp), intent(in) :: mu(0:), c, s, e(:)
    real(dp), intent(out) :: rho(size(e))
    integer :: nmom, n, i
    real(dp) :: g(0:size(mu) - 1), x, th, sum_
    nmom = size(mu)
    do n = 0, nmom - 1
      g(n) = ((nmom - n + 1) * cos(pi * n / (nmom + 1)) + sin(pi * n / (nmom + 1)) / tan(pi / (nmom + 1))) / (nmom + 1)
    end do
    do i = 1, size(e)
      x = (e(i) - c) / s
      if (abs(x) >= 1.0_dp) then; rho(i) = 0.0_dp; cycle; end if
      th = acos(x)
      sum_ = g(0) * mu(0)
      do n = 1, nmom - 1
        sum_ = sum_ + 2.0_dp * g(n) * mu(n) * cos(n * th)
      end do
      rho(i) = sum_ / (pi * sqrt(1.0_dp - x * x) * s)
    end do
  end subroutine kpm_dos

  !> Convolve rho(e) (uniform grid) with a Lorentzian of half-width eta: the DOS seen at Im z = eta.
  subroutine lorentz_broaden(e, rho, eta, rhoeta)
    real(dp), intent(in) :: e(:), rho(:), eta
    real(dp), intent(out) :: rhoeta(size(e))
    integer :: i, j
    real(dp) :: de
    de = e(2) - e(1)
    do i = 1, size(e)
      rhoeta(i) = 0.0_dp
      do j = 1, size(e)
        rhoeta(i) = rhoeta(i) + rho(j) * (eta / pi) / ((e(i) - e(j))**2 + eta**2) * de
      end do
    end do
  end subroutine lorentz_broaden
end module asr_kpm
