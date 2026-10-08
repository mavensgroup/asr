!> Energy-dependent augmented-space recursion (potential-function form of TB-LMTO-ASA).
!>
!> For every complex energy z the recursion runs on  A(z) = P^alpha(z) - S^alpha  with the random diagonal potential
!> function  P^alpha(z) = (z - C)/(Delta + (gamma - alpha)(z - C))  (C, Delta, gamma: gamma-representation parameters of
!> the species on the site).  The physical (orthonormal, gamma-representation) Green function is
!>      G = lambda + mu g^alpha mu,    mu = sqrt(Delta)/(Delta + delta (z - C)),   lambda = delta/(Delta + delta (z - C)),
!> delta = gamma - alpha  (site diagonal), derived from  g^gamma = (1 - delta P) g^alpha (1 - P delta) + delta (1 - P delta).
!> Because n_R is a projector, n_R mu_R^2 = mu_A^2 n_R, so the species-resolved functions need only the average g_RR and the
!> transition element <n_R g_RR> from the complex-symmetric Lanczos chain of A(z) started on |R,{},L>.
module asr_endep
  use asr_kinds, only: dp
  use asr_augsym, only: augsym_op, augsym_ops
  use asr_lmto_realization, only: lmto_real_op, lmto_s_apply
  implicit none
  private
  public :: ed_params_t, ed_work_t, ed_functions, ed_species_green, ed_work_init, ed_ordered_green
  public :: ed_depth_tol, ed_nchain, ed_nlev, cf_species

  !> adaptive depth: from level 8, every 4 levels, stop when the terminated Green function changes by less than ed_depth_tol
  !> (relative); 0 = always use the full depth.  ed_nchain / ed_nlev count chains and levels used (statistics).
  real(dp) :: ed_depth_tol = 0.0_dp
  integer(8) :: ed_nchain = 0, ed_nlev = 0

  !> gamma-representation parameters (9 orbitals) of one species / spin; p = 0: two-parameter P^gamma = (z - C)/Delta,
  !> p > 0: third-order P^gamma = (z - C + p (z - E_nu)^3)/Delta (p = <phidot|phidot>, exact to third order in z - E_nu)
  type :: ed_params_t
    real(dp) :: c(9) = 0.0_dp, del(9) = 0.0_dp, gam(9) = 0.0_dp, p(9) = 0.0_dp, enu(9) = 0.0_dp
    logical :: exact = .false.                   !! use pfx, mux, lmx (exact P^alpha, mu^2, lambda at the current energy, asr_pfun)
    complex(dp) :: pfx(9) = (0.0_dp, 0.0_dp), mux(9) = (0.0_dp, 0.0_dp), lmx(9) = (0.0_dp, 0.0_dp)
  end type ed_params_t

  type :: ed_work_t
    real(dp), allocatable :: re(:), im(:), mre(:), mim(:), sre(:), sim(:)
    complex(dp), allocatable :: u(:), up(:), w(:)
  end type ed_work_t

contains

  subroutine ed_work_init(wk, n)
    type(ed_work_t), intent(out) :: wk
    integer, intent(in) :: n
    allocate (wk%re(n), wk%im(n), wk%mre(n), wk%mim(n), wk%sre(n), wk%sim(n), wk%u(n), wk%up(n), wk%w(n))
  end subroutine ed_work_init

  !> P, mu^2, lambda for the 9 orbitals of a species at energy z (alpha: screening per orbital)
  pure subroutine ed_functions(z, sp, alpha, pf, mu2, lam)
    complex(dp), intent(in) :: z
    type(ed_params_t), intent(in) :: sp
    real(dp), intent(in) :: alpha(9)
    complex(dp), intent(out) :: pf(9), mu2(9), lam(9)
    complex(dp) :: den, nn, nd, ndd, eps
    real(dp) :: dl
    integer :: l
    if (sp%exact) then
      pf = sp%pfx; mu2 = sp%mux; lam = sp%lmx
      return
    end if
    ! P^gamma = N/Delta, N = z - C + p eps^3;  P^alpha = N / (Delta + (gamma - alpha) N);
    ! mu^2 = dP^alpha/dz = N' Delta / den^2,  lambda = -P''/(2 P') = (gamma - alpha) N'/den - N''/(2 N')
    do l = 1, 9
      dl = sp%gam(l) - alpha(l)
      eps = z - sp%enu(l)
      nn = z - sp%c(l) + sp%p(l) * eps**3
      nd = 1.0_dp + 3.0_dp * sp%p(l) * eps**2
      ndd = 6.0_dp * sp%p(l) * eps
      den = sp%del(l) + dl * nn
      pf(l) = nn / den
      mu2(l) = nd * sp%del(l) / den**2
      lam(l) = dl * nd / den - 0.5_dp * ndd / nd
    end do
  end subroutine ed_functions

  !> Species-resolved physical Green functions at energy z for the start orbital l0 (the operator's start orbital
  !> must have been set with augsym_set_start): GA, GB such that  x GA + (1-x) GB = <G_RR>.
  subroutine ed_species_green(op, wk, z, spa, spb, alpha, x, l0, nlev, flip1, md1, sflip, ga, gb, gavg_out)
    type(augsym_op), intent(in) :: op
    type(ed_work_t), intent(inout) :: wk
    complex(dp), intent(in) :: z
    type(ed_params_t), intent(in) :: spa, spb
    real(dp), intent(in) :: alpha(9), x, md1, sflip
    integer, intent(in) :: l0, nlev, flip1
    complex(dp), intent(out) :: ga, gb
    complex(dp), intent(out), optional :: gavg_out
    complex(dp) :: pa(9), pb(9), mua(9), mub(9), lama(9), lamb(9), p0(9), p1(9), a(nlev), b(nlev), cc(nlev)
    complex(dp) :: gavg, fa, gtry, ftry, gprev, fprev
    integer :: k, l, i, nd, m, nn
    nd = size(wk%u); nn = nd / 9
    gprev = (0.0_dp, 0.0_dp); fprev = (0.0_dp, 0.0_dp)
    call ed_functions(z, spa, alpha, pa, mua, lama)
    call ed_functions(z, spb, alpha, pb, mub, lamb)
    p0 = pb; p1 = pa - pb
    wk%u = (0.0_dp, 0.0_dp); wk%up = (0.0_dp, 0.0_dp)
    wk%u(l0) = (1.0_dp, 0.0_dp)
    a = 0.0_dp; b = 0.0_dp; cc = 0.0_dp; m = 0
    do k = 1, nlev
      cc(k) = md1 * wk%u(l0)
      if (flip1 > 0) cc(k) = cc(k) + sflip * wk%u(9 * (flip1 - 1) + l0)
      wk%re = real(wk%u); wk%im = aimag(wk%u)
      call augsym_ops(op, wk%re, wk%mre, wk%sre)
      call augsym_ops(op, wk%im, wk%mim, wk%sim)
      !$omp parallel do private(l)
      do i = 1, nn
        do l = 1, 9
          wk%w(9 * (i - 1) + l) = p0(l) * wk%u(9 * (i - 1) + l) + p1(l) * cmplx(wk%mre(9 * (i - 1) + l), wk%mim(9 * (i - 1) + l), dp) &
                                  - cmplx(wk%sre(9 * (i - 1) + l), wk%sim(9 * (i - 1) + l), dp)
        end do
      end do
      !$omp end parallel do
      a(k) = sum(wk%u * wk%w)
      m = k
      if (k == nlev) exit
      if (ed_depth_tol > 0.0_dp .and. k >= 8 .and. mod(k, 4) == 0) then
        call cf_species(-a, -b, cc, k, (0.0_dp, 0.0_dp), gtry, ftry)
        if (k > 8) then
          if (abs(gtry - gprev) <= ed_depth_tol * abs(gtry) .and. abs(ftry - fprev) <= ed_depth_tol * max(abs(ftry), abs(gtry))) exit
        end if
        gprev = gtry; fprev = ftry
      end if
      wk%w = wk%w - a(k) * wk%u
      if (k > 1) wk%w = wk%w - b(k - 1) * wk%up
      b(k) = sqrt(sum(wk%w * wk%w))
      if (abs(b(k)) < 1.0e-12_dp) exit
      wk%up = wk%u
      wk%u = wk%w / b(k)
    end do
    ed_nchain = ed_nchain + 1; ed_nlev = ed_nlev + m
    call cf_species(-a, -b, cc, m, (0.0_dp, 0.0_dp), gavg, fa)      ! A^-1 = (0 - T)^-1 with T = -A
    if (present(gavg_out)) gavg_out = gavg
    ga = lama(l0) + mua(l0) * fa / x
    gb = lamb(l0) + mub(l0) * (gavg - fa) / (1.0_dp - x)
  end subroutine ed_species_green

  !> Energy-dependent Green function G_RR(z) of ONE site of an ordered crystal (no augmented space): complex-symmetric
  !> Lanczos on A(z) = P_site(z) - S started on |site, l0>, G = lambda + mu^2 [A^-1]_11 with the parameters of the site's
  !> species (op%isa: .true. = species a).  wk must have been initialised with ed_work_init(wk, 9 n).
  subroutine ed_ordered_green(op, wk, z, spa, spb, alpha, site, l0, nlev, g)
    type(lmto_real_op), intent(in) :: op
    type(ed_work_t), intent(inout) :: wk
    complex(dp), intent(in) :: z
    type(ed_params_t), intent(in) :: spa, spb
    real(dp), intent(in) :: alpha(9)
    integer, intent(in) :: site, l0, nlev
    complex(dp), intent(out) :: g
    complex(dp) :: pa(9), pb(9), mua(9), mub(9), lama(9), lamb(9), a(nlev), b(nlev), cc(nlev), gavg, fa, gtry, gprev
    integer :: k, l, i, m, nn
    nn = op%n
    gprev = (0.0_dp, 0.0_dp)
    call ed_functions(z, spa, alpha, pa, mua, lama)
    call ed_functions(z, spb, alpha, pb, mub, lamb)
    wk%u = (0.0_dp, 0.0_dp); wk%up = (0.0_dp, 0.0_dp)
    wk%u(9 * (site - 1) + l0) = (1.0_dp, 0.0_dp)
    a = 0.0_dp; b = 0.0_dp; cc = 0.0_dp; m = 0
    do k = 1, nlev
      call lmto_s_apply(op, wk%u, wk%w)
      !$omp parallel do private(l)
      do i = 1, nn
        do l = 1, 9
          if (op%isa(i)) then
            wk%w(9 * (i - 1) + l) = pa(l) * wk%u(9 * (i - 1) + l) - wk%w(9 * (i - 1) + l)
          else
            wk%w(9 * (i - 1) + l) = pb(l) * wk%u(9 * (i - 1) + l) - wk%w(9 * (i - 1) + l)
          end if
        end do
      end do
      !$omp end parallel do
      a(k) = sum(wk%u * wk%w)
      m = k
      if (k == nlev) exit
      if (ed_depth_tol > 0.0_dp .and. k >= 8 .and. mod(k, 4) == 0) then
        call cf_species(-a, -b, cc, k, (0.0_dp, 0.0_dp), gtry, fa)
        if (k > 8) then
          if (abs(gtry - gprev) <= ed_depth_tol * abs(gtry)) exit
        end if
        gprev = gtry
      end if
      wk%w = wk%w - a(k) * wk%u
      if (k > 1) wk%w = wk%w - b(k - 1) * wk%up
      b(k) = sqrt(sum(wk%w * wk%w))
      if (abs(b(k)) < 1.0e-12_dp) exit
      wk%up = wk%u
      wk%u = wk%w / b(k)
    end do
    ed_nchain = ed_nchain + 1; ed_nlev = ed_nlev + m
    call cf_species(-a, -b, cc, m, (0.0_dp, 0.0_dp), gavg, fa)
    if (op%isa(site)) then
      g = lama(l0) + mua(l0) * gavg
    else
      g = lamb(l0) + mub(l0) * gavg
    end if
  end subroutine ed_ordered_green

  !> g_11 and sum_k c_k [g]_{k1} of the complex-symmetric tridiagonal chain (a, b) with spectral variable z, terminated by a
  !> square-root tail whose band edges are the extreme eigenvalues (real parts) of the finite chain (Beer-Pettifor type): the
  !> tail band then contains the spectral variable wherever the energy lies inside the band, giving a continuum instead of
  !> m isolated poles.
  subroutine cf_species(a, b, c, m, z, g11, f)
    complex(dp), intent(in) :: a(:), b(:), c(:), z
    integer, intent(in) :: m
    complex(dp), intent(out) :: g11, f
    complex(dp) :: gt(m + 1), prod, tail, bb, d, s, g1, g2, aa
    complex(dp) :: tm(m, m), ev(m), vl(1, 1), vr(1, 1), wrk(8 * m + 8)
    real(dp) :: rw(2 * m)
    integer :: j, k, info, ilo, ihi
    external :: zgeev
    if (m >= 2) then
      tm = (0.0_dp, 0.0_dp)
      do j = 1, m
        tm(j, j) = a(j)
        if (j < m) then; tm(j, j + 1) = b(j); tm(j + 1, j) = b(j); end if
      end do
      call zgeev('N', 'N', m, tm, m, ev, vl, 1, vr, 1, wrk, size(wrk), rw, info)
    else
      info = 1
    end if
    if (info == 0) then
      ilo = minloc(real(ev), 1); ihi = maxloc(real(ev), 1)
      aa = 0.5_dp * (ev(ilo) + ev(ihi)); bb = 0.25_dp * (ev(ihi) - ev(ilo))
    else
      aa = a(m); bb = (0.0_dp, 0.0_dp); if (m > 1) bb = b(m - 1)
    end if
    d = z - aa
    if (abs(bb) < 1.0e-14_dp) then
      tail = (0.0_dp, 0.0_dp)
    else
      s = sqrt(d * d - 4.0_dp * bb * bb)
      g1 = (d - s) / (2.0_dp * bb * bb); g2 = (d + s) / (2.0_dp * bb * bb)
      if (abs(bb * g1) <= abs(bb * g2)) then; tail = bb * bb * g1; else; tail = bb * bb * g2; end if
    end if
    gt(m + 1) = (0.0_dp, 0.0_dp)
    gt(m) = 1.0_dp / (z - a(m) - tail)
    do j = m - 1, 1, -1
      gt(j) = 1.0_dp / (z - a(j) - b(j)**2 * gt(j + 1))
    end do
    g11 = gt(1)
    f = (0.0_dp, 0.0_dp); prod = gt(1)
    do k = 1, m
      f = f + c(k) * prod
      if (k < m) prod = prod * b(k) * gt(k + 1)
    end do
  end subroutine cf_species
end module asr_endep
