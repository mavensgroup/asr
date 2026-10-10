!> n-species energy-dependent ASR (AMscf form, see asr_ed) on the n-component augmented space asr_augn: for a central species c
!> (start site frozen to c), H(E) = [e_R(E)] + h' - h' q h',  h' = g + S',  S' = D_c^1/2 S D_c^1/2,
!>   g_s = (D_c/D_s)(C_s - E_nu,s),  q_s = o_s D_s/D_c,  e_s(E) = E - D_c (E - E_nu,s)/D_s   (s = species of the random site),
!> and G_c(E) = <1|(E - H(E))^-1|1> from a complex-symmetric Lanczos chain per energy; one recursion per central species and orbital.
module asr_edn
  use asr_kinds, only: dp
  use asr_augn, only: multi_alloy_t, augn_op, augn_rand_c
  use asr_lmto_ham, only: potpar_t, to_alpha, o_alpha
  use asr_endep, only: ed_depth_tol, ed_nchain, ed_nlev, cf_species
  implicit none
  private
  public :: edn_t, edn_work_t, edn_setup, edn_green, edn_work_init

  integer, parameter :: orb_l(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]

  type :: edn_t
    integer :: ns = 2, central = 1
    real(dp), allocatable :: dc(:), del(:, :), enu(:, :)    !! D_c (9); D_s, E_nu,s (9, ns)
  end type edn_t

  type :: edn_work_t
    real(dp), allocatable :: re(:), im(:), sre(:), sim(:)
    complex(dp), allocatable :: u(:), up(:), w(:), dg(:)
  end type edn_work_t

contains

  subroutine edn_work_init(wk, n)
    type(edn_work_t), intent(out) :: wk
    integer, intent(in) :: n
    allocate (wk%re(n), wk%im(n), wk%sre(n), wk%sim(n), wk%u(n), wk%up(n), wk%w(n), wk%dg(n))
  end subroutine edn_work_init

  !> Operator parameters of `mdl` (9 orbitals, ns species; blocks and concentrations already present) for the central species
  !> `central` from the gamma-representation parameters pp(1:ns) and the screening alpha.  mdl%e = 0: the E_nu term is applied
  !> separately (complex, per energy).
  subroutine edn_setup(mdl, pp, alpha, central, ed)
    type(multi_alloy_t), intent(inout) :: mdl
    type(potpar_t), intent(in) :: pp(:)
    real(dp), intent(in) :: alpha(0:2)
    integer, intent(in) :: central
    type(edn_t), intent(out) :: ed
    type(potpar_t), allocatable :: q(:)
    real(dp) :: o(0:2), dc, ds
    integer :: ns, s, i, ll
    ns = size(pp)
    allocate (q(ns)); ed%ns = ns; ed%central = central
    allocate (ed%dc(9), ed%del(9, ns), ed%enu(9, ns))
    if (.not. allocated(mdl%c)) allocate (mdl%c(9, ns), mdl%d(9, ns), mdl%o(9, ns), mdl%e(9, ns))
    do s = 1, ns; q(s) = to_alpha(pp(s), alpha); end do
    do i = 1, 9
      ll = orb_l(i)
      dc = q(central)%sqrdel(ll)**2
      ed%dc(i) = dc
      do s = 1, ns
        o = o_alpha(q(s), alpha)
        ds = q(s)%sqrdel(ll)**2
        mdl%c(i, s) = dc / ds * (q(s)%c(ll) - q(s)%enu(ll))
        mdl%d(i, s) = abs(q(central)%sqrdel(ll))
        mdl%o(i, s) = o(ll) * ds / dc
        mdl%e(i, s) = 0.0_dp
        ed%del(i, s) = ds; ed%enu(i, s) = q(s)%enu(ll)
      end do
    end do
  end subroutine edn_setup

  !> G_c(E) for the start component iv (9 (node-1) + l0); op carries edn_setup's parameters with the start site frozen to the
  !> central species.  Optionally returns the chain coefficients.
  subroutine edn_green(op, wk, ed, z, iv, nlev, g, acoef, bcoef, mused)
    type(augn_op), intent(in) :: op
    type(edn_work_t), intent(inout) :: wk
    type(edn_t), intent(in) :: ed
    complex(dp), intent(in) :: z
    integer, intent(in) :: iv, nlev
    complex(dp), intent(out) :: g
    complex(dp), intent(out), optional :: acoef(:), bcoef(:)
    integer, intent(out), optional :: mused
    complex(dp) :: ep(9, ed%ns), a(nlev), b(nlev), cc(nlev), gtry, gprev, ftmp, fa
    integer :: k, l, i, m, nn, s
    nn = op%nn
    do s = 1, ed%ns
      do l = 1, 9
        ep(l, s) = z - ed%dc(l) * (z - ed%enu(l, s)) / ed%del(l, s)
      end do
    end do
    gprev = (0.0_dp, 0.0_dp)
    wk%u = (0.0_dp, 0.0_dp); wk%up = (0.0_dp, 0.0_dp)
    wk%u(iv) = (1.0_dp, 0.0_dp)
    a = 0.0_dp; b = 0.0_dp; cc = 0.0_dp; m = 0
    do k = 1, nlev
      wk%re = real(wk%u); wk%im = aimag(wk%u)
      call op%apply(wk%re, wk%sre)
      call op%apply(wk%im, wk%sim)
      call augn_rand_c(op, ep, wk%u, wk%dg)
      !$omp parallel do private(l)
      do i = 1, nn
        do l = 1, 9
          wk%w(9 * (i - 1) + l) = cmplx(wk%sre(9 * (i - 1) + l), wk%sim(9 * (i - 1) + l), dp) + wk%dg(9 * (i - 1) + l)
        end do
      end do
      !$omp end parallel do
      a(k) = sum(wk%u * wk%w)
      m = k
      if (k == nlev) exit
      if (ed_depth_tol > 0.0_dp .and. k >= 8 .and. mod(k, 4) == 0) then
        call cf_species(a, b, cc, k, z, gtry, ftmp)
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
    call cf_species(a, b, cc, m, z, g, fa)
    if (present(acoef)) acoef(1:m) = a(1:m)
    if (present(bcoef)) bcoef(1:m) = b(1:m)
    if (present(mused)) mused = m
  end subroutine edn_green
end module asr_edn
