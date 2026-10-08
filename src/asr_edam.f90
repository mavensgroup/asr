!> Energy-dependent ASR in the literal AMscf form (AMscf ASR/main.f90 ap1..ap13, ASR/hop_endep.f90).
!>
!> For a central species c (the start site is frozen to c, one recursion per species and orbital)
!>     H(E) = E - T (E - E_nu - h) T - T h o h T ,      T_R = (Delta_c / Delta_R)^(1/2)   (site diagonal, per orbital),
!> h = c + d S d (alpha representation; c = C - E_nu, d = sqrt(Delta)),  and   G_c(E) = < R0 | (E - H(E))^-1 | R0 >.
!> Because T d = sqrt(Delta_c) is a constant the structure constants appear only as the non-random
!> S' = Delta_c^1/2 S Delta_c^1/2 and the disorder sits in the site-diagonal terms:
!>     H(E) = [E - Delta_c (E - E_nu,R)/Delta_R] + h' - h' q h',    h' = g + S',  g_R = (Delta_c/Delta_R)(C_R - E_nu,R),
!>     q_R = o_R Delta_R / Delta_c .
!> This is the H^(2)-type operator of asr_augorb/asr_augsym (c -> g, d -> sqrt(Delta_c), o -> q) with the E_nu term replaced by the
!> complex, energy-dependent diagonal e_R(E) = E - Delta_c (E - E_nu,R)/Delta_R.  The real part R = h' - h' q h' is applied to the
!> real and imaginary parts of the Lanczos vector, e_R(E) = e0 + e1 M by M and a complex-symmetric three-term recursion is run
!> for each E.  At infinite depth E - H(E) = T (E - H^(2)) T, so G_c equals the species-resolved EI Green function.
module asr_edam
  use asr_kinds, only: dp
  use asr_operator, only: linop
  use asr_augorb, only: orb_alloy_t
  use asr_lmto_ham, only: potpar_t, to_alpha, o_alpha
  use asr_endep, only: ed_work_t, ed_depth_tol, ed_nchain, ed_nlev, cf_species
  implicit none
  private
  public :: edam_t, edam_setup, edam_green

  integer, parameter :: orb_l(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]

  !> energy-dependent diagonal e_R(E) = E - dc (E - enu_R)/del_R for the two species (per orbital) and the central species flag
  type :: edam_t
    real(dp) :: dc(9) = 1.0_dp, dela(9) = 1.0_dp, delb(9) = 1.0_dp, enua(9) = 0.0_dp, enub(9) = 0.0_dp
    real(dp) :: mc = 1.0_dp                    !! 1: the central species is A, 0: B
  end type edam_t

contains

  !> Fill the operator parameters of `mdl` (blocks, concentration and node-set data already present) for the central species
  !> (central_is_a) from the gamma-representation parameters of A and B and the screening alpha, and the energy-dependent
  !> diagonal `ed`.  mdl%e0 = mdl%e1 = 0: the E_nu term is applied separately (complex, per energy).
  subroutine edam_setup(mdl, ppa, ppb, alpha, central_is_a, ed)
    type(orb_alloy_t), intent(inout) :: mdl
    type(potpar_t), intent(in) :: ppa, ppb
    real(dp), intent(in) :: alpha(0:2)
    logical, intent(in) :: central_is_a
    type(edam_t), intent(out) :: ed
    type(potpar_t) :: qa, qb, qc
    real(dp) :: oa(0:2), ob(0:2), dc, da, db, ga, gb, wa, wb
    integer :: i, ll
    qa = to_alpha(ppa, alpha); qb = to_alpha(ppb, alpha)
    oa = o_alpha(qa, alpha); ob = o_alpha(qb, alpha)
    qc = qb; if (central_is_a) qc = qa
    if (.not. allocated(mdl%c0)) allocate (mdl%c0(9), mdl%c1(9), mdl%d0(9), mdl%d1(9), mdl%o0(9), mdl%o1(9), mdl%e0(9), mdl%e1(9))
    if (.not. allocated(mdl%p0)) allocate (mdl%p0(9), mdl%p1(9))
    do i = 1, 9
      ll = orb_l(i)
      dc = qc%sqrdel(ll)**2; da = qa%sqrdel(ll)**2; db = qb%sqrdel(ll)**2
      ga = dc / da * (qa%c(ll) - qa%enu(ll)); gb = dc / db * (qb%c(ll) - qb%enu(ll))
      wa = oa(ll) * da / dc; wb = ob(ll) * db / dc
      mdl%c0(i) = gb;       mdl%c1(i) = ga - gb
      mdl%d0(i) = abs(qc%sqrdel(ll)); mdl%d1(i) = 0.0_dp
      mdl%o0(i) = wb;       mdl%o1(i) = wa - wb
      mdl%e0(i) = 0.0_dp;   mdl%e1(i) = 0.0_dp
      mdl%p0(i) = 0.0_dp;   mdl%p1(i) = 0.0_dp
      ed%dc(i) = dc; ed%dela(i) = da; ed%delb(i) = db; ed%enua(i) = qa%enu(ll); ed%enub(i) = qb%enu(ll)
    end do
    ed%mc = merge(1.0_dp, 0.0_dp, central_is_a)
  end subroutine edam_setup

  !> G_c(E) for the start orbital l0 (`iv`: index of the start component in the operator's vector, 9 (node-1) + l0).
  !> The operator `op` carries the parameters set by edam_setup, the start site frozen to the central species, and `mapply`.
  !> wk must have been initialised with ed_work_init(wk, op%dim()).
  subroutine edam_green(op, wk, ed, z, iv, nlev, g, acoef, bcoef, mused)
    class(linop), intent(in) :: op
    type(ed_work_t), intent(inout) :: wk
    type(edam_t), intent(in) :: ed
    complex(dp), intent(in) :: z
    integer, intent(in) :: iv, nlev
    complex(dp), intent(out) :: g
    complex(dp), intent(out), optional :: acoef(:), bcoef(:)   !! chain coefficients (diagnostics: exact continued fraction)
    integer, intent(out), optional :: mused
    complex(dp) :: e0(9), e1(9), a(nlev), b(nlev), cc(nlev), gtry, gprev, ftmp, fa, ea, eb
    integer :: k, l, i, m, nn
    nn = op%dim() / 9
    do l = 1, 9
      ea = z - ed%dc(l) * (z - ed%enua(l)) / ed%dela(l)
      eb = z - ed%dc(l) * (z - ed%enub(l)) / ed%delb(l)
      e0(l) = eb; e1(l) = ea - eb
    end do
    gprev = (0.0_dp, 0.0_dp)
    wk%u = (0.0_dp, 0.0_dp); wk%up = (0.0_dp, 0.0_dp)
    wk%u(iv) = (1.0_dp, 0.0_dp)
    a = 0.0_dp; b = 0.0_dp; cc = 0.0_dp; m = 0
    do k = 1, nlev
      wk%re = real(wk%u); wk%im = aimag(wk%u)
      call op%apply(wk%re, wk%sre)                       ! R Re u
      call op%apply(wk%im, wk%sim)                       ! R Im u
      call op%mapply(wk%re, wk%mre)
      call op%mapply(wk%im, wk%mim)
      !$omp parallel do private(l)
      do i = 1, nn
        do l = 1, 9
          wk%w(9 * (i - 1) + l) = cmplx(wk%sre(9 * (i - 1) + l), wk%sim(9 * (i - 1) + l), dp) &
                                  + e0(l) * wk%u(9 * (i - 1) + l) &
                                  + e1(l) * cmplx(wk%mre(9 * (i - 1) + l), wk%mim(9 * (i - 1) + l), dp)
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
  end subroutine edam_green
end module asr_edam
