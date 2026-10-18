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
module asr_ed
  use asr_kinds, only: dp
  use asr_operator, only: linop
  use asr_augorb, only: orb_alloy_t
  use asr_lmto_ham, only: potpar_t, to_alpha, o_alpha
  use asr_endep, only: ed_work_t, ed_depth_tol, ed_nchain, ed_nlev, cf_species
  implicit none
  private
  public :: asr_ed_t, asr_ed_setup, asr_ed_green, asr_ed_medium, asr_ed_medium_k

  integer, parameter :: orb_l(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]

  !> energy-dependent diagonal e_R(E) = E - dc (E - enu_R)/del_R for the two species (per orbital) and the central species flag
  type :: asr_ed_t
    real(dp) :: dc(9) = 1.0_dp, dela(9) = 1.0_dp, delb(9) = 1.0_dp, enua(9) = 0.0_dp, enub(9) = 0.0_dp
    real(dp) :: mc = 1.0_dp                    !! 1: the central species is A, 0: B
  end type asr_ed_t

contains

  !> Fill the operator parameters of `mdl` (blocks, concentration and node-set data already present) for the central species
  !> (central_is_a) from the gamma-representation parameters of A and B and the screening alpha, and the energy-dependent
  !> diagonal `ed`.  mdl%e0 = mdl%e1 = 0: the E_nu term is applied separately (complex, per energy).
  subroutine asr_ed_setup(mdl, ppa, ppb, alpha, central_is_a, ed)
    type(orb_alloy_t), intent(inout) :: mdl
    type(potpar_t), intent(in) :: ppa, ppb
    real(dp), intent(in) :: alpha(0:2)
    logical, intent(in) :: central_is_a
    type(asr_ed_t), intent(out) :: ed
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
  end subroutine asr_ed_setup

  !> Average medium (VCA) of the operator set by asr_ed_setup for the same central species: the random site-diagonal terms g, q
  !> and e_R(E) replaced by their configuration averages (concentration x of A), the hops sqrt(Delta_c) S sqrt(Delta_c) unchanged.
  !> e_R(E) = E - Delta_c (E - E_nu,R)/Delta_R averages to the same form with 1/Delta_v = x/Delta_A + (1-x)/Delta_B and
  !> E_nu,v = Delta_v (x E_nu,A/Delta_A + (1-x) E_nu,B/Delta_B).  mdlm: a copy of mdl with the random parts removed.
  subroutine asr_ed_medium(mdl, ed, x, mdlm, edm)
    type(orb_alloy_t), intent(in) :: mdl
    type(asr_ed_t), intent(in) :: ed
    real(dp), intent(in) :: x
    type(orb_alloy_t), intent(inout) :: mdlm
    type(asr_ed_t), intent(out) :: edm
    if (.not. allocated(mdlm%c0)) allocate (mdlm%c0(9), mdlm%c1(9), mdlm%d0(9), mdlm%d1(9), mdlm%o0(9), mdlm%o1(9), mdlm%e0(9), mdlm%e1(9))
    if (.not. allocated(mdlm%p0)) allocate (mdlm%p0(9), mdlm%p1(9))
    mdlm%c0 = mdl%c0 + x * mdl%c1; mdlm%c1 = 0.0_dp
    mdlm%d0 = mdl%d0 + x * mdl%d1; mdlm%d1 = 0.0_dp
    mdlm%o0 = mdl%o0 + x * mdl%o1; mdlm%o1 = 0.0_dp
    mdlm%e0 = 0.0_dp; mdlm%e1 = 0.0_dp; mdlm%p0 = 0.0_dp; mdlm%p1 = 0.0_dp
    edm%dc = ed%dc
    edm%dela = 1.0_dp / (x / ed%dela + (1.0_dp - x) / ed%delb); edm%delb = edm%dela
    edm%enua = edm%dela * (x * ed%enua / ed%dela + (1.0_dp - x) * ed%enub / ed%delb); edm%enub = edm%enua
    edm%mc = ed%mc
  end subroutine asr_ed_medium

  !> k-space form of the medium at one k:  z - H_med(z, k) = R^1/2 (z - Kt) R^1/2  with R = diag(Delta_c / Delta_v),
  !> Kt = R^-1/2 (h' - h' q h') R^-1/2 + diag(E_nu,v),  h' = g + d S(k) d  (Hermitian and independent of z).  Returns the
  !> eigenvalues lam and the weights wt(L, j) = |U_Lj|^2 / R_L, so that  [z - H_med]^-1_LL = sum_j wt(L, j) / (z - lam(j)).
  subroutine asr_ed_medium_k(mdlm, edm, sk, lam, wt)
    type(orb_alloy_t), intent(in) :: mdlm
    type(asr_ed_t), intent(in) :: edm
    complex(dp), intent(in) :: sk(9, 9)
    real(dp), intent(out) :: lam(9), wt(9, 9)
    complex(dp) :: hp(9, 9), kk(9, 9), cw(36)
    real(dp) :: rw(25), rr(9)
    integer :: i, j, info
    external :: zheev
    do j = 1, 9
      do i = 1, 9
        hp(i, j) = mdlm%d0(i) * sk(i, j) * mdlm%d0(j)
      end do
      hp(j, j) = hp(j, j) + mdlm%c0(j)
    end do
    do j = 1, 9
      do i = 1, 9
        kk(i, j) = hp(i, j) - sum(hp(i, :) * mdlm%o0 * hp(:, j))
      end do
    end do
    rr = edm%dc / edm%dela
    do j = 1, 9
      do i = 1, 9
        kk(i, j) = kk(i, j) / sqrt(rr(i) * rr(j))
      end do
      kk(j, j) = kk(j, j) + edm%enua(j)
    end do
    kk = 0.5_dp * (kk + conjg(transpose(kk)))
    call zheev('V', 'U', 9, kk, 9, lam, cw, size(cw), rw, info)
    if (info /= 0) error stop 'asr_ed_medium_k: zheev failed'
    do j = 1, 9
      wt(:, j) = abs(kk(:, j))**2 / rr
    end do
  end subroutine asr_ed_medium_k

  !> G_c(E) for the start orbital l0 (`iv`: index of the start component in the operator's vector, 9 (node-1) + l0).
  !> The operator `op` carries the parameters set by asr_ed_setup, the start site frozen to the central species, and `mapply`.
  !> wk must have been initialised with ed_work_init(wk, op%dim()).
  subroutine asr_ed_green(op, wk, ed, z, iv, nlev, g, acoef, bcoef, mused)
    class(linop), intent(in) :: op
    type(ed_work_t), intent(inout) :: wk
    type(asr_ed_t), intent(in) :: ed
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
  end subroutine asr_ed_green
end module asr_ed
