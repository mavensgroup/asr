!> Energy-dependent ASR (AMscf form, see asr_ed) on the site-class operator asr_augcls: ordered sublattices, partial disorder
!> and multicomponent alloys.  For a central species c (start site frozen to c, one recursion per (class, species) and orbital)
!>   H(E) = [e_R(E)] + h' - h' q h',  h' = g + S',  S' = D_c^1/2 S D_c^1/2,
!>   g_s = (D_c/D_s)(C_s - E_nu,s),  q_s = o_s D_s/D_c,  e_s(E) = E - D_c (E - E_nu,s)/D_s   (s = species at R, any class),
!> and G_c(E) = <1|(E - H(E))^-1|1>.  T_R = (D_c/D_s)^1/2 is site diagonal for any species layout, so E - H(E) = T (E - H^(2)) T
!> with T = 1 on the central site holds exactly as for the binary.  The chain is the complex-symmetric Lanczos of asr_edn.
module asr_edc
  use asr_kinds, only: dp
  use asr_augcls, only: class_alloy_t, augcls_op, augcls_rand_c, augcls_project_c
  use asr_edn, only: edn_t, edn_work_t, edn_work_init
  use asr_lmto_ham, only: potpar_t, to_alpha, o_alpha
  use asr_endep, only: ed_depth_tol, ed_nchain, ed_nlev, cf_species, ed_project
  use asr_augblk, only: blk_work_t, blk_hamx, blk_project
  implicit none
  private
  public :: edc_setup, edc_green, edc_green_blk, edn_t, edn_work_t, edn_work_init

  integer, parameter :: orb_l(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]

contains

  !> Species parameters of `mdl` (9 orbitals, nsp global species; classes, concentrations and blocks already present) for the
  !> central species `central` from the gamma-representation parameters pp(1:nsp) and the screening alpha.  mdl%e = 0: the E_nu
  !> term is applied separately (complex, per energy).
  subroutine edc_setup(mdl, pp, alpha, central, ed)
    type(class_alloy_t), intent(inout) :: mdl
    type(potpar_t), intent(in) :: pp(:)
    real(dp), intent(in) :: alpha(0:2)
    integer, intent(in) :: central
    type(edn_t), intent(out) :: ed
    type(potpar_t), allocatable :: q(:)
    real(dp) :: o(0:2), dc, ds
    integer :: ns, s, i, ll
    ns = size(pp)
    if (ns /= mdl%nsp) error stop 'edc_setup: number of species parameters /= mdl%nsp'
    allocate (q(ns)); ed%ns = ns; ed%central = central
    allocate (ed%dc(9), ed%del(9, ns), ed%enu(9, ns))
    if (allocated(mdl%c)) deallocate (mdl%c, mdl%d, mdl%o, mdl%e)
    allocate (mdl%c(9, ns), mdl%d(9, ns), mdl%o(9, ns), mdl%e(9, ns))
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
  end subroutine edc_setup

  !> G_c(E) for the start component iv (9 (node-1) + l0; the start node is node 1); op carries edc_setup's parameters with the
  !> start site frozen to the central species and the start orbital set (augcls_set_start).  Optionally the chain coefficients.
  subroutine edc_green(op, wk, ed, z, iv, nlev, g, acoef, bcoef, mused)
    type(augcls_op), intent(in) :: op
    type(edn_work_t), intent(inout) :: wk
    type(edn_t), intent(in) :: ed
    complex(dp), intent(in) :: z
    integer, intent(in) :: iv, nlev
    complex(dp), intent(out) :: g
    complex(dp), intent(out), optional :: acoef(:), bcoef(:)
    integer, intent(out), optional :: mused
    complex(dp) :: ep(9, ed%ns), a(nlev), b(nlev), cc(nlev), gtry, gprev, ftmp, fa
    integer :: k, l, m, s
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
      call augcls_rand_c(op, ep, wk%u, wk%dg)
      wk%w = cmplx(wk%sre, wk%sim, dp) + wk%dg
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
      if (ed_project) call augcls_project_c(op, wk%w)        ! rounding leaves the symmetric subspace and is amplified there
      b(k) = sqrt(sum(wk%w * wk%w))
      if (abs(b(k)) < 1.0e-12_dp) exit
      wk%up = wk%u
      wk%u = wk%w / b(k)
    end do
    !$omp atomic
    ed_nchain = ed_nchain + 1
    !$omp atomic
    ed_nlev = ed_nlev + m
    call cf_species(a, b, cc, m, z, g, fa)
    if (present(acoef)) acoef(1:m) = a(1:m)
    if (present(bcoef)) bcoef(1:m) = b(1:m)
    if (present(mused)) mused = m
  end subroutine edc_green

  !> edc_green for a block of energies zs(1:nb) in lock step (asr_augblk): the same chains, one pass over the operator per step for
  !> all members.  A member leaves the block at its own depth (the test of edc_green); its rows are then zeroed and stay zero, so
  !> its fraction is that of edc_green; mused(j) is its depth.  wk%nb >= nb.
  subroutine edc_green_blk(op, wk, ed, zs, iv, nlev, g, mused)
    type(augcls_op), intent(in) :: op
    type(blk_work_t), intent(inout) :: wk
    type(edn_t), intent(in) :: ed
    complex(dp), intent(in) :: zs(:)
    integer, intent(in) :: iv, nlev
    complex(dp), intent(out) :: g(:)
    integer, intent(out), optional :: mused(:)
    real(dp) :: pa(9, ed%ns), pb(9, ed%ns)
    complex(dp) :: a(nlev, wk%nb), b(nlev, wk%nb), cc(nlev), gtry(wk%nb), gprev(wk%nb), ftmp, fa, acc(wk%nb), zp(wk%nb)
    real(dp) :: keep(wk%nb)
    integer :: k, l, s, nb, j, i, m(wk%nb)
    logical :: act(wk%nb)
    nb = size(zs)
    if (wk%nb < nb) error stop 'edc_green_blk: block work space too small'
    zp(1:nb) = zs; zp(nb + 1:) = zs(nb)                       ! idle members pad a short block (zero rows throughout)
    do s = 1, ed%ns
      do l = 1, 9
        pa(l, s) = 1.0_dp - ed%dc(l) / ed%del(l, s)          ! e_s(z) = z - D_c (z - E_nu,s) / D_s = z pa + pb
        pb(l, s) = ed%dc(l) * ed%enu(l, s) / ed%del(l, s)
      end do
    end do
    associate (u => wk%u, up => wk%up, w => wk%w)
      u = (0.0_dp, 0.0_dp); up = (0.0_dp, 0.0_dp)
      u(1:nb, iv) = (1.0_dp, 0.0_dp)
      a = 0.0_dp; b = 0.0_dp; cc = 0.0_dp; m = nlev; gprev = (0.0_dp, 0.0_dp)
      act = .false.; act(1:nb) = .true.
      do k = 1, nlev
        call blk_hamx(op, wk, pa, pb, zp, u, w)
        acc = (0.0_dp, 0.0_dp)
        do i = 1, size(u, 2); acc = acc + u(:, i) * w(:, i); end do
        do j = 1, wk%nb
          if (act(j)) then; a(k, j) = acc(j); m(j) = k; end if
        end do
        if (k == nlev) exit
        if (ed_depth_tol > 0.0_dp .and. k >= 8 .and. mod(k, 4) == 0) then
          do j = 1, nb
            if (.not. act(j)) cycle
            call cf_species(a(:, j), b(:, j), cc, k, zs(j), gtry(j), ftmp)
            if (k > 8) then
              if (abs(gtry(j) - gprev(j)) <= ed_depth_tol * abs(gtry(j))) act(j) = .false.
            end if
            gprev(j) = gtry(j)
          end do
        end if
        do i = 1, size(u, 2)
          w(:, i) = w(:, i) - a(k, :) * u(:, i)
          if (k > 1) w(:, i) = w(:, i) - b(k - 1, :) * up(:, i)
        end do
        if (ed_project) call blk_project(op, w)
        acc = (0.0_dp, 0.0_dp)
        do i = 1, size(w, 2); acc = acc + w(:, i) * w(:, i); end do
        do j = 1, wk%nb
          if (.not. act(j)) cycle
          b(k, j) = sqrt(acc(j))
          if (abs(b(k, j)) < 1.0e-12_dp) act(j) = .false.
        end do
        if (.not. any(act)) exit
        do j = 1, wk%nb                                       ! members that left: zero rows, they stay zero
          if (act(j)) then
            acc(j) = 1.0_dp / b(k, j); keep(j) = 1.0_dp
          else
            acc(j) = (0.0_dp, 0.0_dp); keep(j) = 0.0_dp
          end if
        end do
        do i = 1, size(u, 2)
          up(:, i) = keep * u(:, i)
          u(:, i) = acc * w(:, i)
        end do
      end do
    end associate
    do j = 1, nb
      call cf_species(a(:, j), b(:, j), cc, m(j), zs(j), g(j), fa)
    end do
    !$omp atomic
    ed_nchain = ed_nchain + nb
    !$omp atomic
    ed_nlev = ed_nlev + sum(m(1:nb))
    if (present(mused)) mused(1:nb) = m
  end subroutine edc_green_blk
end module asr_edc
