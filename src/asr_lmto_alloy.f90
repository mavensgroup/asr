!> Binary-alloy description for the augmented-space recursion from our own LMTO output:
!> periodic multi-shell graph, screened structure-constant bond blocks S^alpha (one COMMON screening for
!> both species, so S does not depend on the local composition), and species parameters mapped onto the
!> random operators  X(n) = X_B + (X_A - X_B) n   (n = 1: species A).
module asr_lmto_alloy
  use asr_kinds, only: dp
  use asr_augspace, only: graph_t
  use asr_augorb, only: orb_alloy_t
  use asr_lattices, only: periodic_graph_shells
  use asr_strconst, only: lattice_cluster, screened_cluster, avg_ws_radius, hardcore_alpha
  use asr_lmto_ham, only: potpar_t, to_alpha, o_alpha
  implicit none
  private
  public :: build_lmto_alloy, set_alloy_params

  integer, parameter :: orb_l(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]

contains

  !> kind: sc|bcc|fcc;  supercell of l conventional cells;  rcl: S^alpha cluster radius (alat);
  !> ppA, ppB in the representation pp%alpha_rep (the SCF returns gamma);  x = concentration of A.
  !> rinv (optional, alat, >= rcl): radius of the cluster on which 1 - S0 alpha is inverted; only bonds within rcl are kept.
  !> alpha_in (optional): screening constants (s, p, d); default hard spheres of radius sigma*w.
  subroutine build_lmto_alloy(kind, l, rcl, sigma, ppA, ppB, x, g, coord, mdl, rinv, alpha_in)
    character(len=*), intent(in) :: kind
    integer, intent(in) :: l
    real(dp), intent(in) :: rcl, sigma, x
    type(potpar_t), intent(in) :: ppA, ppB
    type(graph_t), intent(out) :: g
    integer, allocatable, intent(out) :: coord(:, :)
    type(orb_alloy_t), intent(out) :: mdl
    real(dp), intent(in), optional :: rinv, alpha_in(0:2)
    integer, allocatable :: bvec(:, :)
    real(dp), allocatable :: vec(:, :), blk(:, :, :)
    real(dp) :: alpha(0:2), al4(0:3), w, u(3), bm(9, 9), bp(9, 9)
    integer :: r2max, k, ia, ib, ns

    w = avg_ws_radius(kind)
    al4 = hardcore_alpha(sigma); alpha = al4(0:2)
    if (present(alpha_in)) alpha = alpha_in
    if (present(rinv)) then
      call lattice_cluster(kind, max(rcl, rinv), vec)
    else
      call lattice_cluster(kind, rcl, vec)
    end if
    call screened_cluster(vec, w, alpha, blk)
    r2max = int((2.0_dp * rcl)**2 + 1.0e-9_dp)
    call periodic_graph_shells(kind, l, r2max, g, coord, bvec)
    ns = size(bvec, 2)
    mdl%nl = 9; mdl%x = x
    allocate (mdl%block(9, 9, ns))
    do k = 1, ns
      u = 0.5_dp * real(bvec(:, k), dp)
      ia = find_vec(vec, -u); ib = find_vec(vec, u)
      if (ia == 0 .or. ib == 0) error stop 'build_lmto_alloy: bond vector missing from the screened cluster'
      bm = blk(:, :, ia)                       ! S_{R,R+v} = S^alpha(R - R' = -v)
      bp = blk(:, :, ib)
      mdl%block(:, :, k) = 0.5_dp * (bm + transpose(bp))
    end do
    allocate (mdl%block0(9, 9))
    mdl%block0 = 0.5_dp * (blk(:, :, 1) + transpose(blk(:, :, 1)))       ! on-site S^alpha (origin is cluster site 1)
    call set_alloy_params(mdl, ppA, ppB, alpha)
  end subroutine build_lmto_alloy

  !> (Re)map the species parameters onto the random operators X(n) = X_B + (X_A - X_B) n.
  subroutine set_alloy_params(mdl, ppA, ppB, alpha)
    type(orb_alloy_t), intent(inout) :: mdl
    type(potpar_t), intent(in) :: ppA, ppB
    real(dp), intent(in) :: alpha(0:2)
    type(potpar_t) :: qa, qb
    real(dp) :: oa(0:2), ob(0:2)
    integer :: i, ll
    qa = to_alpha(ppA, alpha); qb = to_alpha(ppB, alpha)
    oa = o_alpha(qa, alpha); ob = o_alpha(qb, alpha)
    if (.not. allocated(mdl%c0)) allocate (mdl%c0(9), mdl%c1(9), mdl%d0(9), mdl%d1(9), mdl%o0(9), mdl%o1(9), mdl%e0(9), mdl%e1(9))
    if (.not. allocated(mdl%p0)) allocate (mdl%p0(9), mdl%p1(9))
    do i = 1, 9
      ll = orb_l(i)
      mdl%c0(i) = qb%c(ll) - qb%enu(ll);   mdl%c1(i) = (qa%c(ll) - qa%enu(ll)) - (qb%c(ll) - qb%enu(ll))
      mdl%d0(i) = qb%sqrdel(ll);           mdl%d1(i) = qa%sqrdel(ll) - qb%sqrdel(ll)
      mdl%o0(i) = ob(ll);                  mdl%o1(i) = oa(ll) - ob(ll)
      mdl%e0(i) = qb%enu(ll);              mdl%e1(i) = qa%enu(ll) - qb%enu(ll)
      mdl%p0(i) = qb%p(ll) - ob(ll)**2;    mdl%p1(i) = (qa%p(ll) - oa(ll)**2) - (qb%p(ll) - ob(ll)**2)   ! p^gamma = p^alpha - o^2
    end do
  end subroutine set_alloy_params

  pure integer function find_vec(vec, u) result(idx)
    real(dp), intent(in) :: vec(:, :), u(3)
    integer :: i
    idx = 0
    do i = 1, size(vec, 2)
      if (maxval(abs(vec(:, i) - u)) < 1.0e-7_dp) then; idx = i; return; end if
    end do
  end function find_vec
end module asr_lmto_alloy
