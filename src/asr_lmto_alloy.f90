!> Binary-alloy description for the augmented-space recursion from our own LMTO output:
!> periodic multi-shell graph, screened structure-constant bond blocks S^alpha (one COMMON screening for
!> both species, so S does not depend on the local composition), and species parameters mapped onto the
!> random operators  X(n) = X_B + (X_A - X_B) n   (n = 1: species A).
!> structure_blocks gives the bond blocks of any crystal structure (asr_crystal_graph): one screened cluster per basis site.
module asr_lmto_alloy
  use asr_kinds, only: dp
  use asr_augspace, only: graph_t
  use asr_lattices, only: periodic_graph_shells
  use asr_strconst, only: lattice_cluster, screened_cluster, avg_ws_radius, hardcore_alpha
  use asr_lmto_ham, only: potpar_t, to_alpha, o_alpha
  use asr_crystal, only: structure_t, ws_radius
  use asr_crystal_graph, only: crystal_graph_t, structure_cluster
  implicit none
  private
  public :: build_lmto_alloy, set_alloy_params, orb_alloy_t, structure_blocks

  integer, parameter :: orb_l(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]

  !> binary model on the graph: species parameters as X0 + X1 n (n = 1: A) and the bond blocks S^alpha.  The class SCF
  !> uses the bond blocks only (block, block0).
  type :: orb_alloy_t
    integer :: nl = 1
    real(dp) :: x = 0.5_dp
    real(dp), allocatable :: c0(:), c1(:), d0(:), d1(:), o0(:), o1(:), e0(:), e1(:)   !! (nl); E_nu = e0 + e1 n
    real(dp), allocatable :: p0(:), p1(:)          !! p^gamma (third-order Hamiltonian only)
    real(dp), allocatable :: block(:, :, :)        !! block(L,L',slot): S_{R,R'} for R' = nbr(R,slot)
    real(dp), allocatable :: block0(:, :)          !! on-site block S_{R,R} (optional; zero if not allocated)
  end type orb_alloy_t

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

  !> Screened structure constants S^alpha of structure s on the slots of its graph cg, per basis site:
  !>     block(:, :, k, b) = S_{R, R + bvec(k)},   block0(:, :, b) = S_{R, R}   for R on basis site b.
  !> On a structure with a basis, S depends on the site and not only on the bond vector (an hcp site has no inversion centre, so
  !> p and d mix on site, with opposite signs on the two sites).  So every basis site gets its own screened cluster of radius
  !> rinv (alat units, at least the longest bond of cg), centred on it.  Equal spheres of the average Wigner-Seitz radius and one
  !> common alpha (s, p, d), as in build_lmto_alloy.  The block of a bond is the mean of the two clusters that hold it,
  !>     block(k, b) = (S_b'(-v) + S_b(v)^T) / 2,   b' = btgt(k, b),   S_c(u) = S_{c + u, c} from the cluster of c,
  !> so that block(k', b') = block(k, b)^T for the reverse bond exactly.  With one basis site this is the formula of
  !> build_lmto_alloy, and the clusters of asr_crystal_graph are those of lattice_cluster: cubic blocks are bit-identical.
  !> Slots that basis site b does not have are zero.
  subroutine structure_blocks(s, cg, rinv, alpha, block, block0)
    type(structure_t), intent(in) :: s
    type(crystal_graph_t), intent(in) :: cg
    real(dp), intent(in) :: rinv, alpha(0:2)
    real(dp), allocatable, intent(out) :: block(:, :, :, :), block0(:, :, :)
    type :: clu_t
      real(dp), allocatable :: vec(:, :), blk(:, :, :)
    end type clu_t
    type(clu_t), allocatable :: cl(:)
    real(dp) :: w, u(3)
    integer :: b, k, t, ia, ib
    w = ws_radius(s)
    allocate (cl(s%nb))
    do b = 1, s%nb
      call structure_cluster(s, b, rinv, cl(b)%vec)
      call screened_cluster(cl(b)%vec, w, alpha, cl(b)%blk)
    end do
    allocate (block(9, 9, cg%nslot, s%nb), block0(9, 9, s%nb)); block = 0.0_dp
    do b = 1, s%nb
      block0(:, :, b) = 0.5_dp * (cl(b)%blk(:, :, 1) + transpose(cl(b)%blk(:, :, 1)))
      do k = 1, cg%nslot
        t = cg%btgt(k, b)
        if (t == 0) cycle
        u = cg%bvec(:, k)
        ia = find_vec(cl(t)%vec, -u); ib = find_vec(cl(b)%vec, u)
        if (ia == 0 .or. ib == 0) error stop 'structure_blocks: a bond is longer than the screened cluster radius rinv'
        block(:, :, k, b) = 0.5_dp * (cl(t)%blk(:, :, ia) + transpose(cl(b)%blk(:, :, ib)))
      end do
    end do
  end subroutine structure_blocks

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
