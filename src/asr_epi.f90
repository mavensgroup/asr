!> Effective pair interactions (EPI) from the augmented-space recursion: the input of the short-range order thermodynamics
!> (tools/sro_alpha).  N species, any shell, any pair of sublattices.
!>
!> Quantity.  Omega(c, s) is the grand potential (frozen potentials, fixed E_F) with species c on site 0, species s on a site j
!> of shell n, and every other site averaged.  The pair table is the mixed second difference against the first species c0, s0
!> of the two sublattices,
!>     U(c, s) = Omega(c, s) - Omega(c0, s) - Omega(c, s0) + Omega(c0, s0).
!> Site terms cancel in it.  The configuration energy of a pair model is E = sum over pairs of U(s_i, s_j), up to terms that
!> canonical exchanges do not change.  On one sublattice the exchange interaction of the unlike pair (p, q) is
!>     V(p, q) = U(p, p) + U(q, q) - U(p, q) - U(q, p)       (V > 0: unlike neighbours are favoured, ordering),
!> and the SRO energy per site is sum_n (z_n / 2) sum_{p < q} x_p x_q V_n(p, q) alpha_n(p, q).
!>
!> Method: thermodynamic integration over the species of site 0.  The H^(2) parameters of site 0 go from c0 to c on a straight
!> path, X(lambda) = (1 - lambda) X_c0 + lambda X_c for X = E_nu, c = C - E_nu, d = sqrt(Delta), o (alpha representation).  Then
!>     dOmega/dlambda = Tr[f(H) dH/dlambda] = -(1/pi) Im sum_j w_j Tr[G(z_j) dH/dlambda]          (Fermi contour z_j, w_j).
!> H = E_nu + h - h o h with h = c + d S d, so dH/dlambda touches site 0 and its bonds.  With the central blocks
!>     X0 = G_00,   X1 = (h G)_00,   X2 = (h G h)_00
!> and the identity (h o h G)_00 = 1 - (z - E_nu) X0 + X1 the trace closes on the diagonal elements of site 0:
!>     Tr[G dH/dlambda] = sum_L { (E_nu' + c') X0 - 2 c' o X1 - o' X2 + 2 kappa [ (z - E_nu - c) X0 + 2 c o X1 - o X2 - 1 ] }_LL,
!> kappa = d'/d, all parameters those of site 0 at lambda.  The constant -2 kappa is real and the same for every partner
!> species.  It is left out.  The expression is linear in G.  So its configuration average is exact in augmented space: the
!> blocks are those of the conditionally averaged G with site 0 on the path species and site j frozen to s.  No logarithm of an
!> averaged determinant is needed (the Lloyd-formula route of the older orbital-peeling work makes that approximation).
!>     U(c, s) = int_0^1 dlambda [ dOmega/dlambda (c0 -> c | s) - dOmega/dlambda (c0 -> c | s0) ]        (Gauss-Legendre).
!> The scaling T_R = (D_c / D_s)^1/2 of asr_edc is 1 on site 0, so X0, X1, X2 are the central blocks of the scaled chains with
!> h -> h' (asr_augcls::augcls_hp_c).  Two chains per start vector as in asr_ed3: from |u> (X0, and X1 by projection on h'|u>)
!> and from h'|u> (X2, and X1 again); the mean X1 is used and the difference measures the tail.  The chains run the full depth:
!> U is a small difference of large numbers, and an adaptive depth would add noise.
!>
!> Symmetry.  The group is the abelian part of the stabiliser of the pair in D4h about z.  Its irreducible representations are
!> one dimensional.  The nine start vectors are symmetry-adapted combinations of the real harmonics of one l (epi_adapted), and
!> the trace over L is the sum over them.
!>
!> Pair exchange.  On one sublattice the exact table is symmetric, U(c, s) = U(s, c).  The two entries come from different
!> chains (path site and frozen site interchanged), so their difference measures the error of the chains and of the path
!> quadrature.  The option sym computes one of them only (epi_skip, epi_mirror): fewer chains, no error measure.
!>
!> Path points (test_epi): two Gauss-Legendre nodes give 0.01 % for species that are close (neighbours in a transition
!> series); one node is not safe (2 % of the largest entry, a factor 2 on a small one); far species need four (0.1 %).
!>
!> Energy blocks.  epi_blocks_blk runs the two chains for a block of contour energies in lock step with the kernels of
!> asr_augblk (one pass over the operator tables per step for all members), as edc_green_blk does in the SCF.  The chains are
!> those of epi_blocks; the results differ at rounding level.
!>
!> Route ed3 (full generalised ASA problem, asr_ed3).  The secular matrix K(E) = E O - H is linear in E, so again
!>     dOmega/dlambda = Sum_i f(e_i) de_i/dlambda = (1/pi) Im sum_j w_j Tr[G(z_j) dK/dlambda],   G = K^-1.
!> Unscaled, K = eps - h + u h + h u - h Q h with eps = z - E_nu, u = o eps, Q = o - p eps (o, p in the representation alpha)
!> and h = c + d S d.  Site 0 on the path changes eps, u, Q, c of site 0 and d of its bonds.  The identity (K G)_00 = 1 removes
!> the terms (h u G)_00 - (h Q h G)_00 that do not close on site 0, and
!>     Tr[G dK/dlambda] = sum_L { -(E_nu' + c') X0 + 2 u' X1 + 2 c' u X0 - Q' X2 - 2 c' Q X1
!>                                + 2 kappa [ 1 - (eps - c) X0 - 2 c u X0 + 2 c Q X1 - Q X2 ] }_LL
!> with u' = o' eps - o E_nu', Q' = o' - p' eps + p E_nu', and X0, X1, X2 the central blocks of asr_ed3 (ed3_blocks: the
!> chains of K' from |a> = |u> + o |v> and |v> = h'|u>).  For p = 0, u -> 0 and Q -> -o this is the H^(2) expression above.
!> epi3_integrand returns minus this trace (the form of epi_integrand), without the constant 2 kappa.
!>
!> Route: H^(2) (route = 'edc') or the full ASA (route = 'ed3').  The screened Coulomb part of the pair table comes from
!> asr_screening::screen_shells.
module asr_epi
  use asr_kinds, only: dp
  use asr_augcls, only: class_alloy_t, augcls_op, augcls_hp_c, augcls_rand_c, augcls_project_c
  use asr_augblk, only: blk_work_t, blk_hamx, blk_project
  use asr_edn, only: edn_t
  use asr_lmto_ham, only: potpar_t, to_alpha, o_alpha
  use asr_endep, only: cf_species, ed_project, ed_nchain, ed_nlev
  use asr_symmetry, only: cubic_site_perms, class_ops
  use asr_lapack, only: dsyev
  use asr_ed3, only: ed3_t, ed3_work_t
  implicit none
  private
  public :: epi_par_t, epi_path_t, epi_work_t, epi_pair_t, epi_table_t
  public :: epi_species, epi_path, epi_setup, epi_integrand, epi_adapted, epi_pair_ops, epi_pair_list
  public :: epi_work_init, epi_start, epi_blocks, epi_blocks_blk, epi_exchange, epi_skip, epi_mirror
  public :: epi3_setup, epi3_start, epi3_integrand

  integer, parameter :: orb_l(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]

  !> Parameters of one species in the screened representation alpha, per orbital (p: route ed3 only)
  type :: epi_par_t
    real(dp) :: e(9) = 0.0_dp, c(9) = 0.0_dp, d(9) = 0.0_dp, o(9) = 0.0_dp   !! E_nu, C - E_nu, |sqrt(Delta)|, o
    real(dp) :: p(9) = 0.0_dp                                               !! p
  end type epi_par_t

  !> site 0 on the path between two species: parameters at lambda and their derivatives
  type :: epi_path_t
    real(dp) :: e(9) = 0.0_dp, c(9) = 0.0_dp, d(9) = 0.0_dp, o(9) = 0.0_dp
    real(dp) :: de(9) = 0.0_dp, dc(9) = 0.0_dp, dk(9) = 0.0_dp, dox(9) = 0.0_dp   !! d/dlambda of E_nu, c, ln d, o
    real(dp) :: p(9) = 0.0_dp, dpx(9) = 0.0_dp                                  !! p and d/dlambda of p (route ed3)
  end type epi_path_t

  type :: epi_work_t
    real(dp), allocatable :: re(:), im(:), sre(:), sim(:)
    complex(dp), allocatable :: u(:), up(:), w(:), dg(:), tv(:)
    complex(dp), allocatable :: s0(:), vv(:)               !! |u> and |v> = h'|u> of the current start vector
    real(dp) :: nv = 0.0_dp                                !! |v|
  end type epi_work_t

  !> one pair type: site 0 (origin) and a representative partner site
  type :: epi_pair_t
    integer :: shell = 0, site = 0, mult = 0, r2 = 0, cls2 = 0, nop = 0
    integer :: vec(3) = 0                                  !! partner - origin, half-lattice units
    integer, allocatable :: orbit(:, :)                    !! orbit(:, 1:mult): the vectors of all sites of the pair type
  end type epi_pair_t

  !> result of one pair type: origin on sublattice kap, partner on kap2
  type :: epi_table_t
    integer :: kap = 0, kap2 = 0, shell = 0, mult = 0, nop = 0, nodes = 0
    integer :: vec(3) = 0                                  !! partner - origin, half-lattice units
    integer, allocatable :: orbit(:, :)                    !! the vectors of all sites of the pair type
    real(dp) :: dist = 0.0_dp                              !! Bohr
    real(dp) :: wr = 0.0_dp                                !! screening kernel W(R) (Ry / e^2)
    real(dp) :: tail = 0.0_dp                              !! largest relative mismatch of X1 between the two chains
    real(dp), allocatable :: uband(:, :), ues(:, :)        !! U(k, k') (Ry): species k of kap on the origin, k' of kap2 on the partner
  end type epi_table_t

contains

  !> H^(2) parameters of a species from its potential parameters pp (any representation) and the screening alpha
  pure function epi_species(pp, alpha) result(par)
    type(potpar_t), intent(in) :: pp
    real(dp), intent(in) :: alpha(0:2)
    type(epi_par_t) :: par
    type(potpar_t) :: q
    real(dp) :: o(0:2)
    integer :: i, l
    q = to_alpha(pp, alpha); o = o_alpha(q, alpha)
    do i = 1, 9
      l = orb_l(i)
      par%e(i) = q%enu(l); par%c(i) = q%c(l) - q%enu(l); par%d(i) = abs(q%sqrdel(l)); par%o(i) = o(l); par%p(i) = q%p(l)
    end do
  end function epi_species

  !> straight path from species a (lambda = 0) to species b (lambda = 1)
  pure function epi_path(a, b, lam) result(pt)
    type(epi_par_t), intent(in) :: a, b
    real(dp), intent(in) :: lam
    type(epi_path_t) :: pt
    pt%e = a%e + lam * (b%e - a%e); pt%de = b%e - a%e
    pt%c = a%c + lam * (b%c - a%c); pt%dc = b%c - a%c
    pt%d = a%d + lam * (b%d - a%d); pt%dk = (b%d - a%d) / pt%d
    pt%o = a%o + lam * (b%o - a%o); pt%dox = b%o - a%o
    pt%p = a%p + lam * (b%p - a%p); pt%dpx = b%p - a%p
  end function epi_path

  !> Operator parameters of `mdl` (species 1 .. ns real, ns + 1 = site 0 on the path; mdl%nsp = ns + 1) in the scaling of
  !> asr_edc with the central Delta of the path point: c = g_s = (D_c / D_s) c_s, d = D_c^1/2, o = q_s = o_s D_s / D_c, e = 0.
  !> ed carries D_c, D_s and E_nu,s for the energy-dependent diagonal e_s(z) = z - D_c (z - E_nu,s) / D_s.
  subroutine epi_setup(mdl, par, pt, ed)
    type(class_alloy_t), intent(inout) :: mdl
    type(epi_par_t), intent(in) :: par(:)
    type(epi_path_t), intent(in) :: pt
    type(edn_t), intent(out) :: ed
    real(dp) :: dc, ds
    integer :: ns, s, i
    ns = size(par)
    if (mdl%nsp /= ns + 1) error stop 'epi_setup: mdl%nsp must be the number of species + 1'
    ed%ns = ns + 1; ed%central = ns + 1
    allocate (ed%dc(9), ed%del(9, ns + 1), ed%enu(9, ns + 1))
    if (allocated(mdl%c)) deallocate (mdl%c, mdl%d, mdl%o, mdl%e)
    allocate (mdl%c(9, ns + 1), mdl%d(9, ns + 1), mdl%o(9, ns + 1), mdl%e(9, ns + 1))
    mdl%e = 0.0_dp
    do i = 1, 9
      dc = pt%d(i)**2
      ed%dc(i) = dc
      do s = 1, ns
        ds = par(s)%d(i)**2
        mdl%c(i, s) = dc / ds * par(s)%c(i); mdl%d(i, s) = pt%d(i); mdl%o(i, s) = par(s)%o(i) * ds / dc
        ed%del(i, s) = ds; ed%enu(i, s) = par(s)%e(i)
      end do
      mdl%c(i, ns + 1) = pt%c(i); mdl%d(i, ns + 1) = pt%d(i); mdl%o(i, ns + 1) = pt%o(i)
      ed%del(i, ns + 1) = dc; ed%enu(i, ns + 1) = pt%e(i)
    end do
  end subroutine epi_setup

  !> [G dH/dlambda]_LL of the module header for one start vector of orbital character l (an orbital index 1..9) from its
  !> central blocks x(0:2) = X0, X1, X2 at z; the constant -2 kappa is left out
  pure function epi_integrand(pt, l, z, x) result(t)
    type(epi_path_t), intent(in) :: pt
    integer, intent(in) :: l
    complex(dp), intent(in) :: z, x(0:)
    complex(dp) :: t
    t = (pt%de(l) + pt%dc(l)) * x(0) - 2.0_dp * pt%dc(l) * pt%o(l) * x(1) - pt%dox(l) * x(2) &
        + 2.0_dp * pt%dk(l) * ((z - pt%e(l) - pt%c(l)) * x(0) + 2.0_dp * pt%c(l) * pt%o(l) * x(1) - pt%o(l) * x(2))
  end function epi_integrand

  !> Route ed3: operator parameters of `mdl` and the coefficients ed of asr_ed3 for species 1 .. ns (par) and site 0 on the path
  !> point pt (species ns + 1, the central species): the scaling T_R = (D_c / D_s)^1/2 with the Delta of the path point, as
  !> asr_ed3::ed3_setup builds it from potential parameters.
  subroutine epi3_setup(mdl, par, pt, ed)
    type(class_alloy_t), intent(inout) :: mdl
    type(epi_par_t), intent(in) :: par(:)
    type(epi_path_t), intent(in) :: pt
    type(ed3_t), intent(out) :: ed
    real(dp) :: dc, ds
    integer :: ns, s, i
    ns = size(par)
    if (mdl%nsp /= ns + 1) error stop 'epi3_setup: mdl%nsp must be the number of species + 1'
    ed%ns = ns + 1; ed%central = ns + 1
    allocate (ed%dc(9), ed%del(9, ns + 1), ed%enu(9, ns + 1), ed%o(9, ns + 1), ed%p(9, ns + 1), ed%g(9, ns + 1))
    if (allocated(mdl%c)) deallocate (mdl%c, mdl%d, mdl%o, mdl%e)
    allocate (mdl%c(9, ns + 1), mdl%d(9, ns + 1), mdl%o(9, ns + 1), mdl%e(9, ns + 1))
    mdl%o = 0.0_dp; mdl%e = 0.0_dp
    do i = 1, 9
      dc = pt%d(i)**2
      ed%dc(i) = dc
      do s = 1, ns
        ds = par(s)%d(i)**2
        ed%del(i, s) = ds; ed%enu(i, s) = par(s)%e(i); ed%o(i, s) = par(s)%o(i); ed%p(i, s) = par(s)%p(i)
        ed%g(i, s) = dc / ds * par(s)%c(i)
      end do
      ed%del(i, ns + 1) = dc; ed%enu(i, ns + 1) = pt%e(i); ed%o(i, ns + 1) = pt%o(i); ed%p(i, ns + 1) = pt%p(i)
      ed%g(i, ns + 1) = pt%c(i)
      mdl%c(i, :) = ed%g(i, :); mdl%d(i, :) = pt%d(i)
    end do
  end subroutine epi3_setup

  !> Route ed3: start vectors of one symmetry-adapted combination v9 of orbital character l (an orbital index 1..9) on the start
  !> node: |u>, |v> = h'|u>, |a> = |u> + o_c |v> (as asr_ed3::ed3_start for one orbital).  op carries epi3_setup's parameters,
  !> the frozen sites and the characters of v9.
  subroutine epi3_start(op, wk, ed, v9, l)
    type(augcls_op), intent(in) :: op
    type(ed3_work_t), intent(inout) :: wk
    type(ed3_t), intent(in) :: ed
    real(dp), intent(in) :: v9(9)
    integer, intent(in) :: l
    wk%s0 = (0.0_dp, 0.0_dp); wk%s0(1:9) = cmplx(v9, 0.0_dp, dp)
    call augcls_hp_c(op, wk%s0, wk%vv)
    if (ed_project) call augcls_project_c(op, wk%vv)
    wk%av = wk%s0 + ed%o(l, ed%central) * wk%vv
    wk%nv = sqrt(real(sum(wk%vv * wk%vv), dp)); wk%na = sqrt(real(sum(wk%av * wk%av), dp))
    if (wk%nv < 1.0e-14_dp .or. wk%na < 1.0e-14_dp) error stop 'epi3_start: zero start vector'
  end subroutine epi3_start

  !> Route ed3: minus Tr[G dK/dlambda]_LL of the module header for one start vector of orbital character l from the central
  !> blocks x(0:2) = X0, X1, X2 of asr_ed3::ed3_blocks at z; the constant 2 kappa is left out.  The sign makes it the
  !> counterpart of epi_integrand: dOmega/dlambda = -(1/pi) Im sum_j w_j (sum of the integrands).
  pure function epi3_integrand(pt, l, z, x) result(t)
    type(epi_path_t), intent(in) :: pt
    integer, intent(in) :: l
    complex(dp), intent(in) :: z, x(0:)
    complex(dp) :: t, eps, u, du, q, dq
    eps = z - pt%e(l)
    u = pt%o(l) * eps;  du = pt%dox(l) * eps - pt%o(l) * pt%de(l)
    q = pt%o(l) - pt%p(l) * eps;  dq = pt%dox(l) - pt%dpx(l) * eps + pt%p(l) * pt%de(l)
    t = -(pt%de(l) + pt%dc(l)) * x(0) + 2.0_dp * du * x(1) + 2.0_dp * pt%dc(l) * u * x(0) - dq * x(2) &
        - 2.0_dp * pt%dc(l) * q * x(1) &
        + 2.0_dp * pt%dk(l) * (-(eps - pt%c(l)) * x(0) - 2.0_dp * pt%c(l) * u * x(0) + 2.0_dp * pt%c(l) * q * x(1) - q * x(2))
    t = -t
  end function epi3_integrand

  !> Exchange interactions of one sublattice from the pair table u(c, s): v(p, q) = u(p,p) + u(q,q) - u(p,q) - u(q,p)
  pure function epi_exchange(u) result(v)
    real(dp), intent(in) :: u(:, :)
    real(dp) :: v(size(u, 1), size(u, 1))
    integer :: p, q
    do q = 1, size(u, 1)
      do p = 1, size(u, 1)
        v(p, q) = u(p, p) + u(q, q) - u(p, q) - u(q, p)
      end do
    end do
  end function epi_exchange

  !> Pair-exchange symmetry on one sublattice: the exact table is symmetric, U(c, s) = U(s, c).  With sym the entries below
  !> the diagonal (1 < s < c) are not computed (epi_skip) and are copied from above it (epi_mirror).  The column s = 1 is
  !> always computed: it is the reference of every row.
  pure logical function epi_skip(sym, same, c, s)
    logical, intent(in) :: sym, same               !! option; origin and partner on the same sublattice
    integer, intent(in) :: c, s
    epi_skip = sym .and. same .and. s > 1 .and. s < c
  end function epi_skip

  pure subroutine epi_mirror(u)
    real(dp), intent(inout) :: u(:, :)
    integer :: c, s
    do c = 3, size(u, 1)
      do s = 2, c - 1
        u(c, s) = u(s, c)
      end do
    end do
  end subroutine epi_mirror

  !> Symmetry-adapted start vectors for an abelian group of involutions with the orbital matrices d(9, 9, nop) (asr_symmetry::
  !> orbital_rep): vec(:, a), a = 1..9, orthonormal, each inside one l, with D(g) vec = chi(g, a) vec, chi = +-1.  lo(a) is an
  !> orbital index of that l (1, 2 or 5).  With the identity alone the vectors are the nine harmonics.
  subroutine epi_adapted(d, vec, chi, lo)
    real(dp), intent(in) :: d(:, :, :)
    real(dp), intent(out) :: vec(9, 9), chi(:, :)
    integer, intent(out) :: lo(9)
    integer, parameter :: i1(3) = [1, 2, 5], i2(3) = [1, 4, 9]
    real(dp) :: m(5, 5), ev(5), work(64), w(9), x
    integer :: nop, g, b, n, info, a, k
    nop = size(d, 3)
    if (size(chi, 1) < nop .or. size(chi, 2) < 9) error stop 'epi_adapted: chi too small'
    do g = 1, nop
      if (maxval(abs(d(:, :, g) - transpose(d(:, :, g)))) > 1.0e-9_dp) &
        error stop 'epi_adapted: an operation is not an involution (the group must be abelian, all elements of order 2)'
    end do
    vec = 0.0_dp
    do b = 1, 3
      n = i2(b) - i1(b) + 1
      m = 0.0_dp
      do g = 1, nop                                     ! generic combination: one eigenvalue per character
        m(1:n, 1:n) = m(1:n, 1:n) + sqrt(real(g, dp) + 0.37_dp) * d(i1(b):i2(b), i1(b):i2(b), g)
      end do
      m(1:n, 1:n) = 0.5_dp * (m(1:n, 1:n) + transpose(m(1:n, 1:n)))
      call dsyev('V', 'U', n, m, 5, ev, work, size(work), info)
      if (info /= 0) error stop 'epi_adapted: dsyev failed'
      vec(i1(b):i2(b), i1(b):i2(b)) = m(1:n, 1:n)
      lo(i1(b):i2(b)) = i1(b)
    end do
    do a = 1, 9
      where (abs(vec(:, a)) < 1.0e-12_dp) vec(:, a) = 0.0_dp
      k = maxloc(abs(vec(:, a)), 1)
      if (vec(k, a) < 0.0_dp) vec(:, a) = -vec(:, a)
      vec(:, a) = vec(:, a) / norm2(vec(:, a))
      do g = 1, nop
        w = matmul(d(:, :, g), vec(:, a)); x = dot_product(vec(:, a), w)
        if (abs(abs(x) - 1.0_dp) > 1.0e-8_dp .or. maxval(abs(w - x * vec(:, a))) > 1.0e-8_dp) &
          error stop 'epi_adapted: a start vector is not a common eigenvector of the group'
        chi(g, a) = sign(1.0_dp, x)
      end do
    end do
  end subroutine epi_adapted

  !> Group of a pair operator: the operations of D4h about z that keep the classes cls (the partner site carries a class of its
  !> own, so they fix it), reduced to an abelian group of involutions.  If the class-preserving operations do not all commute
  !> (partner on the z axis: C4v), the sign changes of the axes among them are kept.  Operation 1 is the identity.
  subroutine epi_pair_ops(coord, mg, cls, perm, mats)
    integer, intent(in) :: coord(:, :), mg, cls(:)
    integer, allocatable, intent(out) :: perm(:, :), mats(:, :, :)
    integer, allocatable :: p0(:, :), m0(:, :, :)
    logical, allocatable :: keep(:)
    logical :: abel
    integer :: n, i, j, k
    integer, parameter :: one(3, 3) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
    call class_ops(coord, mg, cls, p0, m0)
    n = size(m0, 3)
    abel = .true.
    do i = 1, n
      if (any(matmul(m0(:, :, i), m0(:, :, i)) /= one)) abel = .false.
      do j = i + 1, n
        if (any(matmul(m0(:, :, i), m0(:, :, j)) /= matmul(m0(:, :, j), m0(:, :, i)))) abel = .false.
      end do
    end do
    allocate (keep(n)); keep = .true.
    if (.not. abel) then
      do i = 1, n
        keep(i) = m0(1, 2, i) == 0 .and. m0(1, 3, i) == 0 .and. m0(2, 1, i) == 0 .and. m0(2, 3, i) == 0 &
                  .and. m0(3, 1, i) == 0 .and. m0(3, 2, i) == 0
      end do
    end if
    allocate (perm(size(p0, 1), count(keep)), mats(3, 3, count(keep)))
    k = 0
    do i = 1, n
      if (keep(i)) then; k = k + 1; perm(:, k) = p0(:, i); mats(:, :, k) = m0(:, :, i); end if
    end do
    if (any(mats(:, :, 1) /= one)) error stop 'epi_pair_ops: operation 1 is not the identity'
  end subroutine epi_pair_ops

  !> Pair types of the origin (site 1) with the sites of its first nshell coordination shells.  A shell is one distance; it is
  !> split into the orbits of the class-preserving cubic operations (cls: classes of the decorated lattice seen from the
  !> origin).  Each orbit gives one pair type with the orbit size as multiplicity.  Its representative is the member with the
  !> largest group of epi_pair_ops.  msg is empty on success.
  subroutine epi_pair_list(coord, mg, cls, nshell, pairs, msg)
    integer, intent(in) :: coord(:, :), mg, cls(:), nshell
    type(epi_pair_t), allocatable, intent(out) :: pairs(:)
    character(len=:), allocatable, intent(out) :: msg
    integer, allocatable :: v(:, :), r2(:), p48(:, :), m48(:, :, :), cl2(:), pp(:, :), mm(:, :, :), shells(:)
    logical, allocatable :: keep(:), done(:), inorb(:)
    type(epi_pair_t) :: buf(48 * max(nshell, 1)), cand
    integer :: n, i, k, j, np, ish, rcur, rnext
    msg = ''
    n = size(coord, 2)
    allocate (v(3, n), r2(n), done(n), inorb(n), shells(max(nshell, 1)))
    do i = 1, n
      v(:, i) = modulo(coord(:, i) - coord(:, 1) + mg / 2, mg) - mg / 2
      r2(i) = sum(v(:, i)**2)
    end do
    rcur = 0
    do ish = 1, nshell                                   ! the nshell smallest distances
      rnext = minval(r2, mask=r2 > rcur)
      if (4 * rnext >= mg * mg) then
        msg = 'epi: shell '//itoa(ish)//' does not fit into the supercell (increase [asr] supercell)'; allocate (pairs(0)); return
      end if
      shells(ish) = rnext; rcur = rnext
    end do
    call cubic_site_perms(coord, mg, p48, m48)
    allocate (keep(48))
    do k = 1, 48; keep(k) = all(cls(p48(:, k)) == cls); end do
    np = 0; done = .false.
    do ish = 1, nshell
      do i = 2, n
        if (r2(i) /= shells(ish) .or. done(i)) cycle
        inorb = .false.
        do k = 1, 48
          if (keep(k)) inorb(p48(i, k)) = .true.
        end do
        done = done .or. inorb
        cand%nop = 0
        do j = 1, n                                      ! representative: largest pair group
          if (.not. inorb(j)) cycle
          cl2 = cls; cl2(j) = maxval(cls) + 1
          call epi_pair_ops(coord, mg, cl2, pp, mm)
          if (size(pp, 2) > cand%nop) then
            cand%nop = size(pp, 2); cand%site = j; cand%vec = v(:, j)
          end if
        end do
        cand%shell = ish; cand%r2 = shells(ish); cand%mult = count(inorb); cand%cls2 = cls(cand%site)
        if (allocated(cand%orbit)) deallocate (cand%orbit)
        allocate (cand%orbit(3, cand%mult)); k = 0
        do j = 1, n
          if (inorb(j)) then; k = k + 1; cand%orbit(:, k) = v(:, j); end if
        end do
        np = np + 1
        if (np > size(buf)) then; msg = 'epi: too many pair types'; allocate (pairs(0)); return; end if
        buf(np) = cand
      end do
    end do
    pairs = buf(1:np)
  end subroutine epi_pair_list

  subroutine epi_work_init(wk, n)
    type(epi_work_t), intent(out) :: wk
    integer, intent(in) :: n
    allocate (wk%re(n), wk%im(n), wk%sre(n), wk%sim(n), wk%u(n), wk%up(n), wk%w(n), wk%dg(n), wk%tv(n), wk%s0(n), wk%vv(n))
  end subroutine epi_work_init

  !> Start vectors of one symmetry-adapted combination v9 on the start node (node 1): |u> and |v> = h'|u>.  op carries
  !> epi_setup's parameters, the frozen sites, and the characters of v9 in op%chi.
  subroutine epi_start(op, wk, v9)
    type(augcls_op), intent(in) :: op
    type(epi_work_t), intent(inout) :: wk
    real(dp), intent(in) :: v9(9)
    wk%s0 = (0.0_dp, 0.0_dp); wk%s0(1:9) = cmplx(v9, 0.0_dp, dp)
    call augcls_hp_c(op, wk%s0, wk%vv)
    if (ed_project) call augcls_project_c(op, wk%vv)
    wk%nv = sqrt(real(sum(wk%vv * wk%vv), dp))
    if (wk%nv < 1.0e-14_dp) error stop 'epi_start: zero start vector'
  end subroutine epi_start

  !> One complex-symmetric Lanczos chain on H(z) (asr_edc) from wk%u (u^T u = 1), full depth nlev:
  !> g = <u|(z - H(z))^-1|u>,  f = <tv|(z - H(z))^-1|u>, with the tail of cf_species; a closed chain has no tail.
  subroutine chain(op, wk, ed, z, nlev, g, f)
    type(augcls_op), intent(in) :: op
    type(epi_work_t), intent(inout) :: wk
    type(edn_t), intent(in) :: ed
    complex(dp), intent(in) :: z
    integer, intent(in) :: nlev
    complex(dp), intent(out) :: g, f
    complex(dp) :: ep(9, ed%ns), a(nlev), b(nlev), c(nlev)
    integer :: k, l, s, m
    logical :: closed
    do s = 1, ed%ns
      do l = 1, 9
        ep(l, s) = z - ed%dc(l) * (z - ed%enu(l, s)) / ed%del(l, s)
      end do
    end do
    a = (0.0_dp, 0.0_dp); b = (0.0_dp, 0.0_dp); c = (0.0_dp, 0.0_dp)
    wk%up = (0.0_dp, 0.0_dp); m = 0; closed = .false.
    do k = 1, nlev
      c(k) = sum(wk%tv * wk%u)
      wk%re = real(wk%u, dp); wk%im = aimag(wk%u)
      call op%apply(wk%re, wk%sre)
      call op%apply(wk%im, wk%sim)
      call augcls_rand_c(op, ep, wk%u, wk%dg)
      wk%w = cmplx(wk%sre, wk%sim, dp) + wk%dg
      a(k) = sum(wk%u * wk%w)
      m = k
      if (k == nlev) exit
      wk%w = wk%w - a(k) * wk%u
      if (k > 1) wk%w = wk%w - b(k - 1) * wk%up
      if (ed_project) call augcls_project_c(op, wk%w)
      b(k) = sqrt(sum(wk%w * wk%w))
      if (abs(b(k)) < 1.0e-12_dp) then; closed = .true.; exit; end if
      wk%up = wk%u
      wk%u = wk%w / b(k)
    end do
    !$omp atomic
    ed_nchain = ed_nchain + 1
    !$omp atomic
    ed_nlev = ed_nlev + m
    if (closed) then
      call cf_species(a, b, c, m, z, g, f, tail_in='none')
    else
      call cf_species(a, b, c, m, z, g, f)
    end if
  end subroutine chain

  !> Central blocks at z for the start vector of epi_start: x(0) = X0 = <u|G|u>, x(1) = X1 = <u|h'G|u>, x(2) = X2 = <u|h'Gh'|u>,
  !> x(3) = X1 from the |u> chain - X1 from the |v> chain (0 when the chains are exact)
  subroutine epi_blocks(op, wk, ed, z, nlev, x)
    type(augcls_op), intent(in) :: op
    type(epi_work_t), intent(inout) :: wk
    type(edn_t), intent(in) :: ed
    complex(dp), intent(in) :: z
    integer, intent(in) :: nlev
    complex(dp), intent(out) :: x(0:3)
    complex(dp) :: g, fa, fb
    wk%u = wk%s0; wk%tv = wk%vv
    call chain(op, wk, ed, z, nlev, g, fa)
    x(0) = g
    wk%u = wk%vv / wk%nv; wk%tv = wk%s0
    call chain(op, wk, ed, z, nlev, g, fb)
    x(2) = wk%nv**2 * g
    fb = wk%nv * fb
    x(1) = 0.5_dp * (fa + fb)
    x(3) = fa - fb
  end subroutine epi_blocks

  !> chain for a block of energies zs(1:nb) in lock step (asr_augblk): every member starts from s0 (s0^T s0 = 1) and runs the full
  !> depth nlev; g(j), f(j) as chain at zs(j) with the projection vector tv.  A member whose chain closes leaves the block (zero
  !> rows from then on) and gets no tail.  bw%nb >= nb; idle members pad a short block.
  subroutine chain_blk(op, bw, ed, zs, nlev, s0, tv, g, f)
    type(augcls_op), intent(in) :: op
    type(blk_work_t), intent(inout) :: bw
    type(edn_t), intent(in) :: ed
    complex(dp), intent(in) :: zs(:), s0(:), tv(:)
    integer, intent(in) :: nlev
    complex(dp), intent(out) :: g(:), f(:)
    real(dp) :: pa(9, ed%ns), pb(9, ed%ns), keep(bw%nb)
    complex(dp) :: a(nlev, bw%nb), b(nlev, bw%nb), c(nlev, bw%nb), acc(bw%nb), zp(bw%nb)
    integer :: k, l, s, nb, j, i, m(bw%nb)
    logical :: act(bw%nb), closed(bw%nb)
    nb = size(zs)
    zp(1:nb) = zs; zp(nb + 1:) = zs(nb)
    do s = 1, ed%ns
      do l = 1, 9
        pa(l, s) = 1.0_dp - ed%dc(l) / ed%del(l, s)          ! e_s(z) = z - D_c (z - E_nu,s) / D_s = z pa + pb
        pb(l, s) = ed%dc(l) * ed%enu(l, s) / ed%del(l, s)
      end do
    end do
    a = (0.0_dp, 0.0_dp); b = (0.0_dp, 0.0_dp); c = (0.0_dp, 0.0_dp)
    m = nlev; closed = .false.; act = .false.; act(1:nb) = .true.
    associate (u => bw%u, up => bw%up, w => bw%w)
      do i = 1, size(u, 2); u(1:nb, i) = s0(i); u(nb + 1:, i) = (0.0_dp, 0.0_dp); end do
      up = (0.0_dp, 0.0_dp)
      do k = 1, nlev
        acc = (0.0_dp, 0.0_dp)
        do i = 1, size(u, 2); acc = acc + tv(i) * u(:, i); end do
        do j = 1, nb
          if (act(j)) c(k, j) = acc(j)
        end do
        call blk_hamx(op, bw, pa, pb, zp, u, w)
        acc = (0.0_dp, 0.0_dp)
        do i = 1, size(u, 2); acc = acc + u(:, i) * w(:, i); end do
        do j = 1, nb
          if (act(j)) then; a(k, j) = acc(j); m(j) = k; end if
        end do
        if (k == nlev) exit
        do i = 1, size(u, 2)
          w(:, i) = w(:, i) - a(k, :) * u(:, i)
          if (k > 1) w(:, i) = w(:, i) - b(k - 1, :) * up(:, i)
        end do
        if (ed_project) call blk_project(op, w)
        acc = (0.0_dp, 0.0_dp)
        do i = 1, size(w, 2); acc = acc + w(:, i) * w(:, i); end do
        do j = 1, bw%nb
          if (.not. act(j)) cycle
          b(k, j) = sqrt(acc(j))
          if (abs(b(k, j)) < 1.0e-12_dp) then; act(j) = .false.; closed(j) = .true.; end if
        end do
        if (.not. any(act)) exit
        do j = 1, bw%nb                                       ! members that left: zero rows, they stay zero
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
      if (closed(j)) then
        call cf_species(a(:, j), b(:, j), c(:, j), m(j), zs(j), g(j), f(j), tail_in='none')
      else
        call cf_species(a(:, j), b(:, j), c(:, j), m(j), zs(j), g(j), f(j))
      end if
    end do
    !$omp atomic
    ed_nchain = ed_nchain + nb
    !$omp atomic
    ed_nlev = ed_nlev + sum(m(1:nb))
  end subroutine chain_blk

  !> epi_blocks for a block of energies zs(1:nb) in lock step (asr_augblk).  ws holds the start vectors of epi_start, bw the
  !> block work space (bw%nb >= nb); x(0:3, j) as epi_blocks at zs(j).
  subroutine epi_blocks_blk(op, ws, bw, ed, zs, nlev, x)
    type(augcls_op), intent(in) :: op
    type(epi_work_t), intent(in) :: ws
    type(blk_work_t), intent(inout) :: bw
    type(edn_t), intent(in) :: ed
    complex(dp), intent(in) :: zs(:)
    integer, intent(in) :: nlev
    complex(dp), intent(out) :: x(0:, :)
    complex(dp) :: g(size(zs)), fa(size(zs)), fb(size(zs))
    integer :: nb
    nb = size(zs)
    if (bw%nb < nb) error stop 'epi_blocks_blk: block work space too small'
    call chain_blk(op, bw, ed, zs, nlev, ws%s0, ws%vv, g, fa)
    x(0, 1:nb) = g
    call chain_blk(op, bw, ed, zs, nlev, ws%vv / ws%nv, ws%s0, g, fb)
    x(2, 1:nb) = ws%nv**2 * g
    fb = ws%nv * fb
    x(1, 1:nb) = 0.5_dp * (fa + fb)
    x(3, 1:nb) = fa - fb
  end subroutine epi_blocks_blk

  pure function itoa(i) result(s)
    integer, intent(in) :: i
    character(len=:), allocatable :: s
    character(len=16) :: buf
    write (buf, '(i0)') i; s = trim(buf)
  end function itoa
end module asr_epi
