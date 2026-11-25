!> Screened electrostatics of charge fluctuations in a substitutional alloy, single-site form, with first-shell short-range
!> order (derivation below; summary in manual/src/theory/scf_electrostatics.md).
!>
!> Charge-excess model (linear q-V law): every site of the graph lattice responds to its own Madelung potential,
!>     q_i = b_i - (V_i - mu) / a,      V_i = sum_{j /= i} M_ij q_j,      M_ij = 2 / R_ij  (Ry, Bohr; ASA point charges),
!> with one hardness a (Ry per electron) for all sites.  Then q = G b and V = W b with the lattice functions
!>     G(k) = a / (a + M(k)),   W(k) = a M(k) / (a + M(k)) = a (1 - G(k)),
!> M(k) = sum_{R /= 0} (2/R) exp(i k.R) (Ewald).  Shell values: G00 = <G>_BZ, G1 = <G(k) cos(k.R1)>_BZ (same for W).
!> Single-site elimination (one potential per component p = (sublattice, species); first-shell conditional concentrations X):
!>     dq = (G00 I + G1 Ahat) db,  V = (W00 I + W1 Ahat) db   =>   V = F dq,   F = (G00 I + G1 Ahat)^-1 (W00 I + W1 Ahat),
!>     Ahat(p, p') = z(kappa, kappa') (X_s'^(p) - x_s'),   p' = (kappa', s'),   dq_p = q_p - qbar_kappa.
!> No SRO (Ahat = 0): F = (W00/G00) I = -2 beta / R1, the screened-impurity model.  The input is beta; a is solved from
!> beta = -W00 R1 / (2 G00).  beta = 0 is the unscreened limit (a -> infinity: G = 1, W = M, W00 = 0, W1 = 2/R1).  beta has a
!> lattice maximum: a -> a_min = -min_k M(k), where 1 + M/a first vanishes (bcc: a_min R1/2 = 1.7627, the CsCl Madelung constant,
!> at H); G00 stays finite there (integrable singularity), so beta_max = a_min R1/2 (1 - 1/G00) is smaller (bcc about 1.39).
!> diag(w) Ahat is symmetric (w_p = f_kappa x_c; f_kappa z = f_kappa' z', P(c, s) = P(s, c)), so diag(w) F is symmetric and
!> E = 1/2 sum_p w_p dq_p (F dq)_p has dE/dq_p = w_p (F dq)_p.
!>
!> Closure of the first-shell SRO (numerical checks in tools/sro_closure).  The linear form
!> Gamma = G00 + G1 Ahat keeps the correlations of the first shell only.  Its structure factor turns negative at the zone
!> boundary for strong clustering, and G00 + G1 lambda then has a zero (fcc beta 1.0: lambda* = 12.35; beta 1.2: 7.67).  The
!> exact model has no zero: Gamma = <G(k) S(k)> / x with G > 0 and S >= 0, and sum_R G(R) = 0 gives Gamma -> 0 (F -> +infinity)
!> for full segregation.  The code therefore completes the structure factor of every eigenmode of Ahat (eigenvalue lambda,
!> first-shell alpha_1 = lambda / z) with as few correlations beyond the first shell as possible (minimal-tail closure), the
!> configuration model of the ASR itself and of the SQS cells used to test it:
!>   -1 <= lambda <= -1/c1min (fcc -1..3, bcc -1..1): S = 1 + lambda c1(k) >= 0 is realizable; Gamma = G00 + G1 lambda exactly.
!>   lambda > -1/c1min (clustering): S ~ max(0, c1(k) - c0), c0 from c1min upwards; S concentrates at k = 0, where G = 0.
!>   lambda < -1 (ordering):          S ~ max(0, c0 - c1(k)), c0 from 1 downwards; S concentrates where c1 is minimal.
!> c0 follows from alpha_1 = <c1 S> / <S> (bisection; it increases with c0).  Past the last mesh value the order condenses at
!> k = 0 (Gamma -> g_cl0 at alpha_1 = 1, which is 0 on a Bravais lattice) or on the ordering set, where c1 = c1min
!> (Gamma -> g_ord, the average of G over that set; ordering_set).  bcc H and sc R are isolated points: g_ord = G(H), exact
!> for perfect B2 order.  On fcc c1 = -1/3 holds on the whole X-W line, and the first shell does not fix the order along it.
!> The degenerate ground states of the first-shell model fill the line evenly (x = 1/2: any stacking of signed L1_0 planes;
!> their mean exact Gamma is the line average to 4e-16), so g_ord is the line average, between G(X) and G(W).  On a Bravais
!> lattice Gamma is positive, decreasing and continuous on the whole realizable range.
!> Tests (tools/sro_closure): exact model on
!> minimal-tail supercells (inverse Monte Carlo): Gamma within 1 % at beta 1.0 and 1.2 where the linear form is off by 25 %;
!> CuZn alpha -0.148 unchanged.  On thermodynamic (long-tail) supercells the closure is within 4 % at beta 1.0; at beta 1.2 the
!> answer depends on the tails, which the first-shell input does not fix.  No SRO and beta = 0 keep the linear form.
!>
!> Sublattices with n > 1 sites per primitive cell (screen_structure; check: tools/sro_closure/sro_basis_test.py).  The shell is
!> an n x n Hermitian matrix C(k) with the bond phases exp(i k.R) of the full bond vectors R (the gauge of coulomb_k_basis).  Its
!> eigenvalues c1_m(k) and eigenvectors v_m(k) are the branches of a concentration wave.  Branch m has the screening function
!> g_m(k) = v_m^H G_kk(k) v_m, with G = a (a + M(k))^-1 restricted to the sites of the sublattice.  The closure acts on the
!> n Nk branch values without change.  For n = 1 it is the Bravais closure.  The averages <c1_m> = 0, <c1_m^2> = 1/z,
!> <g_m> = G00 and <g_m c1_m> = G1 hold for every n, so the linear range and its edges stay as above.  Gamma stays positive.
!> It need not decrease: dGamma/dalpha_1 = Cov(g, c1) / Var(c1) over the active tail.  In NiAs at c/a 1.6 the anion shell is out
!> of plane.  There Gamma rises from G00 on the ordering side, has a maximum and falls to g_ord < G00 at alpha_1 = -1.  The
!> exact supercell model does the same.  At full clustering (alpha_1 = 1) the clipped structure factor ends evenly on the
!> branches with c1 = 1, so Gamma ends at g_cl0, the average of g over them (clustering_set).  When the shell joins every site
!> of the cell in three dimensions (Bravais lattices, ideal hcp) the only such branch is the uniform mode, and g_cl0 = 0.  On a
!> sublattice smaller than the cell the other sublattices screen (NiAs anions: segregated regions stay neutral through the
!> cation charges).  With a shell in the basal plane any stacking of pure planes has c1 = 1.  There g_cl0 > 0, as the exact
!> supercell model gives.
!>
!> A mode with Gamma = 0 has F = +infinity: its charge deviation must vanish.  V = F dq cannot represent it, but the drive
!> u = db can: dq = Gamma u, V = a (1 - Gamma) u (screen_update).  This form has the same fixed point wherever F is finite.
module asr_screening
  use asr_lapack, only: dgesv, dsyev, zheev, zgesv
  use asr_kinds, only: dp, pi
  use asr_crystal, only: structure_t
  use asr_prototypes, only: prototype_bravais_data, match_cubic_bravais
  implicit none
  private
  public :: screen_t, screen_bravais, screen_lattice, screen_structure, screen_matrix, screen_gamma, screen_update, &
            screen_project, mode_f, coulomb_k, coulomb_k_basis, screen_bravais_shells, screen_shells, screen_structure_shells

  type :: screen_t
    real(dp) :: beta = 0.0_dp            !! input SIM parameter
    real(dp) :: a = huge(1.0_dp)         !! hardness (Ry per electron); huge = unscreened
    real(dp) :: a_min = 0.0_dp           !! stability bound of a
    real(dp) :: beta_max = 0.0_dp        !! largest beta of the lattice
    real(dp) :: r1 = 0.0_dp              !! nearest-neighbour distance (Bohr)
    integer :: z = 0                     !! first-shell coordination of the graph lattice (or correlated sublattice)
    real(dp) :: g00 = 1.0_dp, g1 = 0.0_dp, w00 = 0.0_dp, w1 = 0.0_dp
    real(dp) :: lambda = 1.0_dp          !! linear first-shell weight -a G1 R1 / (2 G00^2) (1 = bare shell)
    logical :: unscreened = .true.
    ! closure of the first-shell SRO (minimal-tail structure factor on the BZ mesh)
    real(dp), allocatable :: gk(:), c1k(:)   !! G(k) and c1(k) on the mesh
    real(dp) :: g_ord = 0.0_dp, c1_ord = -1.0_dp   !! c1 at the ordering set and G averaged over it
    real(dp) :: c0_cl = 1.0_dp, alpha_cl = 1.0_dp, gam_cl = 0.0_dp      !! last clipped point on the mesh, clustering side
    real(dp) :: g_cl0 = 0.0_dp           !! G of complete clustering: g averaged over the c1 = 1 branches (0 on Bravais lattices)
    real(dp) :: c0_ord = -1.0_dp, alpha_ord = -1.0_dp, gam_ord = 0.0_dp !! last clipped point on the mesh, ordering side
  end type screen_t

  abstract interface
    !> Lowest first-shell branch c1(k) at a Cartesian wavevector k (1/Bohr), for ordering_set.
    real(dp) function branch_fn(k)
      import :: dp
      real(dp), intent(in) :: k(3)
    end function branch_fn
  end interface

contains

  !> M(k) = sum_{R /= 0} (2/|R|) cos(k.R) in Ry/e (Ewald; k in 1/Bohr, not a reciprocal lattice vector) for the Bravais lattice
  !> with primitive vectors p (columns, Bohr).  Several k at once (columns of kv).
  subroutine coulomb_k(p, kv, m)
    real(dp), intent(in) :: p(3, 3), kv(:, :)
    real(dp), intent(out) :: m(size(kv, 2))
    real(dp), allocatable :: rv(:, :), wr(:), gv(:, :)
    real(dp) :: b(3, 3), omega, eta, rmax, gmax, v(3), r, q2, s
    integer :: n(3), i1, i2, i3, nr, ng, ik, j
    omega = abs(det3(p))
    b(:, 1) = cross(p(:, 2), p(:, 3)); b(:, 2) = cross(p(:, 3), p(:, 1)); b(:, 3) = cross(p(:, 1), p(:, 2))
    b = 2.0_dp * pi * b / omega
    eta = sqrt(pi) / omega**(1.0_dp / 3.0_dp)
    rmax = 7.0_dp / eta; gmax = 14.0_dp * eta + maxval(norm2(b, dim=1))
    do i1 = 1, 3; n(i1) = ceiling(rmax * norm2(b(:, i1)) / (2.0_dp * pi)) + 1; end do
    allocate (rv(3, product(2 * n + 1)), wr(product(2 * n + 1)))
    nr = 0
    do i1 = -n(1), n(1); do i2 = -n(2), n(2); do i3 = -n(3), n(3)
          v = i1 * p(:, 1) + i2 * p(:, 2) + i3 * p(:, 3); r = norm2(v)
          if (r < 1.0e-10_dp .or. r > rmax) cycle
          nr = nr + 1; rv(:, nr) = v; wr(nr) = erfc(eta * r) / r
        end do; end do; end do
    do i1 = 1, 3; n(i1) = ceiling(gmax * norm2(p(:, i1)) / (2.0_dp * pi)) + 1; end do
    allocate (gv(3, product(2 * n + 1)))
    ng = 0
    do i1 = -n(1), n(1); do i2 = -n(2), n(2); do i3 = -n(3), n(3)
          v = i1 * b(:, 1) + i2 * b(:, 2) + i3 * b(:, 3)
          if (norm2(v) > gmax) cycle
          ng = ng + 1; gv(:, ng) = v
        end do; end do; end do
    do ik = 1, size(kv, 2)
      s = 0.0_dp
      do j = 1, ng
        v = kv(:, ik) + gv(:, j); q2 = dot_product(v, v)
        if (q2 < 1.0e-14_dp) error stop 'coulomb_k: k is a reciprocal lattice vector'
        s = s + exp(-q2 / (4.0_dp * eta**2)) / q2
      end do
      s = 4.0_dp * pi / omega * s
      do j = 1, nr; s = s + cos(dot_product(kv(:, ik), rv(:, j))) * wr(j); end do
      m(ik) = 2.0_dp * (s - 2.0_dp * eta / sqrt(pi))
    end do
  end subroutine coulomb_k

  !> Basis-resolved lattice Coulomb matrix M_bb'(k) in Ry/e (Ewald; k in 1/Bohr, not a reciprocal lattice vector unless
  !> gamma_neutral is true) for a crystal structure with primitive vectors p (columns, Bohr) and Cartesian basis positions
  !> tau(3, nb) (Bohr):
  !>     M_bb'(k) = sum_T (2 / |T + tau_b' - tau_b|) exp(i k . (T + tau_b' - tau_b))   (T = 0 excluded when b = b').
  !> M(k) is Hermitian.  For a single-site Bravais cell (nb = 1, tau = 0) it equals coulomb_k.  eta (optional, 1/Bohr) overrides
  !> the default Ewald splitting parameter.  When gamma_neutral is true and k = 0, the divergent G = 0 term is left out and the
  !> Ewald background term -pi / (Omega eta^2) is kept.  Both are the same for every pair (b, b'), so they change only the
  !> cell-uniform component.  At k = 0 only the neutral block Q^T M Q is defined (Q: orthonormal vectors with zero component
  !> sum); it is the limit k -> 0 of the neutral block of M(k).  Diagonalise it with neutral_eig, which drops the uniform mode.
  subroutine coulomb_k_basis(p, tau, kv, mb, eta, gamma_neutral)
    real(dp), intent(in) :: p(3, 3), tau(:, :), kv(:, :)
    complex(dp), intent(out) :: mb(size(tau, 2), size(tau, 2), size(kv, 2))
    real(dp), intent(in), optional :: eta
    logical, intent(in), optional :: gamma_neutral
    real(dp), allocatable :: tv(:, :), gv(:, :)
    real(dp) :: b(3, 3), omega, et, rmax, gmax, v(3), d(3), r, q2, dtmax, ph
    complex(dp) :: s, sr
    integer :: nb, n(3), i1, i2, i3, nt, ng, ik, j, b1, b2
    logical :: gneut, kzero
    gneut = .false.; if (present(gamma_neutral)) gneut = gamma_neutral
    nb = size(tau, 2)
    omega = abs(det3(p))
    b(:, 1) = cross(p(:, 2), p(:, 3)); b(:, 2) = cross(p(:, 3), p(:, 1)); b(:, 3) = cross(p(:, 1), p(:, 2))
    b = 2.0_dp * pi * b / omega
    et = sqrt(pi) / (omega / real(max(nb, 1), dp))**(1.0_dp / 3.0_dp)
    if (present(eta)) et = eta
    dtmax = 0.0_dp
    do b1 = 1, nb
      do b2 = b1 + 1, nb; dtmax = max(dtmax, norm2(tau(:, b2) - tau(:, b1))); end do
    end do
    rmax = 7.0_dp / et; gmax = 14.0_dp * et + maxval(norm2(b, dim=1))
    do i1 = 1, 3; n(i1) = ceiling((rmax + dtmax) * norm2(b(:, i1)) / (2.0_dp * pi)) + 1; end do
    allocate (tv(3, product(2 * n + 1)))
    nt = 0
    do i1 = -n(1), n(1); do i2 = -n(2), n(2); do i3 = -n(3), n(3)
      v = i1 * p(:, 1) + i2 * p(:, 2) + i3 * p(:, 3)
      if (norm2(v) > rmax + dtmax) cycle
      nt = nt + 1; tv(:, nt) = v
    end do; end do; end do
    do i1 = 1, 3; n(i1) = ceiling(gmax * norm2(p(:, i1)) / (2.0_dp * pi)) + 1; end do
    allocate (gv(3, product(2 * n + 1)))
    ng = 0
    do i1 = -n(1), n(1); do i2 = -n(2), n(2); do i3 = -n(3), n(3)
      v = i1 * b(:, 1) + i2 * b(:, 2) + i3 * b(:, 3)
      if (norm2(v) > gmax) cycle
      ng = ng + 1; gv(:, ng) = v
    end do; end do; end do
    mb = (0.0_dp, 0.0_dp)
    do ik = 1, size(kv, 2)
      kzero = gneut .and. norm2(kv(:, ik)) < 1.0e-12_dp
      do b1 = 1, nb
        do b2 = b1, nb
          d = tau(:, b2) - tau(:, b1)
          s = (0.0_dp, 0.0_dp)
          do j = 1, ng
            v = kv(:, ik) + gv(:, j); q2 = dot_product(v, v)
            if (q2 < 1.0e-14_dp) then
              if (kzero) cycle
              error stop 'coulomb_k_basis: k is a reciprocal lattice vector'
            end if
            ! sum_T erf(eta |T + d|) / |T + d| exp(i k.(T + d)) is a lattice-periodic function of d: its Fourier series has the
            ! phases exp(-i G.d), not exp(i (k + G).d).  With the latter the result depends on eta off the diagonal.
            ph = -dot_product(gv(:, j), d)
            s = s + exp(-q2 / (4.0_dp * et**2)) / q2 * cmplx(cos(ph), sin(ph), dp)
          end do
          s = 4.0_dp * pi / omega * s
          if (kzero) s = s - pi / (omega * et**2)
          sr = (0.0_dp, 0.0_dp)
          do j = 1, nt
            v = tv(:, j) + d; r = norm2(v)
            if (r < 1.0e-10_dp .or. r > rmax) cycle
            ph = dot_product(kv(:, ik), v)
            sr = sr + (erfc(et * r) / r) * cmplx(cos(ph), sin(ph), dp)
          end do
          s = s + sr
          if (b1 == b2) s = cmplx(real(s, dp) - 2.0_dp * et / sqrt(pi), 0.0_dp, dp)
          mb(b1, b2, ik) = 2.0_dp * s
          if (b2 > b1) mb(b2, b1, ik) = conjg(mb(b1, b2, ik))
        end do
      end do
    end do
  end subroutine coulomb_k_basis

  !> Eigenpairs of a basis Coulomb matrix m at k = 0 (coulomb_k_basis with gamma_neutral) on the neutral subspace.  Q is the
  !> Helmert basis: q_j has 1 / sqrt(j (j + 1)) on sites 1 .. j and -j / sqrt(j (j + 1)) on site j + 1 (j = 1 .. nb - 1).  The
  !> nb - 1 neutral eigenpairs come from Q^T m Q, so the uniform mode never enters the eigensolver.  They are returned in
  !> ev(1:nb - 1) and v(:, 1:nb - 1).  Column nb of v is the uniform vector (1, ..., 1) / sqrt(nb) with ev(nb) = 0.  Its
  !> charge costs infinite energy, so the caller must give it a Green function of zero.
  subroutine neutral_eig(m, ev, v)
    complex(dp), intent(in) :: m(:, :)
    real(dp), intent(out) :: ev(size(m, 1))
    complex(dp), intent(out) :: v(size(m, 1), size(m, 1))
    real(dp) :: q(size(m, 1), size(m, 1) - 1), rw(max(1, 3 * size(m, 1)))
    complex(dp) :: a(size(m, 1) - 1, size(m, 1) - 1), cw(max(1, 2 * size(m, 1)))
    integer :: nb, j, info
    nb = size(m, 1)
    q = 0.0_dp
    do j = 1, nb - 1
      q(1:j, j) = 1.0_dp / sqrt(real(j * (j + 1), dp)); q(j + 1, j) = -real(j, dp) / sqrt(real(j * (j + 1), dp))
    end do
    a = matmul(transpose(q), matmul(m, q))
    call zheev('V', 'U', nb - 1, a, nb - 1, ev, cw, size(cw), rw, info)
    if (info /= 0) error stop 'neutral_eig: zheev failed'
    v(:, 1:nb - 1) = matmul(q, a)
    v(:, nb) = 1.0_dp / sqrt(real(nb, dp)); ev(nb) = 0.0_dp
  end subroutine neutral_eig

  !> Wavevectors kc (Cartesian columns, 1/Bohr) that sample the set where the largest eigenvalue of the shell matrix C(k) is 1.
  !> c1 = 1 needs a vector v with exp(i k.R_j) v_r' = v_r on every bond j (r -> r').  That holds exactly when k.L is a multiple
  !> of 2 pi for every closed loop of bonds.  The loop vectors L are lattice translations, and they generate a lattice.  In
  !> fractional reciprocal coordinates f (k = b f) the condition is that f.L is an integer.
  !>   Loops in three directions: the set is finite.  It lies on the grid (1/D) Z^3, D the smallest non-zero |det| of three
  !>   loops, and every grid point with f.L integer for all loops belongs to it (k = 0 always does).
  !>   Loops in two directions only (a shell in the basal plane): the set is lines f0 + s u, u the primitive integer normal of
  !>   the loop plane.  Every line meets the grid (1/D) Z^3 with D = |(L1 x L2).e| (e a unit vector with u.e /= 0); the grid
  !>   points of one line are f0 + (m / D) u.  Each line is sampled at nline midpoints s = (q - 1/2) / nline.
  !> b: reciprocal vectors (columns, 1/Bohr); tau: basis positions (Bohr); sh, sb1, sb2: the shell bonds (structure_sro_shell);
  !> smap: basis site -> index on the sublattice.  msg is set when the loops span fewer than two directions.
  subroutine clustering_set(b, tau, sh, sb1, sb2, smap, nline, kc, msg)
    real(dp), intent(in) :: b(3, 3), tau(:, :), sh(:, :)
    integer, intent(in) :: sb1(:), sb2(:), smap(:), nline
    real(dp), allocatable, intent(out) :: kc(:, :)
    character(len=:), allocatable, intent(out) :: msg
    integer, parameter :: dmax = 64
    integer, allocatable :: cells(:, :), t(:, :), loops(:, :)
    real(dp), allocatable :: fs(:, :)
    logical, allocatable :: got(:), seen(:, :, :)
    real(dp) :: x(3)
    integer :: nsh, nr, nl, np, j, ia, ib, ic, i1, i2, i3, d, iv(3), c(3), u(3), e, m, q
    logical :: changed
    msg = ''
    nsh = size(sh, 2); nr = maxval(smap)
    allocate (cells(3, nsh), t(3, nr), got(nr), loops(3, nsh))
    do j = 1, nsh
      x = matmul(transpose(b), sh(:, j) - (tau(:, sb2(j)) - tau(:, sb1(j)))) / (2.0_dp * pi)
      cells(:, j) = nint(x)
      if (maxval(abs(x - cells(:, j))) > 1.0e-6_dp) then; msg = 'clustering set: a shell bond is not a lattice bond'; return; end if
    end do
    ! cell offsets t of the sublattice sites along a spanning tree of each connected part of the shell graph
    got = .false.; t = 0
    do ia = 1, nr
      if (got(ia)) cycle
      got(ia) = .true.
      changed = .true.
      do while (changed)
        changed = .false.
        do j = 1, nsh
          if (got(smap(sb1(j))) .and. .not. got(smap(sb2(j)))) then
            t(:, smap(sb2(j))) = t(:, smap(sb1(j))) + cells(:, j); got(smap(sb2(j))) = .true.; changed = .true.
          end if
        end do
      end do
    end do
    nl = 0
    do j = 1, nsh
      iv = t(:, smap(sb1(j))) + cells(:, j) - t(:, smap(sb2(j)))
      if (all(iv == 0)) cycle
      if (nl > 0) then
        if (any([(all(loops(:, m) == iv), m=1, nl)])) cycle
      end if
      nl = nl + 1; loops(:, nl) = iv
    end do
    d = 0
    do ia = 1, nl
      do ib = ia + 1, nl
        do ic = ib + 1, nl
          m = abs(idet3(loops(:, ia), loops(:, ib), loops(:, ic)))
          if (m > 0 .and. (d == 0 .or. m < d)) d = m
        end do
      end do
    end do
    if (d > 0) then
      if (d > dmax) then; msg = 'clustering set: loop lattice index too large'; return; end if
      allocate (fs(3, d**3)); np = 0
      do i1 = 0, d - 1; do i2 = 0, d - 1; do i3 = 0, d - 1
        iv = [i1, i2, i3]
        if (all([(modulo(dot_product(iv, loops(:, m)), d) == 0, m=1, nl)])) then
          np = np + 1; fs(:, np) = real(iv, dp) / real(d, dp)
        end if
      end do; end do; end do
    else
      c = 0
      pair: do ia = 1, nl
        do ib = ia + 1, nl
          c = icross(loops(:, ia), loops(:, ib))
          if (any(c /= 0)) exit pair
        end do
      end do pair
      if (all(c == 0)) then
        msg = 'screening: the SRO shell joins the sublattice in fewer than two directions (complete clustering undefined)'
        return
      end if
      u = c / igcd3(c)
      e = minloc(merge(abs(u), huge(1), u /= 0), 1)
      d = abs(c(e))
      if (d > dmax) then; msg = 'clustering set: loop lattice index too large'; return; end if
      allocate (seen(0:d - 1, 0:d - 1, 0:d - 1), fs(3, d**3 * nline)); seen = .false.; np = 0
      do i1 = 0, d - 1; do i2 = 0, d - 1; do i3 = 0, d - 1
        iv = [i1, i2, i3]
        if (seen(i1, i2, i3)) cycle
        if (.not. all([(modulo(dot_product(iv, loops(:, m)), d) == 0, m=1, nl)])) cycle
        do m = 0, d - 1
          c = modulo(iv + m * u, d); seen(c(1), c(2), c(3)) = .true.
        end do
        do q = 1, nline
          np = np + 1; fs(:, np) = real(iv, dp) / real(d, dp) + (real(q, dp) - 0.5_dp) / real(nline, dp) * real(u, dp)
        end do
      end do; end do; end do
    end if
    kc = matmul(b, fs(:, 1:np))
  contains
    pure integer function idet3(a1, a2, a3)
      integer, intent(in) :: a1(3), a2(3), a3(3)
      idet3 = dot_product(a1, icross(a2, a3))
    end function idet3
    pure function icross(a1, a2) result(cr)
      integer, intent(in) :: a1(3), a2(3)
      integer :: cr(3)
      cr = [a1(2) * a2(3) - a1(3) * a2(2), a1(3) * a2(1) - a1(1) * a2(3), a1(1) * a2(2) - a1(2) * a2(1)]
    end function icross
    pure integer function igcd3(v)
      integer, intent(in) :: v(3)
      integer :: g, k, a, bb, tt
      g = 0
      do k = 1, 3
        a = abs(v(k)); bb = g
        do while (bb /= 0); tt = modulo(a, bb); a = bb; bb = tt; end do
        g = a
      end do
      igcd3 = g
    end function igcd3
  end subroutine clustering_set

  !> Wavevectors kc (Cartesian columns, 1/Bohr) with weights wc that sample the set where the lowest first-shell branch
  !> c1min(k) equals c1_ord (the ordering set).  k0 holds the points known to reach c1_ord (special points).  Through each of
  !> them every rational reciprocal direction u (fractional components -2 .. 2, primitive) is tested: if c1min stays at c1_ord
  !> at four generic points along it, the line k0 + s b u is part of the set.  A line is periodic with period b u.  It is
  !> sampled at nline midpoints of weight |b u| / nline, so the measure is uniform per length.  Lines reached from several
  !> points are counted once.  Without lines the set is the points k0, each of weight 1.  b: reciprocal vectors (columns);
  !> p: primitive vectors (columns, Bohr).  msg is set when the minima fill a plane.
  subroutine ordering_set(b, p, k0, c1min, c1_ord, nline, kc, wc, msg)
    real(dp), intent(in) :: b(3, 3), p(3, 3), k0(:, :), c1_ord
    procedure(branch_fn) :: c1min
    integer, intent(in) :: nline
    real(dp), allocatable, intent(out) :: kc(:, :), wc(:)
    character(len=:), allocatable, intent(out) :: msg
    real(dp), parameter :: sg(4) = [0.1234_dp, 0.3141_dp, 0.4567_dp, 0.7071_dp], tol = 1.0e-9_dp
    integer :: dirs(3, 62), nd, i1, i2, i3, ip, id, nl, nloc, j, q, m, iu, jl
    integer, allocatable :: lu(:, :), loc(:)
    real(dp), allocatable :: lf(:, :)
    real(dp) :: f0(3), d(3), sv, bu(3), bu2(3)
    logical :: on, same
    msg = ''
    nd = 0
    do i1 = -2, 2; do i2 = -2, 2; do i3 = -2, 2
      if (i1 == 0 .and. i2 == 0 .and. i3 == 0) cycle
      if (first_nonzero([i1, i2, i3]) < 0) cycle
      if (gcd3([i1, i2, i3]) /= 1) cycle
      nd = nd + 1; dirs(:, nd) = [i1, i2, i3]
    end do; end do; end do
    allocate (lu(3, size(k0, 2) * nd), lf(3, size(k0, 2) * nd), loc(nd)); nl = 0
    do ip = 1, size(k0, 2)
      f0 = matmul(transpose(p), k0(:, ip)) / (2.0_dp * pi)
      nloc = 0
      do id = 1, nd
        bu = matmul(b, real(dirs(:, id), dp))
        on = .true.
        do j = 1, 4
          if (abs(c1min(k0(:, ip) + sg(j) * bu) - c1_ord) > tol) then; on = .false.; exit; end if
        end do
        if (.not. on) cycle
        nloc = nloc + 1; loc(nloc) = id
        same = .false.
        do jl = 1, nl
          if (any(lu(:, jl) /= dirs(:, id))) cycle
          d = f0 - lf(:, jl)
          iu = maxloc(abs(dirs(:, id)), 1)
          do m = -6, 6
            sv = (d(iu) - real(m, dp)) / real(dirs(iu, id), dp)
            if (maxval(abs(d - sv * dirs(:, id) - anint(d - sv * dirs(:, id)))) < 1.0e-8_dp) then; same = .true.; exit; end if
          end do
          if (same) exit
        end do
        if (.not. same) then
          nl = nl + 1; lu(:, nl) = dirs(:, id); lf(:, nl) = f0
        end if
      end do
      do i1 = 1, nloc
        do i2 = i1 + 1, nloc
          bu = matmul(b, real(dirs(:, loc(i1)), dp)); bu2 = matmul(b, real(dirs(:, loc(i2)), dp))
          if (abs(c1min(k0(:, ip) + 0.2371_dp * bu + 0.3617_dp * bu2) - c1_ord) < tol .and. &
              abs(c1min(k0(:, ip) + 0.4113_dp * bu + 0.1729_dp * bu2) - c1_ord) < tol) then
            msg = 'screening: the first-shell ordering minima fill a plane (not supported)'; return
          end if
        end do
      end do
    end do
    if (nl > 0) then
      allocate (kc(3, nl * nline), wc(nl * nline))
      do jl = 1, nl
        bu = matmul(b, real(lu(:, jl), dp))
        do q = 1, nline
          kc(:, (jl - 1) * nline + q) = matmul(b, lf(:, jl)) + (real(q, dp) - 0.5_dp) / real(nline, dp) * bu
          wc((jl - 1) * nline + q) = norm2(bu) / real(nline, dp)
        end do
      end do
    else
      kc = k0
      allocate (wc(size(k0, 2))); wc = 1.0_dp
    end if
  contains
    pure integer function first_nonzero(v)
      integer, intent(in) :: v(3)
      integer :: k
      first_nonzero = 0
      do k = 1, 3
        if (v(k) /= 0) then; first_nonzero = v(k); return; end if
      end do
    end function first_nonzero
    pure integer function gcd3(v)
      integer, intent(in) :: v(3)
      integer :: g, k, a, bb, tt
      g = 0
      do k = 1, 3
        a = abs(v(k)); bb = g
        do while (bb /= 0); tt = modulo(a, bb); a = bb; bb = tt; end do
        g = a
      end do
      gcd3 = g
    end function gcd3
  end subroutine ordering_set

  !> Screened fluctuation electrostatics on any crystal structure st (asr_crystal, alat in Bohr) for SIM parameter beta.
  !> If the basis sites of st form a single Bravais lattice of translations (one site per primitive cell of that lattice: sc, bcc,
  !> fcc, or any cubic decoration cell), the cell is reduced to that Bravais primitive cell and solved with coulomb_k, reproducing
  !> screen_lattice bit for bit.  For a general multi-site structure (hcp, NiAs, explicit cells), the basis-resolved Ewald matrix
  !> M_bb'(k) (coulomb_k_basis) gives G(k) = a (a I + M(k))^-1 and W(k) = a (I - G(k)).  sub_sro (optional): when > 0, the SRO
  !> shell (z, r1 of the shell, c1(k)) is the nearest shell of sublattice sub_sro with itself (e.g. the anion sublattice of NiAs),
  !> while the SIM parameter beta = -W00 R1_nn / (2 G00) still uses the nearest-neighbour distance R1_nn of the structure.
  subroutine screen_structure(st, alat, beta, s, msg, nk, sub_sro)
    type(structure_t), intent(in) :: st
    real(dp), intent(in) :: alat, beta
    type(screen_t), intent(out) :: s
    character(len=:), allocatable, intent(out) :: msg
    integer, intent(in), optional :: nk, sub_sro
    real(dp) :: al_brav, p_pr(3, 3), k_pr(3, 3), p(3, 3), b(3, 3), f(3), lo, hi, mid, bm, r1_nn, r1_sh
    real(dp) :: ev(st%nb), rwork(3 * st%nb), ph
    real(dp), allocatable :: tau(:, :), kv(:, :), ksp(:, :), sh(:, :), ec(:), rw_c(:)
    real(dp), allocatable :: meig(:, :), w0_eig(:, :), c1_br(:, :), ov(:, :, :), sp_c1(:, :)
    complex(dp), allocatable :: mb(:, :, :), msp(:, :, :), cw(:), ck(:, :), cw_c(:)
    integer, allocatable :: sb1(:), sb2(:), sro_sites(:)
    integer :: smap(st%nb), n, nb, i1, i2, i3, ik, j, m, r, idx, nsh, nsp, ksro, nb_sro, info
    complex(dp) :: vz(st%nb, st%nb)
    real(dp), allocatable :: kcl(:, :), wcl(:), kor(:, :)
    logical :: kzero
    msg = ''
    ksro = 0; if (present(sub_sro)) ksro = sub_sro
    if ((ksro == 0 .or. st%nlabel == 1) .and. match_cubic_bravais(st, p_pr, al_brav, k_pr)) then
      call screen_bravais(p_pr * (alat * al_brav), k_pr * (2.0_dp * pi / (alat * al_brav)), beta, s, msg, nk)
      return
    end if
    n = max(12, nint(24.0_dp / real(st%nb, dp)**(1.0_dp / 3.0_dp)))
    if (present(nk)) n = nk
    nb = st%nb
    p = st%avec * alat
    allocate (tau(3, nb)); tau = st%tau * alat
    b(:, 1) = cross(p(:, 2), p(:, 3)); b(:, 2) = cross(p(:, 3), p(:, 1)); b(:, 3) = cross(p(:, 1), p(:, 2))
    b = 2.0_dp * pi * b / det3(p)
    call structure_nn_distance(p, tau, r1_nn)
    call structure_sro_shell(p, tau, st%label, ksro, sh, sb1, sb2, nsh, nb_sro, r1_sh)
    s%z = nsh / max(nb_sro, 1); s%r1 = r1_nn; s%beta = beta
    if (beta < 0.0_dp) then; msg = 'screening: sim_beta must be >= 0'; return; end if
    if (beta == 0.0_dp) then
      s%unscreened = .true.; s%a = huge(1.0_dp); s%g00 = 1.0_dp; s%g1 = 0.0_dp; s%w00 = 0.0_dp; s%w1 = 2.0_dp / r1_sh
      s%lambda = 1.0_dp; return
    end if
    allocate (sro_sites(nb_sro)); smap = 0; r = 0
    do j = 1, nb
      if (ksro == 0 .or. st%label(j) == ksro) then
        r = r + 1; sro_sites(r) = j; smap(j) = r
      end if
    end do
    allocate (kv(3, n**3), mb(nb, nb, n**3), cw(max(1, 2 * nb)), meig(nb, n**3), w0_eig(nb, n**3))
    allocate (ck(nb_sro, nb_sro), ec(nb_sro), rw_c(max(1, 3 * nb_sro)), cw_c(max(1, 2 * nb_sro)))
    allocate (c1_br(nb_sro, n**3), ov(nb_sro, nb, n**3))
    ik = 0
    do i1 = 0, n - 1; do i2 = 0, n - 1; do i3 = 0, n - 1
      ik = ik + 1; f = ([i1, i2, i3] + 0.5_dp) / n; kv(:, ik) = matmul(b, f)
    end do; end do; end do
    call coulomb_k_basis(p, tau, kv, mb)
    s%a_min = 0.0_dp
    do ik = 1, n**3
      call zheev('V', 'U', nb, mb(:, :, ik), nb, ev, cw, size(cw), rwork, info)
      if (info /= 0) error stop 'screen_structure: zheev failed'
      meig(:, ik) = ev
      s%a_min = max(s%a_min, -minval(ev))
      ck = (0.0_dp, 0.0_dp)
      do j = 1, nsh
        ph = dot_product(kv(:, ik), sh(:, j))
        ck(smap(sb1(j)), smap(sb2(j))) = ck(smap(sb1(j)), smap(sb2(j))) + cmplx(cos(ph), sin(ph), dp) / real(max(s%z, 1), dp)
      end do
      call zheev('V', 'U', nb_sro, ck, nb_sro, ec, cw_c, size(cw_c), rw_c, info)
      if (info /= 0) error stop 'screen_structure: shell zheev failed'
      c1_br(:, ik) = ec
      do m = 1, nb
        do r = 1, nb_sro
          ov(r, m, ik) = abs(dot_product(ck(:, r), mb(sro_sites, m, ik)))**2
        end do
        w0_eig(m, ik) = sum(ov(:, m, ik)) / real(nb_sro, dp)
      end do
    end do
    ! high-symmetry points (zone-boundary points, hexagonal K and H, and for nb > 1 also neutral optical modes at k = 0)
    call zone_special_k(b, nb > 1, ksp, nsp)
    allocate (msp(nb, nb, nsp), sp_c1(nb_sro, nsp))
    call coulomb_k_basis(p, tau, ksp, msp, gamma_neutral=.true.)
    do ik = 1, nsp
      kzero = norm2(ksp(:, ik)) < 1.0e-12_dp
      if (kzero) then
        ! k = 0: neutral modes only (the uniform mode carries the G = 0 divergence)
        call neutral_eig(msp(:, :, ik), ev, vz)
        msp(:, :, ik) = vz
        s%a_min = max(s%a_min, -minval(ev(1:nb - 1)))
      else
        call zheev('V', 'U', nb, msp(:, :, ik), nb, ev, cw, size(cw), rwork, info)
        if (info /= 0) error stop 'screen_structure: zheev failed on special k'
        s%a_min = max(s%a_min, -minval(ev))
      end if
      ck = (0.0_dp, 0.0_dp)
      do j = 1, nsh
        ph = dot_product(ksp(:, ik), sh(:, j))
        ck(smap(sb1(j)), smap(sb2(j))) = ck(smap(sb1(j)), smap(sb2(j))) + cmplx(cos(ph), sin(ph), dp) / real(max(s%z, 1), dp)
      end do
      call zheev('V', 'U', nb_sro, ck, nb_sro, ec, cw_c, size(cw_c), rw_c, info)
      if (info /= 0) error stop 'screen_structure: shell zheev failed on special k'
      sp_c1(:, ik) = ec
    end do
    s%beta_max = beta_of_st(s%a_min * (1.0_dp + 1.0e-9_dp))
    if (beta >= s%beta_max) then
      msg = 'screening: sim_beta exceeds the lattice maximum '//trim(rtoa(s%beta_max))//' (charge-transfer instability)'; return
    end if
    lo = log(s%a_min * 1.0e-9_dp); hi = log(1.0e9_dp)
    do j = 1, 200
      mid = 0.5_dp * (lo + hi); bm = beta_of_st(s%a_min + exp(mid))
      if (bm > beta) then; lo = mid; else; hi = mid; end if
      if (hi - lo < 1.0e-13_dp) exit
    end do
    s%a = s%a_min + exp(0.5_dp * (lo + hi)); s%unscreened = .false.
    allocate (s%gk(nb_sro * n**3), s%c1k(nb_sro * n**3))
    idx = 0
    do ik = 1, n**3
      do r = 1, nb_sro
        idx = idx + 1
        s%c1k(idx) = c1_br(r, ik)
        s%gk(idx) = sum(ov(r, :, ik) * (s%a / (s%a + meig(:, ik))))
      end do
    end do
    s%g00 = sum(s%gk) / real(size(s%gk), dp)
    s%g1 = sum(s%gk * s%c1k) / real(size(s%gk), dp)
    s%w00 = s%a * (1.0_dp - s%g00)
    s%w1 = -s%a * s%g1
    s%lambda = -s%a * s%g1 * r1_sh / (2.0_dp * s%g00**2)
    s%c1_ord = min(minval(s%c1k), minval(sp_c1))
    ! G at the ordering point.  The clipped structure factor ends on the branches with c1 = c1_ord.  Where they form lines
    ! (fcc X-W, close-packed and basal-plane shells), any distribution along them has alpha_1 = c1_ord; the degenerate ground
    ! states of the first-shell model fill them evenly (fcc x = 1/2: random stackings, exact).  So g_ord is the mode average of
    ! g over that set (ordering_set), as g_cl0 is on the clustering side.  At an isolated point it is g of that point.
    if (any(sp_c1(1, :) < s%c1_ord + 1.0e-10_dp)) then
      kor = ksp(:, pack([(ik, ik=1, nsp)], sp_c1(1, :) < s%c1_ord + 1.0e-10_dp))
    else
      kor = kv(:, pack([(ik, ik=1, n**3)], c1_br(1, :) < s%c1_ord + 1.0e-10_dp))
    end if
    call ordering_set(b, p, kor, c1_low, s%c1_ord, 24, kcl, wcl, msg)
    if (len(msg) > 0) return
    s%g_ord = branch_green(kcl, wcl, s%c1_ord)
    ! Complete clustering (alpha_1 = 1): the clipped structure factor ends on the branches with c1 = 1.  It spreads evenly over
    ! them, so the limit is g averaged over that set (clustering_set).  The set is k = 0 alone when the shell joins every site
    ! of the cell in three dimensions; then its only branch is the uniform mode and g_cl0 = 0.  On a sublattice smaller than the
    ! cell the other sublattices screen, and with a shell in the basal plane any stacking of pure planes has c1 = 1: g_cl0 > 0.
    call clustering_set(b, tau, sh, sb1, sb2, smap, 24, kcl, msg)
    if (len(msg) > 0) return
    s%g_cl0 = branch_green(kcl, [(1.0_dp, ik=1, size(kcl, 2))], 1.0_dp)
    s%c0_cl = maxval(s%c1k, mask=s%c1k < maxval(s%c1k) - 1.0e-12_dp)
    s%c0_ord = minval(s%c1k, mask=s%c1k > minval(s%c1k) + 1.0e-12_dp)
    call clip_point(s, s%c0_cl, .true., s%alpha_cl, s%gam_cl)
    call clip_point(s, s%c0_ord, .false., s%alpha_ord, s%gam_ord)
  contains
    !> Mode average of g = v^H G_kk v over the branches with c1 = ct at the wavevectors kc, with weights wc (clustering_set,
    !> ordering_set).  At k = 0 the uniform mode has G = 0.
    real(dp) function branch_green(kc, wc, ct) result(gav)
      real(dp), intent(in) :: kc(:, :), wc(:), ct
      complex(dp) :: mk(nb, nb, 1), uk(nb, nb), gs(nb_sro, nb_sro), cc(nb_sro, nb_sro)
      complex(dp) :: cwa(max(1, 2 * nb)), cwb(max(1, 2 * nb_sro))
      real(dp) :: evk(nb), wk(nb), ecc(nb_sro), rwa(max(1, 3 * nb)), rwb(max(1, 3 * nb_sro)), tr, phk, nd
      integer :: ic, r1, r2, jb, infc, mm
      tr = 0.0_dp; nd = 0.0_dp
      do ic = 1, size(kc, 2)
        call coulomb_k_basis(p, tau, kc(:, ic:ic), mk, gamma_neutral=.true.)
        if (norm2(kc(:, ic)) < 1.0e-12_dp) then
          call neutral_eig(mk(:, :, 1), evk, uk)
          wk = s%a / (s%a + evk); wk(nb) = 0.0_dp
        else
          uk = mk(:, :, 1)
          call zheev('V', 'U', nb, uk, nb, evk, cwa, size(cwa), rwa, infc)
          if (infc /= 0) error stop 'screen_structure: zheev failed on the clustering set'
          wk = s%a / (s%a + evk)
        end if
        do r2 = 1, nb_sro
          do r1 = 1, nb_sro
            gs(r1, r2) = sum(uk(sro_sites(r1), :) * wk * conjg(uk(sro_sites(r2), :)))
          end do
        end do
        cc = (0.0_dp, 0.0_dp)
        do jb = 1, nsh
          phk = dot_product(kc(:, ic), sh(:, jb))
          cc(smap(sb1(jb)), smap(sb2(jb))) = cc(smap(sb1(jb)), smap(sb2(jb))) + cmplx(cos(phk), sin(phk), dp) / real(s%z, dp)
        end do
        call zheev('V', 'U', nb_sro, cc, nb_sro, ecc, cwb, size(cwb), rwb, infc)
        if (infc /= 0) error stop 'screen_structure: shell zheev failed on the clustering set'
        do mm = 1, nb_sro
          if (abs(ecc(mm) - ct) < 1.0e-9_dp) then
            tr = tr + wc(ic) * real(dot_product(cc(:, mm), matmul(gs, cc(:, mm))), dp); nd = nd + wc(ic)
          end if
        end do
      end do
      if (nd <= 0.0_dp) error stop 'screen_structure: no branch at the target c1 on the sampled set'
      gav = tr / nd
    end function branch_green

    !> Lowest eigenvalue of the shell matrix C(k).
    real(dp) function c1_low(k)
      real(dp), intent(in) :: k(3)
      complex(dp) :: cl(nb_sro, nb_sro), cwl(max(1, 2 * nb_sro))
      real(dp) :: el(nb_sro), rwl(max(1, 3 * nb_sro)), phl
      integer :: jl, infl
      cl = (0.0_dp, 0.0_dp)
      do jl = 1, nsh
        phl = dot_product(k, sh(:, jl))
        cl(smap(sb1(jl)), smap(sb2(jl))) = cl(smap(sb1(jl)), smap(sb2(jl))) + cmplx(cos(phl), sin(phl), dp) / real(s%z, dp)
      end do
      call zheev('N', 'U', nb_sro, cl, nb_sro, el, cwl, size(cwl), rwl, infl)
      if (infl /= 0) error stop 'screen_structure: shell zheev failed'
      c1_low = minval(el)
    end function c1_low

    real(dp) function beta_of_st(a)
      real(dp), intent(in) :: a
      real(dp) :: g
      g = sum(w0_eig * (a / (a + meig))) / real(n**3, dp)
      beta_of_st = a * (g - 1.0_dp) * s%r1 / (2.0_dp * g)
    end function beta_of_st
  end subroutine screen_structure
  !> Screened fluctuation electrostatics on a single-site Bravais lattice with primitive vectors p (columns, Bohr) and
  !> high-symmetry ordering wavevectors ksp (columns, 1/Bohr).
  subroutine screen_bravais(p, ksp, beta, s, msg, nk)
    real(dp), intent(in) :: p(3, 3), ksp(:, :), beta
    type(screen_t), intent(out) :: s
    character(len=:), allocatable, intent(out) :: msg
    integer, intent(in), optional :: nk
    real(dp) :: b(3, 3), f(3), lo, hi, mid, bm, msp(size(ksp, 2)), c1sp(size(ksp, 2))
    real(dp), allocatable :: kv(:, :), m(:), c1(:), sh(:, :), kor(:, :), wor(:), mor(:)
    integer :: n, i1, i2, i3, ik, j, nsh, nsp
    msg = ''
    n = 28; if (present(nk)) n = nk
    b(:, 1) = cross(p(:, 2), p(:, 3)); b(:, 2) = cross(p(:, 3), p(:, 1)); b(:, 3) = cross(p(:, 1), p(:, 2))
    b = 2.0_dp * pi * b / det3(p)
    call first_shell(p, sh, nsh); s%z = nsh; s%r1 = norm2(sh(:, 1)); s%beta = beta
    if (beta < 0.0_dp) then; msg = 'screening: sim_beta must be >= 0'; return; end if
    if (beta == 0.0_dp) then                    ! unscreened limit
      s%unscreened = .true.; s%a = huge(1.0_dp); s%g00 = 1.0_dp; s%g1 = 0.0_dp; s%w00 = 0.0_dp; s%w1 = 2.0_dp / s%r1
      s%lambda = 1.0_dp; return
    end if
    allocate (kv(3, n**3), m(n**3), c1(n**3))
    ik = 0
    do i1 = 0, n - 1; do i2 = 0, n - 1; do i3 = 0, n - 1
          ik = ik + 1; f = ([i1, i2, i3] + 0.5_dp) / n; kv(:, ik) = matmul(b, f)
        end do; end do; end do
    call coulomb_k(p, kv, m)
    do ik = 1, n**3
      c1(ik) = 0.0_dp
      do j = 1, nsh; c1(ik) = c1(ik) + cos(dot_product(kv(:, ik), sh(:, j))); end do
      c1(ik) = c1(ik) / nsh
    end do
    s%a_min = -minval(m)
    nsp = size(ksp, 2)
    if (nsp > 0) then
      call coulomb_k(p, ksp, msp)
      s%a_min = max(s%a_min, -minval(msp))
    end if
    s%beta_max = beta_of(s%a_min * (1.0_dp + 1.0e-9_dp))
    if (beta >= s%beta_max) then
      msg = 'screening: sim_beta exceeds the lattice maximum '//trim(rtoa(s%beta_max))//' (charge-transfer instability)'; return
    end if
    lo = log(s%a_min * 1.0e-9_dp); hi = log(1.0e9_dp)       ! bisection on log(a - a_min); beta decreases with a
    do j = 1, 200
      mid = 0.5_dp * (lo + hi); bm = beta_of(s%a_min + exp(mid))
      if (bm > beta) then; lo = mid; else; hi = mid; end if
      if (hi - lo < 1.0e-13_dp) exit
    end do
    s%a = s%a_min + exp(0.5_dp * (lo + hi)); s%unscreened = .false.
    s%g00 = sum(s%a / (s%a + m)) / n**3
    s%g1 = sum(s%a / (s%a + m) * c1) / n**3
    s%w00 = s%a * (1.0_dp - s%g00)
    s%w1 = -s%a * s%g1
    s%lambda = -s%a * s%g1 * s%r1 / (2.0_dp * s%g00**2)
    ! closure: lattice functions on the mesh, the ordering set, the last clipped points (next to the extreme
    ! c1 value of the mesh, so that S keeps support)
    s%gk = s%a / (s%a + m); s%c1k = c1
    do j = 1, nsp; c1sp(j) = sum(cos(matmul(ksp(:, j), sh))) / nsh; end do
    s%c1_ord = minval(c1sp)
    ! G at the ordering point: the average of a / (a + M) over the set where c1 = c1_ord (ordering_set; screen_structure).
    ! bcc H and sc R are isolated points, so their value does not change; fcc averages over the X-W lines.
    call ordering_set(b, p, ksp(:, pack([(j, j=1, nsp)], c1sp < s%c1_ord + 1.0e-10_dp)), c1_brav, s%c1_ord, 24, kor, wor, msg)
    if (len(msg) > 0) return
    allocate (mor(size(kor, 2)))
    call coulomb_k(p, kor, mor)
    s%g_ord = sum(wor * (s%a / (s%a + mor))) / sum(wor)
    s%c0_cl = maxval(c1, mask=c1 < maxval(c1) - 1.0e-12_dp)
    s%c0_ord = minval(c1, mask=c1 > minval(c1) + 1.0e-12_dp)
    call clip_point(s, s%c0_cl, .true., s%alpha_cl, s%gam_cl)
    call clip_point(s, s%c0_ord, .false., s%alpha_ord, s%gam_ord)
  contains
    real(dp) function beta_of(a)
      real(dp), intent(in) :: a
      real(dp) :: g
      g = sum(a / (a + m)) / size(m)
      beta_of = a * (g - 1.0_dp) * s%r1 / (2.0_dp * g)
    end function beta_of
    real(dp) function c1_brav(k)
      real(dp), intent(in) :: k(3)
      c1_brav = sum(cos(matmul(k, sh))) / nsh
    end function c1_brav
  end subroutine screen_bravais

  subroutine screen_lattice(lattice, alat, beta, s, msg, nk)
    character(len=*), intent(in) :: lattice
    real(dp), intent(in) :: alat, beta
    type(screen_t), intent(out) :: s
    character(len=:), allocatable, intent(out) :: msg
    integer, intent(in), optional :: nk
    real(dp) :: p_pr(3, 3), k_pr(3, 3)
    call prototype_bravais_data(lattice, p_pr, k_pr, msg)
    if (len(msg) > 0) then; msg = 'screening: '//msg; return; end if
    call screen_bravais(p_pr * alat, k_pr * (2.0_dp * pi / alat), beta, s, msg, nk)
  end subroutine screen_lattice

  !> Lattice functions of the charge-excess model on a Bravais lattice p (columns, Bohr) at lattice vectors rv (columns, Bohr):
  !>     gr = G(R) = <G(k) cos(k.R)>_BZ,     wr = W(R) = a (delta_R0 - G(R)).
  subroutine screen_bravais_shells(p, s, rv, gr, wr, nk)
    real(dp), intent(in) :: p(3, 3), rv(:, :)
    type(screen_t), intent(in) :: s
    real(dp), intent(out) :: gr(:), wr(:)
    integer, intent(in), optional :: nk
    real(dp) :: b(3, 3), f(3)
    real(dp), allocatable :: kv(:, :), m(:)
    integer :: n, i1, i2, i3, ik, j
    logical :: zero
    if (s%unscreened) then
      do j = 1, size(rv, 2)
        zero = norm2(rv(:, j)) < 1.0e-10_dp
        gr(j) = merge(1.0_dp, 0.0_dp, zero)
        wr(j) = 0.0_dp; if (.not. zero) wr(j) = 2.0_dp / norm2(rv(:, j))
      end do
      return
    end if
    n = 28; if (present(nk)) n = nk
    b(:, 1) = cross(p(:, 2), p(:, 3)); b(:, 2) = cross(p(:, 3), p(:, 1)); b(:, 3) = cross(p(:, 1), p(:, 2))
    b = 2.0_dp * pi * b / det3(p)
    allocate (kv(3, n**3), m(n**3))
    ik = 0
    do i1 = 0, n - 1; do i2 = 0, n - 1; do i3 = 0, n - 1
          ik = ik + 1; f = ([i1, i2, i3] + 0.5_dp) / n; kv(:, ik) = matmul(b, f)
        end do; end do; end do
    call coulomb_k(p, kv, m)
    m = s%a / (s%a + m)
    do j = 1, size(rv, 2)
      gr(j) = sum(m * cos(matmul(rv(:, j), kv))) / n**3
      wr(j) = -s%a * gr(j)
      if (norm2(rv(:, j)) < 1.0e-10_dp) wr(j) = s%a * (1.0_dp - gr(j))
    end do
  end subroutine screen_bravais_shells

  subroutine screen_shells(lattice, alat, s, rv, gr, wr, nk)
    character(len=*), intent(in) :: lattice
    real(dp), intent(in) :: alat, rv(:, :)
    type(screen_t), intent(in) :: s
    real(dp), intent(out) :: gr(:), wr(:)
    integer, intent(in), optional :: nk
    real(dp) :: p_pr(3, 3), k_pr(3, 3)
    character(len=:), allocatable :: msg
    call prototype_bravais_data(lattice, p_pr, k_pr, msg)
    if (len(msg) > 0) error stop 'screen_shells: unknown prototype lattice'
    call screen_bravais_shells(p_pr * alat, s, rv, gr, wr, nk)
  end subroutine screen_shells

  !> Lattice functions G(R) and W(R) of the charge-excess model on any crystal structure st at displacement vectors rv(3, npair)
  !> (Bohr) between origin basis site b0 (optional; default 1) and partner basis sites b2(1:npair) (optional; default 1).
  !> For a single-type cubic Bravais structure it uses the matched prototype cell through screen_bravais_shells.
  subroutine screen_structure_shells(st, alat, s, rv, gr, wr, b0, b2, nk)
    type(structure_t), intent(in) :: st
    real(dp), intent(in) :: alat, rv(:, :)
    type(screen_t), intent(in) :: s
    real(dp), intent(out) :: gr(:), wr(:)
    integer, intent(in), optional :: b0, b2(:), nk
    real(dp) :: al_brav, p_pr(3, 3), k_pr(3, 3), p(3, 3), b(3, 3), f(3), ev(st%nb), rwork(3 * st%nb), ph
    real(dp), allocatable :: tau(:, :), kv(:, :)
    complex(dp), allocatable :: mb(:, :, :), cw(:)
    complex(dp) :: gk_el
    integer :: n, nb, i1, i2, i3, ik, j, m, info, b_orig, b_part
    logical :: zero
    if (s%unscreened) then
      do j = 1, size(rv, 2)
        zero = norm2(rv(:, j)) < 1.0e-10_dp
        gr(j) = merge(1.0_dp, 0.0_dp, zero)
        wr(j) = 0.0_dp; if (.not. zero) wr(j) = 2.0_dp / norm2(rv(:, j))
      end do
      return
    end if
    if (match_cubic_bravais(st, p_pr, al_brav, k_pr)) then
      call screen_bravais_shells(p_pr * (alat * al_brav), s, rv, gr, wr, nk)
      return
    end if
    n = max(12, nint(24.0_dp / real(st%nb, dp)**(1.0_dp / 3.0_dp)))
    if (present(nk)) n = nk
    nb = st%nb
    p = st%avec * alat
    allocate (tau(3, nb)); tau = st%tau * alat
    b(:, 1) = cross(p(:, 2), p(:, 3)); b(:, 2) = cross(p(:, 3), p(:, 1)); b(:, 3) = cross(p(:, 1), p(:, 2))
    b = 2.0_dp * pi * b / det3(p)
    allocate (kv(3, n**3), mb(nb, nb, n**3), cw(max(1, 2 * nb)))
    ik = 0
    do i1 = 0, n - 1; do i2 = 0, n - 1; do i3 = 0, n - 1
      ik = ik + 1; f = ([i1, i2, i3] + 0.5_dp) / n; kv(:, ik) = matmul(b, f)
    end do; end do; end do
    call coulomb_k_basis(p, tau, kv, mb)
    b_orig = 1; if (present(b0)) b_orig = b0
    gr = 0.0_dp
    do ik = 1, n**3
      call zheev('V', 'U', nb, mb(:, :, ik), nb, ev, cw, size(cw), rwork, info)
      if (info /= 0) error stop 'screen_structure_shells: zheev failed'
      do j = 1, size(rv, 2)
        b_part = b_orig; if (present(b2)) b_part = b2(j)
        gk_el = (0.0_dp, 0.0_dp)
        do m = 1, nb
          gk_el = gk_el + mb(b_orig, m, ik) * conjg(mb(b_part, m, ik)) * (s%a / (s%a + ev(m)))
        end do
        ph = dot_product(kv(:, ik), rv(:, j))
        gr(j) = gr(j) + real(gk_el * cmplx(cos(ph), -sin(ph), dp), dp)
      end do
    end do
    gr = gr / real(n**3, dp)
    do j = 1, size(rv, 2)
      wr(j) = -s%a * gr(j)
      if (norm2(rv(:, j)) < 1.0e-10_dp) wr(j) = s%a * (1.0_dp - gr(j))
    end do
  end subroutine screen_structure_shells

  !> Screened fluctuation electrostatics for the component matrix ahat (n x n) with the weights w (w_p = f_kappa x_p).
  !> gm = Gamma (dq = Gamma db), wm = a (I - Gamma) (V = wm db), fm = F = a (Gamma^-1 - I) (V = F dq), lam / gam = eigenvalues
  !> of Ahat and of Gamma (one pair per eigenmode), ginv = pseudo-inverse of Gamma (modes with Gamma = 0 left out).
  !> Ahat = 0 (no SRO) and beta = 0 (unscreened): the linear form G00 + G1 Ahat, bit for bit as before.  Otherwise the closure
  !> (screen_gamma) on every eigenmode of the w-symmetric Ahat.  fm with a mode Gamma = 0 (full clustering) stops: use screen_update.
  subroutine screen_matrix(s, ahat, w, fm, gm, wm, lam, gam, ginv)
    type(screen_t), intent(in) :: s
    real(dp), intent(in) :: ahat(:, :), w(:)
    real(dp), intent(out), optional :: fm(:, :), gm(:, :), wm(:, :), lam(:), gam(:), ginv(:, :)
    real(dp) :: g(size(ahat, 1), size(ahat, 1)), f(size(ahat, 1), size(ahat, 1)), y(size(ahat, 1), size(ahat, 1))
    real(dp) :: d(size(ahat, 1)), ev(size(ahat, 1)), ge(size(ahat, 1)), work(3 * size(ahat, 1) + 1)
    integer :: n, i, j, info, ipiv(size(ahat, 1))
    n = size(ahat, 1)
    if (s%unscreened .or. all(ahat == 0.0_dp)) then
      g = s%g1 * ahat; f = s%w1 * ahat
      do i = 1, n; g(i, i) = g(i, i) + s%g00; f(i, i) = f(i, i) + s%w00; end do
      if (present(gm)) gm = g
      if (present(wm)) wm = f
      if (present(lam) .or. present(gam)) then
        call w_sym(ahat, w, d, y)
        call dsyev('N', 'U', n, y, n, ev, work, size(work), info)
        if (info /= 0) error stop 'screen_matrix: dsyev failed'
        if (present(lam)) lam = ev
        if (present(gam)) gam = s%g00 + s%g1 * ev
      end if
      if (present(ginv)) then
        y = g; ginv = 0.0_dp
        do i = 1, n; ginv(i, i) = 1.0_dp; end do
        call dgesv(n, n, y, n, ipiv, ginv, n, info)
        if (info /= 0) error stop 'screen_matrix: singular G00 + G1 Ahat'
      end if
      if (present(fm)) then                          ! as before the closure: F = (G00 + G1 Ahat)^-1 (W00 + W1 Ahat)
        call dgesv(n, n, g, n, ipiv, f, n, info)
        if (info /= 0) error stop 'screen_matrix: singular G00 + G1 Ahat'
        fm = f
      end if
      return
    end if
    ! closure: eigenmodes of D Ahat D^-1 (symmetric, D = diag(sqrt(w))), Gamma = D^-1 Q diag(gamma(lambda)) Q^T D
    call w_sym(ahat, w, d, y)
    call dsyev('V', 'U', n, y, n, ev, work, size(work), info)
    if (info /= 0) error stop 'screen_matrix: dsyev failed'
    do i = 1, n; ge(i) = screen_gamma(s, ev(i)); end do
    if (present(lam)) lam = ev
    if (present(gam)) gam = ge
    if (present(gm)) gm = back(ge)
    if (present(wm)) wm = back(s%a * (1.0_dp - ge))
    if (present(ginv)) ginv = back(merge(1.0_dp / max(ge, tiny(1.0_dp)), 0.0_dp, ge > 1.0e-9_dp * s%g00))
    if (present(fm)) then
      if (any(ge <= 1.0e-9_dp * s%g00)) &
        error stop 'screen_matrix: F is infinite (an SRO mode with Gamma = 0, full clustering); use the drive form screen_update'
      fm = back(s%a * (1.0_dp / ge - 1.0_dp))
    end if
  contains
    function back(v) result(mat)
      real(dp), intent(in) :: v(:)
      real(dp) :: mat(n, n)
      do j = 1, n
        do i = 1, n
          mat(i, j) = sum(y(i, :) * v * y(j, :)) * d(j) / d(i)
        end do
      end do
    end function back
  end subroutine screen_matrix

  !> y = D ahat D^-1 with D = diag(sqrt(w)): symmetric for consistent SRO input (diag(w) Ahat symmetric); checked, then symmetrised.
  subroutine w_sym(ahat, w, d, y)
    real(dp), intent(in) :: ahat(:, :), w(:)
    real(dp), intent(out) :: d(:), y(:, :)
    integer :: i, j
    if (any(w <= 0.0_dp)) error stop 'screen_matrix: component weights must be > 0'
    d = sqrt(w)
    do j = 1, size(w); do i = 1, size(w); y(i, j) = d(i) * ahat(i, j) / d(j); end do; end do
    if (maxval(abs(y - transpose(y))) > 1.0e-6_dp * max(1.0_dp, maxval(abs(y)))) &
      error stop 'screen_matrix: diag(w) Ahat is not symmetric (inconsistent first-shell SRO input)'
    y = 0.5_dp * (y + transpose(y))
  end subroutine w_sym

  !> Gamma of one SRO eigenmode with Ahat eigenvalue lam (first-shell alpha_1 = lam / z): the minimal-tail closure (module
  !> header).  -1 <= lam <= -1/c1(k_o): the linear form G00 + G1 lam, exactly.  Clustering beyond: S ~ max(0, c1 - c0), c0 by
  !> bisection between c1(k_o) and c0_cl, then the condensate at k = 0.  Ordering beyond: S ~ max(0, c0 - c1), c0 between c0_ord
  !> and 1, then the condensate at k_o.  Unscreened (beta = 0): 1.
  function screen_gamma(s, lam) result(gam)
    type(screen_t), intent(in) :: s
    real(dp), intent(in) :: lam
    real(dp) :: gam, al, lo, hi, mid, a1, m
    logical :: up
    integer :: j
    if (s%unscreened) then; gam = 1.0_dp; return; end if
    al = lam / s%z
    if (lam >= -1.0_dp .and. lam <= -1.0_dp / s%c1_ord) then
      gam = s%g00 + s%g1 * lam
      return
    end if
    up = lam > 0.0_dp
    if (up .and. al >= s%alpha_cl) then
      m = min(1.0_dp, (al - s%alpha_cl) / (1.0_dp - s%alpha_cl))
      gam = (1.0_dp - m) * s%gam_cl + m * s%g_cl0
      return
    else if (.not. up .and. al <= s%alpha_ord) then
      ! m > 1 only for alpha_1 below c1(k_o) (fcc: not realizable on one sublattice); G(k_o) >= gam_ord keeps Gamma > 0
      m = (s%alpha_ord - al) / max(s%alpha_ord - s%c1_ord, 1.0e-12_dp)
      gam = (1.0_dp - m) * s%gam_ord + m * s%g_ord
      return
    end if
    if (up) then; lo = s%c1_ord; hi = s%c0_cl; else; lo = s%c0_ord; hi = 1.0_dp; end if
    do j = 1, 200
      mid = 0.5_dp * (lo + hi)
      call clip_point(s, mid, up, a1, gam)
      if (a1 < al) then; lo = mid; else; hi = mid; end if
      if (hi - lo < 1.0e-15_dp) exit
    end do
    call clip_point(s, 0.5_dp * (lo + hi), up, a1, gam)
  end function screen_gamma

  !> F of one SRO mode (Ry/e) from its Ahat eigenvalue lam and its Gamma gam (callers check gam > 0):
  !> a (1/gam - 1); unscreened (beta = 0): W00 + W1 lam = (2/R1) lam, the bare shell charge.
  real(dp) function mode_f(s, lam, gam)
    type(screen_t), intent(in) :: s
    real(dp), intent(in) :: lam, gam
    if (s%unscreened) then
      mode_f = s%w00 + s%w1 * lam
    else
      mode_f = s%a * (1.0_dp / gam - 1.0_dp)
    end if
  end function mode_f

  !> Clipped first-shell structure factor at c0 on the mesh: S ~ max(0, c1 - c0) (up) or max(0, c0 - c1); returns
  !> alpha_1 = <c1 S> / <S> and Gamma = <G S> / <S>.
  subroutine clip_point(s, c0, up, al, gam)
    type(screen_t), intent(in) :: s
    real(dp), intent(in) :: c0
    logical, intent(in) :: up
    real(dp), intent(out) :: al, gam
    real(dp) :: w(size(s%c1k)), sw
    if (up) then; w = max(0.0_dp, s%c1k - c0); else; w = max(0.0_dp, c0 - s%c1k); end if
    sw = sum(w)
    al = sum(s%c1k * w) / sw
    gam = sum(s%gk * w) / sw
  end subroutine clip_point

  !> One SCF step of the SRO electrostatics in the drive u = db of the charge-excess model:
  !>     u_out = P [dq + (I - Gamma) u_in],   v = W u_out,   W = a (I - Gamma)   (gm = Gamma, wm = W from screen_matrix).
  !> At the fixed point u_out = u_in: dq = Gamma u and v = W u, so v = F dq wherever F is finite, and dq = 0 on a mode with
  !> Gamma = 0 (F infinite).  P removes the x-weighted mean of every sublattice (grp: sublattice of each component, x: its
  !> concentration there); dq has no such mean, and Gamma keeps the deviation subspace, so P only stops round-off drift.
  pure subroutine screen_update(gm, wm, x, grp, dq, u_in, u_out, v)
    real(dp), intent(in) :: gm(:, :), wm(:, :), x(:), dq(:), u_in(:)
    integer, intent(in) :: grp(:)
    real(dp), intent(out) :: u_out(:), v(:)
    u_out = dq + u_in - matmul(gm, u_in)
    call screen_project(x, grp, u_out)
    v = matmul(wm, u_out)
  end subroutine screen_update

  !> u -> u minus its x-weighted mean on every sublattice (grp(p): sublattice of component p).
  pure subroutine screen_project(x, grp, u)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: grp(:)
    real(dp), intent(inout) :: u(:)
    integer :: k
    do k = 1, maxval(grp)
      u = u - merge(sum(merge(x * u, 0.0_dp, grp == k)), 0.0_dp, grp == k)
    end do
  end subroutine screen_project

  !> First-shell vectors of the Bravais lattice (columns of sh, Bohr).
  subroutine first_shell(p, sh, nsh)
    real(dp), intent(in) :: p(3, 3)
    real(dp), allocatable, intent(out) :: sh(:, :)
    integer, intent(out) :: nsh
    real(dp) :: v(3), r, rmin, buf(3, 64)
    integer :: i1, i2, i3
    rmin = huge(1.0_dp)
    do i1 = -2, 2; do i2 = -2, 2; do i3 = -2, 2
          v = i1 * p(:, 1) + i2 * p(:, 2) + i3 * p(:, 3); r = norm2(v)
          if (r > 1.0e-10_dp) rmin = min(rmin, r)
        end do; end do; end do
    nsh = 0
    do i1 = -2, 2; do i2 = -2, 2; do i3 = -2, 2
          v = i1 * p(:, 1) + i2 * p(:, 2) + i3 * p(:, 3)
          if (abs(norm2(v) - rmin) < 1.0e-8_dp * rmin) then; nsh = nsh + 1; buf(:, nsh) = v; end if
        end do; end do; end do
    sh = buf(:, 1:nsh)
  end subroutine first_shell

  !> Shortest inter-site distance in the crystal structure (p, tau in Bohr).
  subroutine structure_nn_distance(p, tau, r1)
    real(dp), intent(in) :: p(3, 3), tau(:, :)
    real(dp), intent(out) :: r1
    real(dp) :: v(3), r
    integer :: nb, b1, b2, i1, i2, i3
    nb = size(tau, 2); r1 = huge(1.0_dp)
    do b1 = 1, nb
      do b2 = 1, nb
        do i1 = -2, 2; do i2 = -2, 2; do i3 = -2, 2
          if (b1 == b2 .and. i1 == 0 .and. i2 == 0 .and. i3 == 0) cycle
          v = i1 * p(:, 1) + i2 * p(:, 2) + i3 * p(:, 3) + tau(:, b2) - tau(:, b1)
          r = norm2(v)
          if (r > 1.0e-10_dp) r1 = min(r1, r)
        end do; end do; end do
      end do
    end do
  end subroutine structure_nn_distance

  !> Nearest coordination shell of the structure (ksro = 0) or of sublattice ksro with itself (ksro > 0): returns all directed
  !> bond vectors sh(3, nsh) from basis sites sb1(j) to sb2(j) in the primitive cell, the number of origin basis sites nb_sro,
  !> and the shell radius r1_sh (Bohr).
  subroutine structure_sro_shell(p, tau, label, ksro, sh, sb1, sb2, nsh, nb_sro, r1_sh)
    real(dp), intent(in) :: p(3, 3), tau(:, :)
    integer, intent(in) :: label(:), ksro
    real(dp), allocatable, intent(out) :: sh(:, :)
    integer, allocatable, intent(out) :: sb1(:), sb2(:)
    integer, intent(out) :: nsh, nb_sro
    real(dp), intent(out) :: r1_sh
    real(dp), allocatable :: buf(:, :)
    integer, allocatable :: ib1(:), ib2(:)
    real(dp) :: v(3), r
    integer :: nb, b1, b2, i1, i2, i3
    nb = size(tau, 2)
    r1_sh = huge(1.0_dp); nb_sro = 0
    do b1 = 1, nb
      if (ksro > 0 .and. label(b1) /= ksro) cycle
      nb_sro = nb_sro + 1
      do b2 = 1, nb
        if (ksro > 0 .and. label(b2) /= ksro) cycle
        do i1 = -2, 2; do i2 = -2, 2; do i3 = -2, 2
          if (b1 == b2 .and. i1 == 0 .and. i2 == 0 .and. i3 == 0) cycle
          v = i1 * p(:, 1) + i2 * p(:, 2) + i3 * p(:, 3) + tau(:, b2) - tau(:, b1)
          r = norm2(v)
          if (r > 1.0e-10_dp) r1_sh = min(r1_sh, r)
        end do; end do; end do
      end do
    end do
    allocate (buf(3, 128 * nb), ib1(128 * nb), ib2(128 * nb)); nsh = 0
    do b1 = 1, nb
      if (ksro > 0 .and. label(b1) /= ksro) cycle
      do b2 = 1, nb
        if (ksro > 0 .and. label(b2) /= ksro) cycle
        do i1 = -2, 2; do i2 = -2, 2; do i3 = -2, 2
          v = i1 * p(:, 1) + i2 * p(:, 2) + i3 * p(:, 3) + tau(:, b2) - tau(:, b1)
          if (abs(norm2(v) - r1_sh) < 1.0e-7_dp * r1_sh) then
            nsh = nsh + 1; buf(:, nsh) = v; ib1(nsh) = b1; ib2(nsh) = b2
          end if
        end do; end do; end do
      end do
    end do
    sh = buf(:, 1:nsh); sb1 = ib1(1:nsh); sb2 = ib2(1:nsh)
  end subroutine structure_sro_shell

  !> High-symmetry k points for a general reciprocal cell b (columns, 1/Bohr).  Covers the inversion-inequivalent points of the
  !> quarter-zone grid {0, 1/4, 1/2, 3/4}^3 (zone-boundary points M, L, A, X, R, H, N, plus quarter-zone ordering points such as
  !> fcc W and close-packed (1/4, 1/4, 0)) and the hexagonal 1/3 points K and H for every axis orientation.  Includes k = 0 when
  !> add_gamma is true so multi-site cells sample neutral optical modes.
  subroutine zone_special_k(b, add_gamma, ksp, nsp)
    real(dp), intent(in) :: b(3, 3)
    logical, intent(in) :: add_gamma
    real(dp), allocatable, intent(out) :: ksp(:, :)
    integer, intent(out) :: nsp
    real(dp), parameter :: t1 = 1.0_dp / 3.0_dp, t2 = 2.0_dp / 3.0_dp, h = 0.5_dp
    real(dp) :: f(3, 48), v(3)
    integer :: i1, i2, i3, j1, j2, j3, id1, id2, ax, s2, s3, k
    k = 0
    do i1 = 0, 3; do i2 = 0, 3; do i3 = 0, 3
      if (i1 == 0 .and. i2 == 0 .and. i3 == 0) then
        if (add_gamma) then; k = k + 1; f(:, k) = 0.0_dp; end if
        cycle
      end if
      j1 = modulo(-i1, 4); j2 = modulo(-i2, 4); j3 = modulo(-i3, 4)
      id1 = 16 * i1 + 4 * i2 + i3; id2 = 16 * j1 + 4 * j2 + j3
      if (id1 <= id2) then
        k = k + 1; f(:, k) = 0.25_dp * real([i1, i2, i3], dp)
      end if
    end do; end do; end do
    do ax = 1, 3
      do s2 = 1, 2
        do s3 = 0, 1
          v(ax) = s3 * h
          v(modulo(ax, 3) + 1) = t1
          v(modulo(ax + 1, 3) + 1) = merge(t1, t2, s2 == 1)
          k = k + 1; f(:, k) = v
        end do
      end do
    end do
    nsp = k
    allocate (ksp(3, nsp))
    ksp = matmul(b, f(:, 1:nsp))
  end subroutine zone_special_k

  pure function cross(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
  end function cross

  pure real(dp) function det3(p)
    real(dp), intent(in) :: p(3, 3)
    det3 = dot_product(p(:, 1), cross(p(:, 2), p(:, 3)))
  end function det3

  function rtoa(x) result(s)
    real(dp), intent(in) :: x
    character(len=:), allocatable :: s
    character(len=32) :: buf
    write (buf, '(f8.4)') x; s = trim(adjustl(buf))
  end function rtoa
end module asr_screening
