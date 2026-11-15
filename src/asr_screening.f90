!> Screened electrostatics of charge fluctuations in a substitutional alloy, single-site form, with first-shell short-range
!> order (gate B3; derivation in docs/FINDINGS_LOG.md, "B3 re-derived").
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
!> Closure of the first-shell SRO (docs/FINDINGS_LOG.md, "B3 closure"; scripts in tools/sro_closure).  The linear form
!> Gamma = G00 + G1 Ahat keeps the correlations of the first shell only.  Its structure factor turns negative at the zone
!> boundary for strong clustering, and G00 + G1 lambda then has a zero (fcc beta 1.0: lambda* = 12.35; beta 1.2: 7.67).  The
!> exact model has no zero: Gamma = <G(k) S(k)> / x with G > 0 and S >= 0, and sum_R G(R) = 0 gives Gamma -> 0 (F -> +infinity)
!> for full segregation.  The code therefore completes the structure factor of every eigenmode of Ahat (eigenvalue lambda,
!> first-shell alpha_1 = lambda / z) with as few correlations beyond the first shell as possible (minimal-tail closure), the
!> configuration model of the ASR itself and of the SQS cells of gate B3:
!>   -1 <= lambda <= -1/c1min (fcc -1..3, bcc -1..1): S = 1 + lambda c1(k) >= 0 is realizable; Gamma = G00 + G1 lambda exactly.
!>   lambda > -1/c1min (clustering): S ~ max(0, c1(k) - c0), c0 from c1min upwards; S concentrates at k = 0, where G = 0.
!>   lambda < -1 (ordering):          S ~ max(0, c0 - c1(k)), c0 from 1 downwards; S concentrates at the ordering point k_o.
!> c0 follows from alpha_1 = <c1 S> / <S> (bisection; it increases with c0).  Past the last mesh value the order condenses at
!> k = 0 (Gamma -> 0 at alpha_1 = 1) or at k_o, the point of lowest M where c1 is minimal (bcc H: Gamma -> G(H), exact for
!> perfect B2 order; fcc W: on fcc alpha_1 = -1/3 holds on the whole X-W line and does not fix the order; sc R).  Gamma is positive, decreasing and continuous on the whole realizable range.  Tests (tools/sro_closure): exact model on
!> minimal-tail supercells (inverse Monte Carlo): Gamma within 1 % at beta 1.0 and 1.2 where the linear form is off by 25 %;
!> CuZn alpha -0.148 unchanged.  On thermodynamic (long-tail) supercells the closure is within 4 % at beta 1.0; at beta 1.2 the
!> answer depends on the tails, which the first-shell input does not fix.  No SRO and beta = 0 keep the linear form.
!>
!> A mode with Gamma = 0 has F = +infinity: its charge deviation must vanish.  V = F dq cannot represent it, but the drive
!> u = db can: dq = Gamma u, V = a (1 - Gamma) u (screen_update).  This form has the same fixed point wherever F is finite.
module asr_screening
  use asr_lapack, only: dgesv, dsyev
  use asr_kinds, only: dp, pi
  implicit none
  private
  public :: screen_t, screen_lattice, screen_matrix, screen_gamma, screen_update, screen_project, mode_f, coulomb_k, screen_shells

  type :: screen_t
    real(dp) :: beta = 0.0_dp            !! input SIM parameter
    real(dp) :: a = huge(1.0_dp)         !! hardness (Ry per electron); huge = unscreened
    real(dp) :: a_min = 0.0_dp           !! stability bound of a
    real(dp) :: beta_max = 0.0_dp        !! largest beta of the lattice
    real(dp) :: r1 = 0.0_dp              !! nearest-neighbour distance (Bohr)
    integer :: z = 0                     !! first-shell coordination of the graph lattice
    real(dp) :: g00 = 1.0_dp, g1 = 0.0_dp, w00 = 0.0_dp, w1 = 0.0_dp
    real(dp) :: lambda = 1.0_dp          !! linear first-shell weight -a G1 R1 / (2 G00^2) (1 = bare shell)
    logical :: unscreened = .true.
    ! closure of the first-shell SRO (minimal-tail structure factor on the BZ mesh)
    real(dp), allocatable :: gk(:), c1k(:)   !! G(k) and c1(k) on the mesh
    real(dp) :: g_ord = 0.0_dp, c1_ord = -1.0_dp   !! G and c1 at the ordering point k_o
    real(dp) :: c0_cl = 1.0_dp, alpha_cl = 1.0_dp, gam_cl = 0.0_dp      !! last clipped point on the mesh, clustering side
    real(dp) :: c0_ord = -1.0_dp, alpha_ord = -1.0_dp, gam_ord = 0.0_dp !! last clipped point on the mesh, ordering side
  end type screen_t

contains

  !> Primitive vectors (columns, alat units) of the graph lattice.
  subroutine prim_vectors(lattice, p, msg)
    character(len=*), intent(in) :: lattice
    real(dp), intent(out) :: p(3, 3)
    character(len=:), allocatable, intent(out) :: msg
    msg = ''
    select case (lattice)
    case ('sc'); p = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3]) * 1.0_dp
    case ('bcc'); p = reshape([-1, 1, 1, 1, -1, 1, 1, 1, -1], [3, 3]) * 0.5_dp
    case ('fcc'); p = reshape([0, 1, 1, 1, 0, 1, 1, 1, 0], [3, 3]) * 0.5_dp
    case default; msg = 'screening: unknown graph lattice '//lattice
    end select
  end subroutine prim_vectors

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

  !> Lattice functions of the graph lattice (bcc, fcc, sc) with constant alat (Bohr) for the SIM parameter beta: solves the
  !> hardness a and fills G00, G1, W00, W1, lambda.  nk: BZ mesh per direction (offset mesh; default 28).  msg empty on success.
  subroutine screen_lattice(lattice, alat, beta, s, msg, nk)
    character(len=*), intent(in) :: lattice
    real(dp), intent(in) :: alat, beta
    type(screen_t), intent(out) :: s
    character(len=:), allocatable, intent(out) :: msg
    integer, intent(in), optional :: nk
    real(dp) :: p(3, 3), b(3, 3), f(3), lo, hi, mid, bm, ksp(3, 3), msp(3), c1sp(3)
    real(dp), allocatable :: kv(:, :), m(:), c1(:), sh(:, :)
    integer :: n, i1, i2, i3, ik, j, nsh
    n = 28; if (present(nk)) n = nk
    call prim_vectors(lattice, p, msg); if (len(msg) > 0) return
    p = p * alat
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
    ! the exact minimum of M lies at a high-symmetry point (bcc H, fcc X, sc R) that the offset mesh misses
    ksp = special_k(lattice, alat)
    call coulomb_k(p, ksp, msp)
    s%a_min = max(s%a_min, -minval(msp))
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
    ! closure: lattice functions on the mesh, the ordering point k_o, the last clipped points (next to the extreme
    ! c1 value of the mesh, so that S keeps support)
    s%gk = s%a / (s%a + m); s%c1k = c1
    ! ordering point: lowest M among the special points where c1 is minimal (bcc H, fcc W on the degenerate X-W line, sc R)
    do j = 1, 3; c1sp(j) = sum(cos(matmul(ksp(:, j), sh))) / nsh; end do
    s%c1_ord = minval(c1sp)
    s%g_ord = s%a / (s%a + minval(msp, mask=c1sp < s%c1_ord + 1.0e-10_dp))
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
  end subroutine screen_lattice

  !> Lattice functions of the charge-excess model at lattice vectors rv (columns, Bohr; R = 0 allowed):
  !>     gr = G(R) = <G(k) cos(k.R)>_BZ,     wr = W(R) = a (delta_R0 - G(R)).
  !> W is the pair kernel of the model energy.  The model minimises sum_i a/2 (q_i - b_i)^2 + 1/2 sum_ij q_i M_ij q_j over the
  !> charges; the minimum is E = 1/2 sum_ij b_i W_ij b_j with the drives b (species properties).  So two species with drives
  !> b_p, b_q on sites R apart have the screened Coulomb exchange interaction W(R) (b_p - b_q)^2 (asr_epi).  In the random alloy
  !> the drive deviations are db = dq / G00.  Unscreened (beta = 0): G = delta_R0, W = 2 / R.  Mesh as in screen_lattice.
  subroutine screen_shells(lattice, alat, s, rv, gr, wr, nk)
    character(len=*), intent(in) :: lattice
    real(dp), intent(in) :: alat, rv(:, :)
    type(screen_t), intent(in) :: s
    real(dp), intent(out) :: gr(:), wr(:)
    integer, intent(in), optional :: nk
    real(dp) :: p(3, 3), b(3, 3), f(3)
    real(dp), allocatable :: kv(:, :), m(:)
    character(len=:), allocatable :: msg
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
    call prim_vectors(lattice, p, msg); if (len(msg) > 0) error stop 'screen_shells: unknown graph lattice'
    p = p * alat
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
  end subroutine screen_shells

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
      gam = (1.0_dp - m) * s%gam_cl
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

  !> Zone-boundary high-symmetry points, where the lattice Coulomb kernel takes its lowest values (ordering vectors):
  !> bcc H, N, P; fcc X, W, L; sc R, M, X (Cartesian, 1/Bohr).  The ordering point k_o is the one of lowest M among those with
  !> minimal c1 (bcc H, fcc W, sc R).
  function special_k(lattice, alat) result(kv)
    character(len=*), intent(in) :: lattice
    real(dp), intent(in) :: alat
    real(dp) :: kv(3, 3)
    select case (lattice)
    case ('bcc'); kv = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp], [3, 3])
    case ('fcc'); kv = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp], [3, 3])
    case default; kv = reshape([0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp], [3, 3])
    end select
    kv = kv * 2.0_dp * pi / alat
  end function special_k

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
