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
!> E = 1/2 sum_p w_p dq_p (F dq)_p has dE/dq_p = w_p (F dq)_p; within the model it is the exact configuration average.
module asr_screening
  use asr_kinds, only: dp, pi
  implicit none
  private
  public :: screen_t, screen_lattice, screen_matrix, coulomb_k

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
    real(dp) :: p(3, 3), b(3, 3), f(3), lo, hi, mid, bm
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
    s%a_min = max(s%a_min, -minval(m_special(lattice, alat, p)))
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
  contains
    real(dp) function beta_of(a)
      real(dp), intent(in) :: a
      real(dp) :: g
      g = sum(a / (a + m)) / size(m)
      beta_of = a * (g - 1.0_dp) * s%r1 / (2.0_dp * g)
    end function beta_of
  end subroutine screen_lattice

  !> F = (G00 I + G1 Ahat)^-1 (W00 I + W1 Ahat) for the component matrix ahat (n x n).
  subroutine screen_matrix(s, ahat, fm)
    type(screen_t), intent(in) :: s
    real(dp), intent(in) :: ahat(:, :)
    real(dp), intent(out) :: fm(size(ahat, 1), size(ahat, 1))
    real(dp) :: gm(size(ahat, 1), size(ahat, 1))
    integer :: n, i, info, ipiv(size(ahat, 1))
    n = size(ahat, 1)
    gm = s%g1 * ahat; fm = s%w1 * ahat
    do i = 1, n; gm(i, i) = gm(i, i) + s%g00; fm(i, i) = fm(i, i) + s%w00; end do
    call dgesv(n, n, gm, n, ipiv, fm, n, info)
    if (info /= 0) error stop 'screen_matrix: singular G00 + G1 Ahat'
  end subroutine screen_matrix

  !> M at the zone-boundary high-symmetry points, where the lattice Coulomb kernel takes its lowest values (ordering vectors):
  !> bcc H, N, P; fcc X, W, L; sc R, M, X (Cartesian, 2 pi / alat units).
  function m_special(lattice, alat, p) result(ms)
    character(len=*), intent(in) :: lattice
    real(dp), intent(in) :: alat, p(3, 3)
    real(dp) :: ms(3), kv(3, 3)
    select case (lattice)
    case ('bcc'); kv = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp], [3, 3])
    case ('fcc'); kv = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp], [3, 3])
    case default; kv = reshape([0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp], [3, 3])
    end select
    kv = kv * 2.0_dp * pi / alat
    call coulomb_k(p, kv, ms)
  end function m_special

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
