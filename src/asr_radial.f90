!> Scalar-relativistic (Koelling-Harmon) radial solver on a shifted log mesh.
!> Rydberg units: hbar^2/2m = 1, e^2 = 2, c = 274.0719...
!>
!>   P = r g(r)  (large component),  Q = r (1/2M)(g' )*... defined by
!>   P' =  M Q + P/r
!>   Q' = -Q/r + [ l(l+1)/(M r^2) + V - E ] P ,   M = 1 + (E - V)/c^2
!> (c -> infinity gives the Schroedinger equation, P'' = [l(l+1)/r^2 + V - E] P).
!> Norm:  int (P^2 + Q^2/c^2) dr.   Energy derivatives P_k = d^k P/dE^k (k<=2)
!> are integrated exactly from the variational equations.
module asr_radial
  use asr_kinds, only: dp, pi
  implicit none
  private
  public :: radmesh_t, make_radmesh, integrate_regular, integrate_inward, node_count, bound_state, hankel_logder, norm_pq, c_ry, c_nonrel

  real(dp), parameter :: c_ry = 274.071979_dp       !! speed of light in Rydberg atomic units
  real(dp), parameter :: c_nonrel = 1.0e10_dp

  !> r_i = b (exp(a (i-1)) - 1),  dr/dt = a (r + b)
  type :: radmesh_t
    integer :: n = 0
    real(dp) :: a = 0.0_dp, b = 0.0_dp
    real(dp), allocatable :: r(:), drdt(:)
  end type radmesh_t

contains

  !> Mesh from the origin to rmax.  The first nonzero point is placed well inside
  !> the region r << 2Z/c^2 where the relativistic Coulomb behaviour holds.
  subroutine make_radmesh(z, rmax, a, m, c)
    real(dp), intent(in) :: z, rmax, a
    type(radmesh_t), intent(out) :: m
    real(dp), intent(in), optional :: c
    real(dp) :: cc, b, r2
    integer :: n, i
    cc = c_ry; if (present(c)) cc = c
    r2 = 0.05_dp * 2.0_dp * max(z, 1.0_dp) / (cc * cc)       ! first nonzero radius
    r2 = max(r2, 1.0e-7_dp)
    b = r2 / (exp(a) - 1.0_dp)
    n = ceiling(log(rmax / b + 1.0_dp) / a) + 1
    if (mod(n, 2) == 0) n = n + 1
    m%n = n; m%a = a
    m%b = rmax / (exp(a * (n - 1)) - 1.0_dp)                  ! land exactly on rmax
    ! keep the first-point scale: recompute b so r(n)=rmax; r2 changes slightly
    allocate (m%r(n), m%drdt(n))
    do i = 1, n
      m%r(i) = m%b * (exp(a * (i - 1)) - 1.0_dp)
      m%drdt(i) = a * (m%r(i) + m%b)
    end do
    m%r(1) = 0.0_dp
  end subroutine make_radmesh

  !> Midpoint of r*V on [i,i+1] from the surrounding nodes (cubic).
  pure function rv_mid(rv, i, n) result(u)
    real(dp), intent(in) :: rv(:)
    integer, intent(in) :: i, n
    real(dp) :: u
    if (i >= 2 .and. i + 2 <= n) then
      u = (-rv(i - 1) + 9.0_dp * rv(i) + 9.0_dp * rv(i + 1) - rv(i + 2)) / 16.0_dp
    else if (i + 2 <= n) then
      u = (5.0_dp * rv(i) + 15.0_dp * rv(i + 1) - 5.0_dp * rv(i + 2) + rv(i + 3)) / 16.0_dp
    else
      u = (rv(i - 2) - 5.0_dp * rv(i - 1) + 15.0_dp * rv(i) + 5.0_dp * rv(i + 1)) / 16.0_dp
    end if
  end function rv_mid

  !> Right-hand side of the extended system at (r, rvv = r V): y(1:2,k) = (P_k, Q_k), k=0..nd
  pure subroutine rhs(l, e, c, r, rvv, nd, y, dy)
    integer, intent(in) :: l, nd
    real(dp), intent(in) :: e, c, r, rvv
    real(dp), intent(in) :: y(2, 0:nd)
    real(dp), intent(out) :: dy(2, 0:nd)
    real(dp) :: v, m, a21, a21d, a21dd, ll, c2
    integer :: k
    v = rvv / r; c2 = c * c; ll = real(l * (l + 1), dp)
    m = 1.0_dp + (e - v) / c2
    a21 = ll / (m * r * r) + v - e
    dy(1, 0) = m * y(2, 0) + y(1, 0) / r
    dy(2, 0) = -y(2, 0) / r + a21 * y(1, 0)
    if (nd >= 1) then
      a21d = -ll / (m * m * r * r * c2) - 1.0_dp
      dy(1, 1) = m * y(2, 1) + y(1, 1) / r + y(2, 0) / c2
      dy(2, 1) = -y(2, 1) / r + a21 * y(1, 1) + a21d * y(1, 0)
    end if
    if (nd >= 2) then
      a21dd = 2.0_dp * ll / (m**3 * r * r * c2 * c2)
      dy(1, 2) = m * y(2, 2) + y(1, 2) / r + 2.0_dp * y(2, 1) / c2
      dy(2, 2) = -y(2, 2) / r + a21 * y(1, 2) + 2.0_dp * a21d * y(1, 1) + a21dd * y(1, 0)
    end if
    k = 0
  end subroutine rhs

  !> Integrate the regular solution outward from the origin to point `itop`.
  !> v: true potential on the mesh (incl. -2Z/r).  y(2,0:nd,i) = (P_k,Q_k).
  !> Unnormalised (P ~ r^(l+1) at the origin).
  subroutine integrate_regular(m, v, l, e, z, c, nd, itop, y)
    type(radmesh_t), intent(in) :: m
    real(dp), intent(in) :: v(:), e, z, c
    integer, intent(in) :: l, nd, itop
    real(dp), intent(out) :: y(2, 0:nd, itop)
    real(dp), allocatable :: rv(:)
    real(dp) :: h, r0, rm, um, k1(2, 0:nd), k2(2, 0:nd), k3(2, 0:nd), k4(2, 0:nd), yy(2, 0:nd)
    real(dp) :: s, m0, c1, vv
    integer :: i, n
    n = m%n
    allocate (rv(n))
    rv(2:n) = m%r(2:n) * v(2:n)
    rv(1) = 3.0_dp * rv(2) - 3.0_dp * rv(3) + rv(4)          ! quadratic extrapolation to r=0
    h = 1.0_dp                                   ! step in the mesh index s (dr/ds = a (r+b))
    y = 0.0_dp
    ! ---- start at i = 2 : Frobenius leading terms
    r0 = m%r(2); vv = v(2)
    m0 = 1.0_dp + (e - vv) / (c * c)
    if (m0 > 3.0_dp) then
      ! relativistic Coulomb region: P = r^s, Q = ((s-1) c^2 / (2Z)) r^s
      s = sqrt(real(l * (l + 1), dp) + 1.0_dp - (2.0_dp * z / c)**2)
      y(1, 0, 2) = r0**s
      y(2, 0, 2) = (s - 1.0_dp) * c * c / (2.0_dp * z) * r0**s
    else
      ! non-relativistic start with first Coulomb correction
      c1 = -z / real(l + 1, dp)
      y(1, 0, 2) = r0**(l + 1) * (1.0_dp + c1 * r0)
      y(2, 0, 2) = r0**l * (real(l, dp) + c1 * real(l + 1, dp) * r0) / m0
    end if
    ! derivative components start at zero (energy enters at higher order in r)
    do i = 2, itop - 1
      yy = y(:, :, i)
      call rhs(l, e, c, m%r(i), rv(i), nd, yy, k1)
      k1 = k1 * m%drdt(i)
      rm = m%b * (exp(m%a * (i - 0.5_dp)) - 1.0_dp)
      um = rv_mid(rv, i, n)
      call rhs(l, e, c, rm, um, nd, yy + 0.5_dp * h * k1, k2)
      k2 = k2 * (m%a * (rm + m%b))
      call rhs(l, e, c, rm, um, nd, yy + 0.5_dp * h * k2, k3)
      k3 = k3 * (m%a * (rm + m%b))
      call rhs(l, e, c, m%r(i + 1), rv(i + 1), nd, yy + h * k3, k4)
      k4 = k4 * m%drdt(i + 1)
      y(:, :, i + 1) = yy + h * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4) / 6.0_dp
      ! guard against overflow of the growing solution
      if (abs(y(1, 0, i + 1)) > 1.0e100_dp) y(:, :, i + 1:) = 0.0_dp
      if (abs(y(1, 0, i + 1)) > 1.0e100_dp) exit
    end do
  end subroutine integrate_regular

  !> Integrate inward from index `istart` (decaying start P=1, Q/P=qp) down to index `iend`.
  !> y(2,i) for i = iend..istart (unnormalised, no energy derivatives).
  subroutine integrate_inward(m, v, l, e, c, istart, iend, qp, y)
    type(radmesh_t), intent(in) :: m
    real(dp), intent(in) :: v(:), e, c, qp
    integer, intent(in) :: l, istart, iend
    real(dp), intent(out) :: y(2, iend:istart)
    real(dp), allocatable :: rv(:)
    real(dp) :: rm, um, k1(2, 0:0), k2(2, 0:0), k3(2, 0:0), k4(2, 0:0), yy(2, 0:0)
    integer :: i, n
    n = m%n
    allocate (rv(n)); rv(2:n) = m%r(2:n) * v(2:n); rv(1) = 3.0_dp * rv(2) - 3.0_dp * rv(3) + rv(4)
    y(1, istart) = 1.0_dp; y(2, istart) = qp
    do i = istart, iend + 1, -1
      yy(:, 0) = y(:, i)
      call rhs(l, e, c, m%r(i), rv(i), 0, yy, k1); k1 = k1 * m%drdt(i)
      rm = m%b * (exp(m%a * (i - 1.5_dp)) - 1.0_dp)
      um = rv_mid(rv, i - 1, n)
      call rhs(l, e, c, rm, um, 0, yy - 0.5_dp * k1, k2); k2 = k2 * (m%a * (rm + m%b))
      call rhs(l, e, c, rm, um, 0, yy - 0.5_dp * k2, k3); k3 = k3 * (m%a * (rm + m%b))
      call rhs(l, e, c, m%r(i - 1), rv(i - 1), 0, yy - k3, k4); k4 = k4 * m%drdt(i - 1)
      y(:, i - 1) = yy(:, 0) - (k1(:, 0) + 2.0_dp * k2(:, 0) + 2.0_dp * k3(:, 0) + k4(:, 0)) / 6.0_dp
    end do
  end subroutine integrate_inward

  !> Number of sign changes of P in (0, r_itop].
  function node_count(m, v, l, e, z, c, itop) result(nn)
    type(radmesh_t), intent(in) :: m
    real(dp), intent(in) :: v(:), e, z, c
    integer, intent(in) :: l, itop
    integer :: nn, i
    real(dp), allocatable :: y(:, :, :)
    allocate (y(2, 0:0, itop))
    call integrate_regular(m, v, l, e, z, c, 0, itop, y)
    nn = 0
    do i = 3, itop
      if (y(1, 0, i) == 0.0_dp) cycle
      if (y(1, 0, i) * y(1, 0, i - 1) < 0.0_dp) nn = nn + 1
    end do
  end function node_count

  !> Logarithmic derivative  r P'/P - 1 ... returned as Q/P at radius s for a decaying
  !> Hankel tail in a constant potential vout:  Q/P = D_l / (M s),  D_l = x k_l'/k_l
  function hankel_logder(l, kappa, s, e, vout, c) result(qp)
    integer, intent(in) :: l
    real(dp), intent(in) :: kappa, s, e, vout, c
    real(dp) :: qp, x, k0, k1, k2, dl, mm
    integer :: i
    x = kappa * s
    ! scaled k_l(x) e^x : k0 = 1/x, k1 = 1/x + 1/x^2 (common factor dropped)
    k0 = 1.0_dp / x; k1 = 1.0_dp / x + 1.0_dp / x**2
    if (l == 0) then
      dl = -1.0_dp - x            ! x k0'/k0 with k0 = e^-x / x
    else
      do i = 1, l - 1
        k2 = k0 + (2 * i + 1) / x * k1
        k0 = k1; k1 = k2
      end do
      ! now k0 = k_{l-1}, k1 = k_l
      dl = -x * k0 / k1 - (l + 1)
    end if
    mm = 1.0_dp + (e - vout) / (c * c)
    qp = dl / (mm * s)
  end function hankel_logder

  !> Normalise (P,Q,...) so that int_0^S (P^2 + Q^2/c^2) dr = 1 and make the energy
  !> derivatives those of the normalised family (so <phi|phidot> = 0).
  !> y(2,0:nd,1:itop) modified in place; returns raw norm.
  subroutine norm_pq(m, c, nd, itop, y, nrm)
    type(radmesh_t), intent(in) :: m
    real(dp), intent(in) :: c
    integer, intent(in) :: nd, itop
    real(dp), intent(inout) :: y(2, 0:nd, itop)
    real(dp), intent(out) :: nrm
    real(dp) :: f(itop), n0, n1, n2, s0
    integer :: i
    ! raw integrals via Simpson in t on the first itop points (itop odd)
    f = [(y(1, 0, i)**2 + y(2, 0, i)**2 / c**2, i=1, itop)] * m%drdt(1:itop)
    n0 = simpson_t(f, 1.0_dp, itop)
    nrm = n0
    if (nd >= 1) then
      f = [(2.0_dp * (y(1, 0, i) * y(1, 1, i) + y(2, 0, i) * y(2, 1, i) / c**2), i=1, itop)] * m%drdt(1:itop)
      n1 = simpson_t(f, 1.0_dp, itop)
    else
      n1 = 0.0_dp
    end if
    if (nd >= 2) then
      f = [(2.0_dp * (y(1, 1, i)**2 + y(1, 0, i) * y(1, 2, i) + (y(2, 1, i)**2 + y(2, 0, i) * y(2, 2, i)) / c**2), i=1, itop)] &
          * m%drdt(1:itop)
      n2 = simpson_t(f, 1.0_dp, itop)
    else
      n2 = 0.0_dp
    end if
    s0 = 1.0_dp / sqrt(n0)
    ! phi = y0 N^-1/2;  phidot = y1 N^-1/2 - y0 N' N^-3/2 / 2 ;
    ! phiddot = y2 N^-1/2 - y1 N' N^-3/2 + y0 (3 N'^2/(4 N^2) - N''/(2 N)) N^-1/2
    if (nd >= 2) then
      y(:, 2, :) = y(:, 2, :) * s0 - y(:, 1, :) * n1 * s0**3 + y(:, 0, :) * s0 * (0.75_dp * (n1 / n0)**2 - 0.5_dp * n2 / n0)
    end if
    if (nd >= 1) then
      y(:, 1, :) = y(:, 1, :) * s0 - 0.5_dp * y(:, 0, :) * n1 * s0**3
    end if
    y(:, 0, :) = y(:, 0, :) * s0
  end subroutine norm_pq

  pure function simpson_t(f, h, n) result(s)
    real(dp), intent(in) :: f(:), h
    integer, intent(in) :: n
    real(dp) :: s
    integer :: i
    s = f(1) + f(n)
    do i = 2, n - 1
      if (mod(i, 2) == 0) then; s = s + 4.0_dp * f(i); else; s = s + 2.0_dp * f(i); end if
    end do
    s = s * h / 3.0_dp
  end function simpson_t

  !> Bound state with `nod` nodes in the sphere of mesh index `isph`; potential outside the
  !> sphere is constant = v(isph) (decaying Hankel boundary condition).  Eigenvalue found by
  !> matching the outward and inward solutions at the outer turning point (both stable),
  !> bracketed by node counting.  On output y is normalised on [0, r(isph)].
  subroutine bound_state(m, v, l, nod, z, c, isph, elo, ehi, e, y, ok)
    type(radmesh_t), intent(in) :: m
    real(dp), intent(in) :: v(:), z, c, elo, ehi
    integer, intent(in) :: l, nod, isph
    real(dp), intent(out) :: e
    real(dp), intent(out) :: y(2, 0:0, isph)
    logical, intent(out), optional :: ok     !! .false. if no bracket was found (best estimate returned)
    logical :: found
    real(dp) :: e1, e2, ea, eb, fa, fb, fm, em, vout, s, ehi_eff, ll
    integer :: it, ntot, nc
    real(dp), allocatable :: yo(:, :, :), yi(:, :)
    allocate (yo(2, 0:0, isph), yi(2, 2:isph))
    found = .true.
    s = m%r(isph); vout = v(isph); ll = real(l * (l + 1), dp)
    ehi_eff = min(ehi, vout - 1.0e-8_dp)
    e1 = bisect_nodes(nod)
    e2 = bisect_nodes(nod + 1)
    ea = e1; eb = e2
    call eval(ea, fa, ntot, nc); call eval(eb, fb, ntot, nc)
    if (fa * fb > 0.0_dp) call scan_bracket(ea, eb, fa, fb)
    if (present(ok)) ok = found
    do it = 1, 200
      em = 0.5_dp * (ea + eb)
      call eval(em, fm, ntot, nc)
      if (fa * fm <= 0.0_dp) then; eb = em; fb = fm; else; ea = em; fa = fm; end if
      if (abs(eb - ea) < 1.0e-13_dp * max(1.0_dp, abs(em))) exit
    end do
    e = 0.5_dp * (ea + eb)
    ! final wavefunction: outward to nc, inward beyond, matched in P
    call eval(e, fm, ntot, nc)
    y(:, 0, 1:nc) = yo(:, 0, 1:nc)
    block
      real(dp) :: sc, nrm
      integer :: j
      sc = yo(1, 0, nc) / yi(1, nc)
      if (abs(yi(1, nc)) < 1.0e-300_dp) sc = 1.0_dp
      do j = nc + 1, isph
        y(1, 0, j) = sc * yi(1, j); y(2, 0, j) = sc * yi(2, j)
      end do
      call norm_pq(m, c, 0, isph, y, nrm)
    end block
  contains
    subroutine eval(en, f, nt, nmatch)
      real(dp), intent(in) :: en
      real(dp), intent(out) :: f
      integer, intent(out) :: nt, nmatch
      real(dp) :: kap, qp, fo, fi, cumk, kloc
      integer :: j, it2, istart
      logical :: usewkb
      it2 = 0
      do j = 3, isph
        if (en - v(j) - ll / m%r(j)**2 > 0.0_dp) it2 = j
      end do
      nmatch = min(it2 + 10, isph - 1)
      if (it2 == 0) nmatch = min(isph - 1, 5)
      call integrate_regular(m, v, l, en, z, c, 0, nmatch, yo(:, :, 1:nmatch))
      ! where to start the inward integration: Hankel boundary condition at the sphere edge if the
      ! decay is weak; otherwise a local WKB decaying start where the accumulated decay is ~ e^-35
      cumk = 0.0_dp; istart = isph; usewkb = .false.
      do j = nmatch + 1, isph
        kloc = sqrt(max(v(j) + ll / m%r(j)**2 - en, 0.0_dp))
        cumk = cumk + kloc * (m%r(j) - m%r(j - 1))
        if (cumk > 35.0_dp) then; istart = j; usewkb = .true.; exit; end if
      end do
      if (usewkb) then
        kap = sqrt(max(v(istart) + ll / m%r(istart)**2 - en, 1.0e-12_dp))
        qp = -(kap + 1.0_dp / m%r(istart)) / (1.0_dp + (en - v(istart)) / c**2)
      else
        kap = sqrt(max(vout - en, 1.0e-12_dp))
        qp = hankel_logder(l, kap, s, en, vout, c)
      end if
      yi = 0.0_dp
      call integrate_inward(m, v, l, en, c, istart, nmatch, qp, yi(:, nmatch:istart))
      nt = 0
      do j = 3, nmatch
        if (yo(1, 0, j) * yo(1, 0, j - 1) < 0.0_dp) nt = nt + 1
      end do
      do j = nmatch + 1, isph
        if (yi(1, j) * yi(1, j - 1) < 0.0_dp) nt = nt + 1
      end do
      ! scale-safe: normalise each (P,Q) pair by its own magnitude before forming the cross product
      fo = max(abs(yo(1, 0, nmatch)), abs(yo(2, 0, nmatch)))
      fi = max(abs(yi(1, nmatch)), abs(yi(2, nmatch)))
      f = (yo(2, 0, nmatch) / fo) * (yi(1, nmatch) / fi) - (yi(2, nmatch) / fi) * (yo(1, 0, nmatch) / fo)
    end subroutine eval
    function bisect_nodes(target) result(eb_)
      integer, intent(in) :: target
      real(dp) :: eb_, lo, hi, mid, f
      integer :: k, nt, nmatch
      lo = elo; hi = ehi_eff
      do k = 1, 200
        mid = 0.5_dp * (lo + hi)
        call eval(mid, f, nt, nmatch)
        if (nt >= target) then; hi = mid; else; lo = mid; end if
        if (hi - lo < 1.0e-12_dp) exit
      end do
      eb_ = hi
    end function bisect_nodes
    subroutine scan_bracket(ea_, eb_, fa_, fb_)
      real(dp), intent(inout) :: ea_, eb_, fa_, fb_
      real(dp) :: x0, x1, f0, f1
      integer :: k, nt, nmatch
      x0 = ea_; call eval(x0, f0, nt, nmatch)
      do k = 1, 400
        x1 = ea_ + (eb_ - ea_) * k / 400.0_dp
        call eval(x1, f1, nt, nmatch)
        if (f0 * f1 <= 0.0_dp) then; ea_ = x0; eb_ = x1; fa_ = f0; fb_ = f1; return; end if
        x0 = x1; f0 = f1
      end do
      found = .false.
      ea_ = eb_; fa_ = fb_                       ! best estimate: upper end of the node window
      if (.not. present(ok)) error stop 'bound_state: no sign change of the matching function in the node window'
    end subroutine scan_bracket
  end subroutine bound_state
end module asr_radial
