!> Point-group (D4h, 16 operations) reduced multi-orbital augmented-space operator, for a start orbital among
!> s, p_z, t2g(xy), e_g(3z^2-r^2) at the origin.  Same Hamiltonian as asr_augorb (H = E + h - h o h with random
!> c, d, o, E_nu) but on canonical representatives of the (site, configuration) nodes; each node keeps its 9 orbital
!> components.  Symmetry of the invariant subspace (H commutes with U_g, (U_g psi)(g node) = D(g) psi(node), and
!> the start vector has U_g u1 = chi(g) u1, chi(g) = D(g)[L0,L0]) gives, for a non-representative node f = g^-1 r,
!>      psi(f) = chi(g) D(g)^T psi(r).
!> The vector is stored in normalised coordinates w = sqrt(n_node) psi (n_node = orbit size) so that the standard
!> Lanczos recursion applies.  Links (flip, hop) store the representative index and the operation g.
module asr_augsym
  use asr_kinds, only: dp
  use asr_operator, only: linop
  use asr_augspace, only: graph_t, store_t, store_init, find_or_add, find_state, toggled, is_down
  use asr_augorb, only: orb_alloy_t
  implicit none
  private
  public :: augsym_op, build_augsym_operator, augsym_set_params, augsym_set_start, augsym_ops, augsym_hph

  integer, parameter :: nop = 16

  type, extends(linop) :: augsym_op
    integer :: nn = 0, nl = 9, nslot = 0
    real(dp) :: x = 0.5_dp, s = 0.5_dp, y = 0.5_dp
    real(dp), allocatable :: c0(:), c1(:), d0(:), d1(:), o0(:), o1(:), e0(:), e1(:), block(:, :, :), block0(:, :)
    real(dp), allocatable :: p0(:), p1(:)               !! p^gamma (order 3)
    integer :: order = 2                                !! 2: E + h - hoh;  3: + hohoh - (hphh + hhph)/2
    real(dp), allocatable :: btab(:, :, :, :)           !! btab(l, l2, slot, g): S block times D(g)^T (chi applied at run time)
    integer, allocatable :: flipr(:), flipg(:), hopr(:, :), hopg(:, :), level(:)
    real(dp), allocatable :: md(:), norb(:), sq(:)
    integer :: pi(9, nop), start_l = 1
    real(dp) :: sg(9, nop), chi(nop)
    real(dp), pointer :: wk(:, :) => null()
  contains
    procedure :: apply => augsym_apply
    procedure :: dim => augsym_dim
  end type augsym_op

contains

  pure integer function augsym_dim(self)
    class(augsym_op), intent(in) :: self
    augsym_dim = self%nl * self%nn
  end function augsym_dim

  !> choose the start orbital (1 = s, 3 = p_z, 5 = xy, 7 = 3z^2-r^2) : sets chi(g) = D(g)[L0,L0]
  subroutine augsym_set_start(op, l0, dmat)
    type(augsym_op), intent(inout) :: op
    integer, intent(in) :: l0
    real(dp), intent(in) :: dmat(9, 9, nop)
    integer :: g
    op%start_l = l0
    do g = 1, nop; op%chi(g) = dmat(l0, l0, g); end do
  end subroutine augsym_set_start

  subroutine augsym_set_params(op, mdl)
    type(augsym_op), intent(inout) :: op
    type(orb_alloy_t), intent(in) :: mdl
    op%c0 = mdl%c0; op%c1 = mdl%c1; op%d0 = mdl%d0; op%d1 = mdl%d1; op%o0 = mdl%o0; op%o1 = mdl%o1
    op%e0 = mdl%e0; op%e1 = mdl%e1
    if (allocated(mdl%p0)) then; op%p0 = mdl%p0; op%p1 = mdl%p1; end if
  end subroutine augsym_set_params

  ! ------------------------------------------------------------------ application
  !> (get psi at the node reached by a link) = chi(g) sg(L,g) v(pi(L,g)) from the representative's vector v
  pure subroutine at_link(op, v, r, g, out)
    class(augsym_op), intent(in) :: op
    real(dp), intent(in) :: v(:)
    integer, intent(in) :: r, g
    real(dp), intent(out) :: out(9)
    integer :: l
    do l = 1, 9
      out(l) = op%chi(g) * op%sg(l, g) * v(9 * (r - 1) + op%pi(l, g))
    end do
  end subroutine at_link

  subroutine augsym_apply(self, x, y)
    class(augsym_op), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    integer :: i, l
    if (.not. associated(self%wk)) error stop 'augsym_apply: work space not allocated'
    ! normalised -> amplitude coordinates
    !$omp parallel do private(l)
    do i = 1, self%nn
      do l = 1, 9; self%wk(9 * (i - 1) + l, 7) = x(9 * (i - 1) + l) / self%sq(i); end do
    end do
    !$omp end parallel do
    call apply_h(self, self%wk(:, 7), self%wk(:, 1), self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))
    call apply_rand(self, self%o0, self%o1, self%wk(:, 1), self%wk(:, 2))
    call apply_h(self, self%wk(:, 2), self%wk(:, 3), self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))
    call apply_rand(self, self%e0, self%e1, self%wk(:, 7), self%wk(:, 8))
    if (self%order == 3) then                     ! wk3 <- h o h x - h [ o hw - (p h h x + h p h x)/2 ]
      call apply_h(self, self%wk(:, 1), self%wk(:, 9), self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))      ! h h x
      call apply_rand(self, self%p0, self%p1, self%wk(:, 1), self%wk(:, 10))                              ! p h x
      call apply_h(self, self%wk(:, 10), self%wk(:, 11), self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))    ! h p h x
      call apply_rand(self, self%p0, self%p1, self%wk(:, 9), self%wk(:, 10))                              ! p h h x
      call apply_rand(self, self%o0, self%o1, self%wk(:, 3), self%wk(:, 12))                              ! o hw
      self%wk(:, 12) = self%wk(:, 12) - 0.5_dp * (self%wk(:, 10) + self%wk(:, 11))
      call apply_h(self, self%wk(:, 12), self%wk(:, 9), self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))
      self%wk(:, 3) = self%wk(:, 3) - self%wk(:, 9)
    end if
    !$omp parallel do private(l)
    do i = 1, self%nn
      do l = 1, 9
        y(9 * (i - 1) + l) = self%sq(i) * (self%wk(9 * (i - 1) + l, 8) + self%wk(9 * (i - 1) + l, 1) - self%wk(9 * (i - 1) + l, 3))
      end do
    end do
    !$omp end parallel do
  end subroutine augsym_apply

  !> y = h p h x (normalised coordinates; p = p^gamma): the overlap term of the third-order Hamiltonian, used for the
  !> sphere-charge start vector (1 - h p h / 2)|L>.
  subroutine augsym_hph(self, x, y)
    class(augsym_op), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    integer :: i, l
    if (.not. associated(self%wk)) error stop 'augsym_hph: work space not allocated'
    !$omp parallel do private(l)
    do i = 1, self%nn
      do l = 1, 9; self%wk(9 * (i - 1) + l, 7) = x(9 * (i - 1) + l) / self%sq(i); end do
    end do
    !$omp end parallel do
    call apply_h(self, self%wk(:, 7), self%wk(:, 1), self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))
    call apply_rand(self, self%p0, self%p1, self%wk(:, 1), self%wk(:, 2))
    call apply_h(self, self%wk(:, 2), self%wk(:, 3), self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))
    !$omp parallel do private(l)
    do i = 1, self%nn
      do l = 1, 9; y(9 * (i - 1) + l) = self%sq(i) * self%wk(9 * (i - 1) + l, 3); end do
    end do
    !$omp end parallel do
  end subroutine augsym_hph

  !> Elementary operators on a real vector in normalised coordinates:  mw = M w  (configuration operator
  !> M_R = x + (y-x) P_down + sqrt(xy) T at every node)  and  sw = S w  (structure-constant hop, incl. on-site block).
  subroutine augsym_ops(self, w, mw, sw)
    class(augsym_op), intent(in) :: self
    real(dp), intent(in) :: w(:)
    real(dp), intent(out) :: mw(:), sw(:)
    real(dp), parameter :: z9(9) = 0.0_dp, o9(9) = 1.0_dp
    integer :: i, l
    if (.not. associated(self%wk)) error stop 'augsym_ops: work space not allocated'
    !$omp parallel do private(l)
    do i = 1, self%nn
      do l = 1, 9; self%wk(9 * (i - 1) + l, 7) = w(9 * (i - 1) + l) / self%sq(i); end do
    end do
    !$omp end parallel do
    call apply_rand(self, z9, o9, self%wk(:, 7), self%wk(:, 1))
    call apply_s(self, self%wk(:, 7), self%wk(:, 2))
    !$omp parallel do private(l)
    do i = 1, self%nn
      do l = 1, 9
        mw(9 * (i - 1) + l) = self%sq(i) * self%wk(9 * (i - 1) + l, 1)
        sw(9 * (i - 1) + l) = self%sq(i) * self%wk(9 * (i - 1) + l, 2)
      end do
    end do
    !$omp end parallel do
  end subroutine augsym_ops

  !> t2 = S x  (amplitude coordinates, representatives only)
  subroutine apply_s(self, x, t2)
    class(augsym_op), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: t2(:)
    integer :: i, k, l, l2, g, r
    real(dp) :: acc(9)
    !$omp parallel do private(k, l, l2, g, r, acc)
    do i = 1, self%nn
      acc = 0.0_dp
      do l = 1, 9
        do l2 = 1, 9; acc(l) = acc(l) + self%block0(l, l2) * x(9 * (i - 1) + l2); end do
      end do
      do k = 1, self%nslot
        r = self%hopr(k, i)
        if (r == 0) cycle
        g = self%hopg(k, i)
        do l = 1, 9
          do l2 = 1, 9
            acc(l) = acc(l) + self%chi(g) * self%btab(l, l2, k, g) * x(9 * (r - 1) + l2)
          end do
        end do
      end do
      t2(9 * (i - 1) + 1:9 * i) = acc
    end do
    !$omp end parallel do
  end subroutine apply_s

  !> w = (X0(l) + X1(l) M) x
  subroutine apply_rand(self, x0, x1, x, w)
    class(augsym_op), intent(in) :: self
    real(dp), intent(in) :: x0(:), x1(:), x(:)
    real(dp), intent(out) :: w(:)
    integer :: i, l
    real(dp) :: tf(9), mx
    !$omp parallel do private(l, tf, mx)
    do i = 1, self%nn
      tf = 0.0_dp
      if (self%flipr(i) > 0) call at_link(self, x, self%flipr(i), self%flipg(i), tf)
      do l = 1, 9
        mx = self%md(i) * x(9 * (i - 1) + l) + self%s * tf(l)
        w(9 * (i - 1) + l) = x0(l) * x(9 * (i - 1) + l) + x1(l) * mx
      end do
    end do
    !$omp end parallel do
  end subroutine apply_rand

  !> hx = c x + D S D x
  subroutine apply_h(self, x, hx, t1, t2, t3)
    class(augsym_op), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: hx(:), t1(:), t2(:), t3(:)
    integer :: i, k, l, l2, g, r
    real(dp) :: acc(9), v9(9)
    call apply_rand(self, self%d0, self%d1, x, t1)                     ! t1 = D x   (on representatives)
    !$omp parallel do private(k, l, l2, g, r, acc, v9)
    do i = 1, self%nn                                                  ! t2 = S t1
      acc = 0.0_dp
      do l = 1, 9
        do l2 = 1, 9; acc(l) = acc(l) + self%block0(l, l2) * t1(9 * (i - 1) + l2); end do
      end do
      do k = 1, self%nslot
        r = self%hopr(k, i)
        if (r == 0) cycle
        g = self%hopg(k, i)
        do l = 1, 9
          do l2 = 1, 9
            acc(l) = acc(l) + self%chi(g) * self%btab(l, l2, k, g) * t1(9 * (r - 1) + l2)
          end do
        end do
      end do
      t2(9 * (i - 1) + 1:9 * i) = acc
    end do
    !$omp end parallel do
    call apply_rand(self, self%d0, self%d1, t2, t3)                    ! t3 = D S D x
    call apply_rand(self, self%c0, self%c1, x, hx)                     ! hx = c x
    hx = hx + t3
    v9 = 0.0_dp
  end subroutine apply_h

  ! ------------------------------------------------------------------ canonical nodes
  !> g*(site,cfg) = canonical representative (csite, ccfg); returns that operation g, and the orbit size.
  subroutine canon_node(perm, site, cfg, csite, ccfg, osize, gbest)
    integer, intent(in) :: perm(:, :), site, cfg(:)
    integer, intent(out) :: csite, osize, gbest
    integer, allocatable, intent(out) :: ccfg(:)
    integer :: iop, k, n, s2, ntie
    integer, allocatable :: c2(:)
    logical :: better, equal
    n = size(cfg)
    allocate (c2(n), ccfg(n))
    csite = huge(1); ccfg = huge(1); ntie = 0; gbest = 1
    do iop = 1, nop
      s2 = perm(site, iop)
      do k = 1, n; c2(k) = perm(cfg(k), iop); end do
      call sort_small(c2)
      better = .false.; equal = .false.
      if (s2 < csite) then
        better = .true.
      else if (s2 == csite) then
        equal = .true.
        do k = 1, n
          if (c2(k) < ccfg(k)) then; better = .true.; equal = .false.; exit
          else if (c2(k) > ccfg(k)) then; equal = .false.; exit
          end if
        end do
      end if
      if (better) then
        csite = s2; ccfg = c2; ntie = 1; gbest = iop
      else if (equal) then
        ntie = ntie + 1
      end if
    end do
    osize = nop / ntie
  end subroutine canon_node

  pure subroutine sort_small(a)
    integer, intent(inout) :: a(:)
    integer :: i, j, t
    do i = 2, size(a)
      t = a(i); j = i - 1
      do while (j >= 1)
        if (a(j) <= t) exit
        a(j + 1) = a(j); j = j - 1
      end do
      a(j + 1) = t
    end do
  end subroutine sort_small

  pure subroutine uniq4(a, b, c, d, s, n)
    integer, intent(in) :: a, b, c, d
    integer, allocatable, intent(out) :: s(:)
    integer, intent(out) :: n
    integer :: t(4), v(4), i
    v = [a, b, c, d]; n = 0
    do i = 1, 4
      if (any(t(1:n) == v(i))) cycle
      n = n + 1; t(n) = v(i)
    end do
    allocate (s(n)); s = t(1:n)
  end subroutine uniq4

  pure subroutine uniq3(a, b, c, s, n)
    integer, intent(in) :: a, b, c
    integer, allocatable, intent(out) :: s(:)
    integer, intent(out) :: n
    integer :: t(3)
    n = 1; t(1) = a
    if (b /= a) then; n = n + 1; t(n) = b; end if
    if (c /= a .and. c /= b) then; n = n + 1; t(n) = c; end if
    allocate (s(n)); s = t(1:n)
  end subroutine uniq3

  ! ------------------------------------------------------------------ builder
  !> order 3: three hops per application of H (h o h o h, h p h h), closure over paths R -> R1 -> R2 -> R3.
  subroutine build_augsym_operator(g, mdl, perm, dmat, start, maxlevel, op, maxrank, flippable, order)
    type(graph_t), intent(in) :: g
    type(orb_alloy_t), intent(in) :: mdl
    integer, intent(in) :: perm(:, :), start, maxlevel
    real(dp), intent(in) :: dmat(9, 9, nop)
    type(augsym_op), intent(out) :: op
    integer, intent(in), optional :: maxrank
    logical, intent(in), optional :: flippable(:)
    integer, intent(in), optional :: order
    integer :: ord, k3, r3
    type(store_t) :: st
    logical, allocatable :: canflip(:)
    integer, allocatable :: c0(:), cfg0(:), sites(:), cc(:), tmp(:), ccf(:), c1(:)
    integer :: i, lev, r, r1, r2, k1, k2, nsit, f, j, idx, rkmax, cs, osz, gb, ns, k, l, gg, jj
    logical :: isnew, random_present

    rkmax = huge(1); if (present(maxrank)) rkmax = maxrank
    allocate (canflip(g%n)); canflip = .true.
    if (present(flippable)) canflip = flippable
    random_present = (mdl%x > 0.0_dp .and. mdl%x < 1.0_dp) .and. &
      any([any(mdl%c1 /= 0.0_dp), any(mdl%d1 /= 0.0_dp), any(mdl%o1 /= 0.0_dp), any(mdl%e1 /= 0.0_dp)])
    ord = 2; if (present(order)) ord = order
    if (ord == 3 .and. allocated(mdl%p1)) random_present = random_present .or. &
      ((mdl%x > 0.0_dp .and. mdl%x < 1.0_dp) .and. any(mdl%p1 /= 0.0_dp))
    call store_init(st, 4096)
    allocate (cfg0(0))
    call canon_node(perm, start, cfg0, cs, ccf, osz, gb)
    if (osz /= 1) error stop 'build_augsym_operator: start site is not invariant under the group'
    idx = find_or_add(st, cs, ccf, 0, isnew); st%orb(idx) = 1
    i = 1
    do while (i <= st%n)
      lev = st%level(i)
      if (lev >= maxlevel) then; i = i + 1; cycle; end if
      r = st%site(i)
      c0 = st%data(st%start(i):st%start(i) + st%len(i) - 1)
      do k1 = 0, g%nbr_ptr(r + 1) - g%nbr_ptr(r)
        if (k1 == 0) then; r1 = r; else; r1 = g%nbr(g%nbr_ptr(r) + k1 - 1); end if
        do k2 = 0, g%nbr_ptr(r1 + 1) - g%nbr_ptr(r1)
          if (k2 == 0) then; r2 = r1; else; r2 = g%nbr(g%nbr_ptr(r1) + k2 - 1); end if
          do k3 = 0, merge(g%nbr_ptr(r2 + 1) - g%nbr_ptr(r2), 0, ord == 3)      ! order 3: third hop
          if (k3 == 0) then; r3 = r2; else; r3 = g%nbr(g%nbr_ptr(r2) + k3 - 1); end if
          if (ord == 3) then; call uniq4(r, r1, r2, r3, sites, nsit); else; call uniq3(r, r1, r2, sites, nsit); end if
          do f = 0, merge(2**nsit - 1, 0, random_present)
            cc = c0
            do j = 1, nsit
              if (btest(f, j - 1)) then; tmp = toggled(cc, sites(j)); cc = tmp; end if
            end do
            if (size(cc) > rkmax) cycle
            if (any([(.not. canflip(sites(j)) .and. btest(f, j - 1), j=1, nsit)])) cycle
            call canon_node(perm, r3, cc, cs, ccf, osz, gb)
            idx = find_or_add(st, cs, ccf, lev + 1, isnew)
            if (isnew) st%orb(idx) = osz
          end do
          end do
        end do
      end do
      i = i + 1
    end do

    ! ---- tables
    op%nn = st%n; op%nl = 9; op%x = mdl%x; op%y = 1.0_dp - mdl%x; op%s = sqrt(mdl%x * (1.0_dp - mdl%x))
    ns = maxval(g%nbr_ptr(2:g%n + 1) - g%nbr_ptr(1:g%n)); op%nslot = ns
    allocate (op%flipr(st%n), op%flipg(st%n), op%hopr(ns, st%n), op%hopg(ns, st%n), op%md(st%n), op%level(st%n), op%norb(st%n), op%sq(st%n))
    op%flipr = 0; op%flipg = 1; op%hopr = 0; op%hopg = 1
    op%c0 = mdl%c0; op%c1 = mdl%c1; op%d0 = mdl%d0; op%d1 = mdl%d1; op%o0 = mdl%o0; op%o1 = mdl%o1; op%e0 = mdl%e0; op%e1 = mdl%e1
    allocate (op%p0(9), op%p1(9)); op%p0 = 0.0_dp; op%p1 = 0.0_dp
    if (allocated(mdl%p0)) then; op%p0 = mdl%p0; op%p1 = mdl%p1; end if
    op%order = ord
    op%block = mdl%block
    allocate (op%block0(9, 9)); op%block0 = 0.0_dp
    if (allocated(mdl%block0)) op%block0 = mdl%block0
    ! orbital operation tables: pi(l,g), sg(l,g) with D[pi(l),l] = sg(l)
    do gg = 1, nop
      do l = 1, 9
        op%pi(l, gg) = 0; op%sg(l, gg) = 0.0_dp
        do jj = 1, 9
          if (abs(dmat(jj, l, gg)) > 0.5_dp) then; op%pi(l, gg) = jj; op%sg(l, gg) = dmat(jj, l, gg); end if
        end do
        if (op%pi(l, gg) == 0) error stop 'build_augsym_operator: orbital operation is not monomial'
      end do
    end do
    ! btab(l, l2, slot, g) = sum_{l'} B[l,l'] sg(l',g) delta(l2, pi(l',g))
    allocate (op%btab(9, 9, ns, nop)); op%btab = 0.0_dp
    do gg = 1, nop
      do k = 1, ns
        do l = 1, 9
          do jj = 1, 9
            op%btab(l, op%pi(jj, gg), k, gg) = op%btab(l, op%pi(jj, gg), k, gg) + mdl%block(l, jj, k) * op%sg(jj, gg)
          end do
        end do
      end do
    end do
    op%chi = 1.0_dp
    do i = 1, st%n
      r = st%site(i)
      c0 = st%data(st%start(i):st%start(i) + st%len(i) - 1)
      op%level(i) = st%level(i); op%norb(i) = real(st%orb(i), dp); op%sq(i) = sqrt(op%norb(i))
      op%md(i) = op%x; if (is_down(c0, r)) op%md(i) = op%y
      if (canflip(r)) then
        c1 = toggled(c0, r)
        call canon_node(perm, r, c1, cs, ccf, osz, gb)
        j = find_state(st, cs, ccf)
        if (j > 0) then; op%flipr(i) = j; op%flipg(i) = gb; end if
      end if
      do k = 1, g%nbr_ptr(r + 1) - g%nbr_ptr(r)
        call canon_node(perm, g%nbr(g%nbr_ptr(r) + k - 1), c0, cs, ccf, osz, gb)
        j = find_state(st, cs, ccf)
        if (j > 0) then; op%hopr(k, i) = j; op%hopg(k, i) = gb; end if
      end do
    end do
    allocate (op%wk(9 * st%n, merge(12, 8, ord == 3)))
    call augsym_set_start(op, 1, dmat)
  end subroutine build_augsym_operator
end module asr_augsym
