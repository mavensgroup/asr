!> Augmented-space Hamiltonian for multi-orbital alloys of the TB-LMTO second-order form
!>     H = E_nu + h - h o h ,   h = c + d S d      (c = C - E_nu, d = sqrt(Delta), o)
!> with binary-random per-site, per-orbital parameters X_R(l) = X0(l) + X1(l) n_R, n_R = 1 (A) or 0 (B).
!> In the augmented space  n_R -> M_R = x + (y-x) P_down + sqrt(xy) T  (a projector, M^2 = M), so each
!> random diagonal operator is  X0 + X1 M_R.  H is applied matrix-free as a sequence of the two elementary
!> operators  D (random diagonal, creates flips at the current site)  and  S (hops; bond-block table).
!>
!> Node set: pairs (R, C) (C = set of flipped sites).  It is generated as the closure of the support of H:
!> from (R,C) every path R -> R1 -> R2 (each step a stay or a hop) with any flip subset F of {R,R1,R2},
!> giving (R2, C xor F).  Every intermediate vector of the matrix-free application lies in this set, so
!> the application is exact for source vectors supported on levels < maxlevel.  States with more than
!> `maxrank` flipped sites are dropped (the sites beyond are treated as average medium).
module asr_augorb
  use asr_kinds, only: dp
  use asr_operator, only: linop
  use asr_augspace, only: graph_t, store_t, store_init, find_or_add, find_state, toggled, is_down
  implicit none
  private
  public :: orb_alloy_t, augorb_op, build_augorb_operator, augorb_set_params, augorb_hph, augorb_freeze_site

  type :: orb_alloy_t
    integer :: nl = 1
    real(dp) :: x = 0.5_dp
    real(dp), allocatable :: c0(:), c1(:), d0(:), d1(:), o0(:), o1(:), e0(:), e1(:)   !! (nl); E_nu = e0 + e1 n
    real(dp), allocatable :: p0(:), p1(:)          !! p^gamma (third-order Hamiltonian only)
    real(dp), allocatable :: block(:, :, :)        !! block(L,L',slot): S_{R,R'} for R' = nbr(R,slot)
    real(dp), allocatable :: block0(:, :)          !! on-site block S_{R,R} (optional; zero if not allocated)
  end type orb_alloy_t

  type, extends(linop) :: augorb_op
    integer :: nn = 0                                  !! nodes
    integer :: nl = 1, nslot = 0
    real(dp) :: x = 0.5_dp, s = 0.5_dp, y = 0.5_dp
    real(dp), allocatable :: c0(:), c1(:), d0(:), d1(:), o0(:), o1(:), e0(:), e1(:), block(:, :, :), block0(:, :)
    real(dp), allocatable :: p0(:), p1(:)               !! p^gamma (order 3)
    integer :: order = 2                                !! 2: E + h - hoh;  3: + hohoh - (hphh + hhph)/2
    integer, allocatable :: flip(:), hop(:, :), level(:), site(:)
    real(dp), allocatable :: md(:)                     !! diagonal part of M_R on each node (x or y)
    real(dp), pointer :: wk(:, :) => null()            !! work space (n, 6), allocated once (apply is not thread-safe across calls)
  contains
    procedure :: apply => augorb_apply
    procedure :: dim => augorb_dim
    procedure :: mapply => augorb_mapply
  end type augorb_op

contains

  pure integer function augorb_dim(self)
    class(augorb_op), intent(in) :: self
    augorb_dim = self%nl * self%nn
  end function augorb_dim

  !> y = H x  with  x(L + nl*(node-1)).
  subroutine augorb_apply(self, x, y)
    class(augorb_op), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    integer :: n
    n = self%nl * self%nn
    if (.not. associated(self%wk)) error stop 'augorb_apply: work space not allocated'
    call apply_h(self, x, self%wk(:, 1), self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))     ! hx
    call apply_o(self, self%wk(:, 1), self%wk(:, 2))                                        ! w = o hx
    call apply_h(self, self%wk(:, 2), self%wk(:, 3), self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))   ! hw
    call apply_enu(self, x, y)
    y = y + self%wk(:, 1) - self%wk(:, 3)
    if (self%order == 3) then                     ! + h [ o hw - (p h h x + h p h x)/2 ]
      call apply_h(self, self%wk(:, 1), self%wk(:, 7), self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))     ! h h x
      call apply_rand(self, self%p0, self%p1, self%wk(:, 1), self%wk(:, 8))                              ! p h x
      call apply_h(self, self%wk(:, 8), self%wk(:, 9), self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))     ! h p h x
      call apply_rand(self, self%p0, self%p1, self%wk(:, 7), self%wk(:, 8))                              ! p h h x
      call apply_o(self, self%wk(:, 3), self%wk(:, 10))                                                  ! o hw
      self%wk(:, 10) = self%wk(:, 10) - 0.5_dp * (self%wk(:, 8) + self%wk(:, 9))
      call apply_h(self, self%wk(:, 10), self%wk(:, 7), self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))
      y = y + self%wk(:, 7)
    end if
  end subroutine augorb_apply

  !> y = M x  (the configuration operator n_R -> M_R at every node, no parameters)
  subroutine augorb_mapply(self, x, y)
    class(augorb_op), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    real(dp) :: z9(self%nl), o9(self%nl)
    z9 = 0.0_dp; o9 = 1.0_dp
    call apply_rand(self, z9, o9, x, y)
  end subroutine augorb_mapply

  !> Freeze `site` to a definite species: M_R -> mc (1 = A, 0 = B) on every node of that site.  The site must not have been
  !> flippable when the node set was built (flippable(site) = .false.), so that no node carries a flip of it.
  subroutine augorb_freeze_site(op, site, mc)
    type(augorb_op), intent(inout) :: op
    integer, intent(in) :: site
    real(dp), intent(in) :: mc
    integer :: i
    do i = 1, op%nn
      if (op%site(i) == site) then
        if (op%flip(i) /= 0) error stop 'augorb_freeze_site: the frozen site was flippable'
        op%md(i) = mc
      end if
    end do
  end subroutine augorb_freeze_site

  !> y = h p h x (p = p^gamma): overlap term, for the sphere-charge start vector (1 - h p h / 2)|L>
  subroutine augorb_hph(self, x, y)
    class(augorb_op), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    if (.not. associated(self%wk)) error stop 'augorb_hph: work space not allocated'
    call apply_h(self, x, self%wk(:, 1), self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))
    call apply_rand(self, self%p0, self%p1, self%wk(:, 1), self%wk(:, 2))
    call apply_h(self, self%wk(:, 2), y, self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))
  end subroutine augorb_hph

  subroutine apply_enu(self, x, y)
    class(augorb_op), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    call apply_rand(self, self%e0, self%e1, x, y)
  end subroutine apply_enu

  !> w = (o0 + o1 M) x
  subroutine apply_o(self, x, w)
    class(augorb_op), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: w(:)
    call apply_rand(self, self%o0, self%o1, x, w)
  end subroutine apply_o

  !> w = (X0(l) + X1(l) M) x   (M acts on the configuration part of the node)
  subroutine apply_rand(self, x0, x1, x, w)
    class(augorb_op), intent(in) :: self
    real(dp), intent(in) :: x0(:), x1(:), x(:)
    real(dp), intent(out) :: w(:)
    integer :: i, l, j
    real(dp) :: mx
    !$omp parallel do private(l, j, mx)
    do i = 1, self%nn
      j = self%flip(i)
      do l = 1, self%nl
        mx = self%md(i) * x(self%nl * (i - 1) + l)
        if (j > 0) mx = mx + self%s * x(self%nl * (j - 1) + l)
        w(self%nl * (i - 1) + l) = x0(l) * x(self%nl * (i - 1) + l) + x1(l) * mx
      end do
    end do
    !$omp end parallel do
  end subroutine apply_rand

  !> hx = c x + D S D x
  subroutine apply_h(self, x, hx, t1, t2, t3)
    class(augorb_op), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: hx(:), t1(:), t2(:), t3(:)
    integer :: i, k, j, l, l2, nl
    real(dp) :: acc
    nl = self%nl
    call apply_rand(self, self%d0, self%d1, x, t1)          ! t1 = D x
    !$omp parallel do private(k, j, l, l2, acc)
    do i = 1, self%nn                                       ! t2 = S t1
      do l = 1, nl
        acc = 0.0_dp
        do k = 1, self%nslot
          j = self%hop(k, i)
          if (j == 0) cycle
          do l2 = 1, nl
            acc = acc + self%block(l, l2, k) * t1(nl * (j - 1) + l2)
          end do
        end do
        do l2 = 1, nl                                   ! on-site block S_RR
          acc = acc + self%block0(l, l2) * t1(nl * (i - 1) + l2)
        end do
        t2(nl * (i - 1) + l) = acc
      end do
    end do
    !$omp end parallel do
    call apply_rand(self, self%d0, self%d1, t2, t3)         ! t3 = D S D x
    call apply_rand(self, self%c0, self%c1, x, hx)          ! hx = c x
    hx = hx + t3
  end subroutine apply_h

  !> Build the node set (closure of the support of H) to `maxlevel` and the elementary tables.
  !> order 3: H contains three hops per application (h o h o h, h p h h): the closure follows paths R -> R1 -> R2 -> R3.
  subroutine build_augorb_operator(g, mdl, start, maxlevel, op, maxrank, flippable, order)
    type(graph_t), intent(in) :: g
    type(orb_alloy_t), intent(in) :: mdl
    integer, intent(in) :: start, maxlevel
    type(augorb_op), intent(out) :: op
    integer, intent(in), optional :: maxrank
    logical, intent(in), optional :: flippable(:)     !! sites that may be flipped (others: average medium, never flipped)
    integer, intent(in), optional :: order
    integer :: ord, k3, r3
    logical, allocatable :: canflip(:)
    type(store_t) :: st
    integer, allocatable :: c0(:), cfg0(:), sites(:), cc(:), tmp(:)
    integer :: i, lev, r, r1, r2, k1, k2, nsit, f, nf, rk, idx, j, k, ns
    logical :: isnew
    integer :: rkmax
    logical :: random_present

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
    idx = find_or_add(st, start, cfg0, 0, isnew)
    i = 1
    do while (i <= st%n)
      lev = st%level(i)
      if (lev >= maxlevel) then; i = i + 1; cycle; end if
      r = st%site(i)
      c0 = st%data(st%start(i):st%start(i) + st%len(i) - 1)
      do k1 = 0, g%nbr_ptr(r + 1) - g%nbr_ptr(r)                 ! k1 = 0: stay, else hop slot
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
            idx = find_or_add(st, r3, cc, lev + 1, isnew)
          end do
          end do
        end do
      end do
      i = i + 1
    end do

    ! ---- elementary tables
    op%nn = st%n; op%nl = mdl%nl; op%x = mdl%x; op%y = 1.0_dp - mdl%x; op%s = sqrt(mdl%x * (1.0_dp - mdl%x))
    ns = maxval(g%nbr_ptr(2:g%n + 1) - g%nbr_ptr(1:g%n))
    op%nslot = ns
    allocate (op%flip(st%n), op%hop(ns, st%n), op%md(st%n), op%level(st%n), op%site(st%n))
    op%flip = 0; op%hop = 0
    op%c0 = mdl%c0; op%c1 = mdl%c1; op%d0 = mdl%d0; op%d1 = mdl%d1; op%o0 = mdl%o0; op%o1 = mdl%o1; op%e0 = mdl%e0; op%e1 = mdl%e1
    allocate (op%p0(mdl%nl), op%p1(mdl%nl)); op%p0 = 0.0_dp; op%p1 = 0.0_dp
    if (allocated(mdl%p0)) then; op%p0 = mdl%p0; op%p1 = mdl%p1; end if
    op%order = ord
    op%block = mdl%block
    allocate (op%wk(mdl%nl * st%n, merge(10, 6, ord == 3)))
    allocate (op%block0(mdl%nl, mdl%nl)); op%block0 = 0.0_dp
    if (allocated(mdl%block0)) op%block0 = mdl%block0
    do i = 1, st%n
      r = st%site(i)
      c0 = st%data(st%start(i):st%start(i) + st%len(i) - 1)
      op%level(i) = st%level(i); op%site(i) = r
      op%md(i) = op%x; if (is_down(c0, r)) op%md(i) = op%y
      if (canflip(r)) then
        tmp = toggled(c0, r)
        op%flip(i) = find_state(st, r, tmp)
      end if
      do k = 1, g%nbr_ptr(r + 1) - g%nbr_ptr(r)
        op%hop(k, i) = find_state(st, g%nbr(g%nbr_ptr(r) + k - 1), c0)
      end do
    end do
  end subroutine build_augorb_operator

  !> Refresh the parameter arrays of an existing operator (node structure unchanged).
  subroutine augorb_set_params(op, mdl)
    type(augorb_op), intent(inout) :: op
    type(orb_alloy_t), intent(in) :: mdl
    op%c0 = mdl%c0; op%c1 = mdl%c1; op%d0 = mdl%d0; op%d1 = mdl%d1; op%o0 = mdl%o0; op%o1 = mdl%o1
    op%e0 = mdl%e0; op%e1 = mdl%e1
    if (allocated(mdl%p0)) then; op%p0 = mdl%p0; op%p1 = mdl%p1; end if
  end subroutine augorb_set_params

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
end module asr_augorb
