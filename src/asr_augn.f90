!> n-component augmented space for the TB-LMTO second-order Hamiltonian  H = E_nu + h - h o h,  h = c + d S d,  with a random
!> species s_R in {1..ns} (probabilities x_s) at every site (high-entropy / multicomponent alloys, uncorrelated disorder).
!>
!> Local configuration space of one site = functions on the species set with the inner product <f,g> = sum_s x_s f(s) g(s); an
!> orthonormal basis phi_0 = 1, phi_1 .. phi_{ns-1} (Gram-Schmidt on the species indicators) gives, for the multiplication
!> operator by the indicator of species s, the rank-one matrix  M^(s)_{lk} = x_s phi_l(s) phi_k(s).  A site-diagonal random
!> parameter X_R = sum_s X_s [s_R = s] therefore becomes  W_{lk} = sum_s X_s M^(s)_{lk}  acting on the local state l -> k, and
!> <<f(X)>> = <{}| f(W) |{}> (augmented-space theorem; ns = 2 reproduces asr_augorb with md = x, y and flip amplitude sqrt(xy)).
!> A node is (site R, flipped sites with their states 1..ns-1); a config entry is the integer code  site*ns + state.
!> A site may be declared non-flippable (average medium, W -> X_avg = W_00) or frozen to one species (W -> X_species, used by the
!> energy-dependent route: central site of species c).
module asr_augn
  use asr_kinds, only: dp
  use asr_operator, only: linop
  use asr_augspace, only: graph_t, store_t, store_init, find_or_add, find_state
  implicit none
  private
  public :: multi_alloy_t, augn_op, build_augn_operator, augn_freeze_site, augn_rand, augn_rand_c, augn_basis

  type :: multi_alloy_t
    integer :: nl = 1, ns = 2
    real(dp), allocatable :: conc(:)                     !! x_s
    real(dp), allocatable :: c(:, :), d(:, :), o(:, :), e(:, :)    !! (nl, ns): c = C - E_nu, d = sqrt(Delta), o, E_nu per species
    real(dp), allocatable :: block(:, :, :), block0(:, :)
  end type multi_alloy_t

  type, extends(linop) :: augn_op
    integer :: nn = 0, nl = 1, ns = 2, nslot = 0
    real(dp), allocatable :: conc(:), mt(:, :, :)        !! mt(l+1, k+1, s) = x_s phi_l(s) phi_k(s)
    real(dp), allocatable :: c(:, :), d(:, :), o(:, :), e(:, :), block(:, :, :), block0(:, :)
    integer, allocatable :: flipidx(:, :), hop(:, :), lst(:), site(:), frozen(:), level(:)
    real(dp), pointer :: wk(:, :) => null()
  contains
    procedure :: apply => augn_apply
    procedure :: dim => augn_dim
  end type augn_op

contains

  pure integer function augn_dim(self)
    class(augn_op), intent(in) :: self
    augn_dim = self%nl * self%nn
  end function augn_dim

  !> phi(l, s), l = 0..ns-1 (row l+1), s = 1..ns: orthonormal in <f,g> = sum_s x_s f(s) g(s), phi_0 = 1.
  subroutine augn_basis(conc, phi)
    real(dp), intent(in) :: conc(:)
    real(dp), allocatable, intent(out) :: phi(:, :)
    integer :: ns, l, m
    real(dp) :: ip
    ns = size(conc)
    allocate (phi(ns, ns)); phi = 0.0_dp
    phi(1, :) = 1.0_dp
    do l = 2, ns
      phi(l, :) = 0.0_dp; phi(l, l - 1) = 1.0_dp                         ! indicator of species l-1
      do m = 1, l - 1                                                    ! twice: numerical stability
        ip = sum(conc * phi(m, :) * phi(l, :)); phi(l, :) = phi(l, :) - ip * phi(m, :)
      end do
      do m = 1, l - 1
        ip = sum(conc * phi(m, :) * phi(l, :)); phi(l, :) = phi(l, :) - ip * phi(m, :)
      end do
      phi(l, :) = phi(l, :) / sqrt(sum(conc * phi(l, :)**2))
    end do
  end subroutine augn_basis

  ! ------------------------------------------------------------------ configurations (sorted integer codes site*ns + state)
  pure integer function get_state(cfg, r, ns)
    integer, intent(in) :: cfg(:), r, ns
    integer :: i
    get_state = 0
    do i = 1, size(cfg)
      if (cfg(i) / ns == r) then; get_state = mod(cfg(i), ns); return; end if
    end do
  end function get_state

  pure function set_state(cfg, r, k, ns) result(out)
    integer, intent(in) :: cfg(:), r, k, ns
    integer, allocatable :: out(:)
    integer, allocatable :: t(:)
    integer :: i, n, code
    allocate (t(size(cfg) + 1)); n = 0
    do i = 1, size(cfg)
      if (cfg(i) / ns /= r) then; n = n + 1; t(n) = cfg(i); end if
    end do
    if (k > 0) then
      code = r * ns + k; n = n + 1; t(n) = code
      i = n
      do while (i > 1)                                                   ! insertion sort of the appended entry
        if (t(i - 1) <= t(i)) exit
        code = t(i); t(i) = t(i - 1); t(i - 1) = code; i = i - 1
      end do
    end if
    allocate (out(n)); out = t(1:n)
  end function set_state

  ! ------------------------------------------------------------------ application
  !> y = H x  with  x(L + nl*(node-1)),  H = E_nu + h - h o h.
  subroutine augn_apply(self, x, y)
    class(augn_op), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    if (.not. associated(self%wk)) error stop 'augn_apply: work space not allocated'
    call apply_h(self, x, self%wk(:, 1), self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))
    call apply_rand(self, self%o, self%wk(:, 1), self%wk(:, 2))
    call apply_h(self, self%wk(:, 2), self%wk(:, 3), self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))
    call apply_rand(self, self%e, x, y)
    y = y + self%wk(:, 1) - self%wk(:, 3)
  end subroutine augn_apply

  !> w = X x, X = sum_s X_s(l) [s_R = s] in the augmented space (par(l, s)); a non-flippable site carries the average
  !> X_avg = W_00, a frozen site X_species.
  subroutine apply_rand(self, par, x, w)
    class(augn_op), intent(in) :: self
    real(dp), intent(in) :: par(:, :), x(:)
    real(dp), intent(out) :: w(:)
    real(dp) :: wt(self%nl, self%ns, self%ns)
    integer :: i, lam, l, k, j, s
    real(dp) :: acc
    wt = 0.0_dp
    do lam = 1, self%nl
      do l = 1, self%ns; do k = 1, self%ns
          do s = 1, self%ns; wt(lam, l, k) = wt(lam, l, k) + par(lam, s) * self%mt(l, k, s); end do
      end do; end do
    end do
    !$omp parallel do private(lam, l, k, j, acc)
    do i = 1, self%nn
      do lam = 1, self%nl
        if (self%frozen(i) > 0) then
          w(self%nl * (i - 1) + lam) = par(lam, self%frozen(i)) * x(self%nl * (i - 1) + lam)
          cycle
        end if
        l = self%lst(i) + 1
        acc = 0.0_dp
        do k = 1, self%ns
          j = self%flipidx(k, i)
          if (j == 0) cycle
          acc = acc + wt(lam, l, k) * x(self%nl * (j - 1) + lam)
        end do
        w(self%nl * (i - 1) + lam) = acc
      end do
    end do
    !$omp end parallel do
  end subroutine apply_rand

  !> w = X x with real parameters par(l, s) (public wrapper)
  subroutine augn_rand(self, par, x, w)
    class(augn_op), intent(in) :: self
    real(dp), intent(in) :: par(:, :), x(:)
    real(dp), intent(out) :: w(:)
    call apply_rand(self, par, x, w)
  end subroutine augn_rand

  !> complex parameters and vectors (energy-dependent diagonal e_s(E)): w = X x
  subroutine augn_rand_c(self, par, x, w)
    class(augn_op), intent(in) :: self
    complex(dp), intent(in) :: par(:, :), x(:)
    complex(dp), intent(out) :: w(:)
    complex(dp) :: wt(self%nl, self%ns, self%ns), acc
    integer :: i, lam, l, k, j, s
    wt = (0.0_dp, 0.0_dp)
    do lam = 1, self%nl
      do l = 1, self%ns; do k = 1, self%ns
          do s = 1, self%ns; wt(lam, l, k) = wt(lam, l, k) + par(lam, s) * self%mt(l, k, s); end do
      end do; end do
    end do
    !$omp parallel do private(lam, l, k, j, acc)
    do i = 1, self%nn
      do lam = 1, self%nl
        if (self%frozen(i) > 0) then
          w(self%nl * (i - 1) + lam) = par(lam, self%frozen(i)) * x(self%nl * (i - 1) + lam)
          cycle
        end if
        l = self%lst(i) + 1
        acc = (0.0_dp, 0.0_dp)
        do k = 1, self%ns
          j = self%flipidx(k, i)
          if (j == 0) cycle
          acc = acc + wt(lam, l, k) * x(self%nl * (j - 1) + lam)
        end do
        w(self%nl * (i - 1) + lam) = acc
      end do
    end do
    !$omp end parallel do
  end subroutine augn_rand_c

  !> hx = c x + D S D x
  subroutine apply_h(self, x, hx, t1, t2, t3)
    class(augn_op), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: hx(:), t1(:), t2(:), t3(:)
    integer :: i, k, j, l, l2, nl
    real(dp) :: acc
    nl = self%nl
    call apply_rand(self, self%d, x, t1)                    ! t1 = D x
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
        do l2 = 1, nl
          acc = acc + self%block0(l, l2) * t1(nl * (i - 1) + l2)
        end do
        t2(nl * (i - 1) + l) = acc
      end do
    end do
    !$omp end parallel do
    call apply_rand(self, self%d, t2, t3)                   ! t3 = D S D x
    call apply_rand(self, self%c, x, hx)                    ! hx = c x
    hx = hx + t3
  end subroutine apply_h

  !> Freeze `site` to species `sp`: W -> X_sp on every node of that site.  The site must not have been flippable when the node
  !> set was built (flippable(site) = .false.), so that no node carries a flip of it.
  subroutine augn_freeze_site(op, site, sp)
    type(augn_op), intent(inout) :: op
    integer, intent(in) :: site, sp
    integer :: i
    do i = 1, op%nn
      if (op%site(i) == site) then
        if (op%lst(i) /= 0) error stop 'augn_freeze_site: the frozen site was flippable'
        op%frozen(i) = sp
      end if
    end do
  end subroutine augn_freeze_site

  ! ------------------------------------------------------------------ builder
  !> Node set = closure of the support of H (paths R -> R1 -> R2, every visited site in any of its ns states), rank-limited.
  subroutine build_augn_operator(g, mdl, start, maxlevel, op, maxrank, flippable)
    type(graph_t), intent(in) :: g
    type(multi_alloy_t), intent(in) :: mdl
    integer, intent(in) :: start, maxlevel
    type(augn_op), intent(out) :: op
    integer, intent(in), optional :: maxrank
    logical, intent(in), optional :: flippable(:)
    type(store_t) :: st
    logical, allocatable :: canflip(:)
    integer, allocatable :: c0(:), cfg0(:), sites(:), cc(:), tmp(:)
    real(dp), allocatable :: phi(:, :)
    integer :: ns, i, lev, r, r1, r2, k1, k2, nsit, f, j, rk, idx, k, ncomb, rkmax, dgt, ll
    logical :: isnew, random_present, skip

    ns = mdl%ns
    rkmax = huge(1); if (present(maxrank)) rkmax = maxrank
    allocate (canflip(g%n)); canflip = .true.
    if (present(flippable)) canflip = flippable
    random_present = .false.
    do k = 2, ns
      if (any(mdl%c(:, k) /= mdl%c(:, 1)) .or. any(mdl%d(:, k) /= mdl%d(:, 1)) .or. any(mdl%o(:, k) /= mdl%o(:, 1)) .or. &
          any(mdl%e(:, k) /= mdl%e(:, 1))) random_present = .true.
    end do
    call store_init(st, 4096)
    allocate (cfg0(0))
    idx = find_or_add(st, start, cfg0, 0, isnew)
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
          call uniq3(r, r1, r2, sites, nsit)
          ncomb = 1; if (random_present) ncomb = ns**nsit
          do f = 0, ncomb - 1
            cc = c0; skip = .false.
            do j = 1, nsit
              dgt = mod(f / ns**(j - 1), ns)
              if (dgt > 0 .and. .not. canflip(sites(j))) then; skip = .true.; exit; end if
              if (random_present) then; tmp = set_state(cc, sites(j), dgt, ns); cc = tmp; end if
            end do
            if (skip) cycle
            if (size(cc) > rkmax) cycle
            idx = find_or_add(st, r2, cc, lev + 1, isnew)
          end do
        end do
      end do
      i = i + 1
    end do

    ! ---- tables
    op%nn = st%n; op%nl = mdl%nl; op%ns = ns
    op%conc = mdl%conc
    call augn_basis(mdl%conc, phi)
    allocate (op%mt(ns, ns, ns))
    do k = 1, ns
      do j = 1, ns
        do ll = 1, ns; op%mt(k, j, ll) = mdl%conc(ll) * phi(k, ll) * phi(j, ll); end do
      end do
    end do
    op%nslot = maxval(g%nbr_ptr(2:g%n + 1) - g%nbr_ptr(1:g%n))
    allocate (op%flipidx(ns, st%n), op%hop(op%nslot, st%n), op%lst(st%n), op%site(st%n), op%frozen(st%n), op%level(st%n))
    op%flipidx = 0; op%hop = 0; op%frozen = 0
    op%c = mdl%c; op%d = mdl%d; op%o = mdl%o; op%e = mdl%e
    op%block = mdl%block
    allocate (op%block0(mdl%nl, mdl%nl)); op%block0 = 0.0_dp
    if (allocated(mdl%block0)) op%block0 = mdl%block0
    allocate (op%wk(mdl%nl * st%n, 6))
    do i = 1, st%n
      r = st%site(i)
      c0 = st%data(st%start(i):st%start(i) + st%len(i) - 1)
      op%level(i) = st%level(i); op%site(i) = r
      op%lst(i) = get_state(c0, r, ns)
      do k = 0, ns - 1
        if (k == op%lst(i)) then
          op%flipidx(k + 1, i) = i
        else if (k == 0 .or. canflip(r)) then
          tmp = set_state(c0, r, k, ns)
          op%flipidx(k + 1, i) = find_state(st, r, tmp)
        end if
      end do
      do k = 1, g%nbr_ptr(r + 1) - g%nbr_ptr(r)
        op%hop(k, i) = find_state(st, g%nbr(g%nbr_ptr(r) + k - 1), c0)
      end do
    end do
  end subroutine build_augn_operator

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
end module asr_augn
