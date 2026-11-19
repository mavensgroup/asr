!> Site-class (sublattice) resolved, point-group reduced augmented space for  H = E_nu + h - h o h,  h = c + d S d  (order 2),
!> in which every site class has its own species set and concentrations: ordered sublattices, partial (sublattice) disorder and
!> multicomponent (high-entropy) alloys in one operator.
!>
!> Site R of class kappa carries an independent random species s_R in {spc(1,kappa) .. spc(nsc(kappa),kappa)} (global species
!> indices) with probabilities conc(k, kappa).  Its local configuration space has the orthonormal basis phi_0 = 1, phi_1 ..
!> phi_{nsc-1} (augn_basis on the class concentrations), and a site-diagonal random parameter X_R becomes
!>     W_lk = sum_k' X_{spc(k',kappa)} x_k' phi_l(k') phi_k(k')        (local state l -> k);
!> an ordered class (nsc = 1) never flips and carries X of its species.  A node is (site, flipped sites with their states
!> 1..nsc-1); a config entry is the integer code  site*nsmax + state.  A site may be non-flippable (average medium, W -> W_00) or
!> frozen to one species (energy-dependent route: central site of species c).
!>
!> Group reduction as in asr_augsym: canonical representatives, psi(g^-1 r) = chi(g) D(g)^T psi(r), normalised coordinates
!> w = sqrt(n_node) psi, links store (representative, operation).  The group is supplied by the caller: every operation must map
!> each site onto a site of the same class and fix the start site (asr_symmetry::class_ops), operation 1 must be the identity,
!> and the start orbital must have chi(g) = D(g)[L0,L0] = +-1.  A group of the identity alone gives the unreduced operator.
!> The bond blocks are those of a Bravais graph (one block per neighbour slot, the same for every site): one common S for all
!> species, as in asr_lmto_alloy.
module asr_augcls
  use asr_kinds, only: dp
  use asr_operator, only: linop
  use asr_augspace, only: graph_t, store_t, store_init, find_or_add, find_state
  !$ use omp_lib, only: omp_get_max_threads
  use asr_chain, only: chain_vec_par
  implicit none
  private
  public :: augn_basis, class_alloy_t, augcls_op, build_augcls_operator, augcls_set_start, augcls_set_params, augcls_freeze_site, &
            augcls_rand, augcls_rand_c, augcls_set_class_conc, augcls_hp_c, augcls_project_c, augcls_leak_c, augcls_amp_c, augcls_free

  type :: class_alloy_t
    integer :: nl = 9, nsp = 0, ncls = 0
    integer, allocatable :: nsc(:), spc(:, :)                    !! species per class; global species index spc(k, class)
    real(dp), allocatable :: conc(:, :)                           !! conc(k, class), sum_k = 1
    real(dp), allocatable :: c(:, :), d(:, :), o(:, :), e(:, :)   !! (nl, nsp): c = C - E_nu, d = sqrt(Delta), o, E_nu per species
    real(dp), allocatable :: block(:, :, :), block0(:, :)
  end type class_alloy_t

  type, extends(linop) :: augcls_op
    integer :: nn = 0, nl = 9, nslot = 0, nop = 1, nsmax = 1, ncls = 0, nsp = 0, start_l = 1
    integer, allocatable :: nsc(:), spc(:, :)
    real(dp), allocatable :: mt(:, :, :, :)                       !! mt(l+1, k+1, k', class) = x_k' phi_l(k') phi_k(k')
    real(dp), allocatable :: c(:, :), d(:, :), o(:, :), e(:, :), block0(:, :)
    real(dp), allocatable :: btab(:, :, :, :)                     !! btab(l, l2, slot, g): S block times D(g)^T
    !> The node tables below never change after build_augcls_operator.  They are pointers, so that a copy of the operator (one
    !> per thread in the chain loops of asr_class_scf) shares them instead of duplicating them; only the small parameter arrays,
    !> frozen and the work space wk are per copy.  augcls_free releases them (on the original, never on a copy).
    integer, pointer, contiguous :: flipr(:, :) => null(), flipg(:, :) => null()   !! (nsmax, node): node with the site's own state set to k-1
    integer, pointer, contiguous :: hopr(:, :) => null(), hopg(:, :) => null(), lst(:) => null(), site(:) => null(), &
                                    ncl(:) => null(), level(:) => null()
    integer, allocatable :: frozen(:)
    real(dp), pointer, contiguous :: norb(:) => null(), sq(:) => null()
    integer, allocatable :: pi(:, :)
    real(dp), allocatable :: sg(:, :), chi(:)
    real(dp), allocatable :: csg(:, :)                            !! csg(l, g) = chi(g) sg(l, g), exactly 0 or +-1
    integer, pointer, contiguous :: stab_ptr(:) => null(), stab(:) => null()   !! stabiliser of node i: stab(stab_ptr(i):stab_ptr(i+1)-1)
    real(dp), pointer :: wk(:, :) => null()
  contains
    procedure :: apply => augcls_apply
    procedure :: dim => augcls_dim
  end type augcls_op

contains

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

  pure integer function augcls_dim(self)
    class(augcls_op), intent(in) :: self
    augcls_dim = self%nl * self%nn
  end function augcls_dim

  !> start orbital l0: chi(g) = D(g)[l0,l0], which must be +-1 for every operation of the group
  subroutine augcls_set_start(op, l0, dmat)
    type(augcls_op), intent(inout) :: op
    integer, intent(in) :: l0
    real(dp), intent(in) :: dmat(:, :, :)
    integer :: g
    op%start_l = l0
    do g = 1, op%nop
      if (abs(abs(dmat(l0, l0, g)) - 1.0_dp) > 1.0e-9_dp) error stop 'augcls_set_start: start orbital is not one-dimensional under the group'
      op%chi(g) = sign(1.0_dp, dmat(l0, l0, g))
      op%csg(:, g) = op%chi(g) * op%sg(:, g)
    end do
  end subroutine augcls_set_start

  !> new species parameters (same classes, species and node set)
  subroutine augcls_set_params(op, mdl)
    type(augcls_op), intent(inout) :: op
    type(class_alloy_t), intent(in) :: mdl
    if (mdl%nsp /= op%nsp) error stop 'augcls_set_params: number of species changed'
    op%c = mdl%c; op%d = mdl%d; op%o = mdl%o; op%e = mdl%e
  end subroutine augcls_set_params

  !> Freeze `site` to the global species `sp` (W -> X_sp).  The site must have been non-flippable when the node set was built.
  subroutine augcls_freeze_site(op, site, sp)
    type(augcls_op), intent(inout) :: op
    integer, intent(in) :: site, sp
    integer :: i
    if (.not. any(op%spc(1:op%nsc(op%ncl(1)), op%ncl(1)) == sp) .and. op%site(1) == site) &
      error stop 'augcls_freeze_site: species not in the class of the frozen site'
    do i = 1, op%nn
      if (op%site(i) == site) then
        if (op%lst(i) /= 0 .or. count(op%flipr(:, i) /= 0) > 1) error stop 'augcls_freeze_site: the frozen site was flippable'
        op%frozen(i) = sp
      end if
    end do
  end subroutine augcls_freeze_site

  !> New concentrations of class kap (same species list and node set): short-range order sets the first-shell class of the frozen
  !> start site to the conditional concentrations of the central species.  Species of zero weight are allowed (see class_mt).
  subroutine augcls_set_class_conc(op, kap, conc)
    type(augcls_op), intent(inout) :: op
    integer, intent(in) :: kap
    real(dp), intent(in) :: conc(:)
    if (size(conc) /= op%nsc(kap)) error stop 'augcls_set_class_conc: wrong number of concentrations'
    if (any(conc < 0.0_dp) .or. abs(sum(conc) - 1.0_dp) > 1.0e-8_dp) error stop 'augcls_set_class_conc: concentrations'
    call class_mt(conc, op%mt(:, :, :, kap))
  end subroutine augcls_set_class_conc

  !> mt(l+1, k+1, k') = x_k' phi_l(k') phi_k(k') for one class.  The basis phi_0 .. phi_{n'-1} is built on the n' species of
  !> non-zero weight only (augn_basis is singular otherwise); local states n' .. n-1 then have all mt entries zero, so
  !> they are never coupled to the rest of the node set and drop out of every chain.
  subroutine class_mt(conc, mt)
    real(dp), intent(in) :: conc(:)
    real(dp), intent(out) :: mt(:, :, :)
    real(dp), allocatable :: phi(:, :)
    integer, allocatable :: sup(:)
    integer :: l, k, j, n
    mt = 0.0_dp
    sup = pack([(j, j=1, size(conc))], conc > 0.0_dp)
    n = size(sup)
    call augn_basis(conc(sup), phi)
    do l = 1, n
      do k = 1, n
        do j = 1, n; mt(l, k, sup(j)) = conc(sup(j)) * phi(l, j) * phi(k, j); end do
      end do
    end do
  end subroutine class_mt

  ! ------------------------------------------------------------------ configurations (sorted codes site*nsmax + state)
  pure integer function get_state(cfg, r, nsm)
    integer, intent(in) :: cfg(:), r, nsm
    integer :: i
    get_state = 0
    do i = 1, size(cfg)
      if (cfg(i) / nsm == r) then; get_state = mod(cfg(i), nsm); return; end if
    end do
  end function get_state

  pure function set_state(cfg, r, k, nsm) result(out)
    integer, intent(in) :: cfg(:), r, k, nsm
    integer, allocatable :: out(:)
    integer, allocatable :: t(:)
    integer :: i, n, code
    allocate (t(size(cfg) + 1)); n = 0
    do i = 1, size(cfg)
      if (cfg(i) / nsm /= r) then; n = n + 1; t(n) = cfg(i); end if
    end do
    if (k > 0) then
      code = r * nsm + k; n = n + 1; t(n) = code
      i = n
      do while (i > 1)
        if (t(i - 1) <= t(i)) exit
        code = t(i); t(i) = t(i - 1); t(i - 1) = code; i = i - 1
      end do
    end if
    allocate (out(n)); out = t(1:n)
  end function set_state

  ! ------------------------------------------------------------------ application
  !> psi at the node reached by a link (representative r, operation g), from the amplitude vector v
  pure subroutine at_link(op, v, r, g, out)
    class(augcls_op), intent(in) :: op
    real(dp), intent(in) :: v(:)
    integer, intent(in) :: r, g
    real(dp), intent(out) :: out(9)
    integer :: l
    do l = 1, 9
      out(l) = op%csg(l, g) * v(9 * (r - 1) + op%pi(l, g))
    end do
  end subroutine at_link

  !> Amplitudes (9 orbitals) of the unreduced vector at a node, from the reduced vector v in normalised coordinates.  The node is
  !> given by its link (r, g): r its representative, g the operation that maps it onto r (build_augcls_operator, tnode / tgop).
  !> v holds sqrt(orbit size) times the amplitude at r, so psi(node) = chi(g) D(g)^T v(r) / sqrt(orbit size).  r = 0: zero.
  pure function augcls_amp_c(op, v, r, g) result(amp)
    class(augcls_op), intent(in) :: op
    complex(dp), intent(in) :: v(:)
    integer, intent(in) :: r, g
    complex(dp) :: amp(9)
    if (r <= 0) then; amp = (0.0_dp, 0.0_dp); return; end if
    call at_link_c(op, v, r, g, amp)
    amp = amp / op%sq(r)
  end function augcls_amp_c

  pure subroutine at_link_c(op, v, r, g, out)
    class(augcls_op), intent(in) :: op
    complex(dp), intent(in) :: v(:)
    integer, intent(in) :: r, g
    complex(dp), intent(out) :: out(9)
    integer :: l
    do l = 1, 9
      out(l) = op%csg(l, g) * v(9 * (r - 1) + op%pi(l, g))
    end do
  end subroutine at_link_c

  !> y = H x in normalised coordinates,  H = E_nu + h - h o h
  subroutine augcls_apply(self, x, y)
    class(augcls_op), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    integer :: i, l
    if (.not. associated(self%wk)) error stop 'augcls_apply: work space not allocated'
    !$omp parallel do private(l)
    do i = 1, self%nn
      do l = 1, 9; self%wk(9 * (i - 1) + l, 7) = x(9 * (i - 1) + l) / self%sq(i); end do
    end do
    !$omp end parallel do
    call apply_h(self, self%wk(:, 7), self%wk(:, 1), self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))
    call apply_rand(self, self%o, self%wk(:, 1), self%wk(:, 2))
    call apply_h(self, self%wk(:, 2), self%wk(:, 3), self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))
    call apply_rand(self, self%e, self%wk(:, 7), self%wk(:, 8))
    !$omp parallel do private(l)
    do i = 1, self%nn
      do l = 1, 9
        y(9 * (i - 1) + l) = self%sq(i) * (self%wk(9 * (i - 1) + l, 8) + self%wk(9 * (i - 1) + l, 1) - self%wk(9 * (i - 1) + l, 3))
      end do
    end do
    !$omp end parallel do
  end subroutine augcls_apply

  !> y = h x = (c + D S D) x for complex x in normalised coordinates, one application of S per real part (route ed3; augcls_apply
  !> always runs apply_h twice, also for o = 0).  The term e of apply is not included.
  subroutine augcls_hp_c(self, x, y)
    class(augcls_op), intent(in) :: self
    complex(dp), intent(in) :: x(:)
    complex(dp), intent(out) :: y(:)
    integer :: i
    if (.not. associated(self%wk)) error stop 'augcls_hp_c: work space not allocated'
    !$omp parallel do
    do i = 1, self%nn
      self%wk(9 * (i - 1) + 1:9 * i, 7) = real(x(9 * (i - 1) + 1:9 * i), dp) / self%sq(i)
      self%wk(9 * (i - 1) + 1:9 * i, 8) = aimag(x(9 * (i - 1) + 1:9 * i)) / self%sq(i)
    end do
    !$omp end parallel do
    call apply_h(self, self%wk(:, 7), self%wk(:, 1), self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))
    call apply_h(self, self%wk(:, 8), self%wk(:, 2), self%wk(:, 4), self%wk(:, 5), self%wk(:, 6))
    !$omp parallel do
    do i = 1, self%nn
      y(9 * (i - 1) + 1:9 * i) = self%sq(i) * cmplx(self%wk(9 * (i - 1) + 1:9 * i, 1), self%wk(9 * (i - 1) + 1:9 * i, 2), dp)
    end do
    !$omp end parallel do
  end subroutine augcls_hp_c

  !> Projection of x (normalised coordinates, in place) onto the symmetric subspace V of the reduced operator.  A representative
  !> node r with stabiliser S_r = {g: g r = r} must satisfy psi(r) = chi(g) D(g)^T psi(r) for every g in S_r, so
  !>     (P psi)(r) = 1/|S_r| sum_{g in S_r} chi(g) D(g)^T psi(r).
  !> P is a real symmetric projector.  The operator maps V into V and is symmetric on V only.  Outside V it is some other,
  !> non-symmetric matrix.  A Lanczos chain started in V leaves it by rounding, and the recurrence amplifies that part by
  !> orders of magnitude within tens of steps.  So every chain projects its new vector (asr_edc, asr_ed3).  Nodes with the
  !> trivial stabiliser are not touched; the cost is far below one operator application.
  subroutine augcls_project_c(self, x)
    class(augcls_op), intent(in) :: self
    complex(dp), intent(inout) :: x(:)
    complex(dp) :: acc(9)
    integer :: i, j, l, g, ns
    !$omp parallel do private(j, l, g, ns, acc)
    do i = 1, self%nn
      ns = self%stab_ptr(i + 1) - self%stab_ptr(i)
      if (ns <= 1) cycle
      acc = (0.0_dp, 0.0_dp)
      do j = self%stab_ptr(i), self%stab_ptr(i + 1) - 1
        g = self%stab(j)
        do l = 1, 9
          acc(l) = acc(l) + self%csg(l, g) * x(9 * (i - 1) + self%pi(l, g))
        end do
      end do
      x(9 * (i - 1) + 1:9 * i) = acc / real(ns, dp)
    end do
    !$omp end parallel do
  end subroutine augcls_project_c

  !> |x - P x| / |x|: the part of x outside the symmetric subspace (diagnostic, tests)
  function augcls_leak_c(self, x) result(leak)
    class(augcls_op), intent(in) :: self
    complex(dp), intent(in) :: x(:)
    real(dp) :: leak
    complex(dp), allocatable :: y(:)
    y = x
    call augcls_project_c(self, y)
    leak = sqrt(sum(abs(x - y)**2) / max(sum(abs(x)**2), tiny(1.0_dp)))
  end function augcls_leak_c

  !> local weights wt(l, k, k', class) = sum_s par(l, spc(s)) mt(k, k', s, class)
  pure subroutine local_weights(self, par, wt)
    class(augcls_op), intent(in) :: self
    real(dp), intent(in) :: par(:, :)
    real(dp), intent(out) :: wt(:, :, :, :)
    integer :: lam, k, kk, s, kap
    wt = 0.0_dp
    do kap = 1, self%ncls
      do kk = 1, self%nsc(kap); do k = 1, self%nsc(kap)
        do s = 1, self%nsc(kap)
          do lam = 1, 9
            wt(lam, k, kk, kap) = wt(lam, k, kk, kap) + par(lam, self%spc(s, kap)) * self%mt(k, kk, s, kap)
          end do
        end do
      end do; end do
    end do
  end subroutine local_weights

  !> w = X x (amplitude coordinates), par(l, species)
  subroutine apply_rand(self, par, x, w)
    class(augcls_op), intent(in) :: self
    real(dp), intent(in) :: par(:, :), x(:)
    real(dp), intent(out) :: w(:)
    real(dp) :: wt(9, self%nsmax, self%nsmax, self%ncls), tf(9), acc(9)
    integer :: i, k, kap, l0
    call local_weights(self, par, wt)
    !$omp parallel do private(k, kap, l0, tf, acc)
    do i = 1, self%nn
      if (self%frozen(i) > 0) then
        w(9 * (i - 1) + 1:9 * i) = par(1:9, self%frozen(i)) * x(9 * (i - 1) + 1:9 * i)
        cycle
      end if
      kap = self%ncl(i); l0 = self%lst(i) + 1
      acc = wt(:, l0, l0, kap) * x(9 * (i - 1) + 1:9 * i)
      do k = 1, self%nsc(kap)
        if (k == l0 .or. self%flipr(k, i) == 0) cycle
        call at_link(self, x, self%flipr(k, i), self%flipg(k, i), tf)
        acc = acc + wt(:, l0, k, kap) * tf
      end do
      w(9 * (i - 1) + 1:9 * i) = acc
    end do
    !$omp end parallel do
  end subroutine apply_rand

  !> w = X x in normalised coordinates, real parameters par(l, species) (public)
  subroutine augcls_rand(self, par, x, w)
    class(augcls_op), intent(in) :: self
    real(dp), intent(in) :: par(:, :), x(:)
    real(dp), intent(out) :: w(:)
    integer :: i
    real(dp), allocatable :: t(:), u(:)
    allocate (t(size(x)), u(size(x)))
    do i = 1, self%nn; t(9 * (i - 1) + 1:9 * i) = x(9 * (i - 1) + 1:9 * i) / self%sq(i); end do
    call apply_rand(self, par, t, u)
    do i = 1, self%nn; w(9 * (i - 1) + 1:9 * i) = self%sq(i) * u(9 * (i - 1) + 1:9 * i); end do
  end subroutine augcls_rand

  !> complex parameters and vectors (energy-dependent diagonal e_s(E)), normalised coordinates: w = X x
  subroutine augcls_rand_c(self, par, x, w)
    class(augcls_op), intent(in) :: self
    complex(dp), intent(in) :: par(:, :), x(:)
    complex(dp), intent(out) :: w(:)
    real(dp) :: wr(9, self%nsmax, self%nsmax, self%ncls), wi(9, self%nsmax, self%nsmax, self%ncls)
    complex(dp) :: tf(9), acc(9), xi(9)
    integer :: i, k, kap, l0
    call local_weights(self, real(par), wr)
    call local_weights(self, aimag(par), wi)
    !$omp parallel do private(k, kap, l0, tf, acc, xi)
    do i = 1, self%nn
      xi = x(9 * (i - 1) + 1:9 * i) / self%sq(i)
      if (self%frozen(i) > 0) then
        w(9 * (i - 1) + 1:9 * i) = par(1:9, self%frozen(i)) * x(9 * (i - 1) + 1:9 * i)
        cycle
      end if
      kap = self%ncl(i); l0 = self%lst(i) + 1
      acc = cmplx(wr(:, l0, l0, kap), wi(:, l0, l0, kap), dp) * xi
      do k = 1, self%nsc(kap)
        if (k == l0 .or. self%flipr(k, i) == 0) cycle
        call at_link_c(self, x, self%flipr(k, i), self%flipg(k, i), tf)
        tf = tf / self%sq(self%flipr(k, i))
        acc = acc + cmplx(wr(:, l0, k, kap), wi(:, l0, k, kap), dp) * tf
      end do
      w(9 * (i - 1) + 1:9 * i) = self%sq(i) * acc
    end do
    !$omp end parallel do
  end subroutine augcls_rand_c

  !> hx = c x + D S D x   (amplitude coordinates)
  subroutine apply_h(self, x, hx, t1, t2, t3)
    class(augcls_op), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: hx(:), t1(:), t2(:), t3(:)
    integer :: i, k, l, l2, g, r
    real(dp) :: acc(9), tl
    call apply_rand(self, self%d, x, t1)                               ! t1 = D x
    ! t2 = S t1.  The inner loop runs over the contiguous first index of block0 and btab.  Each acc(l) still sums its terms in
    ! the order l2, then slot, and chi(g) = +-1 multiplies t1 exactly, so the result has the same bits as the l-outer form.
    !$omp parallel do private(k, l, l2, g, r, acc, tl)
    do i = 1, self%nn
      acc = 0.0_dp
      do l2 = 1, 9
        tl = t1(9 * (i - 1) + l2)
        do l = 1, 9; acc(l) = acc(l) + self%block0(l, l2) * tl; end do
      end do
      do k = 1, self%nslot
        r = self%hopr(k, i)
        if (r == 0) cycle
        g = self%hopg(k, i)
        do l2 = 1, 9
          tl = self%chi(g) * t1(9 * (r - 1) + l2)
          do l = 1, 9; acc(l) = acc(l) + self%btab(l, l2, k, g) * tl; end do
        end do
      end do
      t2(9 * (i - 1) + 1:9 * i) = acc
    end do
    !$omp end parallel do
    call apply_rand(self, self%d, t2, t3)                              ! t3 = D S D x
    call apply_rand(self, self%c, x, hx)                               ! hx = c x
    if (chain_vec_par) then
      !$omp parallel do
      do i = 1, size(hx)
        hx(i) = hx(i) + t3(i)
      end do
      !$omp end parallel do
    else
      hx = hx + t3
    end if
  end subroutine apply_h

  ! ------------------------------------------------------------------ canonical nodes
  !> g*(site, cfg) = canonical representative (csite, ccfg) over the nop operations of perm; returns that operation and the
  !> orbit size.  Config codes site*nsm + state map to perm(site)*nsm + state.
  subroutine canon_node(perm, nop, nsm, site, cfg, csite, ccfg, osize, gbest)
    integer, intent(in) :: perm(:, :), nop, nsm, site, cfg(:)
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
      do k = 1, n; c2(k) = perm(cfg(k) / nsm, iop) * nsm + mod(cfg(k), nsm); end do
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

  !> append the record rec to the buffer d (nrec entries used), growing it as needed
  pure subroutine cand_push(d, nrec, rec)
    integer, allocatable, intent(inout) :: d(:)
    integer, intent(inout) :: nrec
    integer, intent(in) :: rec(:)
    integer, allocatable :: t(:)
    if (nrec + size(rec) > size(d)) then
      allocate (t(max(2 * size(d), nrec + size(rec))))
      t(1:nrec) = d(1:nrec)
      call move_alloc(t, d)
    end if
    d(nrec + 1:nrec + size(rec)) = rec
    nrec = nrec + size(rec)
  end subroutine cand_push

  !> release the shared node tables and the work space of an operator from build_augcls_operator.  Call it on the original
  !> only: the copies share the tables (their own wk is released by whoever allocated it).
  subroutine augcls_free(op)
    type(augcls_op), intent(inout) :: op
    if (associated(op%flipr)) deallocate (op%flipr)
    if (associated(op%flipg)) deallocate (op%flipg)
    if (associated(op%hopr)) deallocate (op%hopr)
    if (associated(op%hopg)) deallocate (op%hopg)
    if (associated(op%lst)) deallocate (op%lst)
    if (associated(op%site)) deallocate (op%site)
    if (associated(op%ncl)) deallocate (op%ncl)
    if (associated(op%level)) deallocate (op%level)
    if (associated(op%norb)) deallocate (op%norb)
    if (associated(op%sq)) deallocate (op%sq)
    if (associated(op%stab_ptr)) deallocate (op%stab_ptr)
    if (associated(op%stab)) deallocate (op%stab)
    if (associated(op%wk)) deallocate (op%wk)
    op%nn = 0
  end subroutine augcls_free

  ! ------------------------------------------------------------------ builder
  !> Node set = closure of the support of H (paths R -> R1 -> R2, every visited site in any of its allowed states), rank-limited,
  !> on canonical representatives.  cls(site) = class of every graph site; perm(:, 1:nop), dmat(9, 9, nop): the group.
  !> tsite (optional): target graph sites.  For each target and local state k = 0 .. nsc - 1 (the site alone in state k, every
  !> other site in its reference state) tnode(k + 1, t) returns the representative node and tgop(k + 1, t) the operation that
  !> maps the node onto it (augcls_amp_c).  tnode = 0: the node is not in the node set (or the site cannot flip).  The node store
  !> exists only inside this routine, so the lookup is done here (asr_spectral).
  subroutine build_augcls_operator(g, mdl, cls, perm, dmat, start, maxlevel, op, maxrank, flippable, tsite, tnode, tgop)
    type(graph_t), intent(in) :: g
    type(class_alloy_t), intent(in) :: mdl
    integer, intent(in) :: cls(:), perm(:, :), start, maxlevel
    real(dp), intent(in) :: dmat(:, :, :)
    type(augcls_op), intent(out) :: op
    integer, intent(in), optional :: maxrank
    logical, intent(in), optional :: flippable(:)
    integer, intent(in), optional :: tsite(:)
    integer, allocatable, intent(out), optional :: tnode(:, :), tgop(:, :)
    type(store_t) :: st
    logical, allocatable :: canflip(:)
    integer, allocatable :: c0(:), cfg0(:), sites(:), cc(:), tmp(:), ccf(:), nst(:), c2(:), nstbn(:)
    integer :: nstb, nthr, i2, ii, nrec
    type :: cand_t                                       !! candidates of one node: records (site, orbit size, n, cfg(1:n))
      integer :: n = 0
      integer, allocatable :: d(:)
    end type cand_t
    type(cand_t), allocatable :: cb(:)
    integer :: nop, nsm, i, lev, r, r1, r2, k1, k2, nsit, f, j, idx, rkmax, cs, osz, gb, ns, k, l, gg, jj, kap, ncomb, q, dgt
    logical :: isnew

    nop = size(perm, 2)
    if (size(dmat, 3) < nop) error stop 'build_augcls_operator: dmat has fewer operations than perm'
    if (any(perm(:, 1) /= [(i, i=1, g%n)])) error stop 'build_augcls_operator: operation 1 must be the identity'
    do gg = 1, nop
      if (perm(start, gg) /= start) error stop 'build_augcls_operator: the group does not fix the start site'
      if (any(cls(perm(:, gg)) /= cls)) error stop 'build_augcls_operator: an operation does not preserve the site classes'
    end do
    nsm = maxval(mdl%nsc)
    rkmax = huge(1); if (present(maxrank)) rkmax = maxrank
    allocate (canflip(g%n)); canflip = .true.
    if (present(flippable)) canflip = flippable
    allocate (nst(g%n))
    do r = 1, g%n                                                      ! allowed local states per site
      nst(r) = 1; if (canflip(r)) nst(r) = mdl%nsc(cls(r))
    end do
    call store_init(st, 4096)
    allocate (cfg0(0))
    call canon_node(perm, nop, nsm, start, cfg0, cs, ccf, osz, gb)
    if (osz /= 1) error stop 'build_augcls_operator: start site is not invariant under the group'
    idx = find_or_add(st, cs, ccf, 0, isnew); st%orb(idx) = 1
    ! Breadth-first closure, in chunks of nodes that exist already.  The candidates of a chunk (canonical form of every reachable
    ! configuration) are computed in parallel: canon_node does not depend on the store.  They are then inserted serially in the
    ! order (node, k1, k2, f) of the plain serial loop, so every node gets the same index as there, for any number of threads.
    nthr = 1
    !$ nthr = omp_get_max_threads()
    allocate (cb(max(16, 4 * nthr)))
    i = 1
    do while (i <= st%n)
      i2 = min(st%n, i + size(cb) - 1)
      !$omp parallel do schedule(dynamic) private(lev, r, c0, k1, k2, r1, r2, sites, nsit, ncomb, f, cc, q, j, dgt, tmp, cs, ccf, &
      !$omp&                                       osz, gb, nrec)
      do ii = i, i2
        nrec = 0
        if (.not. allocated(cb(ii - i + 1)%d)) allocate (cb(ii - i + 1)%d(1024))
        lev = st%level(ii)
        if (lev < maxlevel) then
          r = st%site(ii)
          c0 = st%data(st%start(ii):st%start(ii) + st%len(ii) - 1)
          do k1 = 0, g%nbr_ptr(r + 1) - g%nbr_ptr(r)
            if (k1 == 0) then; r1 = r; else; r1 = g%nbr(g%nbr_ptr(r) + k1 - 1); end if
            do k2 = 0, g%nbr_ptr(r1 + 1) - g%nbr_ptr(r1)
              if (k2 == 0) then; r2 = r1; else; r2 = g%nbr(g%nbr_ptr(r1) + k2 - 1); end if
              call uniq3(r, r1, r2, sites, nsit)
              ncomb = product(nst(sites))
              do f = 0, ncomb - 1
                cc = c0; q = f
                do j = 1, nsit
                  if (nst(sites(j)) == 1) cycle
                  dgt = mod(q, nst(sites(j))); q = q / nst(sites(j))
                  tmp = set_state(cc, sites(j), dgt, nsm); cc = tmp
                end do
                if (size(cc) > rkmax) cycle
                call canon_node(perm, nop, nsm, r2, cc, cs, ccf, osz, gb)
                call cand_push(cb(ii - i + 1)%d, nrec, [cs, osz, size(ccf), ccf])
              end do
            end do
          end do
        end if
        cb(ii - i + 1)%n = nrec
      end do
      !$omp end parallel do
      do ii = i, i2                                                      ! serial insertion, in the order of the serial loop
        lev = st%level(ii); k = 0
        do while (k < cb(ii - i + 1)%n)
          cs = cb(ii - i + 1)%d(k + 1); osz = cb(ii - i + 1)%d(k + 2); nsit = cb(ii - i + 1)%d(k + 3)
          idx = find_or_add(st, cs, cb(ii - i + 1)%d(k + 4:k + 3 + nsit), lev + 1, isnew)
          if (isnew) st%orb(idx) = osz
          k = k + 3 + nsit
        end do
      end do
      i = i2 + 1
    end do

    ! ---- tables
    op%nn = st%n; op%nl = 9; op%nop = nop; op%nsmax = nsm; op%ncls = mdl%ncls; op%nsp = mdl%nsp
    op%nsc = mdl%nsc; op%spc = mdl%spc
    allocate (op%mt(nsm, nsm, nsm, mdl%ncls)); op%mt = 0.0_dp
    do kap = 1, mdl%ncls
      call class_mt(mdl%conc(1:mdl%nsc(kap), kap), op%mt(:, :, :, kap))
    end do
    ns = maxval(g%nbr_ptr(2:g%n + 1) - g%nbr_ptr(1:g%n)); op%nslot = ns
    allocate (op%flipr(nsm, st%n), op%flipg(nsm, st%n), op%hopr(ns, st%n), op%hopg(ns, st%n), op%lst(st%n), op%site(st%n), &
              op%ncl(st%n), op%frozen(st%n), op%level(st%n), op%norb(st%n), op%sq(st%n))
    op%flipr = 0; op%flipg = 1; op%hopr = 0; op%hopg = 1; op%frozen = 0
    ! species parameters: copied only if present (callers may set them later with augcls_set_params).  Assigning an unallocated
    ! array is undefined behaviour; ifx leaves the target unallocated, gfortran 14 crashes (test_augcls_blk, ring_case).
    if (allocated(mdl%c)) op%c = mdl%c
    if (allocated(mdl%d)) op%d = mdl%d
    if (allocated(mdl%o)) op%o = mdl%o
    if (allocated(mdl%e)) op%e = mdl%e
    allocate (op%block0(9, 9)); op%block0 = 0.0_dp
    if (allocated(mdl%block0)) op%block0 = mdl%block0
    allocate (op%pi(9, nop), op%sg(9, nop), op%chi(nop), op%csg(9, nop))
    do gg = 1, nop                                                     ! orbital operation tables: D[pi(l),l] = sg(l)
      do l = 1, 9
        op%pi(l, gg) = 0; op%sg(l, gg) = 0.0_dp
        do jj = 1, 9
          if (abs(dmat(jj, l, gg)) > 0.5_dp) then; op%pi(l, gg) = jj; op%sg(l, gg) = dmat(jj, l, gg); end if
        end do
        if (op%pi(l, gg) == 0) error stop 'build_augcls_operator: orbital operation is not monomial'
      end do
    end do
    allocate (op%btab(9, 9, ns, nop)); op%btab = 0.0_dp                ! btab(l, l2, slot, g) = sum_l' B[l,l'] sg(l',g) delta(l2, pi(l',g))
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
    op%csg = op%sg
    ! node tables: stabilisers counted in a first parallel pass (their offsets need a prefix sum), then every column i of the
    ! tables filled in parallel; the store is read only (find_state), so the result does not depend on the threads
    allocate (op%stab_ptr(st%n + 1), nstbn(st%n))
    !$omp parallel do schedule(dynamic, 64) private(r, c0, gg, c2)
    do i = 1, st%n
      r = st%site(i)
      c0 = st%data(st%start(i):st%start(i) + st%len(i) - 1)
      nstbn(i) = 0
      do gg = 1, nop
        if (perm(r, gg) /= r) cycle
        c2 = [(perm(c0(k) / nsm, gg) * nsm + mod(c0(k), nsm), k=1, size(c0))]
        call sort_small(c2)
        if (all(c2 == c0)) nstbn(i) = nstbn(i) + 1
      end do
      if (nstbn(i) * st%orb(i) /= nop) error stop 'build_augcls_operator: stabiliser size /= group order / orbit size'
    end do
    !$omp end parallel do
    op%stab_ptr(1) = 1
    do i = 1, st%n
      op%stab_ptr(i + 1) = op%stab_ptr(i) + nstbn(i)
    end do
    nstb = op%stab_ptr(st%n + 1) - 1
    allocate (op%stab(nstb))
    !$omp parallel do schedule(dynamic, 64) private(r, c0, gg, c2, j, k, tmp, cs, ccf, osz, gb)
    do i = 1, st%n
      r = st%site(i)
      c0 = st%data(st%start(i):st%start(i) + st%len(i) - 1)
      j = op%stab_ptr(i)                                               ! stabiliser: g maps the site and the configuration onto themselves
      do gg = 1, nop
        if (perm(r, gg) /= r) cycle
        c2 = [(perm(c0(k) / nsm, gg) * nsm + mod(c0(k), nsm), k=1, size(c0))]
        call sort_small(c2)
        if (all(c2 == c0)) then; op%stab(j) = gg; j = j + 1; end if
      end do
      op%level(i) = st%level(i); op%norb(i) = real(st%orb(i), dp); op%sq(i) = sqrt(op%norb(i)); op%site(i) = r
      op%ncl(i) = cls(r)
      op%lst(i) = get_state(c0, r, nsm)
      do k = 0, mdl%nsc(cls(r)) - 1
        if (k == op%lst(i)) then
          op%flipr(k + 1, i) = i; op%flipg(k + 1, i) = 1
        else if (k == 0 .or. canflip(r)) then
          tmp = set_state(c0, r, k, nsm)
          call canon_node(perm, nop, nsm, r, tmp, cs, ccf, osz, gb)
          j = find_state(st, cs, ccf)
          if (j > 0) then; op%flipr(k + 1, i) = j; op%flipg(k + 1, i) = gb; end if
        end if
      end do
      do k = 1, g%nbr_ptr(r + 1) - g%nbr_ptr(r)
        call canon_node(perm, nop, nsm, g%nbr(g%nbr_ptr(r) + k - 1), c0, cs, ccf, osz, gb)
        j = find_state(st, cs, ccf)
        if (j > 0) then; op%hopr(k, i) = j; op%hopg(k, i) = gb; end if
      end do
    end do
    !$omp end parallel do
    if (present(tsite)) then                                           ! target nodes (see the header)
      if (.not. (present(tnode) .and. present(tgop))) error stop 'build_augcls_operator: tsite needs tnode and tgop'
      allocate (tnode(nsm, size(tsite)), tgop(nsm, size(tsite))); tnode = 0; tgop = 1
      do i = 1, size(tsite)
        r = tsite(i)
        do k = 0, mdl%nsc(cls(r)) - 1
          if (k > 0 .and. (r == start .or. .not. canflip(r))) cycle
          tmp = set_state(cfg0, r, k, nsm)
          call canon_node(perm, nop, nsm, r, tmp, cs, ccf, osz, gb)
          tnode(k + 1, i) = find_state(st, cs, ccf); tgop(k + 1, i) = gb
        end do
      end do
    end if
    allocate (op%wk(9 * st%n, 8))
    call augcls_set_start(op, 1, dmat)
  end subroutine build_augcls_operator
end module asr_augcls
