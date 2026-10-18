!> Site-class (sublattice) resolved, point-group reduced augmented space for  H = E_nu + h - h o h,  h = c + d S d  (order 2),
!> in which every site class has its own species set and concentrations: ordered sublattices, partial (sublattice) disorder and
!> multicomponent (high-entropy) alloys in one operator.
!>
!> Site R of class kappa carries an independent random species s_R in {spc(1,kappa) .. spc(nsc(kappa),kappa)} (global species
!> indices) with probabilities conc(k, kappa).  Its local configuration space has the orthonormal basis phi_0 = 1, phi_1 ..
!> phi_{nsc-1} (asr_augn::augn_basis on the class concentrations), and a site-diagonal random parameter X_R becomes
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
  use asr_augn, only: augn_basis
  implicit none
  private
  public :: class_alloy_t, augcls_op, build_augcls_operator, augcls_set_start, augcls_freeze_site, augcls_rand, augcls_rand_c

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
    integer, allocatable :: flipr(:, :), flipg(:, :)              !! (nsmax, node): node with the site's own state set to k-1
    integer, allocatable :: hopr(:, :), hopg(:, :), lst(:), site(:), ncl(:), frozen(:), level(:)
    real(dp), allocatable :: norb(:), sq(:)
    integer, allocatable :: pi(:, :)
    real(dp), allocatable :: sg(:, :), chi(:)
    real(dp), pointer :: wk(:, :) => null()
  contains
    procedure :: apply => augcls_apply
    procedure :: dim => augcls_dim
  end type augcls_op

contains

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
      op%chi(g) = dmat(l0, l0, g)
      if (abs(abs(op%chi(g)) - 1.0_dp) > 1.0e-9_dp) error stop 'augcls_set_start: start orbital is not one-dimensional under the group'
    end do
  end subroutine augcls_set_start

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
      out(l) = op%chi(g) * op%sg(l, g) * v(9 * (r - 1) + op%pi(l, g))
    end do
  end subroutine at_link

  pure subroutine at_link_c(op, v, r, g, out)
    class(augcls_op), intent(in) :: op
    complex(dp), intent(in) :: v(:)
    integer, intent(in) :: r, g
    complex(dp), intent(out) :: out(9)
    integer :: l
    do l = 1, 9
      out(l) = op%chi(g) * op%sg(l, g) * v(9 * (r - 1) + op%pi(l, g))
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
    real(dp) :: acc(9)
    call apply_rand(self, self%d, x, t1)                               ! t1 = D x
    !$omp parallel do private(k, l, l2, g, r, acc)
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
    call apply_rand(self, self%d, t2, t3)                              ! t3 = D S D x
    call apply_rand(self, self%c, x, hx)                               ! hx = c x
    hx = hx + t3
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

  ! ------------------------------------------------------------------ builder
  !> Node set = closure of the support of H (paths R -> R1 -> R2, every visited site in any of its allowed states), rank-limited,
  !> on canonical representatives.  cls(site) = class of every graph site; perm(:, 1:nop), dmat(9, 9, nop): the group.
  subroutine build_augcls_operator(g, mdl, cls, perm, dmat, start, maxlevel, op, maxrank, flippable)
    type(graph_t), intent(in) :: g
    type(class_alloy_t), intent(in) :: mdl
    integer, intent(in) :: cls(:), perm(:, :), start, maxlevel
    real(dp), intent(in) :: dmat(:, :, :)
    type(augcls_op), intent(out) :: op
    integer, intent(in), optional :: maxrank
    logical, intent(in), optional :: flippable(:)
    type(store_t) :: st
    logical, allocatable :: canflip(:)
    integer, allocatable :: c0(:), cfg0(:), sites(:), cc(:), tmp(:), ccf(:), nst(:)
    real(dp), allocatable :: phi(:, :)
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
            idx = find_or_add(st, cs, ccf, lev + 1, isnew)
            if (isnew) st%orb(idx) = osz
          end do
        end do
      end do
      i = i + 1
    end do

    ! ---- tables
    op%nn = st%n; op%nl = 9; op%nop = nop; op%nsmax = nsm; op%ncls = mdl%ncls; op%nsp = mdl%nsp
    op%nsc = mdl%nsc; op%spc = mdl%spc
    allocate (op%mt(nsm, nsm, nsm, mdl%ncls)); op%mt = 0.0_dp
    do kap = 1, mdl%ncls
      call augn_basis(mdl%conc(1:mdl%nsc(kap), kap), phi)
      do l = 1, mdl%nsc(kap)
        do k = 1, mdl%nsc(kap)
          do jj = 1, mdl%nsc(kap); op%mt(l, k, jj, kap) = mdl%conc(jj, kap) * phi(l, jj) * phi(k, jj); end do
        end do
      end do
    end do
    ns = maxval(g%nbr_ptr(2:g%n + 1) - g%nbr_ptr(1:g%n)); op%nslot = ns
    allocate (op%flipr(nsm, st%n), op%flipg(nsm, st%n), op%hopr(ns, st%n), op%hopg(ns, st%n), op%lst(st%n), op%site(st%n), &
              op%ncl(st%n), op%frozen(st%n), op%level(st%n), op%norb(st%n), op%sq(st%n))
    op%flipr = 0; op%flipg = 1; op%hopr = 0; op%hopg = 1; op%frozen = 0
    op%c = mdl%c; op%d = mdl%d; op%o = mdl%o; op%e = mdl%e
    allocate (op%block0(9, 9)); op%block0 = 0.0_dp
    if (allocated(mdl%block0)) op%block0 = mdl%block0
    allocate (op%pi(9, nop), op%sg(9, nop), op%chi(nop))
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
    do i = 1, st%n
      r = st%site(i)
      c0 = st%data(st%start(i):st%start(i) + st%len(i) - 1)
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
    allocate (op%wk(9 * st%n, 8))
    call augcls_set_start(op, 1, dmat)
  end subroutine build_augcls_operator
end module asr_augcls
