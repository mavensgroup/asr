!> Augmented space for a binary alloy (single orbital per site).
!>
!> Random site occupation n_R = 1 (A, probability x) or 0 (B, probability y = 1-x).  Each n_R is
!> replaced by the operator  M_R = x + (y-x) P_down_R + sqrt(xy) T_R  acting on the configuration
!> space, whose basis states are labelled by the set C of sites in the "down" (flipped) state; the
!> reference state is C = {} (all up).  Then  <<f(n)>> = <{}| f(M) |{}>  exactly (augmented space theorem).
!>
!> Model:  H = sum_R eps(n_R)|R><R| + sum_<RR'> t(n_R,n_R') (|R><R'| + h.c.),
!>   eps(n) = eps_B + (eps_A - eps_B) n,
!>   t(n,n') = t_BB + (t_AB - t_BB)(n + n') + (t_AA - 2 t_AB + t_BB) n n'.
!> States of the augmented space are pairs (R, C).  The reachable set from a start state is generated
!> breadth-first to a given level and the truncated, exactly symmetric operator is returned as a CSR
!> operator for the Lanczos recursion (exact for the first `maxlevel` recursion levels).
module asr_augspace
  use asr_kinds, only: dp, i8
  use asr_operator, only: csr_op, csr_from_triplets
  implicit none
  private
  public :: graph_t, binary_model_t, build_augmented_operator, build_augmented_operator_sym, aug_state_count
  public :: store_t, store_init, find_or_add, find_state, toggled, is_down, graph_distances

  type :: graph_t
    integer :: n = 0
    integer, allocatable :: nbr_ptr(:), nbr(:)        !! CSR neighbour lists (symmetric graph)
  end type graph_t

  type :: binary_model_t
    real(dp) :: x = 0.5_dp                            !! concentration of A
    real(dp) :: epsa = 0.0_dp, epsb = 0.0_dp
    real(dp) :: tAA = -1.0_dp, tAB = -1.0_dp, tBB = -1.0_dp
  end type binary_model_t

  type :: store_t
    integer :: n = 0, ndata = 0
    integer, allocatable :: site(:), start(:), len(:), data(:), level(:), orb(:)
    integer, allocatable :: table(:)
    integer :: tmask = 0
  end type store_t

contains

  integer function aug_state_count(op)
    type(csr_op), intent(in) :: op
    aug_state_count = op%n
  end function aug_state_count

  ! ---------------------------------------------------------------- hashed state store
  subroutine store_init(st, cap)
    type(store_t), intent(out) :: st
    integer, intent(in) :: cap
    integer :: tsz
    allocate (st%site(cap), st%start(cap), st%len(cap), st%level(cap), st%orb(cap), st%data(8 * cap))
    st%orb = 1
    tsz = 1024
    do while (tsz < 4 * cap); tsz = 2 * tsz; end do
    allocate (st%table(0:tsz - 1)); st%table = 0; st%tmask = tsz - 1
  end subroutine store_init

  pure integer function hash_key(site, cfg)
    integer, intent(in) :: site, cfg(:)
    integer(i8) :: h
    integer :: i
    h = 1469598103_i8 + int(site, i8) * 1099511_i8
    do i = 1, size(cfg)
      h = ieor(h, int(cfg(i), i8) + 1234567_i8)
      h = iand(h * 16777619_i8, 2147483647_i8)
    end do
    hash_key = int(iand(h, 2147483647_i8))
  end function hash_key

  !> Index of state (site, cfg), adding it (with `lev`) if new.
  function find_or_add(st, site, cfg, lev, isnew) result(idx)
    type(store_t), intent(inout) :: st
    integer, intent(in) :: site, cfg(:), lev
    logical, intent(out) :: isnew
    integer :: idx, h, j
    h = iand(hash_key(site, cfg), st%tmask)
    do
      j = st%table(h)
      if (j == 0) exit
      if (st%site(j) == site .and. st%len(j) == size(cfg)) then
        if (all(st%data(st%start(j):st%start(j) + st%len(j) - 1) == cfg)) then
          idx = j; isnew = .false.; return
        end if
      end if
      h = iand(h + 1, st%tmask)
    end do
    ! add
    if (st%n + 1 > size(st%site)) call grow_states(st)
    if (st%ndata + size(cfg) > size(st%data)) call grow_data(st, size(cfg))
    st%n = st%n + 1
    st%site(st%n) = site; st%len(st%n) = size(cfg); st%level(st%n) = lev
    st%start(st%n) = st%ndata + 1
    if (size(cfg) > 0) st%data(st%ndata + 1:st%ndata + size(cfg)) = cfg
    st%ndata = st%ndata + size(cfg)
    st%table(h) = st%n
    idx = st%n; isnew = .true.
    if (4 * st%n > st%tmask) call rehash(st)
  end function find_or_add

  !> Index of an existing state, 0 if absent (never adds).
  function find_state(st, site, cfg) result(idx)
    type(store_t), intent(in) :: st
    integer, intent(in) :: site, cfg(:)
    integer :: idx, h, j
    idx = 0
    h = iand(hash_key(site, cfg), st%tmask)
    do
      j = st%table(h)
      if (j == 0) return
      if (st%site(j) == site .and. st%len(j) == size(cfg)) then
        if (all(st%data(st%start(j):st%start(j) + st%len(j) - 1) == cfg)) then; idx = j; return; end if
      end if
      h = iand(h + 1, st%tmask)
    end do
  end function find_state

  !> Breadth-first graph distance (in neighbour steps) from `start` to every site.
  subroutine graph_distances(g, start, dist)
    type(graph_t), intent(in) :: g
    integer, intent(in) :: start
    integer, allocatable, intent(out) :: dist(:)
    integer, allocatable :: queue(:)
    integer :: qh, qt, r, k, j
    allocate (dist(g%n), queue(g%n)); dist = -1
    dist(start) = 0; queue(1) = start; qh = 1; qt = 1
    do while (qh <= qt)
      r = queue(qh); qh = qh + 1
      do k = g%nbr_ptr(r), g%nbr_ptr(r + 1) - 1
        j = g%nbr(k)
        if (dist(j) < 0) then; dist(j) = dist(r) + 1; qt = qt + 1; queue(qt) = j; end if
      end do
    end do
  end subroutine graph_distances

  subroutine grow_states(st)
    type(store_t), intent(inout) :: st
    integer, allocatable :: t(:)
    integer :: m
    m = 2 * size(st%site)
    call move_alloc_int(st%site, m); call move_alloc_int(st%start, m); call move_alloc_int(st%len, m); call move_alloc_int(st%level, m); call move_alloc_int(st%orb, m)
  end subroutine grow_states

  subroutine move_alloc_int(a, m)
    integer, allocatable, intent(inout) :: a(:)
    integer, intent(in) :: m
    integer, allocatable :: t(:)
    allocate (t(m)); t = 0; t(1:size(a)) = a; call move_alloc(t, a)
  end subroutine move_alloc_int

  subroutine grow_data(st, need)
    type(store_t), intent(inout) :: st
    integer, intent(in) :: need
    integer, allocatable :: t(:)
    allocate (t(2 * size(st%data) + need)); t = 0; t(1:size(st%data)) = st%data
    call move_alloc(t, st%data)
  end subroutine grow_data

  subroutine rehash(st)
    type(store_t), intent(inout) :: st
    integer :: tsz, i, h
    tsz = 2 * (st%tmask + 1)
    deallocate (st%table); allocate (st%table(0:tsz - 1)); st%table = 0; st%tmask = tsz - 1
    do i = 1, st%n
      h = iand(hash_key(st%site(i), st%data(st%start(i):st%start(i) + st%len(i) - 1)), st%tmask)
      do while (st%table(h) /= 0); h = iand(h + 1, st%tmask); end do
      st%table(h) = i
    end do
  end subroutine rehash

  !> cfg with site r toggled (inserted/removed), kept sorted.
  pure function toggled(cfg, r) result(out)
    integer, intent(in) :: cfg(:), r
    integer, allocatable :: out(:)
    integer :: n, k, i
    n = size(cfg)
    k = 0
    do i = 1, n
      if (cfg(i) == r) then; k = i; exit; end if
    end do
    if (k > 0) then
      allocate (out(n - 1)); out(1:k - 1) = cfg(1:k - 1); out(k:) = cfg(k + 1:)
    else
      allocate (out(n + 1))
      k = n + 1
      do i = 1, n
        if (cfg(i) > r) then; k = i; exit; end if
      end do
      out(1:k - 1) = cfg(1:k - 1); out(k) = r; out(k + 1:) = cfg(k:)
    end if
  end function toggled

  pure logical function is_down(cfg, r)
    integer, intent(in) :: cfg(:), r
    is_down = any(cfg == r)
  end function is_down

  ! ---------------------------------------------------------------- the augmented operator
  !> Breadth-first augmented space from |start, {}> to `maxlevel`; returns the symmetric truncated operator
  !> and the index of the start state (1).
  subroutine build_augmented_operator(g, mdl, start, maxlevel, op)
    type(graph_t), intent(in) :: g
    type(binary_model_t), intent(in) :: mdl
    integer, intent(in) :: start, maxlevel
    type(csr_op), intent(out) :: op
    type(store_t) :: st
    integer, allocatable :: ti(:), tj(:), cfg_loc(:)
    real(dp), allocatable :: tv(:)
    integer :: nt, i, r, k, rp, jdx, cur
    real(dp) :: x, y, s, mr, mrp, a, b, de, v
    integer, allocatable :: c0(:), c1(:), c2(:), c3(:)
    logical :: isnew
    integer :: lev

    x = mdl%x; y = 1.0_dp - x; s = sqrt(x * y)
    de = mdl%epsa - mdl%epsb
    a = mdl%tAB - mdl%tBB
    b = mdl%tAA - 2.0_dp * mdl%tAB + mdl%tBB
    call store_init(st, 4096)
    allocate (ti(65536), tj(65536), tv(65536)); nt = 0
    allocate (cfg_loc(0))
    cur = find_or_add(st, start, cfg_loc, 0, isnew)
    i = 1
    do while (i <= st%n)
      lev = st%level(i)
      if (lev >= maxlevel) then; i = i + 1; cycle; end if       ! not expanded
      r = st%site(i)
      c0 = st%data(st%start(i):st%start(i) + st%len(i) - 1)
      mr = x; if (is_down(c0, r)) mr = y
      ! on-site:  eps_B + de M_R
      call addv(i, r, c0, mdl%epsb + de * mr, lev + 1)
      c1 = toggled(c0, r)
      call addv(i, r, c1, de * s, lev + 1)
      ! hopping to neighbours
      do k = g%nbr_ptr(r), g%nbr_ptr(r + 1) - 1
        rp = g%nbr(k)
        mrp = x; if (is_down(c0, rp)) mrp = y
        c1 = toggled(c0, r); c2 = toggled(c0, rp); c3 = toggled(c1, rp)
        call addv(i, rp, c0, mdl%tBB + a * (mr + mrp) + b * mr * mrp, lev + 1)
        call addv(i, rp, c1, a * s + b * s * mrp, lev + 1)
        call addv(i, rp, c2, a * s + b * mr * s, lev + 1)
        call addv(i, rp, c3, b * s * s, lev + 1)
      end do
      i = i + 1
    end do
    ! symmetric completion: entries (j,i) for generated (i,j) whose row state i was not expanded
    block
      integer :: ntold, ri_c, cj_c
      real(dp) :: v_c
      ntold = nt
      do k = 1, ntold
        if (st%level(ti(k)) >= maxlevel .and. ti(k) /= tj(k)) then
          ri_c = tj(k); cj_c = ti(k); v_c = tv(k)            ! copies: put_raw may reallocate ti/tj/tv
          call put_raw(ri_c, cj_c, v_c)
        end if
      end do
    end block
    call csr_from_triplets(st%n, ti(1:nt), tj(1:nt), tv(1:nt), op)
  contains
    !> H|src> contains val |(site,cfg)>: add the target state only if the coupling is nonzero
    subroutine addv(src, tsite, tcfg, val, lv)
      integer, intent(in) :: src, tsite, tcfg(:), lv
      real(dp), intent(in) :: val
      integer :: ridx
      logical :: nw
      if (abs(val) < 1.0e-300_dp) return
      ridx = find_or_add(st, tsite, tcfg, lv, nw)
      call put_raw(ridx, src, val)
    end subroutine addv
    subroutine put_raw(ri, cj, val)
      integer, intent(in) :: ri, cj
      real(dp), intent(in) :: val
      integer, allocatable :: t1(:), t2(:)
      real(dp), allocatable :: t3(:)
      if (nt + 1 > size(ti)) then
        allocate (t1(2 * size(ti)), t2(2 * size(ti)), t3(2 * size(ti)))
        t1(1:nt) = ti(1:nt); t2(1:nt) = tj(1:nt); t3(1:nt) = tv(1:nt)
        call move_alloc(t1, ti); call move_alloc(t2, tj); call move_alloc(t3, tv)
      end if
      nt = nt + 1; ti(nt) = ri; tj(nt) = cj; tv(nt) = val
    end subroutine put_raw
  end subroutine build_augmented_operator

  ! ---------------------------------------------------------------- point-group reduced operator
  !> Canonical representative of the orbit of (site, cfg) under the operations perm(:, iop), the
  !> orbit size, and the canonical form (csite, ccfg).
  subroutine canonical_state(perm, site, cfg, csite, ccfg, osize)
    integer, intent(in) :: perm(:, :), site, cfg(:)
    integer, intent(out) :: csite, osize
    integer, allocatable, intent(out) :: ccfg(:)
    integer :: nop, iop, k, n, s2, ntie
    integer, allocatable :: c2(:)
    logical :: better, equal
    nop = size(perm, 2); n = size(cfg)
    allocate (c2(n), ccfg(n))
    csite = huge(1); ccfg = huge(1); ntie = 0
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
        csite = s2; ccfg = c2; ntie = 1
      else if (equal) then
        ntie = ntie + 1
      end if
    end do
    osize = nop / ntie
  end subroutine canonical_state

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

  !> Same model as build_augmented_operator, in the basis of point-group orbit sums
  !> |O> = |O|^(-1/2) sum_{s in O} |s>.  The start state (start site fixed by the group, empty
  !> configuration) is its own orbit.  H_red[O_i,O_j] = sqrt(|O_j|/|O_i|) sum_{t in O_i} H_{t, s_j}.
  subroutine build_augmented_operator_sym(g, mdl, perm, start, maxlevel, op, maxrank)
    type(graph_t), intent(in) :: g
    type(binary_model_t), intent(in) :: mdl
    integer, intent(in) :: perm(:, :), start, maxlevel
    type(csr_op), intent(out) :: op
    integer, intent(in), optional :: maxrank      !! drop states with more than this many flipped sites (rank truncation)
    integer :: rkmax
    type(store_t) :: st
    integer, allocatable :: ti(:), tj(:), c0(:), c1(:), c2(:), c3(:), cfg_loc(:), cs(:)
    real(dp), allocatable :: tv(:)
    integer :: nt, i, r, k, rp, jdx, lev, csite, osz
    real(dp) :: x, y, s, mr, mrp, a, b, de
    logical :: isnew

    x = mdl%x; y = 1.0_dp - x; s = sqrt(x * y)
    de = mdl%epsa - mdl%epsb
    a = mdl%tAB - mdl%tBB
    b = mdl%tAA - 2.0_dp * mdl%tAB + mdl%tBB
    rkmax = huge(1); if (present(maxrank)) rkmax = maxrank
    call store_init(st, 4096)
    allocate (ti(65536), tj(65536), tv(65536)); nt = 0
    allocate (cfg_loc(0))
    call canonical_state(perm, start, cfg_loc, csite, cs, osz)
    if (osz /= 1) error stop 'build_augmented_operator_sym: start site is not invariant under the group'
    jdx = find_or_add(st, csite, cs, 0, isnew)
    st%orb(jdx) = 1
    i = 1
    do while (i <= st%n)
      lev = st%level(i)
      if (lev >= maxlevel) then; i = i + 1; cycle; end if
      r = st%site(i)
      c0 = st%data(st%start(i):st%start(i) + st%len(i) - 1)
      mr = x; if (is_down(c0, r)) mr = y
      c1 = toggled(c0, r)
      call add(i, r, c0, mdl%epsb + de * mr)
      call add(i, r, c1, de * s)
      do k = g%nbr_ptr(r), g%nbr_ptr(r + 1) - 1
        rp = g%nbr(k)
        mrp = x; if (is_down(c0, rp)) mrp = y
        c2 = toggled(c0, rp); c3 = toggled(c1, rp)
        call add(i, rp, c0, mdl%tBB + a * (mr + mrp) + b * mr * mrp)
        call add(i, rp, c1, a * s + b * s * mrp)
        call add(i, rp, c2, a * s + b * mr * s)
        call add(i, rp, c3, b * s * s)
      end do
      i = i + 1
    end do
    block
      integer :: ntold, ri_c, cj_c
      real(dp) :: v_c
      ntold = nt
      do k = 1, ntold
        if (st%level(ti(k)) >= maxlevel .and. ti(k) /= tj(k)) then
          ri_c = tj(k); cj_c = ti(k); v_c = tv(k)            ! copies: put_raw may reallocate ti/tj/tv
          call put_raw(ri_c, cj_c, v_c)
        end if
      end do
    end block
    call csr_from_triplets(st%n, ti(1:nt), tj(1:nt), tv(1:nt), op)
  contains
    subroutine add(src, tsite, tcfg, val)
      integer, intent(in) :: src, tsite, tcfg(:)
      real(dp), intent(in) :: val
      integer :: cst, osz2, idx2
      integer, allocatable :: ccf(:)
      logical :: nw
      if (abs(val) < 1.0e-300_dp) return
      if (size(tcfg) > rkmax) return                 ! rank truncation: sites beyond are treated as average medium
      call canonical_state(perm, tsite, tcfg, cst, ccf, osz2)
      idx2 = find_or_add(st, cst, ccf, st%level(src) + 1, nw)
      if (nw) st%orb(idx2) = osz2
      call put_raw(idx2, src, val * sqrt(real(st%orb(src), dp) / real(st%orb(idx2), dp)))
    end subroutine add
    subroutine put_raw(ri, cj, val)
      integer, intent(in) :: ri, cj
      real(dp), intent(in) :: val
      integer, allocatable :: t1(:), t2(:)
      real(dp), allocatable :: t3(:)
      if (nt + 1 > size(ti)) then
        allocate (t1(2 * size(ti)), t2(2 * size(ti)), t3(2 * size(ti)))
        t1(1:nt) = ti(1:nt); t2(1:nt) = tj(1:nt); t3(1:nt) = tv(1:nt)
        call move_alloc(t1, ti); call move_alloc(t2, tj); call move_alloc(t3, tv)
      end if
      nt = nt + 1; ti(nt) = ri; tj(nt) = cj; tv(nt) = val
    end subroutine put_raw
  end subroutine build_augmented_operator_sym
end module asr_augspace
