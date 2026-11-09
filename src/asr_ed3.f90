!> Full-ASA energy-dependent recursion (route `ed3`, docs/DERIVATION_ed3.md) on the site-class operator asr_augcls.  For a central
!> species c (start site frozen to c, one recursion per (class, species) and start orbital) the secular matrix of the generalised
!> ASA problem, K(E) = E O - H, scaled with T_R = (D_c / D_s)^1/2 is
!>   K'(E) = T K T = et - h' + u h' + h' u - h' Qt h',   h' = g + S',   S' = D_c^1/2 S D_c^1/2,
!>   et_s = (D_c / D_s) eps_s,  u_s = o_s eps_s,  Qt_s = (D_s / D_c) (o_s - p_s eps_s),  eps_s = E - E_nu,s,  g_s = (D_c / D_s) c_s,
!> with o, p in the representation alpha (p = p^gamma + o^2).  Every coefficient is linear in E, so K' has no poles.  With
!> |u> = |1 L>, |v> = h' |u> and |a> = |u> + o_c |v>, T = 1 on the centre gives the gamma-representation moments directly:
!>   q0 = <a|K'^-1|a>,  q1 = <v|K'^-1|a>,  q2 = <v|K'^-1|v>.
!> Both vectors lie in the reference sector and in the same reduced basis, since h' |u> = g_c |u> + S' |u> is not random.  One
!> complex-symmetric Lanczos chain from each gives q0 and q2 as diagonal elements (error second order in the residual) and q1
!> twice by projection (first order); the mean of the two is used and their difference measures the tail.  So the sphere charge
!> q0 + p^gamma q2 needs no cross term.  docs/DERIVATION_ed3.md section 6 compares this with chains from |u>, S'|u>.
module asr_ed3
  use asr_kinds, only: dp
  use asr_augcls, only: class_alloy_t, augcls_op, augcls_hp_c, augcls_rand_c
  use asr_lmto_ham, only: potpar_t, to_alpha, o_alpha
  use asr_endep, only: ed_depth_tol, ed_nchain, ed_nlev, cf_species
  use asr_augblk, only: blk_work_t, blk_kmul3
  implicit none
  private
  public :: ed3_t, ed3_work_t, ed3_work_init, ed3_setup, ed3_start, ed3_blocks, ed3_kmul, ed3_coefs, ed3_blocks_blk

  integer, parameter :: orb_l(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]

  type :: ed3_t
    integer :: ns = 0, central = 0
    real(dp), allocatable :: dc(:), del(:, :), enu(:, :)   !! D_c (9); D_s, E_nu,s (9, ns)
    real(dp), allocatable :: o(:, :), p(:, :), g(:, :)     !! o_s, p_s (alpha representation), g_s (9, ns)
  end type ed3_t

  type :: ed3_work_t
    complex(dp), allocatable :: u(:), up(:), w(:), t1(:), t2(:), t3(:), t4(:)
    complex(dp), allocatable :: tv(:)                      !! projection vector of the running chain
    complex(dp), allocatable :: s0(:), av(:), vv(:)        !! |u>, |a> = |u> + o |v>, |v> of the current start orbital
    real(dp) :: na = 0.0_dp, nv = 0.0_dp                   !! norms of |a>, |v>
  end type ed3_work_t

contains

  subroutine ed3_work_init(wk, n)
    type(ed3_work_t), intent(out) :: wk
    integer, intent(in) :: n
    allocate (wk%u(n), wk%up(n), wk%w(n), wk%t1(n), wk%t2(n), wk%t3(n), wk%t4(n), wk%tv(n), wk%s0(n), wk%av(n), wk%vv(n))
  end subroutine ed3_work_init

  !> Operator parameters of `mdl` for the central species `central`: c = g_s, d = D_c^1/2, o = e = 0, so that augcls_hp_c is h'.
  !> pp(1:nsp) are the gamma-representation parameters, alpha the screening.
  subroutine ed3_setup(mdl, pp, alpha, central, ed)
    type(class_alloy_t), intent(inout) :: mdl
    type(potpar_t), intent(in) :: pp(:)
    real(dp), intent(in) :: alpha(0:2)
    integer, intent(in) :: central
    type(ed3_t), intent(out) :: ed
    type(potpar_t), allocatable :: q(:)
    real(dp) :: o(0:2), dc, ds
    integer :: ns, s, i, ll
    ns = size(pp)
    if (ns /= mdl%nsp) error stop 'ed3_setup: number of species parameters /= mdl%nsp'
    allocate (q(ns)); ed%ns = ns; ed%central = central
    allocate (ed%dc(9), ed%del(9, ns), ed%enu(9, ns), ed%o(9, ns), ed%p(9, ns), ed%g(9, ns))
    if (allocated(mdl%c)) deallocate (mdl%c, mdl%d, mdl%o, mdl%e)
    allocate (mdl%c(9, ns), mdl%d(9, ns), mdl%o(9, ns), mdl%e(9, ns))
    do s = 1, ns; q(s) = to_alpha(pp(s), alpha); end do
    do i = 1, 9
      ll = orb_l(i)
      dc = q(central)%sqrdel(ll)**2
      ed%dc(i) = dc
      do s = 1, ns
        o = o_alpha(q(s), alpha)
        ds = q(s)%sqrdel(ll)**2
        ed%del(i, s) = ds; ed%enu(i, s) = q(s)%enu(ll); ed%o(i, s) = o(ll); ed%p(i, s) = q(s)%p(ll)
        ed%g(i, s) = dc / ds * (q(s)%c(ll) - q(s)%enu(ll))
        mdl%c(i, s) = ed%g(i, s)
        mdl%d(i, s) = abs(q(central)%sqrdel(ll))
        mdl%o(i, s) = 0.0_dp
        mdl%e(i, s) = 0.0_dp
      end do
    end do
  end subroutine ed3_setup

  !> species coefficients of K'(z)
  pure subroutine ed3_coefs(ed, z, et, uu, qt)
    type(ed3_t), intent(in) :: ed
    complex(dp), intent(in) :: z
    complex(dp), intent(out) :: et(:, :), uu(:, :), qt(:, :)
    complex(dp) :: eps
    integer :: l, s
    do s = 1, ed%ns
      do l = 1, 9
        eps = z - ed%enu(l, s)
        et(l, s) = ed%dc(l) / ed%del(l, s) * eps
        uu(l, s) = ed%o(l, s) * eps
        qt(l, s) = ed%del(l, s) / ed%dc(l) * (ed%o(l, s) - ed%p(l, s) * eps)
      end do
    end do
  end subroutine ed3_coefs

  !> y = K' x = et x - h'x + u h'x + h'(u x - Qt h'x); two applications of S'.  x and y must differ from wk%t1 .. wk%t4.
  subroutine ed3_kmul(op, wk, et, uu, qt, x, y)
    type(augcls_op), intent(in) :: op
    type(ed3_work_t), intent(inout) :: wk
    complex(dp), intent(in) :: et(:, :), uu(:, :), qt(:, :), x(:)
    complex(dp), intent(out) :: y(:)
    call augcls_hp_c(op, x, wk%t1)
    call augcls_rand_c(op, qt, wk%t1, wk%t2)
    call augcls_rand_c(op, uu, x, wk%t3)
    wk%t3 = wk%t3 - wk%t2
    call augcls_hp_c(op, wk%t3, wk%t4)
    call augcls_rand_c(op, et, x, y)
    call augcls_rand_c(op, uu, wk%t1, wk%t2)
    y = y - wk%t1 + wk%t2 + wk%t4
  end subroutine ed3_kmul

  !> start vectors of the orbital iv = 9 (node - 1) + l0 (node 1 is the start node): |u>, |v> = h'|u>, |a> = |u> + o_c |v>.  The
  !> operator builds |v>, so it is in the normalised orbit coordinates of the reduced basis.
  subroutine ed3_start(op, wk, ed, iv)
    type(augcls_op), intent(in) :: op
    type(ed3_work_t), intent(inout) :: wk
    type(ed3_t), intent(in) :: ed
    integer, intent(in) :: iv
    wk%s0 = (0.0_dp, 0.0_dp); wk%s0(iv) = (1.0_dp, 0.0_dp)
    call augcls_hp_c(op, wk%s0, wk%vv)
    wk%av = wk%s0 + ed%o(iv, ed%central) * wk%vv
    wk%nv = sqrt(real(sum(wk%vv * wk%vv), dp)); wk%na = sqrt(real(sum(wk%av * wk%av), dp))
    if (wk%nv < 1.0e-14_dp .or. wk%na < 1.0e-14_dp) error stop 'ed3_start: zero start vector'
  end subroutine ed3_start

  !> One complex-symmetric Lanczos chain on H_eff = z - K' from wk%u (wk%u^T wk%u = 1): g = <start| (z - H_eff)^-1 |start> and
  !> f = <tv| (z - H_eff)^-1 |start>, with the tail of cf_species.  tnorm = |tv| sets the scale of the depth test on f.
  subroutine chain(op, wk, et, uu, qt, z, nlev, tnorm, g, f, m)
    type(augcls_op), intent(in) :: op
    type(ed3_work_t), intent(inout) :: wk
    complex(dp), intent(in) :: et(:, :), uu(:, :), qt(:, :), z
    integer, intent(in) :: nlev
    real(dp), intent(in) :: tnorm
    complex(dp), intent(out) :: g, f
    integer, intent(out) :: m
    complex(dp) :: a(nlev), b(nlev), c(nlev), gtry, ftry, gprev, fprev
    integer :: k
    logical :: closed
    a = (0.0_dp, 0.0_dp); b = (0.0_dp, 0.0_dp); c = (0.0_dp, 0.0_dp)
    gprev = (0.0_dp, 0.0_dp); fprev = (0.0_dp, 0.0_dp)
    wk%up = (0.0_dp, 0.0_dp); m = 0; closed = .false.
    do k = 1, nlev
      c(k) = sum(wk%tv * wk%u)
      call ed3_kmul(op, wk, et, uu, qt, wk%u, wk%w)
      wk%w = z * wk%u - wk%w
      a(k) = sum(wk%u * wk%w)
      m = k
      if (k == nlev) exit
      if (ed_depth_tol > 0.0_dp .and. k >= 8 .and. mod(k, 4) == 0) then
        call cf_species(a, b, c, k, z, gtry, ftry)
        if (k > 8) then
          if (abs(gtry - gprev) <= ed_depth_tol * abs(gtry) .and. &
              abs(ftry - fprev) <= ed_depth_tol * max(abs(ftry), tnorm * abs(gtry))) exit
        end if
        gprev = gtry; fprev = ftry
      end if
      wk%w = wk%w - a(k) * wk%u
      if (k > 1) wk%w = wk%w - b(k - 1) * wk%up
      b(k) = sqrt(sum(wk%w * wk%w))
      if (abs(b(k)) < 1.0e-12_dp) then; closed = .true.; exit; end if
      wk%up = wk%u
      wk%u = wk%w / b(k)
    end do
    !$omp atomic
    ed_nchain = ed_nchain + 1
    !$omp atomic
    ed_nlev = ed_nlev + m
    if (closed) then
      call cf_closed(a, b, c, m, z, g, f)
    else
      call cf_species(a, b, c, m, z, g, f)
    end if
  end subroutine chain

  !> the chain closed (the Krylov space is exhausted): the fraction is exact, no tail
  pure subroutine cf_closed(a, b, c, m, z, g, f)
    complex(dp), intent(in) :: a(:), b(:), c(:), z
    integer, intent(in) :: m
    complex(dp), intent(out) :: g, f
    complex(dp) :: gt(m), prod
    integer :: j, k
    gt(m) = 1.0_dp / (z - a(m))
    do j = m - 1, 1, -1
      gt(j) = 1.0_dp / (z - a(j) - b(j)**2 * gt(j + 1))
    end do
    g = gt(1)
    f = (0.0_dp, 0.0_dp); prod = gt(1)
    do k = 1, m
      f = f + c(k) * prod
      if (k < m) prod = prod * b(k) * gt(k + 1)
    end do
  end subroutine cf_closed

  !> The central blocks at z for the start vectors of ed3_start, from nlev steps of each chain, in the form of bands_core:
  !>   bk(0) = G_11,  bk(1) = (hG)_11,  bk(2) = (hGh)_11,  bk(3) = q1 from the |a> chain - q1 from the |v> chain (0 when exact),
  !> so that q0 = bk(0) + 2 o bk(1) + o^2 bk(2), q1 = bk(1) + o bk(2), q2 = bk(2).  iv is the start component, op carries
  !> ed3_setup's parameters with the start site frozen to the central species.
  subroutine ed3_blocks(op, wk, ed, z, iv, nlev, bk)
    type(augcls_op), intent(in) :: op
    type(ed3_work_t), intent(inout) :: wk
    type(ed3_t), intent(in) :: ed
    complex(dp), intent(in) :: z
    integer, intent(in) :: iv, nlev
    complex(dp), intent(out) :: bk(0:3)
    complex(dp) :: et(9, ed%ns), uu(9, ed%ns), qt(9, ed%ns), gt, fm, q0, q1a, q1b, q2
    real(dp) :: oc
    integer :: m
    oc = ed%o(iv, ed%central)
    call ed3_coefs(ed, z, et, uu, qt)
    wk%u = wk%av / wk%na; wk%tv = wk%vv
    call chain(op, wk, et, uu, qt, z, nlev, wk%nv, gt, fm, m)
    q0 = wk%na**2 * gt; q1a = wk%na * fm
    wk%u = wk%vv / wk%nv; wk%tv = wk%av
    call chain(op, wk, et, uu, qt, z, nlev, wk%na, gt, fm, m)
    q2 = wk%nv**2 * gt; q1b = wk%nv * fm
    fm = 0.5_dp * (q1a + q1b)
    bk(2) = q2
    bk(1) = fm - oc * q2
    bk(0) = q0 - 2.0_dp * oc * fm + oc**2 * q2
    bk(3) = q1a - q1b
  end subroutine ed3_blocks

  !> ed3_blocks for a block of energies zs(1:nb) in lock step (asr_augblk).  ws holds the start vectors of ed3_start, bw the block
  !> work space (bw%nb >= nb); bk(0:3, j) as ed3_blocks at zs(j).
  subroutine ed3_blocks_blk(op, ws, bw, ed, zs, iv, nlev, bk)
    type(augcls_op), intent(in) :: op
    type(ed3_work_t), intent(in) :: ws
    type(blk_work_t), intent(inout) :: bw
    type(ed3_t), intent(in) :: ed
    complex(dp), intent(in) :: zs(:)
    integer, intent(in) :: iv, nlev
    complex(dp), intent(out) :: bk(0:, :)
    real(dp) :: cf(9, ed%ns, 6), oc
    complex(dp) :: gt(size(zs)), fm(size(zs)), q0(size(zs)), q1a(size(zs)), q1b(size(zs)), q2(size(zs))
    integer :: l, s, nb
    nb = size(zs)
    if (bw%nb < nb) error stop 'ed3_blocks_blk: block work space too small'
    oc = ed%o(iv, ed%central)
    do s = 1, ed%ns                                   ! et = z ea + eb, u = z ua + ub, Qt = z qa + qb (ed3_coefs)
      do l = 1, 9
        cf(l, s, 1) = ed%dc(l) / ed%del(l, s); cf(l, s, 2) = -cf(l, s, 1) * ed%enu(l, s)
        cf(l, s, 3) = ed%o(l, s);              cf(l, s, 4) = -ed%o(l, s) * ed%enu(l, s)
        cf(l, s, 5) = -ed%del(l, s) / ed%dc(l) * ed%p(l, s)
        cf(l, s, 6) = ed%del(l, s) / ed%dc(l) * (ed%o(l, s) + ed%p(l, s) * ed%enu(l, s))
      end do
    end do
    call chain_blk(op, bw, cf, zs, nlev, ws%av / ws%na, ws%vv, ws%nv, gt, fm)
    q0 = ws%na**2 * gt; q1a = ws%na * fm
    call chain_blk(op, bw, cf, zs, nlev, ws%vv / ws%nv, ws%av, ws%na, gt, fm)
    q2 = ws%nv**2 * gt; q1b = ws%nv * fm
    fm = 0.5_dp * (q1a + q1b)
    bk(2, 1:nb) = q2
    bk(1, 1:nb) = fm - oc * q2
    bk(0, 1:nb) = q0 - 2.0_dp * oc * fm + oc**2 * q2
    bk(3, 1:nb) = q1a - q1b
  end subroutine ed3_blocks_blk

  !> chain for a block of energies: every member starts from s0 (s0^T s0 = 1); g_j = <s0|(z_j - H_eff)^-1|s0>,
  !> f_j = <tv|(z_j - H_eff)^-1|s0>, depth test and closure as in chain
  subroutine chain_blk(op, bw, cf, zs, nlev, s0, tv, tnorm, g, f)
    type(augcls_op), intent(in) :: op
    type(blk_work_t), intent(inout) :: bw
    real(dp), intent(in) :: cf(:, :, :), tnorm
    complex(dp), intent(in) :: zs(:), s0(:), tv(:)
    integer, intent(in) :: nlev
    complex(dp), intent(out) :: g(:), f(:)
    complex(dp) :: a(nlev, bw%nb), b(nlev, bw%nb), c(nlev, bw%nb), gtry, ftry, gprev(bw%nb), fprev(bw%nb), acc(bw%nb), zp(bw%nb)
    real(dp) :: keep(bw%nb)
    integer :: k, nb, j, i, m(bw%nb)
    logical :: act(bw%nb), closed(bw%nb)
    nb = size(zs)
    zp(1:nb) = zs; zp(nb + 1:) = zs(nb)                       ! idle members pad a short block (zero rows throughout)
    a = (0.0_dp, 0.0_dp); b = (0.0_dp, 0.0_dp); c = (0.0_dp, 0.0_dp)
    gprev = (0.0_dp, 0.0_dp); fprev = (0.0_dp, 0.0_dp)
    m = nlev; closed = .false.; act = .false.; act(1:nb) = .true.
    associate (u => bw%u, up => bw%up, w => bw%w)
      do i = 1, size(u, 2); u(1:nb, i) = s0(i); u(nb + 1:, i) = (0.0_dp, 0.0_dp); end do
      up = (0.0_dp, 0.0_dp)
      do k = 1, nlev
        acc = (0.0_dp, 0.0_dp)
        do i = 1, size(u, 2); acc = acc + tv(i) * u(:, i); end do
        do j = 1, nb
          if (act(j)) c(k, j) = acc(j)
        end do
        call blk_kmul3(op, bw, cf(:, :, 1), cf(:, :, 2), cf(:, :, 3), cf(:, :, 4), cf(:, :, 5), cf(:, :, 6), zp, u, w)
        acc = (0.0_dp, 0.0_dp)
        do i = 1, size(u, 2)
          w(:, i) = zp * u(:, i) - w(:, i)
          acc = acc + u(:, i) * w(:, i)
        end do
        do j = 1, nb
          if (act(j)) then; a(k, j) = acc(j); m(j) = k; end if
        end do
        if (k == nlev) exit
        if (ed_depth_tol > 0.0_dp .and. k >= 8 .and. mod(k, 4) == 0) then
          do j = 1, nb
            if (.not. act(j)) cycle
            call cf_species(a(:, j), b(:, j), c(:, j), k, zs(j), gtry, ftry)
            if (k > 8) then
              if (abs(gtry - gprev(j)) <= ed_depth_tol * abs(gtry) .and. &
                  abs(ftry - fprev(j)) <= ed_depth_tol * max(abs(ftry), tnorm * abs(gtry))) act(j) = .false.
            end if
            gprev(j) = gtry; fprev(j) = ftry
          end do
        end if
        do i = 1, size(u, 2)
          w(:, i) = w(:, i) - a(k, :) * u(:, i)
          if (k > 1) w(:, i) = w(:, i) - b(k - 1, :) * up(:, i)
        end do
        acc = (0.0_dp, 0.0_dp)
        do i = 1, size(w, 2); acc = acc + w(:, i) * w(:, i); end do
        do j = 1, bw%nb
          if (.not. act(j)) cycle
          b(k, j) = sqrt(acc(j))
          if (abs(b(k, j)) < 1.0e-12_dp) then; act(j) = .false.; closed(j) = .true.; end if
        end do
        if (.not. any(act)) exit
        do j = 1, bw%nb
          if (act(j)) then
            acc(j) = 1.0_dp / b(k, j); keep(j) = 1.0_dp
          else
            acc(j) = (0.0_dp, 0.0_dp); keep(j) = 0.0_dp
          end if
        end do
        do i = 1, size(u, 2)
          up(:, i) = keep * u(:, i)
          u(:, i) = acc * w(:, i)
        end do
      end do
    end associate
    do j = 1, nb
      if (closed(j)) then
        call cf_closed(a(:, j), b(:, j), c(:, j), m(j), zs(j), g(j), f(j))
      else
        call cf_species(a(:, j), b(:, j), c(:, j), m(j), zs(j), g(j), f(j))
      end if
    end do
    !$omp atomic
    ed_nchain = ed_nchain + nb
    !$omp atomic
    ed_nlev = ed_nlev + sum(m(1:nb))
  end subroutine chain_blk
end module asr_ed3
