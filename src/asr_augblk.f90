!> Block kernels on the site-class operator asr_augcls for the energy-dependent chains (routes edc and ed3).  A block holds nb
!> complex vectors in the layout x(j, i): j = 1..nb is the member (one contour energy each), i = 1..9 nn the component in the
!> normalised orbit coordinates of augcls_op.  One pass over the operator tables then serves all members: the neighbour indices,
!> the bond blocks and the gathered amplitudes are read once per node instead of once per real vector.
!>
!> The arithmetic is that of augcls_apply, augcls_hp_c and augcls_rand_c; only the order of the products inside one random weight
!> differs (rounding level).  Two exact simplifications:
!>  * a species parameter that is the same for every species (d in both routes, since S' = D_c^1/2 S D_c^1/2) gives W = X 1 on
!>    every node, because the local basis is orthonormal; it is applied as a diagonal product.  Local states of zero weight form
!>    a subspace that no operator couples to the rest, so the start vector never reaches them either way.
!>  * a zero parameter (o in ed3, e in both routes) is skipped.
!> Every energy-dependent species coefficient of both routes is linear in z, X(z) = z A + B, so it is applied as
!> z_j (A x)_j + (B x)_j with real A and B in one pass (blk_rand_lin).
module asr_augblk
  use asr_kinds, only: dp
  use asr_augcls, only: augcls_op
  implicit none
  private
  public :: blk_work_t, blk_work_init, blk_hamx, blk_kmul3, blk_hp, blk_rand_lin, blk_nbmax, blk_project

  !> work space of one thread: u, up, w are the Lanczos vectors of the members, the rest are kernel temporaries.  The kernels
  !> always work on all wk%nb members (whole, contiguous arrays); a shorter block is padded by the caller with idle zero rows.
  type :: blk_work_t
    integer :: nb = 0, n = 0
    complex(dp), allocatable :: u(:, :), up(:, :), w(:, :)
    complex(dp), allocatable :: x0(:, :), h1(:, :), r1(:, :), h2(:, :), e1(:, :), t1(:, :), t2(:, :)
  end type blk_work_t

contains

  !> the largest block size nb <= nreq whose work space (narr vectors of n complex numbers per member) fits into mb megabytes
  pure integer function blk_nbmax(nreq, n, narr, mb)
    integer, intent(in) :: nreq, n, narr
    real(dp), intent(in) :: mb
    blk_nbmax = max(1, min(nreq, int(mb * 1.0e6_dp / (16.0_dp * real(narr, dp) * real(max(n, 1), dp)))))
  end function blk_nbmax

  subroutine blk_work_init(wk, nb, n)
    type(blk_work_t), intent(out) :: wk
    integer, intent(in) :: nb, n
    wk%nb = nb; wk%n = n
    allocate (wk%u(nb, n), wk%up(nb, n), wk%w(nb, n), wk%x0(nb, n), wk%h1(nb, n), wk%r1(nb, n), wk%h2(nb, n), wk%e1(nb, n), &
              wk%t1(nb, n), wk%t2(nb, n))
  end subroutine blk_work_init

  pure logical function uniform(par)
    real(dp), intent(in) :: par(:, :)
    integer :: s
    uniform = .true.
    do s = 2, size(par, 2)
      if (any(par(:, s) /= par(:, 1))) then; uniform = .false.; return; end if
    end do
  end function uniform

  !> local weights wt(l, k, k', class) = sum_s par(l, spc(s)) mt(k, k', s, class) (as asr_augcls)
  pure subroutine weights(op, par, wt)
    type(augcls_op), intent(in) :: op
    real(dp), intent(in) :: par(:, :)
    real(dp), intent(out) :: wt(:, :, :, :)
    integer :: lam, k, kk, s, kap
    wt = 0.0_dp
    do kap = 1, op%ncls
      do kk = 1, op%nsc(kap); do k = 1, op%nsc(kap)
        do s = 1, op%nsc(kap)
          do lam = 1, 9
            wt(lam, k, kk, kap) = wt(lam, k, kk, kap) + par(lam, op%spc(s, kap)) * op%mt(k, kk, s, kap)
          end do
        end do
      end do; end do
    end do
  end subroutine weights

  !> amplitude coordinates x0 = x / sq
  subroutine to_amp(op, x, x0)
    type(augcls_op), intent(in) :: op
    complex(dp), intent(in), contiguous :: x(:, :)
    complex(dp), intent(out), contiguous :: x0(:, :)
    integer :: i, j0, l
    !$omp parallel do private(j0, l)
    do i = 1, op%nn
      j0 = 9 * (i - 1)
      do l = 1, 9; x0(:, j0 + l) = x(:, j0 + l) / op%sq(i); end do
    end do
    !$omp end parallel do
  end subroutine to_amp

  !> w = X x, real species parameter par(l, species), amplitude coordinates
  subroutine rand_r(op, par, x, w)
    type(augcls_op), intent(in) :: op
    real(dp), intent(in) :: par(:, :)
    complex(dp), intent(in), contiguous :: x(:, :)
    complex(dp), intent(out), contiguous :: w(:, :)
    real(dp) :: wt(9, op%nsmax, op%nsmax, op%ncls), f
    integer :: i, k, kap, l0, l, r, g, j0, jr
    if (uniform(par)) then
      !$omp parallel do private(j0, l)
      do i = 1, op%nn
        j0 = 9 * (i - 1)
        do l = 1, 9; w(:, j0 + l) = par(l, 1) * x(:, j0 + l); end do
      end do
      !$omp end parallel do
      return
    end if
    call weights(op, par, wt)
    !$omp parallel do private(k, kap, l0, l, r, g, j0, jr, f)
    do i = 1, op%nn
      j0 = 9 * (i - 1)
      if (op%frozen(i) > 0) then
        do l = 1, 9; w(:, j0 + l) = par(l, op%frozen(i)) * x(:, j0 + l); end do
        cycle
      end if
      kap = op%ncl(i); l0 = op%lst(i) + 1
      do l = 1, 9; w(:, j0 + l) = wt(l, l0, l0, kap) * x(:, j0 + l); end do
      do k = 1, op%nsc(kap)
        if (k == l0 .or. op%flipr(k, i) == 0) cycle
        r = op%flipr(k, i); g = op%flipg(k, i); jr = 9 * (r - 1)
        do l = 1, 9
          f = wt(l, l0, k, kap) * op%chi(g) * op%sg(l, g)
          w(:, j0 + l) = w(:, j0 + l) + f * x(:, jr + op%pi(l, g))
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine rand_r

  !> w_j = z_j (A x)_j + (B x)_j for the species parameters A, B (amplitude coordinates): X(z) = z A + B in one pass
  subroutine rand_lin_amp(op, pa, pb, zs, x, w)
    type(augcls_op), intent(in) :: op
    real(dp), intent(in) :: pa(:, :), pb(:, :)
    complex(dp), intent(in) :: zs(:)
    complex(dp), intent(in), contiguous :: x(:, :)
    complex(dp), intent(out), contiguous :: w(:, :)
    real(dp) :: wa(9, op%nsmax, op%nsmax, op%ncls), wb(9, op%nsmax, op%nsmax, op%ncls), fa, fb
    integer :: i, k, kap, l0, l, r, g, j0, jr, s
    if (uniform(pa) .and. uniform(pb)) then
      !$omp parallel do private(j0, l)
      do i = 1, op%nn
        j0 = 9 * (i - 1)
        do l = 1, 9; w(:, j0 + l) = (zs * pa(l, 1) + pb(l, 1)) * x(:, j0 + l); end do
      end do
      !$omp end parallel do
      return
    end if
    call weights(op, pa, wa); call weights(op, pb, wb)
    !$omp parallel do private(k, kap, l0, l, r, g, j0, jr, fa, fb, s)
    do i = 1, op%nn
      j0 = 9 * (i - 1)
      if (op%frozen(i) > 0) then
        s = op%frozen(i)
        do l = 1, 9; w(:, j0 + l) = (zs * pa(l, s) + pb(l, s)) * x(:, j0 + l); end do
        cycle
      end if
      kap = op%ncl(i); l0 = op%lst(i) + 1
      do l = 1, 9; w(:, j0 + l) = (zs * wa(l, l0, l0, kap) + wb(l, l0, l0, kap)) * x(:, j0 + l); end do
      do k = 1, op%nsc(kap)
        if (k == l0 .or. op%flipr(k, i) == 0) cycle
        r = op%flipr(k, i); g = op%flipg(k, i); jr = 9 * (r - 1)
        do l = 1, 9
          fa = wa(l, l0, k, kap) * op%chi(g) * op%sg(l, g); fb = wb(l, l0, k, kap) * op%chi(g) * op%sg(l, g)
          w(:, j0 + l) = w(:, j0 + l) + (zs * fa + fb) * x(:, jr + op%pi(l, g))
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine rand_lin_amp

  !> t2 = S t1 (amplitude coordinates): on-site block and one 9 x 9 bond block per neighbour slot.  The accumulator is (orbital,
  !> member), so the innermost loop runs over the 9 orbitals of one member, as in asr_augcls; the sums over (slot, l2) keep the
  !> order of asr_augcls.
  subroutine smul(op, t1, t2)
    type(augcls_op), intent(in) :: op
    complex(dp), intent(in), contiguous :: t1(:, :)
    complex(dp), intent(out), contiguous :: t2(:, :)
    complex(dp) :: acc(9, size(t1, 1)), s
    real(dp) :: fc(9)
    integer :: i, k, l, l2, g, r, j0, jr, j, nb
    nb = size(t1, 1)
    !$omp parallel do private(k, l, l2, g, r, j0, jr, j, s, fc, acc)
    do i = 1, op%nn
      j0 = 9 * (i - 1)
      acc = (0.0_dp, 0.0_dp)
      do j = 1, nb
        do l2 = 1, 9
          s = t1(j, j0 + l2)
          acc(:, j) = acc(:, j) + op%block0(:, l2) * s
        end do
      end do
      do k = 1, op%nslot
        r = op%hopr(k, i)
        if (r == 0) cycle
        g = op%hopg(k, i); jr = 9 * (r - 1)
        do l2 = 1, 9
          fc = op%chi(g) * op%btab(:, l2, k, g)
          do j = 1, nb
            acc(:, j) = acc(:, j) + fc * t1(j, jr + l2)
          end do
        end do
      end do
      do l = 1, 9; t2(:, j0 + l) = acc(l, :); end do
    end do
    !$omp end parallel do
  end subroutine smul

  !> hx = c x + D S D x (amplitude coordinates); t1, t2 are temporaries
  subroutine hmul(op, x, hx, t1, t2)
    type(augcls_op), intent(in) :: op
    complex(dp), intent(in), contiguous :: x(:, :)
    complex(dp), intent(out), contiguous :: hx(:, :), t1(:, :), t2(:, :)
    call rand_r(op, op%d, x, t1)
    call smul(op, t1, t2)
    call rand_r(op, op%d, t2, t1)
    call rand_r(op, op%c, x, hx)
    hx = hx + t1
  end subroutine hmul

  !> edc step: y_j = H x_j + X(z_j) x_j with H = h - h o h (+ e if present) and X(z) = z pa + pb (normalised coordinates).  This is
  !> augcls_apply on the real and the imaginary part plus augcls_rand_c with the energy-dependent diagonal of asr_edc.
  subroutine blk_hamx(op, wk, pa, pb, zs, x, y)
    type(augcls_op), intent(in) :: op
    type(blk_work_t), intent(inout) :: wk
    real(dp), intent(in) :: pa(:, :), pb(:, :)
    complex(dp), intent(in) :: zs(:)
    complex(dp), intent(in), contiguous :: x(:, :)
    complex(dp), intent(out), contiguous :: y(:, :)
    integer :: nb, i, j0, l
    nb = size(x, 1)
    if (nb /= wk%nb) error stop 'asr_augblk: the block must have wk%nb members'
    associate (x0 => wk%x0, h1 => wk%h1, r1 => wk%r1, h2 => wk%h2, e1 => wk%e1, &
               t1 => wk%t1, t2 => wk%t2)
      call to_amp(op, x, x0)
      call hmul(op, x0, h1, t1, t2)
      if (any(op%o /= 0.0_dp)) then
        call rand_r(op, op%o, h1, r1)
        call hmul(op, r1, h2, t1, t2)
        h1 = h1 - h2
      end if
      call rand_lin_amp(op, pa, pb, zs, x0, e1)
      if (any(op%e /= 0.0_dp)) then
        call rand_r(op, op%e, x0, r1)
        e1 = e1 + r1
      end if
      !$omp parallel do private(j0, l)
      do i = 1, op%nn
        j0 = 9 * (i - 1)
        do l = 1, 9; y(:, j0 + l) = op%sq(i) * (h1(:, j0 + l) + e1(:, j0 + l)); end do
      end do
      !$omp end parallel do
    end associate
  end subroutine blk_hamx

  !> ed3 step: y_j = K'(z_j) x_j = et x - h'x + u h'x + h'(u x - Qt h'x) (as asr_ed3::ed3_kmul), with every coefficient linear in z:
  !> et = z ea + eb, u = z ua + ub, Qt = z qa + qb (normalised coordinates; op carries h' as in ed3_setup)
  subroutine blk_kmul3(op, wk, ea, eb, ua, ub, qa, qb, zs, x, y)
    type(augcls_op), intent(in) :: op
    type(blk_work_t), intent(inout) :: wk
    real(dp), intent(in) :: ea(:, :), eb(:, :), ua(:, :), ub(:, :), qa(:, :), qb(:, :)
    complex(dp), intent(in) :: zs(:)
    complex(dp), intent(in), contiguous :: x(:, :)
    complex(dp), intent(out), contiguous :: y(:, :)
    integer :: nb, i, j0, l
    nb = size(x, 1)
    if (nb /= wk%nb) error stop 'asr_augblk: the block must have wk%nb members'
    associate (x0 => wk%x0, h1 => wk%h1, r1 => wk%r1, h2 => wk%h2, e1 => wk%e1, &
               t1 => wk%t1, t2 => wk%t2)
      call to_amp(op, x, x0)
      call hmul(op, x0, h1, t1, t2)                     ! h'x
      call rand_lin_amp(op, qa, qb, zs, h1, r1)         ! Qt h'x
      call rand_lin_amp(op, ua, ub, zs, x0, e1)         ! u x
      e1 = e1 - r1
      call hmul(op, e1, h2, t1, t2)                     ! h'(u x - Qt h'x)
      call rand_lin_amp(op, ea, eb, zs, x0, r1)         ! et x
      call rand_lin_amp(op, ua, ub, zs, h1, e1)         ! u h'x
      !$omp parallel do private(j0, l)
      do i = 1, op%nn
        j0 = 9 * (i - 1)
        do l = 1, 9; y(:, j0 + l) = op%sq(i) * (r1(:, j0 + l) - h1(:, j0 + l) + e1(:, j0 + l) + h2(:, j0 + l)); end do
      end do
      !$omp end parallel do
    end associate
  end subroutine blk_kmul3

  !> projection of every member of x onto the symmetric subspace of the reduced operator, in place (as augcls_project_c)
  subroutine blk_project(op, x)
    type(augcls_op), intent(in) :: op
    complex(dp), intent(inout), contiguous :: x(:, :)
    complex(dp) :: acc(size(x, 1), 9)
    integer :: i, j, l, g, ns, j0
    !$omp parallel do private(j, l, g, ns, j0, acc)
    do i = 1, op%nn
      ns = op%stab_ptr(i + 1) - op%stab_ptr(i)
      if (ns <= 1) cycle
      j0 = 9 * (i - 1)
      acc = (0.0_dp, 0.0_dp)
      do j = op%stab_ptr(i), op%stab_ptr(i + 1) - 1
        g = op%stab(j)
        do l = 1, 9
          acc(:, l) = acc(:, l) + op%chi(g) * op%sg(l, g) * x(:, j0 + op%pi(l, g))
        end do
      end do
      do l = 1, 9; x(:, j0 + l) = acc(:, l) / real(ns, dp); end do
    end do
    !$omp end parallel do
  end subroutine blk_project

  !> y = h x = (c + D S D) x, normalised coordinates (as augcls_hp_c)
  subroutine blk_hp(op, wk, x, y)
    type(augcls_op), intent(in) :: op
    type(blk_work_t), intent(inout) :: wk
    complex(dp), intent(in), contiguous :: x(:, :)
    complex(dp), intent(out), contiguous :: y(:, :)
    integer :: nb, i, j0, l
    nb = size(x, 1)
    if (nb /= wk%nb) error stop 'asr_augblk: the block must have wk%nb members'
    associate (x0 => wk%x0, h1 => wk%h1, t1 => wk%t1, t2 => wk%t2)
      call to_amp(op, x, x0)
      call hmul(op, x0, h1, t1, t2)
      do i = 1, op%nn
        j0 = 9 * (i - 1)
        do l = 1, 9; y(:, j0 + l) = op%sq(i) * h1(:, j0 + l); end do
      end do
    end associate
  end subroutine blk_hp

  !> w_j = X(z_j) x_j with X(z) = z pa + pb, normalised coordinates (as augcls_rand_c)
  subroutine blk_rand_lin(op, wk, pa, pb, zs, x, w)
    type(augcls_op), intent(in) :: op
    type(blk_work_t), intent(inout) :: wk
    real(dp), intent(in) :: pa(:, :), pb(:, :)
    complex(dp), intent(in) :: zs(:)
    complex(dp), intent(in), contiguous :: x(:, :)
    complex(dp), intent(out), contiguous :: w(:, :)
    integer :: nb, i, j0, l
    nb = size(x, 1)
    if (nb /= wk%nb) error stop 'asr_augblk: the block must have wk%nb members'
    associate (x0 => wk%x0, e1 => wk%e1)
      call to_amp(op, x, x0)
      call rand_lin_amp(op, pa, pb, zs, x0, e1)
      do i = 1, op%nn
        j0 = 9 * (i - 1)
        do l = 1, 9; w(:, j0 + l) = op%sq(i) * e1(:, j0 + l); end do
      end do
    end associate
  end subroutine blk_rand_lin
end module asr_augblk
