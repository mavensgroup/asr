!> One random realisation of a binary LMTO alloy on a periodic supercell, H = E_nu + h - h o h applied
!> (order 3: + h o h o h - (h p h h + h h p h)/2, the third-order orthonormal ASA Hamiltonian, p = p^gamma)
!> matrix-free (same structure as asr_augorb, but with numerical site occupations).  Ground truth for
!> the augmented-space recursion via KPM / exact diagonalisation of supercells.
module asr_lmto_realization
  use asr_kinds, only: dp
  use asr_operator, only: linop
  use asr_augspace, only: graph_t
  use asr_augorb, only: orb_alloy_t
  implicit none
  private
  public :: lmto_real_op, build_lmto_real, build_lmto_ordered, lmto_s_apply, lmto_hph

  type, extends(linop) :: lmto_real_op
    integer :: n = 0, nl = 9, nslot = 0, order = 2
    logical, allocatable :: isa(:)
    real(dp), allocatable :: c(:, :), d(:, :), o(:, :), e(:, :), p(:, :)     !! (nl, site)
    real(dp), allocatable :: block(:, :, :), block0(:, :)
    integer, allocatable :: nbr(:, :)                               !! (slot, site)
  contains
    procedure :: apply => real_apply
    procedure :: dim => real_dim
  end type lmto_real_op

contains

  pure integer function real_dim(self)
    class(lmto_real_op), intent(in) :: self
    real_dim = self%nl * self%n
  end function real_dim

  subroutine build_lmto_real(g, mdl, seed, op)
    type(graph_t), intent(in) :: g
    type(orb_alloy_t), intent(in) :: mdl
    integer, intent(in) :: seed
    type(lmto_real_op), intent(out) :: op
    real(dp), allocatable :: r(:)
    integer, allocatable :: sd(:)
    logical, allocatable :: isa(:)
    integer :: i, k
    allocate (r(g%n), isa(g%n))
    call random_seed(size=k); allocate (sd(k)); sd = seed + 37 * [(i, i=1, k)]; call random_seed(put=sd)
    call random_number(r); isa = r < mdl%x
    call build_lmto_ordered(g, mdl, isa, op)
  end subroutine build_lmto_real

  !> Operator for a GIVEN occupation pattern isa(site) (.true. = species A parameters, .false. = B): ordered
  !> crystals (B2: A on the even sublattice) and any fixed supercell configuration.  mdl%x is not used.
  subroutine build_lmto_ordered(g, mdl, isa, op)
    type(graph_t), intent(in) :: g
    type(orb_alloy_t), intent(in) :: mdl
    logical, intent(in) :: isa(:)
    type(lmto_real_op), intent(out) :: op
    integer :: i, ns, l
    op%n = g%n; op%nl = mdl%nl
    ns = g%nbr_ptr(2) - g%nbr_ptr(1); op%nslot = ns
    allocate (op%isa(g%n), op%nbr(ns, g%n))
    op%isa = isa
    do i = 1, g%n; op%nbr(:, i) = g%nbr(g%nbr_ptr(i):g%nbr_ptr(i + 1) - 1); end do
    allocate (op%c(mdl%nl, g%n), op%d(mdl%nl, g%n), op%o(mdl%nl, g%n), op%e(mdl%nl, g%n), op%p(mdl%nl, g%n))
    op%p = 0.0_dp
    do i = 1, g%n
      do l = 1, mdl%nl
        if (op%isa(i)) then
          op%c(l, i) = mdl%c0(l) + mdl%c1(l); op%d(l, i) = mdl%d0(l) + mdl%d1(l)
          op%o(l, i) = mdl%o0(l) + mdl%o1(l); op%e(l, i) = mdl%e0(l) + mdl%e1(l)
          if (allocated(mdl%p0)) op%p(l, i) = mdl%p0(l) + mdl%p1(l)
        else
          op%c(l, i) = mdl%c0(l); op%d(l, i) = mdl%d0(l); op%o(l, i) = mdl%o0(l); op%e(l, i) = mdl%e0(l)
          if (allocated(mdl%p0)) op%p(l, i) = mdl%p0(l)
        end if
      end do
    end do
    op%block = mdl%block
    allocate (op%block0(mdl%nl, mdl%nl)); op%block0 = 0.0_dp
    if (allocated(mdl%block0)) op%block0 = mdl%block0
  end subroutine build_lmto_ordered

  !> sx = S x (screened structure-constant hop incl. the on-site block, no D scaling) on a complex vector
  subroutine lmto_s_apply(self, x, sx)
    class(lmto_real_op), intent(in) :: self
    complex(dp), intent(in) :: x(:)
    complex(dp), intent(out) :: sx(:)
    integer :: i, k, j, l, l2, nl
    complex(dp) :: acc
    nl = self%nl
    !$omp parallel do private(l, l2, k, j, acc)
    do i = 1, self%n
      do l = 1, nl
        acc = (0.0_dp, 0.0_dp)
        do l2 = 1, nl; acc = acc + self%block0(l, l2) * x(nl * (i - 1) + l2); end do
        do k = 1, self%nslot
          j = self%nbr(k, i)
          do l2 = 1, nl; acc = acc + self%block(l, l2, k) * x(nl * (j - 1) + l2); end do
        end do
        sx(nl * (i - 1) + l) = acc
      end do
    end do
    !$omp end parallel do
  end subroutine lmto_s_apply

  !> y = h p h x (p = p^gamma): overlap term, for the sphere-charge start vector (1 - h p h / 2)|L>
  subroutine lmto_hph(self, x, y)
    class(lmto_real_op), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    real(dp), allocatable :: t(:)
    integer :: i, l, nl
    nl = self%nl
    allocate (t(size(x)))
    call h_apply(self, x, t)
    !$omp parallel do private(l)
    do i = 1, self%n
      do l = 1, nl; t(nl * (i - 1) + l) = self%p(l, i) * t(nl * (i - 1) + l); end do
    end do
    !$omp end parallel do
    call h_apply(self, t, y)
  end subroutine lmto_hph

  subroutine h_apply(self, x, hx)
    class(lmto_real_op), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: hx(:)
    real(dp), allocatable :: t(:)
    integer :: i, k, j, l, l2, nl
    real(dp) :: acc
    nl = self%nl
    allocate (t(size(x)))
    !$omp parallel do private(l)
    do i = 1, self%n
      do l = 1, nl; t(nl * (i - 1) + l) = self%d(l, i) * x(nl * (i - 1) + l); end do
    end do
    !$omp end parallel do
    !$omp parallel do private(l, l2, k, j, acc)
    do i = 1, self%n
      do l = 1, nl
        acc = 0.0_dp
        do l2 = 1, nl; acc = acc + self%block0(l, l2) * t(nl * (i - 1) + l2); end do
        do k = 1, self%nslot
          j = self%nbr(k, i)
          do l2 = 1, nl; acc = acc + self%block(l, l2, k) * t(nl * (j - 1) + l2); end do
        end do
        hx(nl * (i - 1) + l) = self%c(l, i) * x(nl * (i - 1) + l) + self%d(l, i) * acc
      end do
    end do
    !$omp end parallel do
  end subroutine h_apply

  subroutine real_apply(self, x, y)
    class(lmto_real_op), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    real(dp), allocatable :: hx(:), w(:), hw(:), b(:), r(:), t(:)
    integer :: i, l, nl, j
    nl = self%nl
    allocate (hx(size(x)), w(size(x)), hw(size(x)))
    call h_apply(self, x, hx)
    !$omp parallel do private(l)
    do i = 1, self%n
      do l = 1, nl; w(nl * (i - 1) + l) = self%o(l, i) * hx(nl * (i - 1) + l); end do
    end do
    !$omp end parallel do
    call h_apply(self, w, hw)
    if (self%order == 3) then              ! + h [ o (h o h x) - (p h h x + h p h x)/2 ]  (two applications of h folded into one)
      allocate (b(size(x)), r(size(x)), t(size(x)))
      call h_apply(self, hx, b)                                      ! b = h h x
      !$omp parallel do private(l, j)
      do i = 1, self%n
        do l = 1, nl; j = nl * (i - 1) + l; w(j) = self%p(l, i) * hx(j); end do
      end do
      !$omp end parallel do
      call h_apply(self, w, r)                                       ! r = h p h x
      !$omp parallel do private(l, j)
      do i = 1, self%n
        do l = 1, nl
          j = nl * (i - 1) + l
          t(j) = self%o(l, i) * hw(j) - 0.5_dp * (self%p(l, i) * b(j) + r(j))
        end do
      end do
      !$omp end parallel do
      call h_apply(self, t, b)                                       ! b = h t
      hw = hw - b                                                    ! y = E x + h x - (hw - h t)
    end if
    !$omp parallel do private(l)
    do i = 1, self%n
      do l = 1, nl
        y(nl * (i - 1) + l) = self%e(l, i) * x(nl * (i - 1) + l) + hx(nl * (i - 1) + l) - hw(nl * (i - 1) + l)
      end do
    end do
    !$omp end parallel do
  end subroutine real_apply
end module asr_lmto_realization
