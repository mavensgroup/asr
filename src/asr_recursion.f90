!> Lanczos/Haydock recursion with normalised vectors.
!>   H|u_n> = b_{n-1}|u_{n-1}> + a_n|u_n> + b_n|u_{n+1}>
!> so Haydock's beta_n^2 of the unnormalised form equals b_n**2 here.
module asr_recursion
  use asr_kinds, only: dp
  use asr_operator, only: linop
  implicit none
  private
  public :: lanczos

contains

  !> Run up to nmax steps from v0 (normalised internally).
  !> a(1:n), b(1:n-1) are filled; n returned (< nmax if the Krylov space is
  !> exhausted, i.e. b_n < tiny).  reorth=.true. does full reorthogonalisation
  !> (O(n) stored vectors; for tests and small cases).
  subroutine lanczos(op, v0, nmax, a, b, n, reorth, tiny)
    class(linop), intent(in) :: op
    real(dp), intent(in) :: v0(:)
    integer, intent(in) :: nmax
    real(dp), intent(out) :: a(nmax), b(nmax)
    integer, intent(out) :: n
    logical, intent(in), optional :: reorth
    real(dp), intent(in), optional :: tiny
    real(dp), allocatable :: u(:), up(:), w(:), q(:, :)
    real(dp) :: bn, tol
    logical :: ro
    integer :: i, dim

    ro = .false.; if (present(reorth)) ro = reorth
    tol = 1.0e-12_dp; if (present(tiny)) tol = tiny
    dim = op%dim()
    allocate (u(dim), up(dim), w(dim))
    if (ro) allocate (q(dim, nmax))
    a = 0.0_dp; b = 0.0_dp
    u = v0 / sqrt(dot_product(v0, v0))
    up = 0.0_dp
    bn = 0.0_dp
    n = 0
    do i = 1, nmax
      if (ro) q(:, i) = u
      call op%apply(u, w)
      a(i) = dot_product(u, w)
      w = w - a(i) * u - bn * up
      if (ro) then
        w = w - matmul(q(:, 1:i), matmul(transpose(q(:, 1:i)), w))
        w = w - matmul(q(:, 1:i), matmul(transpose(q(:, 1:i)), w))
      end if
      n = i
      bn = sqrt(dot_product(w, w))
      if (i == nmax) exit
      b(i) = bn
      if (bn < tol) exit
      up = u
      u = w / bn
    end do
    if (n == nmax) b(nmax) = bn   ! b_n of the last step (needed by terminators)
  end subroutine lanczos
end module asr_recursion
