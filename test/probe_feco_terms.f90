!> probe_feco_terms: which term of the orthonormal ASA Hamiltonian controls the B2 FeCo moments?
!> At FIXED self-consistent k-space ASA parameters (standard TB-LMTO screening, as in the recursion route), with
!> K = h (1 + o h)^-1 and the overlap 1 + K p K (p = p^gamma):
!>   full : E_nu + (1 + KpK)^-1/2 K (1 + KpK)^-1/2          (exact ASA; weights = orthonormal-orbital weights)
!>   nop  : E_nu + K                                        (all orders in o, no p: the energy-dependent route)
!>   h2   : E_nu + h - h o h                                (the energy-independent recursion route)
!>   h3o  : h2 + h o h o h
!>   h3   : h3o - (h p h h + h h p h)/2                     (complete third order)
!> One-shot site moments (common E_F, linear tetrahedra) for each, no self-consistency.
!> usage: probe_feco_terms [nk_scf] [nk] [feco|cuzn|fefe]
program probe_feco_terms
  use asr_kinds, only: dp
  use asr_asa
  use asr_atom, only: shell_t
  use asr_lmto_cell
  use asr_lmto_ham, only: potpar_t
  use asr_strconst, only: standard_alpha
  use asr_tetra, only: tetra_weights, tetra_fermi
  implicit none
  type(asa_input_t) :: oi
  type(asa_result_t) :: orr
  type(cell_t) :: cell
  type(cluster_t), allocatable :: cls(:)
  type(potpar_t) :: pps(2, 2)
  character(len=16) :: arg, sys
  real(dp) :: alat_a, nval(2), zat(2), occd(2), occs(2)
  integer :: nsp_run
  character(len=4), parameter :: names(5) = ['full', 'nop ', 'h2  ', 'h3o ', 'h3  ']
  real(dp), allocatable :: alpha_sp(:, :), eb(:, :, :), wt(:, :, :, :), eall(:, :), wall(:, :), wl(:, :, :, :, :), ws(:, :, :, :, :)
  real(dp) :: ql(0:2, 2, 2), qs(0:2, 2, 2)
  complex(dp) :: zv(18), kz(18), om12(18, 18), kmt(18, 18)
  integer :: a, ll
  integer, parameter :: lof(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]
  real(dp) :: bv(3, 3), k(3), ef, m(2), qq(2), nel, ov(18), pv(18), en(18), w(18), rwork(54)
  complex(dp) :: hh(18, 18), kk(18, 18), t(18, 18), hm(18, 18), cw(72)
  integer :: nks, nk, ic, s, l, ix, iy, iz, ik, iv, i, j, info, nst
  external :: zheev

  nks = 16; nk = 24
  if (command_argument_count() >= 1) then; call get_command_argument(1, arg); read (arg, *) nks; end if
  if (command_argument_count() >= 2) then; call get_command_argument(2, arg); read (arg, *) nk; end if
  sys = 'feco'
  if (command_argument_count() >= 3) call get_command_argument(3, sys)
  select case (trim(sys))
  case ('feco'); zat = [26, 27]; nval = [8, 9];   occs = [2, 2]; occd = [6, 7];   alat_a = 2.857_dp;   nsp_run = 2
  case ('cuzn'); zat = [29, 30]; nval = [11, 12]; occs = [1, 2]; occd = [10, 10]; alat_a = 2.92055_dp; nsp_run = 1
  case ('fefe'); zat = [26, 26]; nval = [8, 8];   occs = [2, 2]; occd = [6, 6];   alat_a = 2.8665_dp;  nsp_run = 2
  case default; error stop 'system must be feco, cuzn or fefe'
  end select
  allocate (oi%cls(2))
  do ic = 1, 2; oi%cls(ic)%z = zat(ic); oi%cls(ic)%nval = nint(nval(ic)); end do
  do ic = 1, 2
    oi%cls(ic)%core = core_shells_for([4, 4, 3, 4]); oi%cls(ic)%nvl = [4, 4, 3]
  end do
  do ic = 1, 2; oi%cls(ic)%atom_shells = [oi%cls(ic)%core, shell_t(4, 0, occs(ic)), shell_t(3, 2, occd(ic))]; end do
  oi%alat = alat_a / 0.529177210903_dp; oi%nk = nks; oi%maxit = 200; oi%nspin = nsp_run; oi%zeta0 = 0.1_dp
  oi%beta = 0.2_dp; oi%tol = 1.0e-7_dp
  cell%nat = 2
  cell%avec = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
  allocate (cell%tau(3, 2), cell%isp(2)); cell%tau = 0.0_dp; cell%tau(:, 2) = 0.5_dp; cell%isp = [1, 2]
  call setup_cell(cell)
  do ic = 1, 2; oi%cls(ic)%sws = cell%w * oi%alat; end do
  call asa_scf(oi, cell, orr)
  print '(a,l1,a,f10.6,a,2f9.5)', 'k-space ASA SCF: converged=', orr%converged, '  E_F=', orr%ef, '  moments Fe Co =', &
      orr%cls(1)%moment, orr%cls(2)%moment
  if (nsp_run == 1) print '(a,2f9.5)', 'k-space ASA qnet =', orr%cls(1)%qnet, orr%cls(2)%qnet

  allocate (alpha_sp(0:2, 2), cls(2))
  do ic = 1, 2
    alpha_sp(:, ic) = standard_alpha()
    do s = 1, 2
      do l = 0, 2
        if (nsp_run == 1) orr%cls(ic)%pp(l, 2) = orr%cls(ic)%pp(l, 1)
        pps(ic, s)%enu(l) = orr%cls(ic)%pp(l, s)%enu; pps(ic, s)%c(l) = orr%cls(ic)%pp(l, s)%c
        pps(ic, s)%sqrdel(l) = orr%cls(ic)%pp(l, s)%sqrdel; pps(ic, s)%p(l) = orr%cls(ic)%pp(l, s)%p
        pps(ic, s)%gamma(l) = orr%cls(ic)%pp(l, s)%gamma; pps(ic, s)%alpha_rep(l) = orr%cls(ic)%pp(l, s)%gamma
      end do
    end do
  end do
  do ic = 1, 2; call build_cluster(cell, ic, 4.0_dp * cell%w, alpha_sp, cls(ic)); end do
  bv = recip_vectors(cell%avec)
  nst = nk**3
  allocate (eb(18, nst, 2), wt(18, 2, nst, 2), eall(36, nst), wall(36, nst), wl(18, 0:2, 2, nst, 2), ws(18, 0:2, 2, nst, 2))
  nel = sum(nval)
  do iv = 1, 5
    do s = 1, 2
      !$omp parallel do collapse(3) private(ix, iy, iz, ik, k, hh, ov, pv, en, kk, t, hm, w, cw, rwork, info, i, j, zv, kz, om12, kmt, a, ll)
      do iz = 0, nk - 1
        do iy = 0, nk - 1
          do ix = 0, nk - 1
            ik = 1 + ix + nk * (iy + nk * iz)
            k = (ix + 0.5_dp) / nk * bv(:, 1) + (iy + 0.5_dp) / nk * bv(:, 2) + (iz + 0.5_dp) / nk * bv(:, 3)
            call hsmall_k(cell, cls, pps(:, s), alpha_sp, k, hh, ov, pv, en)
            pv = pv - ov**2                                   ! p^gamma (hsmall_k returns p^alpha = p^gamma + o^2)
            call variant(iv, hh, ov, pv, en, hm, om12, kmt)
            call zheev('V', 'U', 18, hm, 18, w, cw, 72, rwork, info)
            eb(:, ik, s) = w
            do j = 1, 18
              wt(j, 1, ik, s) = sum(abs(hm(1:9, j))**2); wt(j, 2, ik, s) = sum(abs(hm(10:18, j))**2)
              wl(j, :, :, ik, s) = 0.0_dp; ws(j, :, :, ik, s) = 0.0_dp
              if (iv == 5) then                                  ! approximate: z = (1 - h p h / 2) y, |Kz|^2 ~ (e - E_nu)^2 |z|^2
                zv = matmul(om12, hm(:, j))
                zv = hm(:, j) - 0.5_dp * zv + 0.375_dp * matmul(om12, zv)     ! (1 - X/2 + 3X^2/8) y, X = h p h
                do i = 1, 18
                  a = (i - 1) / 9 + 1; ll = lof(mod(i - 1, 9) + 1)
                  wl(j, ll, a, ik, s) = wl(j, ll, a, ik, s) + abs(hm(i, j))**2
                  ws(j, ll, a, ik, s) = ws(j, ll, a, ik, s) + abs(zv(i))**2 * (1.0_dp + pv(i) * (w(j) - en(i))**2)
                end do
              end if
              if (iv == 1) then                                  ! z = O^-1/2 y, sphere charge |z|^2 + p |K z|^2
                zv = matmul(om12, hm(:, j)); kz = matmul(kmt, zv)
                do i = 1, 18
                  a = (i - 1) / 9 + 1; ll = lof(mod(i - 1, 9) + 1)
                  wl(j, ll, a, ik, s) = wl(j, ll, a, ik, s) + abs(hm(i, j))**2
                  ws(j, ll, a, ik, s) = ws(j, ll, a, ik, s) + abs(zv(i))**2 + pv(i) * abs(kz(i))**2
                end do
              end if
            end do
          end do
        end do
      end do
      !$omp end parallel do
      eall((s - 1) * 18 + 1:s * 18, :) = eb(:, :, s)
    end do
    call tetra_fermi(nk, eall, nel, 1.0_dp, ef)
    call tetra_weights(nk, eall, ef, wall)
    m = 0.0_dp; qq = 0.0_dp
    do s = 1, 2
      do ik = 1, nst
        do j = 1, 18
          m = m + merge(1.0_dp, -1.0_dp, s == 1) * wall((s - 1) * 18 + j, ik) * wt(j, :, ik, s)
          qq = qq + wall((s - 1) * 18 + j, ik) * wt(j, :, ik, s)
        end do
      end do
    end do
    print '(a,a4,a,f10.6,a,2f9.5,a,f9.5,a,2f9.4,a,2f9.5)', 'RESULT ', names(iv), '  E_F=', ef, '  m(Fe,Co)=', m, '  cell=', sum(m), &
        '  eig min/max=', minval(eb), maxval(eb), '  qnet=', qq - nval
    if (iv == 1 .or. iv == 5) then
      ql = 0.0_dp; qs = 0.0_dp
      do s = 1, 2
        do ik = 1, nst
          do j = 1, 18
            ql(:, :, s) = ql(:, :, s) + wall((s - 1) * 18 + j, ik) * wl(j, :, :, ik, s)
            qs(:, :, s) = qs(:, :, s) + wall((s - 1) * 18 + j, ik) * ws(j, :, :, ik, s)
          end do
        end do
      end do
      do a = 1, 2
        do s = 1, 2
          print '(a,i1,a,i1,a,3f9.5,a,3f9.5)', 'CHARGE site ', a, ' spin ', s, '  Lowdin s p d', ql(:, a, s), '   sphere s p d', qs(:, a, s)
        end do
      end do
      print '(a,2f9.5,a,2f9.5)', 'CHARGE moments Lowdin', sum(ql(:, :, 1), 1) - sum(ql(:, :, 2), 1), '   sphere', sum(qs(:, :, 1), 1) - sum(qs(:, :, 2), 1)
    end if
  end do
contains
  !> orthonormal Hamiltonian of variant iv (18x18) from h, o, p^gamma, E_nu
  subroutine variant(iv, h, o, p, e, hout, o12, kout)
    integer, intent(in) :: iv
    complex(dp), intent(in) :: h(18, 18)
    real(dp), intent(in) :: o(18), p(18), e(18)
    complex(dp), intent(out) :: hout(18, 18), o12(18, 18), kout(18, 18)
    complex(dp) :: a(18, 18), kmat(18, 18), x(18, 18), u(18, 18), hoh(18, 18), hp(18, 18), hpp(18, 18)
    real(dp) :: lam(18), rw(54)
    complex(dp) :: work(72)
    integer :: ii, jj, ipiv(18), inf
    external :: zgesv, zheev
    do jj = 1, 18; hoh(:, jj) = o * h(:, jj); end do          ! hoh = o h  (row-scaled)
    hoh = matmul(h, hoh)                                      ! h o h
    select case (iv)
    case (1, 2)
      a = hoh * 0.0_dp
      do jj = 1, 18; a(:, jj) = o * h(:, jj); a(jj, jj) = a(jj, jj) + 1.0_dp; end do   ! 1 + o h
      ! K = h (1 + o h)^-1  <=>  K^T solves (1 + o h)^T K^T = h^T
      x = transpose(h); u = transpose(a)
      call zgesv(18, 18, u, 18, ipiv, x, 18, inf)
      kmat = transpose(x)
      kmat = 0.5_dp * (kmat + conjg(transpose(kmat)))
      if (iv == 2) then
        hout = kmat
      else
        do jj = 1, 18; hp(:, jj) = p * kmat(:, jj); end do
        x = matmul(kmat, hp)                                    ! K p K
        do jj = 1, 18; x(jj, jj) = x(jj, jj) + 1.0_dp; end do
        u = x
        call zheev('V', 'U', 18, u, 18, lam, work, 72, rw, inf)
        x = 0.0_dp
        do jj = 1, 18; x(:, jj) = u(:, jj) / sqrt(lam(jj)); end do
        x = matmul(x, conjg(transpose(u)))                      ! (1 + KpK)^-1/2
        hout = matmul(x, matmul(kmat, x))
        o12 = x; kout = kmat
      end if
    case (3)
      hout = h - hoh
    case (4, 5)
      do jj = 1, 18; x(:, jj) = o * hoh(:, jj); end do          ! o h o h
      hout = h - hoh + matmul(h, x)
      if (iv == 5) then
        do jj = 1, 18; hp(:, jj) = p * h(:, jj); end do
        hpp = matmul(h, hp)                                     ! h p h
        hout = hout - 0.5_dp * (matmul(hpp, h) + matmul(h, hpp))
        o12 = hpp
      end if
    end select
    hout = 0.5_dp * (hout + conjg(transpose(hout)))
    do ii = 1, 18; hout(ii, ii) = hout(ii, ii) + e(ii); end do
  end subroutine variant
end program probe_feco_terms
