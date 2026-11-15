!> probe_feco_depth: recursion depth convergence of B2 FeCo at FIXED parameters (no SCF feedback).
!> The self-consistent k-space ASA parameters are fed to the ordered real-space operator (standard screening,
!> third-order H, sphere-charge start vectors (1 - hph/2)|L>); E_F from the electron count (m0 + p m2) and the
!> site moments are compared with the k-space H^(3) value at the same parameters (probe_feco_terms).
!> usage: probe_feco_depth [supercell] [depth ...]
program probe_feco_depth
  use asr_kinds, only: dp, pi
  use asr_asa
  use asr_atom, only: shell_t
  use asr_lmto_cell
  use asr_lmto_ham, only: potpar_t
  use asr_strconst, only: standard_alpha
  use asr_augspace, only: graph_t
  use asr_augorb, only: orb_alloy_t
  use asr_lmto_alloy
  use asr_lmto_realization
  use asr_recursion, only: lanczos
  use asr_greens, only: cf_green, sqrt_tail, tail_params_avg2
  use asr_contour, only: contour_nodes
  implicit none
  type(asa_input_t) :: oi
  type(asa_result_t) :: orr
  type(cell_t) :: cell
  type(potpar_t) :: pps(2, 2)
  type(graph_t) :: g
  type(orb_alloy_t) :: mdl
  type(lmto_real_op) :: rop
  integer, allocatable :: coord(:, :)
  logical, allocatable :: isa(:)
  real(dp), allocatable :: v0(:), hv(:), ta(:), tb(:), ca(:, :, :, :, :)
  integer, parameter :: iorb(4) = [1, 3, 5, 7], lorb(4) = [0, 1, 2, 2], ncont = 80
  real(dp), parameter :: worb(4) = [1.0_dp, 3.0_dp, 3.0_dp, 2.0_dp]
  real(dp) :: nrm(4, 2, 2), lo, hi, ef, nn, mom(2), q(0:2, 0:2, 2, 2), ainf, binf, ebot
  complex(dp) :: zc(ncont), wc(ncont), gg
  integer :: lsup, nlev, nd, ic, s, l, ip, oc, site_b, i, it, k, ia, nargs, m
  integer, allocatable :: mm(:, :, :)
  character(len=16) :: arg

  lsup = 12
  nargs = command_argument_count()
  if (nargs >= 1) then; call get_command_argument(1, arg); read (arg, *) lsup; end if
  allocate (oi%cls(2))
  oi%cls(1)%z = 26.0_dp; oi%cls(1)%nval = 8; oi%cls(2)%z = 27.0_dp; oi%cls(2)%nval = 9
  do ic = 1, 2; oi%cls(ic)%core = core_shells_for([4, 4, 3, 4]); oi%cls(ic)%nvl = [4, 4, 3]; end do
  oi%cls(1)%atom_shells = [oi%cls(1)%core, shell_t(4, 0, 2.0_dp), shell_t(3, 2, 6.0_dp)]
  oi%cls(2)%atom_shells = [oi%cls(2)%core, shell_t(4, 0, 2.0_dp), shell_t(3, 2, 7.0_dp)]
  oi%alat = 2.857_dp / 0.529177210903_dp; oi%nk = 16; oi%maxit = 200; oi%nspin = 2; oi%zeta0 = 0.1_dp
  oi%beta = 0.2_dp; oi%tol = 1.0e-7_dp
  cell%nat = 2
  cell%avec = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
  allocate (cell%tau(3, 2), cell%isp(2)); cell%tau = 0.0_dp; cell%tau(:, 2) = 0.5_dp; cell%isp = [1, 2]
  call setup_cell(cell)
  do ic = 1, 2; oi%cls(ic)%sws = cell%w * oi%alat; end do
  call asa_scf(oi, cell, orr)
  print '(a,f10.6,a,2f9.5)', 'k-space ASA SCF: E_F=', orr%ef, '  moments Fe Co =', orr%cls(1)%moment, orr%cls(2)%moment
  do ic = 1, 2
    do s = 1, 2
      do l = 0, 2
        pps(ic, s)%enu(l) = orr%cls(ic)%pp(l, s)%enu; pps(ic, s)%c(l) = orr%cls(ic)%pp(l, s)%c
        pps(ic, s)%sqrdel(l) = orr%cls(ic)%pp(l, s)%sqrdel; pps(ic, s)%p(l) = orr%cls(ic)%pp(l, s)%p
        pps(ic, s)%gamma(l) = orr%cls(ic)%pp(l, s)%gamma; pps(ic, s)%alpha_rep(l) = orr%cls(ic)%pp(l, s)%gamma
      end do
    end do
  end do
  do ia = 2, max(nargs, 2)
    nlev = 32
    if (nargs >= ia) then; call get_command_argument(ia, arg); read (arg, *) nlev; end if
    allocate (ta(nlev), tb(nlev), ca(2, nlev, 4, 2, 2), mm(4, 2, 2))
    do s = 1, 2
      call build_lmto_alloy('bcc', lsup, 1.2_dp, 0.7_dp, pps(1, s), pps(2, s), 0.5_dp, g, coord, mdl, alpha_in=standard_alpha())
      if (.not. allocated(isa)) then
        allocate (isa(g%n)); isa = modulo(coord(1, :), 2) == 0
        do i = 1, g%n; if (.not. isa(i)) then; site_b = i; exit; end if; end do
        nd = 9 * g%n; allocate (v0(nd), hv(nd))
      end if
      call build_lmto_ordered(g, mdl, isa, rop); rop%order = 3
      do ip = 1, 2
        do oc = 1, 4
          v0 = 0.0_dp; v0(9 * (merge(1, site_b, ip == 1) - 1) + iorb(oc)) = 1.0_dp
          call lmto_hph(rop, v0, hv); v0 = v0 - 0.5_dp * hv
          nrm(oc, s, ip) = dot_product(v0, v0)
          call lanczos(rop, v0, nlev, ta, tb, m)
          mm(oc, s, ip) = m; ca(1, :, oc, s, ip) = ta; ca(2, :, oc, s, ip) = tb
        end do
      end do
    end do
    ebot = minval([(minval(pps(ic, :)%enu(0)), ic=1, 2)]) - 1.2_dp
    lo = -0.5_dp; hi = 0.4_dp
    do it = 1, 50
      ef = 0.5_dp * (lo + hi)
      call occupied(ef)
      nn = 0.0_dp
      do ic = 1, 2; do s = 1, 2; do l = 0, 2
        nn = nn + q(0, l, s, ic) + pps(ic, s)%p(l) * q(2, l, s, ic)
      end do; end do; end do
      if (nn > 17.0_dp) then; hi = ef; else; lo = ef; end if
      if (hi - lo < 1.0e-9_dp) exit
    end do
    do ic = 1, 2
      mom(ic) = 0.0_dp
      do l = 0, 2
        mom(ic) = mom(ic) + q(0, l, 1, ic) + pps(ic, 1)%p(l) * q(2, l, 1, ic) - q(0, l, 2, ic) - pps(ic, 2)%p(l) * q(2, l, 2, ic)
      end do
    end do
    print '(a,i3,a,i3,a,f10.6,a,2f9.5)', 'RESULT supercell', lsup, '  depth', nlev, '  E_F=', ef, '  m(Fe,Co)=', mom
    deallocate (ta, tb, ca, mm)
  end do
contains
  !> sphere-charge moments q(k, l, spin, site) occupied to efv (contour, avg2 square-root tail)
  subroutine occupied(efv)
    real(dp), intent(in) :: efv
    integer :: j, kk
    call contour_nodes(ebot, efv, ncont, zc, wc)
    q = 0.0_dp
    do ic = 1, 2
      do s = 1, 2
        do oc = 1, 4
          call tail_params_avg2(ca(1, :, oc, s, ic), ca(2, :, oc, s, ic), mm(oc, s, ic), ainf, binf)
          do j = 1, ncont
            gg = nrm(oc, s, ic) * cf_green(ca(1, :, oc, s, ic), ca(2, :, oc, s, ic), mm(oc, s, ic), zc(j), sqrt_tail(zc(j), ainf, binf))
            do kk = 0, 2
              q(kk, lorb(oc), s, ic) = q(kk, lorb(oc), s, ic) - worb(oc) * aimag(wc(j) * (zc(j) - pps(ic, s)%enu(lorb(oc)))**kk * gg) / pi
            end do
          end do
        end do
      end do
    end do
  end subroutine occupied
end program probe_feco_depth
