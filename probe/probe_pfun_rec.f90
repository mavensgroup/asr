!> probe_pfun_rec: bcc Fe (k-space ASA potential) on a small periodic supercell; occupied s/p/d counts over the closed
!> contour with EXACT potential functions: (a) dense inverse of P - S, (b) ed_ordered_green recursion; and the same with
!> the parametrised P for reference.
program probe_pfun_rec
  use asr_kinds, only: dp, pi
  use asr_asa
  use asr_atom, only: shell_t
  use asr_lmto_cell
  use asr_radial
  use asr_pfun
  use asr_endep
  use asr_lmto_ham, only: potpar_t
  use asr_augspace, only: graph_t
  use asr_augorb, only: orb_alloy_t
  use asr_lmto_alloy
  use asr_lmto_realization
  use asr_strconst, only: hardcore_alpha
  use asr_contour, only: contour_nodes_closed
  implicit none
  type(asa_input_t) :: oi
  type(asa_result_t) :: orr
  type(cell_t) :: cell
  type(potpar_t) :: pp
  type(graph_t) :: g
  type(orb_alloy_t) :: mdl
  type(lmto_real_op) :: rop
  type(ed_params_t) :: sp
  type(ed_work_t) :: wk
  integer, allocatable :: coord(:, :), ipiv(:)
  logical, allocatable :: isa(:)
  complex(dp), allocatable :: amat(:, :), work(:)
  integer, parameter :: na = 40, nv = 8, lo9(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2], iorb(3) = [1, 3, 5]
  complex(dp) :: zq(na + nv), wq(na + nv), pf(0:2), mu2(0:2), lam(0:2), gr, gd, cnt(3, 2, 2)
  real(dp) :: w, al(0:2), al9(9), eb, ef, al4(0:3)
  integer :: l, j, i, k, kc, ll, nd, info, iex, ist
  allocate (oi%cls(1))
  oi%cls(1)%z = 26.0_dp; oi%cls(1)%nval = 8; oi%cls(1)%core = core_shells_for([4, 4, 3, 4]); oi%cls(1)%nvl = [4, 4, 3]
  oi%cls(1)%atom_shells = [oi%cls(1)%core, shell_t(4, 0, 2.0_dp), shell_t(3, 2, 6.0_dp)]
  oi%alat = 2.8665_dp / 0.529177210903_dp; oi%nk = 8; oi%maxit = 30; oi%nspin = 1
  cell%nat = 1; cell%avec = reshape([-0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp], [3, 3])
  allocate (cell%tau(3, 1), cell%isp(1)); cell%tau = 0.0_dp; cell%isp = 1
  call setup_cell(cell); oi%cls(1)%sws = cell%w * oi%alat
  call asa_scf(oi, cell, orr)
  print '(a,f9.4)', 'k-space ASA E_F ', orr%ef
  w = cell%w * oi%alat
  al4 = hardcore_alpha(0.7_dp); al = al4(0:2)
  do l = 0, 2
    pp%enu(l) = orr%cls(1)%pp(l, 1)%enu; pp%c(l) = orr%cls(1)%pp(l, 1)%c; pp%sqrdel(l) = orr%cls(1)%pp(l, 1)%sqrdel
    pp%p(l) = orr%cls(1)%pp(l, 1)%p; pp%gamma(l) = orr%cls(1)%pp(l, 1)%gamma; pp%alpha_rep(l) = pp%gamma(l)
  end do
  do i = 1, 9
    al9(i) = al(lo9(i)); sp%c(i) = pp%c(lo9(i)); sp%del(i) = pp%sqrdel(lo9(i))**2; sp%gam(i) = pp%gamma(lo9(i)); sp%enu(i) = pp%enu(lo9(i))
  end do
  call build_lmto_alloy('bcc', 3, 1.2_dp, 0.7_dp, pp, pp, 0.5_dp, g, coord, mdl)
  allocate (isa(g%n)); isa = .true.
  call build_lmto_ordered(g, mdl, isa, rop)
  nd = 9 * g%n
  call ed_work_init(wk, nd)
  allocate (amat(nd, nd), ipiv(nd), work(nd))
  eb = minval(pp%enu) - 1.2_dp; ef = orr%ef
  call contour_nodes_closed(eb, ef, 0.01_dp, na, nv, zq, wq)
  cnt = 0.0_dp
  do iex = 1, 2                                     ! 1 = exact, 2 = parametrised
    sp%exact = (iex == 1)
    do j = 1, na + nv
      if (iex == 1) then
        do l = 0, 2
          call pfun_exact(orr%cls(1)%mesh, orr%cls(1)%v(:, 1), l, zq(j), 26.0_dp, c_ry, w, al(l), pf(l), mu2(l), lam(l))
        end do
        do i = 1, 9; sp%pfx(i) = pf(lo9(i)); sp%mux(i) = mu2(lo9(i)); sp%lmx(i) = lam(lo9(i)); end do
      end if
      call ed_functions(zq(j), sp, al9, sp%pfx, sp%mux, sp%lmx)     ! parametrised values (or exact passthrough)
      ! dense
      amat = 0.0_dp
      do i = 1, g%n
        do ll = 1, 9; amat(9 * (i - 1) + ll, 9 * (i - 1) + ll) = sp%pfx(ll); end do
        amat(9 * (i - 1) + 1:9 * i, 9 * (i - 1) + 1:9 * i) = amat(9 * (i - 1) + 1:9 * i, 9 * (i - 1) + 1:9 * i) - mdl%block0
        do k = 1, size(mdl%block, 3)
          kc = g%nbr(g%nbr_ptr(i) + k - 1)
          amat(9 * (i - 1) + 1:9 * i, 9 * (kc - 1) + 1:9 * kc) = amat(9 * (i - 1) + 1:9 * i, 9 * (kc - 1) + 1:9 * kc) - mdl%block(:, :, k)
        end do
      end do
      call zgetrf(nd, nd, amat, nd, ipiv, info)
      do ist = 1, 3
        work = 0.0_dp; work(iorb(ist)) = 1.0_dp
        call zgetrs('N', nd, 1, amat, nd, ipiv, work, nd, info)
        gd = sp%lmx(iorb(ist)) + sp%mux(iorb(ist)) * work(iorb(ist))
        call ed_ordered_green(rop, wk, zq(j), sp, sp, al9, 1, iorb(ist), 24, gr)
        cnt(ist, 1, iex) = cnt(ist, 1, iex) + wq(j) * gd
        cnt(ist, 2, iex) = cnt(ist, 2, iex) - 0.0_dp + wq(j) * gr
      end do
      sp%exact = (iex == 1)
    end do
  end do
  print '(a,3f9.4,a,3f9.4)', 'exact  s p_z t2g occupation: dense', -aimag(cnt(:, 1, 1)) / pi, '   recursion', -aimag(cnt(:, 2, 1)) / pi
  print '(a,3f9.4,a,3f9.4)', 'param  s p_z t2g occupation: dense', -aimag(cnt(:, 1, 2)) / pi, '   recursion', -aimag(cnt(:, 2, 2)) / pi
  print '(a,3f9.4)', 'k-space ASA q0 s p d per orbital (spin-summed/2l+1/2):', (orr%cls(1)%q(0, l, 1) / (2 * (2 * l + 1)), l=0, 2)
end program probe_pfun_rec
