!> End-to-end pure-limit gate of the LMTO augmented-space machinery: fcc Ag with A = B (no disorder) and
!> with x = 1, real-space recursion on the augmented-space operator (graph + screened-S bond blocks + H = E + h - hoh)
!> must give the same recursion coefficients as the ordered k-space Hamiltonian (hk_grid) for the same parameters.
!> Parameters come from our own ASA SCF for Ag.
program test_alloy_ordered
  use asr_kinds, only: dp
  use asr_asa
  use asr_atom, only: shell_t
  use asr_lmto_cell
  use asr_lmto_ham, only: potpar_t, kgrid_op, hk_grid
  use asr_augspace, only: graph_t
  use asr_augorb
  use asr_lmto_alloy
  use asr_strconst, only: lattice_cluster, screened_cluster, avg_ws_radius, hardcore_alpha
  use asr_recursion, only: lanczos
  use asr_testing
  implicit none
  type(asa_input_t) :: inp
  type(asa_result_t) :: res
  type(cell_t) :: cell
  type(potpar_t) :: pp
  type(graph_t) :: g
  type(orb_alloy_t) :: mdl
  type(augorb_op) :: op
  type(kgrid_op) :: kop
  integer, allocatable :: coord(:, :)
  real(dp), allocatable :: v0(:), vk(:), vec(:, :), blk(:, :, :)
  real(dp) :: a1(8), b1(8), a2(8), b2(8), al4(0:3), alpha(0:2), bv(3, 3), dmax(2)
  integer :: l, m, io, ic, iorb(2), nd
  integer, parameter :: depth = 6

  print '(a)', 'test_alloy_ordered'
  allocate (inp%cls(1))
  inp%cls(1)%z = 47.0_dp; inp%cls(1)%nval = 11
  inp%cls(1)%core = core_shells_for([5, 5, 4, 4]); inp%cls(1)%nvl = [5, 5, 4]
  inp%cls(1)%atom_shells = [inp%cls(1)%core, shell_t(5, 0, 1.0_dp), shell_t(4, 2, 10.0_dp)]
  inp%alat = 7.72898293_dp; inp%nk = 10; inp%maxit = 60
  cell%nat = 1
  cell%avec = reshape([0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp], [3, 3])
  allocate (cell%tau(3, 1), cell%isp(1)); cell%tau = 0.0_dp; cell%isp = 1
  call setup_cell(cell); inp%cls(1)%sws = cell%w * inp%alat
  call asa_scf(inp, cell, res)
  call check(res%converged, 'Ag SCF (parameters for the gate)')
  do l = 0, 2
    pp%enu(l) = res%cls(1)%pp(l, 1)%enu; pp%c(l) = res%cls(1)%pp(l, 1)%c; pp%sqrdel(l) = res%cls(1)%pp(l, 1)%sqrdel
    pp%p(l) = res%cls(1)%pp(l, 1)%p; pp%gamma(l) = res%cls(1)%pp(l, 1)%gamma; pp%alpha_rep(l) = res%cls(1)%pp(l, 1)%gamma
  end do

  ! augmented-space operator: A = B = Ag (no disorder), supercell large enough for `depth` levels
  call build_lmto_alloy('fcc', 28, 1.2_dp, 0.7_dp, pp, pp, 0.5_dp, g, coord, mdl)
  call build_augorb_operator(g, mdl, 1, depth, op)
  nd = op%dim()
  print '(a,i0,a,i0,a,i0)', '    sites ', g%n, '   slots ', size(mdl%block, 3), '   nodes ', op%nn

  ! k-space reference operator with the SAME cluster cutoff
  al4 = hardcore_alpha(0.7_dp); alpha = al4(0:2)
  call lattice_cluster('fcc', 1.2_dp, vec)
  call screened_cluster(vec, avg_ws_radius('fcc'), alpha, blk)
  bv = reshape([-1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp], [3, 3])
  call hk_grid(24, bv, vec, blk, pp, alpha, kop)
  allocate (vk(kop%dim()))
  iorb = [1, 7]                                  ! s and an e_g (3z^2-r^2) orbital
  do io = 1, 2
    allocate (v0(nd)); v0 = 0.0_dp; v0(iorb(io)) = 1.0_dp          ! node 1 (origin, no flips), orbital iorb
    call lanczos(op, v0, depth, a1, b1, m, reorth=.true.)
    deallocate (v0)
    vk = 0.0_dp; vk(iorb(io):9 * kop%nk3:9) = 1.0_dp
    call lanczos(kop, vk, depth, a2, b2, m, reorth=.true.)
    print '(a,6f11.5)', '      a real :', a1(1:6)
    print '(a,6f11.5)', '      a kspc :', a2(1:6)
    print '(a,6f11.5)', '      b real :', b1(1:5)
    print '(a,6f11.5)', '      b kspc :', b2(1:5)
    dmax(io) = max(maxval(abs(a1(1:depth) - a2(1:depth))), maxval(abs(b1(1:depth - 1) - b2(1:depth - 1))))
    print '(a,i0,a,es9.2)', '    orbital ', iorb(io), ':  max |coefficient difference| (real space vs k-space) = ', dmax(io)
  end do
  call check_close(dmax(1), 0.0_dp, 1.0e-6_dp, 's orbital: augmented-space recursion == ordered k-space recursion')
  call check_close(dmax(2), 0.0_dp, 1.0e-6_dp, 'e_g orbital: augmented-space recursion == ordered k-space recursion')
  call report()
end program test_alloy_ordered
