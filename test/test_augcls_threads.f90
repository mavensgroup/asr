!> Threads in the class operator (asr_augcls) and its chains:
!>  (B) build_augcls_operator gives the same node set and tables with 1 and with 8 threads (parallel closure with serial
!>      insertion in the serial order; parallel table pass);
!>  (S) a copy of the operator shares the node tables (pointers) and keeps its own work space;
!>  (C) an edc chain outside any parallel region (its inner loops then run on all threads) gives bitwise the same G with 1 and
!>      with 8 threads.
!> Random bcc binary, flippable sites, levels 4, rank 2.
program test_augcls_threads
  use asr_kinds, only: dp
  use asr_augspace, only: graph_t
  use asr_augcls
  use asr_edc
  use asr_lmto_alloy, only: build_lmto_alloy, orb_alloy_t
  use asr_symmetry
  use asr_lmto_ham, only: potpar_t
  use asr_chain, only: ed_depth_tol, ed_tail
  use asr_testing
  !$ use omp_lib, only: omp_set_num_threads, omp_get_max_threads
  implicit none
  integer, parameter :: lsup = 4, nlev = 24
  type(potpar_t) :: pp(2)
  real(dp) :: alpha(0:2)
  type(graph_t) :: g
  type(orb_alloy_t) :: m0
  type(class_alloy_t) :: mdl
  type(augcls_op) :: op1, op8, cp
  type(edn_t) :: en
  type(edn_work_t) :: w1, w8
  integer, allocatable :: coord(:, :), cls(:), perm(:, :), mats(:, :, :)
  real(dp), allocatable :: d(:, :, :)
  logical, allocatable :: fl(:)
  complex(dp) :: g1(3), g8(3), z(3)
  integer :: l, i, nth0, j
  logical :: same

  print '(a)', 'test_augcls_threads'
  nth0 = 1
  !$ nth0 = omp_get_max_threads()
  do l = 0, 2
    pp(1)%enu(l) = 0.1_dp * l; pp(1)%c(l) = 0.3_dp - 0.15_dp * l; pp(1)%sqrdel(l) = 0.30_dp - 0.07_dp * l
    pp(1)%p(l) = 0.9_dp - 0.1_dp * l; pp(1)%gamma(l) = 0.25_dp - 0.1_dp * l
    pp(2)%enu(l) = 0.05_dp * l; pp(2)%c(l) = 0.45_dp - 0.1_dp * l; pp(2)%sqrdel(l) = 0.27_dp - 0.06_dp * l
    pp(2)%p(l) = 0.8_dp - 0.1_dp * l; pp(2)%gamma(l) = 0.2_dp - 0.08_dp * l
  end do
  pp(1)%alpha_rep = pp(1)%gamma; pp(2)%alpha_rep = pp(2)%gamma
  alpha = [0.3485_dp, 0.05303_dp, 0.010714_dp]
  call build_lmto_alloy('bcc', lsup, 1.2_dp, 0.7_dp, pp(1), pp(2), 0.5_dp, g, coord, m0, alpha_in=alpha)
  mdl%nl = 9; mdl%nsp = 2; mdl%ncls = 1
  allocate (mdl%nsc(1), mdl%spc(2, 1), mdl%conc(2, 1))
  mdl%nsc = 2; mdl%spc(:, 1) = [1, 2]; mdl%conc(:, 1) = [0.5_dp, 0.5_dp]
  mdl%block = m0%block; mdl%block0 = m0%block0
  allocate (fl(g%n), cls(g%n)); fl = .true.; fl(1) = .false.; cls = 1
  call class_ops(coord, 2 * lsup, cls, perm, mats)
  call orbital_rep(mats, d)
  call edc_setup(mdl, pp, alpha, 1, en)

  ! ---- (B) the build with 1 and with 8 threads
  !$ call omp_set_num_threads(1)
  call build_augcls_operator(g, mdl, cls, perm, d, 1, 4, op1, maxrank=2, flippable=fl)
  !$ call omp_set_num_threads(8)
  call build_augcls_operator(g, mdl, cls, perm, d, 1, 4, op8, maxrank=2, flippable=fl)
  print '(a,i0,a,i0,a,i0)', '  random bcc binary: nodes ', op1%nn, ', slots ', op1%nslot, ', group order ', size(perm, 2)
  call check(op1%nn > 1000 .and. op1%nn == op8%nn, 'the same number of nodes with 1 and 8 threads')
  same = all(op1%hopr == op8%hopr) .and. all(op1%hopg == op8%hopg) .and. all(op1%flipr == op8%flipr) .and. &
         all(op1%flipg == op8%flipg) .and. all(op1%lst == op8%lst) .and. all(op1%site == op8%site) .and. &
         all(op1%ncl == op8%ncl) .and. all(op1%level == op8%level) .and. all(op1%norb == op8%norb) .and. &
         all(op1%stab_ptr == op8%stab_ptr) .and. all(op1%stab == op8%stab)
  call check(same, 'every node table is the same with 1 and 8 threads (same node order)')

  ! ---- (S) a copy shares the tables
  cp = op1
  call check(associated(cp%hopr, op1%hopr) .and. associated(cp%stab, op1%stab) .and. associated(cp%sq, op1%sq), &
             'a copy of the operator shares the node tables')

  ! ---- (C) chains with 1 and 8 threads (inner parallel loops active)
  call augcls_set_params(op1, mdl); call augcls_freeze_site(op1, 1, 1); call augcls_set_start(op1, 5, d)
  ed_depth_tol = 0.0_dp; ed_tail = 'hull'
  call edn_work_init(w1, op1%dim()); call edn_work_init(w8, op1%dim())
  z = [cmplx(-0.3_dp, 0.05_dp, dp), cmplx(0.1_dp, 0.01_dp, dp), cmplx(0.4_dp, 0.2_dp, dp)]
  !$ call omp_set_num_threads(1)
  do j = 1, 3; call edc_green(op1, w1, en, z(j), 5, nlev, g1(j)); end do
  !$ call omp_set_num_threads(8)
  do j = 1, 3; call edc_green(op1, w8, en, z(j), 5, nlev, g8(j)); end do
  print '(a,3(2es12.4,2x))', '  G(z) ', g1
  call check(all(g1 == g8), 'edc chains: G bitwise the same with 1 and 8 threads')
  !$ call omp_set_num_threads(nth0)
  call augcls_free(op1); call augcls_free(op8)
  call report()
end program test_augcls_threads
