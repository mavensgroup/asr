!> Bloch spectral function through the production routine class_scf (asr_spectral, spectral_final): ordered fcc Cu.
!> For one species the scaled operator is z - H^(2) itself.  The reference is the dense inverse of z - H^(2) on the same node
!> set (unreduced operator, same levels), which is the limit of a long chain (docs/DERIVATION_levels_steps.md, Result S2).  Its
!> diagonal elements G_{0L,RL} on the targets within rmax, unfolded with the Bloch phases, must give the A(k, E) of class_scf
!> (targets, mirror images, symmetrisation, D2h readout, unfolding).  Also: no reciprocity defect.  The node set here is the
!> whole periodic cell (256 sites), so the dense sum over every site is the exact A(k, E) of the cell at the commensurate k of the
!> path, which is non-negative.  The sharp cutoff at rmax gives small negative values in this ordered metal (ringing).
program test_spectral_scf
  use asr_kinds, only: dp, pi
  use asr_driver, only: run_input_t, read_run_input, class_species
  use asr_class_scf, only: class_scf, class_scf_result_t
  use asr_spectral, only: spectral_unfold
  use asr_augspace, only: graph_t
  use asr_augcls
  use asr_edc
  use asr_augorb, only: orb_alloy_t
  use asr_lmto_alloy, only: build_lmto_alloy
  use asr_lmto_ham, only: potpar_t
  use asr_strconst, only: standard_alpha
  use asr_lapack, only: zgesv
  use asr_testing
  implicit none
  type(run_input_t) :: inp
  type(class_scf_result_t) :: r
  character(len=:), allocatable :: msg
  integer :: u
  character(len=*), parameter :: f = 'build/test_spectral_scf_tmp.toml'
  character(len=1), parameter :: nl = new_line('a')
  character(len=*), parameter :: cu = '[structure]'//nl//'lattice = "fcc"'//nl//'alat = 3.61'//nl//'centering = "F"'//nl// &
    '[[species]]'//nl//'symbol = "Cu"'//nl//'z = 29'//nl//'nval = 11'//nl//'konfig = [4, 4, 3, 4]'//nl//'occupations = [1.0, 10.0]'//nl// &
    '[[sublattice]]'//nl//'name = "fcc"'//nl//'positions = [[0.0, 0.0, 0.0]]'//nl//'species = ["Cu"]'//nl// &
    '[scf]'//nl//'spin = 1'//nl//'xc = "vbh"'//nl//'mixing = 0.3'//nl//'max_iterations = 5'//nl//'tolerance = 1.0e3'//nl// &
    '[asr]'//nl//'levels = 3'//nl//'recursion_steps = 8'//nl//'supercell = 4'//nl//'contour_points = 16'//nl// &
    '[spectral]'//nl//'kpath = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.5, 0.0]]'//nl//'labels = ["G", "X", "W"]'//nl// &
    'nk = 2'//nl//'emin = -0.4'//nl//'emax = 0.4'//nl//'npoints = 9'//nl//'eta = 0.03'//nl//'rmax = 1.5'//nl//'steps = 200'//nl
  type(potpar_t) :: pp(1)
  type(graph_t) :: g
  type(orb_alloy_t) :: m0
  type(class_alloy_t) :: mdl
  type(augcls_op) :: opu
  type(edn_t) :: en
  integer, allocatable :: coord(:, :), cls(:), perm1(:, :), tsite(:), tn(:, :), tg(:, :), ipiv(:)
  logical, allocatable :: fl(:)
  real(dp), allocatable :: hm(:, :), x(:), y(:), tdisp(:, :), v(:), aref(:, :), afull(:, :), vall(:), dall(:, :)
  complex(dp), allocatable :: km(:, :), rhs(:, :)
  real(dp) :: alpha(0:2), d1(9, 9, 1), err, amax
  integer :: l, i, j, ik, nd, nt, info, h(3), ne, nkp

  print '(a)', 'test_spectral_scf'
  open (newunit=u, file=f, status='replace'); write (u, '(a)') cu; close (u)
  call read_run_input(f, inp, msg)
  call check(len(msg) == 0, 'fcc Cu input with [spectral] accepted')
  if (len(msg) > 0) then; print '(a)', msg; call report(); stop; end if
  call class_species(inp)
  inp%ci%ndos = 0
  call class_scf(inp%ci, r, verbose=.true.)
  call check(allocated(r%spec_a), 'class_scf returns the spectral function')
  ne = size(r%spec_e); nkp = size(r%spec_dist)
  call check(ne == 9 .and. nkp == 5, 'spectral grid: 9 energies, 5 k points')
  print '(a,i0,a,es10.2)', '  sharp cutoff at rmax = 1.5: points with A < 0: ', r%spec_nneg, ', min A / max A ', r%spec_amin
  call check_close(r%spec_rdef, 0.0_dp, 1.0e-8_dp, 'ordered: G(0L,RL) = G(0L,-RL) (inversion)')

  ! ---- dense reference on the same node set, from the final potential parameters of the run
  do l = 0, 2
    pp(1)%enu(l) = r%comp(1)%pp(l, 1)%enu; pp(1)%c(l) = r%comp(1)%pp(l, 1)%c; pp(1)%sqrdel(l) = r%comp(1)%pp(l, 1)%sqrdel
    pp(1)%p(l) = r%comp(1)%pp(l, 1)%p; pp(1)%gamma(l) = r%comp(1)%pp(l, 1)%gamma
  end do
  pp(1)%alpha_rep = pp(1)%gamma
  alpha = standard_alpha()
  call build_lmto_alloy('fcc', inp%ci%lsuper, inp%ci%rcl, 0.7_dp, pp(1), pp(1), 0.5_dp, g, coord, m0, max(inp%ci%rcl, inp%ci%rinv), alpha)
  mdl%nl = 9; mdl%nsp = 1; mdl%ncls = 1
  allocate (mdl%nsc(1), mdl%spc(1, 1), mdl%conc(1, 1)); mdl%nsc = 1; mdl%spc = 1; mdl%conc = 1.0_dp
  mdl%block = m0%block; mdl%block0 = m0%block0
  call edc_setup(mdl, pp, alpha, 1, en)
  allocate (cls(g%n), fl(g%n)); cls = 1; fl = .true.; fl(1) = .false.
  allocate (perm1(g%n, 1)); perm1(:, 1) = [(i, i=1, g%n)]
  d1 = 0.0_dp; do l = 1, 9; d1(l, l, 1) = 1.0_dp; end do
  nt = 0; allocate (tsite(g%n), tdisp(3, g%n))
  do i = 1, g%n
    h = modulo(coord(:, i) + inp%ci%lsuper, 2 * inp%ci%lsuper) - inp%ci%lsuper
    if (real(sum(h**2), dp) <= 4.0_dp * 1.5_dp**2 + 1.0e-9_dp) then; nt = nt + 1; tsite(nt) = i; tdisp(:, nt) = 0.5_dp * h; end if
  end do
  tsite = tsite(1:nt); tdisp = tdisp(:, 1:nt)
  call check(nt == r%spec_nt, 'reference: the same number of target sites')
  call build_augcls_operator(g, mdl, cls, perm1, d1, 1, inp%ci%nlevels, opu, maxrank=inp%ci%maxrank, flippable=fl, &
                             tsite=tsite, tnode=tn, tgop=tg)
  call augcls_set_params(opu, mdl); call augcls_freeze_site(opu, 1, 1)
  nd = opu%dim()
  print '(a,i0,a,i0)', '  reference: unreduced nodes ', opu%nn, ', dimension ', nd
  allocate (hm(nd, nd), x(nd), y(nd))
  do j = 1, nd
    x = 0.0_dp; x(j) = 1.0_dp
    call opu%apply(x, y)
    hm(:, j) = y
  end do
  do i = 1, opu%nn                                     ! one species: e_s(z) = E_nu, the term edc applies separately
    do l = 1, 9; hm(9 * (i - 1) + l, 9 * (i - 1) + l) = hm(9 * (i - 1) + l, 9 * (i - 1) + l) + en%enu(l, 1); end do
  end do
  allocate (km(nd, nd), rhs(nd, 9), ipiv(nd), v(nt), aref(ne, nkp), afull(ne, nkp), vall(opu%nn), dall(3, opu%nn))
  aref = 0.0_dp; afull = 0.0_dp
  do i = 1, opu%nn                                     ! every node of the cell: displacement of its site (minimum image)
    h = modulo(coord(:, opu%site(i)) + inp%ci%lsuper, 2 * inp%ci%lsuper) - inp%ci%lsuper; dall(:, i) = 0.5_dp * h
  end do
  do j = 1, ne
    km = -hm
    do i = 1, nd; km(i, i) = km(i, i) + cmplx(r%spec_e(j), inp%ci%spectral%eta, dp); end do
    rhs = (0.0_dp, 0.0_dp)
    do l = 1, 9; rhs(l, l) = (1.0_dp, 0.0_dp); end do  ! the origin is node 1
    call zgesv(nd, 9, km, nd, ipiv, rhs, nd, info)
    if (info /= 0) error stop 'test_spectral_scf: zgesv'
    do l = 1, 9
      do i = 1, nt; v(i) = -aimag(rhs(9 * (tn(1, i) - 1) + l, l)) / pi; end do
      do i = 1, opu%nn; vall(i) = -aimag(rhs(9 * (i - 1) + l, l)) / pi; end do
      do ik = 1, nkp
        aref(j, ik) = aref(j, ik) + spectral_unfold(tdisp, v, r%spec_k(:, ik))
        afull(j, ik) = afull(j, ik) + spectral_unfold(dall, vall, r%spec_k(:, ik))
      end do
    end do
  end do
  amax = maxval(abs(aref))
  err = maxval(abs(sum(r%spec_a(:, :, :, 1), dim=3) - aref)) / amax
  print '(a,es10.2,a,f8.3)', '  max |A(class_scf) - A(dense)| / max A = ', err, ',  max A = ', amax
  call check_close(err, 0.0_dp, 1.0e-6_dp, 'A(k, E) of class_scf = dense inverse on the node set (200-step chains)')
  call check(opu%nn == g%n, 'reference: the node set is the whole periodic cell')
  print '(a,es10.2)', '  sum over every site of the cell: min A / max A = ', minval(afull) / maxval(afull)
  call check(minval(afull) > 0.0_dp, 'sum over every site of the cell (exact A of the cell at commensurate k): A > 0')
  call report()
end program test_spectral_scf
