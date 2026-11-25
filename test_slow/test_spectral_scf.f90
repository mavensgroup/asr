!> Bloch spectral function through the production routine class_scf (asr_spectral, spectral_final): ordered fcc Cu.
!> For one species the scaled operator is z - H^(2) itself.  The reference is the dense inverse of z - H^(2) on the same node
!> set (unreduced operator on the same site cluster, same levels; asr_crystal_graph).  Levels 2 with nearest-neighbour bonds keep
!> the node set small (about 380 sites), so chains of 1000 steps exhaust it and must equal the dense inverse.  Its
!> diagonal elements G_{0L,RL} on the targets within rmax, unfolded with the Bloch phases, must give the A(k, E) of class_scf
!> (targets, mirror images, symmetrisation, D2h readout, unfolding).  Also: no reciprocity defect.  The default ball window
!> (asr_spectral) turns the cutoff sum into the exact A convolved in k with a non-negative kernel, so A(k, E) of class_scf must
!> be non-negative (the sharp cutoff rings in this metal).  The bands of H^(2)(k) that class_scf returns depend only on the
!> structure constants and the parameters: on a periodic cell of 4 cubes (asr_lattices), where the node set is the whole cell
!> (256 sites), the dense sum over every site is the exact A(k, E) of the cell at the commensurate k of the path; it must be
!> positive and equal the Lorentzians of the bands.
program test_spectral_scf
  use asr_kinds, only: dp, pi
  use asr_driver, only: run_input_t, read_run_input, class_species
  use asr_class_scf, only: class_scf, class_scf_result_t
  use asr_spectral, only: spectral_unfold, spectral_window
  use asr_augspace, only: graph_t
  use asr_augcls
  use asr_edc
  use asr_lmto_alloy, only: orb_alloy_t, build_lmto_alloy
  use asr_lmto_ham, only: potpar_t
  use asr_strconst, only: standard_alpha
  use asr_lapack, only: zgesv
  use asr_crystal, only: structure_t
  use asr_crystal_symmetry, only: symop_t, space_group
  use asr_crystal_graph, only: crystal_graph_t, crystal_graph
  use asr_lmto_alloy, only: structure_blocks
  use asr_sublattice, only: decoration_structure
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
    '[asr]'//nl//'levels = 2'//nl//'recursion_steps = 8'//nl//'cluster_radius = 0.75'//nl//'contour_points = 16'//nl// &
    '[spectral]'//nl//'kpath = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.5, 0.0]]'//nl//'labels = ["G", "X", "W"]'//nl// &
    'nk = 2'//nl//'emin = -0.4'//nl//'emax = 0.4'//nl//'npoints = 9'//nl//'eta = 0.03'//nl//'rmax = 1.5'//nl//'steps = 1000'//nl
  type(potpar_t) :: pp(1)
  type(graph_t) :: g
  type(orb_alloy_t) :: m0
  type(class_alloy_t) :: mdl
  type(augcls_op) :: opu
  type(edn_t) :: en
  type(structure_t) :: st
  type(symop_t), allocatable :: sops(:)
  type(crystal_graph_t) :: cg
  integer, parameter :: lcell = 4, lvcell = 4         ! periodic cell of the band check: 4 cubes; levels that cover it
  real(dp), allocatable :: blk(:, :, :, :), blk0(:, :, :)
  integer, allocatable :: typ(:)
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
  print '(a,i0,a,es10.2)', '  ball window, rmax = 1.5: points with A < 0: ', r%spec_nneg, ', min A / max A ', r%spec_amin
  call check(trim(inp%ci%spectral%window) == 'ball', 'default Fourier window: ball')
  call check_close(r%spec_rdef, 0.0_dp, 1.0e-8_dp, 'ordered: G(0L,RL) = G(0L,-RL) (inversion)')

  ! ---- dense reference on the same node set, from the final potential parameters of the run
  do l = 0, 2
    pp(1)%enu(l) = r%comp(1)%pp(l, 1)%enu; pp(1)%c(l) = r%comp(1)%pp(l, 1)%c; pp(1)%sqrdel(l) = r%comp(1)%pp(l, 1)%sqrdel
    pp(1)%p(l) = r%comp(1)%pp(l, 1)%p; pp(1)%gamma(l) = r%comp(1)%pp(l, 1)%gamma
  end do
  pp(1)%alpha_rep = pp(1)%gamma
  alpha = standard_alpha()
  ! (1) the node set of class_scf: the site cluster of the structure, unreduced; targets within rmax
  call decoration_structure(inp%ci%dec, st, msg)
  call space_group(st, sops)
  call crystal_graph(st, sops, inp%ci%rcl, inp%ci%nlevels, cg, msg, start=1, rmin=inp%ci%spectral%rmax)
  if (len(msg) > 0) error stop msg
  call structure_blocks(st, cg, max(inp%ci%rcl, inp%ci%rinv), alpha, blk, blk0, typ)
  mdl%nl = 9; mdl%nsp = 1; mdl%ncls = 1
  allocate (mdl%nsc(1), mdl%spc(1, 1), mdl%conc(1, 1)); mdl%nsc = 1; mdl%spc = 1; mdl%conc = 1.0_dp
  mdl%block = blk; mdl%block0 = blk0
  call edc_setup(mdl, pp, alpha, 1, en)
  allocate (cls(cg%n), fl(cg%n)); cls = 1; fl = .true.; fl(1) = .false.
  allocate (perm1(cg%n, 1)); perm1(:, 1) = [(i, i=1, cg%n)]
  d1 = 0.0_dp; do l = 1, 9; d1(l, l, 1) = 1.0_dp; end do
  tsite = pack([(i, i=1, cg%n)], [(4.0_dp * sum(cg%disp(:, i)**2) <= 4.0_dp * 1.5_dp**2 + 1.0e-9_dp, i=1, cg%n)])
  nt = size(tsite); tdisp = cg%disp(:, tsite)
  call check(nt == r%spec_nt, 'reference: the same number of target sites')
  call build_augcls_operator(cg%g, mdl, cls, perm1, d1, 1, inp%ci%nlevels, opu, maxrank=inp%ci%maxrank, flippable=fl, &
                             tsite=tsite, tnode=tn, tgop=tg, stype=typ(cg%basis))
  call augcls_set_params(opu, mdl); call augcls_freeze_site(opu, 1, 1)
  call dense_h(opu, en, hm)
  nd = opu%dim()
  print '(a,i0,a,i0)', '  reference on the cluster: unreduced nodes ', opu%nn, ', dimension ', nd
  allocate (km(nd, nd), rhs(nd, 9), ipiv(nd), v(nt), aref(ne, nkp), afull(ne, nkp))
  aref = 0.0_dp; afull = 0.0_dp
  do j = 1, ne
    km = -hm
    do i = 1, nd; km(i, i) = km(i, i) + cmplx(r%spec_e(j), inp%ci%spectral%eta, dp); end do
    rhs = (0.0_dp, 0.0_dp)
    do l = 1, 9; rhs(l, l) = (1.0_dp, 0.0_dp); end do  ! the origin is node 1
    call zgesv(nd, 9, km, nd, ipiv, rhs, nd, info)
    if (info /= 0) error stop 'test_spectral_scf: zgesv'
    do l = 1, 9
      do i = 1, nt
        v(i) = -aimag(rhs(9 * (tn(1, i) - 1) + l, l)) / pi * spectral_window('ball', norm2(tdisp(:, i)), 1.5_dp)
      end do
      do ik = 1, nkp
        aref(j, ik) = aref(j, ik) + spectral_unfold(tdisp, v, r%spec_k(:, ik))
      end do
    end do
  end do
  call augcls_free(opu); deallocate (km, rhs, ipiv, cls, fl, perm1, hm)
  ! (2) the periodic cell of lcell cubes (asr_lattices): the node set is the whole cell, the sum over every site is exact
  call build_lmto_alloy('fcc', lcell, inp%ci%rcl, 0.7_dp, pp(1), pp(1), 0.5_dp, g, coord, m0, max(inp%ci%rcl, inp%ci%rinv), alpha)
  call augcls_set_blocks(mdl, m0%block, m0%block0)
  allocate (cls(g%n), fl(g%n)); cls = 1; fl = .true.; fl(1) = .false.
  allocate (perm1(g%n, 1)); perm1(:, 1) = [(i, i=1, g%n)]
  call build_augcls_operator(g, mdl, cls, perm1, d1, 1, lvcell, opu, maxrank=inp%ci%maxrank, flippable=fl)
  call augcls_set_params(opu, mdl); call augcls_freeze_site(opu, 1, 1)
  call check(opu%nn == g%n, 'periodic cell: the node set is the whole cell')
  call dense_h(opu, en, hm)
  nd = opu%dim()
  allocate (km(nd, nd), rhs(nd, 9), ipiv(nd), vall(opu%nn), dall(3, opu%nn))
  do i = 1, opu%nn                                     ! every node of the cell: displacement of its site (minimum image)
    h = modulo(coord(:, opu%site(i)) + lcell, 2 * lcell) - lcell; dall(:, i) = 0.5_dp * h
  end do
  do j = 1, ne
    km = -hm
    do i = 1, nd; km(i, i) = km(i, i) + cmplx(r%spec_e(j), inp%ci%spectral%eta, dp); end do
    rhs = (0.0_dp, 0.0_dp)
    do l = 1, 9; rhs(l, l) = (1.0_dp, 0.0_dp); end do
    call zgesv(nd, 9, km, nd, ipiv, rhs, nd, info)
    if (info /= 0) error stop 'test_spectral_scf: zgesv'
    do l = 1, 9
      do i = 1, opu%nn; vall(i) = -aimag(rhs(9 * (i - 1) + l, l)) / pi; end do
      do ik = 1, nkp
        afull(j, ik) = afull(j, ik) + spectral_unfold(dall, vall, r%spec_k(:, ik))
      end do
    end do
  end do
  amax = maxval(abs(aref))
  err = maxval(abs(sum(r%spec_a(:, :, :, 1), dim=3) - aref)) / amax
  print '(a,es10.2,a,f8.3)', '  max |A(class_scf) - A(dense)| / max A = ', err, ',  max A = ', amax
  call check_close(err, 0.0_dp, 1.0e-6_dp, 'A(k, E) of class_scf = dense inverse on the node set (chains exhaust it)')
  call check(minval(sum(r%spec_a(:, :, :, 1), dim=3)) > -1.0e-6_dp * amax, &
             'ball window: A(k, E) >= 0 (exact G of the cell convolved with a non-negative kernel)')
  print '(a,es10.2)', '  sum over every site of the cell: min A / max A = ', minval(afull) / maxval(afull)
  call check(minval(afull) > 0.0_dp, 'sum over every site of the cell (exact A of the cell at commensurate k): A > 0')
  ! the H^(2)(k) bands of class_scf give the same exact A of the cell: A(k, E) = sum_n -Im 1/(E + i eta - e_n(k)) / pi
  call check(allocated(r%spec_bands), 'ordered, one species: class_scf returns the bands of H^(2)(k)')
  err = 0.0_dp
  do ik = 1, nkp
    do j = 1, ne
      err = max(err, abs(afull(j, ik) + sum(aimag(1.0_dp / (cmplx(r%spec_e(j), inp%ci%spectral%eta, dp) - r%spec_bands(:, ik, 1)))) / pi))
    end do
  end do
  print '(a,es10.2)', '  max |A(bands) - A(full cell)| / max A = ', err / maxval(afull)
  call check_close(err / maxval(afull), 0.0_dp, 1.0e-8_dp, 'bands of H^(2)(k) = exact A of the cell (Lorentzians of width eta)')
  call nias_ordered_check()
  call report()
contains

  !> Ordered NiAs (Co on cations, Al on anions; 4 basis sites, 2 sublattices, 3 translation-inversion origin classes).
  !> Checks multi-origin unfolding against (1) the symmetrised dense inverse on the 3 unreduced origin clusters and (2) the exact
  !> periodic 3x3x2 supercell (72 sites) and weighted H^(2)(k) Lorentzians at commensurate wavevectors, including k = (1/3) b_1 +
  !> (1/2) b_3 where no symmetry maps cation 1 to cation 2 while preserving +/- k.
  subroutine nias_ordered_check()
    use asr_spectral, only: spectral_origin_classes
    character(len=*), parameter :: fn = 'build/test_spectral_nias_tmp.toml'
    character(len=*), parameter :: nias_in = &
      '[structure]'//nl//'lattice = "nias"'//nl//'alat = 3.90'//nl//'c_over_a = 1.62'//nl// &
      '[[species]]'//nl//'symbol = "Co"'//nl//'z = 27'//nl//'nval = 9'//nl//'konfig = [4, 4, 3, 4]'//nl//'occupations = [2.0, 7.0]'//nl// &
      '[[species]]'//nl//'symbol = "Al"'//nl//'z = 13'//nl//'nval = 3'//nl//'konfig = [3, 3, 3, 4]'//nl//'occupations = [2.0, 0.0]'//nl// &
      '[[sublattice]]'//nl//'name = "cat"'//nl//'species = ["Co"]'//nl// &
      '[[sublattice]]'//nl//'name = "an"'//nl//'species = ["Al"]'//nl// &
      '[scf]'//nl//'spin = 1'//nl//'xc = "vbh"'//nl//'mixing = 0.2'//nl//'max_iterations = 2'//nl//'tolerance = 1.0e3'//nl// &
      '[asr]'//nl//'levels = 2'//nl//'recursion_steps = 8'//nl//'cluster_radius = 1.05'//nl//'contour_points = 8'//nl// &
      '[spectral]'//nl// &
      'kpath = [[0.0, 0.0, 0.0], [0.3333333333333333, 0.19245008972987526, 0.30864197530864196], '// &
      '[0.3333333333333333, 0.5773502691896258, 0.0], [0.3333333333333333, 0.5773502691896258, 0.30864197530864196]]'//nl// &
      'labels = ["G", "L3", "K", "H"]'//nl//'nk = 1'//nl//'emin = -0.5'//nl//'emax = 0.3'//nl// &
      'npoints = 7'//nl//'eta = 0.04'//nl//'rmax = 1.1'//nl//'steps = 600'//nl
    integer, parameter :: l1 = 3, l2 = 3, l3 = 2, nsup = l1 * l2 * l3 * 4, nds = 9 * nsup
    type(run_input_t) :: ni
    type(class_scf_result_t) :: nr
    type(potpar_t) :: pps(2)
    type(structure_t) :: nst
    type(symop_t), allocatable :: nops(:)
    type(crystal_graph_t) :: ncg(3)
    type(class_alloy_t) :: nmdl
    type(edn_t) :: nen(2)
    type(augcls_op) :: nop
    character(len=:), allocatable :: nmsg
    integer, allocatable :: bcl(:), ob(:), osg(:)
    real(dp), allocatable :: ow(:), nblk(:, :, :, :), nblk0(:, :, :), nhm(:, :), td(:, :, :), vv(:)
    integer, allocatable :: ntyp(:), ncls_(:), nperm(:, :), nts(:), ntn(:, :), ntg(:, :), nipiv(:)
    integer, allocatable :: ntk(:), tbase(:, :), tsign(:, :), tmir(:, :)
    logical, allocatable :: nfl(:)
    complex(dp), allocatable :: nkm(:, :), nrhs(:, :), gdg(:, :, :, :), gdg0(:, :, :, :)
    real(dp), allocatable :: smat(:, :), hmat(:, :), h2mat(:, :)
    real(dp) :: naref(7, 4), nafull(7, 4), nafull_1sub(7, 4), naband(7, 4)
    real(dp) :: rpos(3, nsup), dsup(3, nsup), vsup(nsup), cv(nds), dv(nds), ov(nds), ev0(nds)
    real(dp) :: al(0:2), dm1(9, 9, 1), nerr, nerr1, namax
    integer :: uidx(0:l1 - 1, 0:l2 - 1, 0:l3 - 1, 4), uorig(3)
    integer :: fu, uo, u2o, kap_c, ks_s, ll, ii, jj, kk, b1, b2, ks, i1, i2, i3, j1, j2, j3, u1, u2, ndo, nto, ntmax, inf
    complex(dp) :: zc

    open (newunit=fu, file=fn, status='replace'); write (fu, '(a)') nias_in; close (fu)
    call read_run_input(fn, ni, nmsg)
    call check(len(nmsg) == 0, 'ordered NiAs input with [spectral] accepted')
    if (len(nmsg) > 0) then; print '(a)', nmsg; return; end if
    call class_species(ni)
    ni%ci%ndos = 0
    call class_scf(ni%ci, nr, verbose=.false.)
    call check(allocated(nr%spec_a) .and. size(nr%spec_nodes) == 3, 'ordered NiAs: 3 spectral origin classes')
    call check(allocated(nr%spec_bands) .and. size(nr%spec_bands, 1) == 36, 'ordered NiAs: 36 H^(2)(k) bands')
    nerr = 0.0_dp
    do kk = 1, 4
      nerr = max(nerr, abs(sum(nr%spec_bwgt(:, kk, 1)) - 9.0_dp))
    end do
    call check_close(nerr, 0.0_dp, 1.0e-12_dp, 'ordered NiAs: unfolding weights sum to 9 at every k')

    do kap_c = 1, 2
      do ll = 0, 2
        pps(kap_c)%enu(ll) = nr%comp(kap_c)%pp(ll, 1)%enu
        pps(kap_c)%c(ll) = nr%comp(kap_c)%pp(ll, 1)%c
        pps(kap_c)%sqrdel(ll) = nr%comp(kap_c)%pp(ll, 1)%sqrdel
        pps(kap_c)%p(ll) = nr%comp(kap_c)%pp(ll, 1)%p
        pps(kap_c)%gamma(ll) = nr%comp(kap_c)%pp(ll, 1)%gamma
      end do
      pps(kap_c)%alpha_rep = pps(kap_c)%gamma
    end do
    al = standard_alpha()
    nst = ni%ci%st
    call space_group(nst, nops)
    call spectral_origin_classes(nst, nops, bcl, ob, osg, ow)
    allocate (ntk(3))
    do uo = 1, 3
      call crystal_graph(nst, nops, ni%ci%rcl, ni%ci%nlevels, ncg(uo), nmsg, start=ob(uo), rmin=ni%ci%spectral%rmax)
      if (len(nmsg) > 0) error stop nmsg
      ntk(uo) = count([(4.0_dp * sum(ncg(uo)%disp(:, ii)**2) <= 4.0_dp * 1.1_dp**2 + 1.0e-9_dp, ii=1, ncg(uo)%n)])
    end do
    ntmax = maxval(ntk)
    allocate (td(3, ntmax, 3), tbase(ntmax, 3), tsign(ntmax, 3), tmir(ntmax, 3), gdg(9, ntmax, 7, 3), gdg0(9, ntmax, 7, 3))
    td = 0.0_dp; tbase = 1; tsign = 1; tmir = 0; gdg = (0.0_dp, 0.0_dp)
    call structure_blocks(nst, ncg(1), max(ni%ci%rcl, ni%ci%rinv), al, nblk, nblk0, ntyp)
    nmdl%nl = 9; nmdl%nsp = 2; nmdl%ncls = 2
    allocate (nmdl%nsc(2), nmdl%spc(1, 2), nmdl%conc(1, 2))
    nmdl%nsc = 1; nmdl%spc(1, :) = [1, 2]; nmdl%conc = 1.0_dp
    nmdl%block = nblk; nmdl%block0 = nblk0
    dm1 = 0.0_dp; do ll = 1, 9; dm1(ll, ll, 1) = 1.0_dp; end do

    ! (1) dense reference on each of the 3 unreduced origin clusters
    do uo = 1, 3
      kap_c = nst%label(ob(uo))
      call edc_setup(nmdl, pps, al, kap_c, nen(kap_c))
      allocate (ncls_(ncg(uo)%n), nfl(ncg(uo)%n), nperm(ncg(uo)%n, 1))
      ncls_ = nst%label(ncg(uo)%basis); nfl = .false.
      nperm(:, 1) = [(ii, ii=1, ncg(uo)%n)]
      nts = pack([(ii, ii=1, ncg(uo)%n)], &
                 [(4.0_dp * sum(ncg(uo)%disp(:, ii)**2) <= 4.0_dp * 1.1_dp**2 + 1.0e-9_dp, ii=1, ncg(uo)%n)])
      nto = size(nts)
      td(:, 1:nto, uo) = ncg(uo)%disp(:, nts)
      tbase(1:nto, uo) = bcl(ncg(uo)%basis(nts))
      tsign(1:nto, uo) = osg(ncg(uo)%basis(nts))
      call build_augcls_operator(ncg(uo)%g, nmdl, ncls_, nperm, dm1, 1, ni%ci%nlevels, nop, maxrank=ni%ci%maxrank, &
                                 flippable=nfl, tsite=nts, tnode=ntn, tgop=ntg, stype=ntyp(ncg(uo)%basis))
      call augcls_set_params(nop, nmdl); call augcls_freeze_site(nop, 1, kap_c)
      call dense_h_multi(nop, nen(kap_c), ncls_, nhm)
      ndo = nop%dim()
      allocate (nkm(ndo, ndo), nrhs(ndo, 9), nipiv(ndo))
      do jj = 1, 7
        zc = cmplx(nr%spec_e(jj), ni%ci%spectral%eta, dp)
        nkm = -nhm
        do ii = 1, ndo
          ks_s = ncls_(nop%site((ii - 1) / 9 + 1))
          ll = mod(ii - 1, 9) + 1
          nkm(ii, ii) = nkm(ii, ii) + nen(kap_c)%dc(ll) * (zc - nen(kap_c)%enu(ll, ks_s)) / nen(kap_c)%del(ll, ks_s)
        end do
        nrhs = (0.0_dp, 0.0_dp)
        do ll = 1, 9; nrhs(ll, ll) = (1.0_dp, 0.0_dp); end do
        call zgesv(ndo, 9, nkm, ndo, nipiv, nrhs, ndo, inf)
        if (inf /= 0) error stop 'test_spectral_scf: NiAs cluster zgesv'
        do ll = 1, 9
          do ii = 1, nto
            ks_s = ncls_(nts(ii))
            gdg(ll, ii, jj, uo) = nrhs(9 * (ntn(1, ii) - 1) + ll, ll) * &
                                  sqrt(nen(kap_c)%dc(ll) / nen(kap_c)%del(ll, ks_s))
          end do
        end do
      end do
      call augcls_free(nop)
      deallocate (ncls_, nfl, nperm, nts, ntn, ntg, nhm, nkm, nrhs, nipiv)
    end do
    do uo = 1, 3
      do ii = 1, ntk(uo)
        u2o = tbase(ii, uo); tmir(ii, uo) = ii
        do jj = 1, ntk(u2o)
          if (tbase(jj, u2o) == uo .and. all(abs(td(:, jj, u2o) + tsign(ii, uo) * td(:, ii, uo)) < 1.0e-9_dp)) then
            tmir(ii, uo) = jj; exit
          end if
        end do
      end do
    end do
    gdg0 = gdg
    do uo = 1, 3
      do ii = 1, ntk(uo)
        if (tmir(ii, uo) /= ii .or. tbase(ii, uo) /= uo) then
          gdg(:, ii, :, uo) = 0.5_dp * (gdg0(:, ii, :, uo) + gdg0(:, tmir(ii, uo), :, tbase(ii, uo)))
        end if
      end do
    end do
    naref = 0.0_dp; allocate (vv(ntmax))
    do uo = 1, 3
      nto = ntk(uo)
      do jj = 1, 7
        do ll = 1, 9
          do ii = 1, nto
            vv(ii) = -aimag(gdg(ll, ii, jj, uo)) / pi * spectral_window('ball', norm2(td(:, ii, uo)), 1.1_dp)
          end do
          do kk = 1, 4
            naref(jj, kk) = naref(jj, kk) + ow(uo) * spectral_unfold(td(:, 1:nto, uo), vv(1:nto), nr%spec_k(:, kk))
          end do
        end do
      end do
    end do
    namax = maxval(abs(naref))
    nerr = maxval(abs(sum(nr%spec_a(:, :, :, 1), dim=3) - naref)) / namax
    print '(a,es10.2)', '  NiAs max |A(class_scf) - A(dense 3 origins)| / max A = ', nerr
    call check_close(nerr, 0.0_dp, 1.0e-6_dp, 'ordered NiAs: class_scf = dense inverse on the 3 origin clusters')

    ! (2) exact 3x3x2 periodic supercell (72 sites): unscaled H^(2) = E_nu + h - h o h vs weighted H^(2)(k) Lorentzians
    uo = 0
    do i1 = 0, l1 - 1; do i2 = 0, l2 - 1; do i3 = 0, l3 - 1; do b1 = 1, 4
      uo = uo + 1; uidx(i1, i2, i3, b1) = uo
      rpos(:, uo) = matmul(nst%avec, real([i1, i2, i3], dp)) + nst%tau(:, b1)
      kap_c = nst%label(b1)
      call edc_setup(nmdl, pps, al, kap_c, nen(kap_c))
      cv(9 * (uo - 1) + 1:9 * uo) = nmdl%c(:, kap_c)
      dv(9 * (uo - 1) + 1:9 * uo) = nmdl%d(:, kap_c)
      ov(9 * (uo - 1) + 1:9 * uo) = nmdl%o(:, kap_c)
      ev0(9 * (uo - 1) + 1:9 * uo) = nen(kap_c)%enu(:, kap_c)
    end do; end do; end do; end do
    uorig = [uidx(0, 0, 0, 1), uidx(0, 0, 0, 2), uidx(0, 0, 0, 3)]
    allocate (smat(nds, nds), hmat(nds, nds), h2mat(nds, nds))
    smat = 0.0_dp
    do i1 = 0, l1 - 1; do i2 = 0, l2 - 1; do i3 = 0, l3 - 1; do b1 = 1, 4
      u1 = uidx(i1, i2, i3, b1)
      smat(9 * (u1 - 1) + 1:9 * u1, 9 * (u1 - 1) + 1:9 * u1) = nblk0(:, :, ntyp(b1))
      do ks = 1, ncg(1)%nslot
        b2 = ncg(1)%btgt(ks, b1)
        if (b2 == 0) cycle
        j1 = modulo(i1 + ncg(1)%bdn(1, ks, b1), l1)
        j2 = modulo(i2 + ncg(1)%bdn(2, ks, b1), l2)
        j3 = modulo(i3 + ncg(1)%bdn(3, ks, b1), l3)
        u2 = uidx(j1, j2, j3, b2)
        smat(9 * (u1 - 1) + 1:9 * u1, 9 * (u2 - 1) + 1:9 * u2) = &
          smat(9 * (u1 - 1) + 1:9 * u1, 9 * (u2 - 1) + 1:9 * u2) + nblk(:, :, ks, ntyp(b1))
      end do
    end do; end do; end do; end do
    do jj = 1, nds
      do ii = 1, nds; hmat(ii, jj) = dv(ii) * smat(ii, jj) * dv(jj); end do
      hmat(jj, jj) = hmat(jj, jj) + cv(jj)
    end do
    do jj = 1, nds
      do ii = 1, nds; h2mat(ii, jj) = hmat(ii, jj) - sum(hmat(ii, :) * ov * hmat(:, jj)); end do
      h2mat(jj, jj) = h2mat(jj, jj) + ev0(jj)
    end do
    allocate (nkm(nds, nds), nrhs(nds, 27), nipiv(nds))
    nafull = 0.0_dp; nafull_1sub = 0.0_dp; naband = 0.0_dp
    do jj = 1, 7
      zc = cmplx(nr%spec_e(jj), ni%ci%spectral%eta, dp)
      nkm = -h2mat
      do ii = 1, nds; nkm(ii, ii) = nkm(ii, ii) + zc; end do
      nrhs = (0.0_dp, 0.0_dp)
      do uo = 1, 3
        do ll = 1, 9
          nrhs(9 * (uorig(uo) - 1) + ll, 9 * (uo - 1) + ll) = (1.0_dp, 0.0_dp)
        end do
      end do
      call zgesv(nds, 27, nkm, nds, nipiv, nrhs, nds, inf)
      if (inf /= 0) error stop 'test_spectral_scf: NiAs supercell zgesv'
      do uo = 1, 3
        do ii = 1, nsup; dsup(:, ii) = rpos(:, ii) - rpos(:, uorig(uo)); end do
        do ll = 1, 9
          do ii = 1, nsup
            vsup(ii) = -aimag(nrhs(9 * (ii - 1) + ll, 9 * (uo - 1) + ll)) / pi
          end do
          do kk = 1, 4
            nafull(jj, kk) = nafull(jj, kk) + ow(uo) * spectral_unfold(dsup, vsup, nr%spec_k(:, kk))
            if (uo == 1 .or. uo == 3) then
              nafull_1sub(jj, kk) = nafull_1sub(jj, kk) + 0.5_dp * spectral_unfold(dsup, vsup, nr%spec_k(:, kk))
            end if
          end do
        end do
      end do
      do kk = 1, 4
        naband(jj, kk) = sum(nr%spec_bwgt(:, kk, 1) * (-aimag(1.0_dp / (zc - nr%spec_bands(:, kk, 1))) / pi))
      end do
    end do
    nerr = maxval(abs(naband - nafull)) / maxval(nafull)
    nerr1 = maxval(abs(naband(:, 2) - nafull_1sub(:, 2))) / maxval(nafull)
    print '(a,es10.2,a,es10.2)', '  NiAs max |A(bands) - A(3x3x2 supercell)| / max A = ', nerr, &
      ',  1-origin-per-sublattice error at L3 = ', nerr1
    call check(minval(nafull) > 0.0_dp, 'ordered NiAs periodic supercell: unfolded A(k, E) > 0')
    call check_close(nerr, 0.0_dp, 1.0e-10_dp, 'ordered NiAs: weighted H^(2)(k) Lorentzians = exact unfolded A of supercell')
    call check(nerr1 > 1.0e-2_dp, 'ordered NiAs: single cation origin fails at L3 where k breaks screw and glide')
  end subroutine nias_ordered_check

  !> dense matrix of h' - h' q h' on a multi-component ordered node set (energy-dependent diagonal ep(z) added per energy)
  subroutine dense_h_multi(op, ed, ncls_, hmat)
    type(augcls_op), intent(inout) :: op
    type(edn_t), intent(in) :: ed
    integer, intent(in) :: ncls_(:)
    real(dp), allocatable, intent(out) :: hmat(:, :)
    real(dp), allocatable :: xv(:), yv(:)
    integer :: jj
    if (ed%ns < 0 .or. size(ncls_) < 0) return
    allocate (hmat(op%dim(), op%dim()), xv(op%dim()), yv(op%dim()))
    do jj = 1, op%dim()
      xv = 0.0_dp; xv(jj) = 1.0_dp
      call op%apply(xv, yv)
      hmat(:, jj) = yv
    end do
  end subroutine dense_h_multi

  !> dense matrix of z - H on the node set without z: the unreduced operator plus E_nu (one species: e_s(z) = E_nu, the term edc
  !> applies separately)
  subroutine dense_h(op, ed, hmat)
    type(augcls_op), intent(inout) :: op
    type(edn_t), intent(in) :: ed
    real(dp), allocatable, intent(out) :: hmat(:, :)
    real(dp), allocatable :: xv(:), yv(:)
    integer :: jj, ii, ll
    allocate (hmat(op%dim(), op%dim()), xv(op%dim()), yv(op%dim()))
    do jj = 1, op%dim()
      xv = 0.0_dp; xv(jj) = 1.0_dp
      call op%apply(xv, yv)
      hmat(:, jj) = yv
    end do
    do ii = 1, op%nn
      do ll = 1, 9; hmat(9 * (ii - 1) + ll, 9 * (ii - 1) + ll) = hmat(9 * (ii - 1) + ll, 9 * (ii - 1) + ll) + ed%enu(ll, 1); end do
    end do
  end subroutine dense_h
end program test_spectral_scf
