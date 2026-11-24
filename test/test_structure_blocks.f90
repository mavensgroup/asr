!> Structure constants per basis site (asr_lmto_alloy::structure_blocks) and the class operator with per-type tables
!> (asr_augcls, asr_augblk), step 3 of notes/DESIGN_structures.md.
!>  (A) sc, bcc, fcc and B2 (both basis sites) through the structure: the blocks of build_lmto_alloy, bit for bit.
!>  (B) bcc and B2 class operators on the site cluster of the structure, with its site operations and per-type tables, against
!>      today's operators on a periodic box of 10 cubes (large enough not to wrap at 2 levels): today's reduction group and
!>      orbital matrices bit for bit, the same number of nodes and the same moments for every start orbital.  A box of 4 cubes
!>      wraps around at 2 levels and gives another node set.
!>  (S) every space-group operation g of hcp and NiAs maps the blocks of basis site b onto those of g(b):
!>      block0(g b) = D block0(b) D^T, block(g k, g b) = D block(k, b) D^T; the reverse bond holds the transpose.  hcp: the two
!>      sites are related by inversion, so their on-site p-d blocks are nonzero and opposite, the d-d blocks equal.
!>  (G) gate G2 on small hcp and NiAs maps with disorder on every sublattice: the reduced operator (monomial reduction group)
!>      gives the moments of the unreduced one for every start orbital; partners of an orbital group give equal moments; the
!>      unreduced operator is symmetric; the block kernel (asr_augblk) gives h x of asr_augcls with distinct site types.
program test_structure_blocks
  use asr_kinds, only: dp
  use asr_crystal, only: structure_t
  use asr_prototypes, only: prototype_structure
  use asr_crystal_symmetry, only: symop_t, space_group, site_ops, is_monomial, orbital_matrix, orbital_groups, reduction_ops
  use asr_crystal_graph
  use asr_augspace, only: graph_t, graph_distances
  use asr_augcls
  use asr_augblk, only: blk_work_t, blk_work_init, blk_hp
  use asr_lmto_alloy, only: build_lmto_alloy, orb_alloy_t, structure_blocks
  use asr_lmto_ham, only: potpar_t
  use asr_symmetry, only: class_ops, orbital_rep
  use asr_testing
  implicit none
  real(dp), parameter :: alpha(0:2) = [0.3485_dp, 0.05303_dp, 0.010714_dp]
  type(potpar_t) :: pp(2)
  integer :: l
  print '(a)', 'test_structure_blocks'
  do l = 0, 2
    pp(1)%enu(l) = 0.1_dp * l; pp(1)%c(l) = 0.3_dp - 0.15_dp * l; pp(1)%sqrdel(l) = 0.30_dp - 0.07_dp * l
    pp(1)%p(l) = 0.9_dp - 0.1_dp * l; pp(1)%gamma(l) = 0.25_dp - 0.1_dp * l
    pp(2)%enu(l) = 0.05_dp * l; pp(2)%c(l) = 0.45_dp - 0.1_dp * l; pp(2)%sqrdel(l) = 0.27_dp - 0.06_dp * l
    pp(2)%p(l) = 0.8_dp - 0.1_dp * l; pp(2)%gamma(l) = 0.2_dp - 0.08_dp * l
  end do
  pp(1)%alpha_rep = pp(1)%gamma; pp(2)%alpha_rep = pp(2)%gamma
  call cubic_blocks('sc', 'sc', 3, 1.2_dp, 1.2_dp)
  call cubic_blocks('bcc', 'bcc', 3, 1.2_dp, 1.75_dp)
  call cubic_blocks('fcc', 'fcc', 3, 1.2_dp, 1.5_dp)
  call cubic_blocks('b2', 'bcc', 3, 1.2_dp, 1.75_dp)
  call cubic_operator('bcc')
  call cubic_operator('b2')
  call sym_blocks('hcp', 0.0_dp, 1.05_dp, 1.6_dp)
  call sym_blocks('hcp', 1.7_dp, 1.05_dp, 1.6_dp)
  call sym_blocks('nias', 1.617_dp, 1.0_dp, 1.5_dp)
  call g2_case('hcp', 0.0_dp, 1.05_dp, 1)
  call g2_case('nias', 1.617_dp, 0.9_dp, 1)
  call g2_case('nias', 1.617_dp, 0.9_dp, 3)
  call report()
contains

  subroutine proto(name, ca, s, ops)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: ca
    type(structure_t), intent(out) :: s
    type(symop_t), allocatable, intent(out) :: ops(:)
    character(len=:), allocatable :: msg
    if (ca > 0.0_dp) then
      call prototype_structure(name, s, msg, c_over_a=ca)
    else
      call prototype_structure(name, s, msg)
    end if
    if (len(msg) > 0) error stop msg
    call space_group(s, ops)
  end subroutine proto

  subroutine graph(s, ops, rcl, levels, cg, start)
    type(structure_t), intent(in) :: s
    type(symop_t), intent(in) :: ops(:)
    real(dp), intent(in) :: rcl
    integer, intent(in) :: levels
    type(crystal_graph_t), intent(out) :: cg
    integer, intent(in), optional :: start
    character(len=:), allocatable :: msg
    call crystal_graph(s, ops, rcl, levels, cg, msg, start)
    if (len(msg) > 0) error stop msg
  end subroutine graph

  ! ---------------------------------------------------------------- (A)
  subroutine cubic_blocks(name, kind, l, rcl, rinv)
    character(len=*), intent(in) :: name, kind
    integer, intent(in) :: l
    real(dp), intent(in) :: rcl, rinv
    type(structure_t) :: s
    type(symop_t), allocatable :: ops(:)
    type(crystal_graph_t) :: cg
    type(graph_t) :: g
    type(orb_alloy_t) :: m0
    integer, allocatable :: coord(:, :)
    real(dp), allocatable :: blk(:, :, :, :), blk0(:, :, :)
    integer, allocatable :: typ(:)
    integer :: b
    real(dp) :: e1, e0
    call proto(name, 0.0_dp, s, ops)
    call graph(s, ops, rcl, 1, cg)
    call structure_blocks(s, cg, rinv, alpha, blk, blk0, typ)
    call build_lmto_alloy(kind, l, rcl, 0.7_dp, pp(1), pp(2), 0.5_dp, g, coord, m0, rinv, alpha)
    call check(size(blk, 3) == size(m0%block, 3) .and. size(blk, 4) == 1 .and. all(typ == 1), &
               name//': one table (B2: both sites have the bcc geometry), slots as today')
    if (size(blk, 3) /= size(m0%block, 3)) return
    do b = 1, size(blk, 4)
      e1 = maxval(abs(blk(:, :, :, b) - m0%block)); e0 = maxval(abs(blk0(:, :, b) - m0%block0))
      print '(2x,a,a,i0,a,2es10.2)', name, ' type ', b, ': max |block - today|, |block0 - today| ', e1, e0
      call check(e1 == 0.0_dp .and. e0 == 0.0_dp, name//': blocks of the structure equal build_lmto_alloy bit for bit')
    end do
  end subroutine cubic_blocks

  ! ---------------------------------------------------------------- (B)
  subroutine cubic_operator(name)
    character(len=*), intent(in) :: name
    real(dp), parameter :: rcl = 1.2_dp, rinv = 1.75_dp
    integer, parameter :: nlev = 2, lbig = 10, lsmall = 4, nmom = 8
    type(structure_t) :: s
    type(symop_t), allocatable :: ops(:)
    type(crystal_graph_t) :: cg
    type(class_alloy_t) :: mn
    type(augcls_op) :: opo, opn, ops4
    integer, allocatable :: perm(:, :), mats(:, :, :), gperm(:, :), gidx(:), red(:)
    real(dp), allocatable :: d(:, :, :), rot(:, :, :), dg(:, :, :), blk(:, :, :, :), blk0(:, :, :)
    integer, allocatable :: typ(:)
    logical, allocatable :: fl(:), mono(:)
    real(dp) :: mo_(0:nmom), mn_(0:nmom), sc(0:nmom), err
    integer :: k, grp(9), ngrp, rep(9), wgt(9), nn4
    call proto(name, 0.0_dp, s, ops)
    call graph(s, ops, rcl, nlev, cg)
    call structure_blocks(s, cg, rinv, alpha, blk, blk0, typ)
    ! new: site operations of the cluster, orbital groups and the reduction group from the structure
    call graph_site_ops(cg, s, ops, gperm, rot, gidx)
    allocate (dg(9, 9, size(gperm, 2)), mono(size(gperm, 2)))
    do k = 1, size(gperm, 2)
      call orbital_matrix(rot(:, :, k), dg(:, :, k)); mono(k) = is_monomial(rot(:, :, k))
    end do
    call orbital_groups(dg, mono, grp, ngrp, rep, wgt)
    call reduction_ops(dg, mono, rep(1:ngrp), red)
    call model(mn, s%nb)
    mn%block = blk; mn%block0 = blk0
    allocate (fl(cg%n)); fl = .true.; fl(1) = .false.
    call build_augcls_operator(cg%g, mn, s%label(cg%basis), gperm(:, red), dg(:, :, red), 1, nlev, opn, maxrank=1, flippable=fl, &
                               stype=typ(cg%basis))
    ! today: a periodic box of lbig cubes, classes from the half-lattice coordinates, class_ops
    call torus_operator(s%nb, lbig, nlev, rcl, rinv, opo, perm, mats, d)
    call check(size(red) == size(perm, 2), name//': reduction group of the structure has today''s order')
    if (size(red) /= size(perm, 2)) return
    call check(all(nint(rot(:, :, red)) == mats), name//': reduction group rotations equal class_ops, in order')
    call check(all(dg(:, :, red) == d), name//': orbital matrices equal orbital_rep, bit for bit')
    call check(opn%ntype == 1, name//': one table for every site (bcc geometry)')
    call check(opn%nn == opo%nn, name//': the cluster gives the node set of a box that does not wrap')
    err = 0.0_dp
    do k = 1, ngrp
      call moments(opo, d, rep(k), mo_, sc)
      call moments(opn, dg(:, :, red), rep(k), mn_)
      err = max(err, maxval(abs(mn_ - mo_) / sc))
    end do
    call check(err < 1.0e-13_dp, name//': moments of every start orbital equal those of the box (1e-13)')
    call torus_operator(s%nb, lsmall, nlev, rcl, rinv, ops4)
    nn4 = ops4%nn
    print '(2x,a,a,i0,a,i0,a,i0,a,i0,a,es9.2)', name, ': cluster sites ', cg%n, ', nodes ', opn%nn, ' (box of ', lbig, &
      ' cubes ', opo%nn, '), moments error ', err
    print '(2x,a,a,i0,a,i0)', name, ': box of ', lsmall, ' cubes wraps around at 2 levels, nodes ', nn4
    call check(nn4 /= opn%nn, name//': a box of 4 cubes wraps around at 2 levels (another node set)')
    call augcls_free(opo); call augcls_free(opn); call augcls_free(ops4)
  end subroutine cubic_operator

  !> Today's operator: a periodic bcc box of l cubes (build_lmto_alloy), classes from the half-lattice coordinates (nb = 2: B2),
  !> class_ops and orbital_rep, the model of cubic_operator, nlev levels.
  subroutine torus_operator(nb, l, nlev, rcl, rinv, op, pm, mt, dd)
    integer, intent(in) :: nb, l, nlev
    real(dp), intent(in) :: rcl, rinv
    type(augcls_op), intent(out) :: op
    integer, allocatable, intent(out), optional :: pm(:, :), mt(:, :, :)
    real(dp), allocatable, intent(out), optional :: dd(:, :, :)
    type(graph_t) :: g
    type(orb_alloy_t) :: m0
    type(class_alloy_t) :: mo
    integer, allocatable :: coord(:, :), cls(:), p_(:, :), m_(:, :, :)
    real(dp), allocatable :: d_(:, :, :)
    logical, allocatable :: f_(:)
    integer :: i
    call build_lmto_alloy('bcc', l, rcl, 0.7_dp, pp(1), pp(2), 0.5_dp, g, coord, m0, rinv, alpha)
    allocate (cls(g%n))
    if (nb == 1) then
      cls = 1
    else
      do i = 1, g%n; cls(i) = 1 + modulo(coord(1, i), 2); end do
    end if
    call class_ops(coord, 2 * l, cls, p_, m_)
    call orbital_rep(m_, d_)
    call model(mo, nb)
    call augcls_set_blocks(mo, m0%block, m0%block0)
    allocate (f_(g%n)); f_ = .true.; f_(1) = .false.
    call build_augcls_operator(g, mo, cls, p_, d_, 1, nlev, op, maxrank=1, flippable=f_)
    if (present(pm)) pm = p_
    if (present(mt)) mt = m_
    if (present(dd)) dd = d_
  end subroutine torus_operator

  !> Two species on class 1 (x = 0.35), one species (the second) on class 2 if there are two classes; l-dependent parameters.
  subroutine model(m, ncls)
    type(class_alloy_t), intent(out) :: m
    integer, intent(in) :: ncls
    m%nl = 9; m%nsp = 2; m%ncls = ncls
    allocate (m%nsc(ncls), m%spc(2, ncls), m%conc(2, ncls)); m%spc = 0; m%conc = 0.0_dp
    m%nsc(1) = 2; m%spc(:, 1) = [1, 2]; m%conc(:, 1) = [0.35_dp, 0.65_dp]
    if (ncls == 2) then; m%nsc(2) = 1; m%spc(1, 2) = 2; m%conc(1, 2) = 1.0_dp; end if
    call species(m, 2)
  end subroutine model

  !> l-dependent species parameters c, d, o, e (the same for every m of one l, as in the LMTO model)
  subroutine species(m, nsp)
    type(class_alloy_t), intent(inout) :: m
    integer, intent(in) :: nsp
    integer, parameter :: ol(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2]
    integer :: sp, i
    allocate (m%c(9, nsp), m%d(9, nsp), m%o(9, nsp), m%e(9, nsp))
    do sp = 1, nsp
      do i = 1, 9
        m%c(i, sp) = 0.30_dp - 0.12_dp * ol(i) + 0.07_dp * sp
        m%d(i, sp) = 0.32_dp - 0.06_dp * ol(i) - 0.02_dp * sp
        m%o(i, sp) = 0.04_dp + 0.01_dp * ol(i) * sp
        m%e(i, sp) = 0.05_dp * ol(i) - 0.03_dp * sp
      end do
    end do
  end subroutine species

  ! ---------------------------------------------------------------- (S)
  subroutine sym_blocks(name, ca, rcl, rinv)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: ca, rcl, rinv
    type(structure_t) :: s
    type(symop_t), allocatable :: ops(:)
    type(crystal_graph_t) :: cg
    real(dp), allocatable :: blk(:, :, :, :), blk0(:, :, :)
    integer, allocatable :: typ(:)
    real(dp) :: dm(9, 9), e0, e1, et, par(9, 9), big, pd
    integer :: g, b, b2, k, k2, kr
    integer, parameter :: parity(9) = [1, -1, -1, -1, 1, 1, 1, 1, 1]
    character(len=24) :: tag
    write (tag, '(a,a,f5.3)') name, ' c/a ', ca
    call proto(name, ca, s, ops)
    call graph(s, ops, rcl, 1, cg)
    call structure_blocks(s, cg, rinv, alpha, blk, blk0, typ)
    call check(all(typ == [(b, b=1, s%nb)]), trim(tag)//': no two basis sites related by a translation (one table each)')
    big = maxval(abs(blk))
    e0 = 0.0_dp; e1 = 0.0_dp; et = 0.0_dp
    do g = 1, size(ops)
      call orbital_matrix(ops(g)%rot, dm)
      do b = 1, s%nb
        b2 = ops(g)%perm(b)
        e0 = max(e0, maxval(abs(blk0(:, :, b2) - matmul(dm, matmul(blk0(:, :, b), transpose(dm))))))
        do k = 1, cg%nslot
          k2 = slot_of(cg, matmul(ops(g)%rot, cg%bvec(:, k)))
          if (k2 == 0) error stop 'sym_blocks: rotated bond vector is not a slot'
          if ((cg%btgt(k, b) == 0) .neqv. (cg%btgt(k2, b2) == 0)) error stop 'sym_blocks: slot sets do not match'
          e1 = max(e1, maxval(abs(blk(:, :, k2, b2) - matmul(dm, matmul(blk(:, :, k, b), transpose(dm))))))
        end do
      end do
    end do
    do b = 1, s%nb
      do k = 1, cg%nslot
        if (cg%btgt(k, b) == 0) cycle
        kr = slot_of(cg, -cg%bvec(:, k))
        et = max(et, maxval(abs(blk(:, :, kr, cg%btgt(k, b)) - transpose(blk(:, :, k, b)))))
      end do
    end do
    print '(2x,a,a,3es10.2,a,es9.2)', trim(tag), ': symmetry errors on-site, bonds, reverse bonds ', e0, e1, et, ', max |S| ', big
    call check(e0 < 1.0e-12_dp * big, trim(tag)//': space group maps the on-site blocks of b onto those of g(b)')
    call check(e1 < 1.0e-12_dp * big, trim(tag)//': space group maps the bond blocks of b onto those of g(b)')
    call check(et == 0.0_dp, trim(tag)//': the reverse bond holds the transpose exactly')
    if (name == 'hcp') then
      do k = 1, 9; par(:, k) = 0.0_dp; par(k, k) = real(parity(k), dp); end do
      pd = maxval(abs(blk0(2:4, 5:9, 1)))
      print '(2x,a,a,es10.2)', trim(tag), ': max on-site p-d block of site A ', pd
      call check(pd > 1.0e-4_dp, trim(tag)//': on-site p-d mixing (no inversion centre on the site)')
      call check(maxval(abs(blk0(:, :, 2) - matmul(par, matmul(blk0(:, :, 1), par)))) < 1.0e-12_dp * big, &
                 trim(tag)//': block0(B) = P block0(A) P (inversion)')
      call check(maxval(abs(blk0(2:4, 5:9, 2) + blk0(2:4, 5:9, 1))) < 1.0e-12_dp * big, trim(tag)//': on-site p-d blocks opposite')
      call check(maxval(abs(blk0(5:9, 5:9, 2) - blk0(5:9, 5:9, 1))) < 1.0e-12_dp * big, trim(tag)//': on-site d-d blocks equal')
    end if
  end subroutine sym_blocks

  integer function slot_of(cg, v)
    type(crystal_graph_t), intent(in) :: cg
    real(dp), intent(in) :: v(3)
    integer :: k
    slot_of = 0
    do k = 1, cg%nslot
      if (maxval(abs(cg%bvec(:, k) - v)) < 1.0e-9_dp) then; slot_of = k; return; end if
    end do
  end function slot_of

  ! ---------------------------------------------------------------- (G)
  subroutine g2_case(name, ca, rcl, start)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: ca, rcl
    integer, intent(in) :: start
    integer, parameter :: nmom = 8
    type(structure_t) :: s
    type(symop_t), allocatable :: ops(:)
    type(crystal_graph_t) :: cg
    type(class_alloy_t) :: m
    type(augcls_op) :: opr, opu
    integer, allocatable :: perm(:, :), gidx(:), red(:), cls(:), sidx(:), perm1(:, :)
    real(dp), allocatable :: rot(:, :, :), dg(:, :, :), dcr(:, :, :), blk(:, :, :, :), blk0(:, :, :), x(:), y(:), hy(:), hx(:)
    integer, allocatable :: typ(:)
    logical, allocatable :: fl(:), mono(:), mcr(:)
    real(dp) :: d1(9, 9, 1), mr(0:nmom), mu(0:nmom), sc(0:nmom), err, errp, asym, eb(2)
    real(dp), allocatable :: mall(:, :)
    integer :: k, i, l, grp(9), ngrp, rep(9), wgt(9), cgrp(9), cngrp, crep(9), cwgt(9)
    character(len=24) :: tag
    write (tag, '(a,a,i0)') name, ' start ', start
    call proto(name, ca, s, ops)
    call graph(s, ops, rcl, 2, cg, start)
    call structure_blocks(s, cg, 1.5_dp, alpha, blk, blk0, typ)
    ! orbital groups from the crystal site group, reduction group from the graph operations
    call site_ops(ops, start, sidx)
    allocate (dcr(9, 9, size(sidx)), mcr(size(sidx)))
    do k = 1, size(sidx)
      call orbital_matrix(ops(sidx(k))%rot, dcr(:, :, k)); mcr(k) = is_monomial(ops(sidx(k))%rot)
    end do
    call orbital_groups(dcr, mcr, cgrp, cngrp, crep, cwgt)
    call graph_site_ops(cg, s, ops, perm, rot, gidx)
    allocate (dg(9, 9, size(perm, 2)), mono(size(perm, 2)))
    do k = 1, size(perm, 2)
      call orbital_matrix(rot(:, :, k), dg(:, :, k)); mono(k) = is_monomial(rot(:, :, k))
    end do
    call orbital_groups(dg, mono, grp, ngrp, rep, wgt)
    call check(ngrp == cngrp .and. all(grp == cgrp) .and. all(rep == crep), &
               trim(tag)//': graph operations give the crystal orbital groups')
    call reduction_ops(dg, mono, rep(1:ngrp), red)
    call check(red(1) == 1 .and. size(red) == 4, trim(tag)//': reduction group of order 4, identity first')
    ! disorder on every sublattice: class = label, two species each
    allocate (cls(cg%n)); cls = s%label(cg%basis)
    m%nl = 9; m%ncls = s%nlabel; m%nsp = 2 * s%nlabel
    allocate (m%nsc(m%ncls), m%spc(2, m%ncls), m%conc(2, m%ncls))
    do k = 1, m%ncls
      m%nsc(k) = 2; m%spc(:, k) = [2 * k - 1, 2 * k]; m%conc(:, k) = [0.3_dp + 0.1_dp * k, 0.7_dp - 0.1_dp * k]
    end do
    call species(m, m%nsp)
    m%block = blk; m%block0 = blk0
    allocate (fl(cg%n)); fl = .true.; fl(1) = .false.
    allocate (perm1(cg%n, 1)); perm1(:, 1) = [(i, i=1, cg%n)]
    d1 = 0.0_dp; do l = 1, 9; d1(l, l, 1) = 1.0_dp; end do
    call build_augcls_operator(cg%g, m, cls, perm(:, red), dg(:, :, red), 1, 2, opr, maxrank=1, flippable=fl, stype=typ(cg%basis))
    call build_augcls_operator(cg%g, m, cls, perm1, d1, 1, 2, opu, maxrank=1, flippable=fl, stype=typ(cg%basis))
    print '(2x,a,a,i0,a,i0,a,i0,a,i0)', trim(tag), ': sites ', cg%n, ', slots ', cg%nslot, ', nodes reduced ', opr%nn, &
      ', unreduced ', opu%nn
    ! moments of every start orbital, reduced against unreduced
    allocate (mall(0:nmom, 9))
    err = 0.0_dp
    do l = 1, 9
      call moments(opu, d1, l, mu, sc); mall(:, l) = mu
      if (rep(grp(l)) /= l) cycle
      call moments(opr, dg(:, :, red), l, mr)
      err = max(err, maxval(abs(mr - mu) / sc))
    end do
    errp = 0.0_dp
    do l = 1, 9
      call moments(opu, d1, rep(grp(l)), mu, sc)
      errp = max(errp, maxval(abs(mall(:, l) - mu) / sc))
    end do
    print '(2x,a,a,2es10.2)', trim(tag), ': moments reduced - unreduced, partner - representative ', err, errp
    call check(err < 1.0e-12_dp, trim(tag)//': G2, reduced moments equal unreduced moments (1e-12)')
    call check(errp < 1.0e-12_dp, trim(tag)//': partners of an orbital group have the moments of its representative')
    ! symmetry of the unreduced operator
    allocate (x(opu%dim()), y(opu%dim()), hx(opu%dim()), hy(opu%dim()))
    call lcg(x); call lcg(y)
    call opu%apply(x, hx); call opu%apply(y, hy)
    asym = abs(dot_product(y, hx) - dot_product(hy, x)) / (norm2(hx) * norm2(y))
    call check(asym < 1.0e-13_dp, trim(tag)//': the unreduced operator is symmetric')
    call block_kernel(opu, eb(1)); call block_kernel(opr, eb(2))
    print '(2x,a,a,es9.2,a,2es10.2)', trim(tag), ': asymmetry ', asym, ', block kernel - augcls (unreduced, reduced) ', eb
    call check(maxval(eb) < 1.0e-14_dp, trim(tag)//': asr_augblk h x equals asr_augcls h x with distinct site types')
    call augcls_free(opr); call augcls_free(opu)
  end subroutine g2_case

  !> Relative difference of h x from asr_augblk (a block of three complex vectors) and from asr_augcls (one vector at a time).
  subroutine block_kernel(op, e)
    type(augcls_op), intent(inout) :: op
    real(dp), intent(out) :: e
    type(blk_work_t) :: bw
    complex(dp), allocatable :: zx(:, :), zb(:, :), zc(:)
    real(dp), allocatable :: x(:)
    integer :: k
    allocate (zx(3, op%dim()), zb(3, op%dim()), zc(op%dim()), x(op%dim()))
    do k = 1, 3
      call lcg(x); zx(k, :) = cmplx(x, 0.0_dp, dp); call lcg(x); zx(k, :) = zx(k, :) + cmplx(0.0_dp, x, dp)
    end do
    call blk_work_init(bw, 3, op%dim())
    call blk_hp(op, bw, zx, zb)
    e = 0.0_dp
    do k = 1, 3
      call augcls_hp_c(op, zx(k, :), zc)
      e = max(e, maxval(abs(zb(k, :) - zc)) / maxval(abs(zc)))
    end do
  end subroutine block_kernel

  !> mu(n) = <e_l0 | H^n | e_l0> on the start node; sc(n) = |H^n e_l0| (scale of the error)
  subroutine moments(op, dm, l0, mu, sc)
    type(augcls_op), intent(inout) :: op
    real(dp), intent(in) :: dm(:, :, :)
    integer, intent(in) :: l0
    real(dp), intent(out) :: mu(0:)
    real(dp), intent(out), optional :: sc(0:)
    real(dp), allocatable :: v(:), w(:)
    integer :: n
    call augcls_set_start(op, l0, dm)
    allocate (v(op%dim()), w(op%dim())); v = 0.0_dp; v(l0) = 1.0_dp
    mu(0) = 1.0_dp; if (present(sc)) sc(0) = 1.0_dp
    do n = 1, ubound(mu, 1)
      call op%apply(v, w)
      mu(n) = w(l0); if (present(sc)) sc(n) = norm2(w)
      v = w
    end do
  end subroutine moments

  subroutine lcg(v)
    real(dp), intent(out) :: v(:)
    integer(8), save :: st = 20261124_8
    integer :: k
    do k = 1, size(v)
      st = mod(st * 6364136223846793005_8 + 1442695040888963407_8, 9223372036854775807_8)
      v(k) = real(mod(abs(st), 1000003_8), dp) / 1000003.0_dp - 0.5_dp
    end do
  end subroutine lcg
end program test_structure_blocks
