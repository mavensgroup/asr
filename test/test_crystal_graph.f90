!> Real-space graph of a crystal structure (asr_crystal_graph), step 2 of notes/DESIGN_structures.md: a finite cluster about the
!> start site, sized from the recursion depth.
!>  (C) sc, bcc, fcc: the bond vectors of asr_lattices::periodic_graph_shells in order; the sites are those of
!>      asr_strconst::lattice_cluster for the same radius, in the same order; the 48 site operations in the order of
!>      asr_symmetry::cubic_site_perms.  B2: the bcc graph with the two basis sites alternating.
!>  (H) hcp: 12 nearest neighbours per site at the ideal c/a, 6 + 6 otherwise; NiAs: the cation has 6 anions and then 2
!>      cations as its nearest shells, the anion 6 cations.
!>  (S) every slot of every site: exact labels, neighbour - site = bvec(k) exactly, entries 0 (no such bond) and -1 (outside the
!>      cluster) where they belong; the graph is symmetric.
!>  (R) size: every site within 2 levels - 1 bonds of site 1 has all its neighbours (the closure of asr_augcls never meets -1);
!>      rmin enlarges the cluster.
!>  (G) every site operation maps neighbours onto neighbours with the rotated bond vector and rotates the displacements.
program test_crystal_graph
  use asr_kinds, only: dp
  use asr_crystal, only: structure_t
  use asr_prototypes, only: prototype_structure
  use asr_crystal_symmetry, only: symop_t, space_group, is_monomial
  use asr_crystal_graph
  use asr_augspace, only: graph_t, graph_distances
  use asr_lattices, only: periodic_graph, periodic_graph_shells
  use asr_strconst, only: lattice_cluster
  use asr_symmetry, only: cubic_site_perms
  use asr_testing
  implicit none
  real(dp), parameter :: ca_ideal = 1.632993161855452_dp
  print '(a)', 'test_crystal_graph'
  call cubic_case('sc', 1.2_dp, 3)
  call cubic_case('bcc', 1.2_dp, 2)
  call cubic_case('bcc', 1.45_dp, 2)
  call cubic_case('fcc', 1.2_dp, 2)
  call cubic_case('fcc', 0.75_dp, 3)
  call b2_case(1.2_dp, 2)
  call hcp_case(ca_ideal, [12, 12], 0)
  call hcp_case(1.7_dp, [6, 6], 6)
  call nias_case()
  call rmin_case()
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

  subroutine graph(tag, s, ops, rcl, levels, cg, start, rmin)
    character(len=*), intent(in) :: tag
    type(structure_t), intent(in) :: s
    type(symop_t), intent(in) :: ops(:)
    real(dp), intent(in) :: rcl
    integer, intent(in) :: levels
    type(crystal_graph_t), intent(out) :: cg
    integer, intent(in), optional :: start
    real(dp), intent(in), optional :: rmin
    character(len=:), allocatable :: msg
    call crystal_graph(s, ops, rcl, levels, cg, msg, start, rmin)
    call check(len(msg) == 0, tag//': graph builds')
    if (len(msg) > 0) error stop msg
  end subroutine graph

  subroutine cubic_case(kind, rcl, levels)
    character(len=*), intent(in) :: kind
    real(dp), intent(in) :: rcl
    integer, intent(in) :: levels
    type(structure_t) :: s
    type(symop_t), allocatable :: ops(:)
    type(crystal_graph_t) :: cg
    type(graph_t) :: g
    integer, allocatable :: coord(:, :), bv(:, :), p48(:, :), m48(:, :, :), perm(:, :), idx(:)
    real(dp), allocatable :: rot(:, :, :), vec(:, :)
    character(len=40) :: tag
    write (tag, '(a,a,f4.2,a,i0)') kind, ' rcl=', rcl, ' levels=', levels
    call proto(kind, 0.0_dp, s, ops)
    call graph(trim(tag), s, ops, rcl, levels, cg)
    call periodic_graph_shells(kind, 4, int((2.0_dp * rcl)**2 + 1.0e-9_dp), g, coord, bv)
    call check(cg%nslot == size(bv, 2), trim(tag)//': number of slots as asr_lattices')
    if (cg%nslot == size(bv, 2)) call check(all(cg%bvec == 0.5_dp * real(bv, dp)), &
                                            trim(tag)//': bond vectors as asr_lattices, in order')
    call check(abs(cg%radius - 2.0_dp * levels * maxval(norm2(cg%bvec, 1))) < 1.0e-14_dp, trim(tag)//': radius 2 levels bmax')
    call lattice_cluster(kind, cg%radius, vec)
    call check(size(vec, 2) == cg%n, trim(tag)//': the sites of lattice_cluster within the radius')
    if (size(vec, 2) == cg%n) call check(all(vec == cg%disp), trim(tag)//': same displacements in the same order, bit for bit')
    call check(all(cg%g%nbr /= 0), trim(tag)//': Bravais lattice: no entry 0')
    call periodic_graph(kind, 2, g, coord)                           ! only for the order of the 48 rotations
    call cubic_site_perms(coord, 4, p48, m48)
    call graph_site_ops(cg, s, ops, perm, rot, idx)
    call check(size(perm, 2) == 48, trim(tag)//': 48 site operations')
    if (size(perm, 2) == 48) call check(all(nint(rot) == m48) .and. maxval(abs(rot - real(m48, dp))) == 0.0_dp, &
                                        trim(tag)//': rotations in the order of cubic_site_perms')
    call slot_checks(s, cg, trim(tag))
    call size_checks(cg, trim(tag))
    call op_checks(s, cg, ops, trim(tag))
  end subroutine cubic_case

  subroutine b2_case(rcl, levels)
    real(dp), intent(in) :: rcl
    integer, intent(in) :: levels
    type(structure_t) :: s, sb
    type(symop_t), allocatable :: ops(:), opsb(:)
    type(crystal_graph_t) :: cg, cgb
    integer, allocatable :: perm(:, :), idx(:), permb(:, :), idxb(:)
    real(dp), allocatable :: rot(:, :, :), rotb(:, :, :)
    integer :: b0, i
    call proto('bcc', 0.0_dp, sb, opsb)
    call graph('bcc', sb, opsb, rcl, levels, cgb)
    call graph_site_ops(cgb, sb, opsb, permb, rotb, idxb)
    do b0 = 1, 2
      call proto('b2', 0.0_dp, s, ops)
      call graph('b2', s, ops, rcl, levels, cg, start=b0)
      call check(cg%n == cgb%n .and. all(cg%disp == cgb%disp) .and. all(cg%g%nbr == cgb%g%nbr), &
                 'b2: the bcc graph, site by site')
      call check(all([(cg%basis(i) == 1 + modulo(nint(2.0_dp * cg%disp(1, i)) + b0 - 1, 2), i=1, cg%n)]), &
                 'b2: basis sites alternate, start on b0')
      call check(cg%norbit == 2, 'b2: two symmetry-distinct sites')
      call graph_site_ops(cg, s, ops, perm, rot, idx)
      call check(size(perm, 2) == 48, 'b2: 48 site operations')
      if (size(perm, 2) == 48) call check(all(perm == permb), 'b2: site permutations of bcc')
      call slot_checks(s, cg, 'b2')
      call size_checks(cg, 'b2')
      call op_checks(s, cg, ops, 'b2')
    end do
  end subroutine b2_case

  subroutine hcp_case(ca, nshell, nsecond)
    real(dp), intent(in) :: ca
    integer, intent(in) :: nshell(2), nsecond
    type(structure_t) :: s
    type(symop_t), allocatable :: ops(:)
    type(crystal_graph_t) :: cg
    integer, allocatable :: perm(:, :), idx(:)
    real(dp), allocatable :: rot(:, :, :)
    integer :: b, k, i, n1, n2, nmono
    real(dp) :: d, d1
    character(len=24) :: tag
    write (tag, '(a,f5.3)') 'hcp c/a=', ca
    call proto('hcp', ca, s, ops)
    call graph(trim(tag), s, ops, 1.05_dp, 2, cg)
    call check(cg%norbit == 1, trim(tag)//': one symmetry-distinct site')
    do b = 1, 2
      d1 = huge(1.0_dp)
      do k = 1, cg%nslot
        if (cg%btgt(k, b) > 0) d1 = min(d1, norm2(cg%bvec(:, k)))
      end do
      n1 = 0; n2 = 0
      do k = 1, cg%nslot
        if (cg%btgt(k, b) == 0) cycle
        d = norm2(cg%bvec(:, k))
        if (abs(d - d1) < 1.0e-9_dp) then; n1 = n1 + 1; else; n2 = n2 + 1; end if
      end do
      call check(n1 == nshell(b) .and. n2 == nsecond, trim(tag)//': nearest shells of each site')
      if (n1 /= nshell(b) .or. n2 /= nsecond) print '(2x,a,i0,a,2i4)', 'site ', b, ' shells ', n1, n2
    end do
    call check(cg%nslot == 18, trim(tag)//': union of the bond vectors of A and B (6 in plane, 6 + 6 out of plane)')
    call check(all(count(cg%btgt > 0, 1) == 12), trim(tag)//': 12 of the 18 slots filled on each site')
    call graph_site_ops(cg, s, ops, perm, rot, idx)
    nmono = count([(is_monomial(rot(:, :, i)), i=1, size(rot, 3))])
    call check(size(perm, 2) == 12 .and. nmono == 4, trim(tag)//': site group D3h (12), 4 monomial')
    call check(all([(is_monomial(rot(:, :, i)), i=1, nmono)]), trim(tag)//': monomial operations first')
    call slot_checks(s, cg, trim(tag))
    call size_checks(cg, trim(tag))
    call op_checks(s, cg, ops, trim(tag))
  end subroutine hcp_case

  subroutine nias_case()
    type(structure_t) :: s
    type(symop_t), allocatable :: ops(:)
    type(crystal_graph_t) :: cg
    integer, allocatable :: perm(:, :), idx(:)
    real(dp), allocatable :: rot(:, :, :)
    integer :: b, k, nte, nmn, nother
    real(dp) :: dmin
    do b = 1, 4
      call proto('nias', 1.617_dp, s, ops)
      call graph('nias', s, ops, 0.9_dp, 2, cg, start=b)
      call check(cg%norbit == 2, 'nias: two symmetry-distinct sites')
      dmin = minval([(norm2(cg%bvec(:, k)), k=1, cg%nslot)], mask=cg%btgt(:, b) > 0)
      nte = 0; nmn = 0; nother = 0
      do k = 1, cg%nslot
        if (cg%btgt(k, b) == 0) cycle
        if (abs(norm2(cg%bvec(:, k)) - dmin) < 1.0e-9_dp .and. s%label(cg%btgt(k, b)) /= s%label(b)) then
          nte = nte + 1
        else if (s%label(cg%btgt(k, b)) == s%label(b)) then
          nmn = nmn + 1
        else
          nother = nother + 1
        end if
      end do
      if (s%label(b) == 1) then
        call check(nte == 6 .and. nmn == 2 .and. nother == 0, 'nias: cation shells 6 anions, then 2 cations')
      else
        call check(nte == 6 .and. nmn == 0 .and. nother == 0, 'nias: anion shell 6 cations')
      end if
      call graph_site_ops(cg, s, ops, perm, rot, idx)
      call check(size(perm, 2) == 12, 'nias: site group of order 12 (D3d cation, D3h anion)')
      call slot_checks(s, cg, 'nias')
      call size_checks(cg, 'nias')
      call op_checks(s, cg, ops, 'nias')
    end do
  end subroutine nias_case

  !> rmin larger than 2 levels bmax sets the radius; levels 0 gives the start site alone.
  subroutine rmin_case()
    type(structure_t) :: s
    type(symop_t), allocatable :: ops(:)
    type(crystal_graph_t) :: cg
    real(dp), allocatable :: vec(:, :)
    call proto('fcc', 0.0_dp, s, ops)
    call graph('fcc rmin', s, ops, 0.75_dp, 1, cg, rmin=3.1_dp)
    call lattice_cluster('fcc', 3.1_dp, vec)
    call check(abs(cg%radius - 3.1_dp) < 1.0e-14_dp .and. cg%n == size(vec, 2), 'rmin sets the cluster radius')
    call graph('fcc levels 0', s, ops, 0.75_dp, 0, cg)
    call check(cg%n == 1 .and. all(cg%g%nbr == -1), 'levels 0: the start site alone, every neighbour outside')
  end subroutine rmin_case

  !> Slot consistency: exact labels, Cartesian vectors, entries 0 and -1, symmetry of the graph.
  subroutine slot_checks(s, cg, tag)
    type(structure_t), intent(in) :: s
    type(crystal_graph_t), intent(in) :: cg
    character(len=*), intent(in) :: tag
    integer :: i, j, k, k2, b, t
    logical :: lab_ok, cart_ok, zero_ok, sym_ok
    real(dp) :: v(3)
    lab_ok = .true.; cart_ok = .true.; zero_ok = .true.; sym_ok = .true.
    do k = 1, cg%nslot
      do b = 1, cg%nb
        if (cg%btgt(k, b) == 0) cycle
        v = matmul(s%avec, real(cg%bdn(:, k, b), dp) + s%frac(:, cg%btgt(k, b)) - s%frac(:, b))
        if (maxval(abs(v - cg%bvec(:, k))) > 1.0e-12_dp) cart_ok = .false.
      end do
    end do
    do i = 1, cg%n
      b = cg%basis(i)
      if (graph_site(cg, cg%cell(:, i), b) /= i) lab_ok = .false.
      v = matmul(s%avec, real(cg%cell(:, i), dp) + s%frac(:, b) - s%frac(:, cg%b0))
      if (maxval(abs(v - cg%disp(:, i))) > 1.0e-12_dp) cart_ok = .false.
      do k = 1, cg%nslot
        j = cg%g%nbr(cg%g%nbr_ptr(i) + k - 1)
        if ((j == 0) .neqv. (cg%btgt(k, b) == 0)) zero_ok = .false.
        if (j == 0) cycle
        t = graph_site(cg, cg%cell(:, i) + cg%bdn(:, k, b), cg%btgt(k, b))
        if (j == -1) then
          if (t /= 0 .or. norm2(cg%disp(:, i) + cg%bvec(:, k)) <= cg%radius + rtol) zero_ok = .false.
          cycle
        end if
        if (j /= t .or. cg%basis(j) /= cg%btgt(k, b)) lab_ok = .false.
        if (maxval(abs(cg%disp(:, j) - cg%disp(:, i) - cg%bvec(:, k))) > 1.0e-12_dp) cart_ok = .false.
        k2 = find_slot(cg, -cg%bvec(:, k))
        if (k2 == 0) then; sym_ok = .false.; cycle; end if
        if (cg%g%nbr(cg%g%nbr_ptr(j) + k2 - 1) /= i) sym_ok = .false.
      end do
    end do
    call check(lab_ok, tag//': neighbour = site + slot offset, exact labels')
    call check(cart_ok, tag//': neighbour - site = bvec(k), displacements exact')
    call check(zero_ok, tag//': entries 0 where the basis site has no such bond, -1 exactly beyond the radius')
    call check(sym_ok, tag//': the graph is symmetric (slot -v leads back)')
  end subroutine slot_checks

  !> Every site within 2 levels - 1 bonds of site 1 has all its neighbours in the cluster.
  subroutine size_checks(cg, tag)
    type(crystal_graph_t), intent(in) :: cg
    character(len=*), intent(in) :: tag
    integer, allocatable :: dist(:)
    integer :: i
    logical :: ok
    call graph_distances(cg%g, 1, dist)
    ok = .true.
    do i = 1, cg%n
      if (dist(i) < 0 .or. dist(i) > 2 * cg%levels - 1) cycle
      if (any(cg%g%nbr(cg%g%nbr_ptr(i):cg%g%nbr_ptr(i + 1) - 1) < 0)) ok = .false.
    end do
    call check(ok, tag//': sites within 2 levels - 1 bonds have all their neighbours')
    call check(any(cg%g%nbr == -1), tag//': the cluster is finite (entries -1 at the edge)')
  end subroutine size_checks

  !> Site operations: neighbours onto neighbours with the rotated bond vector (entries 0 and -1 too), displacements rotated.
  subroutine op_checks(s, cg, ops, tag)
    type(structure_t), intent(in) :: s
    type(crystal_graph_t), intent(in) :: cg
    type(symop_t), intent(in) :: ops(:)
    character(len=*), intent(in) :: tag
    integer, allocatable :: perm(:, :), idx(:), kmap(:)
    real(dp), allocatable :: rot(:, :, :)
    integer :: g, i, k, j, j2
    logical :: nb_ok, disp_ok, map_ok
    call graph_site_ops(cg, s, ops, perm, rot, idx)
    nb_ok = .true.; disp_ok = .true.; map_ok = .true.
    allocate (kmap(cg%nslot))
    do g = 1, size(perm, 2)
      do k = 1, cg%nslot
        kmap(k) = find_slot(cg, matmul(rot(:, :, g), cg%bvec(:, k)))
      end do
      if (any(kmap == 0)) then; map_ok = .false.; cycle; end if
      do i = 1, cg%n
        if (maxval(abs(cg%disp(:, perm(i, g)) - matmul(rot(:, :, g), cg%disp(:, i)))) > 1.0e-12_dp) disp_ok = .false.
        do k = 1, cg%nslot
          j = cg%g%nbr(cg%g%nbr_ptr(i) + k - 1)
          j2 = cg%g%nbr(cg%g%nbr_ptr(perm(i, g)) + kmap(k) - 1)
          if (j <= 0) then
            if (j2 /= j) nb_ok = .false.
          else
            if (j2 /= perm(j, g)) nb_ok = .false.
          end if
        end do
      end do
    end do
    call check(map_ok, tag//': every site operation maps the slot set onto itself')
    call check(nb_ok, tag//': site operations map neighbours onto neighbours with the rotated bond vector')
    call check(disp_ok, tag//': site operations rotate the displacements')
    call check(all(perm(1, :) == 1), tag//': site operations fix site 1')
  end subroutine op_checks

  integer function find_slot(cg, v)
    type(crystal_graph_t), intent(in) :: cg
    real(dp), intent(in) :: v(3)
    integer :: k
    find_slot = 0
    do k = 1, cg%nslot
      if (maxval(abs(cg%bvec(:, k) - v)) < 1.0e-9_dp) then; find_slot = k; return; end if
    end do
  end function find_slot
end program test_crystal_graph
