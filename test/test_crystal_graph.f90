!> Periodic graph of a crystal structure (asr_crystal_graph), step 2 of notes/DESIGN_structures.md.
!>  (C) sc, bcc, fcc on L conventional cubes: the graph, the bond vectors, the minimum-image displacements and the 48 site
!>      operations are those of asr_lattices::periodic_graph_shells and asr_symmetry::cubic_site_perms, site by site.  B2 is
!>      the bcc graph with the two basis sites alternating.
!>  (H) hcp: 12 nearest neighbours per site at the ideal c/a, 6 + 6 otherwise; NiAs: the cation has 6 anions and then 2
!>      cations as its nearest shells, the anion 6 cations.
!>  (S) every slot of every site: neighbour - site = bvec(k) as a vector of the supercell lattice (exact labels and Cartesian);
!>      the graph is symmetric; a supercell that is too small is flagged.
!>  (G) every site operation maps neighbours onto neighbours with the rotated bond vector and keeps the displacement lengths;
!>      the displacements are the shortest images also on skewed supercells.
program test_crystal_graph
  use asr_kinds, only: dp
  use asr_crystal, only: structure_t
  use asr_prototypes, only: prototype_structure
  use asr_crystal_symmetry, only: symop_t, space_group, is_monomial
  use asr_crystal_graph
  use asr_augspace, only: graph_t
  use asr_lattices, only: periodic_graph_shells
  use asr_symmetry, only: cubic_site_perms
  use asr_testing
  implicit none
  real(dp), parameter :: ca_ideal = 1.632993161855452_dp
  print '(a)', 'test_crystal_graph'
  call cubic_case('sc', 3, 1.2_dp)
  call cubic_case('bcc', 4, 1.2_dp)
  call cubic_case('bcc', 3, 1.45_dp)
  call cubic_case('fcc', 3, 1.2_dp)
  call cubic_case('fcc', 4, 0.75_dp)
  call b2_case(4, 1.2_dp)
  call alias_case()
  call hcp_case(ca_ideal, [12, 12], 0)
  call hcp_case(1.7_dp, [6, 6], 6)
  call nias_case()
  call skew_case()
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

  subroutine graph(tag, s, ops, msc, rcl, cg, start)
    character(len=*), intent(in) :: tag
    type(structure_t), intent(in) :: s
    type(symop_t), intent(in) :: ops(:)
    integer, intent(in) :: msc(3, 3)
    real(dp), intent(in) :: rcl
    type(crystal_graph_t), intent(out) :: cg
    integer, intent(in), optional :: start
    character(len=:), allocatable :: msg
    call crystal_graph(s, ops, msc, rcl, cg, msg, start)
    call check(len(msg) == 0, tag//': graph builds')
    if (len(msg) > 0) print '(2x,a)', msg
  end subroutine graph

  function cubes(s, l) result(msc)
    type(structure_t), intent(in) :: s
    integer, intent(in) :: l
    integer :: msc(3, 3), i
    real(dp) :: c(3, 3)
    character(len=:), allocatable :: msg
    c = 0.0_dp; do i = 1, 3; c(i, i) = real(l, dp); end do
    call supercell_from_cart(s, c, msc, msg)
    if (len(msg) > 0) error stop msg
  end function cubes

  subroutine cubic_case(kind, l, rcl)
    character(len=*), intent(in) :: kind
    integer, intent(in) :: l
    real(dp), intent(in) :: rcl
    type(structure_t) :: s
    type(symop_t), allocatable :: ops(:)
    type(crystal_graph_t) :: cg
    type(graph_t) :: g
    integer, allocatable :: coord(:, :), bv(:, :), p48(:, :), m48(:, :, :), perm(:, :), idx(:)
    real(dp), allocatable :: rot(:, :, :)
    character(len=40) :: tag
    integer :: mg, i
    real(dp) :: err
    write (tag, '(a,a,i0,a,f4.2)') kind, ' L=', l, ' rcl=', rcl
    call proto(kind, 0.0_dp, s, ops)
    call graph(trim(tag), s, ops, cubes(s, l), rcl, cg)
    call periodic_graph_shells(kind, l, int((2.0_dp * rcl)**2 + 1.0e-9_dp), g, coord, bv)
    mg = 2 * l
    call check(cg%n == g%n, trim(tag)//': number of sites as asr_lattices')
    call check(cg%nslot == size(bv, 2), trim(tag)//': number of slots as asr_lattices')
    if (cg%n /= g%n .or. cg%nslot /= size(bv, 2)) return
    call check(all(cg%bvec == 0.5_dp * real(bv, dp)), trim(tag)//': bond vectors as asr_lattices, in order')
    call check(all(cg%g%nbr_ptr == g%nbr_ptr) .and. all(cg%g%nbr == g%nbr), trim(tag)//': neighbour lists as asr_lattices, site by site')
    err = 0.0_dp
    do i = 1, cg%n
      err = max(err, maxval(abs(cg%disp(:, i) - 0.5_dp * real(modulo(coord(:, i) + mg / 2, mg) - mg / 2, dp))))
    end do
    call check(err == 0.0_dp, trim(tag)//': minimum-image displacements as the half-lattice box')
    call cubic_site_perms(coord, mg, p48, m48)
    call graph_site_ops(cg, s, ops, perm, rot, idx)
    call check(size(perm, 2) == 48, trim(tag)//': 48 site operations')
    if (size(perm, 2) /= 48) return
    call check(all(perm == p48), trim(tag)//': site permutations as cubic_site_perms, in order')
    call check(all(nint(rot) == m48) .and. maxval(abs(rot - real(m48, dp))) == 0.0_dp, trim(tag)//': rotations as cubic_site_perms')
    call check(.not. cg%aliased, trim(tag)//': not aliased')
    call slot_checks(s, cg, trim(tag))
    call op_checks(s, cg, ops, trim(tag))
  end subroutine cubic_case

  subroutine b2_case(l, rcl)
    integer, intent(in) :: l
    real(dp), intent(in) :: rcl
    type(structure_t) :: s
    type(symop_t), allocatable :: ops(:)
    type(crystal_graph_t) :: cg
    type(graph_t) :: g
    integer, allocatable :: coord(:, :), bv(:, :), perm(:, :), idx(:), p48(:, :), m48(:, :, :)
    real(dp), allocatable :: rot(:, :, :)
    integer :: b0, i
    do b0 = 1, 2
      call proto('b2', 0.0_dp, s, ops)
      call graph('b2', s, ops, cubes(s, l), rcl, cg, start=b0)
      call periodic_graph_shells('bcc', l, int((2.0_dp * rcl)**2 + 1.0e-9_dp), g, coord, bv)
      call check(cg%n == g%n .and. all(cg%g%nbr == g%nbr), 'b2: the bcc graph, site by site')
      call check(all([(cg%basis(i) == 1 + modulo(coord(1, i) + b0 - 1, 2), i=1, cg%n)]), 'b2: basis sites alternate, start on b0')
      call check(cg%norbit == 2, 'b2: two symmetry-distinct sites')
      call cubic_site_perms(coord, 2 * l, p48, m48)
      call graph_site_ops(cg, s, ops, perm, rot, idx)
      call check(size(perm, 2) == 48, 'b2: 48 site operations')
      if (size(perm, 2) == 48) call check(all(perm == p48), 'b2: site permutations as cubic_site_perms')
      call slot_checks(s, cg, 'b2')
      call op_checks(s, cg, ops, 'b2')
    end do
  end subroutine b2_case

  !> bcc on 2 cubes with the second shell: (1, 0, 0) and (-1, 0, 0) are one site.
  subroutine alias_case()
    type(structure_t) :: s
    type(symop_t), allocatable :: ops(:)
    type(crystal_graph_t) :: cg
    call proto('bcc', 0.0_dp, s, ops)
    call graph('bcc L=2', s, ops, cubes(s, 2), 1.2_dp, cg)
    call check(cg%aliased, 'bcc L=2 with two shells is flagged as aliased')
    call graph('bcc L=2', s, ops, cubes(s, 2), 0.9_dp, cg)
    call check(.not. cg%aliased, 'bcc L=2 with one shell is not aliased')
  end subroutine alias_case

  subroutine hcp_case(ca, nshell, nsecond)
    real(dp), intent(in) :: ca
    integer, intent(in) :: nshell(2), nsecond
    type(structure_t) :: s
    type(symop_t), allocatable :: ops(:)
    type(crystal_graph_t) :: cg
    integer, allocatable :: perm(:, :), idx(:)
    real(dp), allocatable :: rot(:, :, :)
    integer :: msc(3, 3), b, k, i, n1, n2, nmono
    real(dp) :: d, d1
    character(len=24) :: tag
    write (tag, '(a,f5.3)') 'hcp c/a=', ca
    msc = 0; msc(1, 1) = 6; msc(2, 2) = 6; msc(3, 3) = 4
    call proto('hcp', ca, s, ops)
    call graph(trim(tag), s, ops, msc, 1.05_dp, cg)
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
    call check(size(perm, 2) == 12 .and. nmono == 4, trim(tag)//': site group D3h (12) on the graph, 4 monomial')
    call check(all([(is_monomial(rot(:, :, i)), i=1, nmono)]), trim(tag)//': monomial operations first')
    call slot_checks(s, cg, trim(tag))
    call op_checks(s, cg, ops, trim(tag))
  end subroutine hcp_case

  subroutine nias_case()
    type(structure_t) :: s
    type(symop_t), allocatable :: ops(:)
    type(crystal_graph_t) :: cg
    integer, allocatable :: perm(:, :), idx(:)
    real(dp), allocatable :: rot(:, :, :)
    integer :: msc(3, 3), b, k, nte, nmn, nother
    real(dp) :: dmin
    msc = 0; msc(1, 1) = 5; msc(2, 2) = 5; msc(3, 3) = 4
    do b = 1, 4
      call proto('nias', 1.617_dp, s, ops)
      call graph('nias', s, ops, msc, 0.9_dp, cg, start=b)
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
      call check(size(perm, 2) == 12, 'nias: site group of order 12 on the graph (D3d cation, D3h anion)')
      call slot_checks(s, cg, 'nias')
      call op_checks(s, cg, ops, 'nias')
    end do
  end subroutine nias_case

  !> A skewed hexagonal supercell (non-diagonal msc): site count, slot consistency, shortest images by brute force.
  subroutine skew_case()
    type(structure_t) :: s
    type(symop_t), allocatable :: ops(:)
    type(crystal_graph_t) :: cg
    integer :: msc(3, 3), i, e1, e2, e3
    real(dp) :: p(3), dmin, v(3)
    logical :: ok
    msc = reshape([4, 2, 0, -2, 3, 0, 0, 0, 3], [3, 3])
    call proto('nias', 1.617_dp, s, ops)
    call graph('skewed nias', s, ops, msc, 0.9_dp, cg)
    call check(cg%n == 4 * 16 * 3, 'skewed supercell: |det msc| cells of 4 sites')
    ok = .true.
    do i = 1, cg%n
      p = real(cg%cell(:, i), dp) + s%frac(:, cg%basis(i)) - s%frac(:, cg%b0)
      dmin = huge(1.0_dp)
      do e3 = -3, 3; do e2 = -3, 3; do e1 = -3, 3
        v = matmul(s%avec, p + real(matmul(msc, [e1, e2, e3]), dp))
        dmin = min(dmin, norm2(v))
      end do; end do; end do
      if (norm2(cg%disp(:, i)) > dmin + 1.0e-12_dp) ok = .false.
      v = matmul(s%avec, p) - cg%disp(:, i)                          ! disp is an image of the site
      if (.not. lattice_vector(s, cg, v)) ok = .false.
    end do
    call check(ok, 'skewed supercell: displacements are the shortest images')
    call check(all(cg%disp(:, 1) == 0.0_dp), 'skewed supercell: site 1 is the origin')
    call slot_checks(s, cg, 'skewed nias')
    call op_checks(s, cg, ops, 'skewed nias')
  end subroutine skew_case

  !> Is the Cartesian vector v a vector of the supercell lattice?
  logical function lattice_vector(s, cg, v)
    type(structure_t), intent(in) :: s
    type(crystal_graph_t), intent(in) :: cg
    real(dp), intent(in) :: v(3)
    real(dp) :: f(3), fs(3)
    f = matmul(transpose(s%bvec), v)                                 ! primitive coordinates
    lattice_vector = all(abs(f - anint(f)) < 1.0e-9_dp)
    if (.not. lattice_vector) return
    fs = matmul(inv_int(cg%msc), f)
    lattice_vector = all(abs(fs - anint(fs)) < 1.0e-9_dp)
  end function lattice_vector

  function inv_int(m) result(a)
    integer, intent(in) :: m(3, 3)
    real(dp) :: a(3, 3), r(3, 3), d
    r = real(m, dp)
    a(1, 1) = r(2, 2) * r(3, 3) - r(2, 3) * r(3, 2); a(1, 2) = r(1, 3) * r(3, 2) - r(1, 2) * r(3, 3)
    a(1, 3) = r(1, 2) * r(2, 3) - r(1, 3) * r(2, 2); a(2, 1) = r(2, 3) * r(3, 1) - r(2, 1) * r(3, 3)
    a(2, 2) = r(1, 1) * r(3, 3) - r(1, 3) * r(3, 1); a(2, 3) = r(1, 3) * r(2, 1) - r(1, 1) * r(2, 3)
    a(3, 1) = r(2, 1) * r(3, 2) - r(2, 2) * r(3, 1); a(3, 2) = r(1, 2) * r(3, 1) - r(1, 1) * r(3, 2)
    a(3, 3) = r(1, 1) * r(2, 2) - r(1, 2) * r(2, 1)
    d = r(1, 1) * a(1, 1) + r(1, 2) * a(2, 1) + r(1, 3) * a(3, 1)
    a = a / d
  end function inv_int

  !> Slot consistency: exact labels, Cartesian vectors, zero slots, symmetry of the graph.
  subroutine slot_checks(s, cg, tag)
    type(structure_t), intent(in) :: s
    type(crystal_graph_t), intent(in) :: cg
    character(len=*), intent(in) :: tag
    integer :: i, j, k, k2, b
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
      do k = 1, cg%nslot
        j = cg%g%nbr(cg%g%nbr_ptr(i) + k - 1)
        if ((j == 0) .neqv. (cg%btgt(k, b) == 0)) zero_ok = .false.
        if (j == 0) cycle
        if (cg%basis(j) /= cg%btgt(k, b)) lab_ok = .false.
        if (any(graph_reduce_cell(cg, cg%cell(:, j) - cg%cell(:, i) - cg%bdn(:, k, b)) /= 0)) lab_ok = .false.
        if (.not. lattice_vector(s, cg, cg%disp(:, j) - cg%disp(:, i) - cg%bvec(:, k))) cart_ok = .false.
        k2 = find_slot(cg, -cg%bvec(:, k))
        if (k2 == 0) then; sym_ok = .false.; cycle; end if
        if (cg%g%nbr(cg%g%nbr_ptr(j) + k2 - 1) /= i) sym_ok = .false.
      end do
    end do
    call check(lab_ok, tag//': neighbour - site = slot offset, exact labels')
    call check(cart_ok, tag//': neighbour - site = bvec(k) modulo the supercell, Cartesian')
    call check(zero_ok, tag//': zero slots exactly where the basis site has no such bond')
    call check(sym_ok, tag//': the graph is symmetric (slot -v leads back)')
  end subroutine slot_checks

  !> Site operations: neighbours onto neighbours with the rotated bond vector, displacement lengths kept.
  subroutine op_checks(s, cg, ops, tag)
    type(structure_t), intent(in) :: s
    type(crystal_graph_t), intent(in) :: cg
    type(symop_t), intent(in) :: ops(:)
    character(len=*), intent(in) :: tag
    integer, allocatable :: perm(:, :), idx(:), kmap(:)
    real(dp), allocatable :: rot(:, :, :)
    integer :: g, i, k, j, j2
    logical :: nb_ok, len_ok, map_ok
    call graph_site_ops(cg, s, ops, perm, rot, idx)
    nb_ok = .true.; len_ok = .true.; map_ok = .true.
    allocate (kmap(cg%nslot))
    do g = 1, size(perm, 2)
      do k = 1, cg%nslot
        kmap(k) = find_slot(cg, matmul(rot(:, :, g), cg%bvec(:, k)))
      end do
      if (any(kmap == 0)) then; map_ok = .false.; cycle; end if
      do i = 1, cg%n
        if (abs(norm2(cg%disp(:, perm(i, g))) - norm2(cg%disp(:, i))) > 1.0e-9_dp) len_ok = .false.
        do k = 1, cg%nslot
          j = cg%g%nbr(cg%g%nbr_ptr(i) + k - 1)
          j2 = cg%g%nbr(cg%g%nbr_ptr(perm(i, g)) + kmap(k) - 1)
          if (j == 0) then
            if (j2 /= 0) nb_ok = .false.
          else
            if (j2 /= perm(j, g)) nb_ok = .false.
          end if
        end do
      end do
    end do
    call check(map_ok, tag//': every site operation maps the slot set onto itself')
    call check(nb_ok, tag//': site operations map neighbours onto neighbours with the rotated bond vector')
    call check(len_ok, tag//': site operations keep the displacement lengths')
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
