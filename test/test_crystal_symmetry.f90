!> Symmetry of crystal structures found from the structure alone (asr_crystal, asr_prototypes, asr_crystal_symmetry).
!> Gate G1 of notes/DESIGN_structures.md for every prototype: space-group (point part) and site-group orders as in the International
!> Tables; the integer action on sites reproduces the Cartesian one; orbital matrices are orthogonal, also for threefold axes;
!> orbital groups, representatives and the monomial reduction group.
program test_crystal_symmetry
  use asr_kinds, only: dp
  use asr_crystal, only: structure_t, to_cart
  use asr_prototypes, only: prototype_structure
  use asr_crystal_symmetry, only: symop_t, space_group, site_ops, is_monomial, orbital_matrix, orbital_groups, reduction_ops
  use asr_testing
  implicit none
  integer, parameter :: cubic_grp(9) = [1, 2, 2, 2, 3, 3, 4, 3, 4]    ! s | p y z x | xy yz 3z2 xz x2-y2
  integer, parameter :: hex_grp(9) = [1, 2, 3, 2, 4, 5, 6, 5, 4]      ! s | p_xy, p_z | d_xy,x2-y2 | d_yz,xz | d_3z2
  print '(a)', 'test_crystal_symmetry'
  ! name, c/a, space-group size, per-site site-group order, monomial site order, number of orbital groups, reduction order
  call case_('sc', 0.0_dp, 48, [48], [48], 4, 16, cubic_grp, [1, 3, 5, 7])
  call case_('bcc', 0.0_dp, 48, [48], [48], 4, 16, cubic_grp, [1, 3, 5, 7])
  call case_('fcc', 0.0_dp, 48, [48], [48], 4, 16, cubic_grp, [1, 3, 5, 7])
  call case_('b2', 0.0_dp, 48, [48, 48], [48, 48], 4, 16, cubic_grp, [1, 3, 5, 7])
  call case_('hex', 1.6_dp, 24, [24], [8], 6, 8, hex_grp, [0])
  call case_('hcp', 0.0_dp, 24, [12, 12], [4, 4], 6, 4, hex_grp, [0])
  call case_('nias', 1.617_dp, 24, [12, 12, 12, 12], [4, 4, 4, 4], 6, 4, hex_grp, [0])
  call report()
contains

  subroutine case_(name, ca, nsg, nsite, nmono, ngrp_ref, nred_ref, grp_ref, rep_ref)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: ca
    integer, intent(in) :: nsg, nsite(:), nmono(:), ngrp_ref, nred_ref, grp_ref(9), rep_ref(:)
    type(structure_t) :: s
    type(symop_t), allocatable :: ops(:)
    character(len=:), allocatable :: msg
    integer, allocatable :: idx(:), red(:)
    real(dp), allocatable :: d(:, :, :)
    logical, allocatable :: mono(:)
    integer :: b, k, grp(9), ngrp, rep(9), weight(9)
    real(dp) :: err_act, err_orth, r(3)
    if (ca > 0.0_dp) then
      call prototype_structure(name, s, msg, c_over_a=ca)
    else
      call prototype_structure(name, s, msg)
    end if
    call check(len(msg) == 0, name//': prototype builds')
    call space_group(s, ops)
    call check(size(ops) == nsg, name//': space-group operations')
    if (size(ops) /= nsg) print '(2x,a,i0,a,i0)', 'found ', size(ops), ' expected ', nsg
    ! exact integer action on sites against the Cartesian action
    err_act = 0.0_dp
    do k = 1, size(ops)
      do b = 1, s%nb
        r = matmul(ops(k)%rot, s%tau(:, b)) + to_cart(s, ops(k)%t) - s%tau(:, ops(k)%perm(b)) - to_cart(s, real(ops(k)%shift(:, b), dp))
        err_act = max(err_act, maxval(abs(r)))
      end do
    end do
    call check(err_act < 1.0e-10_dp, name//': integer site action equals the Cartesian one')
    do b = 1, s%nb
      call site_ops(ops, b, idx)
      call check(size(idx) == nsite(b), name//': site-group order')
      if (size(idx) /= nsite(b)) print '(2x,a,i0,a,i0,a,i0)', 'site ', b, ' found ', size(idx), ' expected ', nsite(b)
      allocate (d(9, 9, size(idx)), mono(size(idx)))
      err_orth = 0.0_dp
      do k = 1, size(idx)
        call orbital_matrix(ops(idx(k))%rot, d(:, :, k))
        mono(k) = is_monomial(ops(idx(k))%rot)
        err_orth = max(err_orth, maxval(abs(matmul(transpose(d(:, :, k)), d(:, :, k)) - eye9())))
      end do
      call check(err_orth < 1.0e-10_dp, name//': orbital matrices are orthogonal')
      call check(count(mono) == nmono(b), name//': monomial site operations')
      if (count(mono) /= nmono(b)) print '(2x,a,i0,a,i0)', 'found ', count(mono), ' expected ', nmono(b)
      call orbital_groups(d, mono, grp, ngrp, rep, weight)
      call check(ngrp == ngrp_ref, name//': number of orbital groups')
      call check(all(grp == grp_ref), name//': orbital groups')
      if (any(grp /= grp_ref)) print '(2x,a,9i2)', 'groups ', grp
      call check(sum(weight(1:ngrp)) == 9, name//': weights sum to nine orbitals')
      if (rep_ref(1) > 0) call check(all(rep(1:ngrp) == rep_ref), name//': representatives as today (s, p_z, d_xy, d_3z2)')
      call reduction_ops(d, mono, rep(1:ngrp), red)
      call check(size(red) == nred_ref, name//': reduction group order')
      if (size(red) /= nred_ref) print '(2x,a,i0,a,i0)', 'found ', size(red), ' expected ', nred_ref
      deallocate (d, mono)
    end do
  end subroutine case_

  pure function eye9() result(e)
    real(dp) :: e(9, 9)
    integer :: i
    e = 0.0_dp
    do i = 1, 9; e(i, i) = 1.0_dp; end do
  end function eye9
end program test_crystal_symmetry
