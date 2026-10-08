!> Collinear magnetism generated from scratch (no tabulated input): bcc Fe, fcc Ni, fcc Cu.
!> Gates are literature LDA/ASA values: Fe ~2.2 muB, Ni ~0.6 muB, Cu non-magnetic (a seeded
!> polarisation must relax to zero).  Moments depend on k-mesh/smearing at the 2-3% level.
program test_asa_mag
  use asr_kinds, only: dp
  use asr_asa
  use asr_atom, only: shell_t
  use asr_lmto_cell
  use asr_testing
  implicit none
  print '(a)', 'test_asa_mag'
  call run('Fe', 26.0_dp, 8, 2.8665_dp, [-0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp], 2, 6, 2.15_dp, 2.30_dp)
  call run('Ni', 28.0_dp, 10, 3.524_dp, [0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp], 2, 8, 0.58_dp, 0.68_dp)
  call run('Cu', 29.0_dp, 11, 3.61_dp, [0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp], 1, 10, -0.002_dp, 0.002_dp)
  call report()
contains
  subroutine run(name, z, nval, a_ang, avec, n4s, n3d, mlo, mhi)
    character(*), intent(in) :: name
    real(dp), intent(in) :: z, a_ang, avec(9), mlo, mhi
    integer, intent(in) :: nval, n4s, n3d
    type(asa_input_t) :: inp
    type(asa_result_t) :: res
    type(cell_t) :: cell
    allocate (inp%cls(1))
    inp%cls(1)%z = z; inp%cls(1)%nval = nval
    inp%cls(1)%core = core_shells_for([4, 4, 3, 4]); inp%cls(1)%nvl = [4, 4, 3]
    inp%cls(1)%atom_shells = [inp%cls(1)%core, shell_t(4, 0, real(n4s, dp)), shell_t(3, 2, real(n3d, dp))]
    inp%alat = a_ang / 0.529177_dp
    inp%nspin = 2; inp%zeta0 = 0.4_dp; inp%nk = 12; inp%smear = 0.01_dp; inp%beta = 0.25_dp; inp%maxit = 200; inp%tol = 1.0e-5_dp
    cell%nat = 1; cell%avec = reshape(avec, [3, 3])
    allocate (cell%tau(3, 1), cell%isp(1)); cell%tau = 0.0_dp; cell%isp = 1
    call setup_cell(cell)
    inp%cls(1)%sws = cell%w * inp%alat
    call asa_scf(inp, cell, res)
    print '(a,a,a,l1,a,i0,a,f8.4,a,es9.2)', '    ', name, ' converged=', res%converged, ' it=', res%iter, '  moment=', res%cls(1)%moment, '  qnet=', res%cls(1)%qnet
    call check(res%converged, name//' SCF converged')
    call check(abs(res%cls(1)%moment) > mlo .and. abs(res%cls(1)%moment) < mhi, name//' |moment| in literature range')
    call check_close(res%cls(1)%qnet, 0.0_dp, 1.0e-4_dp, name//' sphere neutral')
  end subroutine run
end program test_asa_mag
