!> libxc binding: unpolarised == polarised at zeta=0; polarised exchange scaling; sign of spin splitting.
program test_xc
  use asr_kinds, only: dp
  use asr_xc
  use asr_testing
  implicit none
  type(xc_lda_t) :: f
  real(dp) :: rho(3), e0(3), v0(3), ep(3), vu(3), vd(3), rup(3), rdn(3)
  print '(a)', 'test_xc'
  rho = [0.001_dp, 0.05_dp, 1.0_dp]
  call xc_init(f, 'vbh')
  call xc_eval(f, rho, e0, v0)
  rup = 0.5_dp * rho; rdn = 0.5_dp * rho
  call xc_eval_pol(f, rup, rdn, ep, vu, vd)
  call check_close(maxval(abs(ep - e0)), 0.0_dp, 1.0e-10_dp, 'polarised exc at zeta=0 == unpolarised')
  call check_close(maxval(abs(vu - v0)), 0.0_dp, 1.0e-10_dp, 'polarised vxc at zeta=0 == unpolarised')
  rup = 0.7_dp * rho; rdn = 0.3_dp * rho
  call xc_eval_pol(f, rup, rdn, ep, vu, vd)
  call check(all(vu < vd), 'majority-spin xc potential is deeper (more attractive)')
  call xc_finalize(f)
  call report()
end program test_xc
