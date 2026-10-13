!> libxc binding: unpolarised == polarised at zeta=0; polarised exchange scaling; sign of spin splitting.
program test_xc
  use asr_kinds, only: dp
  use asr_xc
  use asr_testing
  use asr_radial, only: radmesh_t, make_radmesh
  use asr_atom, only: cumint
  implicit none
  type(xc_lda_t) :: f
  real(dp) :: rho(3), e0(3), v0(3), ep(3), vu(3), vd(3), rup(3), rdn(3)
  type(radmesh_t) :: m
  real(dp), allocatable :: n0(:), dn(:), na(:), nb(:), ex(:), vx(:), vx2(:), w(:), cum(:)
  real(dp) :: ep_, em_, dfun, lin, h
  integer :: k
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
  ! PBE (GGA) on a radial density: the potential is the functional derivative,
  ! dE[n + h dn]/dh = int 4 pi r^2 v dn  (tests the divergence term -(1/r^2) d/dr[r^2 2 vsigma n'])
  call make_radmesh(10.0_dp, 12.0_dp, 0.012_dp, m)
  allocate (n0(m%n), dn(m%n), na(m%n), nb(m%n), ex(m%n), vx(m%n), vx2(m%n), w(m%n), cum(m%n))
  n0 = 2.0_dp * exp(-1.3_dp * m%r) + 0.3_dp * exp(-0.4_dp * m%r**2)
  dn = m%r**2 * exp(-0.8_dp * m%r**2)
  w = 4.0_dp * 3.14159265358979324_dp * m%r**2
  h = 1.0e-4_dp
  call xc_init(f, 'pbe')
  call xc_eval(f, n0 + h * dn, ex, vx, m%r); na = w * (n0 + h * dn) * ex; call cumint(m, na, cum); ep_ = cum(m%n)
  call xc_eval(f, n0 - h * dn, ex, vx, m%r); na = w * (n0 - h * dn) * ex; call cumint(m, na, cum); em_ = cum(m%n)
  dfun = (ep_ - em_) / (2.0_dp * h)
  call xc_eval(f, n0, ex, vx, m%r); na = w * vx * dn; call cumint(m, na, cum); lin = cum(m%n)
  print '(a,2es16.8)', '    PBE  dE/dh, int v dn: ', dfun, lin
  call check_close(lin, dfun, 1.0e-5_dp * abs(dfun), 'PBE potential == functional derivative (unpolarised)')
  ! polarised: perturb the up density only
  na = 0.6_dp * n0; nb = 0.4_dp * n0
  call xc_eval_pol(f, na + h * dn, nb, ex, vx, vx2, m%r); w = 4.0_dp * 3.14159265358979324_dp * m%r**2
  cum = w * (na + h * dn + nb) * ex; call cumint(m, cum, ex); ep_ = ex(m%n)
  call xc_eval_pol(f, na - h * dn, nb, ex, vx, vx2, m%r)
  cum = w * (na - h * dn + nb) * ex; call cumint(m, cum, ex); em_ = ex(m%n)
  dfun = (ep_ - em_) / (2.0_dp * h)
  call xc_eval_pol(f, na, nb, ex, vx, vx2, m%r); cum = w * vx * dn; call cumint(m, cum, ex); lin = ex(m%n)
  print '(a,2es16.8)', '    PBE pol dE/dh_up, int v_up dn: ', dfun, lin
  call check_close(lin, dfun, 1.0e-5_dp * abs(dfun), 'PBE spin-up potential == functional derivative (polarised)')
  call xc_eval_pol(f, 0.5_dp * n0, 0.5_dp * n0, ex, vx, vx2, m%r)
  call xc_eval(f, n0, cum, w, m%r)
  call check_close(maxval(abs(vx(2:) - w(2:)) / (1.0_dp + abs(w(2:)))), 0.0_dp, 1.0e-8_dp, 'PBE polarised vxc at zeta=0 == unpolarised')
  call xc_finalize(f)
  call report()
end program test_xc
