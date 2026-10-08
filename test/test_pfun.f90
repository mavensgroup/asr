!> Exact potential function at complex energy (asr_pfun) for a bcc Fe ASA potential:
!> (1) on the real axis P^0 equals the real-axis radial solver (asr_radial) value;
!> (2) at E_nu, alpha = gamma: P^gamma = (E_nu - C)/Delta, dP/dz = 1/Delta (definition of the gamma parameters);
!> (5) the single-site term of the density is real on the real axis and carries no charge;
!> (3) mu^2 = dP/dz is the complex derivative (Cauchy-Riemann: same along the imaginary direction);
!> (4) the sphere-normalised phi^2 integrates to one.
program test_pfun
  use asr_kinds, only: dp
  use asr_asa
  use asr_atom, only: shell_t
  use asr_lmto_cell
  use asr_radial
  use asr_pfun
  use asr_testing
  implicit none
  type(asa_input_t) :: oi
  type(asa_result_t) :: orr
  type(cell_t) :: cell
  real(dp), allocatable :: y(:, :, :)
  complex(dp), allocatable :: phi2(:)
  complex(dp) :: pf, mu2, lam, pfp, pfm, du, dummy1, dummy2, z
  real(dp) :: w, s, d, mm, p0r, e, e1, e2, e3, e4, e5, e6, e7, nrm
  integer :: l, n, i
  print '(a)', 'test_pfun'
  allocate (oi%cls(1))
  oi%cls(1)%z = 26.0_dp; oi%cls(1)%nval = 8; oi%cls(1)%core = core_shells_for([4, 4, 3, 4]); oi%cls(1)%nvl = [4, 4, 3]
  oi%cls(1)%atom_shells = [oi%cls(1)%core, shell_t(4, 0, 2.0_dp), shell_t(3, 2, 6.0_dp)]
  oi%alat = 2.8665_dp / 0.529177210903_dp; oi%nk = 8; oi%maxit = 4; oi%nspin = 1
  cell%nat = 1; cell%avec = reshape([-0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp], [3, 3])
  allocate (cell%tau(3, 1), cell%isp(1)); cell%tau = 0.0_dp; cell%isp = 1
  call setup_cell(cell); oi%cls(1)%sws = cell%w * oi%alat
  call asa_scf(oi, cell, orr)
  w = cell%w * oi%alat
  associate (ms => orr%cls(1)%mesh, v => orr%cls(1)%v(:, 1))
    n = ms%n; s = ms%r(n)
    allocate (y(2, 0:0, n), phi2(n))
    e1 = 0.0_dp; e2 = 0.0_dp; e3 = 0.0_dp; e4 = 0.0_dp; e5 = 0.0_dp; e6 = 0.0_dp; e7 = 0.0_dp
    do l = 0, 2
      associate (pp => orr%cls(1)%pp(l, 1))
        ! (1) real axis, alpha = 0
        e = pp%enu + 0.13_dp
        call integrate_regular(ms, v, l, e, 26.0_dp, c_ry, 0, n, y)
        mm = 1.0_dp + (e - v(n)) / c_ry**2
        d = mm * y(2, 0, n) * s / y(1, 0, n)
        p0r = 2.0_dp * (2 * l + 1) * (w / s)**(2 * l + 1) * (d + l + 1) / (d - l)
        call pfun_exact(ms, v, l, cmplx(e, 0.0_dp, dp), 26.0_dp, c_ry, w, 0.0_dp, pf, mu2, lam)
        e1 = max(e1, abs(pf - p0r) / abs(p0r))
        ! (2) gamma parameters at E_nu
        call pfun_exact(ms, v, l, cmplx(pp%enu, 0.0_dp, dp), 26.0_dp, c_ry, w, pp%gamma, pf, mu2, lam, phi2)
        e2 = max(e2, abs(pf - (pp%enu - pp%c) / pp%sqrdel**2) / abs((pp%enu - pp%c) / pp%sqrdel**2))
        e3 = max(e3, abs(mu2 * pp%sqrdel**2 - 1.0_dp))     ! lambda = C'/C is not zero at E_nu (only the 2-parameter -P''/2P' was)
        ! (5) single-site term: no charge, real on the real axis
        block
          complex(dp) :: ts(n), ph2(n), a1, a2, a3
          real(dp) :: tq
          integer :: k
          call pfun_exact(ms, v, l, cmplx(pp%enu + 0.07_dp, 0.0_dp, dp), 26.0_dp, c_ry, w, 0.02_dp, a1, a2, a3, ph2, ts)
          e6 = max(e6, maxval(abs(aimag(ts))) / maxval(abs(ts)))
          call pfun_exact(ms, v, l, cmplx(pp%enu - 0.1_dp, 0.2_dp, dp), 26.0_dp, c_ry, w, 0.02_dp, a1, a2, a3, ph2, ts)
          tq = 0.0_dp
          do k = 1, n
            tq = tq + merge(1.0_dp, merge(4.0_dp, 2.0_dp, mod(k, 2) == 0), k == 1 .or. k == n) * abs(ts(k)) * ms%drdt(k)
          end do
          e7 = max(e7, abs(sum(ts * ms%drdt * [(merge(1.0_dp, merge(4.0_dp, 2.0_dp, mod(k, 2) == 0), k == 1 .or. k == n), k=1, n)])) / tq)
        end block
        ! (4) normalisation
        nrm = 0.0_dp
        do i = 1, n
          nrm = nrm + merge(1.0_dp, merge(4.0_dp, 2.0_dp, mod(i, 2) == 0), i == 1 .or. i == n) * real(phi2(i)) * ms%drdt(i)
        end do
        e5 = max(e5, abs(nrm / 3.0_dp - 1.0_dp))
        ! (3) complex derivative along the imaginary axis
        z = cmplx(pp%enu - 0.05_dp, 0.08_dp, dp)
        call pfun_exact(ms, v, l, z, 26.0_dp, c_ry, w, 0.02_dp, pf, mu2, lam)
        call pfun_exact(ms, v, l, z + (0.0_dp, 1.0e-4_dp), 26.0_dp, c_ry, w, 0.02_dp, pfp, dummy1, dummy2)
        call pfun_exact(ms, v, l, z - (0.0_dp, 1.0e-4_dp), 26.0_dp, c_ry, w, 0.02_dp, pfm, dummy1, dummy2)
        du = (pfp - pfm) / (0.0_dp, 2.0e-4_dp)
        e4 = max(e4, abs(du - mu2) / abs(mu2))
      end associate
    end do
  end associate
  print '(a,es9.2,a,es9.2,a,es9.2,a,es9.2,a,es9.2)', '    real-axis P0 ', e1, '  P^g(E_nu) ', e2, '  mu2,lambda at E_nu ', e3, '  CR ', e4, '  norm ', e5
  call check(e1 < 1.0e-8_dp, 'complex solver == real-axis radial solver on the real axis')
  ! potpar's C, Delta rest on the Wronskian-normalisation identity, satisfied only to the Simpson accuracy of the norm (~1e-4);
  ! the exact route does not use it, so the two agree to that level (same 1e-4 in probe_pfun with converged potentials)
  call check(e2 < 3.0e-4_dp, 'P^gamma(E_nu) = (E_nu - C)/Delta (to the potential-parameter accuracy)')
  call check(e3 < 3.0e-4_dp, 'dP^gamma/dz = 1/Delta at E_nu (to the potential-parameter accuracy)')
  call check(e4 < 1.0e-6_dp, 'mu^2 is the complex derivative (Cauchy-Riemann)')
  call check(e5 < 1.0e-8_dp, 'phi^2 normalised in the sphere')
  print '(a,es9.2,a,es9.2)', '    single-site term: Im on real axis ', e6, '   charge ', e7
  call check(e6 < 1.0e-12_dp, 'single-site term is real on the real axis')
  call check(e7 < 1.0e-10_dp, 'single-site term carries no charge')
  call report()
end program test_pfun
