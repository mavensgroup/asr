!> probe_pfun_cmp: exact (asr_pfun) vs parametrised (asr_endep::ed_functions) P^alpha, mu^2, lambda at complex energies.
program probe_pfun_cmp
  use asr_kinds, only: dp
  use asr_asa
  use asr_atom, only: shell_t
  use asr_lmto_cell
  use asr_radial
  use asr_pfun
  use asr_endep
  implicit none
  type(asa_input_t) :: oi
  type(asa_result_t) :: orr
  type(cell_t) :: cell
  type(ed_params_t) :: sp
  complex(dp) :: pf, mu2, lam, z, pp9(9), mm9(9), ll9(9)
  real(dp) :: w, al(0:2), al9(9)
  integer :: l, k, i
  real(dp), parameter :: dre(4) = [-0.6_dp, -0.2_dp, 0.0_dp, 0.2_dp], dim(4) = [0.5_dp, 0.2_dp, 0.05_dp, 0.1_dp]
  integer, parameter :: lo9(9) = [0, 1, 1, 1, 2, 2, 2, 2, 2], rep(0:2) = [1, 2, 5]
  allocate (oi%cls(1))
  oi%cls(1)%z = 26.0_dp; oi%cls(1)%nval = 8; oi%cls(1)%core = core_shells_for([4, 4, 3, 4]); oi%cls(1)%nvl = [4, 4, 3]
  oi%cls(1)%atom_shells = [oi%cls(1)%core, shell_t(4, 0, 2.0_dp), shell_t(3, 2, 6.0_dp)]
  oi%alat = 2.8665_dp / 0.529177210903_dp; oi%nk = 8; oi%maxit = 30; oi%nspin = 1
  cell%nat = 1; cell%avec = reshape([-0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp], [3, 3])
  allocate (cell%tau(3, 1), cell%isp(1)); cell%tau = 0.0_dp; cell%isp = 1
  call setup_cell(cell); oi%cls(1)%sws = cell%w * oi%alat
  call asa_scf(oi, cell, orr)
  w = cell%w * oi%alat
  al = [0.35_dp, 0.0571667_dp, 0.0168070_dp]          ! hard-core sigma 0.7
  do i = 1, 9
    al9(i) = al(lo9(i)); sp%c(i) = orr%cls(1)%pp(lo9(i), 1)%c; sp%del(i) = orr%cls(1)%pp(lo9(i), 1)%sqrdel**2
    sp%gam(i) = orr%cls(1)%pp(lo9(i), 1)%gamma; sp%enu(i) = orr%cls(1)%pp(lo9(i), 1)%enu
  end do
  do k = 1, 4
    z = cmplx(orr%cls(1)%pp(2, 1)%enu + dre(k), dim(k), dp)
    call ed_functions(z, sp, al9, pp9, mm9, ll9)
    do l = 0, 2
      call pfun_exact(orr%cls(1)%mesh, orr%cls(1)%v(:, 1), l, z, 26.0_dp, c_ry, w, al(l), pf, mu2, lam)
      print '(a,2f7.3,a,i1,a,2es11.3,a,2es11.3,a,2es11.3,a,2es11.3,a,2es11.3,a,2es11.3)', 'z=', z, ' l=', l, &
          ' P ex', pf, ' par', pp9(rep(l)), '  mu2 ex', mu2, ' par', mm9(rep(l)), '  lam ex', lam, ' par', ll9(rep(l))
    end do
  end do
  ! single-site (S = 0) count: G = lambda + mu^2 / P over the closed contour E_nu - 1.2 .. E_nu + 0.3 (one d level at C expected)
  block
    use asr_contour, only: contour_nodes_closed
    use asr_kinds, only: pi
    integer, parameter :: na = 60, nv = 8
    complex(dp) :: zq(na + nv), wq(na + nv), gsum(0:2), gpar(0:2), pfq
    real(dp) :: eb, ef
    integer :: j
    eb = orr%cls(1)%pp(2, 1)%enu - 1.2_dp; ef = orr%cls(1)%pp(2, 1)%enu + 0.3_dp
    call contour_nodes_closed(eb, ef, 0.01_dp, na, nv, zq, wq)
    gsum = 0.0_dp; gpar = 0.0_dp
    do j = 1, na + nv
      call ed_functions(zq(j), sp, al9, pp9, mm9, ll9)
      do l = 0, 2
        call pfun_exact(orr%cls(1)%mesh, orr%cls(1)%v(:, 1), l, zq(j), 26.0_dp, c_ry, w, al(l), pf, mu2, lam)
        gsum(l) = gsum(l) + wq(j) * (lam + mu2 / pf)
        gpar(l) = gpar(l) + wq(j) * (ll9(rep(l)) + mm9(rep(l)) / pp9(rep(l)))
      end do
    end do
    block
      complex(dp), allocatable :: ph(:)
      real(dp), allocatable :: rho(:)
      real(dp) :: q
      integer :: nn, ii
      nn = orr%cls(1)%mesh%n
      allocate (ph(nn), rho(nn)); rho = 0.0_dp
      do j = 1, na + nv
        call pfun_exact(orr%cls(1)%mesh, orr%cls(1)%v(:, 1), 2, zq(j), 26.0_dp, c_ry, w, al(2), pfq, mu2, lam, ph)
        rho = rho - aimag(wq(j) * (lam + mu2 / pfq) * ph) / pi
      end do
      q = 0.0_dp
      do ii = 1, nn
        q = q + merge(1.0_dp, merge(4.0_dp, 2.0_dp, mod(ii, 2) == 0), ii == 1 .or. ii == nn) * rho(ii) * orr%cls(1)%mesh%drdt(ii)
      end do
      block
        complex(dp), allocatable :: ts(:), phd(:)
        real(dp), allocatable :: rho2(:), ref(:)
        complex(dp) :: pa, pb, pm, d1, d2
        real(dp) :: ea, eb2, em, e1t, e2t
        integer :: it2
        allocate (ts(nn), phd(nn), rho2(nn), ref(nn)); rho2 = 0.0_dp
        do j = 1, na + nv
          call pfun_exact(orr%cls(1)%mesh, orr%cls(1)%v(:, 1), 2, zq(j), 26.0_dp, c_ry, w, al(2), pfq, mu2, lam, ph, ts)
          rho2 = rho2 - aimag(wq(j) * ((lam + mu2 / pfq) * ph + ts)) / pi
        end do
        ! reference: phi^2 at the real zero of P_d (the single bound level in the window)
        ea = orr%cls(1)%pp(2, 1)%c - 0.1_dp; eb2 = orr%cls(1)%pp(2, 1)%c + 0.1_dp
        do it2 = 1, 80
          em = 0.5_dp * (ea + eb2)
          call pfun_exact(orr%cls(1)%mesh, orr%cls(1)%v(:, 1), 2, cmplx(ea, 0.0_dp, dp), 26.0_dp, c_ry, w, al(2), pa, d1, d2)
          call pfun_exact(orr%cls(1)%mesh, orr%cls(1)%v(:, 1), 2, cmplx(em, 0.0_dp, dp), 26.0_dp, c_ry, w, al(2), pm, d1, d2)
          if (real(pa) * real(pm) <= 0.0_dp) then; eb2 = em; else; ea = em; end if
        end do
        call pfun_exact(orr%cls(1)%mesh, orr%cls(1)%v(:, 1), 2, cmplx(em, 0.0_dp, dp), 26.0_dp, c_ry, w, al(2), pa, d1, d2, phd, ts)
        ref = real(phd)
        e1t = maxval(abs(rho - ref)) / maxval(abs(ref)); e2t = maxval(abs(rho2 - ref)) / maxval(abs(ref))
        print '(a,f9.5,a,es10.2,a,es10.2,a,es10.2)', 'single-site d level', em, ': max|rho - phi^2(E_d)|/max: without T', e1t, &
            '   with T', e2t, '   max|Im T(E_d)|', maxval(abs(aimag(ts)))
      end block
      print '(a,f10.5,a,es10.2)', 'single-site d: integral of -1/pi Im sum w G phi^2 =', q / 3.0_dp, '   max rho', maxval(abs(rho))
    end block
    print '(a,3f10.4)', 'single-site count exact s p d:', -aimag(gsum) / pi
    print '(a,3f10.4)', 'single-site count param s p d:', -aimag(gpar) / pi
    print '(a,3f9.4)', 'C (s p d):', (orr%cls(1)%pp(l, 1)%c, l=0, 2)
    print '(a,2f9.4)', 'window:', eb, ef
  end block
end program probe_pfun_cmp
