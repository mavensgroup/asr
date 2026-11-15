!> probe_pfun: exact ASA potential function P^gamma_l(E) of self-consistent bcc Fe (majority spin) from the radial
!> logarithmic derivative,  P^0 = (w/S)^(2l+1) 2(2l+1) (D + l + 1)/(D - l),  P^gamma = P^0 / (1 - gamma P^0),
!> compared with the parametrisations used by the energy-dependent recursion:
!>   two-parameter   (E - C)/Delta
!>   third-order     (E - C + p (E - E_nu)^3)/Delta
!> on real energies E_nu - 1 ... E_nu + 1 Ry.
program probe_pfun
  use asr_kinds, only: dp
  use asr_asa
  use asr_atom, only: shell_t
  use asr_lmto_cell
  use asr_radial
  implicit none
  type(asa_input_t) :: oi
  type(asa_result_t) :: orr
  type(cell_t) :: cell
  real(dp), allocatable :: y(:, :, :)
  real(dp) :: e, s, w, d, mm, p0, pg, p2, p3, eps, c, avw
  integer :: l, ie, n, sp

  allocate (oi%cls(1))
  oi%cls(1)%z = 26.0_dp; oi%cls(1)%nval = 8; oi%cls(1)%core = core_shells_for([4, 4, 3, 4]); oi%cls(1)%nvl = [4, 4, 3]
  oi%cls(1)%atom_shells = [oi%cls(1)%core, shell_t(4, 0, 2.0_dp), shell_t(3, 2, 6.0_dp)]
  oi%alat = 2.8665_dp / 0.529177210903_dp; oi%nk = 12; oi%maxit = 200; oi%nspin = 2; oi%zeta0 = 0.4_dp; oi%tol = 1.0e-6_dp
  cell%nat = 1; cell%avec = reshape([-0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp], [3, 3])
  allocate (cell%tau(3, 1), cell%isp(1)); cell%tau = 0.0_dp; cell%isp = 1
  call setup_cell(cell); oi%cls(1)%sws = cell%w * oi%alat
  call asa_scf(oi, cell, orr)
  print '(a,f9.4)', 'bcc Fe moment ', orr%cls(1)%moment
  c = c_ry; avw = cell%w * oi%alat
  n = orr%cls(1)%mesh%n; s = orr%cls(1)%mesh%r(n); w = avw
  allocate (y(2, 0:0, n))
  sp = 1; if (orr%cls(1)%moment < 0.0_dp) sp = 2                ! majority spin
  do l = 0, 2
    associate (pp => orr%cls(1)%pp(l, sp))
      print '(a,i1,a,f9.4,a,f9.4,a,f9.5,a,f9.4,a,f9.5)', 'l=', l, '  E_nu=', pp%enu, '  C=', pp%c, '  sqrt(Delta)=', pp%sqrdel, &
          '  p=', pp%p, '  gamma=', pp%gamma
      print '(a)', '   E-E_nu    exact P^g    2-param     3rd-order   (relative errors)'
      do ie = -10, 10
        eps = 0.1_dp * ie; e = pp%enu + eps
        call integrate_regular(orr%cls(1)%mesh, orr%cls(1)%v(:, sp), l, e, 26.0_dp, c, 0, n, y)
        mm = 1.0_dp + (e - orr%cls(1)%v(n, sp)) / (c * c)
        d = mm * y(2, 0, n) * s / y(1, 0, n)
        p0 = (w / s)**(2 * l + 1) * 2.0_dp * (2 * l + 1) * (d + l + 1) / (d - l)
        pg = p0 / (1.0_dp - pp%gamma * p0)
        p2 = (e - pp%c) / pp%sqrdel**2
        p3 = (e - pp%c + pp%p * eps**3) / pp%sqrdel**2
        print '(f8.2,3es13.4,2f10.4)', eps, pg, p2, p3, (p2 - pg) / max(abs(pg), 1.0e-3_dp), (p3 - pg) / max(abs(pg), 1.0e-3_dp)
      end do
    end associate
  end do
end program probe_pfun
