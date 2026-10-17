!> Self-consistent TB-LMTO-ASA for ordered cells: any number of atoms, any number of inequivalent
!> classes (own sphere radius, Z, core, valence), collinear spin polarisation.  Own radial solver,
!> relaxed core in the sphere (Hankel tail at S), neutral-sphere Hartree (+ optional Madelung shift),
!> libxc, first-order ASA bands, moments -> charge, E_nu to band centre, potential mixing.
!> Every parameter is generated here; nothing tabulated.
module asr_asa
  use asr_kinds, only: dp, pi
  use asr_radial
  use asr_xc
  use asr_atom, only: cumint, poisson_hartree, shell_t, solve_atom, atom_result_t
  use asr_potpar, only: potpar_l_t, potpar_at_enu, clamp_enu
  use asr_lmto_ham, only: potpar_t, to_alpha, o_alpha
  use asr_lmto_cell
  use asr_ewald, only: madelung_phi
  use asr_mixing, only: anderson_t, anderson_init, anderson_step
  use asr_tetra, only: tetra_weights, tetra_fermi
  implicit none
  private
  public :: asa_class_t, asa_input_t, asa_class_result_t, asa_result_t, asa_scf, core_shells_for, asa_dos, atomic_start, bands_moments

  type :: asa_class_t
    real(dp) :: z = 0.0_dp
    real(dp) :: sws = 0.0_dp                !! sphere radius (Bohr)
    integer :: nval = 0                     !! valence electrons per atom of this class
    integer :: nvl(0:2) = [4, 4, 3]         !! valence principal quantum numbers (s,p,d)
    type(shell_t), allocatable :: core(:)
    type(shell_t), allocatable :: atom_shells(:)   !! neutral free atom, all shells (starting potential)
    real(dp) :: spin_sign = 1.0_dp          !! sign of the starting spin polarisation (zeta0) of this class
  end type asa_class_t

  type :: asa_input_t
    type(asa_class_t), allocatable :: cls(:)
    real(dp) :: alat = 0.0_dp               !! lattice constant (Bohr)
    character(len=8) :: xc = 'vbh'
    real(dp) :: sigma = 0.7_dp              !! hard-core screening
    integer :: nspin = 1
    integer :: nk = 12
    character(len=8) :: bz = 'tetra'        !! BZ integration: 'tetra' (linear tetrahedra) | 'gauss'
    real(dp) :: smear = 0.01_dp             !! Gaussian smearing (Ry), used when bz = 'gauss'
    real(dp) :: beta = 0.3_dp               !! potential mixing
    integer :: maxit = 100
    real(dp) :: tol = 1.0e-6_dp
    logical :: relativistic = .true.
    real(dp) :: zeta0 = 0.0_dp              !! initial spin polarisation of the atomic start (nspin = 2)
    integer :: mix_history = 8              !! Anderson mixing history
    real(dp) :: dvmax = 0.15_dp             !! safety clamp: max change of the potential per iteration, per point (Ry)
    real(dp) :: denumax = 0.08_dp           !! safety clamp: max change of E_nu per iteration (Ry)
    logical :: second_order = .false.       !! diagnostic: orthogonal H = E_nu + h - h o h (the recursion Hamiltonian) instead of the generalised problem
    logical :: trace = .false.              !! print per iteration, class and spin: E_nu, C, sqrt(Delta)_d, q0, V(S), Madelung shift
  end type asa_input_t

  type :: asa_class_result_t
    type(radmesh_t) :: mesh
    real(dp), allocatable :: v(:, :)        !! (n, spin) true potential
    real(dp), allocatable :: nrad(:, :)     !! (n, spin) radial charge 4 pi r^2 rho
    type(potpar_l_t) :: pp(0:2, 2)          !! gamma-rep parameters at converged E_nu, (l, spin)
    real(dp) :: q(0:2, 0:2, 2) = 0.0_dp     !! gamma-rep moments (k, l, spin) per atom
    real(dp) :: enu(0:2, 2) = 0.0_dp
    real(dp) :: moment = 0.0_dp             !! spin moment inside the sphere (mu_B per atom)
    real(dp) :: qnet = 0.0_dp               !! electrons in sphere - Z
  end type asa_class_result_t

  type :: asa_result_t
    type(asa_class_result_t), allocatable :: cls(:)
    real(dp) :: ef = 0.0_dp, dv = 0.0_dp
    real(dp) :: moment_cell = 0.0_dp
    integer :: iter = 0
    logical :: converged = .false.
    !! SCF history, rows 1:iter (dV in Ry; moment in mu_B and qnet in electrons per class)
    real(dp), allocatable :: hist_dv(:), hist_ef(:), hist_moment(:, :), hist_qnet(:, :)
  end type asa_result_t

contains

  function core_shells_for(nvl) result(sh)
    integer, intent(in) :: nvl(0:3)
    type(shell_t), allocatable :: sh(:)
    type(shell_t) :: tmp(40)
    integer :: l, n, k
    k = 0
    do l = 0, 3
      do n = l + 1, nvl(l) - 1
        k = k + 1; tmp(k) = shell_t(n, l, real(2 * (2 * l + 1), dp))
      end do
    end do
    allocate (sh(k)); sh = tmp(1:k)
  end function core_shells_for

  subroutine asa_scf(inp, cell, res, verbose)
    type(asa_input_t), intent(in) :: inp
    type(cell_t), intent(in) :: cell
    type(asa_result_t), intent(out) :: res
    logical, intent(in), optional :: verbose

    type(xc_lda_t) :: xf
    type(cluster_t), allocatable :: cls(:)
    type(potpar_t), allocatable :: pps(:, :)          ! (class, spin)
    type(potpar_l_t) :: ppl
    real(dp), allocatable :: alpha_sp(:, :), vnew(:, :), exc(:), vxc(:, :), rho(:, :), vh(:), ntot(:), tmp(:), cum(:)
    real(dp), allocatable :: y(:, :, :, :, :), yc(:, :, :), qat(:, :, :, :), qcl(:, :, :, :), enu(:, :, :), vmad(:)
    real(dp), allocatable :: nval(:, :, :), ncore(:, :, :), vold(:, :, :)
    integer, allocatable :: natc(:)
    real(dp), allocatable :: phimad(:, :), qatom(:), vmad_at(:), xin(:), xout(:), xnext(:), wmix(:)
    integer, allocatable :: offx(:)
    integer :: nx, ii
    type(anderson_t) :: mixer
    real(dp) :: c, avw, ef, nelec, dv, dvc, e, ratio, qcore_s, qval_s
    integer :: ncl, nsp, it, ic, s, l, is, k, n, ia, nod
    logical :: vb, core_ok
    type(asa_class_result_t), allocatable :: cr(:)

    vb = .false.; if (present(verbose)) vb = verbose
    ncl = size(inp%cls); nsp = inp%nspin
    c = c_nonrel; if (inp%relativistic) c = c_ry
    avw = cell%w * inp%alat
    allocate (cr(ncl), alpha_sp(0:2, ncl), pps(ncl, nsp), natc(ncl), enu(0:2, nsp, ncl), vmad(ncl))
    natc = 0
    do ia = 1, cell%nat; natc(cell%isp(ia)) = natc(cell%isp(ia)) + 1; end do
    nelec = 0.0_dp
    do ia = 1, cell%nat; nelec = nelec + real(inp%cls(cell%isp(ia))%nval, dp); end do
    call xc_init(xf, inp%xc)
    allocate (res%hist_dv(inp%maxit), res%hist_ef(inp%maxit), res%hist_moment(inp%maxit, ncl), res%hist_qnet(inp%maxit, ncl))
    res%hist_dv = 0.0_dp; res%hist_ef = 0.0_dp; res%hist_moment = 0.0_dp; res%hist_qnet = 0.0_dp
    vmad = 0.0_dp
    allocate (phimad(cell%nat, cell%nat), qatom(cell%nat), vmad_at(cell%nat))
    call madelung_phi(cell%avec, cell%tau, inp%alat, phimad)

    do ic = 1, ncl
      call make_radmesh(inp%cls(ic)%z, inp%cls(ic)%sws, 0.012_dp, cr(ic)%mesh, c)
      n = cr(ic)%mesh%n
      allocate (cr(ic)%v(n, nsp), cr(ic)%nrad(n, nsp))
      alpha_sp(:, ic) = alpha_from_sigma(inp%sigma, inp%cls(ic)%sws / avw)
    end do
    allocate (cls(cell%nat))
    do ia = 1, cell%nat
      call build_cluster(cell, ia, 4.0_dp * cell%w, alpha_sp, cls(ia))
    end do
    allocate (offx(ncl + 1)); offx(1) = 0
    do ic = 1, ncl; offx(ic + 1) = offx(ic) + cr(ic)%mesh%n * nsp + 3 * nsp; end do
    nx = offx(ncl + 1)
    allocate (xin(nx), xout(nx), xnext(nx))
    call anderson_init(mixer, nx, inp%mix_history)

    ! ---- starting potentials: superposed free atoms, spin-split by zeta0
    do ic = 1, ncl
      n = cr(ic)%mesh%n
      allocate (rho(n, 2), exc(n), vxc(n, 2))
      call atomic_start(inp, ic, cr(ic)%mesh, c, xf, nsp, cr(ic)%v, enu(:, :, ic))
      deallocate (rho, exc, vxc)
    end do

    do it = 1, inp%maxit
      do ic = 1, ncl
        n = cr(ic)%mesh%n
        do s = 1, nsp
          xin(offx(ic) + (s - 1) * n + 2:offx(ic) + s * n) = cr(ic)%v(2:n, s) + 2.0_dp * inp%cls(ic)%z / cr(ic)%mesh%r(2:n)
          xin(offx(ic) + (s - 1) * n + 1) = xin(offx(ic) + (s - 1) * n + 2)
          xin(offx(ic) + n * nsp + 3 * (s - 1) + 1:offx(ic) + n * nsp + 3 * s) = enu(:, s, ic)
        end do
      end do
      ! (1) potential parameters at E_nu, per class / spin / l
      do ic = 1, ncl
        n = cr(ic)%mesh%n
        if (.not. allocated(y)) allocate (y(2, 0:2, maxval(cr(:)%mesh%n), 0:2, 2 * ncl))
        do s = 1, nsp
          do l = 0, 2
            call potpar_at_enu(cr(ic)%mesh, cr(ic)%v(:, s), l, enu(l, s, ic), inp%cls(ic)%z, c, avw, ppl, &
                               y(:, :, 1:n, l, 2 * (ic - 1) + s))
            cr(ic)%pp(l, s) = ppl
            pps(ic, s)%enu(l) = ppl%enu; pps(ic, s)%c(l) = ppl%c; pps(ic, s)%sqrdel(l) = ppl%sqrdel
            pps(ic, s)%p(l) = ppl%p; pps(ic, s)%gamma(l) = ppl%gamma; pps(ic, s)%alpha_rep(l) = ppl%gamma
          end do
        end do
      end do
      ! (2) bands, common Fermi level, moments per atom -> per class
      if (allocated(qat)) deallocate (qat)
      allocate (qat(0:2, 0:2, cell%nat, nsp))
      call bands_moments(cell, cls, pps, alpha_sp, inp%nk, nsp, nelec, inp%smear, inp%bz, ef, qat, inp%second_order)
      do ic = 1, ncl
        cr(ic)%q = 0.0_dp
        do ia = 1, cell%nat
          if (cell%isp(ia) == ic) cr(ic)%q(:, :, 1:nsp) = cr(ic)%q(:, :, 1:nsp) + qat(:, :, ia, 1:nsp) / real(natc(ic), dp)
        end do
      end do
      ! (3) charge densities per class / spin (valence from moments, core relaxed in V_s)
      if (.not. allocated(vnew)) allocate (vnew(maxval(cr(:)%mesh%n), nsp))
      dv = 0.0_dp
      do ic = 1, ncl
        n = cr(ic)%mesh%n
        allocate (nval(n, nsp, 1), ncore(n, nsp, 1), tmp(n), cum(n), yc(2, 0:0, n))
        nval = 0.0_dp; ncore = 0.0_dp
        do s = 1, nsp
          do l = 0, 2
            k = 2 * (ic - 1) + s
            nval(:, s, 1) = nval(:, s, 1) &
              + cr(ic)%q(0, l, s) * (y(1, 0, 1:n, l, k)**2 + y(2, 0, 1:n, l, k)**2 / c**2) &
              + 2.0_dp * cr(ic)%q(1, l, s) * (y(1, 0, 1:n, l, k) * y(1, 1, 1:n, l, k) + y(2, 0, 1:n, l, k) * y(2, 1, 1:n, l, k) / c**2) &
              + cr(ic)%q(2, l, s) * (y(1, 1, 1:n, l, k)**2 + y(2, 1, 1:n, l, k)**2 / c**2)
          end do
          do is = 1, size(inp%cls(ic)%core)
            nod = inp%cls(ic)%core(is)%n - inp%cls(ic)%core(is)%l - 1
            call bound_state(cr(ic)%mesh, cr(ic)%v(:, s), inp%cls(ic)%core(is)%l, nod, inp%cls(ic)%z, c, n, &
                             -1.6_dp * inp%cls(ic)%z**2 - 20.0_dp, cr(ic)%v(n, s) - 1.0e-6_dp, e, yc, core_ok)
            if (.not. core_ok .and. vb) print '(a,i0,a,i0)', '    warning: core state not bracketed (n,l)=', inp%cls(ic)%core(is)%n, ',', inp%cls(ic)%core(is)%l
            ncore(:, s, 1) = ncore(:, s, 1) + inp%cls(ic)%core(is)%occ / real(nsp, dp) * (yc(1, 0, :)**2 + yc(2, 0, :)**2 / c**2)
          end do
        end do
        cr(ic)%nrad = nval(:, :, 1) + ncore(:, :, 1)
        tmp = sum(cr(ic)%nrad, dim=2); call cumint(cr(ic)%mesh, tmp, cum)
        cr(ic)%qnet = cum(n) - inp%cls(ic)%z
        if (nsp == 2) then
          tmp = cr(ic)%nrad(:, 1); call cumint(cr(ic)%mesh, tmp, cum); qval_s = cum(n)
          tmp = cr(ic)%nrad(:, 2); call cumint(cr(ic)%mesh, tmp, cum); qcore_s = cum(n)
          cr(ic)%moment = qval_s - qcore_s
        end if
        deallocate (nval, ncore, tmp, cum, yc)
      end do
      ! (4) Madelung shifts from the net sphere charges (electrons minus Z): V_i = 2 sum_j Phi_ij qnet_j (Ry)
      do ia = 1, cell%nat; qatom(ia) = cr(cell%isp(ia))%qnet; end do
      vmad_at = 2.0_dp * matmul(phimad, qatom)
      vmad = 0.0_dp
      do ia = 1, cell%nat; vmad(cell%isp(ia)) = vmad(cell%isp(ia)) + vmad_at(ia) / real(natc(cell%isp(ia)), dp); end do
      ! (5) Hartree (neutral-sphere form) + Madelung shift + xc -> new potential; mix; E_nu update
      do ic = 1, ncl
        n = cr(ic)%mesh%n
        allocate (ntot(n), vh(n), exc(n), vxc(n, 2), rho(n, 2))
        ntot = sum(cr(ic)%nrad, dim=2)
        call poisson_hartree(cr(ic)%mesh, ntot, vh)
        do k = 2, n; rho(k, :) = cr(ic)%nrad(k, :) / (4.0_dp * pi * cr(ic)%mesh%r(k)**2); end do
        rho(1, :) = rho(2, :)
        if (nsp == 1) then
          call xc_eval(xf, rho(:, 1), exc, vxc(:, 1), cr(ic)%mesh%r); vxc(:, 2) = vxc(:, 1)
        else
          call xc_eval_pol(xf, rho(:, 1), rho(:, 2), exc, vxc(:, 1), vxc(:, 2), cr(ic)%mesh%r)
        end if
        do s = 1, nsp
          vnew(2:n, s) = -2.0_dp * inp%cls(ic)%z / cr(ic)%mesh%r(2:n) + vh(2:n) + vxc(2:n, s) + vmad(ic)
          vnew(1, s) = vnew(2, s)
          dvc = maxval(abs(vnew(2:n, s) - cr(ic)%v(2:n, s)) * cr(ic)%mesh%r(2:n) / (1.0_dp + cr(ic)%mesh%r(2:n)))
          dv = max(dv, dvc)
          if (inp%trace) print '(a,i4,i3,i2,a,2f10.5,a,f10.5,a,i0)', 'trace_v', it, ic, s, ' V_in(S) V_out(S)', cr(ic)%v(n, s), vnew(n, s), &
            '  max |dV|', dvc, ' at i ', maxloc(abs(vnew(2:n, s) - cr(ic)%v(2:n, s)) * cr(ic)%mesh%r(2:n) / (1.0_dp + cr(ic)%mesh%r(2:n)), 1) + 1
          xout(offx(ic) + (s - 1) * n + 2:offx(ic) + s * n) = vnew(2:n, s) + 2.0_dp * inp%cls(ic)%z / cr(ic)%mesh%r(2:n)
          xout(offx(ic) + (s - 1) * n + 1) = xout(offx(ic) + (s - 1) * n + 2)
          do l = 0, 2
            ratio = cr(ic)%q(1, l, s) / max(cr(ic)%q(0, l, s), 1.0e-6_dp)
            xout(offx(ic) + n * nsp + 3 * (s - 1) + l + 1) = enu(l, s, ic) + 0.7_dp * ratio
          end do
        end do
        deallocate (ntot, vh, exc, vxc, rho)
      end do
      if (inp%trace) then
        do ic = 1, ncl
          do s = 1, nsp
            print '(a,i4,i3,i2,a,f9.5,a,3f9.5,a,3f9.5,a,f8.5,a,3f8.4,a,f9.5,a,f9.5)', 'trace', it, ic, s, ' EF', ef, &
              ' Enu', pps(ic, s)%enu, ' C', pps(ic, s)%c, ' sqD_d', pps(ic, s)%sqrdel(2), ' q0', cr(ic)%q(0, :, s), &
              ' V(S)', cr(ic)%v(cr(ic)%mesh%n, s), ' vmad', vmad(ic)
          end do
        end do
      end if
      ! (6) Anderson step on the full state (potentials and E_nu), with safety clamps
      ! Anderson metric = the dV measure (r/(1+r))^2 on the potentials, 1 on E_nu: unweighted, the log mesh's many points near the
      ! nucleus dominated the fit and the step near S went against its residual (bcc Fe: E_F to +1.5 Ry in the first iterations)
      if (.not. allocated(wmix)) then
        allocate (wmix(nx)); wmix = 1.0_dp
        do ic = 1, ncl
          n = cr(ic)%mesh%n
          do s = 1, nsp
            wmix(offx(ic) + (s - 1) * n + 1:offx(ic) + s * n) = (cr(ic)%mesh%r(1:n) / (1.0_dp + cr(ic)%mesh%r(1:n)))**2
          end do
        end do
      end if
      call anderson_step(mixer, xin, xout, inp%beta, xnext, wmix)
      do ic = 1, ncl
        n = cr(ic)%mesh%n
        do s = 1, nsp
          do ii = 1, n
            xnext(offx(ic) + (s - 1) * n + ii) = xin(offx(ic) + (s - 1) * n + ii) &
              + max(-inp%dvmax, min(inp%dvmax, xnext(offx(ic) + (s - 1) * n + ii) - xin(offx(ic) + (s - 1) * n + ii)))
          end do
          cr(ic)%v(2:n, s) = xnext(offx(ic) + (s - 1) * n + 2:offx(ic) + s * n) - 2.0_dp * inp%cls(ic)%z / cr(ic)%mesh%r(2:n)
          cr(ic)%v(1, s) = cr(ic)%v(2, s)
          do l = 0, 2
            ii = offx(ic) + n * nsp + 3 * (s - 1) + l + 1
            enu(l, s, ic) = xin(ii) + max(-inp%denumax, min(inp%denumax, xnext(ii) - xin(ii)))
            if (clamp_enu(cr(ic)%mesh, cr(ic)%v(:, s), l, inp%cls(ic)%nvl(l) - l - 1, inp%cls(ic)%z, c, enu(l, s, ic))) then
              if (vb) print '(a,i0,a,i0,a,i0,a,f9.5)', '    E_nu outside its band window: class ', ic, ' spin ', s, ' l ', l, ' -> ', enu(l, s, ic)
            end if
          end do
        end do
      end do
      res%iter = it; res%dv = dv; res%ef = ef
      res%hist_dv(it) = dv; res%hist_ef(it) = ef
      do ic = 1, ncl; res%hist_moment(it, ic) = cr(ic)%moment; res%hist_qnet(it, ic) = cr(ic)%qnet; end do
      if (vb) then
        write (*, '(a,i4,a,es10.3,a,f10.6)', advance='no') 'it', it, '  dV=', dv, '  EF=', ef
        write (*, '(a)', advance='no') '  m='
        do ic = 1, ncl; write (*, '(f9.4)', advance='no') cr(ic)%moment; end do
        write (*, '(a)', advance='no') '  qnet='
        do ic = 1, ncl; write (*, '(f9.5)', advance='no') cr(ic)%qnet; end do
        write (*, *)
      end if
      if (dv < inp%tol .and. it > 3) then; res%converged = .true.; exit; end if
    end do
    call xc_finalize(xf)
    ! final parameters at converged potential and E_nu
    do ic = 1, ncl
      n = cr(ic)%mesh%n
      do s = 1, nsp
        do l = 0, 2
          call potpar_at_enu(cr(ic)%mesh, cr(ic)%v(:, s), l, enu(l, s, ic), inp%cls(ic)%z, c, avw, cr(ic)%pp(l, s))
        end do
      end do
      cr(ic)%enu = enu(:, :, ic)
    end do
    res%cls = cr
    res%moment_cell = 0.0_dp
    do ia = 1, cell%nat; res%moment_cell = res%moment_cell + cr(cell%isp(ia))%moment; end do
  end subroutine asa_scf

  !> Free-atom (all electrons) potential interpolated onto the class sphere mesh; spin-split with zeta0;
  !> E_nu from the atomic valence eigenvalues.
  subroutine atomic_start(inp, ic, m, c, xf, nsp, v, enu)
    type(asa_input_t), intent(in) :: inp
    integer, intent(in) :: ic, nsp
    type(radmesh_t), intent(in) :: m
    real(dp), intent(in) :: c
    type(xc_lda_t), intent(in) :: xf
    real(dp), intent(out) :: v(:, :), enu(0:, :)
    type(atom_result_t) :: at
    real(dp), allocatable :: rho(:), e0(:), v0(:), eu(:), vu(:), vd(:)
    integer :: k, n, ish, s
    n = m%n
    enu = -0.4_dp
    if (.not. allocated(inp%cls(ic)%atom_shells)) then
      do k = 2, n
        v(k, 1) = -2.0_dp * inp%cls(ic)%z / m%r(k) * exp(-m%r(k) * inp%cls(ic)%z**(1.0_dp / 3.0_dp) / 2.0_dp) - 0.5_dp
      end do
      v(1, 1) = v(2, 1)
      if (size(v, 2) >= 2) v(:, 2) = v(:, 1)
      return
    end if
    call solve_atom(inp%cls(ic)%z, inp%cls(ic)%atom_shells, inp%xc, 30.0_dp, 0.01_dp, inp%relativistic, at)
    allocate (rho(n), e0(n), v0(n), eu(n), vu(n), vd(n))
    block
      real(dp), allocatable :: usm(:)
      allocate (usm(at%mesh%n))
      usm(2:) = at%v(2:) + 2.0_dp * inp%cls(ic)%z / at%mesh%r(2:); usm(1) = usm(2)
      do k = 2, n
        v(k, 1) = interp_cubic(at%mesh%r, usm, m%r(k)) - 2.0_dp * inp%cls(ic)%z / m%r(k)
        rho(k) = interp_cubic(at%mesh%r, at%rho, m%r(k))
      end do
    end block
    v(1, 1) = v(2, 1); rho(1) = rho(2)
    if (nsp == 2 .and. size(v, 2) >= 2) then
      call xc_eval(xf, rho, e0, v0, m%r)
      call xc_eval_pol(xf, 0.5_dp * (1.0_dp + inp%cls(ic)%spin_sign * inp%zeta0) * rho, &
                       0.5_dp * (1.0_dp - inp%cls(ic)%spin_sign * inp%zeta0) * rho, eu, vu, vd, m%r)
      v(:, 2) = v(:, 1) - v0 + vd
      v(:, 1) = v(:, 1) - v0 + vu
    else if (size(v, 2) >= 2) then
      v(:, 2) = v(:, 1)
    end if
    do ish = 1, size(inp%cls(ic)%atom_shells)
      if (inp%cls(ic)%atom_shells(ish)%l <= 2 .and. inp%cls(ic)%atom_shells(ish)%occ > 0.0_dp) then
        if (inp%cls(ic)%atom_shells(ish)%n == inp%cls(ic)%nvl(inp%cls(ic)%atom_shells(ish)%l)) then
          do s = 1, size(enu, 2); enu(inp%cls(ic)%atom_shells(ish)%l, s) = at%eps(ish); end do
        end if
      end if
    end do
  end subroutine atomic_start

  pure function interp_cubic(r, f, r0) result(y)
    real(dp), intent(in) :: r(:), f(:), r0
    real(dp) :: y, w
    integer :: i, j, k, lo, n
    n = size(r)
    lo = 1; i = n
    do while (i - lo > 1)
      k = (lo + i) / 2
      if (r(k) > r0) then; i = k; else; lo = k; end if
    end do
    lo = max(1, min(lo - 1, n - 3))
    y = 0.0_dp
    do j = lo, lo + 3
      w = 1.0_dp
      do k = lo, lo + 3
        if (k /= j) w = w * (r0 - r(k)) / (r(j) - r(k))
      end do
      y = y + w * f(j)
    end do
  end function interp_cubic

  pure integer function orb_l_of(i)
    integer, intent(in) :: i
    if (i == 1) then; orb_l_of = 0
    else if (i <= 4) then; orb_l_of = 1
    else; orb_l_of = 2
    end if
  end function orb_l_of

  !> First-order ASA generalised eigenproblem on the k-mesh, for every spin.
  !>   O = 1 + (o_i+o_j) h + h p h ;  H = diag(E_nu) + (w5_i+w5_j) h + h w6 h ;  H z = e O z
  !> eb(state, spin), state = nd*ik + ib (ik 0-based);  xm(0:2, 0:2, (state-1)*nat + atom, spin) holds
  !> per-(atom,l) sums of |z|^2, Re z*(hz), |hz|^2.
  subroutine bands_core(c, cls, pps, alpha_sp, nk, nsp, eb, xm, h2)
    type(cell_t), intent(in) :: c
    type(cluster_t), intent(in) :: cls(:)
    type(potpar_t), intent(in) :: pps(:, :)
    real(dp), intent(in) :: alpha_sp(0:, :)
    integer, intent(in) :: nk, nsp
    real(dp), allocatable, intent(out) :: eb(:, :), xm(:, :, :, :)
    logical, intent(in), optional :: h2           !! orthogonal second-order H; xm then holds |z|^2 (e-E_nu)^k
    real(dp), allocatable :: ov(:), pv(:), en(:), rwork(:)
    logical :: so
    complex(dp), allocatable :: hh(:, :), hm(:, :), om(:, :), cw(:), zc(:), cc(:)
    real(dp) :: bv(3, 3), k(3)
    integer :: ix, iy, iz, ik, ib, i, j, l, a, nst, nd, info, s, js
    external :: zhegv, zheev
    so = .false.; if (present(h2)) so = h2
    nd = 9 * c%nat
    bv = recip_vectors(c%avec)
    nst = nk**3 * nd
    allocate (eb(nst, nsp), xm(0:2, 0:2, nst * c%nat, nsp))
    xm = 0.0_dp
    do s = 1, nsp
      !$omp parallel private(hh, hm, om, ov, pv, en, rwork, cw, zc, cc, ix, iy, iz, ik, k, i, j, l, a, ib, info, js)
      allocate (hh(nd, nd), hm(nd, nd), om(nd, nd), ov(nd), pv(nd), en(nd), rwork(3 * nd), cw(4 * nd), zc(nd), cc(nd))
      !$omp do collapse(2)
      do iz = 0, nk - 1
        do iy = 0, nk - 1
          do ix = 0, nk - 1
            ik = ix + nk * (iy + nk * iz)
            k = (ix + 0.5_dp) / nk * bv(:, 1) + (iy + 0.5_dp) / nk * bv(:, 2) + (iz + 0.5_dp) / nk * bv(:, 3)
            call hsmall_k(c, cls, pps(:, s), alpha_sp, k, hh, ov, pv, en)
            if (so) then                          ! H = E_nu + h - h o h, orthonormal basis
              do j = 1, nd
                do i = 1, nd
                  hm(i, j) = hh(i, j) - sum(hh(i, :) * ov * hh(:, j))
                end do
                hm(j, j) = hm(j, j) + en(j)
              end do
              hm = 0.5_dp * (hm + conjg(transpose(hm)))
              call zheev('V', 'U', nd, hm, nd, eb(nd * ik + 1:nd * ik + nd, s), cw, 4 * nd, rwork, info)
              if (info /= 0) error stop 'bands_core: zheev failed'
              do ib = 1, nd
                do i = 1, nd
                  a = (i - 1) / 9 + 1
                  l = orb_l_of(mod(i - 1, 9) + 1)
                  js = (nd * ik + ib - 1) * c%nat + a
                  xm(0, l, js, s) = xm(0, l, js, s) + abs(hm(i, ib))**2
                  xm(1, l, js, s) = xm(1, l, js, s) + (eb(nd * ik + ib, s) - en(i)) * abs(hm(i, ib))**2
                  xm(2, l, js, s) = xm(2, l, js, s) + (eb(nd * ik + ib, s) - en(i))**2 * abs(hm(i, ib))**2
                end do
              end do
              cycle
            end if
            do j = 1, nd
              do i = 1, nd
                om(i, j) = (ov(i) + ov(j)) * hh(i, j) + sum(hh(i, :) * pv * hh(:, j))
                hm(i, j) = (0.5_dp + en(i) * ov(i) + 0.5_dp + en(j) * ov(j)) * hh(i, j) &
                           + sum(hh(i, :) * (ov + pv * en) * hh(:, j))
              end do
              om(j, j) = om(j, j) + 1.0_dp
              hm(j, j) = hm(j, j) + en(j)
            end do
            call zhegv(1, 'V', 'U', nd, hm, nd, om, nd, eb(nd * ik + 1:nd * ik + nd, s), cw, 4 * nd, rwork, info)
            if (info /= 0) error stop 'bands_core: zhegv failed'
            do ib = 1, nd
              zc = hm(:, ib)
              cc = matmul(hh, zc)
              do i = 1, nd
                a = (i - 1) / 9 + 1
                l = orb_l_of(mod(i - 1, 9) + 1)
                js = (nd * ik + ib - 1) * c%nat + a
                xm(0, l, js, s) = xm(0, l, js, s) + abs(zc(i))**2
                xm(1, l, js, s) = xm(1, l, js, s) + real(conjg(zc(i)) * cc(i), dp)
                xm(2, l, js, s) = xm(2, l, js, s) + abs(cc(i))**2
              end do
            end do
          end do
        end do
      end do
      !$omp end do
      deallocate (hh, hm, om, ov, pv, en, rwork, cw, zc, cc)
      !$omp end parallel
    end do
  end subroutine bands_core

  !> Fermi level and per-atom gamma-representation moments qat(k, l, atom, spin).
  subroutine bands_moments(c, cls, pps, alpha_sp, nk, nsp, nelec, sigma, bz, ef, qat, h2)
    type(cell_t), intent(in) :: c
    type(cluster_t), intent(in) :: cls(:)
    type(potpar_t), intent(in) :: pps(:, :)
    real(dp), intent(in) :: alpha_sp(0:, :), nelec, sigma
    character(len=*), intent(in) :: bz
    integer, intent(in) :: nk, nsp
    real(dp), intent(out) :: ef, qat(0:2, 0:2, c%nat, nsp)
    logical, intent(in), optional :: h2           !! orthogonal second-order H: moments taken directly (no o transformation)
    real(dp), allocatable :: eb(:, :), xm(:, :, :, :), eall(:, :), wall(:, :)
    real(dp) :: wk, f, lo, hi, nn, ol(0:2), fs
    integer :: i, l, a, nst, nd, it, s, ic, js
    type(potpar_t) :: q
    nd = 9 * c%nat
    nst = nk**3 * nd
    call bands_core(c, cls, pps, alpha_sp, nk, nsp, eb, xm, h2)
    wk = 1.0_dp / real(nk**3, dp)
    fs = 2.0_dp / real(nsp, dp)                    ! spin degeneracy per state
    if (trim(bz) == 'tetra') then
      allocate (eall(nd * nsp, nk**3), wall(nd * nsp, nk**3))
      do s = 1, nsp
        eall((s - 1) * nd + 1:s * nd, :) = reshape(eb(:, s), [nd, nk**3])
      end do
      call tetra_fermi(nk, eall, nelec, fs, ef)
      call tetra_weights(nk, eall, ef, wall)
    else
      lo = minval(eb) - 1.0_dp; hi = maxval(eb) + 1.0_dp
      do it = 1, 200
        ef = 0.5_dp * (lo + hi)
        nn = 0.0_dp
        do s = 1, nsp
          do i = 1, nst
            nn = nn + fs * wk * 0.5_dp * erfc((eb(i, s) - ef) / sigma)
          end do
        end do
        if (nn > nelec) then; hi = ef; else; lo = ef; end if
        if (hi - lo < 1.0e-12_dp) exit
      end do
      ef = 0.5_dp * (lo + hi)
    end if
    qat = 0.0_dp
    do s = 1, nsp
      do i = 1, nst
        if (trim(bz) == 'tetra') then
          f = fs * wall((s - 1) * nd + 1 + modulo(i - 1, nd), 1 + (i - 1) / nd)
        else
          f = fs * wk * 0.5_dp * erfc((eb(i, s) - ef) / sigma)
        end if
        do a = 1, c%nat
          ic = c%isp(a)
          q = to_alpha(pps(ic, s), alpha_sp(:, ic)); ol = o_alpha(q, alpha_sp(:, ic))
          if (present(h2)) then; if (h2) ol = 0.0_dp; end if
          js = (i - 1) * c%nat + a
          do l = 0, 2
            qat(0, l, a, s) = qat(0, l, a, s) + f * (xm(0, l, js, s) + 2.0_dp * ol(l) * xm(1, l, js, s) + ol(l)**2 * xm(2, l, js, s))
            qat(1, l, a, s) = qat(1, l, a, s) + f * (xm(1, l, js, s) + ol(l) * xm(2, l, js, s))
            qat(2, l, a, s) = qat(2, l, a, s) + f * xm(2, l, js, s)
          end do
        end do
      end do
    end do
  end subroutine bands_moments

  !> (class-averaged ... per atom) l-resolved DOS (states/Ry/atom/spin) of a converged result on
  !> a Gaussian-broadened k-mesh.  State weight per (atom,l) is the normalised ASA weight
  !> (xm0 + 2 o xm1 + p^alpha xm2)/norm, so each state carries unit weight in total.
  subroutine asa_dos(inp, cell, res, nk, egrid, sigma, dos)
    type(asa_input_t), intent(in) :: inp
    type(cell_t), intent(in) :: cell
    type(asa_result_t), intent(in) :: res
    integer, intent(in) :: nk
    real(dp), intent(in) :: egrid(:), sigma
    real(dp), intent(out) :: dos(size(egrid), 0:2, cell%nat, inp%nspin)
    type(cluster_t), allocatable :: cls(:)
    type(potpar_t), allocatable :: pps(:, :)
    real(dp), allocatable :: alpha_sp(:, :), eb(:, :), xm(:, :, :, :), olc(:, :, :), pa(:, :, :)
    real(dp) :: avw, wk, tot, w(0:2), x, g, fs
    integer :: ncl, nsp, ic, s, l, ia, i, nd, nst, js, ie
    type(potpar_t) :: q
    ncl = size(inp%cls); nsp = inp%nspin; avw = cell%w * inp%alat
    allocate (alpha_sp(0:2, ncl), pps(ncl, nsp), olc(0:2, ncl, nsp), pa(0:2, ncl, nsp))
    do ic = 1, ncl
      alpha_sp(:, ic) = alpha_from_sigma(inp%sigma, inp%cls(ic)%sws / avw)
      do s = 1, nsp
        do l = 0, 2
          pps(ic, s)%enu(l) = res%cls(ic)%pp(l, s)%enu; pps(ic, s)%c(l) = res%cls(ic)%pp(l, s)%c
          pps(ic, s)%sqrdel(l) = res%cls(ic)%pp(l, s)%sqrdel; pps(ic, s)%p(l) = res%cls(ic)%pp(l, s)%p
          pps(ic, s)%gamma(l) = res%cls(ic)%pp(l, s)%gamma; pps(ic, s)%alpha_rep(l) = res%cls(ic)%pp(l, s)%gamma
        end do
        q = to_alpha(pps(ic, s), alpha_sp(:, ic))
        olc(:, ic, s) = o_alpha(q, alpha_sp(:, ic)); pa(:, ic, s) = q%p
      end do
    end do
    allocate (cls(cell%nat))
    do ia = 1, cell%nat; call build_cluster(cell, ia, 4.0_dp * cell%w, alpha_sp, cls(ia)); end do
    call bands_core(cell, cls, pps, alpha_sp, nk, nsp, eb, xm)
    nd = 9 * cell%nat; nst = nk**3 * nd; wk = 1.0_dp / real(nk**3, dp)
    fs = 1.0_dp                                    ! per spin channel
    dos = 0.0_dp
    do s = 1, nsp
      do i = 1, nst
        tot = 0.0_dp
        do ia = 1, cell%nat
          ic = cell%isp(ia); js = (i - 1) * cell%nat + ia
          do l = 0, 2
            tot = tot + xm(0, l, js, s) + 2.0_dp * olc(l, ic, s) * xm(1, l, js, s) + pa(l, ic, s) * xm(2, l, js, s)
          end do
        end do
        do ia = 1, cell%nat
          ic = cell%isp(ia); js = (i - 1) * cell%nat + ia
          do l = 0, 2
            w(l) = (xm(0, l, js, s) + 2.0_dp * olc(l, ic, s) * xm(1, l, js, s) + pa(l, ic, s) * xm(2, l, js, s)) / tot
          end do
          do ie = 1, size(egrid)
            x = (egrid(ie) - eb(i, s)) / sigma
            if (abs(x) > 6.0_dp) cycle
            g = wk * exp(-0.5_dp * x * x) / (sigma * sqrt(2.0_dp * pi))
            dos(ie, :, ia, s) = dos(ie, :, ia, s) + g * w
          end do
        end do
      end do
    end do
  end subroutine asa_dos
end module asr_asa
