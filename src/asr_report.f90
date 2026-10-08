!> Human-readable run report (OUTCAR / pw.out style): banner, input echo, structure, SCF history, converged
!> per-class quantities, DOS summary, warnings, timing.  Plain text, fixed section headers ("=== name ===") so it can be
!> grepped; every number is in the units stated in the header of its section.
module asr_report
  use asr_kinds, only: dp, pi
  use asr_asa, only: asa_class_result_t
  implicit none
  private
  public :: rep_echo_input, rep_open, rep_close, rep_section, rep_line, rep_structure, rep_scf_history, rep_classes, rep_dos_at_ef, rep_wall

  character(len=*), parameter :: program_version = 'asr_new 0.1 (own TB-LMTO-ASA + augmented-space recursion)'
  integer, parameter :: dp_wall = selected_int_kind(18)

contains

  !> Wall-clock seconds since an arbitrary origin (differences only).
  function rep_wall() result(t)
    real(dp) :: t
    integer(dp_wall) :: c, r
    call system_clock(c, r)
    t = real(c, dp) / real(r, dp)
  end function rep_wall

  !> Open the report and write the banner.  nthreads from the caller (omp_get_max_threads).
  subroutine rep_open(u, file, title, infile, nthreads)
    integer, intent(out) :: u
    character(len=*), intent(in) :: file, title, infile
    integer, intent(in) :: nthreads
    character(len=8) :: d
    character(len=10) :: t
    open (newunit=u, file=file, status='replace', action='write')
    call date_and_time(date=d, time=t)
    write (u, '(a)') repeat('=', 78)
    write (u, '(a)') ' '//program_version
    write (u, '(a)') ' '//title
    write (u, '(a)') repeat('=', 78)
    write (u, '(a,a,a,a,a,a,a,a,a,a)') ' run started  ', d(1:4), '-', d(5:6), '-', d(7:8), '  ', t(1:2)//':'//t(3:4)//':'//t(5:6)
    write (u, '(a,a)') ' input file   ', trim(infile)
    write (u, '(a,i0)') ' OpenMP threads ', nthreads
    write (u, '(a)') ' units: lengths Bohr (Angstrom where stated), energies Rydberg, charges electrons, moments mu_B'
  end subroutine rep_open

  subroutine rep_close(u, converged, total_seconds)
    integer, intent(in) :: u
    logical, intent(in) :: converged
    real(dp), intent(in) :: total_seconds
    call rep_section(u, 'end of run')
    write (u, '(a,f10.2,a)') ' total wall time   ', total_seconds, ' s'
    if (converged) then
      write (u, '(a)') ' status            SCF CONVERGED'
    else
      write (u, '(a)') ' status            SCF NOT CONVERGED  (results below are not self-consistent)'
    end if
    close (u)
  end subroutine rep_close

  !> Copy the TOML input verbatim into the report.
  subroutine rep_echo_input(u, infile)
    integer, intent(in) :: u
    character(len=*), intent(in) :: infile
    character(len=512) :: line
    integer :: v, ios
    call rep_section(u, 'input (verbatim)')
    open (newunit=v, file=infile, status='old', action='read', iostat=ios)
    if (ios /= 0) then; write (u, '(a)') ' (input file not readable)'; return; end if
    do
      read (v, '(a)', iostat=ios) line
      if (ios /= 0) exit
      write (u, '(a)') '   | '//trim(line)
    end do
    close (v)
  end subroutine rep_echo_input

  subroutine rep_section(u, name)
    integer, intent(in) :: u
    character(len=*), intent(in) :: name
    write (u, '(a)') ''
    write (u, '(a)') '=== '//name//' '//repeat('=', max(4, 70 - len(name)))
  end subroutine rep_section

  subroutine rep_line(u, text)
    integer, intent(in) :: u
    character(len=*), intent(in) :: text
    write (u, '(a)') ' '//text
  end subroutine rep_line

  !> Structure: lattice, volume, sphere radii and filling.
  subroutine rep_structure(u, lattice, alat, avec, tau, isp, symbols, sws)
    integer, intent(in) :: u
    character(len=*), intent(in) :: lattice
    real(dp), intent(in) :: alat, avec(3, 3), tau(:, :), sws(:)
    integer, intent(in) :: isp(:)
    character(len=*), intent(in) :: symbols(:)
    real(dp), parameter :: ang = 0.529177210903_dp
    real(dp) :: det, vol, vsph
    integer :: ia
    det = avec(1, 1) * (avec(2, 2) * avec(3, 3) - avec(3, 2) * avec(2, 3)) &
        - avec(1, 2) * (avec(2, 1) * avec(3, 3) - avec(3, 1) * avec(2, 3)) &
        + avec(1, 3) * (avec(2, 1) * avec(3, 2) - avec(3, 1) * avec(2, 2))
    vol = abs(det) * alat**3
    call rep_section(u, 'structure')
    write (u, '(a,a)') ' lattice               ', trim(lattice)
    write (u, '(a,f12.6,a,f12.6,a)') ' lattice constant       ', alat, ' Bohr  = ', alat * ang, ' Angstrom'
    write (u, '(a)') ' primitive vectors (units of the lattice constant):'
    do ia = 1, 3
      write (u, '(a,i1,a,3f10.5)') '   a', ia, '  ', avec(:, ia)
    end do
    write (u, '(a,f14.5,a,f12.5,a)') ' cell volume            ', vol, ' Bohr^3 = ', vol * ang**3, ' Angstrom^3'
    write (u, '(a,i0)') ' atoms per cell         ', size(isp)
    write (u, '(a,f12.5,a)') ' volume per atom        ', vol * ang**3 / real(size(isp), dp), ' Angstrom^3'
    write (u, '(a)') ' atom  species   class        position (alat units)           S (Bohr)'
    vsph = 0.0_dp
    do ia = 1, size(isp)
      write (u, '(i5,2x,a8,i7,3x,3f10.5,f14.5)') ia, symbols(isp(ia)), isp(ia), tau(:, ia), sws(isp(ia))
      vsph = vsph + 4.0_dp * pi / 3.0_dp * sws(isp(ia))**3
    end do
    write (u, '(a,f9.5)') ' sphere volume / cell volume  ', vsph / vol
  end subroutine rep_structure

  !> SCF iteration table.
  subroutine rep_scf_history(u, symbols, niter, dv, ef, moment, qnet, tol, converged)
    integer, intent(in) :: u, niter
    character(len=*), intent(in) :: symbols(:)
    real(dp), intent(in) :: dv(:), ef(:), moment(:, :), qnet(:, :), tol
    logical, intent(in) :: converged
    integer :: it, ic, ncl
    ncl = size(symbols)
    call rep_section(u, 'SCF history')
    write (u, '(a,es9.2,a)') ' convergence criterion: max change of r*V/(1+r) between iterations < ', tol, ' Ry'
    write (u, '(a)', advance='no') '  iter      dV(Ry)     E_F(Ry)'
    do ic = 1, ncl; write (u, '(a12)', advance='no') 'm:'//trim(symbols(ic)); end do
    do ic = 1, ncl; write (u, '(a13)', advance='no') 'qnet:'//trim(symbols(ic)); end do
    write (u, '(a)') ''
    do it = 1, niter
      write (u, '(i6,es12.3,f12.6)', advance='no') it, dv(it), ef(it)
      do ic = 1, ncl; write (u, '(f12.5)', advance='no') moment(it, ic); end do
      do ic = 1, ncl; write (u, '(f13.5)', advance='no') qnet(it, ic); end do
      write (u, '(a)') ''
    end do
    if (converged) then
      write (u, '(a,i0,a)') ' converged after ', niter, ' iterations'
    else
      write (u, '(a,i0,a)') ' WARNING: not converged after ', niter, ' iterations'
    end if
  end subroutine rep_scf_history

  !> Converged per-class results: charges, moments, E_nu and potential parameters (gamma rep), s/p/d occupations.
  subroutine rep_classes(u, symbols, z, nval, sws, nat_class, cr, nspin, ef, moment_cell, label)
    integer, intent(in) :: u, nspin
    character(len=*), intent(in) :: symbols(:), label
    real(dp), intent(in) :: z(:), sws(:), ef, moment_cell
    integer, intent(in) :: nval(:), nat_class(:)
    type(asa_class_result_t), intent(in) :: cr(:)
    integer :: ic, s, l, ncl
    real(dp) :: qsum, qtot
    character(len=1), parameter :: lname(0:2) = ['s', 'p', 'd']
    ncl = size(symbols)
    call rep_section(u, 'converged results'//label)
    write (u, '(a,f14.8,a)') ' Fermi energy                ', ef, ' Ry'
    qsum = 0.0_dp
    do ic = 1, ncl
      qsum = qsum + real(nat_class(ic), dp) * cr(ic)%qnet
    end do
    write (u, '(a,es12.3,a)') ' net charge of the cell       ', qsum, ' e   (sum of class qnet x atoms; 0 = neutral)'
    if (nspin == 2) write (u, '(a,f12.5,a)') ' spin moment per cell         ', moment_cell, ' mu_B'
    do ic = 1, ncl
      write (u, '(a)') ''
      write (u, '(a,i0,a,a,a,f6.1,a,i0,a,f10.5,a,i0)') ' class ', ic, '  ', trim(symbols(ic)), '   Z=', z(ic), '  valence=', nval(ic), &
        '  S=', sws(ic), ' Bohr   atoms/cell=', nat_class(ic)
      write (u, '(a,f11.5,a)') '   qnet (electrons in sphere - Z)  ', cr(ic)%qnet, ' e'
      if (nspin == 2) write (u, '(a,f11.5,a)') '   spin moment in sphere           ', cr(ic)%moment, ' mu_B'
      write (u, '(a)') '   spin  l      E_nu         C        sqrt(Delta)       p         gamma        occ(q0)       q1         q2'
      do s = 1, nspin
        do l = 0, 2
          write (u, '(i7,2x,a1,f12.6,f12.6,f13.6,f12.6,f12.6,f12.5,f11.5,f11.5)') s, lname(l), cr(ic)%enu(l, s), cr(ic)%pp(l, s)%c, &
            cr(ic)%pp(l, s)%sqrdel, cr(ic)%pp(l, s)%p, cr(ic)%pp(l, s)%gamma, cr(ic)%q(0, l, s), cr(ic)%q(1, l, s), cr(ic)%q(2, l, s)
        end do
      end do
      qtot = 0.0_dp
      do s = 1, nspin; qtot = qtot + sum(cr(ic)%q(0, :, s)); end do
      write (u, '(a,f10.5,a)') '   valence electrons in sphere (sum q0) ', qtot, ' e'
      do s = 1, nspin
        write (u, '(a,i0,a,3f10.5)') '   occupation spin ', s, '  s p d:', cr(ic)%q(0, 0, s), cr(ic)%q(0, 1, s), cr(ic)%q(0, 2, s)
      end do
    end do
    write (u, '(a)') ''
    write (u, '(a)') ' (C, sqrt(Delta), p, gamma, q0..q2 are in the gamma representation; E_nu, C in Ry)'
  end subroutine rep_classes

  !> DOS summary: N(E_F) per species and l (states/Ry/atom/spin), interpolated linearly from the DOS grid.
  subroutine rep_dos_at_ef(u, symbols, egrid, dos, ef, smear, nspin, note)
    integer, intent(in) :: u, nspin
    character(len=*), intent(in) :: symbols(:), note
    real(dp), intent(in) :: egrid(:), dos(:, 0:, :, :), ef, smear
    integer :: ie, ia, s, l
    real(dp) :: w, n(0:2), ntot
    call rep_section(u, 'density of states at E_F')
    write (u, '(a,f9.4,a,f9.4,a,i0,a,f8.4,a)') ' window ', egrid(1), ' .. ', egrid(size(egrid)), ' Ry, ', size(egrid), ' points, broadening ', smear, ' Ry'
    write (u, '(a)') ' '//note
    if (ef < egrid(1) .or. ef > egrid(size(egrid))) then
      write (u, '(a)') ' WARNING: E_F is outside the DOS window'
      return
    end if
    ie = 1
    do while (ie < size(egrid) - 1 .and. egrid(ie + 1) < ef); ie = ie + 1; end do
    w = (ef - egrid(ie)) / (egrid(ie + 1) - egrid(ie))
    ntot = 0.0_dp
    write (u, '(a)') '  species  spin       s          p          d        total   (states/Ry/atom)'
    do ia = 1, size(symbols)
      do s = 1, nspin
        do l = 0, 2
          n(l) = (1.0_dp - w) * dos(ie, l, ia, s) + w * dos(ie + 1, l, ia, s)
        end do
        write (u, '(2x,a8,i5,4f11.4)') symbols(ia), s, n, sum(n)
        ntot = ntot + sum(n)
      end do
    end do
    write (u, '(a,f11.4,a)') ' sum over listed species and spins  ', ntot, '  (not concentration-weighted)'
  end subroutine rep_dos_at_ef

end module asr_report
