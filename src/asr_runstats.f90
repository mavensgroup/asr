!> Run statistics that the SCF routes keep for the report: wall and CPU time per part of the run, and the names of the
!> total-energy terms (the index order of eterm / ecomp in asr_class_scf and asr_asa).
!>
!> Timers are exclusive buckets.  A part that calls another part books the inner part in the inner bucket: the DOS on the
!> real axis calls the potential parameters and the recursion chains, so its time lies in tb_params and tb_dos (the chains of
!> the DOS).  The time outside every bucket is the total minus their sum.  Timers are started and stopped outside OpenMP
!> parallel regions only.  The wall time is that of this MPI rank.  The CPU time is the process CPU time of this rank
!> (cpu_time: the sum over its threads), so CPU / wall is at most the number of threads.  Busy-waiting threads count as CPU.
!>
!> work(bucket) counts the work units of a bucket.  The recursion buckets count chain evaluations: one chain at one energy
!> (start orbital x component x spin x energy), summed over all MPI ranks.
module asr_runstats
  use asr_kinds, only: dp, i8
  implicit none
  private
  public :: run_timer_t, timer_start, timer_stop, wall_seconds, cpu_seconds
  public :: tb_setup, tb_params, tb_chains, tb_below, tb_dos, tb_bands, tb_density, tb_mix, tb_io, n_tb, tb_name
  public :: et_band, et_core, et_dc, et_hartree, et_xc, et_mad, et_sim, n_et, et_name

  !> timer buckets
  integer, parameter :: tb_setup = 1, &   !! graph, class operators, screening model, k-space clusters (once)
                        tb_params = 2, &  !! radial solutions and potential parameters at E_nu
                        tb_chains = 3, &  !! recursion chains on the occupation contour (E_F search included)
                        tb_below = 4, &   !! recursion chains of the check below the contour bottom
                        tb_dos = 5, &     !! recursion chains of the DOS on the real axis
                        tb_bands = 6, &   !! k-space: bands, Brillouin-zone integration and the Fermi level
                        tb_density = 7, & !! densities, core states, Hartree, xc, electrostatics, total energy
                        tb_mix = 8, &     !! Anderson mixing, step limits and E_nu windows
                        tb_io = 9, &      !! restart file
                        n_tb = 9
  character(len=*), parameter :: tb_name(n_tb) = [character(len=36) :: &
    'setup (graph, operators, screening)', 'potential parameters', 'recursion chains: SCF contour', &
    'recursion chains: bottom check', 'recursion chains: DOS', 'k-space bands and BZ integration', &
    'density, core, Hartree, xc, energy', 'mixing and E_nu windows', 'restart file']

  !> total-energy terms (Ry).  Kohn-Sham form: output density in the input potential.  Their sum is etot (to rounding).
  !> v~_in = v_in + 2Z/r: the nuclear terms -2Z int n / r of the kinetic energy and of the sphere electrostatics cancel, so
  !> neither the double counting nor the Hartree term holds them.
  integer, parameter :: et_band = 1, &    !! valence band energy sum_l (q1 + E_nu Q_l), Q_l the sphere charge of the channel
                        et_core = 2, &    !! sum of the core eigenvalues (relaxed core)
                        et_dc = 3, &      !! double counting -int n v~_in
                        et_hartree = 4, & !! sphere Hartree energy 1/2 int n V_H[n] of the electrons
                        et_xc = 5, &      !! exchange-correlation energy int n eps_xc[n]
                        et_mad = 6, &     !! Madelung energy 1/2 q V_M (sublattice-average charges; k-space: the cell charges)
                        et_sim = 7, &     !! screened energy of the charge deviations 1/2 q V_F (SIM; with SRO the drive form)
                        n_et = 7
  character(len=*), parameter :: et_name(n_et) = [character(len=40) :: &
    'valence band energy', 'core eigenvalue sum', 'double counting -int n (v_in + 2Z/r)', 'sphere Hartree 1/2 int n V_H', &
    'exchange-correlation int n eps_xc', 'Madelung 1/2 q V_M', 'screened deviations 1/2 q V_F']

  type :: run_timer_t
    real(dp) :: wall(n_tb) = 0.0_dp, cpu(n_tb) = 0.0_dp   !! seconds per bucket
    integer(i8) :: work(n_tb) = 0_i8                      !! work units per bucket (see the module header)
    integer :: ncall(n_tb) = 0                            !! completed start/stop pairs per bucket
    real(dp) :: wall_total = 0.0_dp, cpu_total = 0.0_dp   !! the whole SCF call
    real(dp) :: w0(n_tb) = 0.0_dp, c0(n_tb) = 0.0_dp      !! start stamps (internal)
  end type run_timer_t

contains

  !> Wall-clock seconds since an arbitrary origin (differences only).
  function wall_seconds() result(t)
    real(dp) :: t
    integer(i8) :: c, r
    call system_clock(c, r)
    t = real(c, dp) / real(max(r, 1_i8), dp)
  end function wall_seconds

  !> Process CPU seconds of this rank (all threads).
  function cpu_seconds() result(t)
    real(dp) :: t
    call cpu_time(t)
  end function cpu_seconds

  subroutine timer_start(tm, ib)
    type(run_timer_t), intent(inout) :: tm
    integer, intent(in) :: ib
    tm%w0(ib) = wall_seconds(); tm%c0(ib) = cpu_seconds()
  end subroutine timer_start

  !> Stop bucket ib and add the elapsed time; work (optional) adds work units.
  subroutine timer_stop(tm, ib, work)
    type(run_timer_t), intent(inout) :: tm
    integer, intent(in) :: ib
    integer(i8), intent(in), optional :: work
    tm%wall(ib) = tm%wall(ib) + max(0.0_dp, wall_seconds() - tm%w0(ib))
    tm%cpu(ib) = tm%cpu(ib) + max(0.0_dp, cpu_seconds() - tm%c0(ib))
    tm%ncall(ib) = tm%ncall(ib) + 1
    if (present(work)) tm%work(ib) = tm%work(ib) + work
  end subroutine timer_stop

end module asr_runstats
