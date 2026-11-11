!> Guards of the class SCF (docs/REVIEW_fable.md section 6): a result is valid only if (1) E_nu is not held at its band window
!> in the last iteration (clamp_enu, asr_potpar::enu_window_margin) and (2) no valence weight lies below the bottom of the
!> occupation contour (class_scf_result_t%wbelow).  Small fcc Cu runs through the production routine class_scf.
program test_scf_guards
  use asr_kinds, only: dp
  use asr_driver, only: run_input_t, read_run_input, class_species
  use asr_class_scf, only: class_scf, class_scf_result_t, wbelow_tol
  use asr_potpar, only: enu_window_margin
  use asr_testing
  implicit none
  type(run_input_t) :: inp
  type(class_scf_result_t) :: r0, r1, r2, r3
  character(len=:), allocatable :: msg
  integer :: u
  character(len=*), parameter :: f = 'build/test_scf_guards_tmp.toml'
  character(len=1), parameter :: nl = new_line('a')
  character(len=*), parameter :: cu = '[structure]'//nl//'lattice = "fcc"'//nl//'alat = 3.61'//nl//'centering = "F"'//nl// &
    '[[species]]'//nl//'symbol = "Cu"'//nl//'z = 29'//nl//'nval = 11'//nl//'konfig = [4, 4, 3, 4]'//nl//'occupations = [1.0, 10.0]'//nl// &
    '[[sublattice]]'//nl//'name = "fcc"'//nl//'positions = [[0.0, 0.0, 0.0]]'//nl//'species = ["Cu"]'//nl// &
    '[scf]'//nl//'spin = 1'//nl//'xc = "vbh"'//nl//'mixing = 0.3'//nl//'max_iterations = 5'//nl//'tolerance = 1.0e3'//nl// &
    '[asr]'//nl//'levels = 2'//nl//'recursion_steps = 8'//nl//'supercell = 4'//nl//'contour_points = 16'//nl

  print '(a)', 'test_scf_guards'
  open (newunit=u, file=f, status='replace'); write (u, '(a)') cu; close (u)
  call read_run_input(f, inp, msg)
  call check(len(msg) == 0, 'fcc Cu input accepted')
  if (len(msg) > 0) then; print '(a)', msg; call report(); stop; end if
  call class_species(inp)
  inp%ci%ndos = 0

  ! (0) production settings: the tolerance is huge, so the run stops at iteration 4 on dV alone
  call class_scf(inp%ci, r0)
  print '(a,i0,a,i0,a,i0,a,2es10.2)', '    default: iterations ', r0%iter, '  clamp last/run ', r0%nclamp, ' / ', r0%nclamp_run, &
    '  weight below (first, final) ', r0%wbelow1, r0%wbelow
  call check(r0%iter == 4 .and. r0%converged, 'default margin: stops at iteration 4, reported converged')
  call check(r0%nclamp == 0 .and. r0%nclamp_run == 0, 'default margin: the E_nu clamp never fires')
  call check(r0%below_ok .and. abs(r0%wbelow) < 1.0e-5_dp .and. abs(r0%wbelow1) < 1.0e-5_dp, &
             'no valence weight below the contour bottom (< 1e-5 e, first and final potential)')

  ! (1) the margin is inert while the clamp does not fire: a looser margin (silent a fortiori) gives the same bits
  enu_window_margin = 0.01_dp
  call class_scf(inp%ci, r1)
  call check(r1%nclamp_run == 0 .and. r1%ef == r0%ef .and. r1%etot == r0%etot .and. r1%comp(1)%qnet == r0%comp(1)%qnet, &
             'margin 0.01, clamp silent: E_F, total energy and charge bit-identical to margin 0.05')

  ! (2) a margin that pushes E_nu out of its fixed point: the same dV test must not report convergence
  enu_window_margin = 0.999_dp
  call class_scf(inp%ci, r2)
  print '(a,i0,a,i0,a,l1,a,es10.2)', '    margin 0.999: clamp last/run ', r2%nclamp, ' / ', r2%nclamp_run, '  converged ', r2%converged, &
    '  E_F shift ', r2%ef - r0%ef
  call check(r2%nclamp > 0, 'margin 0.999: the clamp is active in the last iteration')
  call check(.not. r2%converged, 'active clamp in the last iteration: not reported as converged')
  call check(abs(r2%ef - r0%ef) > 1.0e-4_dp, 'test sensitivity: the clamped run has a different E_F')
  enu_window_margin = 0.05_dp

  ! (3) a contour that starts inside the band: the missing weight is measured and flagged
  inp%ci%ebot_offset = 0.05_dp
  call class_scf(inp%ci, r3)
  print '(a,2es10.2)', '    contour bottom 0.05 Ry below the lowest E_nu: weight below (first, final) ', r3%wbelow1, r3%wbelow
  call check(abs(r3%wbelow) > wbelow_tol .and. .not. r3%below_ok, 'contour bottom inside the band: missing weight found and flagged')
  call report()
end program test_scf_guards
