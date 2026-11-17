!> Configuration-averaged Bloch spectral function A(k, E) of the energy-dependent ASR (route edc).  Derivation and review:
!> docs/DERIVATION_spectral_argon.md (Q1, Q2 Form 2B, Q3b, Q7 and the review section).
!>
!> Method.  The off-site elements come from the chains that the SCF runs anyway.  With the centre frozen to species c the scaled
!> operator gives G^(2) = T (E - H(E))^-1 T, and T = 1 on the centre.  So
!>     <G^(2)_{0L, RL'}> = sum_c x_c  <0L| (E - H^(c)(E))^-1 | T_R R L'>,
!> where T_R |R L'> = sum_k tau_k |R L', R in local state k> is the scaling of the target site in augmented space (dressing
!> weights tau_k = sum_s x_s T_s(L') phi_k(s), T_s = (D_c / D_s)^1/2).  The chain from |0 L> determines the vector everywhere,
!> so one chain per centre species, start orbital and energy gives every target site (asr_edc::edc_green_proj).  The spectral
!> chains use the operator reduced with D2h (the diagonal operations of the site group), where all nine real orbitals are one
!> dimensional, so every start orbital is allowed and no rotation between orbitals is needed.
!>
!> Symmetrisation.  The truncated augmented space does not give G_{0L,RL'} = G_{0L',-R L} exactly (two different chains, the
!> second from the origin of the partner sublattice).  The code averages the two and prints their largest difference.
!>
!> Unfolding.  A(k, E) is given on the Brillouin zone of the graph lattice (bcc, fcc, sc) per graph site and per spin:
!>     A(k, E) = sum_kappa f_kappa sum_L sum_R cos(2 pi k . R) ( -1/pi Im <G_{0 kappa L, R L}(E + i eta)> ),
!> R the displacement in units of alat, k in units of 2 pi / alat (Cartesian).  It is the trace over L of <k L| A |k L> with Bloch
!> vectors of the graph lattice, so it is non-negative for the exact average, and its BZ average is the DOS per site.  For an
!> ordered decoration it is the unfolded band structure.  It needs every sublattice to be one lattice of translations
!> (random alloys, B2, L2_1; not L1_2, whose face sites are related by rotations only): spectral_translations_ok.
!> The truncation of the Fourier sum at rmax is sharp (no window; review of the derivation).
module asr_spectral
  use asr_kinds, only: dp, pi
  use asr_sublattice, only: decoration_t
  use asr_augn, only: augn_basis
  implicit none
  private
  public :: spectral_input_t, spectral_kpath, spectral_translations_ok, spectral_unfold, spectral_dressing

  type :: spectral_input_t
    logical :: on = .false.
    real(dp), allocatable :: kpts(:, :)                !! corners of the k path (3, n), units 2 pi / alat, Cartesian
    character(len=16), allocatable :: labels(:)        !! their labels
    integer :: nk = 30                                 !! points per segment
    integer :: ne = 201                                !! energies
    integer :: nsteps = 0                              !! chain length (0: recursion_steps of the SCF)
    real(dp) :: emin = -0.9_dp, emax = 0.5_dp          !! energy window (Ry)
    real(dp) :: eta = 0.01_dp                          !! imaginary part of the energy (Ry)
    real(dp) :: rmax = 2.0_dp                          !! largest displacement in the Fourier sum (units of alat)
  end type spectral_input_t

contains

  !> k points along the path kpts(:, 1) -> kpts(:, 2) -> ...: nk points per segment, every corner once.  dist: the path length
  !> (units 2 pi / alat); icorner(j): the index of corner j in kv.
  subroutine spectral_kpath(kpts, nk, kv, dist, icorner)
    real(dp), intent(in) :: kpts(:, :)
    integer, intent(in) :: nk
    real(dp), allocatable, intent(out) :: kv(:, :), dist(:)
    integer, allocatable, intent(out) :: icorner(:)
    integer :: nc, ns, j, i, n
    nc = size(kpts, 2); ns = max(nc - 1, 0)
    n = ns * nk + 1
    allocate (kv(3, n), dist(n), icorner(nc))
    kv(:, 1) = kpts(:, 1); dist(1) = 0.0_dp; icorner(1) = 1; n = 1
    do j = 1, ns
      do i = 1, nk
        n = n + 1
        kv(:, n) = kpts(:, j) + (kpts(:, j + 1) - kpts(:, j)) * real(i, dp) / real(nk, dp)
        dist(n) = dist(n - 1) + norm2(kv(:, n) - kv(:, n - 1))
      end do
      icorner(j + 1) = n
    end do
  end subroutine spectral_kpath

  !> .true. if every class of the decoration is one lattice of the translations that keep the decoration: then the Bloch sum
  !> from one origin per class is the class average (asr_spectral header).  Translations: graph vectors v between two sites of
  !> the cell with tab(h + v) = tab(h) for every site h of the cell.  Each class must hold as many sites per cell as there are
  !> such translations.
  logical function spectral_translations_ok(dec) result(ok)
    type(decoration_t), intent(in) :: dec
    integer :: i, j, nt, v(3), h(3), m
    logical :: keep
    m = dec%m; nt = 0
    do i = 1, size(dec%pcls)
      v = dec%pts(:, i) - dec%pts(:, 1)
      keep = .true.
      do j = 1, size(dec%pcls)
        h = modulo(dec%pts(:, j) + v, m)
        if (dec%tab(h(1), h(2), h(3)) /= dec%pcls(j)) then; keep = .false.; exit; end if
      end do
      if (keep) nt = nt + 1
    end do
    ok = .true.
    do i = 1, dec%ncls
      if (count(dec%pcls == i) /= nt) ok = .false.
    end do
  end function spectral_translations_ok

  !> sum_t cos(2 pi k . d_t) v_t: the unfolded Bloch sum of real-space values v_t at the displacements d_t (units alat).
  pure real(dp) function spectral_unfold(disp, v, k) result(a)
    real(dp), intent(in) :: disp(:, :), v(:), k(3)
    integer :: t
    a = 0.0_dp
    do t = 1, size(v)
      a = a + cos(2.0_dp * pi * dot_product(k, disp(:, t))) * v(t)
    end do
  end function spectral_unfold

  !> Dressing weights tau(l, k + 1) = sum_s x_s T_s(l) phi_k(s), T_s(l) = (dc(l) / del(l, s))^1/2, for a target site whose class
  !> has the global species spc(1:n) with concentrations conc(1:n).  The local basis phi is that of asr_augcls::class_mt (built
  !> on the species of non-zero weight; states beyond them get zero).  frozen: the target is the frozen centre (tau = 1 for the
  !> reference state).
  subroutine spectral_dressing(dc, del, spc, conc, frozen, tau)
    real(dp), intent(in) :: dc(:), del(:, :), conc(:)
    integer, intent(in) :: spc(:)
    logical, intent(in) :: frozen
    real(dp), intent(out) :: tau(:, :)
    real(dp), allocatable :: phi(:, :)
    integer, allocatable :: sup(:)
    integer :: l, k, j, n
    tau = 0.0_dp
    if (frozen) then; tau(:, 1) = 1.0_dp; return; end if
    sup = pack([(j, j=1, size(conc))], conc > 0.0_dp)
    n = size(sup)
    call augn_basis(conc(sup), phi)
    do k = 1, n
      do l = 1, size(dc)
        do j = 1, n
          tau(l, k) = tau(l, k) + conc(sup(j)) * sqrt(dc(l) / del(l, spc(sup(j)))) * phi(k, j)
        end do
      end do
    end do
  end subroutine spectral_dressing
end module asr_spectral
