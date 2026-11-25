!> Configuration-averaged Bloch spectral function A(k, E) of the energy-dependent ASR (route edc).  Derivation:
!> manual/src/theory/spectral.md.
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
!> Unfolding.  A(k, E) is given per crystal site and per spin on Cartesian wavevectors k (units 2 pi / alat):
!>     A(k, E) = sum_u w_u sum_L sum_R cos(2 pi k . R) ( -1/pi Im <G_{0_u L, R L}(E + i eta)> ),
!> R the displacement in units of alat, w_u the site fraction of origin class u.  Basis sites of a sublattice fall into one origin
!> class when a pure translation or an inversion of the structure maps one onto the other (cos(2 pi k . R) is even in R and
!> D_LL(-I)^2 = 1).  Sites related only by rotations, screws or glides (such as NiAs cations and L1_2 face sites) form separate
!> origin classes whose unfolded sums are averaged with w_u.  Every sublattice must still be a single orbit of the space group so
!> its sites share one single-site potential: spectral_structure_ok.
!>
!> Window.  The Fourier sum stops at rmax.  A sharp stop rings: for an ordered metal with eta 0.01 Ry, A dips to -30 % of its
!> peak.  The default window is the ball window: w(R) is the overlap volume of two
!> balls of diameter rmax at distance R, divided by the volume of one ball, w(R) = (1 - R/rmax)^2 (1 + R/(2 rmax)) for R < rmax.
!> Its Fourier transform is |FT of the ball|^2 >= 0.  Samples of such a function on a lattice keep a non-negative Fourier series
!> (Poisson summation: the series is the periodic sum of the continuum transform).  So the windowed A is the exact A convolved in
!> k with a non-negative kernel of width about 1 / rmax: A >= 0 wherever the real-space elements are exact, at every k.  w(0) = 1,
!> so the BZ average (the on-site DOS) and the energy integral at each k are unchanged.  As rmax grows, w -> 1 and A -> the exact A.
module asr_spectral
  use asr_kinds, only: dp, pi
  use asr_crystal, only: structure_t
  use asr_crystal_symmetry, only: symop_t, space_group
  use asr_sublattice, only: decoration_t, decoration_structure
  use asr_augcls, only: augn_basis
  implicit none
  private
  public :: spectral_input_t, spectral_kpath, spectral_translations_ok, spectral_structure_ok, spectral_origin_classes, &
            spectral_unfold, spectral_dressing, spectral_window

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
    character(len=8) :: window = 'ball'                !! 'ball' (positive, module header) or 'sharp' (w = 1 up to rmax)
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

  !> .true. if every sublattice of the structure s is a single space-group orbit under the space group of s (with labels).
  !> Each sublattice holds one set of component potentials in class_scf, so every site of the sublattice must be symmetry
  !> equivalent to its representative site.  Structures where sites of a sublattice are related by rotations, screws or glides
  !> (such as NiAs and L1_2) are accepted and unfolded over all translation-inversion origin classes.
  logical function spectral_structure_ok(s) result(ok)
    type(structure_t), intent(in) :: s
    type(symop_t), allocatable :: ops(:)
    integer :: kap, b0, b, g
    logical :: hit
    call space_group(s, ops)
    ok = .true.
    do kap = 1, s%nlabel
      b0 = findloc(s%label, kap, 1)
      if (b0 == 0) cycle
      do b = 1, s%nb
        if (s%label(b) /= kap) cycle
        hit = .false.
        do g = 1, size(ops)
          if (ops(g)%perm(b0) == b) then; hit = .true.; exit; end if
        end do
        if (.not. hit) then; ok = .false.; return; end if
      end do
    end do
  end function spectral_structure_ok

  !> Partition the basis sites 1 .. s%nb into equivalence classes under pure translations (rot = +I) and inversions (rot = -I)
  !> in ops (space_group(s, ops)).  Returns bclass(1:nb) in 1 .. norig, representative basis sites orig_b(1:norig) ordered by
  !> sublattice and basis index, orig_sign(1:nb) (+1 if a pure translation maps b to orig_b(bclass(b)), -1 if an inversion does),
  !> and site fractions orig_w(1:norig) summing to 1.
  subroutine spectral_origin_classes(s, ops, bclass, orig_b, orig_sign, orig_w)
    type(structure_t), intent(in) :: s
    type(symop_t), intent(in) :: ops(:)
    integer, allocatable, intent(out) :: bclass(:), orig_b(:), orig_sign(:)
    real(dp), allocatable, intent(out) :: orig_w(:)
    integer :: ob(s%nb), kap, b, b2, g, norig
    allocate (bclass(s%nb), orig_sign(s%nb))
    bclass = 0; orig_sign = 1; norig = 0
    do kap = 1, s%nlabel
      do b = 1, s%nb
        if (s%label(b) /= kap .or. bclass(b) /= 0) cycle
        norig = norig + 1
        ob(norig) = b; bclass(b) = norig; orig_sign(b) = 1
        do b2 = b + 1, s%nb
          if (s%label(b2) /= kap .or. bclass(b2) /= 0) cycle
          do g = 1, size(ops)
            if (ops(g)%perm(b2) == b .and. maxval(abs(ops(g)%rot - eye3())) < 1.0e-7_dp) then
              bclass(b2) = norig; orig_sign(b2) = 1; exit
            end if
          end do
        end do
        do b2 = b + 1, s%nb
          if (s%label(b2) /= kap .or. bclass(b2) /= 0) cycle
          do g = 1, size(ops)
            if (ops(g)%perm(b2) == b .and. maxval(abs(ops(g)%rot + eye3())) < 1.0e-7_dp) then
              bclass(b2) = norig; orig_sign(b2) = -1; exit
            end if
          end do
        end do
      end do
    end do
    orig_b = ob(1:norig)
    orig_w = [(real(count(bclass == g), dp) / real(s%nb, dp), g=1, norig)]
  contains
    pure function eye3() result(e)
      real(dp) :: e(3, 3)
      integer :: i
      e = 0.0_dp; do i = 1, 3; e(i, i) = 1.0_dp; end do
    end function eye3
  end subroutine spectral_origin_classes

  !> .true. if every class of the decoration is one lattice of the translations that keep the decoration: then the Bloch sum
  !> from one origin per class is the class average (asr_spectral header).
  logical function spectral_translations_ok(dec) result(ok)
    type(decoration_t), intent(in) :: dec
    type(structure_t) :: s
    character(len=:), allocatable :: msg
    integer :: i, j, nt, v(3), h(3), m
    logical :: keep
    if (allocated(dec%tab)) then
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
    else
      call decoration_structure(dec, s, msg)
      ok = spectral_structure_ok(s)
    end if
  end function spectral_translations_ok

  !> weight of a displacement of length r (units alat) in the Fourier sum: 'ball' (module header) or 'sharp'
  pure real(dp) function spectral_window(kind, r, rmax) result(w)
    character(len=*), intent(in) :: kind
    real(dp), intent(in) :: r, rmax
    real(dp) :: d
    d = r / rmax
    if (d > 1.0_dp + 1.0e-12_dp) then
      w = 0.0_dp
    else if (trim(kind) == 'sharp') then
      w = 1.0_dp
    else
      d = min(d, 1.0_dp)
      w = (1.0_dp - d)**2 * (1.0_dp + 0.5_dp * d)
    end if
  end function spectral_window

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
