!> Terminator unit tests + the gapped-band stress test.
program test_terminators
  use asr_kinds, only: dp, pi
  use asr_operator, only: csr_op
  use asr_lattice_toy, only: cubic_cluster
  use asr_recursion, only: lanczos
  use asr_greens, only: sqrt_tail, periodic_g
  use asr_kspace, only: kspace_green_b2
  use asr_terminator
  use asr_testing, only: check
  use asr_endep, only: cf_species, ed_tail
  implicit none
  integer, parameter :: ne = 81
  integer :: n, k, centre, m, ie
  real(dp) :: a(60), b(60), ainf, binf, egrid(ne), e_avg, e_per, d
  real(dp), allocatable :: v(:)
  complex(dp) :: z, g1, g2, gref(ne), g
  type(csr_op) :: op
  type(terminator_t) :: t
  integer :: ni
  logical :: ok
  complex(dp) :: ca(6), cb(6), cc(6), fdum
  ! test 7: two energy-dependent Cr e_g chains (Fe66Ni14Cr20 fcc levels 4, 24 steps, last main-line contour node, scan energies
  ! -0.070915 / -0.070885 Ry of the 2026-10-30 chain dump); columns Re a, Im a, Re b, Im b
  real(dp), parameter :: cr1(4, 24) = reshape([ &
    -3.33122997700464329e-02_dp, 0.00000000000000000e+00_dp, 2.70339555656921693e-01_dp, 0.00000000000000000e+00_dp, &
    6.75072018375181826e-01_dp, -4.07445143372574953e-03_dp, 7.16164049263389835e-01_dp, 1.14024258075473214e-02_dp, &
    4.06668138190933504e-01_dp, -2.47341884740947433e-02_dp, 4.15664263725236793e-01_dp, 1.15821666175592846e-02_dp, &
    6.88657654055072954e-01_dp, -3.45877843937521315e-02_dp, 6.89862167617009803e-01_dp, 1.47683471278086519e-02_dp, &
    6.20654515566229303e-01_dp, -1.46111060613956425e-02_dp, 5.05954297529322750e-01_dp, 1.61937217610468552e-02_dp, &
    5.49198220147881666e-01_dp, -3.93618477486828711e-02_dp, 6.48394243412637516e-01_dp, 3.69004790075446029e-03_dp, &
    7.60694418275561168e-01_dp, 2.19437840533611090e-03_dp, 5.53909598640540013e-01_dp, 8.96764562700077573e-03_dp, &
    4.88907554095473540e-01_dp, 2.18437347959869440e-03_dp, 6.14846054530356789e-01_dp, -1.43912983587100411e-02_dp, &
    7.38783739494219827e-01_dp, 3.65454756475329032e-02_dp, 5.95156227152380146e-01_dp, -2.63550008883733783e-02_dp, &
    5.21617174910467996e-01_dp, 9.03499756950708194e-02_dp, 6.38690254787860123e-01_dp, -4.71604193136712785e-02_dp, &
    4.83354078739727910e-01_dp, 7.11316414019380339e-02_dp, 7.09550819652795428e-01_dp, 3.78722126507435126e-03_dp, &
    4.72073557913776198e-01_dp, -7.95950280263672399e-02_dp, 6.42156902030795562e-01_dp, 5.55600153210076900e-02_dp, &
    5.56459546769577096e-01_dp, -1.00151301713672772e-01_dp, 6.01369894754317680e-01_dp, 3.97625264383723726e-02_dp, &
    6.33838885610439107e-01_dp, -4.94625333908548862e-02_dp, 5.99592980474140669e-01_dp, 1.17319382488059575e-02_dp, &
    6.79255310308420013e-01_dp, 1.09063796395689628e-02_dp, 5.71932919191747491e-01_dp, -2.09532206142636138e-02_dp, &
    6.22883758672819532e-01_dp, 1.01671394065303747e-01_dp, 6.14084362982432741e-01_dp, -8.13750322806767595e-02_dp, &
    4.94236829924825649e-01_dp, 2.25687250131880013e-01_dp, 7.33242165052063410e-01_dp, -9.09440944262238649e-02_dp, &
    1.51070009754013052e-01_dp, -9.38378672686341064e-02_dp, 6.93818970102222710e-01_dp, 1.43858697912973399e-01_dp, &
    6.85142331460038112e-01_dp, -2.26592831808112460e-01_dp, 6.23244772875802711e-01_dp, 7.73612105190774058e-02_dp, &
    6.78578633147831733e-01_dp, -7.66459811854565826e-02_dp, 5.40625401316164544e-01_dp, 2.36770437467367957e-02_dp, &
    6.95683189903854915e-01_dp, -1.33769545046196169e-02_dp, 5.53552862613546748e-01_dp, -1.05647155292150213e-02_dp, &
    7.87606361209532535e-01_dp, 7.50418057840613606e-02_dp, 4.45638634049426496e-01_dp, -1.29433302111628484e-01_dp, &
    8.75608211656740965e-01_dp, 1.01873008612916327e+00_dp, 1.28527845528295903e+00_dp, -6.48000841064881095e-01_dp, &
    -7.42764733190824256e-01_dp, -8.80382693048855414e-01_dp, 0.00000000000000000e+00_dp, 0.00000000000000000e+00_dp], [4, 24])
  real(dp), parameter :: cr2(4, 24) = reshape([ &
    -3.33122997700464329e-02_dp, 0.00000000000000000e+00_dp, 2.70339555656921693e-01_dp, 0.00000000000000000e+00_dp, &
    6.75070804159054072e-01_dp, -4.07445143372574953e-03_dp, 7.16167448770739901e-01_dp, 1.14025192465645412e-02_dp, &
    4.06660761208871901e-01_dp, -2.47348751569104082e-02_dp, 4.15667722698309317e-01_dp, 1.15831402674129846e-02_dp, &
    6.88647280344467205e-01_dp, -3.45922171239722090e-02_dp, 6.89866592876719809e-01_dp, 1.47704590860994224e-02_dp, &
    6.20650349097869980e-01_dp, -1.46139418786301409e-02_dp, 5.05959038837746400e-01_dp, 1.61966007668465153e-02_dp, &
    5.49185922294442408e-01_dp, -3.93696713002834042e-02_dp, 6.48395810798234273e-01_dp, 3.69197403259828236e-03_dp, &
    7.60694402764079403e-01_dp, 2.18825680790949882e-03_dp, 5.53912984210284010e-01_dp, 8.97511059856990340e-03_dp, &
    4.88906727850456191e-01_dp, 2.16936922941306535e-03_dp, 6.14840742143402896e-01_dp, -1.43871086799508101e-02_dp, &
    7.38801863288650384e-01_dp, 3.65289807024149621e-02_dp, 5.95139462676501751e-01_dp, -2.63416194742365009e-02_dp, &
    5.21679462916404768e-01_dp, 9.03424297052562580e-02_dp, 6.38654109818242777e-01_dp, -4.71820798901987984e-02_dp, &
    4.83393579792765815e-01_dp, 7.12337645301274885e-02_dp, 7.09575811981096649e-01_dp, 3.73084287807427676e-03_dp, &
    4.71964721217190852e-01_dp, -7.95514410671090549e-02_dp, 6.42205202700369293e-01_dp, 5.55736412966769122e-02_dp, &
    5.56405678242677726e-01_dp, -1.00208605844399848e-01_dp, 6.01381134885095037e-01_dp, 3.97911649313058904e-02_dp, &
    6.33831596696031840e-01_dp, -4.95128995788042189e-02_dp, 5.99592653447139856e-01_dp, 1.17566606027181986e-02_dp, &
    6.79264600372840310e-01_dp, 1.08470573611063933e-02_dp, 5.71920640239631184e-01_dp, -2.09102160694130956e-02_dp, &
    6.22967841745487250e-01_dp, 1.01557896068829132e-01_dp, 6.13980945416682355e-01_dp, -8.13340764549757705e-02_dp, &
    4.94584372082812473e-01_dp, 2.25823674170973820e-01_dp, 7.33154731064137422e-01_dp, -9.12034699727222592e-02_dp, &
    1.50662454808811214e-01_dp, -9.34477633958859799e-02_dp, 6.94038881420415366e-01_dp, 1.43945840332219382e-01_dp, &
    6.84999462408643178e-01_dp, -2.26845942341028783e-01_dp, 6.23249774374431409e-01_dp, 7.74450404443587370e-02_dp, &
    6.78613588348078056e-01_dp, -7.67531840139244886e-02_dp, 5.40622135568460616e-01_dp, 2.37183282150916185e-02_dp, &
    6.95659194664182912e-01_dp, -1.34765033011578182e-02_dp, 5.53576367387244472e-01_dp, -1.04841123686454738e-02_dp, &
    7.87543270537403028e-01_dp, 7.47057409973977127e-02_dp, 4.45516616044494940e-01_dp, -1.28956038759317854e-01_dp, &
    8.79522300924381439e-01_dp, 1.01723332039711778e+00_dp, 1.28313745737377172e+00_dp, -6.51408065828392058e-01_dp, &
    -7.47591743472729053e-01_dp, -8.78137955125679492e-01_dp, 0.00000000000000000e+00_dp, 0.00000000000000000e+00_dp], [4, 24])
  complex(dp), parameter :: zc1 = cmplx(-8.30731261151762546e-02_dp, 1.00530964914873386e-01_dp, dp), &
                            zc2 = cmplx(-8.30431671572485763e-02_dp, 1.00530964914873386e-01_dp, dp)
  complex(dp) :: a1(24), b1(24), a2(24), b2(24), c0(24)
  real(dp) :: rh, rb

  print '(a)', 'test_terminators'
  ! 1. p=1 periodic tail == sqrt terminator
  z = cmplx(1.3_dp, 0.2_dp, dp)
  g1 = periodic_g(z, [0.4_dp], [1.7_dp]); g2 = sqrt_tail(z, 0.4_dp, 1.7_dp) / 1.7_dp**2
  call check_c(abs(g1 - g2), 1.0e-13_dp, 'periodic_g(p=1) == sqrt terminator')
  ! 2. period-2 tail is its own fixed point: g = 1/(z-a1-b1^2/(z-a2-b2^2 g))
  g1 = periodic_g(z, [0.5_dp, -0.5_dp], [1.2_dp, 0.9_dp])
  g2 = 1.0_dp / (z - 0.5_dp - 1.2_dp**2 / (z + 0.5_dp - 0.9_dp**2 * g1))
  call check_c(abs(g1 - g2), 1.0e-13_dp, 'period-2 tail is a fixed point')
  call check_c(max(0.0_dp, aimag(g1)), 1.0e-14_dp, 'period-2 tail on physical sheet (Im g<0)')
  ! 3. Beer-Pettifor on a uniform chain recovers (0,1) as depth grows
  a = 0.0_dp; b = 1.0_dp
  do n = 8, 56, 16
    call bp_asymptote(a, b, n, ainf, binf)
    print '(a,i3,a,es10.3,a,es10.3)', '    BP uniform chain n=', n, '  |A|=', abs(ainf), '  |B-1|=', abs(binf - 1.0_dp)
  end do
  call check(abs(binf - 1.0_dp) < 1.0e-2_dp .and. abs(ainf) < 1.0e-12_dp, 'BP -> (0,1) on uniform chain, n=56')
  ! 4. Luchini-Nex joins linearly between last computed and asymptote
  a = 0.0_dp; b = 1.0_dp; a(8) = 0.4_dp; b(8) = 1.4_dp
  ni = 4
  t = make_terminator('sqrt-avg2', 'luchini-nex', a, b, 8, ni)
  call check(t%n == 12, 'LN extends the chain by ninterp levels')
  call check_c(abs(t%a(12) - (a(8) + 4.0_dp / 5.0_dp * (t%ainf - a(8)))), 1.0e-14_dp, 'LN a interpolation')
  call check_c(abs(t%b(10) - (b(8) + 2.0_dp / 5.0_dp * (t%binf - b(8)))), 1.0e-14_dp, 'LN b interpolation')
  ! 5. gapped B2 compound: period-2 tail must beat every one-band tail
  do ie = 1, ne
    egrid(ie) = -9.0_dp + 18.0_dp * (ie - 1) / real(ne - 1, dp)
  end do
  d = 3.0_dp
  gref = kspace_green_b2(1.0_dp, d, egrid, 0.1_dp, 96)
  call cubic_cluster('bcc', 2 * 32 + 6, 1.0_dp, op, centre, b2_delta=d)
  allocate (v(op%n)); v = 0.0_dp; v(centre) = 1.0_dp
  call lanczos(op, v, 32, a, b, m)
  e_avg = 0.0_dp; e_per = 0.0_dp
  do k = 1, 2
    if (k == 1) then; t = make_terminator('sqrt-avg2', 'butt', a, b, 32, 0)
    else; t = make_terminator('periodic2', 'butt', a, b, 32, 0); end if
    do ie = 1, ne
      z = cmplx(egrid(ie), 0.1_dp, dp)
      g = terminated_green(t, z)
      if (k == 1) then; e_avg = max(e_avg, abs(aimag(g - gref(ie))) / pi)
      else; e_per = max(e_per, abs(aimag(g - gref(ie))) / pi); end if
    end do
  end do
  print '(a,es10.3,a,es10.3)', '    B2 gapped, depth 32: sqrt-avg2 err=', e_avg, '   periodic2 err=', e_per
  ok = e_per < 1.0e-3_dp .and. e_per * 20.0_dp < e_avg
  call check(ok, 'B2: periodic2 < 1e-3 and >=20x better than one-band tail')
  ! 6. complex-symmetric chain tail (asr_endep%cf_species, ed_tail = 'bp') on a real chain == real Beer-Pettifor (bp_asymptote)
  a(1:6) = [0.10_dp, -0.05_dp, 0.20_dp, 0.00_dp, 0.15_dp, 0.05_dp]; b(1:6) = [0.90_dp, 1.10_dp, 0.95_dp, 1.05_dp, 1.00_dp, 1.00_dp]
  call bp_asymptote(a, b, 6, ainf, binf)
  ca = a(1:6); cb = b(1:6); cc = (0.0_dp, 0.0_dp)
  ed_tail = 'bp'
  d = 0.0_dp
  do ie = 1, 9
    z = cmplx(-2.4_dp + 0.6_dp * (ie - 1), 0.05_dp, dp)
    call cf_species(ca, cb, cc, 6, z, g1, fdum)
    g2 = 1.0_dp / (z - a(6) - sqrt_tail(z, ainf, binf))
    do k = 5, 1, -1
      g2 = 1.0_dp / (z - a(k) - b(k)**2 * g2)
    end do
    d = max(d, abs(g1 - g2))
  end do
  ed_tail = 'hull'
  call check_c(d, 1.0e-10_dp, 'cf_species bp tail == bp_asymptote on a real chain')
  ! 7. Nevanlinna bound.  G_c(z) = <1|(z - H^(2))^-1|1> exactly (E - H(E) = T (E - H^(2)) T, T = 1 on the centre, H^(2) real symmetric),
  !    so |dG/dz| <= -Im G / Im z for Im z > 0, whatever the system.  A terminated chain whose G breaks this has a spurious pole near z.
  !    hull obeys it on the Cr e_g chains above (ratio ~0.4); bp exceeds it ~1000x (its fraction has a pole 0.01 Ry from z), which is
  !    what stalled the levels-4 SCF.
  a1 = cmplx(cr1(1, :), cr1(2, :), dp); b1 = cmplx(cr1(3, :), cr1(4, :), dp)
  a2 = cmplx(cr2(1, :), cr2(2, :), dp); b2 = cmplx(cr2(3, :), cr2(4, :), dp); c0 = (0.0_dp, 0.0_dp)
  rh = nevanlinna_ratio('hull'); rb = nevanlinna_ratio('bp')
  print '(a,f8.3,a,f10.1)', '    Cr e_g chain, |dG/dz| / (-Im G / Im z): hull', rh, '   bp', rb
  call check(rh < 1.0_dp, 'hull tail obeys the Nevanlinna bound |dG/dz| <= -Im G / Im z (energy-dependent Cr e_g chain)')
  call check(rb > 1.0_dp, 'test sensitivity: the bp tail breaks the bound on the same chain')
  call report_all()
contains
  !> |G(z2) - G(z1)| / |z2 - z1| divided by min(-Im G) / Im z for the two chains of test 7, terminated with the given tail
  real(dp) function nevanlinna_ratio(tail)
    character(*), intent(in) :: tail
    complex(dp) :: gz1, gz2, fz
    call cf_species(a1, b1, c0, 24, zc1, gz1, fz, tail_in=tail)
    call cf_species(a2, b2, c0, 24, zc2, gz2, fz, tail_in=tail)
    nevanlinna_ratio = abs(gz2 - gz1) / abs(zc2 - zc1) / (min(-aimag(gz1), -aimag(gz2)) / aimag(zc1))
  end function nevanlinna_ratio
  subroutine check_c(val, tol, name)
    real(dp), intent(in) :: val, tol
    character(*), intent(in) :: name
    call check(val <= tol, name)
  end subroutine check_c
  subroutine report_all()
    use asr_testing, only: report
    call report()
  end subroutine report_all
end program test_terminators
