/* Canonical Monte Carlo (species exchange inside one sublattice, Metropolis) for the pair model of asr_epi.
 * Kernel of sro_alpha.py.  Build: gcc -O2 -shared -fPIC -o epi_mc.so epi_mc.c -lm
 *
 * Sites 0 .. nsite-1.  spec[i] = global species index of site i (0 .. nsp-1).
 * Neighbour entries of site i: k = ptr[i] .. ptr[i+1]-1, partner nbr[k], energy table tab[k].
 * eps[(t*nsp + a)*nsp + b] = energy in units of kT of a bond whose table is t, with species a on the site that owns the entry
 * and b on the partner.  Every bond appears in the lists of both of its sites, with tables that are transposes of each other.
 * So the energy change of a site is the sum over its own entries, and a bond between the two exchanged sites cancels.
 * sub[i] = sublattice of site i; sl[sp[s] .. sp[s+1]-1] = sites of sublattice s.  An exchange picks a site and a second
 * site of the same sublattice.  Returns the number of accepted exchanges.  One sweep = nsite attempts.
 */
#include <math.h>
#include <stdint.h>

static uint64_t st[2];
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t rnd(void) {                       /* xoroshiro128+ */
  uint64_t s0 = st[0], s1 = st[1], r = s0 + s1;
  s1 ^= s0; st[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); st[1] = rotl(s1, 37);
  return r;
}
static inline double u01(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }
static uint64_t splitmix(uint64_t *x) {
  uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

static inline double local(long i, int a, const int8_t *spec, const int64_t *ptr, const int64_t *nbr, const int32_t *tab,
                           int nsp, const double *eps, long skip) {
  double e = 0.0;
  for (long k = ptr[i]; k < ptr[i + 1]; k++) {
    if (nbr[k] == skip) continue;
    e += eps[((long)tab[k] * nsp + a) * nsp + spec[nbr[k]]];
  }
  return e;
}

long epi_mc(long nsite, int8_t *spec, const int64_t *ptr, const int64_t *nbr, const int32_t *tab, int nsp, const double *eps,
            const int32_t *sub, const int64_t *sp, const int64_t *sl, long nsweep, uint64_t seed) {
  long acc = 0;
  st[0] = splitmix(&seed); st[1] = splitmix(&seed);
  for (long sw = 0; sw < nsweep; sw++) {
    for (long t = 0; t < nsite; t++) {
      long i = (long)(u01() * nsite);
      int s = sub[i];
      long j = sl[sp[s] + (long)(u01() * (sp[s + 1] - sp[s]))];
      int si = spec[i], sj = spec[j];
      if (si == sj) continue;
      /* the bond i-j itself keeps its energy (its tables are transposes): leave it out on both sides */
      double de = local(i, sj, spec, ptr, nbr, tab, nsp, eps, j) - local(i, si, spec, ptr, nbr, tab, nsp, eps, j)
                + local(j, si, spec, ptr, nbr, tab, nsp, eps, i) - local(j, sj, spec, ptr, nbr, tab, nsp, eps, i);
      if (de <= 0.0 || u01() < exp(-de)) { spec[i] = (int8_t)sj; spec[j] = (int8_t)si; acc++; }
    }
  }
  return acc;
}

/* total energy in units of kT: every bond once (half of the sum over the entries of both ends) */
double epi_energy(long nsite, const int8_t *spec, const int64_t *ptr, const int64_t *nbr, const int32_t *tab, int nsp,
                  const double *eps) {
  double e = 0.0;
  for (long i = 0; i < nsite; i++) e += local(i, spec[i], spec, ptr, nbr, tab, nsp, eps, -1);
  return 0.5 * e;
}
