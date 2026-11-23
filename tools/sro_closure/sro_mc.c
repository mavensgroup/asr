/* Canonical Monte Carlo (species exchange, Metropolis) on a lattice embedded in a periodic simple-cubic grid of n^3 points.
 * Kernel of sro_exact_test.py (notes/FINDINGS_LOG.md, "B3 closure").  Build: gcc -O2 -shared -fPIC -o sro_mc.so sro_mc.c
 *
 * spec[(ix*n + iy)*n + iz] = species (0 .. nsp-1) of a lattice site, -1 on grid points that are not lattice sites.
 * sites[nsite]            = grid indices of the lattice sites.
 * nb[3*nnb]               = first-shell neighbour offsets in grid units.
 * eps[nsp*nsp]            = symmetric pair energies of first-shell pairs in units of kT.
 * local /= 0: exchange with a random first-shell neighbour; local = 0: exchange with a random site anywhere.
 * Returns the number of accepted exchanges.  One sweep = nsite attempts.
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
static inline long nbr(int n, long i, const int32_t *d) {
  int iz = (int)(i % n), iy = (int)((i / n) % n), ix = (int)(i / ((long)n * n));
  ix = (ix + d[0] + n) % n; iy = (iy + d[1] + n) % n; iz = (iz + d[2] + n) % n;
  return ((long)ix * n + iy) * n + iz;
}

long mc(int n, int8_t *spec, long nsite, const int64_t *sites, int nnb, const int32_t *nb, int nsp, const double *eps,
        long nsweep, uint64_t seed, int local) {
  long acc = 0;
  st[0] = splitmix(&seed); st[1] = splitmix(&seed);
  for (long sw = 0; sw < nsweep; sw++) {
    for (long t = 0; t < nsite; t++) {
      long i = sites[(long)(u01() * nsite)], j;
      if (local) j = nbr(n, i, nb + 3 * (int)(u01() * nnb)); else j = sites[(long)(u01() * nsite)];
      int si = spec[i], sj = spec[j];
      if (si == sj) continue;
      double de = 0.0;
      for (int k = 0; k < nnb; k++) {
        long a = nbr(n, i, nb + 3 * k), b = nbr(n, j, nb + 3 * k);
        if (a != j) de += eps[sj * nsp + spec[a]] - eps[si * nsp + spec[a]];
        if (b != i) de += eps[si * nsp + spec[b]] - eps[sj * nsp + spec[b]];
      }
      if (de <= 0.0 || u01() < exp(-de)) { spec[i] = (int8_t)sj; spec[j] = (int8_t)si; acc++; }
    }
  }
  return acc;
}

/* Same with pair energies on several shells: nb[3*nnb] are the offsets of shells 1..nshell, shell[nnb] the shell (0-based) of
 * each offset, eps[(s*nsp + a)*nsp + b] the energy (kT) of a pair (a, b) at distance shell s.  Used by the inverse Monte Carlo
 * that builds configurations with a given first-shell alpha and as little order as possible on the further shells. */
long mc_shells(int n, int8_t *spec, long nsite, const int64_t *sites, int nnb, const int32_t *nb, const int32_t *shell,
               int nsp, const double *eps, long nsweep, uint64_t seed) {
  long acc = 0;
  st[0] = splitmix(&seed); st[1] = splitmix(&seed);
  for (long sw = 0; sw < nsweep; sw++) {
    for (long t = 0; t < nsite; t++) {
      long i = sites[(long)(u01() * nsite)], j = sites[(long)(u01() * nsite)];
      int si = spec[i], sj = spec[j];
      if (si == sj) continue;
      double de = 0.0;
      for (int k = 0; k < nnb; k++) {
        const double *e = eps + (long)shell[k] * nsp * nsp;
        long a = nbr(n, i, nb + 3 * k), b = nbr(n, j, nb + 3 * k);
        if (a != j) de += e[sj * nsp + spec[a]] - e[si * nsp + spec[a]];
        if (b != i) de += e[si * nsp + spec[b]] - e[sj * nsp + spec[b]];
      }
      if (de <= 0.0 || u01() < exp(-de)) { spec[i] = (int8_t)sj; spec[j] = (int8_t)si; acc++; }
    }
  }
  return acc;
}
