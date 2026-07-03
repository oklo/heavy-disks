// Quiet-start IC of the pinned LB94 basic-state family (not the noisy grid snapshot):
//
//   Sigma(R) = 48.4 (R/100AU)^-2.51 exp[-(54.2AU/R) - (R/267AU)^2] g/cm^2
//   c_s(R)   = 0.534 (R/100AU)^-1/4 km/s,   M_* = 0.340 Msun (collisionless, gravity-only)
//
// Placement designed for azimuthal quietness:
//   * constant-n rings (n azimuthal sites, EQUALLY spaced, random phase per ring):
//     every m < n vanishes per ring to roundoff, and -- because n is the same for all
//     rings -- every particle in a ring sees an identical neighbor geometry, so the SPH
//     forces are azimuthally symmetric too;
//   * ring radii at equal-mass quantiles of the family's cumulative M(R);
//   * z sampled from a Gaussian with H = (1/H_tid + 1/H_sg)^-1, assigned in same-z
//     opposite-azimuth pairs so the m=1 bending seed cancels as well;
//   * v_phi initialized from the enclosed-mass rotation (the driver rebalances against
//     the code's own measured forces).
// Code units: L = 100 AU, V = 1 km/s, G = 1 (as the rest of lb94/sph).
#include "sph.hpp"
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

using namespace sph;

int main(int argc, char** argv) {
  const int N = (argc > 1) ? std::atoi(argv[1]) : 25000;
  const int nring = (argc > 2) ? std::atoi(argv[2]) : 125;   // azimuthal sites per ring (even)

  // cgs / code units
  const double AU = 1.496e13, G_cgs = 6.674e-8, Msun = 1.989e33;
  const double Lu = 100 * AU, Vu = 1e5, Mu = Vu * Vu * Lu / G_cgs;
  const double SigU = Mu / (Lu * Lu);                        // = 100.2 g/cm^2

  // the family (code units)
  const double Sig0 = 48.4 / SigU, pP = 2.51, Rin = 0.542, Rout = 2.67;
  const double cs0 = 0.534e5 / Vu, Mstar = 0.340 * Msun / Mu;
  auto Sigma = [&](double R) {
    return Sig0 * std::pow(R, -pP) * std::exp(-Rin / R - (R / Rout) * (R / Rout));
  };
  auto cs2_of = [&](double R) { return cs0 * cs0 * std::pow(R, -0.5); };

  // cumulative disk mass on a fine grid (0.02 -> 4.0 = 2 -> 400 AU)
  const int NG = 40000;
  std::vector<double> Rg(NG), Mc(NG);
  double acc = 0;
  for (int i = 0; i < NG; ++i) {
    Rg[i] = 0.02 + (4.0 - 0.02) * (i + 0.5) / NG;
    acc += 2.0 * PI * Rg[i] * Sigma(Rg[i]) * (4.0 - 0.02) / NG;
    Mc[i] = acc;
  }
  double Mdisk = acc;
  int nrings = N / nring;
  double mpart = Mdisk / (nrings * nring);

  std::mt19937 rng(777);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  std::normal_distribution<double> Gz(0.0, 1.0);

  std::vector<Particle> P;
  Particle star{}; star.heavy = true; star.gas = false;
  star.m = Mstar; star.h = 0.11;                            // 11 AU gravitational softening
  P.push_back(star);

  int ig = 0;
  for (int k = 0; k < nrings; ++k) {
    double Mtarget = (k + 0.5) * Mdisk / nrings;            // ring at the mass quantile
    while (ig < NG - 1 && Mc[ig] < Mtarget) ++ig;
    double Rk = Rg[ig];
    // scale height: tidal + self-gravity combined
    double Menc = Mstar + Mc[ig];
    double Om = std::sqrt(Menc / (Rk * Rk * Rk));
    double cs = std::sqrt(cs2_of(Rk));
    double Htid = cs / Om, Hsg = cs * cs / (PI * Sigma(Rk));
    double H = 1.0 / (1.0 / Htid + 1.0 / Hsg);
    double vphi = std::sqrt(std::max(Menc / Rk, 0.0));      // initial guess; driver rebalances
    double phase = U(rng);
    // same-z opposite-azimuth pairs: draw nring/2 z values, give each to j and j+nring/2
    std::vector<double> zs(nring);
    for (int j = 0; j < nring / 2; ++j) {
      double z = H * std::max(-3.0, std::min(3.0, Gz(rng)));
      zs[j] = z; zs[j + nring / 2] = z;
    }
    for (int j = 0; j < nring; ++j) {
      double phi = 2.0 * PI * (j + phase) / nring;
      Particle q{};
      q.x = Rk * std::cos(phi); q.y = Rk * std::sin(phi); q.z = zs[j];
      q.vx = -vphi * std::sin(phi); q.vy = vphi * std::cos(phi); q.vz = 0;
      q.m = mpart; q.h = 0.05; q.cs2_iso = cs2_of(Rk);
      P.push_back(q);
    }
  }

  FILE* out = fopen("lb94/sph/quiet_ic.dat", "w");
  fprintf(out, "# QUIET-START IC of the basic-state family (L=100AU, V=1km/s, G=1)\n");
  fprintf(out, "# N=%zu  nrings=%d x n=%d  M_star=%.4f  M_disk=%.4f (code)  seed=777\n",
          P.size(), nrings, nring, Mstar, Mdisk);
  fprintf(out, "# x y z vx vy vz m cs2 heavy\n");
  for (auto& q : P)
    fprintf(out, "%.17e %.17e %.17e %.17e %.17e %.17e %.10e %.10e %d\n",
            q.x, q.y, q.z, q.vx, q.vy, q.vz, q.m, q.cs2_iso, q.heavy ? 1 : 0);
  fclose(out);

  // placement quietness check: mass-weighted azimuthal moments of the gas
  for (int m = 1; m <= 4; ++m) {
    double cr = 0, ci = 0, mt = 0;
    for (auto& q : P) {
      if (q.heavy) continue;
      double phi = std::atan2(q.y, q.x);
      cr += q.m * std::cos(m * phi); ci += q.m * std::sin(m * phi); mt += q.m;
    }
    printf("  placement |c%d| = %.3e\n", m, std::sqrt(cr * cr + ci * ci) / mt);
  }
  printf("wrote lb94/sph/quiet_ic.dat  (%d rings x %d, m_part=%.3e, M_d/M_*=%.3f)\n",
         nrings, nring, mpart, Mdisk / Mstar);
  return 0;
}
