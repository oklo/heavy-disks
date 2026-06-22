// Build the 3D SPH initial condition from the 2D collapse disk (LB94 pp.340-342):
//   - read the finest-grid (R,Z) state rho, v_phi, T from ybl_grid_*.dat (cgs);
//   - the central zone + first annulus + sink core M_c become ONE massive SPH particle;
//   - distribute N equal-mass particles over the disk cells per the density, each at a
//     jittered (R,Z), random azimuth, mirrored to +/-Z, with the cell's rotation and a
//     locally-isothermal cs^2 = k T / (mu m_H);
//   - write the particle IC in LB94 code units (length 100 AU, velocity 1 km/s, G=1).
#include "sph.hpp"
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

using namespace sph;

int main(int argc, char** argv) {
  const char* gridfile = (argc > 1) ? argv[1] : "lb94/ybl_grid_a0.000_40kyr.dat";
  const int N = (argc > 2) ? std::atoi(argv[2]) : 25000;
  const int i_core = 2;                                  // cells i<i_core -> central particle

  // cgs constants and LB94 code units
  const double AU = 1.496e13, G_cgs = 6.674e-8, Msun = 1.989e33;
  const double kB = 1.380649e-16, mH = 1.6726e-24, mu = 2.34;
  const double Lu = 100 * AU, Vu = 1e5, Mu = Vu * Vu * Lu / G_cgs;   // G=1 in code units
  const double Tu = Lu / Vu, yr = 3.156e7;

  FILE* fp = fopen(gridfile, "r");
  if (!fp) { fprintf(stderr, "cannot open %s\n", gridfile); return 1; }
  int NR, NZ; double dR, dZ, Mc, Jc, t_yr, alpha;
  char line[512];
  fgets(line, sizeof(line), fp);
  sscanf(line, "# NR=%d NZ=%d dR_cm=%lf dZ_cm=%lf Mc_g=%lf Jc=%lf t_yr=%lf alpha=%lf",
         &NR, &NZ, &dR, &dZ, &Mc, &Jc, &t_yr, &alpha);
  fgets(line, sizeof(line), fp);                         // column header
  std::vector<double> R(NR * NZ), Z(NR * NZ), rho(NR * NZ), vphi(NR * NZ), T(NR * NZ);
  for (int k = 0; k < NR * NZ; ++k)
    fscanf(fp, "%lf %lf %lf %lf %lf", &R[k], &Z[k], &rho[k], &vphi[k], &T[k]);
  fclose(fp);
  auto idx = [&](int i, int j) { return i * NZ + j; };
  auto vol = [&](int i) { return 2.0 * (2.0 * PI * R[idx(i, 0)] * dR * dZ); };  // full ring (+/-Z)

  // (1) central massive particle: sink core + the inner i<i_core gas
  double core_m = Mc, core_mT = 0;
  for (int i = 0; i < i_core; ++i)
    for (int j = 0; j < NZ; ++j) { double mm = rho[idx(i, j)] * vol(i); core_m += mm; core_mT += mm * T[idx(i, j)]; }
  double core_T = core_mT > 0 ? core_mT / (core_m - Mc + 1e-300) : 20.0;

  // (2) disk mass outside the core, and the equal particle mass
  double Mdisk = 0;
  for (int i = i_core; i < NR; ++i)
    for (int j = 0; j < NZ; ++j) Mdisk += rho[idx(i, j)] * vol(i);
  double mpart = Mdisk / N;

  std::mt19937 rng(2024);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  std::vector<Particle> P;
  Particle core{}; core.heavy = true; core.m = core_m / Mu; core.h = 2.0 * dR / Lu;
  core.cs2_iso = (kB * core_T / (mu * mH)) / (Vu * Vu);
  P.push_back(core);

  // (3) systematic sampling: accumulate fractional particle counts cell by cell
  double acc = 0;
  for (int i = i_core; i < NR; ++i)
    for (int j = 0; j < NZ; ++j) {
      acc += rho[idx(i, j)] * vol(i) / mpart;
      double cs2 = (kB * std::max(T[idx(i, j)], 3.0) / (mu * mH)) / (Vu * Vu);
      while (acc >= 1.0) {
        acc -= 1.0;
        double Rp = R[idx(i, j)] + (U(rng) - 0.5) * dR;
        double Zp = (Z[idx(i, j)] + (U(rng) - 0.5) * dZ) * (U(rng) < 0.5 ? 1.0 : -1.0);
        double phi = 2.0 * PI * U(rng), cphi = std::cos(phi), sphi = std::sin(phi);
        Particle q{};
        q.x = Rp * cphi / Lu; q.y = Rp * sphi / Lu; q.z = Zp / Lu;
        q.vx = -vphi[idx(i, j)] * sphi / Vu; q.vy = vphi[idx(i, j)] * cphi / Vu; q.vz = 0;
        q.m = mpart / Mu; q.h = 2.0 * dR / Lu; q.cs2_iso = cs2;
        P.push_back(q);
      }
    }

  // (4) write the IC + diagnostics
  FILE* out = fopen("lb94/sph/disk_ic.dat", "w");
  fprintf(out, "# LB94 SPH IC (code units: L=100AU, V=1km/s, G=1, Mu=%.4e g=%.4f Msun, Tu=%.0f yr)\n",
          Mu, Mu / Msun, Tu / yr);
  fprintf(out, "# N=%zu  from %s (t=%.0f yr, alpha=%.3f)\n", P.size(), gridfile, t_yr, alpha);
  fprintf(out, "# x y z vx vy vz m cs2 heavy\n");
  for (auto& q : P)
    fprintf(out, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %d\n",
            q.x, q.y, q.z, q.vx, q.vy, q.vz, q.m, q.cs2_iso, q.heavy ? 1 : 0);
  fclose(out);

  printf("Disk SPH IC from %s (t=%.0f yr, alpha=%.3f):\n", gridfile, t_yr, alpha);
  printf("  code units: Mu=%.4f Msun, Lu=100 AU, Vu=1 km/s, Tu=%.0f yr\n", Mu / Msun, Tu / yr);
  printf("  central particle: M=%.4f Msun (= M_c %.3f + inner gas), h=%.1f AU, T=%.0f K\n",
         core_m / Msun, Mc / Msun, 2 * dR / AU, core_T);
  printf("  disk: M=%.4f Msun, %d particles, m_part=%.3e Msun (%.2e Mu)\n",
         Mdisk / Msun, (int)P.size() - 1, mpart / Msun, mpart / Mu);
  printf("  wrote lb94/sph/disk_ic.dat\n");
  return 0;
}
