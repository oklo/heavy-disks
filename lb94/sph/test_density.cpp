// Validate the SPH kernel + adaptive-h density estimate. On a lattice (or relaxed glass) the
// estimate must recover the known mass density rho_true = M/L^3 for interior particles. A
// random (Poisson) distribution shows the usual SPH sampling bias -- which is exactly the
// Monte-Carlo noise LB94 use to seed the spiral modes.
#include "sph.hpp"
#include <cstdio>
#include <random>

using namespace sph;

int main() {
  const double L = 1.0, M = 1.0, rho_true = M / (L * L * L);
  std::mt19937 rng(12345);
  std::uniform_real_distribution<double> U(0.0, L);

  auto run = [&](bool grid, const char* label) {
    std::vector<Particle> p;
    if (grid) {
      int ns = 16;                                   // 16^3 = 4096 particles on a cubic lattice
      double d = L / ns;
      for (int i = 0; i < ns; ++i)
        for (int j = 0; j < ns; ++j)
          for (int k = 0; k < ns; ++k)
            p.push_back({(i + 0.5) * d, (j + 0.5) * d, (k + 0.5) * d});
    } else {
      p.resize(4000);
      for (auto& q : p) { q.x = U(rng); q.y = U(rng); q.z = U(rng); }
    }
    double h0 = L * std::cbrt(50.0 / (4.18879 * p.size()));
    for (auto& q : p) { q.m = M / p.size(); q.h = h0; }
    density_bruteforce(p, /*N_target=*/50);
    double sum = 0, sum2 = 0, hmean = 0; int ni = 0;
    for (auto& q : p)
      if (q.x > 0.25 && q.x < 0.75 && q.y > 0.25 && q.y < 0.75 && q.z > 0.25 && q.z < 0.75) {
        double e = q.rho / rho_true; sum += e; sum2 += e * e; hmean += q.h; ++ni;
      }
    double mean = sum / ni, rms = std::sqrt(sum2 / ni - mean * mean);
    printf("  %-8s N=%5zu  interior=%4d  <h>=%.4f  <rho/rho_true>=%.4f  rms=%4.1f%%  %s\n",
           label, p.size(), ni, hmean / ni, mean, 100 * rms,
           grid ? ((std::fabs(mean - 1.0) < 0.03) ? "PASS" : "CHECK") : "(Poisson bias, expected)");
  };
  printf("SPH density test (target 50 neighbors, rho_true=%.1f):\n", rho_true);
  run(true, "lattice");
  run(false, "random");
  return 0;
}
