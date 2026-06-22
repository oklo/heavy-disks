// Validate the block-timestep integrator on the 3D Sedov blast: with the post-blast central
// particles on deep rungs and the ambient on shallow rungs, it must still give the Sedov
// scaling R_s ~ t^0.4 (matching the global-timestep test's 0.398).
#include "block.hpp"
#include <cstdio>
#include <vector>

using namespace sph;

int main() {
  const int ns = 24;
  const double L = 1.0, rho0 = 1.0, E0 = 1.0, gamma = 5.0 / 3.0;
  const double d = L / ns, m = rho0 * L * L * L / (ns * ns * ns);
  std::vector<Particle> p;
  for (int i = 0; i < ns; ++i)
    for (int j = 0; j < ns; ++j)
      for (int k = 0; k < ns; ++k) {
        Particle q{};
        q.x = (i + 0.5) * d - 0.5 * L; q.y = (j + 0.5) * d - 0.5 * L; q.z = (k + 0.5) * d - 0.5 * L;
        q.m = m; q.h = 1.5 * d; q.u = 1e-4;
        p.push_back(q);
      }
  double rdep = 2.0 * d, Mdep = 0;
  for (auto& q : p) if (q.x * q.x + q.y * q.y + q.z * q.z < rdep * rdep) Mdep += q.m;
  for (auto& q : p) if (q.x * q.x + q.y * q.y + q.z * q.z < rdep * rdep) q.u += E0 / Mdep;

  SPHParams par; par.gamma = gamma; par.isothermal = false; par.selfgrav = false;
  Tree tree; tree.P = &p;
  BlockStepper bs; bs.Rmax = 8; bs.dt0 = 1.0e-3;   // fixed dt0 -> modest rung spread (no limiter)
  bs.init(p, tree, par);
  printf("3D Sedov, BLOCK timesteps (N=%zu): dt0=%.2e, dt_min=%.2e (Rmax=%d)\n",
         p.size(), bs.dt0, bs.dt_min(), bs.Rmax);

  std::vector<std::pair<double, double>> track;
  double samples[] = {0.02, 0.03, 0.04, 0.05}; int si = 0;
  auto cb = [&](double t) {
    if (si < 4 && t >= samples[si] - 1e-9) {
      double rmax = 0, rhomax = 0;
      for (auto& q : p) {
        double r = std::sqrt(q.x * q.x + q.y * q.y + q.z * q.z);
        if (q.rho > rhomax) { rhomax = q.rho; rmax = r; }
      }
      double Rsed = 1.15 * std::pow(E0 * t * t / rho0, 0.2);
      printf("  t=%.3f  R_s=%.4f  Sedov=%.4f  ratio=%.3f\n", t, rmax, Rsed, rmax / Rsed);
      track.push_back({t, rmax}); ++si;
    }
  };
  bs.run(p, tree, par, 0.05, 0.0025, cb);

  int nn = track.size(); double sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (auto& q : track) { double x = std::log(q.first), y = std::log(q.second);
    sx += x; sy += y; sxx += x * x; sxy += x * y; }
  double slope = (nn * sxy - sx * sy) / (nn * sxx - sx * sx);
  printf("d log R_s/d log t = %.3f  (Sedov 0.400; global-step test gave 0.398)\n", slope);
  return 0;
}
