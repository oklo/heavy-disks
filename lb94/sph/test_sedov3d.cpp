// 3D Sedov-Taylor blast: deposit energy E0 at the center of a uniform particle lattice and
// evolve adiabatically (gamma=5/3, no self-gravity). The spherical shock expands as
// R_s = xi0 (E0 t^2 / rho0)^(1/5) ~ t^0.4, with xi0 ~ 1.15 for gamma=5/3.
#include "integrate.hpp"
#include <cstdio>

using namespace sph;

int main() {
  const int ns = 24;                                     // 24^3 = 13824 particles
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
  // deposit E0 in the central particles (r < 2d), weighted by mass
  double rdep = 2.0 * d, Mdep = 0;
  for (auto& q : p) if (q.x * q.x + q.y * q.y + q.z * q.z < rdep * rdep) Mdep += q.m;
  for (auto& q : p) if (q.x * q.x + q.y * q.y + q.z * q.z < rdep * rdep) q.u += E0 / Mdep;

  SPHParams par; par.gamma = gamma; par.isothermal = false; par.selfgrav = false;
  Tree tree; tree.P = &p;
  compute_forces(p, tree, par);

  auto shock_radius = [&]() {                            // radius of the peak density (shell)
    double rmax = 0, rhomax = 0;
    for (auto& q : p) {
      double r = std::sqrt(q.x * q.x + q.y * q.y + q.z * q.z);
      if (q.rho > rhomax) { rhomax = q.rho; rmax = r; }
    }
    return rmax;
  };

  printf("3D Sedov (N=%zu, gamma=%.3f):\n   t       R_s(sim)  R_s(Sedov)  ratio\n", p.size(), gamma);
  double t = 0; int nstep = 0;
  double samples[] = {0.02, 0.03, 0.04, 0.05}; int ns_i = 0;
  std::vector<std::pair<double, double>> track;
  while (t < 0.055 && ns_i < 4 && nstep < 100000) {
    double dt = std::min(timestep(p, par), samples[ns_i] - t);
    if (dt <= 0) dt = 1e-6;
    leapfrog_step(p, tree, par, dt); t += dt; ++nstep;
    if (t >= samples[ns_i] - 1e-9) {
      double Rs = shock_radius(), Rsed = 1.15 * std::pow(E0 * t * t / rho0, 0.2);
      printf("  %.3f    %.4f    %.4f    %.3f\n", t, Rs, Rsed, Rs / Rsed);
      track.push_back({t, Rs}); ++ns_i;
    }
  }
  double sx = 0, sy = 0, sxx = 0, sxy = 0; int nn = track.size();
  for (auto& q : track) { double x = std::log(q.first), y = std::log(q.second);
    sx += x; sy += y; sxx += x * x; sxy += x * y; }
  double slope = (nn * sxy - sx * sy) / (nn * sxx - sx * sx);
  printf("d log R_s / d log t = %.3f   (Sedov: 0.400)   steps=%d\n", slope, nstep);
  return 0;
}
