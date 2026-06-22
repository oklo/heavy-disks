// Time integration: kick-drift-kick leapfrog (HK89 eq. 2.35) with the Courant + acceleration
// timestep criteria (eq. 2.37). Global timestep here; individual/block timesteps are a later
// optimization. The internal energy is advanced together with the velocity (trapezoidal).
#ifndef SPH_INTEGRATE_HPP
#define SPH_INTEGRATE_HPP
#include <vector>
#include <cmath>
#include "forces.hpp"

namespace sph {

inline double timestep(const std::vector<Particle>& p, const SPHParams& par,
                       double C_cour = 0.3, double C_acc = 0.3) {
  double dt = 1e300;
  for (const auto& q : p) {
    double vsig = q.cs + 1.2 * (par.alpha * q.cs + par.beta * q.h * std::max(-q.divv, 0.0));
    if (vsig > 0) dt = std::min(dt, C_cour * q.h / vsig);
    double a = std::sqrt(q.ax * q.ax + q.ay * q.ay + q.az * q.az);
    if (a > 0) dt = std::min(dt, C_acc * std::sqrt(q.h / a));
  }
  return dt;
}

// one KDK leapfrog step; forces must be current on entry, and are current on exit
inline void leapfrog_step(std::vector<Particle>& p, Tree& tree, const SPHParams& par, double dt) {
  int n = p.size();
  for (int i = 0; i < n; ++i) {                          // kick (half) + drift
    Particle& q = p[i];
    q.vx += 0.5 * dt * q.ax; q.vy += 0.5 * dt * q.ay; q.vz += 0.5 * dt * q.az;
    if (!par.isothermal) q.u = std::max(q.u + 0.5 * dt * q.dudt, 1e-300);
    q.x += dt * q.vx; q.y += dt * q.vy; q.z += dt * q.vz;
  }
  compute_forces(p, tree, par);                          // a^{n+1}, dudt^{n+1} at v^{n+1/2}
  for (int i = 0; i < n; ++i) {                          // kick (half)
    Particle& q = p[i];
    q.vx += 0.5 * dt * q.ax; q.vy += 0.5 * dt * q.ay; q.vz += 0.5 * dt * q.az;
    if (!par.isothermal) q.u = std::max(q.u + 0.5 * dt * q.dudt, 1e-300);
  }
}

}  // namespace sph
#endif
