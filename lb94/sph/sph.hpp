// Core SPH state and estimators (Hernquist & Katz 1989). Neighbor finding is brute-force
// here; it is replaced by the Barnes-Hut tree (tree.hpp) for production. Units are arbitrary
// (the driver picks G, length, velocity scales).
#ifndef SPH_SPH_HPP
#define SPH_SPH_HPP
#include <vector>
#include <cmath>
#include "kernel.hpp"

namespace sph {

struct Particle {
  double x, y, z;          // position
  double vx, vy, vz;       // velocity
  double ax, ay, az;       // acceleration
  double m, h, rho;        // mass, smoothing length, density
  double u, dudt;          // specific internal energy and its rate (adiabatic EOS)
  double cs2_iso = 0;      // assigned isothermal sound speed^2 (locally-isothermal EOS, LB94)
  double P, cs;            // pressure, sound speed
  double divv;             // velocity divergence (for the scalar AV / dh/dt)
  bool heavy = false;      // the central massive ("core") particle
  bool gas = true;         // false: collisionless (gravity-only) particle -- e.g. the star
  int rung = 0;            // block-timestep rung: dt = dt0 / 2^rung
  long long ti_end = 0;    // integer time (units of dt_min) at which this particle's step ends
};

inline double dist2(const Particle& a, const Particle& b) {
  double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
  return dx * dx + dy * dy + dz * dz;
}

// Adapt each h_i so the number of neighbors within 2h_i approaches N_target (HK89 eq. 2.17),
// then evaluate the symmetrized density (eq. 2.16). Brute force O(N^2): test/reference only.
inline void density_bruteforce(std::vector<Particle>& p, int N_target, int n_iter = 4) {
  int n = p.size();
  for (int it = 0; it < n_iter; ++it) {
    // (1) update h from the current neighbor counts
    for (int i = 0; i < n; ++i) {
      int cnt = 0;
      double tol2 = 4.0 * p[i].h * p[i].h;   // (2h)^2
      for (int j = 0; j < n; ++j)
        if (dist2(p[i], p[j]) < tol2) ++cnt;
      double ratio = std::cbrt((double)N_target / std::max(cnt, 1));
      p[i].h *= 0.5 * (1.0 + ratio);
    }
  }
  // (2) symmetrized density
  for (int i = 0; i < n; ++i) {
    double rho = 0.0;
    for (int j = 0; j < n; ++j) {
      double r = std::sqrt(dist2(p[i], p[j]));
      if (r < 2.0 * p[i].h || r < 2.0 * p[j].h)
        rho += p[j].m * 0.5 * (W(r, p[i].h) + W(r, p[j].h));
    }
    p[i].rho = rho;
  }
}

}  // namespace sph
#endif
