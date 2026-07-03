// SPH forces (Hernquist & Katz 1989): adaptive-h symmetrized density (2.16), EOS, the
// symmetrized pressure-gradient acceleration with Monaghan-Gingold artificial viscosity
// (2.20-2.23), and the thermal-energy rate (2.29). The per-particle kernels (density_one,
// hydro_one) let the block-timestep integrator recompute forces for an active SUBSET while
// the full-set wrappers (compute_forces) serve the global integrator and the tests.
#ifndef SPH_FORCES_HPP
#define SPH_FORCES_HPP
#include <vector>
#include <cmath>
#include "tree.hpp"

namespace sph {

struct SPHParams {
  double alpha = 1.0, beta = 2.0;     // Monaghan-Gingold artificial viscosity coefficients
  double eta2_fac = 0.01;             // eta^2 = eta2_fac * hbar^2 (softens mu_ij at r->0)
  double gamma = 5.0 / 3.0;           // adiabatic index (adiabatic EOS)
  bool isothermal = false;            // true: P = rho cs2_iso (locally isothermal, LB94)
  bool selfgrav = true;               // include tree gravity
  int N_neigh = 50;                   // target neighbor number
};

// ---- per-particle kernels (tree must be current) ----

// iterate h_i toward N_neigh via the gather count (HK89 eq. 2.17); a few passes suffice when
// h is already close (per step); returns the final neighbor count
inline int adapt_h_one(std::vector<Particle>& p, Tree& tree, const SPHParams& par, int i, int iters) {
  if (!p[i].gas) return 0;                       // collisionless: h is a fixed softening
  int cnt = par.N_neigh;
  for (int it = 0; it < iters; ++it) {
    cnt = tree.gather_count(i, 2.0 * p[i].h);
    double ratio = std::cbrt((double)par.N_neigh / std::max(cnt, 1));
    p[i].h *= 0.5 * (1.0 + std::min(std::max(ratio, 0.5), 1.5));
    if (std::abs(cnt - par.N_neigh) <= 1) break;
  }
  return cnt;
}

// symmetrized density of particle i (eq. 2.16), including the self term m_i W(0,h_i)
inline void density_one(std::vector<Particle>& p, Tree& tree, int i, std::vector<int>& nb) {
  if (!p[i].gas) { p[i].rho = p[i].m * W(0.0, p[i].h); return; }
  double rho = p[i].m * W(0.0, p[i].h);
  tree.sym_neighbors(i, nb);                     // full symmetric pair set (exact eq 2.16)
  for (int j : nb) {
    if (!p[j].gas) continue;                     // collisionless particles carry no SPH mass
    double r = std::sqrt(dist2(p[i], p[j]));
    rho += p[j].m * 0.5 * (W(r, p[i].h) + W(r, p[j].h));
  }
  p[i].rho = rho;
}

inline void eos_one(Particle& q, const SPHParams& par) {
  if (!q.gas) { q.P = 0.0; q.cs = 0.0; return; }
  if (par.isothermal) { q.P = q.rho * q.cs2_iso; q.cs = std::sqrt(q.cs2_iso); }
  else { q.P = (par.gamma - 1.0) * q.rho * q.u; q.cs = std::sqrt(par.gamma * std::max(q.P, 0.0) / q.rho); }
}

// pressure + artificial-viscosity acceleration, du/dt, divv, and gravity for particle i
inline void hydro_one(std::vector<Particle>& p, Tree& tree, const SPHParams& par, int i,
                      std::vector<int>& nb) {
  Particle& pi = p[i];
  if (!pi.gas) {                                 // collisionless: gravity only
    double gx = 0, gy = 0, gz = 0;
    if (par.selfgrav) tree.accel(i, gx, gy, gz);
    pi.ax = gx; pi.ay = gy; pi.az = gz; pi.dudt = 0; pi.divv = 0;
    return;
  }
  double ax = 0, ay = 0, az = 0, dudt = 0, divv = 0;
  double Pi_over = pi.P / (pi.rho * pi.rho);
  tree.sym_neighbors(i, nb);                     // symmetric pair set -> reciprocal forces
  for (int j : nb) {
    Particle& pj = p[j];
    if (!pj.gas) continue;                       // no SPH coupling to collisionless particles
    double dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
    double r2 = dx * dx + dy * dy + dz * dz, r = std::sqrt(r2);
    if (r == 0) continue;
    double dvx = pi.vx - pj.vx, dvy = pi.vy - pj.vy, dvz = pi.vz - pj.vz;
    double gw = 0.5 * (gradW(r, pi.h) + gradW(r, pj.h));
    double gwx = gw * dx, gwy = gw * dy, gwz = gw * dz;
    double Pij = 0.0, vdotr = dvx * dx + dvy * dy + dvz * dz;
    if (vdotr < 0.0) {
      double hbar = 0.5 * (pi.h + pj.h), cbar = 0.5 * (pi.cs + pj.cs), rhobar = 0.5 * (pi.rho + pj.rho);
      double mu = hbar * vdotr / (r2 + par.eta2_fac * hbar * hbar);
      Pij = (-par.alpha * cbar * mu + par.beta * mu * mu) / rhobar;
    }
    double fac = Pi_over + pj.P / (pj.rho * pj.rho) + Pij;
    ax -= pj.m * fac * gwx; ay -= pj.m * fac * gwy; az -= pj.m * fac * gwz;
    double vdotgw = dvx * gwx + dvy * gwy + dvz * gwz;
    dudt += pj.m * (Pi_over + 0.5 * Pij) * vdotgw;
    divv -= pj.m * vdotgw;
  }
  if (par.selfgrav) { double gx, gy, gz; tree.accel(i, gx, gy, gz); ax += gx; ay += gy; az += gz; }
  pi.ax = ax; pi.ay = ay; pi.az = az; pi.dudt = dudt; pi.divv = divv / pi.rho;
}

// ---- full-set wrappers (global integrator + tests) ----

inline void adapt_h(std::vector<Particle>& p, Tree& tree, const SPHParams& par, int iters = 30) {
  int n = p.size();
  for (int it = 0; it < iters; ++it) {
    tree.build();
    int nbad = 0;
#pragma omp parallel for reduction(+ : nbad) schedule(dynamic, 256)
    for (int i = 0; i < n; ++i) {
      int cnt = tree.gather_count(i, 2.0 * p[i].h);
      if (std::abs(cnt - par.N_neigh) > 1) ++nbad;
      double ratio = std::cbrt((double)par.N_neigh / std::max(cnt, 1));
      p[i].h *= 0.5 * (1.0 + std::min(std::max(ratio, 0.5), 1.5));
    }
    if (nbad == 0) break;
  }
  tree.build();
}

inline void compute_forces(std::vector<Particle>& p, Tree& tree, const SPHParams& par) {
  int n = p.size();
  adapt_h(p, tree, par, 6);
#pragma omp parallel for schedule(dynamic, 256)
  for (int i = 0; i < n; ++i) { std::vector<int> nb; density_one(p, tree, i, nb); }
#pragma omp parallel for schedule(dynamic, 256)
  for (int i = 0; i < n; ++i) eos_one(p[i], par);
#pragma omp parallel for schedule(dynamic, 256)
  for (int i = 0; i < n; ++i) { std::vector<int> nb; hydro_one(p, tree, par, i, nb); }
}

// ---- active-subset evaluation (block-timestep integrator); tree must be pre-built ----
inline void compute_forces_active(std::vector<Particle>& p, Tree& tree, const SPHParams& par,
                                  const std::vector<int>& active) {
  int na = active.size();
#pragma omp parallel for schedule(dynamic, 64)
  for (int a = 0; a < na; ++a) adapt_h_one(p, tree, par, active[a], 4);
#pragma omp parallel for schedule(dynamic, 64)
  for (int a = 0; a < na; ++a) { std::vector<int> nb; density_one(p, tree, active[a], nb); eos_one(p[active[a]], par); }
#pragma omp parallel for schedule(dynamic, 64)
  for (int a = 0; a < na; ++a) { std::vector<int> nb; hydro_one(p, tree, par, active[a], nb); }
}

}  // namespace sph
#endif
