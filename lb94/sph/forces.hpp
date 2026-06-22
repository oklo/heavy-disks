// SPH forces (Hernquist & Katz 1989): adaptive-h symmetrized density (2.16), EOS, the
// symmetrized pressure-gradient acceleration with Monaghan-Gingold artificial viscosity
// (2.20-2.23), and the thermal-energy rate (2.29). EOS is either ideal-gas adiabatic
// (P=(gamma-1) rho u, for tests like Sedov) or locally isothermal (P = rho cs2_iso, LB94).
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
  int N_neigh = 50, h_iter = 4;       // target neighbor number, h iterations per converged call
  double h_tol = 0.2;                  // accept |N_i - N_neigh| < h_tol*N_neigh
};

// adaptive smoothing lengths via the gather neighbor count (HK89 eq. 2.17). Converges h_i to
// within h_tol of the target; per-step (h changes little) a few iterations suffice.
inline void adapt_h(std::vector<Particle>& p, Tree& tree, const SPHParams& par) {
  int n = p.size();
  for (int it = 0; it < par.h_iter; ++it) {
    tree.build();
    int nbad = 0;
#pragma omp parallel for reduction(+ : nbad) schedule(dynamic, 256)
    for (int i = 0; i < n; ++i) {
      int cnt = tree.gather_count(i, 2.0 * p[i].h);
      if (std::abs(cnt - par.N_neigh) > par.h_tol * par.N_neigh) ++nbad;
      double ratio = std::cbrt((double)par.N_neigh / std::max(cnt, 1));
      p[i].h *= 0.5 * (1.0 + std::min(std::max(ratio, 0.5), 1.5));   // wider, faster step
    }
    if (nbad == 0) break;
  }
  tree.build();   // final tree with converged h_max
}

// symmetrized density (eq. 2.16), including the self term m_i W(0,h_i)
inline void density(std::vector<Particle>& p, Tree& tree) {
  int n = p.size();
#pragma omp parallel for schedule(dynamic, 256)
  for (int i = 0; i < n; ++i) {
    std::vector<int> nb;
    double rho = p[i].m * W(0.0, p[i].h);
    tree.neighbors(i, nb);
    for (int j : nb) {
      double r = std::sqrt(dist2(p[i], p[j]));
      rho += p[j].m * 0.5 * (W(r, p[i].h) + W(r, p[j].h));
    }
    p[i].rho = rho;
  }
}

inline void equation_of_state(std::vector<Particle>& p, const SPHParams& par) {
  for (auto& q : p) {
    if (par.isothermal) { q.P = q.rho * q.cs2_iso; q.cs = std::sqrt(q.cs2_iso); }
    else { q.P = (par.gamma - 1.0) * q.rho * q.u; q.cs = std::sqrt(par.gamma * std::max(q.P, 0.0) / q.rho); }
  }
}

// pressure-gradient + artificial-viscosity acceleration and du/dt (eqs. 2.20-2.23, 2.29)
inline void hydro_forces(std::vector<Particle>& p, Tree& tree, const SPHParams& par) {
  int n = p.size();
#pragma omp parallel for schedule(dynamic, 256)
  for (int i = 0; i < n; ++i) {
    std::vector<int> nb;
    double ax = 0, ay = 0, az = 0, dudt = 0, divv = 0;
    Particle& pi = p[i];
    double Pi_over = pi.P / (pi.rho * pi.rho);
    tree.neighbors(i, nb);
    for (int j : nb) {
      Particle& pj = p[j];
      double dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
      double r2 = dx * dx + dy * dy + dz * dz, r = std::sqrt(r2);
      if (r == 0) continue;
      double dvx = pi.vx - pj.vx, dvy = pi.vy - pj.vy, dvz = pi.vz - pj.vz;
      // symmetrized kernel gradient: grad_i Wbar = 0.5[gradW(h_i)+gradW(h_j)] * (r_i - r_j)
      double gw = 0.5 * (gradW(r, pi.h) + gradW(r, pj.h));
      double gwx = gw * dx, gwy = gw * dy, gwz = gw * dz;
      // Monaghan-Gingold artificial viscosity (2.22-2.23)
      double Pij = 0.0, vdotr = dvx * dx + dvy * dy + dvz * dz;
      if (vdotr < 0.0) {
        double hbar = 0.5 * (pi.h + pj.h), cbar = 0.5 * (pi.cs + pj.cs), rhobar = 0.5 * (pi.rho + pj.rho);
        double mu = hbar * vdotr / (r2 + par.eta2_fac * hbar * hbar);
        Pij = (-par.alpha * cbar * mu + par.beta * mu * mu) / rhobar;
      }
      double fac = Pi_over + pj.P / (pj.rho * pj.rho) + Pij;
      ax -= pj.m * fac * gwx; ay -= pj.m * fac * gwy; az -= pj.m * fac * gwz;
      // thermal energy: du/dt = sum m_j (P_i/rho_i^2 + 1/2 Pij) v_ij . grad_i Wbar
      double vdotgw = dvx * gwx + dvy * gwy + dvz * gwz;
      dudt += pj.m * (Pi_over + 0.5 * Pij) * vdotgw;
      divv -= pj.m * vdotgw;                              // rho_i div v = -sum m_j v_ij.gradW
    }
    pi.ax = ax; pi.ay = ay; pi.az = az;
    pi.dudt = dudt; pi.divv = divv / pi.rho;
  }
  if (par.selfgrav)
#pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < n; ++i) {
      double gx, gy, gz; tree.accel(i, gx, gy, gz);
      p[i].ax += gx; p[i].ay += gy; p[i].az += gz;
    }
}

// full force evaluation: h, density, EOS, hydro+gravity accelerations and du/dt
inline void compute_forces(std::vector<Particle>& p, Tree& tree, const SPHParams& par) {
  adapt_h(p, tree, par);
  density(p, tree);
  equation_of_state(p, par);
  hydro_forces(p, tree, par);
}

}  // namespace sph
#endif
