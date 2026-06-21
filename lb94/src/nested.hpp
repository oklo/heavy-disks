// Explicit Nested Grids (ENG) self-gravity (Yorke & Kaisig style): n self-similar
// grids sharing the origin, each NR x NZ, level l covering [0, L/2^l]^2 (factor-2
// refinement).  Density is restricted fine->coarse (conservative), then the potential
// is solved coarse->fine, each interior grid taking its outer-boundary potential by
// interpolation from its parent.
#ifndef LB94_NESTED_HPP
#define LB94_NESTED_HPP

#include <vector>
#include <cmath>
#include "poisson.hpp"
#include "hydro.hpp"

namespace lb94 {

// bilinear interpolation of cell-centred field f on grid g at physical point (Rt,Zt)
inline double bilerp(const Poisson& g, const std::vector<double>& f, double Rt, double Zt) {
  double fi = Rt / g.dR - 0.5, fj = Zt / g.dZ - 0.5;
  int i0 = (int)std::floor(fi), j0 = (int)std::floor(fj);
  double a = fi - i0, b = fj - j0;
  auto get = [&](int i, int j) {
    i = std::min(std::max(i, 0), g.NR - 1); j = std::min(std::max(j, 0), g.NZ - 1);
    return f[i * g.NZ + j];
  };
  return (1 - a) * (1 - b) * get(i0, j0) + a * (1 - b) * get(i0 + 1, j0) +
         (1 - a) * b * get(i0, j0 + 1) + a * b * get(i0 + 1, j0 + 1);
}

struct NestedGravity {
  int nlev, NR, NZ;
  std::vector<Poisson> lev;                       // 0 = coarsest ... nlev-1 = finest
  std::vector<std::vector<double>> rho, Phi;

  NestedGravity(int nlev_, int NR_, int NZ_, double L, double G = 6.67430e-8)
      : nlev(nlev_), NR(NR_), NZ(NZ_) {
    for (int l = 0; l < nlev; ++l) {
      double sz = L / std::pow(2.0, l);
      lev.emplace_back(NR_, NZ_, sz, sz, G);
    }
    rho.assign(nlev, std::vector<double>((size_t)NR * NZ, 0.0));
    Phi.assign(nlev, std::vector<double>((size_t)NR * NZ, 0.0));
  }

  // conservative density restriction over the overlap (inner quarter of each parent)
  void restrict_density() {
    for (int l = nlev - 1; l >= 1; --l) {
      Poisson& f = lev[l];
      for (int ip = 0; ip < NR / 2; ++ip)
        for (int jp = 0; jp < NZ / 2; ++jp) {
          double mass = 0, vol = 0;
          for (int a = 0; a < 2; ++a)
            for (int b = 0; b < 2; ++b) {
              int i = 2 * ip + a, j = 2 * jp + b;
              double vf = f.R[i] * f.dR * f.dZ;   // ~ ring volume / 2pi
              mass += rho[l][i * NZ + j] * vf; vol += vf;
            }
          rho[l - 1][ip * NZ + jp] = mass / vol;
        }
    }
  }

  void solve() {
    restrict_density();
    lev[0].solve(rho[0], Phi[0]);                 // coarsest: isolated multipole boundary
    for (int l = 1; l < nlev; ++l) {
      Poisson& c = lev[l];
      std::vector<double> gR(NZ), gZ(NR);
      for (int j = 0; j < NZ; ++j) gR[j] = bilerp(lev[l - 1], Phi[l - 1], c.R[NR - 1] + c.dR, c.Z[j]);
      for (int i = 0; i < NR; ++i) gZ[i] = bilerp(lev[l - 1], Phi[l - 1], c.R[i], c.Z[NZ - 1] + c.dZ);
      c.solve_dirichlet(rho[l], Phi[l], gR, gZ);
    }
  }
};

// ---- Explicit Nested Grids hydrodynamics (synchronous time stepping) ----
struct NestedHydro {
  int nlev, NR, NZ;
  double L, G;
  std::vector<Hydro> lev;                          // 0 coarsest ... nlev-1 finest

  NestedHydro(int nlev_, int NR_, int NZ_, double L_, double G_,
              double cs2, double rho_crit, double gam)
      : nlev(nlev_), NR(NR_), NZ(NZ_), L(L_), G(G_) {
    lev.reserve(nlev);
    for (int l = 0; l < nlev; ++l) {
      double sz = L / std::pow(2.0, l);
      lev.emplace_back(NR_, NZ_, sz, sz, G_, cs2, rho_crit, gam);
    }
  }
  Hydro& finest() { return lev[nlev - 1]; }

  // conservative restriction of all fluid variables fine->coarse over the overlap
  void restrict_fluid() {
    for (int l = nlev - 1; l >= 1; --l) {
      Hydro& f = lev[l]; Hydro& p = lev[l - 1];
      for (int ip = 0; ip < NR / 2; ++ip)
        for (int jp = 0; jp < NZ / 2; ++jp) {
          double vol = 0, mr = 0, mu = 0, mw = 0, ma = 0, me = 0;
          for (int a = 0; a < 2; ++a)
            for (int b = 0; b < 2; ++b) {
              int c = (2 * ip + a) * NZ + (2 * jp + b);
              double vf = f.R[2 * ip + a] * f.dR * f.dZ;
              vol += vf; mr += f.rho[c] * vf; mu += f.sR[c] * vf;
              mw += f.sZ[c] * vf; ma += f.A[c] * vf; me += f.e[c] * vf;
            }
          int pc = ip * NZ + jp;
          p.rho[pc] = mr / vol; p.sR[pc] = mu / vol; p.sZ[pc] = mw / vol; p.A[pc] = ma / vol;
          p.e[pc] = me / vol;
        }
    }
  }

  // fill each fine grid's outer R,Z ghost state from its parent (conservative inflow BC)
  void prolong_boundary() {
    for (int l = 1; l < nlev; ++l) {
      Hydro& c = lev[l]; Hydro& p = lev[l - 1];
      const Poisson& pg = p.poisson;
      std::vector<double>* pf[4] = {&p.rho, &p.sR, &p.sZ, &p.A};
      c.ghostBC = true;
      for (int k = 0; k < 4; ++k) {
        for (int j = 0; j < NZ; ++j) c.ghostR[k][j] = bilerp(pg, *pf[k], c.R[NR - 1] + c.dR, c.Z[j]);
        for (int i = 0; i < NR; ++i) c.ghostZ[k][i] = bilerp(pg, *pf[k], c.R[i], c.Z[NZ - 1] + c.dZ);
      }
    }
  }

  void solve_gravity() {
    lev[0].poisson.solve(lev[0].rho, lev[0].Phi);
    for (int l = 1; l < nlev; ++l) {
      Hydro& c = lev[l]; const Poisson& pg = lev[l - 1].poisson;
      std::vector<double>& pPhi = lev[l - 1].Phi;
      std::vector<double> gR(NZ), gZ(NR);
      for (int j = 0; j < NZ; ++j) gR[j] = bilerp(pg, pPhi, c.R[NR - 1] + c.dR, c.Z[j]);
      for (int i = 0; i < NR; ++i) gZ[i] = bilerp(pg, pPhi, c.R[i], c.Z[NZ - 1] + c.dZ);
      c.poisson.solve_dirichlet(c.rho, c.Phi, gR, gZ);
    }
  }

  double timestep(double cfl) {
    double dt = 1e300;
    for (auto& h : lev) dt = std::min(dt, h.timestep(cfl));
    return dt;
  }

  void step(double dt) {
    double Mc = finest().M_c, Jc = finest().J_c;        // shared central core
    for (auto& h : lev) { h.M_c = Mc; h.J_c = Jc; }
    restrict_fluid();                                   // coarse reflects fine
    solve_gravity();                                    // nested potential
    prolong_boundary();                                 // fine boundaries from parent
    for (auto& h : lev) h.step(dt);                     // finest also accretes -> M_c
  }
};

}  // namespace lb94
#endif
