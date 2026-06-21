// 2D axisymmetric cylindrical (R,Z) self-gravitating hydrodynamics.
// State (cell-centred, conservative): rho, sR=rho v_R, sZ=rho v_Z, A=rho R v_phi.
// Operator split: gravity/pressure/centrifugal source -> directionally-split
// upwind transport.  Barotropic EOS (isothermal Stage 1a; barotropic switch ready).
// Axis (R=0) and midplane (Z=0) reflect; outer R,Z edges outflow.
#ifndef LB94_HYDRO_HPP
#define LB94_HYDRO_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include "poisson.hpp"

namespace lb94 {

struct Hydro {
  int NR, NZ;
  double dR, dZ, G;
  std::vector<double> R, Z;
  std::vector<double> rho, sR, sZ, A, Phi;
  std::vector<double> e;              // internal energy density (when use_energy)
  bool use_energy = false;           // false: barotropic EOS; true: ideal gas with energy eq
  double emin = 1e-300;
  // nested-grid ghost state at the outer R and Z boundaries ([rho,sR,sZ,A] each)
  bool ghostBC = false;
  std::vector<std::vector<double>> ghostR{4}, ghostZ{4};
  Poisson poisson;
  // barotropic EOS: P = cs2 * rho  (rho < rho_crit);  P = cs2*rho_crit*(rho/rho_crit)^gam above
  double cs2, rho_crit, gam;
  double floor = 1e-20;
  double rho_active = 1e-3;            // below this a cell is "vacuum": carries no momentum
  inline double vdiv(double d) const { return std::max(d, rho_active); }

  Hydro(int NR_, int NZ_, double Rmax, double Zmax, double G_ = 6.67430e-8,
        double cs2_ = 1.0, double rho_crit_ = 1e30, double gam_ = 1.4)
      : NR(NR_), NZ(NZ_), dR(Rmax / NR_), dZ(Zmax / NZ_), G(G_),
        poisson(NR_, NZ_, Rmax, Zmax, G_), cs2(cs2_), rho_crit(rho_crit_), gam(gam_) {
    R = poisson.R; Z = poisson.Z;
    size_t n = (size_t)NR * NZ;
    rho.assign(n, floor); sR.assign(n, 0); sZ.assign(n, 0); A.assign(n, 0); Phi.assign(n, 0);
    e.assign(n, 0.0);
    for (int k = 0; k < 4; ++k) { ghostR[k].assign(NZ, 0.0); ghostZ[k].assign(NR, 0.0); }
  }
  double Cq = 1.4;                    // von Neumann artificial-viscosity coefficient (R85)
  // central unresolved core/sink (BYRT90): within R_sink the density is capped at
  // rho_ceiling (a resolved value, Jeans length > few cells) and the excess mass+angular
  // momentum is accreted onto the core M_c, J_c.  Avoids the under-resolved Truelove runaway.
  double R_sink = 0.0, rho_ceiling = 0.0, M_c = 0.0, J_c = 0.0, Mdot = 0.0;
  double rho_max = 1e300;             // hard density clamp (stabilizes under-resolved coarse
                                      // overlap cells, which are overwritten by restriction)
  inline int idx(int i, int j) const { return i * NZ + j; }
  static inline double vanleer(double a, double b) {
    return (a * b > 0.0) ? 2.0 * a * b / (a + b) : 0.0;
  }
  inline double pres(double d) const {
    return (d < rho_crit) ? cs2 * d : cs2 * rho_crit * std::pow(d / rho_crit, gam);
  }
  inline double csound(double d) const {
    return (d < rho_crit) ? std::sqrt(cs2) : std::sqrt(gam * pres(d) / d);
  }
  // pressure / sound speed per cell (ideal gas P=(gam-1)e when use_energy, else barotropic)
  inline double Pcell(int c) const { return use_energy ? (gam - 1.0) * e[c] : pres(rho[c]); }
  inline double cs_cell(int c) const {
    return use_energy ? std::sqrt(gam * std::max((gam - 1.0) * e[c], 0.0) / vdiv(rho[c]))
                      : csound(rho[c]);
  }

  double timestep(double cfl) const {
    double dtinv = 1e-300;
    for (int i = 0; i < NR; ++i)
      for (int j = 0; j < NZ; ++j) {
        double d = rho[idx(i, j)];
        if (d < rho_active) continue;                 // ignore vacuum cells
        double cs = cs_cell(idx(i, j));
        double vr = std::fabs(sR[idx(i, j)] / d), vz = std::fabs(sZ[idx(i, j)] / d);
        double Om = std::fabs(A[idx(i, j)] / (d * R[i] * R[i]));   // angular frequency v_phi/R
        double wff = std::sqrt(4.0 * PI * G * d);                  // local free-fall rate
        dtinv = std::max(dtinv, (vr + cs) / dR + (vz + cs) / dZ + Om + wff);
      }
    return cfl / dtinv;
  }

  void update_gravity() { poisson.solve(rho, Phi); }

  // von Neumann artificial-viscosity pressure (compression only): Q = Cq rho (L div v)^2
  void avisc_pressure(std::vector<double>& Q) const {
    double L = std::min(dR, dZ);
    Q.assign((size_t)NR * NZ, 0.0);
    for (int i = 0; i < NR; ++i)
      for (int j = 0; j < NZ; ++j) {
        double d = rho[idx(i, j)];
        if (d < rho_active) continue;
        int ip = std::min(i + 1, NR - 1), im = std::max(i - 1, 0);
        int jp = std::min(j + 1, NZ - 1), jm = std::max(j - 1, 0);
        double vRp = sR[idx(ip, j)] / vdiv(rho[idx(ip, j)]), vRm = sR[idx(im, j)] / vdiv(rho[idx(im, j)]);
        double vZp = sZ[idx(i, jp)] / vdiv(rho[idx(i, jp)]), vZm = sZ[idx(i, jm)] / vdiv(rho[idx(i, jm)]);
        double vR = sR[idx(i, j)] / vdiv(d);
        double divv = (vRp - vRm) / ((ip - im) * dR) + vR / R[i] + (vZp - vZm) / ((jp - jm) * dZ);
        if (divv < 0) Q[idx(i, j)] = Cq * d * (L * divv) * (L * divv);
      }
  }

  // --- source: pressure(+artificial viscosity) gradient, gravity, centrifugal,
  //     and (use_energy) the compression + artificial-viscosity heating of e ---
  void source(double dt) {
    std::vector<double> nsR = sR, nsZ = sZ, ne = e, Q;
    avisc_pressure(Q);
    auto Peff = [&](int i, int j) { return Pcell(idx(i, j)) + Q[idx(i, j)]; };
    for (int i = 0; i < NR; ++i)
      for (int j = 0; j < NZ; ++j) {
        int c = idx(i, j);
        double d = std::max(rho[c], floor);
        int ip = std::min(i + 1, NR - 1), im = std::max(i - 1, 0);
        int jp = std::min(j + 1, NZ - 1), jm = std::max(j - 1, 0);
        double dPdR = (Peff(ip, j) - Peff(im, j)) / ((ip - im) * dR);
        double dPdZ = (Peff(i, jp) - Peff(i, jm)) / ((jp - jm) * dZ);
        double dPhidR = (Phi[idx(ip, j)] - Phi[idx(im, j)]) / ((ip - im) * dR);
        double dPhidZ = (Phi[idx(i, jp)] - Phi[idx(i, jm)]) / ((jp - jm) * dZ);
        double vphi = A[c] / (vdiv(d) * R[i]);                 // v_phi = A/(rho R)
        // central core point-mass gravity (Plummer-softened at R_sink)
        double gcR = 0.0, gcZ = 0.0;
        if (M_c > 0) {
          double s = std::pow(R[i] * R[i] + Z[j] * Z[j] + R_sink * R_sink, 1.5);
          gcR = -G * M_c * R[i] / s; gcZ = -G * M_c * Z[j] / s;
        }
        nsR[c] += dt * (-dPdR - d * dPhidR + d * gcR + d * vphi * vphi / R[i]);
        nsZ[c] += dt * (-dPdZ - d * dPhidZ + d * gcZ);
        if (use_energy) {                                     // de = -(P+q) div(v) dt
          double vRp = sR[idx(ip, j)] / vdiv(rho[idx(ip, j)]), vRm = sR[idx(im, j)] / vdiv(rho[idx(im, j)]);
          double vZp = sZ[idx(i, jp)] / vdiv(rho[idx(i, jp)]), vZm = sZ[idx(i, jm)] / vdiv(rho[idx(i, jm)]);
          double divv = (vRp - vRm) / ((ip - im) * dR) + (sR[c] / vdiv(d)) / R[i]
                      + (vZp - vZm) / ((jp - jm) * dZ);
          ne[c] = std::max(e[c] - (Pcell(c) + Q[c]) * divv * dt, emin);
        }
      }
    sR.swap(nsR); sZ.swap(nsZ);
    if (use_energy) e.swap(ne);
  }

  // density-based sink: any cell above rho_ceiling (the Truelove-unresolved gas, i.e. the
  // central object + any inner-edge pile-up) is capped and its excess mass + angular
  // momentum accreted onto the core.  No hard radius edge -> no edge pile-up.
  void accrete(double dt) {
    if (rho_ceiling <= 0) return;
    double dM = 0.0, dJ = 0.0;
    for (int i = 0; i < NR; ++i)
      for (int j = 0; j < NZ; ++j) {
        int c = idx(i, j);
        if (rho[c] <= rho_ceiling) continue;
        double f = rho_ceiling / rho[c];                       // keep this fraction
        double vol = 2.0 * (2.0 * PI * R[i] * dR * dZ);        // full ring (mirror Z<0)
        dM += rho[c] * (1.0 - f) * vol; dJ += A[c] * (1.0 - f) * vol;
        rho[c] = rho_ceiling; sR[c] *= f; sZ[c] *= f; A[c] *= f;   // keep velocities
      }
    M_c += dM; J_c += dJ; Mdot = dM / dt;
  }

  // van Leer 2nd-order conservative transport (R then Z), face velocities fixed per sweep
  void transport(double dt) {
    std::vector<double>* q[5] = {&rho, &sR, &sZ, &A, &e};
    int nq = use_energy ? 5 : 4;
    auto ghR = [&](int k, int j) { return k < 4 ? ghostR[k][j] : e[idx(NR - 1, j)]; };
    auto ghZ = [&](int k, int i) { return k < 4 ? ghostZ[k][i] : e[idx(i, NZ - 1)]; };

    // R-sweep:  d_t q + (1/R) d_R(R q v_R) = 0
    std::vector<double> vc((size_t)NR * NZ);
    for (int c = 0; c < NR * NZ; ++c) vc[c] = sR[c] / vdiv(rho[c]);
    std::vector<double> vgR(NZ, 0.0);                      // ghost face velocity (outer R)
    if (ghostBC)
      for (int j = 0; j < NZ; ++j)
        vgR[j] = 0.5 * (vc[idx(NR - 1, j)] + ghostR[1][j] / vdiv(ghostR[0][j]));
    for (int k = 0; k < nq; ++k) {
      std::vector<double>& Q = *q[k];
      std::vector<double> nQ = Q;
      for (int j = 0; j < NZ; ++j) {
        std::vector<double> F(NR + 1, 0.0);
        for (int i = 0; i < NR - 1; ++i) {
          double vf = 0.5 * (vc[idx(i, j)] + vc[idx(i + 1, j)]);
          double qf;
          if (vf > 0) {
            double sl = (i > 0) ? vanleer(Q[idx(i + 1, j)] - Q[idx(i, j)],
                                          Q[idx(i, j)] - Q[idx(i - 1, j)]) : 0.0;
            qf = Q[idx(i, j)] + 0.5 * sl * (1.0 - vf * dt / dR);
          } else {
            double sl = (i + 1 < NR - 1) ? vanleer(Q[idx(i + 2, j)] - Q[idx(i + 1, j)],
                                                   Q[idx(i + 1, j)] - Q[idx(i, j)]) : 0.0;
            qf = Q[idx(i + 1, j)] - 0.5 * sl * (1.0 + vf * dt / dR);
          }
          F[i + 1] = (i + 1) * dR * vf * qf;                  // R_{i+1/2} * v * q
        }
        if (ghostBC) {
          double vf = vgR[j];
          double qf = (vf > 0) ? Q[idx(NR - 1, j)] : ghR(k, j);      // inflow from parent
          F[NR] = NR * dR * vf * qf;
        } else
          F[NR] = NR * dR * std::max(vc[idx(NR - 1, j)], 0.0) * Q[idx(NR - 1, j)];
        for (int i = 0; i < NR; ++i)
          nQ[idx(i, j)] -= dt / (R[i] * dR) * (F[i + 1] - F[i]);
      }
      Q.swap(nQ);
    }
    // Z-sweep:  d_t q + d_Z(q v_Z) = 0
    for (int c = 0; c < NR * NZ; ++c) vc[c] = sZ[c] / vdiv(rho[c]);
    std::vector<double> vgZ(NR, 0.0);                      // ghost face velocity (outer Z)
    if (ghostBC)
      for (int i = 0; i < NR; ++i)
        vgZ[i] = 0.5 * (vc[idx(i, NZ - 1)] + ghostZ[2][i] / vdiv(ghostZ[0][i]));
    for (int k = 0; k < nq; ++k) {
      std::vector<double>& Q = *q[k];
      std::vector<double> nQ = Q;
      for (int i = 0; i < NR; ++i) {
        std::vector<double> F(NZ + 1, 0.0);
        for (int j = 0; j < NZ - 1; ++j) {
          double vf = 0.5 * (vc[idx(i, j)] + vc[idx(i, j + 1)]);
          double qf;
          if (vf > 0) {
            double sl = (j > 0) ? vanleer(Q[idx(i, j + 1)] - Q[idx(i, j)],
                                          Q[idx(i, j)] - Q[idx(i, j - 1)]) : 0.0;
            qf = Q[idx(i, j)] + 0.5 * sl * (1.0 - vf * dt / dZ);
          } else {
            double sl = (j + 1 < NZ - 1) ? vanleer(Q[idx(i, j + 2)] - Q[idx(i, j + 1)],
                                                   Q[idx(i, j + 1)] - Q[idx(i, j)]) : 0.0;
            qf = Q[idx(i, j + 1)] - 0.5 * sl * (1.0 + vf * dt / dZ);
          }
          F[j + 1] = vf * qf;
        }
        F[0] = 0.0;                                            // midplane reflect: no flux
        if (ghostBC) {
          double vf = vgZ[i];
          double qf = (vf > 0) ? Q[idx(i, NZ - 1)] : ghZ(k, i);
          F[NZ] = vf * qf;
        } else
          F[NZ] = std::max(vc[idx(i, NZ - 1)], 0.0) * Q[idx(i, NZ - 1)];  // outflow
        for (int j = 0; j < NZ; ++j)
          nQ[idx(i, j)] -= dt / dZ * (F[j + 1] - F[j]);
      }
      Q.swap(nQ);
    }
    // floor/clamp density and clear momentum in vacuum cells
    for (int c = 0; c < NR * NZ; ++c) {
      if (rho[c] < floor) rho[c] = floor;
      if (rho[c] > rho_max) {                              // clamp (keep velocity)
        double f = rho_max / rho[c]; sR[c] *= f; sZ[c] *= f; A[c] *= f; rho[c] = rho_max;
      }
      if (rho[c] < rho_active) { sR[c] = sZ[c] = A[c] = 0.0; }
      if (use_energy && e[c] < emin) e[c] = emin;
    }
  }

  void step(double dt) { source(dt); transport(dt); accrete(dt); }
};

}  // namespace lb94
#endif
