// YBL93 standard-case collapse in physical (cgs) units, on Explicit Nested Grids.
//   1 Msun, sphere R=4e16 cm (~2700 AU), rho ~ 1/r, isothermal 20 K,
//   uniform Omega=4.4e-13 s^-1.  Compare to YBL93: ~0.45 Msun central object at
//   ~8e4 yr, disk Keplerian radius of the outer element ~250 AU.
#include "nested.hpp"
#include <cstdio>
#include <cmath>

using namespace lb94;

int main() {
  // cgs constants
  const double G = 6.67430e-8, Msun = 1.989e33, AU = 1.496e13, yr = 3.156e7;
  const double kB = 1.380649e-16, mH = 1.6726e-24, mu = 2.34;

  // YBL standard case
  const double Mtot = 1.0 * Msun, Rcloud = 4.0e16;       // 1 Msun, ~2700 AU
  const double Tgas = 20.0, Omega0 = 4.4e-13;
  const double cs2 = kB * Tgas / (mu * mH);              // isothermal sound speed^2
  const double rho_crit = 1.0e-13, gam = 1.4;            // barotropic (first core)
  const int NR = 60, NZ = 60, nlev = 4;
  const double L = 2.0 * Rcloud;                         // box >> cloud: cloud well inside the
                                                         // coarsest grid, disk inside the finest

  const double rho_mean = Mtot / (4.0 / 3.0 * PI * Rcloud * Rcloud * Rcloud);
  const double t_ff = std::sqrt(3.0 * PI / (32.0 * G * rho_mean));
  const double Ccoef = Mtot / (2.0 * PI * Rcloud * Rcloud);   // rho = Ccoef / r
  const double j_out = Omega0 * Rcloud * Rcloud;
  const double Rkep = j_out * j_out / (G * Mtot);

  NestedHydro nh(nlev, NR, NZ, L, G, cs2, rho_crit, gam);
  nh.lev[0].closed_outer = true;        // the bound cloud cannot lose mass through the box edge
  double dRfine = nh.finest().dR;
  for (auto& h : nh.lev) {
    h.R_sink = 7.0 * dRfine; h.rho_active = rho_mean * 1e-4;
    h.vmax = 1.0e7;                                      // 100 km/s velocity ceiling (no mass loss)
    h.use_energy = true; h.radiation = true; h.Tamb = Tgas;   // energy eq + radiative cooling
    h.floor = rho_mean * 1e-3;                           // ~500x below cloud density: stable van Leer
                                                         // advection (much lower -> negative-rho
                                                         // overshoot at the steep cloud/vacuum edge)
  }
  nh.finest().rho_ceiling = 5.0e-14;                     // Truelove-resolved on the finest grid
  nh.finest().t_drain = 300.0 * yr;                      // drain the unresolved central object
  for (int l = 0; l < nlev; ++l) {
    Hydro& h = nh.lev[l];
    for (int i = 0; i < NR; ++i)
      for (int j = 0; j < NZ; ++j) {
        int c = h.idx(i, j);
        double r = std::sqrt(h.R[i] * h.R[i] + h.Z[j] * h.Z[j]);
        if (r < Rcloud) {
          double rr = std::max(r, 0.5 * dRfine);
          h.rho[c] = Ccoef / rr;
          h.A[c] = h.rho[c] * h.R[i] * (Omega0 * h.R[i]);
          h.e[c] = h.rho[c] * cs2 / (gam - 1.0);          // isothermal 20 K internal energy
        } else { h.rho[c] = h.floor; h.e[c] = h.floor * cs2 / (gam - 1.0); }
      }
  }

  auto gas_mass = [&]() { return nh.composite_mass(); };   // finest-resolved, not lagged
  double m0 = gas_mass();
  printf("YBL standard case:  M=%.3f Msun  Rcloud=%.0f AU  cs=%.3f km/s  t_ff=%.0f yr\n",
         m0 / Msun, Rcloud / AU, std::sqrt(cs2) / 1e5, t_ff / yr);
  printf("  finest dR=%.2f AU (nlev=%d)  outer Keplerian radius=%.0f AU\n",
         dRfine / AU, nlev, Rkep / AU);
  printf("\n t(yr)    t/t_ff   M_core(Msun)  M_disk(Msun)  cons    Rout(AU)  T_disk(K)\n");

  double t = 0, nextlog = 0; int nstep = 0;
  double tmax = 6.2e4 * yr;                              // the LB94 SPH-initial epoch (~60 kyr)
  while (t < tmax && nstep < 400000) {
    double dt = std::min(nh.timestep(0.3), 2.0e-3 * t_ff);
    nh.step(dt); t += dt; ++nstep;     // refluxing (in step) conserves; no renormalization
    {
      int bl = -1, bc = -1;
      for (int l = 0; l < nlev && bl < 0; ++l) {
        Hydro& h = nh.lev[l];
        for (int c = 0; c < NR * NZ; ++c)
          if (!std::isfinite(h.rho[c]) || !std::isfinite(h.sR[c]) || !std::isfinite(h.sZ[c]) ||
              !std::isfinite(h.A[c]) || !std::isfinite(h.Phi[c])) { bl = l; bc = c; break; }
      }
      if (bl >= 0) {
        Hydro& h = nh.lev[bl]; int i = bc / NZ, j = bc % NZ;
        fprintf(stderr, "first NaN t=%.0fyr step=%d: LEVEL %d cell(%d,%d) R=%.2fAU Z=%.2fAU\n"
                "   rho=%.3e sR=%.3e sZ=%.3e A=%.3e Phi=%.3e  (dR_lev=%.1fAU, M_c=%.3e)\n",
                t / yr, nstep, bl, i, j, h.R[i] / AU, h.Z[j] / AU,
                h.rho[bc], h.sR[bc], h.sZ[bc], h.A[bc], h.Phi[bc], h.dR / AU, h.M_c);
        break;
      }
    }
    if (!std::isfinite(nh.finest().rho[0])) { printf("  [blowup t=%.0f yr]\n", t / yr); break; }
    if (t >= nextlog) {
      // disk outer radius + midplane temperature at the densest disk cell
      Hydro& f = nh.finest(); double Rout = 0, rpk = 0, Tdisk = 0;
      double cv = kB / ((gam - 1.0) * mu * mH);
      for (int i = 0; i < NR; ++i) {
        double d = f.rho[f.idx(i, 0)];
        if (d > 10.0 * rho_mean) Rout = f.R[i];
        if (f.R[i] > f.R_sink && d > rpk) { rpk = d; Tdisk = f.e[f.idx(i, 0)] / (d * cv); }
      }
      double floored = 0;                               // fictitious mass added by the vacuum floor
      for (auto& h : nh.lev) floored += h.dbg_mass_floored;
      printf(" %.2e  %.3f   %.4f       %.4f      %.4f   %.0f      %.0f\n",
             t / yr, t / t_ff, nh.finest().M_c / Msun, gas_mass() / Msun,
             (gas_mass() + nh.finest().M_c - floored) / m0, Rout / AU, Tdisk);
      fflush(stdout); nextlog += 0.5e4 * yr;
    }
  }
  Hydro& f = nh.finest();
  FILE* fp = fopen("lb94/ybl_rho.dat", "w");
  for (int i = 0; i < NR; ++i)
    for (int j = 0; j < NZ; ++j)
      fprintf(fp, "%.4e %.4e %.6e\n", f.R[i] / AU, f.Z[j] / AU, f.rho[f.idx(i, j)]);
  fclose(fp);
  double floored = 0;
  for (auto& h : nh.lev) floored += h.dbg_mass_floored;
  printf("\nfinal: M_core=%.3f Msun (%.0f%%) at t=%.0f yr,  steps=%d\n",
         nh.finest().M_c / Msun, 100 * nh.finest().M_c / m0, t / yr, nstep);
  printf("       cumulative mass added by density floor = %.4f Msun (%.2f%% of m0)\n",
         floored / Msun, 100 * floored / m0);
  return 0;
}
