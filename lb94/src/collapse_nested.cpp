// Rotating collapse on Explicit Nested Grids: same physics as collapse.cpp but the
// inner region is resolved by nested grids (factor 2 each), so the disk/core sit on
// the finest grid.  Validates the nested hydro (restriction/prolongation + nested gravity).
#include "nested.hpp"
#include <cstdio>
#include <cmath>

using namespace lb94;

int main() {
  const double G = 1.0, rho0 = 1.0, a = 1.0, L = 1.5;       // coarsest box [0,1.5]^2
  const int NR = 80, NZ = 80, nlev = 3;                     // finest covers [0,0.375]
  const double Omega0 = 0.8, cs2 = 0.05;
  const double M = 4.0 / 3.0 * PI * rho0 * a * a * a;
  const double t_ff = std::sqrt(3.0 * PI / (32.0 * G * rho0));
  const double Rcf = Omega0 * Omega0 * a * a * a * a / (G * M);

  NestedHydro nh(nlev, NR, NZ, L, G, cs2, /*rho_crit=*/30.0, /*gam=*/5.0 / 3.0);
  for (auto& h : nh.lev) { h.R_sink = 0.08; h.rho_active = 1e-3; }
  nh.finest().rho_ceiling = 30.0;                           // only the finest accretes
  for (int l = 0; l < nlev; ++l) {
    Hydro& h = nh.lev[l];
    for (int i = 0; i < NR; ++i)
      for (int j = 0; j < NZ; ++j) {
        int c = h.idx(i, j);
        if (std::sqrt(h.R[i] * h.R[i] + h.Z[j] * h.Z[j]) < a) {
          h.rho[c] = rho0; h.A[c] = rho0 * h.R[i] * (Omega0 * h.R[i]);
        } else h.rho[c] = h.floor;
      }
  }

  auto gas_mass = [&]() {                                   // from the coarsest grid (full domain)
    Hydro& h = nh.lev[0]; double m = 0;
    for (int i = 0; i < NR; ++i)
      for (int j = 0; j < NZ; ++j) m += h.rho[h.idx(i, j)] * 2 * (2 * PI * h.R[i] * h.dR * h.dZ);
    return m;
  };
  double m0 = gas_mass();
  printf("M=%.4f  t_ff=%.4f  R_cf=%.3f  nlev=%d  finest dR=%.4f (single-grid dR=%.4f)\n",
         m0, t_ff, Rcf, nlev, nh.finest().dR, L / NR);
  printf("\n t/t_ff   M_core    M_gas   total/M0   disk R_peak(finest)\n");

  double t = 0, nextlog = 0; int nstep = 0;
  while (t < 2.0 * t_ff && nstep < 400000) {
    double dt = std::min(nh.timestep(0.3), 0.004 * t_ff);
    nh.step(dt); t += dt; ++nstep;
    if (!std::isfinite(nh.finest().rho[0])) { printf("  [blowup t/tff=%.3f]\n", t / t_ff); break; }
    if (t >= nextlog) {
      Hydro& f = nh.finest(); int ipk = 0; double rpk = 0;
      for (int i = 0; i < NR; ++i)
        if (f.R[i] > f.R_sink && f.rho[f.idx(i, 0)] > rpk) { rpk = f.rho[f.idx(i, 0)]; ipk = i; }
      printf(" %.3f   %.4f   %.4f   %.5f    %.4f\n",
             t / t_ff, nh.finest().M_c, gas_mass(), (gas_mass() + nh.finest().M_c) / m0, f.R[ipk]);
      fflush(stdout); nextlog += 0.25 * t_ff;
    }
  }
  // dump finest-grid density field
  Hydro& f = nh.finest();
  FILE* fp = fopen("lb94/collapse_nested_rho.dat", "w");
  for (int i = 0; i < NR; ++i)
    for (int j = 0; j < NZ; ++j)
      fprintf(fp, "%.5f %.5f %.6e\n", f.R[i], f.Z[j], f.rho[f.idx(i, j)]);
  fclose(fp);
  printf("\nfinal: M_core=%.4f (%.0f%%)  steps=%d  wrote lb94/collapse_nested_rho.dat\n",
         nh.finest().M_c, 100 * nh.finest().M_c / m0, nstep);
  return 0;
}
