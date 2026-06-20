// Rotating collapse of a cold cloud -> centrifugally supported protostellar disk.
// Uniform sphere (rho0, radius a) in solid-body rotation Omega0, G=1 units.
// Low-angular-momentum polar gas accretes onto the central sink (M_c); the rest
// settles into a rotating midplane disk near the centrifugal radius.
#include "hydro.hpp"
#include <cstdio>
#include <cmath>

using namespace lb94;

int main() {
  const double G = 1.0, rho0 = 1.0, a = 1.0, Rmax = 1.5, Zmax = 1.5;
  const int N = 120;
  const double Omega0 = 0.8, cs2 = 0.05;
  const double M = 4.0 / 3.0 * PI * rho0 * a * a * a;
  const double t_ff = std::sqrt(3.0 * PI / (32.0 * G * rho0));
  const double Rcf = Omega0 * Omega0 * a * a * a * a / (G * M);   // outer-shell centrifugal radius

  // barotropic EOS (Larson): isothermal until rho_crit, then adiabatic -> pressure-supported
  // central core that halts the collapse runaway.
  Hydro h(N, N, Rmax, Zmax, G, cs2, /*rho_crit=*/30.0, /*gam=*/5.0 / 3.0);
  h.rho_active = 1e-3;
  h.R_sink = 0.08; h.rho_ceiling = 30.0;        // density-ceiling sink (resolved core zone)
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) {
      int c = h.idx(i, j);
      if (std::sqrt(h.R[i] * h.R[i] + h.Z[j] * h.Z[j]) < a) {
        h.rho[c] = rho0;
        h.A[c] = rho0 * h.R[i] * (Omega0 * h.R[i]);            // A = rho R v_phi, v_phi=Omega0 R
      } else h.rho[c] = h.floor;
    }

  auto gas_mass = [&]() {
    double m = 0;
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j) m += h.rho[h.idx(i, j)] * 2 * (2 * PI * h.R[i] * h.dR * h.dZ);
    return m;
  };
  double m0 = gas_mass();
  printf("M=%.4f  t_ff=%.4f  centrifugal radius (outer shell)=%.3f  R_sink=%.3f\n",
         m0, t_ff, Rcf, h.R_sink);
  printf("\n t/t_ff   M_core    M_gas   total/M0   disk R_peak(midplane)\n");

  double t = 0;
  int nstep = 0;
  double tmax = 2.5 * t_ff;
  double nextlog = 0.0;
  while (t < tmax && nstep < 200000) {
    h.poisson.solve(h.rho, h.Phi, 3000, 1e-5);
    double dt = std::min(h.timestep(0.3), 0.004 * t_ff);
    h.step(dt); t += dt; ++nstep;
    {
      double rmax = 0, vmax = 0; int ic = -1;
      for (int c = 0; c < N * N && ic < 0; ++c)
        if (!std::isfinite(h.rho[c]) || !std::isfinite(h.sR[c]) || !std::isfinite(h.sZ[c]) ||
            !std::isfinite(h.A[c]) || h.rho[c] < 0) ic = c;
      for (int c = 0; c < N * N; ++c) {
        rmax = std::max(rmax, h.rho[c]);
        vmax = std::max(vmax, std::fabs(h.sR[c]) / std::max(h.rho[c], h.rho_active));
      }
      if (nstep % 25 == 0 || ic >= 0)
        fprintf(stderr, "  step %d t/tff=%.3f dt/tff=%.1e rho_max=%.2e v_max=%.2e M_c=%.4f%s\n",
                nstep, t / t_ff, dt / t_ff, rmax, vmax, h.M_c, ic >= 0 ? "  <-- BAD!" : "");
      if (ic >= 0) {
        int i = ic / N, j = ic % N;
        fprintf(stderr, "  first bad cell %d (i=%d,j=%d) R=%.3f Z=%.3f: rho=%.3e sR=%.3e A=%.3e Phi=%.3e\n",
                ic, i, j, h.R[i], h.Z[j], h.rho[ic], h.sR[ic], h.A[ic], h.Phi[ic]);
        break;
      }
    }
    if (t >= nextlog) {
      // disk peak: max midplane (j=0) gas density and its radius
      int ipk = 0; double rpk = 0;
      for (int i = 0; i < N; ++i)
        if (h.R[i] > h.R_sink && h.rho[h.idx(i, 0)] > rpk) { rpk = h.rho[h.idx(i, 0)]; ipk = i; }
      printf(" %.3f   %.4f   %.4f   %.5f    %.3f\n",
             t / t_ff, h.M_c, gas_mass(), (gas_mass() + h.M_c) / m0, h.R[ipk]);
      fflush(stdout);
      nextlog += 0.25 * t_ff;
    }
  }

  // dump final density field and midplane profile
  FILE* f = fopen("lb94/collapse_rho.dat", "w");
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      fprintf(f, "%.4f %.4f %.6e\n", h.R[i], h.Z[j], h.rho[h.idx(i, j)]);
  fclose(f);
  FILE* g = fopen("lb94/collapse_midplane.dat", "w");
  for (int i = 0; i < N; ++i) {
    double d = h.rho[h.idx(i, 0)];
    double vphi = (d > h.rho_active) ? h.A[h.idx(i, 0)] / (d * h.R[i]) : 0.0;
    double vkep = std::sqrt(G * h.M_c / h.R[i]);               // core-only Keplerian
    fprintf(g, "%.4f %.6e %.4f %.4f\n", h.R[i], d, vphi, vkep);
  }
  fclose(g);
  printf("\nfinal: M_core=%.4f (%.0f%% of M)  steps=%d\n", h.M_c, 100 * h.M_c / m0, nstep);
  printf("wrote lb94/collapse_rho.dat, lb94/collapse_midplane.dat\n");
  return 0;
}
