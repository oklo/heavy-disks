// Sedov-Taylor blast wave: validate the energy equation (adiabatic shock).
// Deposit energy E0 at the origin in a cold uniform medium; the spherical blast
// expands as R_S = xi0 (E0 t^2 / rho0)^(1/5) ~ t^0.4 (3D, gamma=5/3).
// (The axisymmetric R,Z grid with rotation about the axis represents 3D.)
#include "hydro.hpp"
#include <cstdio>
#include <cmath>

using namespace lb94;

int main() {
  const double rho0 = 1.0, E0 = 1.0, gamma = 5.0 / 3.0, L = 1.0;
  const int N = 160;
  Hydro h(N, N, L, L, /*G=*/0.0, /*cs2=*/0.0, /*rho_crit=*/1e300, /*gam=*/gamma);
  h.use_energy = true;
  h.rho_active = 1e-6;

  const double r_blast = 0.06;
  double Vb = 0.0;
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      if (std::sqrt(h.R[i] * h.R[i] + h.Z[j] * h.Z[j]) < r_blast)
        Vb += 2.0 * (2.0 * PI * h.R[i] * h.dR * h.dZ);     // 3D volume of blast region
  double e_amb = 1e-5;
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) {
      int c = h.idx(i, j);
      h.rho[c] = rho0;
      h.e[c] = (std::sqrt(h.R[i] * h.R[i] + h.Z[j] * h.Z[j]) < r_blast) ? E0 / Vb : e_amb;
    }

  auto shock_radius = [&]() {                              // largest midplane R with rho>1.3 rho0
    double Rs = 0;
    for (int i = 0; i < N; ++i)
      if (h.rho[h.idx(i, 0)] > 1.3 * rho0) Rs = h.R[i];
    return Rs;
  };
  double xi0 = 1.15;                                       // Sedov constant for gamma=5/3
  printf("Sedov blast: E0=%.1f rho0=%.1f gamma=%.3f, blast region r<%.3f (V=%.4f)\n",
         E0, rho0, gamma, r_blast, Vb);
  printf("\n   t       R_S(sim)   R_S(Sedov)   ratio\n");

  double t = 0; int nstep = 0;
  double samples[] = {0.03, 0.05, 0.08, 0.11};
  int ns = 0;
  std::vector<std::pair<double, double>> track;
  while (t < 0.12 && ns < 4 && nstep < 200000) {
    double dt = h.timestep(0.3);
    if (t + dt > samples[ns]) dt = samples[ns] - t;
    h.step(dt); t += dt; ++nstep;
    if (t >= samples[ns] - 1e-12) {
      double Rs = shock_radius();
      double Rsed = xi0 * std::pow(E0 * t * t / rho0, 0.2);
      printf("  %.3f    %.4f      %.4f     %.3f\n", t, Rs, Rsed, Rs / Rsed);
      track.push_back({t, Rs}); ns++;
    }
  }
  // fit d log R_S / d log t  (target 0.4)
  int n = track.size();
  double sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (auto& p : track) { double x = std::log(p.first), y = std::log(p.second);
    sx += x; sy += y; sxx += x * x; sxy += x * y; }
  double slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
  printf("\nd log R_S / d log t = %.3f   (Sedov: 0.400; YBL got 0.39)\n", slope);
  return 0;
}
