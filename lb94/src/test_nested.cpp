// Validate nested-grid self-gravity: a small dense sphere resolved by the inner grid.
// Compare the central potential from (a) coarse grid only and (b) the nested solve to
// the analytic uniform-sphere value, showing the nesting sharpens the centre.
#include "nested.hpp"
#include <cstdio>
#include <cmath>

int main() {
  using namespace lb94;
  const double G = 1.0, rho0 = 1.0, a = 0.5, L = 4.0;       // sphere a=0.5 fits the finest grid
  const int NR = 64, NZ = 64, nlev = 3;
  NestedGravity ng(nlev, NR, NZ, L, G);

  // sample the sphere on every level
  for (int l = 0; l < nlev; ++l)
    for (int i = 0; i < NR; ++i)
      for (int j = 0; j < NZ; ++j)
        ng.rho[l][i * NZ + j] =
            (std::sqrt(ng.lev[l].R[i] * ng.lev[l].R[i] + ng.lev[l].Z[j] * ng.lev[l].Z[j]) < a)
                ? rho0 : 0.0;
  ng.solve();

  double M = 4.0 / 3.0 * PI * rho0 * a * a * a;
  double phi_an = -2.0 * PI * G * rho0 * a * a;             // analytic centre
  // coarse-only reference (single grid over the whole box)
  Poisson coarse(NR, NZ, L, L, G);
  std::vector<double> rc((size_t)NR * NZ, 0.0), pc;
  for (int i = 0; i < NR; ++i)
    for (int j = 0; j < NZ; ++j)
      if (std::sqrt(coarse.R[i] * coarse.R[i] + coarse.Z[j] * coarse.Z[j]) < a)
        rc[i * NZ + j] = rho0;
  coarse.solve(rc, pc);

  printf("sphere a=%.2f, M=%.4f, finest resolution = a/%.0f cells\n",
         a, M, a / ng.lev[nlev - 1].dR);
  printf("analytic central Phi          = %.5f\n", phi_an);
  printf("coarse-grid (a/%.1f cells)     = %.5f   (err %.1f%%)\n",
         a / coarse.dR, pc[0], 100 * std::fabs(pc[0] - phi_an) / std::fabs(phi_an));
  printf("nested finest grid central    = %.5f   (err %.1f%%)\n",
         ng.Phi[nlev - 1][0], 100 * std::fabs(ng.Phi[nlev - 1][0] - phi_an) / std::fabs(phi_an));
  return 0;
}
