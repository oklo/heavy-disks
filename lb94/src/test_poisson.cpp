// Validate the cylindrical Poisson solver against a uniform sphere (analytic):
//   Phi_in(r)  = -2 pi G rho0 (a^2 - r^2/3),   Phi_out(r) = -G M / r,  M = 4/3 pi rho0 a^3.
#include "poisson.hpp"
#include <cstdio>
#include <cmath>

int main() {
  using namespace lb94;
  const double G = 1.0, rho0 = 1.0, a = 1.0;        // dimensionless units
  const double Rmax = 4.0, Zmax = 4.0;              // box >> sphere so multipole BC is clean
  int N = 128;
  Poisson p(N, N, Rmax, Zmax, G);

  std::vector<double> rho((size_t)N * N, 0.0), Phi;
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      if (std::sqrt(p.R[i] * p.R[i] + p.Z[j] * p.Z[j]) < a) rho[p.idx(i, j)] = rho0;

  int iters = p.solve(rho, Phi, 50000, 1e-9);
  double M = 4.0 / 3.0 * PI * rho0 * a * a * a;

  auto analytic = [&](double r) {
    return (r < a) ? -2.0 * PI * G * rho0 * (a * a - r * r / 3.0) : -G * M / r;
  };
  // relative error sampled across the grid
  double l2 = 0, ref = 0, maxrel = 0;
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) {
      double r = std::sqrt(p.R[i] * p.R[i] + p.Z[j] * p.Z[j]);
      double an = analytic(r), nu = Phi[p.idx(i, j)];
      l2 += (nu - an) * (nu - an); ref += an * an;
      maxrel = std::max(maxrel, std::fabs(nu - an) / std::fabs(an));
    }
  printf("SOR iters=%d\n", iters);
  printf("center  Phi: numeric=%.6f  analytic=%.6f\n", Phi[p.idx(0, 0)], analytic(0.0));
  printf("at r~a   Phi: numeric=%.6f  analytic=%.6f\n",
         Phi[p.idx(N / 4, 0)], analytic(p.R[N / 4]));
  printf("L2 relative error = %.3e   max pointwise rel = %.3e\n",
         std::sqrt(l2 / ref), maxrel);
  return 0;
}
