// Free-fall collapse of a cold uniform sphere vs the analytic homologous solution.
//   shell:  r = r0 cos^2(b),   t = sqrt(3/(8 pi G rho0)) (b + 1/2 sin 2b)
//   central density:  rho_c(t)/rho0 = sec^6(b),   t_ff = sqrt(3 pi/(32 G rho0)).
#include "hydro.hpp"
#include <cstdio>
#include <cmath>

using namespace lb94;

// invert t(b) for b by Newton
static double beta_of_t(double t, double K) {
  double b = 1.0;
  for (int k = 0; k < 60; ++k) {
    double f = K * (b + 0.5 * std::sin(2 * b)) - t;
    double df = K * (1.0 + std::cos(2 * b));
    b -= f / df;
    b = std::min(std::max(b, 0.0), 0.5 * PI - 1e-6);
  }
  return b;
}

int main() {
  const double G = 1.0, rho0 = 1.0, a = 1.0, Rmax = 2.0, Zmax = 2.0;
  const int N = 80;
  const double t_ff = std::sqrt(3.0 * PI / (32.0 * G * rho0));
  const double K = std::sqrt(3.0 / (8.0 * PI * G * rho0));

  Hydro h(N, N, Rmax, Zmax, G, /*cs2=*/1e-3);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      h.rho[h.idx(i, j)] =
          (std::sqrt(h.R[i] * h.R[i] + h.Z[j] * h.Z[j]) < a) ? rho0 : h.floor;

  auto mass = [&]() {
    double m = 0;
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j) m += h.rho[h.idx(i, j)] * 2 * PI * h.R[i] * h.dR * h.dZ;
    return 2 * m;                                              // mirror Z<0
  };
  double m0 = mass();
  printf("t_ff = %.4f,  initial mass = %.4f (analytic %.4f)\n",
         t_ff, m0, 4.0 / 3.0 * PI * rho0 * a * a * a);
  printf("\n t/t_ff   rho_c(sim)   rho_c(theory)   ratio    mass/mass0\n");

  double t = 0;
  double checkpoints[] = {0.3, 0.5, 0.7};
  int nc = 0, nstep = 0;
  while (t < 0.71 * t_ff && nc < 3 && nstep < 30000) {
    h.poisson.solve(h.rho, h.Phi, 3000, 1e-5);                // warm-started, looser tol
    double dt = std::min(h.timestep(0.3), 0.005 * t_ff);
    if (t + dt > checkpoints[nc] * t_ff) dt = checkpoints[nc] * t_ff - t;
    h.step(dt);
    t += dt; ++nstep;
    if (nstep % 200 == 0) fprintf(stderr, "  step %d  t/t_ff=%.3f  rho_c=%.2f\n",
                                  nstep, t / t_ff, h.rho[h.idx(0, 0)]);
    if (t >= checkpoints[nc] * t_ff - 1e-12) {
      double b = beta_of_t(t, K);
      double th = rho0 * std::pow(1.0 / std::cos(b), 6);
      double sim = h.rho[h.idx(0, 0)];
      printf(" %.3f   %9.3f   %11.3f   %.3f    %.5f\n",
             t / t_ff, sim, th, sim / th, mass() / m0);
      fflush(stdout);
      nc++;
    }
  }
  return 0;
}
