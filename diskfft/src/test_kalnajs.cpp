// Validate the Kalnajs FFT self-gravity solver against a direct O(N^2) convolution.
#include <cstdio>
#include <cmath>
#include <vector>
#include "kalnajs.hpp"

using namespace diskfft;

int main() {
  int Nu = 32, Nphi = 32;
  double Rin = 0.25, Rout = 1.0, eps = 0.02;
  KalnajsSolver solver(Nu, Nphi, Rin, Rout, eps);

  // test surface density: Gaussian ring (R0=0.45,w2=0.05) + small m=2 perturbation
  std::vector<double> Sigma(Nu * Nphi);
  for (int i = 0; i < Nu; ++i) {
    double R = solver.R[i];
    double S0 = 0.372 * std::exp(-(R - 0.45) * (R - 0.45) / 0.05);
    for (int j = 0; j < Nphi; ++j) {
      double phi = (j + 0.5) * solver.dphi;
      Sigma[i * Nphi + j] = S0 * (1.0 + 0.3 * std::cos(2.0 * phi));
    }
  }

  std::vector<double> Phi_fft, Phi_dir;
  solver.potential(Sigma, Phi_fft);
  solver.potential_direct(Sigma, Phi_dir);

  double maxabs = 0.0, maxdiff = 0.0;
  for (int k = 0; k < Nu * Nphi; ++k) {
    maxabs = std::max(maxabs, std::fabs(Phi_dir[k]));
    maxdiff = std::max(maxdiff, std::fabs(Phi_fft[k] - Phi_dir[k]));
  }
  std::printf("Kalnajs FFT vs direct convolution (Nu=%d Nphi=%d eps=%.3f):\n", Nu, Nphi, eps);
  std::printf("  max|Phi| = %.6e   max|FFT-direct| = %.3e   rel = %.3e\n",
              maxabs, maxdiff, maxdiff / maxabs);
  // sample a few values
  std::printf("  Phi[mid,0]   fft=%.6e dir=%.6e\n",
              Phi_fft[(Nu / 2) * Nphi + 0], Phi_dir[(Nu / 2) * Nphi + 0]);
  std::printf("  Phi[mid,Nphi/2] fft=%.6e dir=%.6e\n",
              Phi_fft[(Nu / 2) * Nphi + Nphi / 2], Phi_dir[(Nu / 2) * Nphi + Nphi / 2]);
  return (maxdiff / maxabs < 1e-10) ? 0 : 1;
}
