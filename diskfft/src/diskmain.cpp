// Driver: LKA98 heavy disk with the from-scratch van-Leer hydro + Kalnajs FFT gravity.
//   ./diskdisk Nu Nphi eps tlim selfgrav(0/1) amp
// Outputs t, |C0| (mass), |C2|, arg(C2) -> growth rate & pattern speed of the m=2 spiral.
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "hydro.hpp"
#include "kalnajs.hpp"

using namespace diskfft;

int main(int argc, char **argv) {
  int Nu = (argc > 1) ? atoi(argv[1]) : 128;
  int Nphi = (argc > 2) ? atoi(argv[2]) : 128;
  double eps = (argc > 3) ? atof(argv[3]) : 0.02;
  double tlim = (argc > 4) ? atof(argv[4]) : 8.0;
  int selfgrav = (argc > 5) ? atoi(argv[5]) : 1;
  double amp = (argc > 6) ? atof(argv[6]) : 1e-3;
  int nthreads = (argc > 7) ? atoi(argv[7]) : 1;
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif

  double Rin = 0.25, Rout = 1.0, R0 = 0.45, w2 = 0.05, S0 = 0.372;
  double Kpoly = 0.25, gamma_gas = 2.0, gmstar = 0.6, G = 1.0;

  DiskHydro hyd(Nu, Nphi, Rin, Rout, Kpoly, gamma_gas, gmstar);
  KalnajsSolver grav(Nu, Nphi, Rin, Rout, eps, G, nthreads);

  auto Sigma0 = [&](double R) { return S0 * std::exp(-(R - R0) * (R - R0) / w2); };

  // axisymmetric disk radial acceleration (for the self-gravitating equilibrium)
  std::vector<double> gR_eq(Nu * Nphi, 0.0), gp_eq, Phi_eq, Sig_axi(Nu * Nphi);
  for (int i = 0; i < Nu; ++i)
    for (int j = 0; j < Nphi; ++j) Sig_axi[i * Nphi + j] = Sigma0(hyd.R[i]);
  std::vector<double> gRdisk(Nu, 0.0);
  if (selfgrav) {
    grav.potential(Sig_axi, Phi_eq);
    grav.forces(Phi_eq, gR_eq, gp_eq);
    for (int i = 0; i < Nu; ++i) {           // phi-average (already axisymmetric)
      double s = 0; for (int j = 0; j < Nphi; ++j) s += gR_eq[i * Nphi + j];
      gRdisk[i] = s / Nphi;
    }
  }

  // initial condition: equilibrium + m=2 density seed
  for (int i = 0; i < Nu; ++i) {
    double R = hyd.R[i];
    double dSig = Sigma0(R) * (-2.0 * (R - R0) / w2);
    double dPdR = Kpoly * gamma_gas * std::pow(Sigma0(R), gamma_gas - 1.0) * dSig;
    double v2 = gmstar / R + (R / Sigma0(R)) * dPdR - R * gRdisk[i];
    double vphi = (v2 > 0) ? std::sqrt(v2) : 0.0;
    for (int j = 0; j < Nphi; ++j) {
      double phi = (j + 0.5) * hyd.dphi;
      hyd.Sig[i * Nphi + j] = Sigma0(R) * (1.0 + amp * std::cos(2.0 * phi));
      hyd.vR[i * Nphi + j] = 0.0;
      hyd.vp[i * Nphi + j] = vphi;
    }
  }

  // total disk mass for normalisation
  double mD = 0.0;
  for (int i = 0; i < Nu; ++i) {
    double dA = hyd.R[i] * (hyd.Rf[i + 1] - hyd.Rf[i]) * hyd.dphi;
    for (int j = 0; j < Nphi; ++j) mD += hyd.Sig[i * Nphi + j] * dA;
  }

  auto Cm = [&](int m, double &mag, double &phase) {
    double re = 0, im = 0;
    for (int i = 0; i < Nu; ++i) {
      double dA = hyd.R[i] * (hyd.Rf[i + 1] - hyd.Rf[i]) * hyd.dphi;
      for (int j = 0; j < Nphi; ++j) {
        double phi = (j + 0.5) * hyd.dphi;
        re += hyd.Sig[i * Nphi + j] * std::cos(m * phi) * dA;
        im -= hyd.Sig[i * Nphi + j] * std::sin(m * phi) * dA;
      }
    }
    mag = std::sqrt(re * re + im * im) / mD;
    phase = std::atan2(im, re);
  };

  FILE *fp = fopen("diskfft_modes.dat", "w");
  fprintf(fp, "# t  |C0|  |C2|  arg(C2)\n");
  double t = 0.0, dtout = 0.25, tnext = 0.0;
  std::vector<double> Phi, gR, gp;
  double c0, c2, ph2, c0p, ph0;
  int step = 0;
  while (t < tlim) {
    if (t >= tnext - 1e-12) {
      Cm(0, c0, ph0); Cm(2, c2, ph2);
      printf("t=%7.3f  |C0|=%.5f  |C2|=%.6e\n", t, c0, c2);
      fprintf(fp, "%.5f %.6e %.6e %.6e\n", t, c0, c2, ph2);
      fflush(fp);
      tnext += dtout;
    }
    double dt = hyd.timestep(0.3);
    if (t + dt > tlim) dt = tlim - t;
    if (selfgrav) {
      grav.potential(hyd.Sig, Phi);
      grav.forces(Phi, gR, gp);
      hyd.step(dt, &gR, &gp);
    } else {
      hyd.step(dt);
    }
    t += dt; ++step;
  }
  Cm(2, c2, ph2);
  printf("final t=%.3f  |C2|=%.6e  steps=%d\n", t, c2, step);
  fclose(fp);
  return 0;
}
