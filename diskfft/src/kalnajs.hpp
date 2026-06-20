// Kalnajs (1971) logarithmic-spiral FFT solver for razor-thin disk self-gravity.
// FFTW real-to-complex implementation (optionally multi-threaded).
//
// In u = ln R the softened razor-thin potential is a convolution:
//   V(u,phi) = -G du dphi  sum_{u',phi'} S(u',phi') K(u-u', phi-phi'),
//   S = R^{3/2} Sigma,  V = R^{1/2} Phi,
//   K(du,dphi) = 1 / sqrt(2 cosh(du) - 2 cos(dphi) + eps).
// Circular in phi (period 2pi), linear in u (zero-padded to M=2*Nu).  Real fields,
// so r2c/c2r transforms are used (half the work of a full complex FFT).
#ifndef DISKFFT_KALNAJS_HPP
#define DISKFFT_KALNAJS_HPP

#include <vector>
#include <cmath>
#include <fftw3.h>

namespace diskfft {

class KalnajsSolver {
 public:
  int Nu, Nphi, M, Nc;      // M=2*Nu padded radial size; Nc=Nphi/2+1 (r2c last dim)
  double umin, umax, du, dphi, eps, G;
  std::vector<double> R, sqrtR, R32;

  KalnajsSolver(int Nu_, int Nphi_, double Rin, double Rout, double eps_,
                double G_ = 1.0, int nthreads = 1)
      : Nu(Nu_), Nphi(Nphi_), M(2 * Nu_), Nc(Nphi_ / 2 + 1), eps(eps_), G(G_) {
    umin = std::log(Rin); umax = std::log(Rout);
    du = (umax - umin) / Nu; dphi = 2.0 * M_PI / Nphi;
    R.resize(Nu); sqrtR.resize(Nu); R32.resize(Nu);
    for (int i = 0; i < Nu; ++i) {
      R[i] = std::exp(umin + (i + 0.5) * du);
      sqrtR[i] = std::sqrt(R[i]); R32[i] = R[i] * sqrtR[i];
    }
    if (nthreads > 1) { fftw_init_threads(); fftw_plan_with_nthreads(nthreads); }
    rbuf = fftw_alloc_real(M * Nphi);
    cbuf = fftw_alloc_complex(M * Nc);
    Khat = fftw_alloc_complex(M * Nc);
    plan_fwd = fftw_plan_dft_r2c_2d(M, Nphi, rbuf, cbuf, FFTW_MEASURE);
    plan_inv = fftw_plan_dft_c2r_2d(M, Nphi, cbuf, rbuf, FFTW_MEASURE);
    build_kernel();
  }
  ~KalnajsSolver() {
    fftw_destroy_plan(plan_fwd); fftw_destroy_plan(plan_inv);
    fftw_free(rbuf); fftw_free(cbuf); fftw_free(Khat);
  }

  void potential(const std::vector<double> &Sigma, std::vector<double> &Phi) {
    for (int k = 0; k < M * Nphi; ++k) rbuf[k] = 0.0;
    for (int i = 0; i < Nu; ++i)
      for (int j = 0; j < Nphi; ++j) rbuf[i * Nphi + j] = R32[i] * Sigma[i * Nphi + j];
    fftw_execute(plan_fwd);                       // rbuf -> cbuf = Shat
    for (int k = 0; k < M * Nc; ++k) {            // cbuf *= Khat (complex)
      double a = cbuf[k][0], b = cbuf[k][1], c = Khat[k][0], d = Khat[k][1];
      cbuf[k][0] = a * c - b * d; cbuf[k][1] = a * d + b * c;
    }
    fftw_execute(plan_inv);                        // cbuf -> rbuf = V (unnormalised)
    double pref = -G * du * dphi / (static_cast<double>(M) * Nphi);
    Phi.assign(Nu * Nphi, 0.0);
    for (int i = 0; i < Nu; ++i)
      for (int j = 0; j < Nphi; ++j)
        Phi[i * Nphi + j] = pref * rbuf[i * Nphi + j] / sqrtR[i];
  }

  void forces(const std::vector<double> &Phi, std::vector<double> &gR,
              std::vector<double> &gphi) {
    gR.assign(Nu * Nphi, 0.0); gphi.assign(Nu * Nphi, 0.0);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < Nu; ++i)
      for (int j = 0; j < Nphi; ++j) {
        int jm = (j - 1 + Nphi) % Nphi, jp = (j + 1) % Nphi;
        double dPhidphi = (Phi[i * Nphi + jp] - Phi[i * Nphi + jm]) / (2.0 * dphi);
        double dPhidu;
        if (i == 0) dPhidu = (Phi[(i + 1) * Nphi + j] - Phi[i * Nphi + j]) / du;
        else if (i == Nu - 1) dPhidu = (Phi[i * Nphi + j] - Phi[(i - 1) * Nphi + j]) / du;
        else dPhidu = (Phi[(i + 1) * Nphi + j] - Phi[(i - 1) * Nphi + j]) / (2.0 * du);
        gR[i * Nphi + j] = -dPhidu / R[i];
        gphi[i * Nphi + j] = -dPhidphi / R[i];
      }
  }

  // direct O(N^2) convolution (validation reference)
  void potential_direct(const std::vector<double> &Sigma, std::vector<double> &Phi) {
    Phi.assign(Nu * Nphi, 0.0);
    double pref = -G * du * dphi;
    for (int i = 0; i < Nu; ++i)
      for (int j = 0; j < Nphi; ++j) {
        double V = 0.0;
        for (int ip = 0; ip < Nu; ++ip) {
          double Du = (i - ip) * du;
          for (int jp = 0; jp < Nphi; ++jp) {
            double Dphi = (j - jp) * dphi;
            V += R32[ip] * Sigma[ip * Nphi + jp]
                 / std::sqrt(2.0 * std::cosh(Du) - 2.0 * std::cos(Dphi) + eps);
          }
        }
        Phi[i * Nphi + j] = pref * V / sqrtR[i];
      }
  }

 private:
  double *rbuf; fftw_complex *cbuf, *Khat;
  fftw_plan plan_fwd, plan_inv;

  void build_kernel() {
    for (int p = 0; p < M; ++p) {
      int ioff = (p <= M / 2) ? p : p - M;
      double Du = ioff * du;
      for (int q = 0; q < Nphi; ++q) {
        double Dphi = q * dphi;
        rbuf[p * Nphi + q] = 1.0 / std::sqrt(2.0 * std::cosh(Du) - 2.0 * std::cos(Dphi) + eps);
      }
    }
    fftw_execute(plan_fwd);                        // rbuf -> cbuf = Khat
    for (int k = 0; k < M * Nc; ++k) { Khat[k][0] = cbuf[k][0]; Khat[k][1] = cbuf[k][1]; }
  }
};

}  // namespace diskfft
#endif
