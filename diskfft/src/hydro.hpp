// 2D razor-thin polar (R,phi) hydro: operator-split 2nd-order van-Leer (MUSCL)
// advection + barotropic polytropic EOS P=K Sigma^gamma + source terms.
// Transport restructured as race-free flux-then-gather; OpenMP-parallel throughout.
#ifndef DISKFFT_HYDRO_HPP
#define DISKFFT_HYDRO_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace diskfft {

inline double vanleer(double a, double b) {
  return (a * b > 0.0) ? 2.0 * a * b / (a + b) : 0.0;
}

struct DiskHydro {
  int Nu, Nphi;
  double umin, umax, du, dphi, Kpoly, gamma_gas, gmstar;
  std::vector<double> R, Rf, Vol;       // centres (Nu), faces (Nu+1), cell volume (Nu)
  std::vector<double> Sig, vR, vp;      // primitive state (Nu*Nphi)
  // pre-allocated per-step work buffers (avoid malloc churn / improve scaling)
  std::vector<double> wA, wB, mR, mp, Fm, FmR, FmP;

  DiskHydro(int Nu_, int Nphi_, double Rin, double Rout,
            double K_, double g_, double GMstar_)
      : Nu(Nu_), Nphi(Nphi_), Kpoly(K_), gamma_gas(g_), gmstar(GMstar_) {
    umin = std::log(Rin); umax = std::log(Rout);
    du = (umax - umin) / Nu; dphi = 2.0 * M_PI / Nphi;
    R.resize(Nu); Rf.resize(Nu + 1); Vol.resize(Nu);
    for (int i = 0; i <= Nu; ++i) Rf[i] = std::exp(umin + i * du);
    for (int i = 0; i < Nu; ++i) {
      R[i] = std::exp(umin + (i + 0.5) * du);
      Vol[i] = R[i] * (Rf[i + 1] - Rf[i]) * dphi;
    }
    Sig.assign(Nu * Nphi, 0.0); vR.assign(Nu * Nphi, 0.0); vp.assign(Nu * Nphi, 0.0);
    wA.assign(Nu * Nphi, 0.0); wB.assign(Nu * Nphi, 0.0);
    mR.assign(Nu * Nphi, 0.0); mp.assign(Nu * Nphi, 0.0);
    Fm.assign((Nu + 1) * Nphi, 0.0); FmR.assign((Nu + 1) * Nphi, 0.0); FmP.assign((Nu + 1) * Nphi, 0.0);
  }
  inline int idx(int i, int j) const { return i * Nphi + j; }
  inline double cs2(double s) const { return gamma_gas * Kpoly * std::pow(s, gamma_gas - 1.0); }
  inline double pres(double s) const { return Kpoly * std::pow(s, gamma_gas); }

  double timestep(double cfl) const {
    double dtinv = 1e-30;
    #pragma omp parallel for reduction(max:dtinv) schedule(static)
    for (int i = 0; i < Nu; ++i) {
      double dR = Rf[i + 1] - Rf[i];
      for (int j = 0; j < Nphi; ++j) {
        double cs = std::sqrt(cs2(std::max(Sig[idx(i, j)], 1e-12)));
        double sR = (std::fabs(vR[idx(i, j)]) + cs) / dR;
        double sp = (std::fabs(vp[idx(i, j)]) + cs) / (R[i] * dphi);
        if (sR + sp > dtinv) dtinv = sR + sp;
      }
    }
    return cfl / dtinv;
  }

  void step(double dt, const std::vector<double> *gR_sg = nullptr,
            const std::vector<double> *gp_sg = nullptr) {
    source(dt, gR_sg, gp_sg);
    transportR(dt);
    transportPhi(dt);
    #pragma omp parallel for schedule(static)
    for (int k = 0; k < Nu * Nphi; ++k) Sig[k] = std::max(Sig[k], 1e-10);
  }

  // --- source step: pressure grad, gravity, geometric (centrifugal/coriolis) ---
  void source(double dt, const std::vector<double> *gR_sg, const std::vector<double> *gp_sg) {
    std::vector<double> &nvR = wA, &nvp = wB;
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < Nu; ++i) {
      int im = std::max(i - 1, 0), ip = std::min(i + 1, Nu - 1);
      double dRc = R[ip] - R[im];
      for (int j = 0; j < Nphi; ++j) {
        int c = idx(i, j);
        int jm = idx(i, (j - 1 + Nphi) % Nphi), jp = idx(i, (j + 1) % Nphi);
        double s = Sig[c];
        double dPdR = (pres(Sig[idx(ip, j)]) - pres(Sig[idx(im, j)])) / dRc;
        double dPdphi = (pres(Sig[jp]) - pres(Sig[jm])) / (2.0 * dphi);
        double gR = -gmstar / (R[i] * R[i]);
        double gp = 0.0;
        if (gR_sg) gR += (*gR_sg)[c];
        if (gp_sg) gp += (*gp_sg)[c];
        double aR = -dPdR / s + gR + vp[c] * vp[c] / R[i];
        double ap = -(1.0 / (s * R[i])) * dPdphi + gp - vR[c] * vp[c] / R[i];
        nvR[c] = vR[c] + dt * aR;
        nvp[c] = vp[c] + dt * ap;
      }
    }
    vR.swap(nvR); vp.swap(nvp);
  }

  // --- radial transport: flux at faces (race-free), then gather into cells ---
  void transportR(double dt) {
    #pragma omp parallel for schedule(static)
    for (int k = 0; k < Nu * Nphi; ++k) { mR[k] = Sig[k] * vR[k]; mp[k] = Sig[k] * vp[k]; }
    // reflecting boundary faces f=0 and f=Nu -> zero flux (buffer is shared w/ transportPhi)
    for (int j = 0; j < Nphi; ++j) {
      Fm[j] = FmR[j] = FmP[j] = 0.0;
      Fm[Nu * Nphi + j] = FmR[Nu * Nphi + j] = FmP[Nu * Nphi + j] = 0.0;
    }
    auto slopeR = [&](const std::vector<double> &a, int i, int j) {
      if (i <= 0 || i >= Nu - 1) return 0.0;
      return vanleer(a[idx(i + 1, j)] - a[idx(i, j)], a[idx(i, j)] - a[idx(i - 1, j)]);
    };
    #pragma omp parallel for schedule(static)
    for (int f = 1; f < Nu; ++f) {                 // interior radial faces (between f-1, f)
      double Af = Rf[f] * dphi;
      double dRl = Rf[f] - Rf[f - 1], dRr = Rf[f + 1] - Rf[f];
      for (int j = 0; j < Nphi; ++j) {
        double vf = 0.5 * (vR[idx(f - 1, j)] + vR[idx(f, j)]);
        double sf, vRf, vpf;
        if (vf > 0) {
          sf  = Sig[idx(f - 1, j)] + 0.5 * slopeR(Sig, f - 1, j) * (1.0 - vf * dt / dRl);
          vRf = vR[idx(f - 1, j)]  + 0.5 * slopeR(vR, f - 1, j) * (1.0 - vf * dt / dRl);
          vpf = vp[idx(f - 1, j)]  + 0.5 * slopeR(vp, f - 1, j) * (1.0 - vf * dt / dRl);
        } else {
          sf  = Sig[idx(f, j)] - 0.5 * slopeR(Sig, f, j) * (1.0 + vf * dt / dRr);
          vRf = vR[idx(f, j)]  - 0.5 * slopeR(vR, f, j) * (1.0 + vf * dt / dRr);
          vpf = vp[idx(f, j)]  - 0.5 * slopeR(vp, f, j) * (1.0 + vf * dt / dRr);
        }
        double fm = vf * sf * Af;
        Fm[f * Nphi + j] = fm; FmR[f * Nphi + j] = fm * vRf; FmP[f * Nphi + j] = fm * vpf;
      }
    }
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < Nu; ++i)
      for (int j = 0; j < Nphi; ++j) {
        int c = idx(i, j);
        double ns = Sig[c] + dt * (Fm[i * Nphi + j] - Fm[(i + 1) * Nphi + j]) / Vol[i];
        double nmR = mR[c] + dt * (FmR[i * Nphi + j] - FmR[(i + 1) * Nphi + j]) / Vol[i];
        double nmP = mp[c] + dt * (FmP[i * Nphi + j] - FmP[(i + 1) * Nphi + j]) / Vol[i];
        Sig[c] = ns; vR[c] = nmR / std::max(ns, 1e-10); vp[c] = nmP / std::max(ns, 1e-10);
      }
  }

  // --- azimuthal transport (periodic): flux at faces, gather into cells ---
  void transportPhi(double dt) {
    #pragma omp parallel for schedule(static)
    for (int k = 0; k < Nu * Nphi; ++k) { mR[k] = Sig[k] * vR[k]; mp[k] = Sig[k] * vp[k]; }
    auto slopeP = [&](const std::vector<double> &a, int i, int j) {
      int jm = (j - 1 + Nphi) % Nphi, jp = (j + 1) % Nphi;
      return vanleer(a[idx(i, jp)] - a[idx(i, j)], a[idx(i, j)] - a[idx(i, jm)]);
    };
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < Nu; ++i) {
      double dl = R[i] * dphi, Aface = Rf[i + 1] - Rf[i];
      for (int j = 0; j < Nphi; ++j) {             // face f at j-1/2 (between j-1, j)
        int jm = (j - 1 + Nphi) % Nphi;
        double vf = 0.5 * (vp[idx(i, jm)] + vp[idx(i, j)]);
        double sf, vRf, vpf;
        if (vf > 0) {
          sf  = Sig[idx(i, jm)] + 0.5 * slopeP(Sig, i, jm) * (1.0 - vf * dt / dl);
          vRf = vR[idx(i, jm)]  + 0.5 * slopeP(vR, i, jm) * (1.0 - vf * dt / dl);
          vpf = vp[idx(i, jm)]  + 0.5 * slopeP(vp, i, jm) * (1.0 - vf * dt / dl);
        } else {
          sf  = Sig[idx(i, j)] - 0.5 * slopeP(Sig, i, j) * (1.0 + vf * dt / dl);
          vRf = vR[idx(i, j)]  - 0.5 * slopeP(vR, i, j) * (1.0 + vf * dt / dl);
          vpf = vp[idx(i, j)]  - 0.5 * slopeP(vp, i, j) * (1.0 + vf * dt / dl);
        }
        double fm = vf * sf * Aface;
        Fm[idx(i, j)] = fm; FmR[idx(i, j)] = fm * vRf; FmP[idx(i, j)] = fm * vpf;
      }
    }
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < Nu; ++i)
      for (int j = 0; j < Nphi; ++j) {
        int c = idx(i, j); int jp = idx(i, (j + 1) % Nphi);
        double ns = Sig[c] + dt * (Fm[c] - Fm[jp]) / Vol[i];
        double nmR = mR[c] + dt * (FmR[c] - FmR[jp]) / Vol[i];
        double nmP = mp[c] + dt * (FmP[c] - FmP[jp]) / Vol[i];
        Sig[c] = ns; vR[c] = nmR / std::max(ns, 1e-10); vp[c] = nmP / std::max(ns, 1e-10);
      }
  }
};

}  // namespace diskfft
#endif
