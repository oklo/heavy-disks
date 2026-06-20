// 2D axisymmetric cylindrical (R,Z) self-gravity: del^2 Phi = 4 pi G rho.
//   (1/R) d_R(R d_R Phi) + d_ZZ Phi = 4 pi G rho
// Single quadrant R in [0,Rmax], Z in [0,Zmax]; axis (R=0) and midplane (Z=0) are
// reflecting; the outer R,Z edges use an isolated multipole (l=0,2,4) boundary.
// Solved by red-black SOR. Cell-centered, uniform dR,dZ.
#ifndef LB94_POISSON_HPP
#define LB94_POISSON_HPP

#include <vector>
#include <cmath>
#include <algorithm>

namespace lb94 {

constexpr double PI = 3.14159265358979323846;

struct Poisson {
  int NR, NZ;
  double dR, dZ, G;
  std::vector<double> R, Z;            // cell centres

  Poisson(int NR_, int NZ_, double Rmax, double Zmax, double G_ = 6.67430e-8)
      : NR(NR_), NZ(NZ_), dR(Rmax / NR_), dZ(Zmax / NZ_), G(G_) {
    R.resize(NR); Z.resize(NZ);
    for (int i = 0; i < NR; ++i) R[i] = (i + 0.5) * dR;
    for (int j = 0; j < NZ; ++j) Z[j] = (j + 0.5) * dZ;
  }
  inline int idx(int i, int j) const { return i * NZ + j; }

  static inline double P2(double x) { return 0.5 * (3 * x * x - 1); }
  static inline double P4(double x) { return 0.125 * (35 * x * x * x * x - 30 * x * x + 3); }

  // even axisymmetric multipole moments of the (Z>=0, mirrored) mass distribution
  void multipoles(const std::vector<double>& rho, double& M0, double& M2, double& M4) const {
    M0 = M2 = M4 = 0.0;
    for (int i = 0; i < NR; ++i)
      for (int j = 0; j < NZ; ++j) {
        double dm = rho[idx(i, j)] * 2.0 * PI * R[i] * dR * dZ;   // ring mass (Z>=0)
        double r = std::sqrt(R[i] * R[i] + Z[j] * Z[j]);
        double x = Z[j] / r;
        M0 += dm; M2 += dm * r * r * P2(x); M4 += dm * r * r * r * r * P4(x);
      }
    M0 *= 2; M2 *= 2; M4 *= 2;                                    // mirror Z<0
  }
  double bphi(double Rb, double Zb, double M0, double M2, double M4) const {
    double r = std::sqrt(Rb * Rb + Zb * Zb), x = Zb / r;
    return -G * (M0 / r + M2 * P2(x) / (r * r * r) + M4 * P4(x) / (r * r * r * r * r));
  }

  // Solve for Phi (size NR*NZ). Reuses Phi as the initial guess if correctly sized.
  int solve(const std::vector<double>& rho, std::vector<double>& Phi,
            int maxiter = 20000, double tol = 1e-7) const {
    double M0, M2, M4; multipoles(rho, M0, M2, M4);
    double scale = G * M0 / R[NR - 1];
    if (Phi.size() != (size_t)NR * NZ) Phi.assign((size_t)NR * NZ, -scale);
    const double omega = 2.0 / (1.0 + std::sin(PI / std::max(NR, NZ)));   // ~optimal SOR
    const double src = 4.0 * PI * G;
    int it = 0;
    for (; it < maxiter; ++it) {
      double maxres = 0.0;
      for (int color = 0; color < 2; ++color)
        for (int i = 0; i < NR; ++i) {
          double cE = (i + 1) * dR / (R[i] * dR * dR);            // east face R_{i+1/2}
          double cW = i * dR / (R[i] * dR * dR);                  // west face R_{i-1/2} (=0 at axis)
          double cN = 1.0 / (dZ * dZ), cS = cN, cD = cE + cW + cN + cS;
          for (int j = (i + color) & 1; j < NZ; j += 2) {
            double pE = (i < NR - 1) ? Phi[idx(i + 1, j)] : bphi(R[i] + dR, Z[j], M0, M2, M4);
            double pW = (i > 0)      ? Phi[idx(i - 1, j)] : Phi[idx(i, j)];     // axis reflect
            double pN = (j < NZ - 1) ? Phi[idx(i, j + 1)] : bphi(R[i], Z[j] + dZ, M0, M2, M4);
            double pS = (j > 0)      ? Phi[idx(i, j - 1)] : Phi[idx(i, j)];     // midplane mirror
            double rhs = cE * pE + cW * pW + cN * pN + cS * pS - src * rho[idx(i, j)];
            double res = rhs / cD - Phi[idx(i, j)];
            Phi[idx(i, j)] += omega * res;
            maxres = std::max(maxres, std::fabs(res));
          }
        }
      if (maxres < tol * std::fabs(scale)) break;
    }
    return it;
  }

  // Solve with given outer-boundary potentials (for nested interior grids):
  //   gR[j] = Phi in the ghost just outside R=Rmax;  gZ[i] = ghost just outside Z=Zmax.
  // Axis (R=0) and midplane (Z=0) stay reflecting.
  int solve_dirichlet(const std::vector<double>& rho, std::vector<double>& Phi,
                      const std::vector<double>& gR, const std::vector<double>& gZ,
                      int maxiter = 20000, double tol = 1e-7) const {
    double scale = 0.0;
    for (double v : gR) scale = std::max(scale, std::fabs(v));
    if (scale == 0) scale = 1.0;
    if (Phi.size() != (size_t)NR * NZ) Phi.assign((size_t)NR * NZ, -scale);
    const double omega = 2.0 / (1.0 + std::sin(PI / std::max(NR, NZ)));
    const double src = 4.0 * PI * G;
    int it = 0;
    for (; it < maxiter; ++it) {
      double maxres = 0.0;
      for (int color = 0; color < 2; ++color)
        for (int i = 0; i < NR; ++i) {
          double cE = (i + 1) * dR / (R[i] * dR * dR), cW = i * dR / (R[i] * dR * dR);
          double cN = 1.0 / (dZ * dZ), cS = cN, cD = cE + cW + cN + cS;
          for (int j = (i + color) & 1; j < NZ; j += 2) {
            double pE = (i < NR - 1) ? Phi[idx(i + 1, j)] : gR[j];
            double pW = (i > 0)      ? Phi[idx(i - 1, j)] : Phi[idx(i, j)];
            double pN = (j < NZ - 1) ? Phi[idx(i, j + 1)] : gZ[i];
            double pS = (j > 0)      ? Phi[idx(i, j - 1)] : Phi[idx(i, j)];
            double rhs = cE * pE + cW * pW + cN * pN + cS * pS - src * rho[idx(i, j)];
            double res = rhs / cD - Phi[idx(i, j)];
            Phi[idx(i, j)] += omega * res;
            maxres = std::max(maxres, std::fabs(res));
          }
        }
      if (maxres < tol * std::fabs(scale)) break;
    }
    return it;
  }
};

}  // namespace lb94
#endif
