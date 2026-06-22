// 3D cubic-spline SPH kernel (Monaghan & Lattanzio 1985; Hernquist & Katz 1989, eq. 2.8).
// W(r,h) = (1/pi h^3) w(q),  q = r/h,  with compact support r < 2h:
//   w(q) = 1 - 3/2 q^2 + 3/4 q^3      0 <= q < 1
//   w(q) = 1/4 (2 - q)^3              1 <= q < 2
//   w(q) = 0                          q >= 2
// Also the spline-softened gravity (potential and force) consistent with this kernel,
// from the integral of the enclosed kernel mass (HK89 use the same spline for softening).
#ifndef SPH_KERNEL_HPP
#define SPH_KERNEL_HPP
#include <cmath>

namespace sph {

constexpr double PI = 3.14159265358979323846;

// W(r,h)
inline double W(double r, double h) {
  double q = r / h, c = 1.0 / (PI * h * h * h);
  if (q < 1.0)      return c * (1.0 - 1.5 * q * q + 0.75 * q * q * q);
  else if (q < 2.0) { double a = 2.0 - q; return c * 0.25 * a * a * a; }
  return 0.0;
}

// scalar gradient factor gW such that grad_i W(|r_i-r_j|,h) = gW * (r_i - r_j).
// gW = (1/r) dW/dr.  dW/dr = (1/pi h^4) w'(q):
//   w'(q) = -3q + 9/4 q^2     0<=q<1 ;   w'(q) = -3/4 (2-q)^2   1<=q<2
inline double gradW(double r, double h) {
  if (r <= 0.0) return 0.0;
  double q = r / h, c = 1.0 / (PI * h * h * h * h);
  double dwdr;
  if (q < 1.0)      dwdr = c * (-3.0 * q + 2.25 * q * q);
  else if (q < 2.0) { double a = 2.0 - q; dwdr = c * (-0.75 * a * a); }
  else return 0.0;
  return dwdr / r;     // multiply by the vector (r_i - r_j) to get grad_i W
}

// dW/dh at fixed r (needed for the "grad-h" terms if used; provided for completeness)
inline double dWdh(double r, double h) {
  double q = r / h;
  if (q >= 2.0) return 0.0;
  double c = 1.0 / (PI * h * h * h * h);
  // dW/dh = -(3/h) W - (q/h) dW/dq * (1/pi h^3)... computed directly:
  double w, dwdq;
  if (q < 1.0) { w = 1.0 - 1.5 * q * q + 0.75 * q * q * q; dwdq = -3.0 * q + 2.25 * q * q; }
  else { double a = 2.0 - q; w = 0.25 * a * a * a; dwdq = -0.75 * a * a; }
  return -c * (3.0 * w + q * dwdq);
}

// --- spline-softened gravity (softening length eps plays the role of h) ---
// Force law |g|(r) = G M f(r,eps)/r^2 with f the fraction of the kernel mass within r;
// returns the "g/r" factor so that the acceleration is gfac * (r_i - r_j) (points from i to j
// is handled by the caller's sign).  Newtonian for r >= 2 eps.
inline double grav_factor(double r, double eps) {
  if (r >= 2.0 * eps) return 1.0 / (r * r * r);          // Newtonian: g = -GM r_hat / r^2
  double q = r / eps, q2 = q * q, q3 = q2 * q;
  double fr;   // enclosed-mass fraction M(<r)/M for the cubic spline (Hernquist & Katz 1989)
  if (q < 1.0)
    fr = (4.0 / 3.0) * q3 - 1.2 * q3 * q2 + 0.5 * q3 * q3;          // (4/3)q^3 -6/5 q^5 +1/2 q^6
  else
    fr = (8.0 / 3.0) * q3 - 3.0 * q3 * q + 1.2 * q3 * q2 - (1.0 / 6.0) * q3 * q3 - 1.0 / 15.0;
  return fr / (eps * eps * eps * q3);   // = [M(<r)/M] / r^3, so g = -GM * gfac * r_vec
}

}  // namespace sph
#endif
