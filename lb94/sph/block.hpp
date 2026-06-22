// Power-of-2 block-timestep KDK leapfrog (the individual-timestep scheme HK89/LB94 use).
// Each particle is on a rung r with dt = dt0/2^r; the system advances in integer units of
// dt_min = dt0/2^Rmax. Per sub-step: drift ALL particles to the next sync time, rebuild the
// tree, recompute forces for the ACTIVE set only (those whose step ends now), then their
// end-of-step + start-of-next half-kicks. Rung changes obey block alignment: a particle may
// only lengthen its step (lower rung) at a time that is a boundary of the new rung.
#ifndef SPH_BLOCK_HPP
#define SPH_BLOCK_HPP
#include <vector>
#include <cmath>
#include "forces.hpp"

namespace sph {

struct BlockStepper {
  double dt0 = 0.0;                 // rung-0 (max) timestep; 0 -> auto from the initial state
  int Rmax = 24;                    // deepest rung
  double C_cour = 0.25, C_acc = 0.25;
  long long ti = 0;                 // current integer time (units of dt_min)

  double dt_min() const { return dt0 / (double)(1LL << Rmax); }

  double dt_desired(const Particle& q, const SPHParams& par) const {
    double dt = dt0;
    double vsig = q.cs + 1.2 * (par.alpha * q.cs + par.beta * q.h * std::max(-q.divv, 0.0));
    if (vsig > 0) dt = std::min(dt, C_cour * q.h / vsig);
    double a = std::sqrt(q.ax * q.ax + q.ay * q.ay + q.az * q.az);
    if (a > 0) dt = std::min(dt, C_acc * std::sqrt(q.h / a));
    return dt;
  }
  int rung_for(double dt) const {
    if (dt >= dt0) return 0;
    int r = (int)std::ceil(std::log2(dt0 / dt));
    return std::min(std::max(r, 0), Rmax);
  }
  // deepest rung whose step boundary lands on integer time t (alignment floor for new rungs)
  int rung_floor(long long t) const {
    if (t == 0) return 0;
    int tz = __builtin_ctzll((unsigned long long)t);
    return std::max(0, Rmax - tz);
  }
  static void half_kick(Particle& q, double dth, bool iso) {
    q.vx += dth * q.ax; q.vy += dth * q.ay; q.vz += dth * q.az;
    if (!iso) q.u = std::max(q.u + dth * q.dudt, 1e-300);
  }

  // initial forces, rung assignment, and first half-kick
  void init(std::vector<Particle>& p, Tree& tree, const SPHParams& par) {
    int n = p.size();
    tree.build();
    compute_forces(p, tree, par);
    if (dt0 <= 0.0) {                                   // auto: dt0 = the largest desired step
      double dmax = 1e-300;
      for (int i = 0; i < n; ++i) dmax = std::max(dmax, dt_desired(p[i], par));
      dt0 = dmax;
    }
    for (int i = 0; i < n; ++i) {
      p[i].rung = rung_for(dt_desired(p[i], par));
      double dti = dt0 / (double)(1LL << p[i].rung);
      half_kick(p[i], 0.5 * dti, par.isothermal);
      p[i].ti_end = (1LL << (Rmax - p[i].rung));
    }
    ti = 0;
  }

  // advance to T_total (code units); cb(t) is called at multiples of out_dt
  template <class CB>
  void run(std::vector<Particle>& p, Tree& tree, const SPHParams& par,
           double T_total, double out_dt, CB cb) {
    int n = p.size();
    long long ti_total = (long long)std::llround(T_total / dt_min());
    long long ti_out = 0, out_step = std::max(1LL, (long long)std::llround(out_dt / dt_min()));
    std::vector<int> active;
    while (ti < ti_total) {
      int rmax_inuse = 0;
      for (int i = 0; i < n; ++i) rmax_inuse = std::max(rmax_inuse, p[i].rung);
      long long step = (1LL << (Rmax - rmax_inuse));
      long long ti_next = std::min(ti + step, ti_total);
      double dt_drift = (ti_next - ti) * dt_min();
#pragma omp parallel for schedule(static)
      for (int i = 0; i < n; ++i) {
        p[i].x += dt_drift * p[i].vx; p[i].y += dt_drift * p[i].vy; p[i].z += dt_drift * p[i].vz;
      }
      ti = ti_next;
      active.clear();
      for (int i = 0; i < n; ++i) if (p[i].ti_end <= ti) active.push_back(i);
      tree.build();
      compute_forces_active(p, tree, par, active);
      int rfloor = rung_floor(ti);
      for (int a : active) {
        Particle& q = p[a];
        double dti_old = dt0 / (double)(1LL << q.rung);
        half_kick(q, 0.5 * dti_old, par.isothermal);                // close the finished step
        int rnew = std::max(rung_for(dt_desired(q, par)), rfloor);  // accuracy + alignment
        q.rung = rnew;
        double dti_new = dt0 / (double)(1LL << rnew);
        half_kick(q, 0.5 * dti_new, par.isothermal);                // open the next step
        q.ti_end = ti + (1LL << (Rmax - rnew));
      }
      if (ti >= ti_out) { cb(ti * dt_min()); ti_out += out_step; }
    }
  }
};

}  // namespace sph
#endif
