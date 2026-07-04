// Quiet-start null experiment: evolve the ring-stratified realization of the basic-state
// family and watch the azimuthal moments. Linear-theory prediction (LINEAR_PREDICTIONS.md):
// |c_m| FLATLINES (no growth above gamma ~ 0.04/T_unit). The Poisson-IC runs grew m=1 at
// gamma ~ 3-6/T_unit -- if that were a real mode of the disk, it must appear here too,
// growing from the (tree-noise-level) seed.
//
// Setup differences vs treesph.cpp: the star is COLLISIONLESS (gravity only, no SPH
// coupling); per-ring v_phi is REBALANCED against the code's own measured radial forces
// before release; no absorbing boundary (the family taper is self-contained); theta=0.4;
// fine-cadence output of c_m amplitude AND phase, plus per-annulus m=1.
#include "block.hpp"
#include <cstdio>
#include <cmath>
#include <complex>
#include <chrono>
#include <map>
#include <string>

using namespace sph;

int main(int argc, char** argv) {
  const char* icfile = (argc > 1) ? argv[1] : "lb94/sph/quiet_ic.dat";
  const double t_end = (argc > 2) ? std::atof(argv[2]) : 3.0;
  const int m_seed = (argc > 3) ? std::atoi(argv[3]) : 0;       // 0 = unseeded control
  const double eps_seed = (argc > 4) ? std::atof(argv[4]) : 0.0;

  std::vector<Particle> p;
  FILE* fp = fopen(icfile, "r");
  if (!fp) { fprintf(stderr, "cannot open %s\n", icfile); return 1; }
  char line[512];
  while (fgets(line, sizeof(line), fp)) {
    if (line[0] == '#') continue;
    Particle q{}; int heavy;
    if (sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %d", &q.x, &q.y, &q.z,
               &q.vx, &q.vy, &q.vz, &q.m, &q.cs2_iso, &heavy) == 9) {
      q.heavy = heavy; q.gas = !heavy; if (!q.gas) q.h = 0.11; else q.h = 0.05;
      p.push_back(q);
    }
  }
  fclose(fp);
  int n = p.size(), star_i = -1;
  double Mdisk = 0;
  for (int i = 0; i < n; ++i) { if (p[i].heavy) star_i = i; else Mdisk += p[i].m; }

  SPHParams par; par.isothermal = true; par.selfgrav = true; par.N_neigh = 50;
  Tree tree; tree.P = &p; tree.theta = 0.4;                 // tighter: lower force-noise seed
  adapt_h(p, tree, par, 30);
  compute_forces(p, tree, par);

  // ---- per-ring v_phi rebalance against the measured forces (2 iterations) ----
  // rings are identified by their exact cylindrical radius (identical per ring)
  std::map<long long, std::vector<int>> rings;
  for (int i = 0; i < n; ++i) {
    if (p[i].heavy) continue;
    double R = std::hypot(p[i].x, p[i].y);
    rings[(long long)std::llround(R * 1e12)].push_back(i);
  }
  for (int it = 0; it < 2; ++it) {
    int nneg = 0;
    for (auto& [key, mem] : rings) {
      double aR = 0, R = std::hypot(p[mem[0]].x, p[mem[0]].y);
      for (int i : mem) aR += (p[i].ax * p[i].x + p[i].ay * p[i].y) / R;
      aR /= mem.size();
      double v2 = -aR * R;
      if (v2 <= 0) { ++nneg; v2 = 0; }
      double vph = std::sqrt(v2);
      for (int i : mem) {
        double phi = std::atan2(p[i].y, p[i].x);
        p[i].vx = -vph * std::sin(phi); p[i].vy = vph * std::cos(phi); p[i].vz = 0;
      }
    }
    compute_forces(p, tree, par);
    if (nneg) printf("  rebalance iter %d: %d rings with outward net force (v set 0)\n", it, nneg);
  }

  // ---- diagnostics ----
  auto modes = [&](double cm[5], double ph[5]) {
    double hx = p[star_i].x, hy = p[star_i].y;
    std::complex<double> C[5] = {};
    for (auto& q : p) {
      if (q.heavy) continue;
      double phi = std::atan2(q.y - hy, q.x - hx);
      for (int m = 1; m <= 4; ++m) C[m] += q.m * std::exp(std::complex<double>(0, m * phi));
    }
    for (int m = 1; m <= 4; ++m) { cm[m] = std::abs(C[m]) / Mdisk; ph[m] = std::arg(C[m]); }
  };
  double cm0[5], ph0[5]; modes(cm0, ph0);
  printf("post-rebalance |c1..c4| = %.2e %.2e %.2e %.2e   (seed floor)\n",
         cm0[1], cm0[2], cm0[3], cm0[4]);

  // ---- coherent seeding: Lagrangian radial displacement R -> R (1 + eps cos(m phi)) ----
  if (m_seed > 0 && eps_seed > 0.0) {
    for (auto& q : p) {
      if (!q.gas) continue;
      double R = std::hypot(q.x, q.y), phi = std::atan2(q.y, q.x);
      double f = 1.0 + eps_seed * std::cos(m_seed * phi);
      q.x = R * f * std::cos(phi); q.y = R * f * std::sin(phi);   // velocities unchanged
    }
    modes(cm0, ph0);
    printf("seeded m=%d at eps=%.4f -> |c1..c4| = %.2e %.2e %.2e %.2e\n",
           m_seed, eps_seed, cm0[1], cm0[2], cm0[3], cm0[4]);
  }

  char fmname[128];
  std::snprintf(fmname, sizeof(fmname), "lb94/sph/quiet_modes_m%d_e%.4f.dat", m_seed, eps_seed);
  FILE* fm = fopen(fmname, "w");
  fprintf(fm, "# t  |c1| |c2| |c3| |c4|  phase1  |c1|(30-70AU) |c1|(70-120) |c1|(120-180)\n");
  auto t0 = std::chrono::steady_clock::now();
  BlockStepper bs; bs.dt0 = 0.05; bs.Rmax = 12;
  bs.init(p, tree, par);
  auto cb = [&](double t) {
    double cm[5], ph[5]; modes(cm, ph);
    double hx = p[star_i].x, hy = p[star_i].y;
    double ann[3] = {0, 0, 0}; double annM[3] = {0, 0, 0};
    std::complex<double> CA[3] = {};
    for (auto& q : p) {
      if (q.heavy) continue;
      double R = std::hypot(q.x - hx, q.y - hy);
      int b = (R > 0.3 && R < 0.7) ? 0 : (R > 0.7 && R < 1.2) ? 1 : (R > 1.2 && R < 1.8) ? 2 : -1;
      if (b < 0) continue;
      double phi = std::atan2(q.y - hy, q.x - hx);
      CA[b] += q.m * std::exp(std::complex<double>(0, phi)); annM[b] += q.m;
    }
    for (int b = 0; b < 3; ++b) ann[b] = std::abs(CA[b]) / std::max(annM[b], 1e-30);
    fprintf(fm, "%.4f %.6e %.6e %.6e %.6e %.4f %.6e %.6e %.6e\n",
            t, cm[1], cm[2], cm[3], cm[4], ph[1], ann[0], ann[1], ann[2]);
    fflush(fm);
    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    // momentum-conservation diagnostic: |P_tot| relative to the disk's typical momentum
    double px = 0, py = 0, pz = 0, mref = 0, comx = 0, comy = 0, Mtot = 0;
    for (auto& q : p) {
      px += q.m * q.vx; py += q.m * q.vy; pz += q.m * q.vz;
      mref += q.m * std::sqrt(q.vx * q.vx + q.vy * q.vy);
      comx += q.m * q.x; comy += q.m * q.y; Mtot += q.m;
    }
    double Prel = std::sqrt(px * px + py * py + pz * pz) / std::max(mref, 1e-300);
    printf("  t=%.2f  |c1|=%.2e |c2|=%.2e |c3|=%.2e  star_r=%.4f  |P|/P_ref=%.1e com=%.1e (%.0f s)\n",
           t, cm[1], cm[2], cm[3], std::hypot(p[star_i].x, p[star_i].y), Prel,
           std::hypot(comx, comy) / Mtot, wall);
    fflush(stdout);
  };
  bs.run(p, tree, par, t_end, 0.05, cb);
  fclose(fm);
  printf("done: t=%.1f T_unit\n", t_end);
  return 0;
}
