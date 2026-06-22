// LB94 3D TreeSPH driver: read the disk IC (from the 2D collapse), evolve the locally-
// isothermal, self-gravitating disk + central heavy particle, with an absorbing boundary at
// R = 230 AU, and track the azimuthal Fourier mode amplitudes |c_m|(t) (LB94 eq. 7) that
// measure the growth of the m=1/2 spiral instabilities. Code units: L=100 AU, V=1 km/s, G=1.
#include "block.hpp"
#include <cstdio>
#include <cmath>
#include <complex>
#include <chrono>
#include <string>

using namespace sph;

int main(int argc, char** argv) {
  const char* icfile = (argc > 1) ? argv[1] : "lb94/sph/disk_ic.dat";
  const double t_end = (argc > 2) ? std::atof(argv[2]) : 10.0;   // in T_unit (LB94 ran 10)
  const double R_edge = 2.30;                                    // 230 AU absorbing boundary

  // read the IC (code units)
  std::vector<Particle> p;
  FILE* fp = fopen(icfile, "r");
  if (!fp) { fprintf(stderr, "cannot open %s\n", icfile); return 1; }
  char line[512];
  while (fgets(line, sizeof(line), fp)) {
    if (line[0] == '#') continue;
    Particle q{}; int heavy;
    if (sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %d", &q.x, &q.y, &q.z,
               &q.vx, &q.vy, &q.vz, &q.m, &q.cs2_iso, &heavy) == 9) {
      q.heavy = heavy; q.h = 0.04; p.push_back(q);   // ~4 AU initial guess (refined below)
    }
  }
  fclose(fp);
  int heavy_i = -1;
  double Mdisk = 0;
  for (int i = 0; i < (int)p.size(); ++i) { if (p[i].heavy) heavy_i = i; else Mdisk += p[i].m; }

  SPHParams par; par.isothermal = true; par.selfgrav = true; par.N_neigh = 50;
  Tree tree; tree.P = &p; tree.theta = 0.7;
  adapt_h(p, tree, par, 30);                                      // converge h from the IC guess
  compute_forces(p, tree, par);

  // azimuthal Fourier amplitudes |c_m| = |sum_p m_p e^{i m phi}| / M_disk, in the heavy frame
  auto modes = [&](double cm[5]) {
    std::complex<double> C[5] = {};
    double hx = p[heavy_i].x, hy = p[heavy_i].y;
    for (auto& q : p) {
      if (q.heavy) continue;
      double phi = std::atan2(q.y - hy, q.x - hx);
      for (int m = 0; m <= 4; ++m) C[m] += q.m * std::exp(std::complex<double>(0, m * phi));
    }
    for (int m = 0; m <= 4; ++m) cm[m] = std::abs(C[m]) / Mdisk;
  };

  // initial surface density and Toomre Q in cylindrical-radius bins
  printf("LB94 TreeSPH: N=%zu (disk M=%.3f, heavy M=%.3f), R_edge=%.0f AU, evolve %.1f T_unit\n",
         p.size(), Mdisk, p[heavy_i].m, R_edge * 100, t_end);
  {
    const int NB = 20; double Rmax = R_edge;
    double mb[NB] = {}, csb[NB] = {}, omb[NB] = {}; int cnt[NB] = {};
    double hx = p[heavy_i].x, hy = p[heavy_i].y;
    for (auto& q : p) {
      if (q.heavy) continue;
      double R = std::hypot(q.x - hx, q.y - hy); int b = (int)(R / Rmax * NB);
      if (b < 0 || b >= NB) continue;
      double vphi = (-(q.vx) * (q.y - hy) + (q.vy) * (q.x - hx)) / std::max(R, 1e-9);
      mb[b] += q.m; csb[b] += std::sqrt(q.cs2_iso); omb[b] += vphi / std::max(R, 1e-9); ++cnt[b];
    }
    printf("  R(AU)   Sigma     Q(Toomre)\n");
    for (int b = 0; b < NB; ++b) {
      if (cnt[b] < 5) continue;
      double R = (b + 0.5) / NB * Rmax, dR = Rmax / NB;
      double Sig = mb[b] / (2 * PI * R * dR);            // code-unit surface density
      double cs = csb[b] / cnt[b], Om = omb[b] / cnt[b];
      double Q = cs * Om / (PI * Sig);                  // kappa ~ Omega (near-Keplerian), G=1
      printf("  %5.0f  %8.3f   %6.2f\n", R * 100, Sig, Q);
    }
  }

  double cm0[5]; modes(cm0);
  printf("initial modes |c1..c4| = %.3f %.3f %.3f %.3f\n", cm0[1], cm0[2], cm0[3], cm0[4]);

  // evolve with the block-timestep integrator; mode file tagged by N, snapshots optional
  bool do_snap = (argc > 3 && std::string(argv[3]) == "snap");
  char fmname[64]; std::snprintf(fmname, sizeof(fmname), "lb94/sph/modes_N%05zu.dat", p.size());
  FILE* fm = fopen(fmname, "w");
  fprintf(fm, "# t(T_unit)  |c1| |c2| |c3| |c4|   N=%zu\n", p.size());
  auto t0 = std::chrono::steady_clock::now();
  BlockStepper bs; bs.dt0 = 0.05; bs.Rmax = 12;
  bs.init(p, tree, par);
  printf("block stepper: dt0=%.3f, dt_min=%.2e (Rmax=%d)\n", bs.dt0, bs.dt_min(), bs.Rmax);
  bs.rung_histogram(p);
  auto cb = [&](double t) {
    // absorbing boundary: damp outward radial velocity beyond R_edge
    double hx = p[heavy_i].x, hy = p[heavy_i].y;
    for (auto& q : p) {
      if (q.heavy) continue;
      double R = std::hypot(q.x - hx, q.y - hy);
      if (R > R_edge) {
        double vr = (q.vx * (q.x - hx) + q.vy * (q.y - hy)) / std::max(R, 1e-9);
        if (vr > 0) { q.vx -= vr * (q.x - hx) / R; q.vy -= vr * (q.y - hy) / R; }
      }
    }
    double cm[5]; modes(cm);
    fprintf(fm, "%.4f %.4f %.4f %.4f %.4f\n", t, cm[1], cm[2], cm[3], cm[4]); fflush(fm);
    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    printf("  t=%.2f  |c1|=%.3f |c2|=%.3f |c3|=%.3f  (%.0f s)\n", t, cm[1], cm[2], cm[3], wall);
    fflush(stdout);
    // face-on snapshot (x,y,rho in the heavy-particle frame) for overdensity maps (LB94 Fig 5)
    if (!do_snap) return;
    static int isnap = 0;
    char sf[64]; std::snprintf(sf, sizeof(sf), "lb94/sph/snap_%02d.dat", isnap++);
    FILE* fs = fopen(sf, "w");
    fprintf(fs, "# t=%.3f  x y z rho  (heavy frame)\n", t);
    double hz = p[heavy_i].z;
    for (auto& q : p)
      if (!q.heavy) fprintf(fs, "%.4f %.4f %.4f %.4e\n", q.x - hx, q.y - hy, q.z - hz, q.rho);
    fclose(fs);
  };
  bs.run(p, tree, par, t_end, 0.1, cb);
  fclose(fm);
  printf("done: t=%.1f T_unit\n", t_end);
  return 0;
}
