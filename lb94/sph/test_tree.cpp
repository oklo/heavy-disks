// Validate the Barnes-Hut tree: (a) gravity vs a direct O(N^2) sum with the same spline
// softening -- RMS force error should fall with theta and improve with the quadrupole;
// (b) the symmetric neighbor search must exactly match a brute-force list.
#include "tree.hpp"
#include <cstdio>
#include <random>
#include <set>

using namespace sph;

int main() {
  const int N = 3000;
  const double eps = 0.05;
  std::mt19937 rng(7);
  std::normal_distribution<double> Gauss(0.0, 0.3);
  std::vector<Particle> p(N);
  for (auto& q : p) {                                  // a centrally-concentrated blob
    q.x = Gauss(rng); q.y = Gauss(rng); q.z = Gauss(rng);
    q.m = 1.0 / N; q.h = eps;
  }

  // direct N-body acceleration (same spline softening) as the reference
  auto direct = [&](int i, double& ax, double& ay, double& az) {
    ax = ay = az = 0;
    for (int j = 0; j < N; ++j) {
      if (j == i) continue;
      double dx = p[i].x - p[j].x, dy = p[i].y - p[j].y, dz = p[i].z - p[j].z;
      double r = std::sqrt(dx * dx + dy * dy + dz * dz), gf = grav_factor(r, eps);
      ax += -p[j].m * dx * gf; ay += -p[j].m * dy * gf; az += -p[j].m * dz * gf;   // G=1
    }
  };

  printf("Tree gravity vs direct (N=%d, G=1):\n", N);
  printf("  theta   monopole-only      monopole+quad   (RMS |da|/|a|)\n");
  for (double th : {1.0, 0.7, 0.5, 0.3}) {
    double err[2];
    for (int quad = 0; quad < 2; ++quad) {
      Tree T; T.P = &p; T.G = 1.0; T.theta = th; T.use_quad = quad; T.build();
      double s2 = 0, a2 = 0;
      for (int i = 0; i < N; i += 3) {                 // sample every 3rd particle
        double dax, day, daz, tx, ty, tz;
        direct(i, dax, day, daz);
        T.accel(i, tx, ty, tz);
        double ex = tx - dax, ey = ty - day, ez = tz - daz;
        s2 += ex * ex + ey * ey + ez * ez; a2 += dax * dax + day * day + daz * daz;
      }
      err[quad] = std::sqrt(s2 / a2);
    }
    printf("  %.2f    %10.2e        %10.2e\n", th, err[0], err[1]);
  }

  // neighbor search vs brute force (adaptive h first)
  density_bruteforce(p, 50);                            // sets h per particle
  Tree T; T.P = &p; T.build();
  int mism = 0, checked = 0;
  std::vector<int> nb;
  for (int i = 0; i < N; i += 50) {
    T.neighbors(i, nb);
    std::set<int> tree_set(nb.begin(), nb.end()), bf_set;
    for (int j = 0; j < N; ++j) {
      if (j == i) continue;
      double r = std::sqrt(dist2(p[i], p[j]));
      if (r < 2 * p[i].h || r < 2 * p[j].h) bf_set.insert(j);
    }
    if (tree_set != bf_set) ++mism;
    ++checked;
  }
  printf("Neighbor search: %d/%d sampled particles match brute force  %s\n",
         checked - mism, checked, mism == 0 ? "PASS" : "FAIL");
  return 0;
}
