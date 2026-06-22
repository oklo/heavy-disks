// Barnes-Hut octree (Hernquist & Katz 1989): gravity with monopole + quadrupole moments and
// the s/d < theta opening criterion, spline-softened at small separations; plus a symmetric
// (gather + scatter) nearest-neighbor search used by the SPH estimators. Top-down build.
#ifndef SPH_TREE_HPP
#define SPH_TREE_HPP
#include <vector>
#include <cmath>
#include <algorithm>
#include "sph.hpp"

namespace sph {

struct Tree {
  std::vector<Particle>* P = nullptr;
  double theta = 0.7;            // opening angle (HK89 use ~0.7-1.0)
  double G = 1.0;
  bool use_quad = true;          // include the quadrupole moment (HK89)

  struct Node {
    double cx, cy, cz, half;     // geometric center and half-side of the cubic cell
    double mass, comx, comy, comz;
    double qxx, qyy, qzz, qxy, qxz, qyz;   // traceless quadrupole about the com
    double hmax;                 // max smoothing length in the subtree (for neighbor search)
    int child[8];                // child node indices, -1 if empty
    int pidx;                    // particle index if a single-particle leaf, else -1
  };
  std::vector<Node> nodes;
  int root = -1;

  void build() {
    nodes.clear();
    int n = P->size();
    if (n == 0) { root = -1; return; }
    double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
    for (auto& q : *P) {
      lo[0] = std::min(lo[0], q.x); hi[0] = std::max(hi[0], q.x);
      lo[1] = std::min(lo[1], q.y); hi[1] = std::max(hi[1], q.y);
      lo[2] = std::min(lo[2], q.z); hi[2] = std::max(hi[2], q.z);
    }
    double cx = 0.5 * (lo[0] + hi[0]), cy = 0.5 * (lo[1] + hi[1]), cz = 0.5 * (lo[2] + hi[2]);
    double half = 0.5 * std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2], 1e-30}) * 1.0001;
    std::vector<int> idx(n);
    for (int i = 0; i < n; ++i) idx[i] = i;
    nodes.reserve(2 * n);
    root = build_node(idx, cx, cy, cz, half);
  }

  int build_node(std::vector<int>& idx, double cx, double cy, double cz, double half) {
    int ni = nodes.size();
    nodes.emplace_back();
    Node nd{}; nd.cx = cx; nd.cy = cy; nd.cz = cz; nd.half = half; nd.pidx = -1;
    for (int k = 0; k < 8; ++k) nd.child[k] = -1;
    if (idx.size() == 1) {
      const Particle& q = (*P)[idx[0]];
      nd.pidx = idx[0]; nd.mass = q.m; nd.comx = q.x; nd.comy = q.y; nd.comz = q.z; nd.hmax = q.h;
      nodes[ni] = nd; return ni;
    }
    std::vector<int> oct[8];
    for (int p : idx) {
      const Particle& q = (*P)[p];
      int o = (q.x > cx) | ((q.y > cy) << 1) | ((q.z > cz) << 2);
      oct[o].push_back(p);
    }
    double q4 = 0.5 * half;
    int children[8];
    for (int o = 0; o < 8; ++o) {
      if (oct[o].empty()) { children[o] = -1; continue; }
      double ox = cx + ((o & 1) ? q4 : -q4), oy = cy + ((o & 2) ? q4 : -q4), oz = cz + ((o & 4) ? q4 : -q4);
      children[o] = build_node(oct[o], ox, oy, oz, q4);   // recursion may realloc `nodes`
    }
    // accumulate mass, com, hmax from children (indices are realloc-safe)
    double mass = 0, mx = 0, my = 0, mz = 0, hmax = 0;
    for (int o = 0; o < 8; ++o)
      if (children[o] >= 0) {
        Node& c = nodes[children[o]];
        mass += c.mass; mx += c.mass * c.comx; my += c.mass * c.comy; mz += c.mass * c.comz;
        hmax = std::max(hmax, c.hmax);
      }
    double comx = mx / mass, comy = my / mass, comz = mz / mass;
    // quadrupole about this com via the parallel-axis theorem
    double qxx = 0, qyy = 0, qzz = 0, qxy = 0, qxz = 0, qyz = 0;
    for (int o = 0; o < 8; ++o)
      if (children[o] >= 0) {
        Node& c = nodes[children[o]];
        double Dx = c.comx - comx, Dy = c.comy - comy, Dz = c.comz - comz, D2 = Dx * Dx + Dy * Dy + Dz * Dz;
        qxx += c.qxx + c.mass * (3 * Dx * Dx - D2);
        qyy += c.qyy + c.mass * (3 * Dy * Dy - D2);
        qzz += c.qzz + c.mass * (3 * Dz * Dz - D2);
        qxy += c.qxy + c.mass * (3 * Dx * Dy);
        qxz += c.qxz + c.mass * (3 * Dx * Dz);
        qyz += c.qyz + c.mass * (3 * Dy * Dz);
      }
    nd.mass = mass; nd.comx = comx; nd.comy = comy; nd.comz = comz; nd.hmax = hmax;
    nd.qxx = qxx; nd.qyy = qyy; nd.qzz = qzz; nd.qxy = qxy; nd.qxz = qxz; nd.qyz = qyz;
    for (int o = 0; o < 8; ++o) nd.child[o] = children[o];
    nodes[ni] = nd; return ni;
  }

  // gravitational acceleration on particle i (spline-softened with eps = h_i)
  void accel(int i, double& ax, double& ay, double& az) {
    ax = ay = az = 0.0;
    if (root >= 0) walk_grav(root, i, (*P)[i].h, ax, ay, az);
  }

  void walk_grav(int ni, int i, double eps, double& ax, double& ay, double& az) {
    Node& nd = nodes[ni];
    if (nd.mass == 0.0 || nd.pidx == i) return;     // empty, or i's own leaf
    double dx = (*P)[i].x - nd.comx, dy = (*P)[i].y - nd.comy, dz = (*P)[i].z - nd.comz;
    double r2 = dx * dx + dy * dy + dz * dz, r = std::sqrt(r2);
    double size = 2.0 * nd.half;
    if (nd.pidx >= 0 || size < theta * r) {          // accept: leaf, or s/d < theta
      double gf = grav_factor(r, eps);               // = [M(<r)/M]/r^3 (Newtonian for r>2eps)
      ax += -G * nd.mass * dx * gf; ay += -G * nd.mass * dy * gf; az += -G * nd.mass * dz * gf;
      if (use_quad && nd.pidx < 0 && r > 0) {         // quadrupole correction (dr = r_i - com)
        double r5 = r2 * r2 * r, r7 = r5 * r2;
        double Qx = nd.qxx * dx + nd.qxy * dy + nd.qxz * dz;
        double Qy = nd.qxy * dx + nd.qyy * dy + nd.qyz * dz;
        double Qz = nd.qxz * dx + nd.qyz * dy + nd.qzz * dz;
        double drQdr = dx * Qx + dy * Qy + dz * Qz;
        double c2 = 2.5 * drQdr / r7;
        ax += G * (Qx / r5 - c2 * dx); ay += G * (Qy / r5 - c2 * dy); az += G * (Qz / r5 - c2 * dz);
      }
      return;
    }
    for (int o = 0; o < 8; ++o)
      if (nd.child[o] >= 0) walk_grav(nd.child[o], i, eps, ax, ay, az);
  }

  // count of particles j!=i within radius R of i (gather; used for the h iteration)
  int gather_count(int i, double R) {
    int cnt = 0;
    if (root >= 0) gwalk(root, i, R, cnt);
    return cnt;
  }
  void gwalk(int ni, int i, double R, int& cnt) {
    Node& nd = nodes[ni];
    if (nd.mass == 0.0) return;
    double dx = (*P)[i].x - nd.cx, dy = (*P)[i].y - nd.cy, dz = (*P)[i].z - nd.cz;
    double d = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (d - nd.half * 1.7320508 > R) return;          // node entirely outside R
    if (nd.pidx >= 0) {
      if (nd.pidx != i && std::sqrt(dist2((*P)[i], (*P)[nd.pidx])) < R) ++cnt;
      return;
    }
    for (int o = 0; o < 8; ++o)
      if (nd.child[o] >= 0) gwalk(nd.child[o], i, R, cnt);
  }

  // symmetric neighbor list of i: all j with r_ij < 2 h_i OR r_ij < 2 h_j
  void neighbors(int i, std::vector<int>& out) {
    out.clear();
    if (root >= 0) walk_nbr(root, i, out);
  }

  void walk_nbr(int ni, int i, std::vector<int>& out) {
    Node& nd = nodes[ni];
    if (nd.mass == 0.0) return;
    double dx = (*P)[i].x - nd.cx, dy = (*P)[i].y - nd.cy, dz = (*P)[i].z - nd.cz;
    double d = std::sqrt(dx * dx + dy * dy + dz * dz);
    double reach = 2.0 * std::max((*P)[i].h, nd.hmax) + nd.half * 1.7320508;  // + cell half-diagonal
    if (d > reach) return;
    if (nd.pidx >= 0) {
      int j = nd.pidx;
      if (j != i) {
        double r = std::sqrt(dist2((*P)[i], (*P)[j]));
        if (r < 2.0 * (*P)[i].h || r < 2.0 * (*P)[j].h) out.push_back(j);
      }
      return;
    }
    for (int o = 0; o < 8; ++o)
      if (nd.child[o] >= 0) walk_nbr(nd.child[o], i, out);
  }
};

}  // namespace sph
#endif
