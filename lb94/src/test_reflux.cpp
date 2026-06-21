// Isolate which physics op breaks mass conservation: uniform inflow on a 2-level grid,
// adding features one at a time. composite_mass + M_c must be conserved in every case.
#include "nested.hpp"
#include <cstdio>

using namespace lb94;

double run(const char* name, bool gravity, bool energy, bool cooling, bool sink, bool vcap,
           bool closed = false) {
  const int NR = 16, NZ = 16, nlev = 2;
  const double L = 2.0;
  NestedHydro nh(nlev, NR, NZ, L, gravity ? 1.0 : 0.0, /*cs2=*/energy ? 0.1 : 0.0,
                 /*rho_crit=*/1e300, /*gam=*/1.4);
  nh.lev[0].closed_outer = closed;
  for (auto& h : nh.lev) {
    h.use_energy = energy; h.radiation = cooling; h.Tamb = 10; h.mu_gas = 1; h.kap0 = 1e-3;
    if (vcap) h.vmax = 0.5;
    for (int c = 0; c < NR * NZ; ++c) {
      h.rho[c] = 1.0; h.sR[c] = -0.3; h.sZ[c] = -0.3;
      if (energy) h.e[c] = 1.0;
    }
  }
  if (sink) { nh.finest().rho_ceiling = 1.5; nh.finest().R_sink = 0.3; nh.finest().t_drain = 0.1; }

  auto total = [&]() { return nh.composite_mass() + nh.finest().M_c; };
  double m0 = total(), dt = 0.02;
  for (int n = 0; n < 30; ++n) nh.step(dt);
  double drift = total() / m0 - 1.0;
  printf("  %-28s drift = %+.3e\n", name, drift);
  return drift;
}

int main() {
  printf("conservation of (composite gas + M_c) after 30 steps of uniform inflow:\n");
  run("transport only", false, false, false, false, false);
  run("+ gravity", true, false, false, false, false);
  run("+ energy (open bndry)", false, true, false, false, false);
  run("+ energy (closed bndry)", false, true, false, false, false, /*closed=*/true);
  run("+ sink", false, false, false, true, false);
  run("+ velocity ceiling", false, false, false, false, true);
  run("everything (closed)", true, true, true, true, true, /*closed=*/true);
  return 0;
}
