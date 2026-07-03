# The quiet-start null experiment — verdict

Setup (`sph/ic_quiet.cpp`, `sph/quiet.cpp`): ring-stratified realization of the pinned
basic-state family — constant-n rings (n=125, equal azimuthal spacing, random ring phases;
placement |c_m| ~ 5e-17), same-z opposite-azimuth pairs (no bending seed), collisionless
gravity-only star, per-ring v_phi rebalanced against the code's own measured forces,
theta = 0.4, no absorbing boundary.

## Two conservation bugs found by the null run (both fixed)

The first quiet run erupted in m=1 at gamma_eff ~ 30-50/T_unit from machine-level seeds —
far too fast to be physics. Diagnosis:
1. **Receiver-only gravitational softening** (eps = h_receiver): star-gas pairs near the
   center interacted with non-reciprocal forces — a net momentum pump at the exact place
   the forces are largest, i.e. an artificial dipole (m=1) driver. Fixed: symmetric pair
   softening eps_ij = (h_i+h_j)/2 at tree leaf level (tree.hpp).
2. **Gather-only pair forces**: pairs with 2h_j > r > 2h_i were felt by one member only.
   Momentum drifted coherently (|P| and the COM grew linearly). Fixed: symmetric neighbor
   lists (r < 2 max(h_i,h_j)) for density and forces (tree.hpp::sym_neighbors, forces.hpp).
Sedov re-validated after both fixes (d log R_s/d log t = 0.398).

## The verdict (fixed code, same physics, two ICs)

| | quiet start (smooth family) | Poisson-mapped grid IC |
|---|---|---|
| seed |c1| | ~1e-4 (tree/SPH noise floor) | ~7e-3 (sqrt-N Poisson) |
| t = 0.6 | ~2e-2 (floor) | **0.43** (erupting) |
| t = 3.0 | **7e-2, flat/diffusive** (mostly COM random-walk + neutral slosh) | (disrupts by ~1) |

**The quiet start does not grow for 3 T_unit — matching the linear-theory prediction**
(no robust growing mode of the smooth basic state; LINEAR_PREDICTIONS.md). **The
Poisson-mapped IC still erupts on the original schedule with the conservative code** —
the violent m=1 is not a code artifact, but it is also not a linear mode: if a uniform
exponential amplifier at the observed effective rate (~15/T_unit locally) existed, the
quiet run's 1e-4 seeds would have reached unity by t ~ 0.6. They did not in 3 T_unit.

## Interpretation

The system is **linearly stable but nonlinearly (subcritically) unstable**: the violent
m=1 requires finite-amplitude input. The Poisson realization supplies it — ~1e-2 mode
seeds and ~14% local density noise (at N_neigh = 50), plus the grid disk's sharp
plateau/cliff structure, the heavy *gas* central particle, and the unbalanced start.
This sharpens the reading of the 1994-style SPH stage: Monte-Carlo ICs at 25k particles
carry 2-3% noise — above threshold by construction. The onset scaling with 1/sqrt(N)
(M1_ANALYSIS.md) fits naturally: the trigger is the noise *amplitude* entering a
subcritical amplifier, not the seed of a linear eigenmode.

## The discriminating next pair

1. **Seeded quiet start**: coherent m=1 at 0.25/0.5/1/2% on the quiet IC — locate the
   subcritical threshold; check whether m=2 seeding also triggers or only m=1.
2. **Poisson-sampled smooth family** (same noise level as the grid IC, but no sharp
   features, collisionless star, rebalance): separates noise amplitude from the other
   realization differences (disk structure / heavy gas particle / imbalance).
