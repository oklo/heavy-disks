# LB94 collapse code — end-to-end audit / self-referee

Audit of the 2D axisymmetric radiation-hydro nested-grid collapse (`lb94/src/`) that
produces the LB94 (Laughlin & Bodenheimer 1994) protostellar-disk initial model, following
YBL93 (Yorke, Bodenheimer & Laughlin 1993) for the method. Written as an adversarial
self-review: what is solid, what is approximate, what is wrong.

## Verdict

The numerics that were *claimed* solid are solid: the Poisson solver, the conservative
van Leer transport, and the Berger–Colella refluxing are correct (the last is exact to
machine precision in isolation, `test_reflux.cpp`). Mass conservation is **0.02%** (physical)
over 60 kyr. At 5.6 AU resolution the collapse reproduces the YBL/LB94 structure to the
digit: a **0.681 M⊙** central object (LB94 0.68) plus a retained 0.33 M⊙ disk to 253 AU
(R_kep 250) at ~20 K.

The remaining weaknesses are in **physics fidelity** (radiation, the central sink), not the
core discretization. The two original majors — the velocity ceiling and resolution — are now
resolved (see below). None of the rest invalidate the result; the gray-cooling simplification
and the sink parameters are the items to understand before the SPH hand-off.

## Solid / validated
- **Poisson** (`poisson.hpp`): cylindrical FD + red-black SOR, isolated multipole (l=0,2,4)
  boundary; validated against the analytic uniform sphere. Nested Dirichlet solve takes its
  outer boundary from the parent potential.
- **Transport** (`hydro.hpp::transport`): van Leer 2nd-order, flux-conservative form; the
  reflux unit test confirms machine-precision conservation under pure advection.
- **Refluxing** (`nested.hpp::reflux`): analytically derived for the cylindrical factor-2
  interface and verified exact; the earlier apparent "residual" was a physics leak, not a
  refluxing error.
- **Energy equation**: Sedov–Taylor blast gives d log R_s/d log t = 0.396 (Sedov 0.40; YBL 0.39).
- **Conservation**: closed outer boundary + floor-mass accounting → 0.03% over 60 kyr.

## Issues, by severity

### Major (could shift the core mass / disk at the ~10% level)
1. **Velocity ceiling — RESOLVED.** Was an artificial 100 km/s cap on |v_R|,|v_Z|. At 5.6 AU
   resolution a run with the ceiling OFF gives a bit-for-bit identical trajectory (same
   39237 steps, same core mass) — i.e. the sink already bounds the center and the ceiling
   never triggered. It is now off (`vmax = 1e300`); confirmed not masking under-resolution.
2. **Radiative cooling is a crude gray relaxation** (`radiative_cooling`), not the
   flux-limited diffusion LB94/YBL used. e relaxes toward ρ c_v T_amb on an optical-depth-
   suppressed time; the opacity scale height uses the Jeans length, not c_s/Ω. The disk lands
   at ~20 K (right ballpark) but the radial/vertical *thermal structure* — which sets H, Q,
   and fragmentation susceptibility — is only approximate. This is the largest physics
   simplification relative to the paper.
3. **Central sink has free parameters** (`rho_ceiling = 5e-14`, `t_drain = 300 yr`,
   `R_sink = 7 dR`). These set how fast mass leaves the grid for the core and therefore the
   core/disk split directly. `rho_ceiling` must stay above the physical disk density (it does
   here — the disk is retained), but the core mass is parameter-sensitive and not derived.
4. **Resolution — RESOLVED.** Was 11 AU (nlev=4). Now 5.6 AU (nlev=5): core mass 0.681 M⊙
   (LB94 0.68), disk R_out 253 AU (R_kep 250), conservation 0.02%, floor artifact 0.10%.
   Confirms the 0.73 vs 0.68 gap was inner-disk/accretion under-resolution.

### Minor (real but small)
5. **First-order operator splitting + lagged gravity**: `step()` does Lie splitting
   (source → transport, full dt each) and solves gravity once per step on the *start*-of-step
   density. Both are first-order in time; a Strang split + a mid-step gravity re-solve would
   be 2nd order. For a smooth collapse the phase error is small but present.
6. **One-sided gradients at the reflecting axis/midplane** (`source`): for j=0 (and i=0,
   and the outer faces) the pressure/gravity gradient uses `(P[1]-P[0])/dZ`, which for the
   reflect boundary is ~2× the symmetric value `(P[1]-P[0])/(2dZ)`. Over-compresses the single
   midplane/axis cell layer. One-cell effect; minor but a genuine discretization error.
7. **Artificial viscosity has no linear (von Neumann–Landshoff) term** — only the quadratic
   `Cq ρ (L div v)²`. Some post-shock ringing; fine for this smooth-ish collapse.
8. **Reflux density floor not tracked** (`nested.hpp:163-166`): the post-reflux `rho<floor`
   clamp can add mass that is *not* included in `dbg_mass_floored`, so it is not subtracted in
   the conservation metric. Negligible at the current 0.03% level but an accounting gap.
9. **Synchronous time stepping**: all levels share the global-min dt (the original ENG
   subcycles the fine grids). Correct, just inefficient — the coarse grids are over-stepped.
10. **Multipole boundary truncated at l=4**: adequate because the coarsest boundary sits 2×
    the cloud radius away and interior grids inherit the parent potential, but a flattened
    disk near a boundary would want higher l.

### Dead code to remove
- `NestedHydro::conserve_mass` (the old renormalization stopgap, superseded by refluxing).
- `freeze_overlap` path and the `debug_mass` accumulators (`dloss_*`, `dbg_ghost_in`,
  `dbg_reflux_added`, `dloss_gap`) — diagnostics; keep `dbg_mass_floored` (now used).

## Recommended next steps
1. ~~Push resolution to ~5.5 AU and retire the velocity ceiling~~ — DONE: core mass 0.681 M⊙;
   ceiling confirmed inactive and switched off.
2. Sensitivity scan over (`rho_ceiling`, `t_drain`) to bound the core/disk split.
3. Track the reflux floor in the conservation accounting (close issue 8).
4. Then the LB94 payoff: viscous (α≈0.01) vs inviscid runs → Fig. 1 (M_core(t)) and Fig. 2
   (Σ, T, j profiles).
