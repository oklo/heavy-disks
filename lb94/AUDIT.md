# LB94 collapse code — end-to-end audit / self-referee

Audit of the 2D axisymmetric radiation-hydro nested-grid collapse (`lb94/src/`) that
produces the LB94 (Laughlin & Bodenheimer 1994) protostellar-disk initial model, following
YBL93 (Yorke, Bodenheimer & Laughlin 1993) for the method. Written as an adversarial
self-review: what is solid, what is approximate, what is wrong.

## Verdict

The numerics that were *claimed* solid are solid: the Poisson solver, the conservative
van Leer transport, and the Berger–Colella refluxing are correct (the last is exact to
machine precision in isolation, `test_reflux.cpp`). Mass conservation is now **0.03%**
(physical) over 60 kyr. The collapse reproduces the YBL/LB94 structure: a ~0.73 M⊙ central
object plus a retained ~0.30 M⊙ disk to ~280 AU at ~20 K, vs LB94's 0.68 M⊙.

The weaknesses are in the **physics fidelity** (radiation, the central sink, stabilizers),
not the core discretization. None invalidate the qualitative result; several could move the
core mass at the ~10% level and should be understood before the SPH hand-off.

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
1. **Velocity ceiling `vmax = 100 km/s`** (`hydro.hpp:279`). An artificial cap on |v_R|,|v_Z|
   that bounds the central runaway. It is non-physical drag; real free-fall onto a 0.7 M⊙
   core exceeds 100 km/s inside ~few AU. It acts mostly inside R_sink (already being drained),
   but it is a stabilizer masking under-resolution rather than resolving it. **Right fix:**
   finer central resolution + a better sink so the ceiling never binds; then remove it.
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
4. **Resolution = 11 AU** (nlev=4, traded for a box 2× the cloud). LB94 resolved finer. The
   0.73 vs 0.68 M⊙ gap is plausibly in part under-resolution of the inner disk/accretion.

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
1. Push resolution to ~5.5 AU and retire the velocity ceiling; re-check core mass → 0.68.
2. Sensitivity scan over (`rho_ceiling`, `t_drain`) to bound the core/disk split.
3. Track the reflux floor in the conservation accounting (close issue 8).
4. Then the LB94 payoff: viscous (α≈0.01) vs inviscid runs → Fig. 1 (M_core(t)) and Fig. 2
   (Σ, T, j profiles).
