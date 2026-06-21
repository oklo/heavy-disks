# Independent code-review prompt (hand to GPT-5.5 / Codex)

You are an expert computational astrophysicist and numerical-methods referee. Perform an
**independent, adversarial review** of a from-scratch reproduction of the protostellar-disk
initial model from **Laughlin & Bodenheimer 1994 (ApJ 436, 335)**, built on the collapse
method of **Yorke, Bodenheimer & Laughlin 1993 (ApJ 411, 274)** ("YBL93"). Do not assume the
code is correct. Your job is to find errors, unjustified approximations, and conservation or
consistency bugs — and to judge whether the result faithfully reproduces the target physics.

## What the code is supposed to do
A 2D axisymmetric (cylindrical R,Z) self-gravitating radiation-hydrodynamics collapse of a
1 M⊙, ~2700 AU, ρ∝1/r, 20 K, uniformly rotating (Ω=4.4e-13 s⁻¹) cloud, on **Explicit Nested
Grids** (static factor-2 mesh refinement sharing the origin, Berger–Colella style). It should
collapse in ~one free-fall time (~24 kyr here, run to ~60 kyr) to a central protostellar
object plus a rotationally-supported disk — the model LB94 then hand to a 3D SPH code.
Target numbers (YBL/LB94): central object ~0.5–0.7 M⊙ by ~60–80 kyr, disk Keplerian radius
of the outermost element ~250 AU, disk T~20 K, mass conservation at the sub-percent level.

## Source layout (all in `lb94/src/`, C++17, header-only)
- `poisson.hpp` — cylindrical Poisson (∇²Φ=4πGρ) via red-black SOR; isolated multipole
  (l=0,2,4) outer boundary; `solve_dirichlet` for nested interior grids (boundary potential
  supplied from the parent).
- `hydro.hpp` — single-grid hydro: conservative state (ρ, ρv_R, ρv_Z, A=ρRv_φ, e). van Leer
  2nd-order directionally-split transport; operator-split source (pressure + self-gravity +
  central point-mass + centrifugal); von Neumann artificial viscosity; ideal-gas energy
  equation with PdV + AV heating; simplified gray radiative cooling; optional Shakura–Sunyaev
  α-viscosity; a density-cap + radius-drain central sink onto (M_c, J_c); stabilizers (density
  floor, density clamp, velocity ceiling).
- `nested.hpp` — the nested-grid driver: conservative fine→coarse restriction, coarse→fine
  potential solve, parent→child ghost prolongation, **Berger–Colella refluxing**, synchronous
  (single global-dt) time stepping; `composite_mass()` (finest-resolved total).
- `collapse_ybl.cpp` — the YBL standard-case driver in cgs.
- Validation drivers: `test_poisson.cpp` (vs analytic sphere), `test_freefall.cpp`,
  `test_nested.cpp`, `test_sedov.cpp` (Sedov–Taylor, should give d log R_s/d log t ≈ 0.4),
  `test_reflux.cpp` (isolated mass-conservation harness for the refluxing).

Build any driver with: `g++ -O2 -std=c++17 -Isrc src/<driver>.cpp -o <bin>`.

## What to scrutinize (in priority order)
1. **Conservation**: is mass/momentum/angular-momentum/energy actually conserved by the
   transport + restriction + prolongation + refluxing chain? Derive the cylindrical factor-2
   refluxing relations yourself and check them against `nested.hpp::reflux`. Check the
   restriction and ghost prolongation are conservative and mutually consistent.
2. **The central sink** (`hydro.hpp::accrete`): is the mass/angular-momentum removed to the
   core exactly what leaves the gas? Are the free parameters (rho_ceiling, t_drain, R_sink)
   physically justified, and could they bias the core/disk mass split?
3. **Stabilizers**: the velocity ceiling, density clamp, and density floor. Which of these are
   physically benign and which could alter the dynamics or hide under-resolution? Quantify.
4. **Self-gravity**: correctness of the cylindrical Poisson discretization (axis/midplane
   reflecting BCs, the multipole boundary), the nested Dirichlet coupling, and the time-level
   at which gravity is evaluated relative to the source/transport (splitting order).
5. **Discretization at boundaries**: the reflecting axis (R=0) and midplane (Z=0) — are the
   gradients (pressure, gravity) and fluxes symmetric/consistent there?
6. **Thermal physics**: the gray cooling vs proper radiation transport — what does the
   simplification cost for the disk's thermal/vertical structure, and hence the LB94 result?
7. **Operator-splitting / time integration order** and the synchronous-vs-subcycled tradeoff.
8. **Truelove resolution** of the collapsing core: is the Jeans length resolved, or is the
   sink/ceiling compensating for under-resolution?

## Deliverable
A referee report: (a) confirmed-correct components; (b) bugs (with file:line and a proposed
fix); (c) approximations ranked by their likely effect on the core mass and disk structure;
(d) a short list of the highest-value changes to improve fidelity to LB94/YBL. Run the
validation drivers if you can and report the numbers. Be specific and quantitative; cite the
paper equations where relevant. Do not soft-pedal problems.
