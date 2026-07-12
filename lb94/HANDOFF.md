# HANDOFF — LB94 reproduction & the m=1 investigation (state as of 2026-07-12)

Read this first on re-up. It is the single entry point; per-topic detail lives in the
docs it points to. The repo is `oklo/heavy-disks` (private); everything below is
committed through the seeded-batch data.

## 0. One paragraph of context

The `heavy-disks` repo reproduces GL's trilogy ARS89 → LB94 → LKA98. ARS89 and LKA98 are
done. LB94 is done end-to-end (2D nested-grid radiation-hydro collapse → 3D TreeSPH →
resolution audit) and the project has since evolved into its most interesting phase: a
Toomre-style controlled-experiment program interrogating whether the violent m=1 spiral
in the SPH stage is disk physics or a property of noisy realizations. Current verdict:
**the smooth disk is linearly stable (eigenproblem + quiet start agree); the violence is
subcritical/finite-amplitude**; the seeded-quiet-start ladder (run overnight 2026-07-03,
data in repo, ANALYSIS PENDING) shows a slow (~0.5–1/T_unit) m=1 amplifier awakened by
any finite seed regardless of azimuthal number.

## 1. The scientific narrative (chronological, with doc pointers)

1. **2D collapse** (`src/`, `PLAN.md`, `AUDIT.md`): nested-grid radiation-hydro,
   Berger–Colella refluxing exact, physical conservation 0.02%, irradiated cooling
   (L_acc ≈ 25 L⊙ → disk T gradient 235→50 K) + Truelove-limited sink. Reproduces LB94
   Fig 1/2 well (`COMPARISON.md`, `cmp_fig*.png`); inviscid core 0.47 vs LB94 0.42 M⊙.
   The SPH IC epoch is the 40 kyr **viscous** disk (min Q = 1.24; the inviscid disk's
   Q = 0.68 was a replication error, fixed).
2. **TreeSPH** (`sph/`): HK89-faithful (kernel, BH tree + quadrupole, Monaghan AV,
   block timesteps). Validations: lattice density 0.14%, tree vs direct, 3D Sedov 0.398.
3. **Resolution audit** (`sph/plot_audit.py`, `fig_audit.png`): peak |c1| falls
   monotonically with N (0.97→0.39 for 6.25k→50k). First interpretation (converging to
   stability) was WRONG — see 4.
4. **M1_ANALYSIS.md**: noise floor = Poisson sqrt(pi/4N) exactly; onset shifts by
   ln(sqrt(N/N0))/gamma — fixed-time amplitude is a *seed* artifact. ARS89's note in
   proof read: their half-integer SLING modes (incl. the canonical eigenvalue) are
   stabilized at full edge correction (c=1) — do not attribute anything to SLING.
   Swing: X(m=1) ≈ 3–4 (weak), X(m=2) ≈ 1.5–2 (strong) yet m=1 dominates → not swing.
5. **Toomre program** (`TOOMRE_PROGRAM.md`, `paper/Toomre81.pdf`): controlled basic
   state + quiet starts + one-knob surgery. Basic state pinned (`fit_basic_state.py`):
   Sigma = 48.4 (R/100AU)^-2.51 exp[-(54.2AU/R)-(R/267AU)^2] g/cm², c_s = 0.534
   (R/100AU)^-1/4 km/s (q=1/2 = irradiation), M_*=0.340, M_d=0.374 M⊙. Predicts (not
   fits) v_phi to 6.4% (needs eps=H/2 thickness), Q to 7.8%, Q_min 1.35.
   Taper indices are NOT selected by the data (12% ring-structure rms floor) — they are
   dials. Family Q_min at 56 AU vs grid 32 AU (the smoothed 19–31 AU plateau/cliff is
   the first candidate "groove" feature to reintroduce deliberately).
6. **Linear eigenproblem** (`lineig.py`, `lineig_check.py`, `LINEAR_PREDICTIONS.md`,
   `fig_lineig.png`): ars89 machinery generalized to the family (regression EXACT after
   matching eta_eq=0.05). **No robust growing mode**: m=1 gamma ≤ 0.04/T_unit (outer-
   taper residents; die with domain/thickness); one m=2 inner-cavity mode 0.83/T_unit
   (dies with r_in move or thickness); argument-principle winding over gamma=0.5–5/Tu
   is exactly ZERO. Eigenvalues are exquisitely sensitive to equilibrium details
   (0.13 shift from a 0.05→0.1 softening change) — the note-in-proof lesson recurs.
7. **Quiet start** (`sph/ic_quiet.cpp`, `sph/quiet.cpp`, `QUIET_START.md`): ring-
   stratified (constant-n rings, placement |c_m| ~ 5e-17), collisionless star
   (Particle::gas flag), per-ring v_phi rebalance. The null run EXPOSED TWO
   CONSERVATION BUGS (both fixed, Sedov revalidated): receiver-only grav softening
   (non-reciprocal star–gas forces = artificial dipole pump) → symmetric leaf pair
   softening; gather-only pair forces (one-sided pressure for 2h_j > r > 2h_i) →
   `tree.sym_neighbors` for density+forces. FIXED-code verdict: quiet start stable
   3 T_unit; Poisson grid-IC still erupts on schedule (0.43 by t=0.6) → **linearly
   stable, subcritically unstable**.

## 2. OPEN THREAD: the seeded quiet start (data in repo, analysis pending)

Overnight batch (2026-07-03), N = 100k quiet IC (400 rings × 250; runtime floor
|c1| ~ 1e-5 at t=0.2, ~15× below 25k), radial-displacement seeds
R → R(1+eps cos m phi), 4 T_unit each. Data: `sph/quiet_modes_m{M}_e{EPS}.dat`
(t, |c1..c4|, phase1, |c1| in three annuli). Stdout logs were lost to /tmp cleanup
(lesson: write logs into the repo or ~ next time).

Quick-look (|c1| at t=1/2/4, peak):

| run | t=1 | t=2 | t=4 | peak |
|---|---|---|---|---|
| control (m=0) | 1.9e-2 | 1.1e-2 | 1.7e-2 | 3.8e-2 |
| m=1, 0.25% | 2.4e-2 | 3.0e-2 | **1.40e-1** | 1.55e-1 |
| m=1, 0.5% | 1.3e-2 | 4.6e-3 | 1.29e-1 | 1.39e-1 |
| m=1, 1% | 2.2e-2 | 4.0e-2 | 9.9e-2 | 9.9e-2 |
| m=1, 2% | 7.9e-3 | 2.4e-2 | 7.4e-2 | 7.4e-2 |
| m=2, 1% | 1.8e-2 | 6.4e-2 | 1.59e-1 | 1.59e-1 |
| m=2, 2% | 1.1e-2 | 6.4e-2 | 8.5e-2 | 8.7e-2 |

Provisional reading (NOT yet analyzed properly): every seeded run — any amplitude
(≥0.25%), any m — pulls away from the control after t ≈ 1–2 and grows at ~0.5–1/T_unit
toward |c1| ~ 0.1; the response is m=1 regardless of the seeded m (|c2| peaks stay
≤1.5e-2); the ordering is INVERSE in eps (0.25% highest at t=4). So: a slow m=1
amplifier awakened by any finite disturbance; seed amplitude appears to set phase-mixing
history, not the outcome. Rate is ~10–20× the linear prediction, ~5–10× below the
Poisson-IC eruption. Two diagnostic caveats for the analysis: (a) the global |c_m| is
R-blind, so radial-displacement seeds register only after dynamics phase-mixes them
(first-output |c1| = 2.14e-3·eps, linear — one substep of response, not the seed);
(b) the control's own floor (1–4e-2) is COM-slosh contaminated — consider measuring
about the system COM and/or using the per-annulus columns and phase coherence.

**First tasks on re-up:**
1. Analyze the seven mode files properly: growth-rate fits per run in a common window,
   per-annulus amplitudes (cols 7–9), phase1 coherence/pattern speed, seeded-run minus
   control comparisons. Figure + verdict: is the 0.5–1/T_unit amplifier real, what is
   its Omega_p, and does it match anything (e.g. the razor-thin m=2 0.83/Tu? a slow
   eccentric mode the thick linear theory marginally stabilized?).
2. If real: rerun 1–2 cases with per-annulus PHASE output and longer baseline (8 Tu),
   and a repeat control with different ring-phase seed (floor statistics).
3. Then the remaining discriminators: pinned star (kills the reflex channel — cleanest
   mechanism test for the slow amplifier), Poisson-sampled SMOOTH family (noise
   amplitude vs grid-disk structure), thickness/active-fraction dials.
4. Independent thread: the external GPT-5.5/Codex review of the 2D code
   (`REVIEW_PROMPT.md`) was sent out ~2026-06-22 — fold in its findings when back.

## 3. Practical gotchas

- **Builds** (macOS, Apple clang): plain: `g++ -O2 -std=c++17 -Ilb94/sph <src> -o <bin>`;
  OpenMP: add `-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
  -L/opt/homebrew/opt/libomp/lib -lomp` (14 cores; `OMP_NUM_THREADS=14`). No `timeout`
  cmd on macOS. clang has no native `-fopenmp`.
- **Python venv**: /tmp is periodically cleaned — recreate:
  `python3 -m venv /tmp/lb94venv && /tmp/lb94venv/bin/pip install numpy matplotlib
  scipy pymupdf`. System python3 has no numpy. (Consider making a persistent venv in
  the repo next time.) `np.trapz` → `np.trapezoid` (numpy 2.x).
- **Long-run logs**: write into the repo or ~, not /tmp (we lost the seed-batch logs).
- **Papers** (gitignored, local): `paper/ARS89.pdf` (note in proof = p. 976 = PDF p.18),
  `paper/Toomre81.pdf` (ADS scan 1981seng.proc..111T), `1989ApJS...70..419H.pdf` (HK89),
  `LB94.pdf`, `YBL93.pdf`, `R85.pdf`, `BYRT90.pdf`, page scans in `lb94/paper/` etc.
  Read equations from page IMAGES (render via pymupdf), never the text layer.
- **Units**: SPH code units L=100 AU, V=1 km/s, G=1 → Mu=0.1127 M⊙, Tu(SPH)=474 yr.
  Linear code: G=M_*=1, R_u=100 AU → t_u=272.9 yr; conv to SPH T_unit = ×1.7369.
  ARS89 convention omega = m·Omega_p − i·gamma (growing ⇒ Im<0).
- **ICs**: `sph/quiet_ic_100k.dat` (400×250, seed 777) and `sph/ic_quiet.cpp` to
  regenerate; the Poisson IC generator is `sph/ic_disk.cpp` reading
  `ybl_grid_a0.010_40kyr.dat` (regenerate via `lb94/collapse_ybl 0.01` if the .dat
  files are missing — gitignored). `disk_ic.dat` currently = 25k viscous Poisson.
- **Key SPH knobs**: theta=0.4 for quiet work (0.7 production), N_neigh=50,
  BlockStepper dt0=0.05 Rmax=12. The star: `heavy=1` in IC files; quiet.cpp maps
  heavy → gas=false (collisionless); treesph.cpp keeps heavy as gas (LB94-faithful).

## 4. Key numbers (memorize-level)

- LB94 targets (from the page): viscous 0.602 M⊙ within 209 R⊙, J_c=2.31e52, β=0.439;
  Fig 1 inviscid core 0.42, viscous 0.60; Fig 2 at 20 kyr, T 205→33 K over 11–265 AU.
- Our 2D @ 40 kyr viscous: min Q 1.24–1.35 at ~30–56 AU; M_d/M_* = 1.10.
- SPH Poisson eruption: onset t≈0.4–0.7 Tu, effective gamma 3–6/Tu (up to ~15 locally),
  m=1 ≫ m=2; peak |c1| falls with N as seeds shrink (0.97→0.39, 6.25k→50k).
- Linear: m=1 ≤ 0.04/Tu (fragile), m=2 0.83/Tu (cavity artifact), thick: nothing.
- Quiet control: |c1| floor 1–4e-2 (25k: 2–7e-2) over 3–4 Tu, no growth.
- Seeded (pending analysis): all seeds → m=1 at ~0.5–1/Tu toward |c1| ~ 0.1 by t=4.
