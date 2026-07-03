# Reading Toomre (1981) — and what he would have done with our 2D→3D handoff

Source: Toomre 1981, "What amplifies the spirals?", in *The Structure and Evolution of
Normal Galaxies*, p. 111 (ADS scan `1981seng.proc..111T`; local copy `paper/Toomre81.pdf`).
Context: GL's observation that the paper's power comes from *exquisite control over initial
conditions* — the opposite pole from the Poisson-sampled 25k-particle ICs of LB94.

## 1. The methodological anatomy of the paper

What actually makes it work, experiment by experiment:

1. **Smooth, analytically specified basic states.** Zang's V=const disk (with taper
   functions and a center cutout whose *index* is a dial), and the Gaussian disk of §4.
   Never a "realistic" messy equilibrium; always one whose linear theory is computable.
2. **Deterministic, coherent perturbations — never noise-as-signal.** Imposed transient
   tidal masses (Figs 1–2), an initialized packet of leading waves ("dust to ashes",
   Fig 8), a passing M51 companion (Figs 13–15). Amplitude is a dial (1%, 2%, 4% rows) so
   *linearity itself is tested* — he explicitly notes where the 4% case starts to cheat.
3. **Null and kinematic controls.** The same forcing with self-gravity off ("kinematic"
   rows in Figs 1–2) isolates exactly what self-gravity adds. Quiet-start N-body wherever
   particles are used (Zang & Hohl; the 20,000-body M51 runs), with "clock particles" as
   references.
4. **One-knob-at-a-time mechanism surgery.** The Fig 11 experiment: growth rates and
   pattern speeds of the Gaussian-disk modes A–F as a function of the *fraction f of the
   disk mass that is self-gravitating ("active")*, the rest frozen into the axisymmetric
   potential. And the feedback-loop ablations: insert an ILR barrier / make the center
   stiff or immobile — the bar modes die. Mechanism is identified by ablation, not
   attribution.
5. **Wave-anatomy diagnostics.** Corotation circles drawn on every mode; the 90°
   interference-node spacing used to decompose a standing pattern into trailing + leading
   traveling parts; growth rate and pattern speed always reported as functions of the dial.
6. **The local amplifier quantified.** Swing amplification in the shearing sheet, gain as
   a function of Q and X, checked three independent ways (GLB, LSK, Zang; his Fig 7).
   The global modes of §4 are then *explained* as feedback cycles of this local transient
   amplifier — the "Maclaurin analogy is fallacious" passage traces mode growth not to
   edges but to the loop through the center.

The paper's deep point for us: **growing "modes" are usually feedback cycles of a transient
amplifier, and the way to understand them is to measure the amplifier's gain and then cut
the loop in every place it could close.**

## 2. Would he have mapped the 2D viscous IC into particles? No — not first.

The Monte-Carlo mapping (Poisson-sample the grid, let shot noise seed whatever grows) is
the exact opposite of every design choice above: the seed is uncontrolled, broadband,
resolution-dependent (∝1/√N — as our audit measured), and entangled with the start-up
transient of an IC that is not in discrete equilibrium. Toomre would have treated the 2D
multi-grid output as *defining an axisymmetric basic state*, then studied that state:

**Step 0 — idealize the basic state.** Fit Σ(R), c_s(R), Ω(R), M_d/M_c from the 40 kyr
viscous disk with a smooth analytic family. The 2D output picks the family member; the
family (not the noisy snapshot) is the object of study. Knobs it exposes: taper shape,
inner cutout, Q(R) normalization, disk/star mass ratio.

**Step 1 — linear theory of that state.** We already own the machinery (the ARS89
eigenvalue code from the first leg of the trilogy): recompute the m=1 and m=2 global modes
of *our* disk profile with a smooth taper (no truncated edge → no c-ambiguity), the
indirect term included, softening ↔ finite thickness. Deliverables: γ_m, Ω_p,
eigenfunctions — predictions, in advance, for the particle experiments.

**Step 2 — a quiet-start 3D realization that flatlines.** Rings of equally spaced
particles in azimuth (stratified/deterministic sampling of the same Σ, not Poisson), each
ring's phase chosen to cancel low-m to machine precision; vertical structure hydrostatic
*in the SPH's own discrete forces* (relax axisymmetrically first — this also eliminates
the inner-hole refill transient we identified). The null experiment: verify |c_m| stays
at ~10⁻⁴ for several T_unit. Until the null runs clean, no instability claim is
meaningful — this is his kinematic-row discipline.

**Step 3 — controlled seeding, one mode at a time.**
- Impose a coherent m=1 at 0.5%, 1%, 2% → γ₁, eigenfunction, and the amplitude-scaling
  linearity check.
- Same for m=2 → γ₂. (Swing predicts m=2 favored, X₂≈1.5–2 vs X₁≈3–4; the eccentric-mode
  picture predicts m=1 wins. Clean discrimination.)
- Launch a *leading* wave packet and measure the swing gain directly against the X,Q
  prediction — the Fig 8 experiment in our disk.

**Step 4 — mechanism surgery (the loop-cutting).** For m=1 in a near-Keplerian disk there
is effectively no ILR (κ≈Ω ⇒ Ω−κ/1 ≈ 0), so the feedback loop must close through either
the *center* (the star's reflex motion — the indirect term) or the *outer turning region*.
Hence:
- **Pin the star** (or blend: star free ↔ star on rails). A true eccentric mode dies or
  slows drastically; a local GI/swing response doesn't care. This is the single most
  decisive experiment — the N-body version of ARS89's c-dial, with no bookkeeping
  ambiguity.
- **Active-fraction dial** (Toomre Fig 11): gravity g = f·g_live + (1−f)·g_frozen-axisym;
  measure γ₁(f).
- **Inner stiffness barrier**: heat the inner disk (raise Q inside ~30 AU) to block wave
  transmission through the center region; his bar-stabilization trick, transposed.
- **Outer control**: move/soften the 230 AU absorbing ring (anti-edge-feedback check).

**Step 5 — only now, the mapped IC.** Re-run the Poisson-mapped IC as a *validation* that
the messy realization grows the same eigenmode (same γ, Ω_p, eigenfunction) from its
uncontrolled seed — demoting the 1994-style IC from "the experiment" to "the consistency
check." Our onset-scaling result (delay = ln√(N/N₀)/γ) already shows this demotion is
justified: the noise only sets the *starting amplitude*.

**On "3D" specifically:** he'd treat dimensionality as a dial too. Finite thickness is,
in his own framing, a softening that *reduces* the swing gain — so: same basic state run
(a) razor-thin-2D-like (small h/eps), (b) full 3D thickness; compare γ₁. Plus 3D adds
genuinely new channels worth separating: bending/buckling (out-of-plane m=1) and vertical
resonances — measurable only because the in-plane experiment is under control first.

## 3. Diagnostics to build (all cheap)

- c_m(R, t) **amplitude and phase** at fine cadence (~0.02 T_unit) → rigid-pattern test,
  corotation radius, radial eigenfunction, trailing/leading decomposition (interference
  node spacing à la his §4).
- Growth rate and pattern speed *as functions of the dial* (f, star mobility, Q(R),
  thickness, taper) — every result a curve, not a single run.
- Clock particles / reference circles in all visualizations.

## 3.5 Step 0 done — the pinned basic state (fit_basic_state.py, fig_basic_state.png)

    Sigma(R) = 48.4 (R/100AU)^(-2.51) exp[ -(54.2 AU / R) - (R / 267 AU)^2 ]  g/cm^2
    c_s(R)   = 0.534 (R/100AU)^(-1/4)  km/s        (q = 1/2 fixed; T(100 AU) = 81 K)
    M_star   = 0.340 Msun,  M_disk = 0.374 Msun  ->  M_d/M_star = 1.10

- Fit window 18–280 AU (inside 18 AU is the sink cliff, whose role the SPH heavy particle
  takes over). rms(ln Sigma) = 12% — this is the disk's real ring/bump structure, and it is
  a *floor*: all 16 integer taper-index pairs (n,m) in {1..4}^2 fit within 0.116–0.156, so
  **the data cannot select the taper shape**. Canonical choice (n=1, m=2) by
  interpretability: p matches the directly measured 60–180 AU slope (2.57), R_out matches
  the visible outer edge, Gaussian outer taper. Taper indices are therefore a *dial* for
  the stability analysis, not a measurement.
- **Predictions (not fitted):** self-consistent rotation curve (star + softened razor-thin
  disk quadrature + midplane pressure term) matches the grid v_phi to **6.4% rms** with
  thickness-mimicking softening eps = H/2 (10% if razor-thin — the thickness correction is
  real at H/R ~ 0.25); Q(R) to **7.8% rms**, Q_min = 1.35 (grid: 1.28).
- Honest wart: the family puts Q_min at 56 AU; the grid's is at 32 AU, sitting on the
  19–31 AU plateau/cliff that the smooth family deliberately rounds. That plateau edge is
  the first candidate "feature" to reintroduce as a parameterized groove later (§2 item 5).
  Validated domain of the family: 30–220 AU.

## 4. Suggested order of attack

1. Quiet-start machinery + null run (Step 2) — everything else depends on it.
2. Seeded m=1 vs m=2 growth rates (Step 3) — one clean discriminating number pair.
3. Pinned vs free star (Step 4) — the mechanism verdict.
4. Linear eigenproblem of the basic state (Step 1, in parallel — reuses ars89 code).
5. Active-fraction and thickness dials; then the mapped-IC validation last.
