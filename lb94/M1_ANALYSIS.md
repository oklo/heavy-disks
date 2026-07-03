# Is the m=1 real? — the ARS89 note in proof, and the resolution progression re-examined

Two questions (GL, 2026-07-03): (1) does the ARS89 note added in proof invalidate the SLING
headline result we had been citing as the mechanism? (2) In the resolution audit, is the
falling peak-|c1|-at-fixed-time a *seed* effect (a genuine linear mode started from smaller
Poisson noise) or a start-up transient? Analysis: `sph/analyze_m1.py` → `sph/fig_m1_analysis.png`.

## 1. The ARS89 note added in proof (ApJ 347, 959, p. 976)

Summary of the note: a sharply truncated power-law disk contributes extra *edge* terms to
both the direct potential (their eq. N1) and the indirect potential (eq. N2), because an
elliptical distortion of the edge carries zeroth-order surface density into new area. N1
formally diverges at the edge (needs softening); the edge term in the *indirect* potential is
multiplied by an uncertainty factor c. Their reruns:

- **c = 0** (edge terms off): recovers the paper's results — except that small softening
  introduces a new growing "edge mode" (m=1 analog of the Toomre-1981 Maclaurin m=2 case).
- **c = 1** (full correction): the paper's half-integer m=1 modes — the headline SLING
  modes, including the canonical eigenvalue we reproduced — **"become completely
  stabilized,"** replaced by only *weakly* growing whole-integer modes.
- **c = 2** (arbitrary): whole-integer modes grow at rates comparable to the original ones.

Reading: **the specific published SLING growth rates are not trustworthy** — the instability
they computed lives or dies by the bookkeeping of a sharp outer edge in the indirect
potential, and the physically complete treatment (c=1) kills it. Their closing sentence
retreats to the hope that eccentric disturbances survive "under a sufficiently wide variety
of circumstances" — a hope, not a result.

Implications for us:
- Stop attributing the SPH m=1 to "the ARS89/SLING mechanism" (as earlier commits did).
  SLING requires the sharp-edge feedback loop; our disk tapers smoothly, so that mechanism
  isn't even available in our system.
- Structural advantage of the SPH realization worth noting: **the entire c-ambiguity is
  moot in the particle simulation** — the indirect potential is included automatically and
  exactly (the central heavy particle moves freely under momentum-conserving forces), and
  there is no truncated edge. Whatever grows in the SPH does not depend on edge bookkeeping.

## 2. The resolution progression: seed effect, not stabilization

### Noise floor = Poisson expectation
The pre-growth |c1| (t ≤ 0.3) matches E|c1| = sqrt(pi/4N) at every N:

| N | measured floor | sqrt(pi/4N) |
|---|---|---|
| 6,250 | 0.0148 | 0.0112 |
| 12,500 | 0.0084 | 0.0079 |
| 25,000 | 0.0064 | 0.0056 |
| 50,000 | 0.0039 | 0.0040 |

The seed is pure particle noise, scaling as 1/sqrt(N).

### Onset delay matches the common-growth-rate prediction
If a genuine mode with N-independent gamma grows from a 1/sqrt(N) seed, the time to reach a
fixed amplitude shifts by ln(sqrt(N/N0))/gamma. With the fitted mean gamma ≈ 4.4/T_unit:

| N0 → N | predicted delay | observed |
|---|---|---|
| 6.25k → 25k | +0.16 | +0.18 |
| 6.25k → 50k | +0.23 | +0.22 |

(The 6.25k→12.5k pair is within the 0.1-T_unit sampling noise.) Fitted gamma: 5.7 (25k),
3.2 (50k) — scatter from sparse sampling, no systematic N-trend.

**Conclusion: the falling peak-|c1|-at-fixed-time in the audit was a seed-amplitude effect,
not convergence toward stability.** The earlier audit commit over-claimed ("converges toward
marginally-stable") — corrected here. What *is* genuinely resolution-dependent is the
nonlinear outcome: at N ≤ 12.5k the disk disrupts with a broad spectrum (m=1–4 all ~0.9,
noise-triggered clumping); at 50k it remains a cleaner m=1-dominated spiral.

### So: was the handoff disk actually m=1 unstable?
As realized in this SPH system — **yes, genuinely dynamically unstable**: sustained
exponential growth at gamma ≈ 3–6 /T_unit (e-fold ≈ 100–150 yr) over ≥4 e-folds, long after
start-up settles, with the noise-seed scaling a true linear instability predicts. Caveats
that keep "as realized" doing real work:
- **Start-up transient exists but is localized**: the inner hole (cells subsumed into the
  heavy particle) refills by t ≈ 0.2 (particle count at R≈13 AU: 31 → 1921); the bulk disk
  rearranges only 2–6% (median). A transient *trigger* can't be fully excluded — but a
  one-shot amplifier would give a plateau, not sustained exponential growth, and a
  systematic (N-independent) trigger would break the observed onset scaling.
- The mode emerges **inside-out** (inner annulus leads at t=0.4; by t=1.0 the amplitude
  peaks at 120–180 AU) — natural for shortest-dynamical-time-first, but it means the linear
  eigenfunction and any inner-edge influence need finer diagnostics.
- Pattern speeds (0.1-T_unit snapshot cadence, noisy): mid/outer annuli give
  Omega_p ≈ 2.2/1.6 rad per T_unit → corotation ≈ 90–120 AU, and the outer pattern rotates
  *faster than the local material* — wave-like, not material winding. Inner annulus
  measurement unreliable (saturates early). Suggestive of a coherent global pattern; not yet
  a clean demonstration.

## 3. Swing amplification: insufficient for m=1 (confirming the 1994 analysis)

Toomre's swing parameter X_m = kappa^2 R/(2 pi G Sigma m) for the viscous IC:
X(m=1) ≈ 2.8–4.2 across 25–220 AU (marginal-to-weak; swing needs X ≲ 3), while
X(m=2) ≈ 1.4–2.1 (squarely in the strong-swing regime), with Q ≈ 1.4–1.8. **If the growth
were swing amplification of noise, m=2 should dominate. We observe m=1 ≫ m=2 in the linear
phase** — so swing is insufficient for the m=1, exactly as the 1994 analysis concluded.

The natural surviving interpretation: a **global eccentric (m=1) instability of a
comparable-mass disk–star system** (M_disk/M_c ≈ 1 here), in which the star's reflex motion
(the indirect term — automatic in the SPH) participates. This is the class of mechanism the
ARS89 note's closing sentence hoped would persist — but established here without any sharp
edge, and without their c-ambiguity.

## 4. Candidate next experiments (for discussion)

1. **Quiet start**: azimuthally phase-cancelled particle placement to suppress the Poisson
   m=1 seed by ~10–100x. Prediction if the mode is a real linear instability: same gamma,
   onset delayed by ln(seed ratio)/gamma. A start-up-transient driver would break this.
2. **Controlled seeding**: inject a coherent m=1 (then, separately, m=2) at fixed small
   amplitude on a quiet start → clean gamma_1, gamma_2 and the eigenfunction; directly tests
   the swing prediction (gamma_2 > gamma_1?) against the eccentric-mode prediction.
3. **Pin the central particle** (suppress its motion): removes the indirect/reflex channel.
   A true eccentric mode should slow or die; local GI/swing would not care. This is the
   N-body analog of ARS89's c-experiment, minus the bookkeeping ambiguity.
4. **Inner-hole control**: gentler core–disk interface (smaller subsumed region, or an
   axisymmetrized relaxation phase before releasing m>0) to remove the refill transient.
5. **Fine-cadence mode output**: c_m(R, t) amplitude *and phase* every ~0.02 T_unit →
   rigid-pattern test, corotation radius, radial eigenfunction.
6. **Fixed-amplitude resolution comparison + seed ensembles; extend 50k (add 100k) into
   saturation** — the remaining "was the 1994 outcome resolution-limited" endgame, done at
   matched amplitude rather than matched time.
7. **Outer-boundary control**: move/soften the 230 AU absorbing ring (anti-SLING check).
