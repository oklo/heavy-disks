# Linear global modes of the LB94 basic state — predictions for the SPH experiments

Machinery: `lineig.py` — the ars89 eigenvalue code generalized from pure power-law Sigma to
the pinned tapered family (V-bracket power-law constants promoted to local log-derivative
fields; sigma0-ratio weights in the Poisson and indirect matrices; S-bracket unchanged since
c_s is an exact q=1/2 power law). Checks: `lineig_check.py`, figure `fig_lineig.png`.

**Regression:** with tapers off and matched equilibrium softening, the generalized matrix
equals ars89's to 1.8e-15 and the canonical eigenvalue agrees to 3e-13. (The first apparent
failure was a softening-default mismatch, eta_eq 0.05 vs 0.1 — itself a reminder of how
sensitive these m=1 growth rates are to small equilibrium changes.)

Units: G = M_* = 1, R_u = 100 AU, t_u = 272.9 yr; conversions quoted per SPH T_unit
(474 yr), conv = 1.737. ARS89 convention omega = m Omega_p - i gamma.

## Results (razor-thin, eta_eq = 0.1, N = 700; domain 4–350 AU)

| mode | omega (code) | gamma [/T_unit] | Omega_p [/T_unit] | R_cor |
|---|---|---|---|---|
| m=1 (leading) | 1.899 − 0.022i | **0.04** | 3.30 | 75 AU |
| m=1 (next three) | 2.98/3.51/4.03 − 0.01–0.02i | 0.02–0.03 | 5.2–7.0 | 42–53 AU |
| m=2 (only root) | 6.792 − 0.479i | **0.83** | 5.90 | 48 AU |

**Argument-principle root count** (exhaustive within the window): winding number over
Re omega ∈ [0.15, 4] (Omega_p down to corotation ≈ 40 AU), gamma_code ∈ [0.3, 3]
(i.e. gamma = 0.5–5 per T_unit): **0.000 — no fast m=1 root exists.** The sigma_min hunt
missed nothing.

**Sensitivity of the found modes:**

| variation | m=1 leading | m=2 |
|---|---|---|
| N 700 → 1000 | unchanged | unchanged |
| r_in 4 → 6 AU | unchanged | **gamma → 0** |
| r_out 350 → 450 AU | **gamma → 0** | unchanged |
| thickness (softened kernel, eta = 0.25) | **gamma → 0** | **gamma → 0** |

And the eigenfunctions (fig_lineig.png) say why: the slow m=1 modes are oscillatory wave
packets in the **outer taper** (amplitude peaking at 300+ AU); the m=2 mode is localized at
the **inner boundary**. Both families are boundary/taper residents.

## The finding

**The smooth LB94 basic state supports no robust growing linear mode, m=1 or m=2.** The
razor-thin residual growth (gamma ≤ 0.04/T_unit for m=1; 0.83/T_unit for one inner-cavity
m=2) is boundary-fed and vanishes with finite thickness or domain changes — the same kind
of edge-bookkeeping fragility the ARS89 note in proof exposed, here diagnosed by ablation.
The SPH-measured m=1 growth (gamma ≈ 3–6 per T_unit, e-fold ~120 yr) exceeds anything the
basic state's linear theory supports by a factor of ~100, and no linear mode matches the
SPH's measured pattern speed (1.6–2.2 /T_unit, corotation 90–120 AU) either.

This **revises the M1_ANALYSIS.md conclusion** ("a genuine linear instability of the disk"):
the onset-scaling evidence there showed amplitude-proportional amplification of sqrt(N)-
scaled noise with a common effective rate — but the smooth disk provides no eigenmode to do
it. The live interpretations are now:

1. **Toomre's own §5 picture: recurrent swing amplification of fresh particle noise.** Not
   an eigenmode; a noise-*collective* amplifier that still scales with seed amplitude
   (consistent with the onset scaling) but has no smooth-disk counterpart. At our X_1 ≈ 3–4
   the single-pass m=1 swing gain is modest, but recycling + finite-amplitude noise
   (density fluctuations ~ 1/sqrt(N_neigh) ≈ 14% locally!) could plausibly sustain it.
2. **Realization artifacts**: the inner-hole refill transient; the heavy central particle
   — which is a *gas* particle participating in SPH pressure forces with 0.35 Msun on one
   smoothing kernel — stirring the inner disk.
3. A physics channel absent from the linear model: the Lagrangian-isothermal SPH EOS
   (particles carry their birth c_s) vs the Eulerian-isothermal linear theory; or genuine
   3D (bending) coupling. Noted, but neither plausibly supplies a factor ~100 — and
   thickness corrections push the *wrong way* (they stabilize).

## Predictions for the SPH program (falsifiable, in advance)

- **Quiet start, unseeded:** |c_m| flatlines for many T_unit. Any visible growth within
  ~5 T_unit contradicts the smooth-disk linear theory.
- **Seeded m=1 (0.5–2%):** no exponential growth; at most neutral pattern rotation near
  Omega_p ≈ 3.3/T_unit. Growth at ~4/T_unit would falsify the linear result outright.
- **Seeded m=2:** for the thick disk, nothing robust (the razor-thin 0.83/T_unit cavity
  mode requires the disk to extend sharply to <5 AU, which the real system does not).
- **Poisson-mapped IC (the 1994-style run):** if the violent m=1 persists there while the
  quiet start flatlines, the 1994 SPH instability was a property of the noisy realization
  — seeded and sustained by particle noise — not of the underlying disk.

Caveats kept in view: the linear model is 2D (razor-thin or softened), Eulerian-locally-
isothermal, point-mass star + exact indirect term; the winding window covered Omega_p up to
4 code (corotation ≥ 40 AU) and gamma up to 5/T_unit — the SPH-measured values sit well
inside. Extending the count window is cheap diligence if wanted.
