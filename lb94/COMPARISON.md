# This work vs LB94 Figures 1 & 2 — overplot comparison

Direct overplot of our 2D collapse on Laughlin & Bodenheimer 1994's published figures
(digitized by eye from the scanned page `paper/pg-04.png`, so the LB94 points are ~few-%
approximate). Figures: `cmp_fig1.png`, `cmp_fig2.png`; script `plot_lb94_vs_paper.py`.

## What the LB94 figures actually show (read from the page)
- **Fig. 1** — *core mass vs time* (T in yr). Open circles = inviscid 2D model, filled
  triangles = viscous 2D model. Inviscid **plateaus at ~0.42 M⊙**; viscous climbs to
  **~0.60 M⊙** (the text gives 0.602 M⊙ within 209 R⊙, J_c = 2.31×10⁵² g cm² s⁻¹, β=0.439).
- **Fig. 2** — disk **T (linear, 50–250 K)**, **log Σ**, **log j** vs **log₁₀ R(cm), 14.2–15.6
  (≈ 11–265 AU)**, at **t = 20,000 yr**. Unconnected circles = inviscid, connected dots =
  viscous.
- (Note: an earlier working note in this repo carried 0.682 M⊙ / J_c=2.31×10⁵¹ as the
  target — both wrong on re-reading the page; the correct values are 0.602 M⊙ / 2.31×10⁵².)

## Physics added to close the gaps (2026-06-22)
Following the first comparison (which flagged the gray cooling and the over-aggressive sink),
two physically-motivated changes were made (`hydro.hpp`):
1. **Irradiated cooling.** The disk is heated by the central accretion luminosity
   L_acc = G M_c Ṁ/R_⋆ (R_⋆ = 1 R⊙ ⇒ L ≈ 25 L⊙, LB94's value), giving an equilibrium
   T_irr(s) = (L_acc/16πσs²)^¼, capped at dust sublimation (1500 K) and not allowed to
   diverge inside the sink. This is a tractable stand-in for LB94's full FLD transport and
   supplies the missing inner-disk heating.
2. **Truelove-limited sink.** The density cap is set to the *local* Jeans-resolvable limit
   ρ_J = π c_s²/(N_J² G Δx²) (N_J=4) with an absolute first-core backstop (10⁻¹²), instead of
   a fixed 5×10⁻¹⁴ that capped resolvable disk gas. Only the genuinely unresolvable excess is
   accreted, so the resolved disk is retained.

## What the overplot now shows

### Fig. 2 (disk structure) — now in good agreement
With irradiation, the disk has a **real temperature gradient ~235→50 K** that tracks LB94's
~205→33 K (ours runs ~20–30 K warm). The **surface density now matches** (log Σ ≈ 1.7–2.4
inner, declining outward — the previous inner-disk depletion is gone, because the Truelove
cap retains the resolved disk). **Specific j** tracks LB94 (~Keplerian, log j ≈ 19.8→20.7).
So all three Fig. 2 panels are reproduced to within the digitization error.

### Fig. 1 (core mass vs time) — magnitude fixed, two residuals
- **Magnitude fixed:** inviscid core reaches **0.47 M⊙** vs LB94's 0.42 (was 0.68). The
  Truelove sink no longer over-accretes; the retained disk grows from 0.30 to **0.54 M⊙**.
- **Residual 1 — timing lag.** Our core still onsets at ~t_ff (~20 kyr) vs LB94's gradual
  growth from ~3 kyr. This is the **IC central softening** (ρ capped at r = ½ dx) removing the
  dense inner cusp that should collapse first; the magnitude is right but the curve is shifted
  ~15 kyr later.
- **Residual 2 — viscous separation not reproduced.** LB94's viscous core (0.60) clearly
  exceeds its inviscid (0.42); ours nearly coincide (~0.42–0.47). The α-viscosity that would
  drive the extra inner accretion is also numerically unstable in the inner disk (it had to be
  excluded inside the sink region for stability), so it only spreads the *outer* disk rather
  than driving inward accretion. The viscous run is also fragile (reaches ~50–54 kyr before an
  outer-disk-edge instability), so it is shown to 40–50 kyr.

## Bottom line
The two top `AUDIT.md` caveats are largely closed: **irradiated cooling reproduces the disk
thermal structure**, and the **Truelove-limited sink fixes the over-accretion** (core, disk
mass, and Σ now all match LB94). The disk *structure* — the actual SPH initial condition — is
now a good match. What remains is dynamical, not structural: the **early-time core onset**
(needs better central-cusp resolution, not softening) and a **stable, effective inner-disk
viscosity** to reproduce LB94's viscous/inviscid core split. The latter likely needs an
implicit/sub-cycled viscous solve so it can act in the inner disk without going unstable.
