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

## What the overplot shows

### Fig. 1 (core mass vs time)
1. **Timing lag.** LB94's core grows *gradually from ~3 kyr* and is largely in place by
   ~20 kyr (< 1 t_ff). Ours stays ~0 until ~20 kyr (≈ t_ff) then rises steeply. This is
   physical for a ρ∝1/r cloud (the dense inner cusp has a short local t_ff and should
   collapse first) — our **central IC softening** (ρ capped at r = ½ dx) plus the sink
   activation threshold removes that early inner collapse, delaying core onset.
2. **Over-accretion.** Our core reaches **0.68 M⊙** — well above LB94's inviscid 0.42 and
   even above the viscous 0.60. The **central sink is too aggressive** (density cap +
   drain), moving disk material to the core that LB94 retains.
3. **No viscous separation.** LB94 shows a clear inviscid/viscous split (0.42 vs 0.60);
   ours nearly coincide. With our over-accretion the inviscid core is already near the
   ceiling, leaving little room for the (weak, α=0.01) viscosity to add more.

### Fig. 2 (disk structure)
LB94's disk is formed by 20 kyr; ours forms ~30–40 kyr later (Fig. 1), so we overplot our
disk-formed epoch (60 kyr) and, as a thin dotted line, our 20 kyr state (still infalling).
1. **Temperature — the largest discrepancy.** LB94: a real **30–205 K gradient**
   (optically-thick inner disk heated by compression/accretion). Ours: **isothermal ~20 K**.
   This is the direct cost of the **gray cooling-relaxation** vs LB94's flux-limited
   radiation transport (audit issue #2) — the inner disk never heats up.
2. **Surface density.** Our **outer disk (~100–250 AU) matches LB94** (log Σ ≈ 1); our
   **inner disk (< 60 AU) is depleted** (log Σ dips to −1) — the same sink over-accretion as
   Fig. 1. LB94's Σ is higher and smoother overall (more disk mass retained).
3. **Specific angular momentum.** Reasonable agreement in the outer disk (both → log j ≈
   20.7–20.8, ~Keplerian); ours is modestly low in the inner disk. At 20 kyr our j is ~10¹⁸
   (un-spun-up infalling envelope) — the clearest signature of the timing lag.

## Bottom line
The **rotation/angular-momentum structure of the outer disk is reproduced well**; the two
real failures are (a) the **thermal structure** (no inner-disk heating — gray cooling), and
(b) the **central sink over-accreting** (too-massive core, depleted inner disk, washed-out
viscous separation, delayed onset compounded by the IC central softening). Both were the
top-ranked physics caveats in `AUDIT.md` (#2 cooling, #3 sink) — this comparison confirms
that prioritization and quantifies the gaps.
