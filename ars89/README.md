# ARS89 — m=1 SLING gravitational instability of a nearly-Keplerian disk

A clean, vectorized Python reproduction of the linear normal-mode eigenvalue code of

> **Adams, F. C., Ruden, S. P., & Shu, F. H. 1989, ApJ, 347, 959**,
> *"Eccentric Gravitational Instabilities in Nearly Keplerian Disks"*
> ([1989ApJ...347..959A](https://ui.adsabs.harvard.edu/abs/1989ApJ...347..959A)).

The code solves the `m = 1` ("eccentric") gravitational instability of a razor-thin,
nearly-Keplerian, self-gravitating disk around a central star — the **SLING** mechanism,
in which the reflex motion of the central star (the indirect, `m = 1` potential term)
feeds back on the disk and destabilizes a one-armed mode.

## Result

For the paper's canonical model (`p = 3/2`, `q = 1/2`, `R_D/R_* = 1e4`, `M_D = M_*`,
`Q_* = 10`) the code recovers the lowest-order growing mode:

| quantity | this code | ARS89 (Fig. 3) |
|---|---|---|
| eigenfrequency `ω` (units `Ω_D`) | `4.22 − 0.242 i` | `4.26 − 0.232 i` |
| pattern speed `Ω_p = Re ω` | `4.22` | `4.26` |
| growth rate `γ = −Im ω` | `0.242` | `0.232` |
| corotation radius `R_CR` | `0.453 R_D` | `0.452 R_D` |

Pattern speed agrees to **~1%**, growth rate to **~4%**, and the corotation radius to 0.2%.
The eigenfunction `S(r) = σ₁/σ₀` reproduces the Fig. 3 morphology: ~10 regular radial
oscillations across the inner disk with a 180° phase offset between `Re S` and `Im S` (the
alternating "bananas"), giving way to smoother structure peaking near the outer edge.

**Independent energy check (ARS89 eq 33).** Substituting the eigenfunction into the modal
energy budget — Reynolds stress + acoustic flux + direct/indirect gravitational work —
yields a growth rate `γ_energy` that matches the eigenvalue `γ` to **< 0.3 %** (converging
to < 0.05 % at high resolution). This reproduces ARS89's qualitative findings: the Reynolds
stress is negative through most of the disk (the disturbance feeds the shear), the net
acoustic flux nearly vanishes, and the indirect (SLING) term is a locally significant
contributor in the resonant cavity.

```
python -m ars89.validate            # solve, print ω + energy balance, write fig3/fig5 PNGs
python -m ars89.validate --N 2000   # finer grid
pytest -q                           # validation + unit tests
```

`Ω_D ≡ (G M_*/R_D³)^{1/2}` is the frequency unit; the code uses `G = M_* = R_D = 1`, so
eigenvalues come out directly in units of `Ω_D`.

## Method

The reproduction follows ARS89's own formulation and numerical recipe (their Appendix B).

**Governing equations.** Perturbations `∝ exp[i(mθ − ωt)]` with `ω = mΩ_p − iγ` (eq 6).
The linearized continuity + momentum equations (7a–c) are reduced, by eliminating the
velocity perturbations, to a single second-order ODE (eq 21):

```
L(h₁ + ψ₁ + ψ̃₁) + C h₁ = 0,     L ≡ d²/dr² + A d/dr + B        (eqs 8, 9)
```

with coefficients `A` (10a), `B` (10b), `C = −κ²(1−ν²)/a₀²` (10c) and `ν ≡ (ω−mΩ)/κ` (10d).
Here `h₁ = a₀² σ₁/σ₀` is the enthalpy perturbation, `ψ₁` the perturbed self-gravity, and
`ψ̃₁` the indirect (SLING) potential.

**Self-gravity** (razor-thin Poisson integral, eqs 7d/7e):
`ψ₁ = −2πG ∫ K_m(r,ρ) σ₁(ρ) dρ`, with
`K_m(x) = (x/π) ∫₀^π cos(mα)(1+x²−2x cosα)^{−1/2} dα`, `x = ρ/r`. Implemented in closed form
via complete elliptic integrals (`K₀(x) = (x/π)(2/(1+x))K(k)`,
`K₁(x) = (1/π)[(1+x²)/(1+x) K(k) − (1+x) E(k)]`, `k² = 4x/(1+x)²`). The integrable
logarithmic singularity at `x = 1` is handled exactly as in ARS89 (Appendix B): the two
cells straddling `ρ = r` are integrated by representing the smooth factor as a cubic and
integrating the singular kernel analytically (cubic product integration in `log x`).

**Energy analysis** (eq 33, `ars89/energy.py`): the velocity perturbations `u₁, v₁` are
reconstructed from the momentum equations (7b,7c), and the modal energy budget — Reynolds
stress, acoustic flux divergence, and direct/indirect gravitational work — is evaluated
against `2γ⟨ℰ⟩`. This provides an eigenfunction-level accuracy check independent of the
eigenvalue solve.

**Indirect / SLING term** (eqs 18b, 20b): `ψ̃₁(r) = ω² R₀ r` with
`R₀ = π/(M_*+M_D) ∫ ρ² σ₁ dρ` — a global, rank-1 coupling of every radius to the disk's
`m = 1` mass moment, i.e. the stellar reflex motion.

**Equilibrium** (Section III, Appendix A): power-law surface density
`σ₀ = σ_*(R_*/r)^p` (23a,b), sound speed `a₀ ∝ r^{−q/2}`, and the rotation curve
`Ω² = GM_*/r³ + (1/r)dψ_disk/dr + (1/r)dh₀/dr` (5a) built from Kepler (A2) + disk
self-gravity (A3, softened with `η = 0.1`, A5) + pressure (A4); epicyclic frequency
`κ² = (1/r³) d/dr[(r²Ω)²]` (5b).

**Boundary conditions:** inner `u₁ = 0` at `R_*` → Robin condition on `(ψ₁+h₁)` (eq 16);
outer Lagrangian `Δσ = 0` at `R_D` (eq 14a, with `𝒟 = κ² − (ω−mΩ)²`, 14b).

**Discretization** (Appendix B): geometric/log radial grid (`r_{j+1} = f r_j`, eq B1);
log-derivative finite-difference matrices `D¹, D²` (B4); the Poisson integral becomes a
dense matrix `J` (B2/27); the assembled complex matrix `M(ω)` (B6) has the boundary
conditions in rows 1 and `N`.

**Eigenvalue solver.** `M(ω)` is *non-linear* in `ω` (via `ν, ν²` in the coefficients and
`ω²` in the indirect term). The growing SLING mode is a genuine zero of `det M(ω)` but is
**not** a smallest-singular-value minimum — it is shielded by the dense set of neutral
real-axis cavity modes. Its existence near the ARS89 value is confirmed by an
argument-principle winding number of 1, and it is found by Newton iteration on the
determinant's logarithmic derivative,

```
ω ← ω − 1 / (d/dω log det M) = ω − 1 / tr(M⁻¹ dM/dω),
```

seeded at the ARS89 Fig. 3 value. The null right-singular vector at the solution is the
eigenfunction `S(r)`.

## Module layout

| file | role |
|---|---|
| `ars89/model.py` | `DiskModel` — parameters, unit system, equilibrium profiles (23, 24) |
| `ars89/discretize.py` | log-radial grid (B1) + `D¹, D²` log-derivative matrices (B4) |
| `ars89/kernel.py` | self-gravity kernel `K_m` (7e) → Poisson matrix `J` (27) + indirect operator (B5) |
| `ars89/physics.py` | rotation curve (A1–A5), `κ` (5b), `ω`-dependent coefficients `A,B,C` (10) |
| `ars89/eigensolve.py` | assemble `M(ω)` (B6) + BCs (14a,16); Newton det-root solver |
| `ars89/energy.py` | velocity reconstruction + modal energy budget (eq 33) |
| `ars89/validate.py` | run canonical case, compare to Fig. 3 & 5, plot `S(r)` and the budget |
| `tests/` | eigenvalue + energy-balance validation, kernel/derivative/equilibrium unit tests |

## Fidelity notes

* **Rotation-curve softening (`η`).** ARS89 use the *unsoftened* disk gravity (A2–A4) in the
  interior and soften (`η = 0.1`) only at the edges (A5). Their unsoftened interior curve is
  the faithful target; a direct comparison shows it matches the `η ≈ 0.04–0.05` member of the
  softened-everywhere family near the mode peak (where `Ω` is ~1% above the `η = 0.1` curve).
  The default `η = 0.05` therefore reproduces ARS89's interior rotation curve while remaining
  numerically robust at all resolutions. A fully unsoftened option
  (`build_equilibrium(..., unsoftened_rotation=True)`, with the A5 edge handling) is also
  provided but is more fragile near the disk edge. The eigenvalue tracks `η` smoothly:
  `η = 0.10 → 4.08 − 0.21 i`, `η = 0.05 → 4.22 − 0.24 i`, approaching the ARS89 value.
* **Poisson diagonal.** The singular `x = 1` cell uses ARS89's cubic-spline-over-singularity
  treatment (cubic product integration across the two adjacent cells).
* The `paper/` directory holds the scanned article and the page renders used to transcribe
  the equations.

## Reference

Adams, Ruden & Shu 1989, ApJ 347, 959. Companion theory: Shu, Tremaine, Adams & Ruden
1990, ApJ 358, 495 (1990ApJ...358..495S).
