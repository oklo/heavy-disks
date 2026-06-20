# LKA98 — transcribed & verified equations (Laughlin, Korchagin & Adams 1998, ApJ 504, 945)

Source: `paper/Laughlin_1998_ApJ_504_945.pdf` (clean text layer + page renders `paper/page-NN.png`).
Convention (eq 21): perturbation `f = f0(r) + f1(r) exp[i(ωt − mφ)]`, with (eq 22)
`ω = m Ωp − iγ1` (Ωp = pattern speed, γ1 = growth rate, m = azimuthal number).
Reference model is m = 2.

## Governing equations (§2)
- Continuity / radial mom / azimuthal mom in polar coords (1)–(3).
- Softened Poisson (4): `Ψ = −G ∫_Rin^RD Σ(r')r'dr' ∫_0^2π dφ'/√(r²+r'²−2rr'cosφ'+g²(r))`.
- Polytropic enthalpy (5): `h = [γp/(γp−1)] K Σ^{γp−1}`; sound speed (12): `cs² = K γp Σ0^{γp−1}`.
- Index relation (13): `γp = 3 − 2/γ`; reference model uses `γ = γp = 2` (stiff EOS).

## Reference ("standard") disk (§3, eq 14)
- `Σ0(r) = S0 exp[−(r−R0)²/w²]`, R0 = 0.45 (peak), Rin = 0.25, RD = 1.0.
- S0 = 0.372 ⇒ disk mass mD = 0.4, central star m* = 0.6 (Mtot = 1). G = 1.
- K = 0.25, w² = 0.05.  Toomre Q = cs κ/(πGΣ0), min Q = 1.27 at r = 0.504 RD.
- Heavy disk: mD/m* = 0.67 (vs ARS89 nearly-Keplerian).

## Linear eigenvalue problem (§4) — structurally identical to ARS89 (7a–c), no SLING term
Let `A ≡ γ1 + i m Ωp − i m Ω = i(ω − mΩ)`, `ν ≡ (ω − mΩ)/κ` (so `A = iκν`),
`W ≡ (cs²/Σ0) σ1 + Ψ1` (enthalpy + potential perturbation).
- (23) continuity: `A σ1 + (1/r) d/dr(r Σ0 u1) − (im/r) Σ0 v1 = 0`
- (24) radial mom: `A u1 − 2Ω v1 + dW/dr = 0`
- (25) azimuthal mom: `A v1 + (κ²/2Ω) u1 − (im/r) W = 0`

### Direct-integration (shooting) ODEs (26),(27) — INDEPENDENTLY RE-DERIVED, both correct
- (26): `dσ1/dr = [2ΩΣ0 m/(cs² ν κ r)] W − [A Σ0 (ν²−1)/(cs² ν²)] u1
                 − (Σ0/cs²) d/dr(cs²/Σ0) σ1 − (Σ0/cs²) dΨ1/dr`
- (27): `du1/dr = (im/(νr))[(m/(rκ)) W + (A/(2νΩ)) u1] − (A/Σ0) σ1 − (1/(Σ0 r)) d/dr(rΣ0) u1`
Verified by eliminating v1 from (23)–(25) and using A = iκν; every term matches.

### Boundary conditions
Reflective at both edges: `u1 = 0 at r = Rin` and `u1 = 0 at r = RD`.

### Solution methods (paper compares both)
1. Matrix method (à la ARS): reduce to one integro-diff eq in σ1, discretize, solve det = 0.
2. Direct integration: RK from inner BC (u1 = 0, σ1 arbitrary), 2-D Newton–Raphson on
   (γ1, Ωp) to hit u1(RD) = 0, iterate the potential Ψ1 to convergence.

## Equilibrium rotation curve (Appendix C)
- (C1): `r Ω²(r) = ∂/∂r(Ψ + Ψ* + h0)`  (radial balance, u = 0).
- (C2): stellar + pressure contributions; Ψ* = −G m*/r ⇒ ∂Ψ*/∂r = G m*/r².
- (C3): `dh0/dr = −(2 S0 K γp (r−R0)/w²) exp[−(γp−1)(r−R0)²/w²]` (written for γp = 2).
- (C4): disk self-gravity contribution `r ∂Ψ0/∂r` from the softened kernel.

## Softening & self-gravity kernel (Appendix C)
- Softened axisymmetric kernel (C5): `K(r',r) = (r'/πr) ∫_0^π dφ/√(1+(r'/r)²−2(r'/r)cosφ+η²)`,
  with η(r) = g(r)/r dimensionless.
- Softening (C6): `η(r)² = 0.01·[(r−Rin)/(RD−Rin)]^6 + 1e−4`  (η≈0.01 interior → ≈0.1 outer edge).
- Potential-derivative via elliptic integrals (C7), Binney & Tremaine (1987) eq (2-146):
  `r ∂Ψ0/∂r = (1/√r) ∫_Rin^RD [F(π/2,ξ) − (1/4)(ξ²/(1−ξ²))(r'/r − r/r' + η²/(rr')) E(π/2,ξ)] ξ Σ0(r') √r' dr'`.
- Perturbed potential Ψ1 (linear): m-th harmonic of the softened Green's function — the
  softened analogue of ARS89's K_m (i.e. cos(mφ) in the numerator, +η² in the denominator).

### ✗ CONFIRMED TYPO — (C8) elliptic modulus
Printed: `ξ² = 4rr'/(r + r'² + η²(r))`.  Dimensionally inconsistent (r + r'²), missing the
2rr' cross-term. The correct (self-derived, = BT 2-146) modulus is
    `ξ² = 4rr'/[(r+r')² + η²(r) r²] = 4rr'/[r² + 2rr' + r'² + η²(r) r²]`.
VERIFIED (2026-06-19): with the correct modulus the equilibrium reproduces the paper's own
stated fact Q_min = 1.273 at r = 0.508 (paper: 1.27 at 0.504) and disk mass = 0.400; the
literal (C8) gives ξ² > 1 → NaN (no valid equilibrium). So (C8) is a typeset corruption of
`(r+r')²` plus a dropped `r²` on the softening term.  Implemented in lka98/selfgravity.py
(modulus="correct" vs "paper_c8").

## LINEAR RESULT (this code, matrix method, 2026-06-19) — DISCREPANCY TO RESOLVE
With the (C6) Appendix-C softening, my solver finds the dominant growing m=2 mode at
**Ωp = 1.573, γ1 = 0.724** (ω = 3.146 − 0.724i), N-converged (N=400→700 identical) and
robust to softening scale (0.1×–3× → Ωp 1.59→1.56). Found via argument-principle contour
localization of the det-zero (Newton alone drifts to the neutral real-axis modes, as in ARS89).

This matches the paper's *FFT-softening / simulation* value (Ωp=1.59, γ1=0.74) to ~1–2%,
but sits ~8%/10% BELOW the paper's quoted *Appendix-C* matrix value (Ωp=1.70, γ1=0.804).
Softening magnitude does NOT bridge the gap (even near-unsoftened gives ~1.60). All internal
checks pass (equilibrium Q_min, kernel vs quad to 1e-14, rotation curve cross-check to 0.00%,
m=0 kernel ↔ equilibrium consistency). Open question: does the page-7 "1.70/0.804 with the
Appendix-C softening" actually use lighter softening than (C6) — i.e. is (C6) [described as the
"sizable" softening for the 3rd-order calc] the same softening that gave 1.70? My (C6) result
coincides with the sim value, suggesting (C6) → ~1.59, not 1.70. To adjudicate: (a) Appendix-A
energy check (self-consistency), (b) Athena hydro growth rate (independent ground truth; paper's
sim gives γ≈0.74, matching my 0.72).

## SOFTENING SWEEP (2026-06-19) — cannot recover 1.70/0.804
Sweeping a constant softening eps in BOTH equilibrium and perturbation, contour-localized,
N-converged: eps=3e-3 -> Op=1.48/g=0.73 ; eps=1e-4 -> 1.57/0.77 ; eps=3e-5 -> 1.583/0.762.
Monotonic, asymptotes at Op~1.58, g~0.76 as eps->0 (the unsoftened "proper thin-disk" limit),
with Q_min->1.27. The C6 softening gives 1.573/0.724 (= the FFT/simulation value 1.59/0.74).
=> Reducing softening does NOT reach the paper's quoted matrix value 1.70/0.804; my formulation
caps ~6-7% below it. Both my matrix result and the paper's OWN simulation (FFT-softened, gamma~0.74)
sit at ~1.58/0.75. The quoted 1.70/0.804 is the outlier; to adjudicate, implement the independent
direct-integration method (paper compares matrix vs direct) and run the Athena hydro (ground truth).

## DIRECT-INTEGRATION CROSS-CHECK (2026-06-19) — confirms 1.574, not 1.70
Implemented the shooting method (lka98/directint.py): RK45 integration of (26,27) from Rin,
inner Newton on omega for u1(RD)=0 at fixed Psi1, outer loop revising Psi1=P sigma1 (seeded
from the matrix eigenfunction, per the paper's "refine a known mode" usage). It shares neither
the derivative operator nor the BC rows with the matrix solver.
RESULT: shooting -> Omega_p=1.5745, gamma1=0.7239; matrix -> 1.5744, 0.7235 (agree to 4 sig figs).
Crucially, seeding the shooting AT the paper value (omega=3.40-0.804i) flows DOWNHILL to 1.5745 —
there is no self-consistent mode at 1.70/0.804. Two independent methods + winding number (one root)
+ softening sweep (asymptote 1.58) all converge on Omega_p~1.57-1.58, gamma1~0.72-0.76, = the paper's
OWN simulation value (1.59/0.74). The quoted matrix/direct value 1.70/0.804 is not reproducible;
most consistent with a numerical issue in the original 1998 code or an undisclosed config parameter.

## ENERGY CHECK (Appendix A, eq A1) — PASSES (2026-06-19)
(A1) [Broadbent & Moore 1979], for m>1: gamma sigma0 <u1^2+v1^2+h1^2/cs^2> = -r(dOmega/dr)sigma0<u1 v1>
 - (1/r)d/dr(r sigma0<h1 u1>) - sigma0<u1.grad Psi1>. Single gamma is correct (the 1/2 in the energy
density and the 2 in 2*gamma cancel). Implemented in lka98/energy.py (same azimuthal-conjugation care
as ARS89: the (im/r) factor in the gravity-work term is NOT conjugated, only Psi1_hat).
RESULT: gamma_energy = gamma1 to 0.05% at N=800 (0.22% at N=400, converging). Gravity work dominates the
budget (gravitational instability), Reynolds smaller, acoustic ~0 — physically sensible. THIRD independent
confirmation of Omega_p~1.574, gamma1~0.724 (after matrix == shooting). No typo found in (A1).

## REFERENCE-CODE LINK (ucolick ~laugh/oxide/codes) — note
The page has Fortran codes (fewbody, fluid.iso/fluid 1D Lagrangian hydro, von_3daniken 3D Eulerian hydro,
integrator) but NOT a 2D disk code with van Leer advection + Kalnajs FFT. The original V2D.F is not posted
there; a van-Leer/Kalnajs C++ reimplementation would be from-scratch (standard methods).

## VALIDATION TARGETS (dominant m=2 spiral, standard disk)
- Matrix method:        Ωp = 1.70,  γ1 = 0.804
- Direct integration:   Ωp = 1.71,  γ1 = 0.794   (Appendix-C softening)
- FFT-mimicking softening (matches the hydro sim): Ωp = 1.59, γ1 = 0.74  (from LKA / Fig 1)

## Still to transcribe (for later tasks)
- Appendix A: energy analysis (cross-check on γ1, like ARS89 eq 33).
- §5: first-order mode-mode coupling (2nd-order transport).
- §6: third-order saturation + higher-order energy constraint.
- Appendix B: boundary conditions for the nonlinear problem (eqs B3–B14).

## Diagnostics from the simulations (§3.1)
- Fourier amplitude (15): `|a_m(r,t)| = |(1/2π)∫ Σ e^{−imφ} dφ|`; phase (16); pattern speed (17).
- Global amplitude (18): `|C_m| = |(1/mD)∫∫ Σ(r,φ) r dr e^{−imφ} dφ|` (note 1/mD normalization).
- m=0 deformation (19); global growth rate (20): `γ_m = d/dt ln|C_m|`.
- Paper's hydro: "V2D.F" (Laughlin 1994), 2nd-order van Leer + FFT Poisson, 256×256. We use Athena++.
