# LKA98 weakly-nonlinear theory — algebra check

## Setup (gamma_p = 2 simplification)
For the reference model gamma_p=2, the enthalpy is EXACTLY linear in sigma:
h = K gamma_p (Sigma0+sigma)^{gamma_p-1} - ... = (cs^2/Sigma0) sigma, with cs^2=K gamma_p Sigma0.
So the pressure term (cs^2/(gamma_p-1))[(1+sigma/Sigma0)^{gamma_p-1}-1] = cs^2 sigma/Sigma0 = h~,
i.e. NO pressure nonlinearity. The only quadratic nonlinearities are advective + geometric:
  (9)  radial:    v~^2/r - u~ d_r u~ - (v~/r) d_phi u~        [v~^2/r is the centrifugal/geometric]
  (10) azimuthal: -v~ u~/r - u~ d_r v~ - (v~/r) d_phi v~       [-v~u~/r is the curvature term — CONFIRMED present]
  (11) continuity:-(1/r) d_r(r sigma~ u~) - (1/r) d_phi(sigma~ v~)

## §5 SECOND-ORDER MODE COUPLING (eqs 28-34) — VERIFIED (2026-06-19)
Ansatz (28): f~ = f~0 + e^{c1 t} (1/2)[ f1 e^{i m(Omega_p t - phi)} + c.c.], m=2, c1=gamma1.
The m=2 self-interaction (quadratic, ~ e^{2 c1 t}) forces an axisymmetric (m=0) response:
  (29) d_t u0 - 2 Omega v0 + d_r[(cs^2/Sigma0) sigma0 + Phi0] = NL1(r) e^{2c1 t}
  (31) d_t v0 + (r Omega' + 2 Omega) u0                       = NL2(r) e^{2c1 t}
  (33) d_t sigma0 + (1/r) d_r(r Sigma0 u0)                    = NL3(r) e^{2c1 t}

Independently re-derived the forcing functions (substitute ansatz into 9-11, take m=0 / azimuthal avg):
  NL1 = v1 v1*/(2r) - (1/4) d_r(u1 u1*) + (i m/4r)(v1* u1 - v1 u1*)          -> MATCHES eq (30)
  NL2 = -(1/4r)(u1 v1* + v1 u1*) - (1/4)(u1 d_r v1* + u1* d_r v1)            -> MATCHES eq (32)
  NL3 = -(1/r) d_r[ (r/4)(sigma1 u1* + sigma1* u1) ]                        -> MATCHES eq (34)
All three match term-by-term. NL1's last term is the one the authors corrected from LKA97
("opposite sign ... error discovered after LKA went to press") — the LKA98 sign is the correct one
(my derivation confirms it). NL2's first term requires the -v~u~/r curvature term in (10) (confirmed
present in the printed equation). NL3 = -div of the m=0 mass flux (1/2)Re(sigma1 u1*) — the transport.

Status: §5 algebra CORRECT. Optional numerical check: measure the m=0 response (sigma0, mass flux)
in the from-scratch hydro and compare to integrating (29)-(34) with the linear eigenfunctions (the
paper's Fig 11 shows "nearly perfect agreement").

## §6 THIRD-ORDER / SATURATION (eqs 35-44) — in progress
Full ansatz (35): f~ = f~0(r,t) + e^{c1 t}[ (1/2)(f1 e^{i m chi}+c.c.) + (1/2)(f2 e^{2i m chi}+c.c.) ],
chi = Omega_p t - phi, m=m1=2; f1,f2 now time-dependent amplitudes (not fixed eigenfunctions).
Nine coupled equations: m=0 (36-38), m=2 (39-41, the cubic saturation feedback), m=4=2m1 (42-44).

CHECKED so far (all CORRECT by independent derivation):
- (36)-(38) m=0: identical to the second-order forcing (29-34); verified.
- (42) m=4 continuity: forcing = [-(1/2r) d_r(r sigma1 u1) + (i m1/r) sigma1 v1] e^{2c1 t}.
  My derivation (m=4 = e^{2i chi} part of -(1/r)d_r(r sigma~u~) - (1/r)d_phi(sigma~v~), with d_phi e^{2i chi}
  = -2 i m1) gives exactly +(i m1/r) sigma1 v1 -> MATCHES the image.

IMPORTANT METHOD NOTE: the pdftotext layer GARBLES SIGNS in these equations. Two apparent
discrepancies ((38) and (42) v1-term signs) were BOTH text-layer artifacts; the rendered images
confirm the equations are correct. => verify every sign against the page image, never the text layer.

## FULL SYMBOLIC VERIFICATION (2026-06-19) — ALL EQUATIONS CORRECT
Used sympy (lka98/derive_nonlinear.py) to mechanically substitute the three-mode ansatz (35)
into the nonlinear eqs (9-11) for gamma_p=2 and collect the m=0, m=2, m=4 harmonics
(d/dphi = -i m E d/dE; forcing for mode n = [n!=0 ? 2 : 1] x coeff of E^n). Compared each
result against the PAGE IMAGES (not the text layer). EVERY forcing term matches term-by-term:

  eq   mode      what                                            result
  ---  --------  ----------------------------------------------  ------
  30   m=0 rad   NL1                                             MATCH
  32   m=0 azim  NL2 (needs the -v~u~/r curvature term)          MATCH
  34   m=0 cont  NL3                                             MATCH
  39   m=2 cont  sigma0 u1 + sigma1 u0 + 1/2 u2 sigma1* + ...    MATCH (all 1/2's, conjugates)
  40   m=2 rad   ...- d_r(u0 u1 + 1/2 u2 u1*) + (im/r)(...-1/2 v2 u1*)  MATCH (incl. -1/2 sign)
  41   m=2 azim  geometric + u d_r v + im, all 1/2's/conjugates  MATCH
  42   m=4 cont  -(1/2r) d_r(r sigma1 u1) + (im/r) sigma1 v1     MATCH
  43   m=4 rad   (1/2r) v1^2 - 1/2 u1 d_r u1 + (im/2r) v1 u1     MATCH
  44   m=4 azim  -(1/2r) u1 v1 - 1/2 u1 d_r v1 + (im/2r) v1^2    MATCH

CONCLUSION: the LKA98 weakly-nonlinear governing equations (29-44) are ALGEBRAICALLY CORRECT in
their entirety. No errors found. The one historical slip (LKA97's NL1 sign) was already corrected
in LKA98, and the sympy independently confirms the corrected sign. (Linear LHS operators of 39-44
reduce to the verified linear operators 23-25 with m1 -> 2 m1 for the m=4 set.)
WARNING repeated: the pdftotext layer garbles signs/conjugates; (38) and (42) looked wrong in the
text but were correct in the images. Always check the page image.

## SATURATION via the nine-equation system (2026-06-19) — SOLVED
Two approaches in lka98/saturation.py (lab-frame 3-harmonic truncation of the 2D hydro,
dz_m/dt = L_m z_m + N_m(z), L_m the linearised operator per harmonic, N pseudo-spectral):

(A) Direct IVP, IMEX (implicit L, explicit pseudo-spectral N).  Linear operators validated
   (L2 leading eigenvalue = i*omega = 0.722+3.15j, exact). Linear growth + 2nd-order driving
   correct (m=2 ~ e^{gamma t}, m=0/m=4 ~ e^{2 gamma t}). BUT blows up at t~6 at the inner edge
   BEFORE saturating: the disk is MULTIPLY unstable (secondary m=2 mode ~0.4, m=4 self-mode ~0.38),
   and m=4 (driven at 2 gamma > gamma) overtakes m=2 with no cascade sink in the truncation.
   Extending to m=0..8 + viscosity delayed but did not cure it.  => This is exactly the
   instability LKA98 hit in 1998 ("susceptible to numerical instabilities"); the severe
   truncation + edge-concentrated modes make the raw IVP genuinely ill-behaved.

(B) Stuart-Landau / eigenmode (Galerkin) reduction -- THE CLEAN ANSWER.  Project on the dominant
   m=2 eigenmode phi2 (s2=gamma1+i 2 Om_p); slave the 2nd-order responses
     zeta0=(2 gamma1 - L0)^-1 N0(phi2,phi2*),  zeta4=(2 s2 - L4)^-1 N4(phi2,phi2),
   form the cubic back-reaction G=N2(zeta0,phi2)+N2(zeta4,phi2*), project on the adjoint:
     dA/dt = s2 A - beta |A|^2 A,   beta = -<psi2,G>/<psi2,phi2>.
   RESULT:  beta = 312 + 174 j.  Re(beta) > 0  => SUPERCRITICAL => the m=2 mode SATURATES
   (the back-reaction is stabilising -- the physics the raw IVP couldn't reach).
     saturated m=2 contrast |C2|/|C0|:  theory 0.27  vs  hydro 0.24   (~10%)
     pattern speed at saturation:       theory 1.57->1.37  vs  hydro 1.46->1.34  (both DROP ~10%)
   The Im(beta)>0 phase shift = the pattern-speed decrease the paper reports (Fig 1).
   Robust to viscosity/resolution.  Driver: lka98/run_saturation.py -> lka98/saturation.png.

CONCLUSION: the nine-equation system DOES predict saturation; the 1998 IVP instability was
numerical (severe truncation), not physical. The weakly-nonlinear amplitude equation recovers
the saturation amplitude AND the pattern-speed shift, both matching the full hydro.

(Algebra of the 9 equations: fully verified above, all correct.)
