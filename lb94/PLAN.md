# LB94 — Stage 1: 2D axisymmetric collapse to a protostellar disk

Goal: reproduce, from scratch in C++, the 2D radiation-hydro collapse that produces the
protostellar disk handed to the 3D SPH code in **Laughlin & Bodenheimer 1994 (ApJ 436, 335)**.
The collapse calculation is the YBL93 standard case (**Yorke, Bodenheimer & Laughlin 1993,
ApJ 411, 274**) plus the LB94 addition of physical (alpha) viscosity. Numerical method = the
Explicit Nested Grids (ENG) scheme of **Yorke & Kaisig 1994** on the **Różyczka 1985** 2D
Eulerian hydro base.

## Physical model (Black & Bodenheimer 1975 eqs 1-4; cylindrical R,Z; axisymmetric)
2D grid in (R, Z), single quadrant, mirror symmetry about Z=0, axis at R=0. Three velocity
components (v_R, v_Z, v_phi) — axisymmetric but rotating. State per cell:
rho, rho v_R, rho v_Z, angular momentum density A = rho R v_phi, (internal energy e for RT).
- continuity:   d_t rho + (1/R) d_R(R rho v_R) + d_Z(rho v_Z) = 0
- R-momentum:   d_t(rho v_R) + div(rho v_R v_pol) = -d_R P - rho d_R Phi + A^2/(rho R^3)   [centrifugal]
- Z-momentum:   d_t(rho v_Z) + div(rho v_Z v_pol) = -d_Z P - rho d_Z Phi
- ang. mom.:    d_t A + div(A v_pol) = 0      (specific j = R v_phi advected; no phi-forces)
- Poisson:      (1/R) d_R(R d_R Phi) + d_ZZ Phi = 4 pi G rho
- EOS:          Stage 1a barotropic (Larson: isothermal c_s^2 rho below rho_crit ~1e-13,
                then adiabatic gamma=7/5); Stage 1c full ideal-gas + dissociation + RT energy eq.
- shocks:       von Neumann-Richtmyer artificial viscosity (Różyczka 1985).
- central sink: unresolved zone interior to ~5 AU. Mass flows IN not OUT; track M_c, J_c, Mdot;
                accretion luminosity from a Maclaurin-spheroid model (BYRT90); held at 25 Lsun
                after 6000 yr in the standard case.

## ENG nested grids (Yorke & Kaisig; static Berger-Colella-style AMR, factor 2)
- n self-similar grids, each 60x60, fixed in time, centered on the core; grid l+1 is half the
  linear size of grid l (refinement factor 2). Standard case n=4 -> innermost res 5.55 AU,
  outer R_max=Z_max=4e16 cm (~2700 AU).
- time subcycling: inner grid takes 2 steps of dt_l, then outer takes 1 step of dt_{l-1}=2 dt_l,
  recursively. dt_n = min_k (1/2^k) dt_k^CFL. INNER grid integrated FIRST (reverse of
  Berger-Colella) — needed for radiation-transfer stability.
- BCs: Z=0 and R=0 as unnested (reflect/axis); each interior grid's outer boundary filled by
  bilinear interpolation from its parent. Restriction: outer cell = conservative sum of its 4
  child cells. Berger-Colella flux correction at coarse/fine interfaces -> exact mass conservation.
- ENG applied to BOTH explicit (hydro) and implicit (Poisson ADI, FLD radiation) substeps.

## Standard case (target to reproduce)
1 Msun, sphere R=4e16 cm, rho ~ r^-1, isothermal 20 K, Omega=4.4e-13 s^-1, J=2.8e53,
j_outer=7e20 cm^2/s, thermal/|grav|=0.37, rot/|grav|=0.01. Followed ~8e4 yr ->
0.45 Msun central object + ~equilibrium disk. (Keplerian radius of outer element ~250 AU.)

## Staged roadmap
- [1a] single-grid 2D cylindrical isothermal self-gravitating hydro + central sink.
       Tests: free-fall time of a uniform sphere; rotating-collapse disk formation.
- [1b] ENG nested grids. Test: Sedov blast wave across grid boundaries (d log R_S/d log t -> 0.4;
       YBL got 0.39). Mass conservation < 0.4%.
- [1c] radiative transfer (flux-limited diffusion, implicit) + opacities (Pollack et al 1985 dust,
       Alexander 1975 molecular). Or barotropic stand-in first.
- [1d] physical alpha-viscosity eta = alpha c_s H rho (LB94 addition); reproduce the LB94 disk.
- Stage 2: hand the 2D (R,Z) disk to a from-scratch 3D SPH code; m=1/m=2 spiral instabilities.
- Stage 3: resolution audit (Truelove 1997 / Bate & Burkert 1997 criteria postdate LB94).

## Validation targets
- free-fall collapse: central density runaway on ~t_ff = sqrt(3 pi/(32 G rho)).
- Sedov blast: R_S ~ t^0.4 across nested-grid boundaries.
- standard-case disk: ~equilibrium disk, central object ~0.45 Msun at 8e4 yr.
- mass conservation to <~0.4% over the run.

## References (in lb94/paper/, gitignored)
YBL93 (ApJ 411,274), LB94 (ApJ 436,335). Incoming: Różyczka 1985 (A&A 143,59),
Bodenheimer/Yorke/Różyczka/Tohline 1990 (ApJ 355,651), Pollack/McKay/Christofferson 1985,
Alexander 1975 (ApJS 29,363). Yorke & Kaisig 1994 (Comput.Phys.Comm.) — abstract only;
method reconstructed from the YBL93 description above + general AMR practice.
