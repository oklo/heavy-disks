# Athena++ LKA98 disk — setup notes

Athena++ install: `/Users/greglaughlin/src/athena` (C++ Athena++, Stone et al. 2020).

## Problem generator
`src/pgen/lka98_disk.cpp` — 2D razor-thin cylindrical (R,phi) disk:
- IC: Gaussian Sigma0 = S0 exp[-(R-R0)^2/w2] (S0=0.372, R0=0.45, w2=0.05), Rin=0.25, RD=1.0.
- Polytropic via adiabatic EOS gamma=2 with uniform entropy (P = K Sigma^2, K=0.25) ≡ barotropic.
- Central star gravity m*=0.6 applied as an EnrollUserExplicitSourceFunction (-GM/R^2 radial).
- Equilibrium rotation balances ONLY star + pressure (no disk self-gravity yet):
  vphi^2 = GM/R - 4 K R Sigma0 (R-R0)/w2.
- m=2 density seed (amp=1e-4 cos 2phi).
- Reflecting radial BCs (= LKA98 reflective u1=0 edges), periodic in phi.

## Build & run
```
cd /Users/greglaughlin/src/athena
python3 configure.py --prob lka98_disk --coord cylindrical --eos adiabatic --flux hllc
make -j4
cd runs/lka98 && ../../bin/athena -i ../../inputs/hydro/athinput.lka98
```
Input: `inputs/hydro/athinput.lka98` (nx1=200 R, nx2=256 phi, nx3=1; VTK + hst output).

## Status (2026-06-19)
Runs stably to t=5 (4971 cycles, 21 s, 1.2e7 zone-cyc/s). Mass conserved exactly (0.39996),
radial KE ~1e-5 (equilibrium holds). NON-self-gravitating reference confirmed.

## TODO
- Build the razor-thin self-gravity (Kalnajs 1971 log-spiral FFT) — NOT native in Athena++
  (its FFT gravity is 3D Poisson). Add as a custom source term; reuse the same solver in the
  from-scratch van-Leer hydro code. With self-gravity, switch IC rotation to the full
  (self-gravity-included) equilibrium curve and expect the m=2 spiral to grow at gamma~0.74.
- Fourier diagnostics (eqs 15-20): |C_m(t)|, pattern speed, growth rate from the VTK outputs.
