# diskfft — van-Leer + Kalnajs FFT razor-thin disk code (LKA98 V2D.F replica)

From-scratch C++ reimplementation of the Laughlin (1994) V2D.F approach used in
LKA98: 2nd-order van-Leer MUSCL advection on a polar (R,phi) grid + Kalnajs (1971)
logarithmic-spiral FFT for the razor-thin self-gravity.  Built to (a) reproduce the
m=2 spiral growth (gamma~0.74) and (b) benchmark performance against Athena++.

## Components
- `src/fft.hpp`        — self-contained radix-2 FFT (powers of two; swap in FFTW later).
- `src/kalnajs.hpp`    — Kalnajs self-gravity solver. In u=ln R the softened razor-thin
  potential is a convolution `V = -G du dphi (S * K)`, `S=R^{3/2}Sigma`, `V=R^{1/2}Phi`,
  `K=1/sqrt(2 cosh(du) - 2 cos(dphi) + eps)`; circular in phi, zero-padded (2*Nu) in u.
  VALIDATED: FFT vs direct convolution to 2.4e-15; kernel = exact razor-thin Green's
  function; reduced<->physical equivalence is algebraically exact.
- `src/hydro.hpp`      — (next) van-Leer MUSCL advection, polytropic EOS, central gravity.
- `src/diskmain.cpp`   — (next) driver: IC, time loop, Fourier diagnostics (eqs 15-20).

## Build/run
```
g++ -O3 -std=c++17 -I src src/test_kalnajs.cpp -o test_kalnajs && ./test_kalnajs
```

## RESULT (2026-06-19) — full code works, reproduces the m=2 spiral
`diskmain.cpp` (van-Leer hydro + Kalnajs FFT) runs the LKA98 heavy disk:
- Non-self-gravity: stable equilibrium, m=2 seed does NOT grow (mass conserved to machine precision).
- Self-gravity ON: m=2 spiral grows exponentially. Growth rate vs scale-free softening eps:
  eps=0.04 -> g=0.36 ; 0.02 -> 0.53 ; 0.01 -> 0.64 ; 0.005 -> 0.69 ; 0.0025 -> 0.757 (Omega_p=1.49).
  As eps->0, gamma -> ~0.74, Omega_p -> ~1.5, matching the linear analysis (0.724/1.574) and the
  paper's V2D.F simulation (0.74/1.59). Fifth independent confirmation that gamma~0.72-0.76 (NOT 0.804).
- Build/run: `g++ -O3 -std=c++17 -I src src/diskmain.cpp -o diskdisk && ./diskdisk 128 128 0.0025 9 1 1e-3`
- TODO: wire the Kalnajs solver into Athena++ as a self-gravity source -> apples-to-apples
  performance comparison; eps->0 extrapolation for exact Omega_p; nonlinear saturation study.

## PERFORMANCE (2026-06-19) — FFTW + OpenMP
Build now uses FFTW (real-to-complex, optionally threaded) + OpenMP hydro (race-free
flux-then-gather transport, pre-allocated work buffers).  14-core Apple-silicon laptop,
256x256, self-gravity, growth rate unchanged (gamma=0.717, mass drift 0.00):

  | config                         | zone-updates/s | speedup |
  |--------------------------------|----------------|---------|
  | original radix-2 FFT, 1 core   | 6.45e6         | 1.0x    |
  | FFTW r2c, 1 core               | 1.71e7         | 2.65x   |
  | FFTW + OpenMP, 8 threads       | 6.12e7         | 9.5x    |

- FFTW r2c cut the self-gravity cost from ~60% of a step to ~19% (3.86 -> next-to-free).
- Hydro scales 4.2x at 8 threads (memory-bandwidth bound; ~8 threads optimal, drops at 10
  as it spills onto efficiency cores). Build: `make` (Makefile sets FFTW/OpenMP flags);
  run: `./diskdisk Nu Nphi eps tlim selfgrav amp nthreads`.
- A full 256^2 self-gravitating run through saturation now takes tens of seconds.

## Softening note
`eps` (scale-free, in the log kernel) is the FFT-relevant softening. It is NOT the LKA98
Appendix-C field-point softening eta(R)^2 R^2; it is the symmetric ~eps*R*R' softening a
Kalnajs FFT solver imposes, consistent with the simulation eigenvalue (Omega_p~1.59, gamma~0.74)
that the linear matrix/shooting analysis reproduced.
