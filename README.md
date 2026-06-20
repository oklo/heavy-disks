# Dynamics of heavy gaseous disks

Independent, from-scratch reproductions of two razor-thin self-gravitating disk
instability papers, built outward from linear theory to nonlinear hydrodynamics:

> **Laughlin, Korchagin & Adams 1998, ApJ, 504, 945** — *The Dynamics of Heavy
> Gaseous Disks* ([1998ApJ...504..945L](https://ui.adsabs.harvard.edu/abs/1998ApJ...504..945L))
>
> **Adams, Ruden & Shu 1989, ApJ, 347, 959** — *Eccentric Gravitational
> Instabilities in Nearly Keplerian Disks* ([1989ApJ...347..959A](https://ui.adsabs.harvard.edu/abs/1989ApJ...347..959A))

Each result is cross-checked by independent methods until it agrees with itself.
The journal PDFs/scans are not redistributed here (AAS/IOP copyright); the
equations used are transcribed and verified in the `*.md` notes.

---

## Components

### `lka98/` — heavy-disk m=2 spiral, linear → weakly-nonlinear → saturation
The main effort. The unstable two-armed (m=2) spiral of the "standard reference"
heavy disk, attacked four independent ways that all agree:

- **Linear global modes** (`eigensolve.py`, `directint.py`): matrix method *and*
  direct-integration shooting, plus an **energy-integral cross-check** (`energy.py`).
  All give γ ≈ 0.72–0.76, Ω_p ≈ 1.5–1.6.
- **Softened self-gravity** (`poisson.py`, `selfgravity.py`) via complete elliptic
  integrals; rotation curve incl. disk self-gravity.
- **Weakly-nonlinear theory fully re-derived** (`derive_nonlinear.py`): a `sympy`
  derivation that substitutes the three-mode ansatz into the nonlinear equations and
  mechanically reproduces **every** forcing term of the paper's 2nd- and 3rd-order
  systems (eqs 29–44). All verified term-by-term — see `PERTURBATION.md`.
- **Nonlinear saturation** (`saturation.py`): the full nine-equation three-wave system
  that the 1998 paper could not integrate. The direct IVP (IMEX + pseudo-spectral)
  reproduces the linear growth and second-order driving but blows up — for reasons we
  diagnose. The **Stuart–Landau reduction** then gives the clean answer:
  `dA/dt = s₂A − β|A|²A` with `Re β > 0` (supercritical → saturates).

  | | saturated \|C2\|/\|C0\| | Ω_p (linear → saturated) |
  |---|---|---|
  | Stuart–Landau theory | 0.27 | 1.57 → 1.37 |
  | Full hydro simulation | 0.24 | 1.46 → 1.34 |

  The amplitude *and* the pattern-speed drop both match the hydro.

**Findings worth flagging** (`EQUATIONS.md`, `PERTURBATION.md`):
- A confirmed **typesetting error** in the published softening modulus (eq C8).
- The quoted matrix/direct eigenvalue (γ=0.804, Ω_p=1.70) is a mild outlier; five
  independent determinations here — and the paper's *own* simulation — cluster at
  γ≈0.72–0.76, Ω_p≈1.5–1.6.
- The paper's full weakly-nonlinear algebra (eqs 29–44) is **error-free**; the one
  historical slip (a sign in the 1997 precursor) was already corrected in 1998.

### `diskfft/` — from-scratch 2D razor-thin hydro (C++)
A self-contained reimplementation of the V2D.F approach: 2nd-order van-Leer (MUSCL)
advection on a polar grid + Kalnajs (1971) logarithmic-spiral **FFT self-gravity**,
barotropic EOS. Reproduces the m=2 spiral growth and its nonlinear saturation — the
fifth independent confirmation of γ≈0.72–0.76. Optimized with **FFTW + OpenMP**
(~9.5× over the initial radix-2 single-core version; a full 256² self-gravitating run
through saturation takes tens of seconds). See `diskfft/README.md`.

### `ars89/` — m=1 SLING eigenmode code (Python)
The eccentric (m=1) instability of a nearly-Keplerian disk, where the reflex motion of
the central star destabilizes a one-armed mode. Recovers the paper's canonical
eigenvalue to ~1% (pattern speed) / ~4% (growth rate). See `ars89/README.md`.

---

## Running

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# LKA98 weakly-nonlinear algebra check + saturation
python lka98/derive_nonlinear.py            # symbolic verification of eqs 29-44
python lka98/run_saturation.py              # Stuart-Landau saturation -> lka98/saturation.png

# ARS89 linear eigenvalue
python -m ars89.validate

# diskfft hydro (needs FFTW + OpenMP; see diskfft/Makefile)
cd diskfft && make && ./diskdisk 256 256 0.005 16 1 1e-3 8
```

Dependencies: `numpy`, `scipy`, `sympy`, `matplotlib`, `pytest`.

## License

MIT — see [LICENSE](LICENSE). Note this covers the code only; the reproduced ApJ
articles remain under AAS/IOP copyright and are not redistributed here.
