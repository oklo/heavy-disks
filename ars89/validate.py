"""Validation against ARS89 Fig. 3 (canonical m=1 mode).

Run as a module:

    python -m ars89.validate              # solve + print + save figure
    python -m ars89.validate --N 2000     # finer grid

Canonical model (ARS89 Section IV / Fig. 3 caption):
    p = 3/2, q = 1/2, R_D/R_* = 1e4, M_D = M_*, Q_* = 10.
ARS89 report the lowest-order eigenfrequency omega = 4.26 - 0.232 i (units Omega_D);
the pattern speed Omega_p = Re(omega) sets corotation at R_CR = 0.452 R_D.
"""

from __future__ import annotations

import argparse
import numpy as np

from .model import DiskModel
from .eigensolve import EigenProblem, solve_mode, ARS89_FIG3_EIGENVALUE
from .energy import energy_budget


def run(N: int = 1200, eta: float = 0.05, make_plot: bool = True, plot_path: str = "fig3_eigenfunction.png"):
    model = DiskModel(p=1.5, q=0.5, Rd_over_Rstar=1.0e4, Md_over_Mstar=1.0, Qstar=10.0, m=1)
    problem = EigenProblem.build(model, N=N, eta=eta)

    omega, S, info = solve_mode(problem, omega_guess=ARS89_FIG3_EIGENVALUE)

    target = ARS89_FIG3_EIGENVALUE
    print("ARS89 reproduction — canonical m=1 mode (p=3/2, q=1/2, R_D/R_*=1e4, M_D=M_*, Q_*=10)")
    print(f"  grid: N = {N}, log-radial, eta = {eta}")
    print(f"  eigenfrequency  omega = {omega.real:+.4f} {omega.imag:+.4f} i   (units Omega_D)")
    print(f"  ARS89 Fig. 3    omega = {target.real:+.4f} {target.imag:+.4f} i")
    print(f"  pattern speed   Omega_p = {omega.real:.4f}   (ARS89: {target.real:.4f})")
    print(f"  growth rate     gamma   = {-omega.imag:.4f}   (ARS89: {-target.imag:.4f})")
    re_err = abs(omega.real - target.real) / abs(target.real)
    im_err = abs(omega.imag - target.imag) / abs(target.imag)
    print(f"  relative error: Re {100*re_err:.1f}%, Im {100*im_err:.1f}%")
    print(f"  solver: nit={info['nit']}, residual sigma_min={info['residual']:.2e}, "
          f"converged={info['converged']}")

    # corotation radius from the computed pattern speed
    r = problem.grid.r
    icr = int(np.argmin(np.abs(problem.eq.Omega - omega.real)))
    print(f"  corotation radius R_CR = {r[icr]:.3f} R_D   (ARS89: 0.452 R_D)")

    # --- energy analysis (ARS89 eq 33): independent growth-rate check ---
    eb = energy_budget(problem, omega, S)
    print("  energy budget (eq 33):")
    print(f"    gamma (eigenvalue) = {eb.gamma_eigenvalue:.5f}")
    print(f"    gamma (energy)     = {eb.gamma_energy:.5f}   "
          f"[LHS=RHS to {100*eb.relative_error:.2f}%]")

    if make_plot:
        _plot_eigenfunction(problem, omega, S, plot_path)
        print(f"  eigenfunction figure written to {plot_path}")
        _plot_energy(problem, eb, "fig5_energy.png")
        print("  energy-budget figure written to fig5_energy.png")

    return model, problem, omega, S, info, eb


def _plot_energy(problem, eb, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    r = problem.grid.r
    Rstar = problem.model.Rstar
    x = np.log10(r / Rstar)
    # ARS89 Fig. 5 normalisation: work per unit area divided by sigma_0 Omega^2 r^2
    norm = problem.eq.sigma0 * problem.eq.Omega**2 * r**2

    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.plot(x, eb.reynolds / norm, "-", lw=1.0, label="Reynolds")
    ax.plot(x, eb.acoustic / norm, ":", lw=1.2, label="acoustic flux")
    ax.plot(x, eb.grav_direct / norm, "--", lw=1.0, label="gravity (direct)")
    ax.plot(x, eb.grav_indirect / norm, "-.", lw=1.2, label="gravity (indirect/SLING)")
    ax.set_xlabel(r"$\log_{10}(r/R_*)$")
    ax.set_ylabel(r"work$\,/\,\sigma_0\Omega^2 r^2$")
    ax.set_title(
        "ARS89 Fig. 5 reproduction — modal energy budget  "
        rf"($\gamma_{{E}}={eb.gamma_energy:.3f}$ vs $\gamma={eb.gamma_eigenvalue:.3f}$)"
    )
    ax.legend(loc="upper left", fontsize=8)
    ax.axhline(0, color="0.7", lw=0.6)
    # the work terms are concentrated in the outer resonant cavity; scale to them
    cavity = np.log10(r / Rstar) > 2.4
    pk = np.max(np.abs(eb.grav_direct[cavity] / norm[cavity]))
    ax.set_ylim(-1.2 * pk, 1.2 * pk)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _plot_eigenfunction(problem, omega, S, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    r = problem.grid.r
    Rstar = problem.model.Rstar
    x = np.log10(r / Rstar)

    # rotate global phase so the real part carries the bulk of the amplitude
    phase = np.angle(S[np.argmax(np.abs(S))])
    Srot = S * np.exp(-1j * phase)

    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.plot(x, np.real(Srot), "-", lw=1.1, color="k", label="Re S(r)")
    ax.plot(x, np.imag(Srot), "--", lw=1.1, color="C3", label="Im S(r)")
    ax.set_xlabel(r"$\log_{10}(r/R_*)$")
    ax.set_ylabel(r"$S(r)=\sigma_1/\sigma_0$  (arb. units)")
    ax.set_title(
        f"ARS89 Fig. 3 reproduction — m=1 mode,  "
        rf"$\omega={omega.real:.3f}{omega.imag:+.3f}i\,\Omega_D$"
    )
    ax.legend(loc="upper left")
    ax.axhline(0, color="0.7", lw=0.6)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Reproduce ARS89 Fig. 3 eigenvalue/eigenfunction.")
    ap.add_argument("--N", type=int, default=1200, help="number of radial grid points")
    ap.add_argument("--eta", type=float, default=0.05, help="rotation-curve softening (ARS89 A5)")
    ap.add_argument("--no-plot", action="store_true", help="skip the eigenfunction figure")
    args = ap.parse_args()
    run(N=args.N, eta=args.eta, make_plot=not args.no_plot)


if __name__ == "__main__":
    main()
