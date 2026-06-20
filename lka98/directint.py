"""Direct-integration (shooting) method for the linear m=2 mode (LKA98 §4).

Independent of the matrix method: integrate the reduced first-order ODEs (26,27)
for (sigma1, u1) from the inner boundary with u1(Rin)=0, sigma1(Rin)=1, using a
trial eigenfrequency omega and a trial potential Psi1(r); recompute Psi1 = P sigma1
(softened m-harmonic) and iterate to self-consistency; then Newton-iterate omega so
that the reflective outer condition u1(RD)=0 is satisfied.

This shares the governing equations and the Poisson kernel with eigensolve.py but
uses an entirely different numerical scheme (adaptive RK + fixed-point potential
iteration + complex Newton), so agreement is a meaningful cross-check.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from .model import DiskModel
from .selfgravity import rotation_curve
from .poisson import build_poisson_matrix


class ShootingProblem:
    def __init__(self, model: DiskModel, N: int = 800, modulus: str = "correct"):
        self.model = model
        self.r = np.linspace(model.Rin, model.RD, N)
        Om, kap, _, _ = rotation_curve(model, self.r, modulus=modulus)
        self.Om = CubicSpline(self.r, Om)
        self.kap = CubicSpline(self.r, kap)
        self.P = build_poisson_matrix(model, self.r, model.m)

    # --- equilibrium coefficient pieces at arbitrary r --------------------

    def _eq(self, r):
        m = self.model
        Sig = m.S0 * np.exp(-((r - m.R0) ** 2) / m.w2)
        dSig = Sig * (-2.0 * (r - m.R0) / m.w2)
        cs2 = m.K * m.gamma_p * Sig ** (m.gamma_p - 1.0)
        # d/dr(cs2/Sig) = K gamma_p (gamma_p-2) Sig^{gamma_p-3} dSig
        dcs2_over_Sig = m.K * m.gamma_p * (m.gamma_p - 2.0) * Sig ** (m.gamma_p - 3.0) * dSig
        drSig = Sig + r * dSig                      # d/dr(r Sig)
        return Sig, cs2, dcs2_over_Sig, drSig

    def rhs(self, r, y, omega, Psi1, dPsi1):
        """ODE RHS for y=[sigma1, u1] at radius r (eqs 26,27).

        Psi1, dPsi1 are callables (splines of the current potential iterate).
        """
        m = self.model.m
        sigma1, u1 = y
        Om = float(self.Om(r)); kap = float(self.kap(r))
        Sig, cs2, dcs2_over_Sig, drSig = self._eq(r)
        nu = (omega - m * Om) / kap
        A = 1j * kap * nu
        P1 = Psi1(r); dP1 = dPsi1(r)
        W = (cs2 / Sig) * sigma1 + P1

        a11 = 2.0 * Om * m / (nu * kap * r) - (Sig / cs2) * dcs2_over_Sig
        a12 = -A * Sig * (nu**2 - 1.0) / (cs2 * nu**2)
        b1 = 2.0 * Om * Sig * m / (cs2 * nu * kap * r)
        c1 = -Sig / cs2
        a21 = 1j * m**2 * cs2 / (nu * r**2 * kap * Sig) - A / Sig
        a22 = 1j * m * A / (2.0 * nu**2 * Om * r) - drSig / (Sig * r)
        b2 = 1j * m**2 / (nu * r**2 * kap)

        dsigma = a11 * sigma1 + a12 * u1 + b1 * P1 + c1 * dP1
        du = a21 * sigma1 + a22 * u1 + b2 * P1
        return [dsigma, du]

    def shoot(self, omega, Psi1_vals):
        """One integration pass: returns sigma1(grid), u1(grid), u1(RD).

        Psi1_vals is the current potential on the grid; interpolated (with its
        derivative) for the RK substeps.
        """
        Psi1 = CubicSpline(self.r, Psi1_vals)
        dPsi1 = Psi1.derivative()
        sol = solve_ivp(
            self.rhs, (self.r[0], self.r[-1]), [1.0 + 0j, 0.0 + 0j],
            t_eval=self.r, args=(omega, Psi1, dPsi1),
            method="RK45", rtol=1e-8, atol=1e-10, dense_output=False,
        )
        sigma1 = sol.y[0]
        u1 = sol.y[1]
        return sigma1, u1, u1[-1]

    def newton_omega(self, omega, Psi1_vals, tol=1e-9, max_iter=40, h=1e-5):
        """Inner loop: Newton on omega so u1(RD)=0 at FIXED potential Psi1_vals."""
        for _ in range(max_iter):
            _, _, f0 = self.shoot(omega, Psi1_vals)
            if abs(f0) < tol:
                break
            _, _, fp = self.shoot(omega + h, Psi1_vals)
            _, _, fm = self.shoot(omega - h, Psi1_vals)
            step = -f0 / ((fp - fm) / (2.0 * h))
            if abs(step) > 0.3:
                step *= 0.3 / abs(step)
            omega = omega + step
            if abs(step) < tol:
                break
        return omega


def solve_shooting(problem: ShootingProblem, omega_guess: complex, Psi1_init,
                   n_outer: int = 30, relax: float = 1.0, tol: float = 1e-7):
    """Outer loop: Newton on omega at fixed Psi1, then revise Psi1 = P sigma1.

    Psi1_init is the seed potential on the grid (e.g. from the matrix eigenfunction,
    normalised so sigma1(Rin)=1).  Following LKA98, the direct-integration method
    refines a mode already located by the matrix method; it cross-checks the
    discretisation/BC handling since it shares neither with the matrix solver.
    """
    Psi1 = np.asarray(Psi1_init, dtype=complex).copy()
    omega = complex(omega_guess)
    history = []
    for outer in range(n_outer):
        omega = problem.newton_omega(omega, Psi1)
        sigma1, u1, _ = problem.shoot(omega, Psi1)
        Psi1_new = problem.P @ sigma1
        dP = np.max(np.abs(Psi1_new - Psi1)) / (np.max(np.abs(Psi1_new)) + 1e-30)
        history.append((omega, dP))
        Psi1 = (1.0 - relax) * Psi1 + relax * Psi1_new
        if dP < tol:
            break
    sigma1, u1, res = problem.shoot(omega, Psi1)
    m = problem.model.m
    info = {"n_outer": outer + 1, "rel_dPsi": dP, "residual_uRD": abs(res),
            "Omega_p": omega.real / m, "gamma1": -omega.imag, "history": history}
    return omega, sigma1, u1, info
