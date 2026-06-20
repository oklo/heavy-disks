"""Linear energy-integral cross-check (LKA98 Appendix A, eq A1).

For a self-gravitating mode with m>1 the linearized equations (23-25) imply
(Broadbent & Moore 1979):

    gamma sigma0 <u1^2 + v1^2 + h1^2/cs^2>
        = -r (dOmega/dr) sigma0 <u1 v1>          (Reynolds stress)
          - (1/r) d/dr( r sigma0 <h1 u1> )        (acoustic flux)
          - sigma0 <u1 . grad Psi1>,              (gravitational work)

with <.> the azimuthal average.  For modal amplitudes (hatted), the azimuthal
average of a product is <X Y> = (1/2) Re[X_hat Y_hat*]; the common 1/2 cancels in

    gamma_energy = Integral(RHS) / Integral(sigma0 (|u1|^2+|v1|^2+|h1|^2/cs^2)),

which must equal the eigenvalue growth rate gamma1.  As in the ARS89 check, the
azimuthal gravitational-work term is Re[v1_hat (i m/r) Psi1_hat*] -- the (i m/r)
factor is NOT conjugated, only Psi1_hat.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .eigensolve import LinearProblem


@dataclass
class EnergyBudget:
    r: np.ndarray
    energy_density: np.ndarray
    reynolds: np.ndarray
    acoustic: np.ndarray
    gravity: np.ndarray
    u1: np.ndarray
    v1: np.ndarray
    gamma_energy: float
    gamma_eigenvalue: float

    @property
    def relative_error(self) -> float:
        return abs(self.gamma_energy - self.gamma_eigenvalue) / abs(self.gamma_eigenvalue)


def energy_budget(problem: LinearProblem, omega: complex, sigma1: np.ndarray,
                  u1: np.ndarray) -> EnergyBudget:
    """Evaluate the LKA98 linear energy budget (A1) for a computed mode."""
    m = problem.model.m
    r = problem.r
    Om, kap = problem.Omega, problem.kappa
    cs2, Sig = problem.cs2, problem.Sigma0
    D, P = problem.D, problem.P

    nu = (omega - m * Om) / kap
    A = 1j * kap * nu
    h1 = (cs2 / Sig) * sigma1
    Psi1 = P @ sigma1
    W = h1 + Psi1
    # azimuthal momentum (25): A v1 + (kappa^2/2Omega) u1 - (i m/r) W = 0
    v1 = ((1j * m / r) * W - (kap**2 / (2.0 * Om)) * u1) / A

    dOmega = D @ Om
    dPsi1 = D @ Psi1

    # energy density sigma0 (|u1|^2+|v1|^2+|h1|^2/cs^2) (1/2 dropped; cancels in ratio)
    energy = Sig * (np.abs(u1) ** 2 + np.abs(v1) ** 2 + np.abs(h1) ** 2 / cs2)

    reynolds = -r * dOmega * Sig * np.real(u1 * np.conj(v1))
    flux = r * Sig * np.real(h1 * np.conj(u1))
    acoustic = -(1.0 / r) * (D @ flux)
    gravity = -Sig * np.real(u1 * np.conj(dPsi1) + v1 * (1j * m / r) * np.conj(Psi1))

    dA = 2.0 * np.pi * r
    num = np.trapezoid((reynolds + acoustic + gravity) * dA, r)
    den = np.trapezoid(energy * dA, r)
    gamma_energy = num / den

    return EnergyBudget(
        r=r, energy_density=energy, reynolds=reynolds, acoustic=acoustic,
        gravity=gravity, u1=u1, v1=v1,
        gamma_energy=float(np.real(gamma_energy)),
        gamma_eigenvalue=float(-omega.imag),
    )
