r"""Energy analysis of a computed mode (ARS89 eq 33, Section IV.b).

For a normal mode the azimuthally-averaged disturbance energy per unit area

    <E> = (1/2) sigma_0 < u_1^2 + v_1^2 + h_1^2/a_0^2 >          (kinetic + acoustic)

grows as d<E>/dt = 2 gamma <E>, and ARS89 (eq 33) write its local budget as the
sum of three work terms,

    d<E>/dt = -r (dOmega/dr) sigma_0 <u_1 v_1>                  (Reynolds stress)
              - (1/r) d/dr ( r sigma_0 <h_1 u_1> )             (acoustic flux)
              - sigma_0 < u_1 . grad(psi_1 + psitilde_1) >,    (gravitational work)

where the gravitational term splits into the direct self-gravity (psi_1) and the
indirect / SLING (psitilde_1) contributions.  Equating the volume integrals of
the two sides gives an *independent* estimate of the growth rate,

    gamma_energy = Integral(RHS work) / ( 2 Integral(<E>) ),

which should match the eigenvalue gamma = -Im(omega).  This is ARS89's accuracy
check (their Fig. 5 shows the four work terms; the indirect term is comparable in
magnitude to the sum of the other three -- the hallmark of the SLING mechanism).

Azimuthal averages of mode quadratic quantities use
    <Re[A e^{i phi}] Re[B e^{i phi}]> = (1/2) Re[A B*],
so every term below carries the modal amplitudes (hatted) and a factor 1/2; the
common growth factor e^{2 gamma t} cancels between the two sides, and the overall
normalisation of the eigenfunction cancels in the ratio.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .eigensolve import EigenProblem


@dataclass
class EnergyBudget:
    """Radial profiles of the energy density and the four work terms."""

    r: np.ndarray
    energy_density: np.ndarray      # <E>(r)  (spatial part)
    reynolds: np.ndarray            # T1(r)
    acoustic: np.ndarray            # T2(r)
    grav_direct: np.ndarray         # T3 from psi_1
    grav_indirect: np.ndarray       # T3 from psitilde_1
    u1: np.ndarray                  # radial velocity amplitude
    v1: np.ndarray                  # azimuthal velocity amplitude
    gamma_energy: float             # growth rate from the energy budget
    gamma_eigenvalue: float         # -Im(omega)

    @property
    def work_total(self) -> np.ndarray:
        return self.reynolds + self.acoustic + self.grav_direct + self.grav_indirect

    @property
    def relative_error(self) -> float:
        return abs(self.gamma_energy - self.gamma_eigenvalue) / abs(self.gamma_eigenvalue)


def perturbation_potentials(problem: EigenProblem, omega: complex, S: np.ndarray):
    """Return (h1, psi1_direct, psi1_indirect) modal amplitudes for eigenfunction S.

    h1 = a_0^2 S,  psi1 = 2 pi G sigma_0 r (J S)  (eq 25b),
    psitilde1(r) = omega^2 R0 r with R0 = pi/(M_*+M_D) Integral(rho^2 sigma_1 drho) (eq 18b).
    """
    model = problem.model
    r = problem.grid.r
    eq = problem.eq

    h1 = eq.a0sq * S
    psi1 = 2.0 * np.pi * model.G * eq.sigma0 * r * (problem.J @ S)

    # indirect potential: R0 from the m=1 mass moment of the perturbation
    f = np.exp(problem.grid.log_f)
    w = r * 0.5 * (f - 1.0 / f)
    w = w.copy()
    w[0] *= 0.5
    w[-1] *= 0.5
    sigma1 = eq.sigma0 * S
    R0 = np.pi / model.Mtot * np.sum(r**2 * sigma1 * w)
    psi1_ind = omega**2 * R0 * r
    return h1, psi1, psi1_ind


def velocities(problem: EigenProblem, omega: complex, S: np.ndarray):
    """Radial/azimuthal velocity amplitudes from the momentum equations (7b,7c).

    With forcing W = psi_1 + psitilde_1 + h_1 and sigma_tilde = omega - m Omega,
    D = kappa^2 - sigma_tilde^2 (eq 14b):

        u_1 = (1/D)[ -i sigma_tilde W' + (2 i m Omega / r) W ]
        v_1 = (1/D)[ (kappa^2/2Omega) W' - (m sigma_tilde / r) W ].
    """
    model = problem.model
    m = model.m
    r = problem.grid.r
    eq = problem.eq
    ops = problem.ops

    h1, psi1, psi1_ind = perturbation_potentials(problem, omega, S)
    W = psi1 + psi1_ind + h1
    dW = ops.ddr(W)

    sig = omega - m * eq.Omega
    D = eq.kappa**2 - sig**2
    u1 = (-1j * sig * dW + (2j * m * eq.Omega / r) * W) / D
    v1 = ((eq.kappa**2 / (2.0 * eq.Omega)) * dW - (m * sig / r) * W) / D
    return u1, v1, (h1, psi1, psi1_ind)


def energy_budget(problem: EigenProblem, omega: complex, S: np.ndarray) -> EnergyBudget:
    """Evaluate the ARS89 energy budget (eq 33) for a computed mode."""
    model = problem.model
    m = model.m
    r = problem.grid.r
    eq = problem.eq
    ops = problem.ops

    u1, v1, (h1, psi1, psi1_ind) = velocities(problem, omega, S)
    sigma0 = eq.sigma0

    # energy density <E>(r) (spatial part), eq for <E>
    energy = 0.25 * sigma0 * (np.abs(u1) ** 2 + np.abs(v1) ** 2 + np.abs(h1) ** 2 / eq.a0sq)

    # Reynolds stress: -(1/2) r Omega' sigma_0 Re[u1 v1*]
    reynolds = -0.5 * r * eq.dOmega * sigma0 * np.real(u1 * np.conj(v1))

    # acoustic flux: -(1/r) d/dr [ r sigma_0 (1/2) Re[h1 u1*] ]
    flux = r * sigma0 * 0.5 * np.real(h1 * np.conj(u1))
    acoustic = -(1.0 / r) * ops.ddr(flux)

    # gravitational work: -(1/2) sigma_0 Re[ u1 conj(dpsi/dr) + v1 (i m/r) conj(psi) ].
    # The azimuthal average of v1 (1/r) d/dtheta psi yields (1/2) Re[v_hat (i m/r) psi_hat*];
    # the conjugate applies to psi_hat only, not to the (i m/r) factor.
    def grav_work(psi):
        dpsi = ops.ddr(psi)
        return -0.5 * sigma0 * np.real(
            u1 * np.conj(dpsi) + v1 * (1j * m / r) * np.conj(psi)
        )

    grav_direct = grav_work(psi1)
    grav_indirect = grav_work(psi1_ind)

    # global balance: gamma_energy = Integral(work) / (2 Integral(<E>)), dA = 2 pi r dr
    dA = 2.0 * np.pi * r
    work_total = reynolds + acoustic + grav_direct + grav_indirect
    num = np.trapezoid(work_total * dA, r)
    den = 2.0 * np.trapezoid(energy * dA, r)
    gamma_energy = num / den

    return EnergyBudget(
        r=r,
        energy_density=energy,
        reynolds=reynolds,
        acoustic=acoustic,
        grav_direct=grav_direct,
        grav_indirect=grav_indirect,
        u1=u1,
        v1=v1,
        gamma_energy=float(gamma_energy),
        gamma_eigenvalue=float(-omega.imag),
    )
