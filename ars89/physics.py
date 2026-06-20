"""Equilibrium rotation curve and the reduced-ODE coefficients (ARS89).

The equilibrium (Appendix A) builds Omega(r) from three contributions (eq 5a/A1):

    Omega^2(r) = G M_*/r^3 + (1/r) dpsi_disk/dr + (1/r) dh_0/dr
               = Omega_K^2  + Omega_disk^2          + Omega_press^2          (A2,A3,A4)

* Kepler (A2)   : Omega_K^2 = G M_*/r^3.
* Pressure (A4) : (1/r) dh_0/dr = -p a_0^2 / r^2.
* Disk (A3)     : self-gravity, computed from the axisymmetric potential with
                  the eta=0.1 softened kernel (A5); a small, nearly-Keplerian
                  correction except near the outer edge.

The epicyclic frequency is kappa^2 = (1/r^3) d/dr[(r^2 Omega)^2] (eq 5b).

The reduced 2nd-order ODE L(h_1 + psi_1 + psitilde_1) + C h_1 = 0 (eq 21) has
L = d^2/dr^2 + A d/dr + B (eq 9) with coefficients (eqs 10a-d, m=1 written out):

    nu = (omega - m Omega)/kappa
    A  = d/dr log[ sigma_0 r / (kappa^2 (1-nu^2)) ]
       = sigma_0'/sigma_0 + 1/r - 2 kappa'/kappa + 2 nu nu'/(1-nu^2)
    B  = -m^2/r^2
         - (4 m/r^2)(Omega/kappa)(r/(1-nu^2)) nu'
         + (2 m/(r nu))(Omega/kappa) d/dr log(kappa^2/(Omega sigma_0))
    C  = -kappa^2 (1-nu^2)/a_0^2
    nu'= -m Omega'/kappa - nu kappa'/kappa

A, B, C, nu depend on the (complex) eigenfrequency omega; the equilibrium fields
Omega, kappa, sigma_0 and their radial derivatives do not and are precomputed.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .model import DiskModel
from .discretize import LogGrid, DiffOperators
from .kernel import (
    axisym_softened_integral,
    axisym_unsoftened_integral,
    cubic_singular_band_weights,
)


@dataclass
class Equilibrium:
    """Precomputed equilibrium fields on the radial grid."""

    model: DiskModel
    grid: LogGrid
    ops: DiffOperators

    sigma0: np.ndarray
    dsigma0: np.ndarray          # d sigma_0 / dr
    a0sq: np.ndarray
    Omega: np.ndarray
    dOmega: np.ndarray           # d Omega / dr
    kappa: np.ndarray
    dkappa: np.ndarray           # d kappa / dr
    Omega2: np.ndarray
    kappa2: np.ndarray

    # --- omega-dependent coefficient fields -------------------------------

    def nu(self, omega: complex) -> np.ndarray:
        m = self.model.m
        return (omega - m * self.Omega) / self.kappa

    def dnu(self, omega: complex) -> np.ndarray:
        """nu' = -m Omega'/kappa - nu kappa'/kappa."""
        m = self.model.m
        return -m * self.dOmega / self.kappa - self.nu(omega) * self.dkappa / self.kappa

    def coefficients(self, omega: complex):
        """Return (A, B, C) coefficient arrays for the operator L (eqs 10a-c)."""
        m = self.model.m
        r = self.grid.r
        Omega, kappa = self.Omega, self.kappa
        nu = self.nu(omega)
        dnu = self.dnu(omega)
        one_minus_nu2 = 1.0 - nu**2

        dlog_sigma0 = self.dsigma0 / self.sigma0      # sigma_0'/sigma_0
        dlog_kappa = self.dkappa / self.kappa          # kappa'/kappa
        dlog_Omega = self.dOmega / self.Omega          # Omega'/Omega

        A = (
            dlog_sigma0
            + 1.0 / r
            - 2.0 * dlog_kappa
            + 2.0 * nu * dnu / one_minus_nu2
        )

        # d/dr log(kappa^2/(Omega sigma_0)) = 2 kappa'/kappa - Omega'/Omega - sigma_0'/sigma_0
        dlog_ratio = 2.0 * dlog_kappa - dlog_Omega - dlog_sigma0
        B = (
            -(m**2) / r**2
            - (4.0 * m / r**2) * (Omega / kappa) * (r / one_minus_nu2) * dnu
            + (2.0 * m / (r * nu)) * (Omega / kappa) * dlog_ratio
        )

        C = -kappa**2 * one_minus_nu2 / self.a0sq
        return A, B, C


def _disk_self_gravity_potential(model: DiskModel, grid: LogGrid, eta: float) -> np.ndarray:
    r"""Softened axisymmetric disk potential psi_disk(r_i) (ARS89 A3 with A5).

    psi_disk(r) = -(G/r) \int_{R_*}^{R_D} sigma_0(rho) rho Phi(rho/r) drho,
    Phi(x) = \int_0^{2pi} dphi/sqrt(1 + x^2 - 2 x cos phi + eta^2).
    Trapezoidal in rho (smooth integrand for eta>0).
    """
    r = grid.r
    rho = grid.r
    f = np.exp(grid.log_f)

    sigma0_rho = model.sigma0(rho)
    w = rho * 0.5 * (f - 1.0 / f)
    w = w.copy()
    w[0] *= 0.5
    w[-1] *= 0.5

    X = rho[None, :] / r[:, None]
    Phi = axisym_softened_integral(X, eta)
    integrand = sigma0_rho[None, :] * rho[None, :] * Phi * w[None, :]
    return -(model.G / r) * integrand.sum(axis=1)


def _disk_self_gravity_potential_unsoftened(model: DiskModel, grid: LogGrid) -> np.ndarray:
    r"""Unsoftened axisymmetric disk potential (ARS89 A3, interior treatment).

    psi_disk(r) = -G r \int g_axi(s) sigma_0(r e^s) ds,
    g_axi(s) = e^{2s} Phi(e^s),  Phi(x) = \int_0^{2pi} dphi/sqrt(1+x^2-2x cos phi)
                                        = (4/(1+x)) K(k).
    Phi is log-divergent at x=1 (rho=r); the two cells straddling the singularity
    are integrated with the same cubic product-integration weights used for the
    perturbation Poisson kernel (ARS89 Appendix B), the rest by trapezoid in rho.
    """
    r = grid.r
    rho = grid.r
    f = np.exp(grid.log_f)
    log_f = grid.log_f
    N = grid.N

    sigma0_rho = model.sigma0(rho)
    w = rho * 0.5 * (f - 1.0 / f)
    w = w.copy()
    w[0] *= 0.5
    w[-1] *= 0.5

    X = rho[None, :] / r[:, None]
    Phi = axisym_unsoftened_integral(X)
    np.fill_diagonal(Phi, 0.0)                    # singular cells added below
    integrand = sigma0_rho[None, :] * rho[None, :] * Phi * w[None, :]
    psi = -(model.G / r) * integrand.sum(axis=1)

    # cubic-spline-over-singularity correction.  Field interpolated is
    # sigma_0(rho); the prefactor for row i is -G r_i and g_axi carries the
    # kernel times the rho-Jacobian: psi_i += -G r_i sum_d Wcub[d] sigma_0(r_{i+d}).
    g_axi = lambda s: float(axisym_unsoftened_integral(np.exp(s))) * np.exp(2.0 * s)
    Wcub = cubic_singular_band_weights(g_axi, log_f)
    # trapezoidal singular-cell parts already in psi at offsets +-1 (subtract):
    # node i+-1 weight from the singular cell = sigma_0 rho Phi (cell half-width)
    Wtrap = {
        -1: float(axisym_unsoftened_integral(1.0 / f)) * (1.0 / f) * (1.0 - 1.0 / f) / 2.0,
        +1: float(axisym_unsoftened_integral(f)) * f * (f - 1.0) / 2.0,
    }
    # convert the rho-trapezoid (which multiplied sigma_0(rho_j) rho_j Phi w_j and
    # the -(G/r_i) prefactor) into the same per-offset form. In rho-units the
    # left/right singular-cell trapezoidal contribution to psi_i is
    #   -(G/r_i) sigma_0(r_{i+d}) * [rho-cell weight] with rho=r_i x, giving
    #   -(G r_i) * [x Phi(x) (cell width in x)/2] sigma_0  ==  -(G r_i) Wtrap[d] sigma_0.
    delta = {
        -2: -Wcub[-2],
        -1: -Wcub[-1] + Wtrap[-1],
        0: -Wcub[0],
        +1: -Wcub[+1] + Wtrap[+1],
        +2: -Wcub[+2],
    }
    idx = np.arange(N)
    for d, val in delta.items():
        i = idx[max(0, -d):N - max(0, d)]
        psi[i] += -(model.G * r[i]) * val * sigma0_rho[i + d]
    return psi


def _smooth_outer_rotation(model, grid, ops, eta, Omega2_unsoftened):
    """ARS89 (A5) outer smoothing of the rotation curve.

    The unsoftened disk self-gravity integral (A3) is the faithful interior
    rotation curve, but near the disk edges the rho=r singularity sits at the
    truncated integration boundary and differentiating the edge-truncated
    potential becomes noisy (kappa^2 dips below zero in the last few percent of
    the radius).  ARS89 soften (eta=0.1) only to obtain finite endpoint values
    and bridge to them near the edge.

    We realise this by blending the unsoftened interior curve smoothly into the
    eta=0.1 softened curve over a narrow window just outside the mode peak (the
    two differ by <1% there): Omega = (1-w) Omega_unsoft + w Omega_soft, with w a
    smoothstep from 0 at 0.88 R_D to 1 at 0.94 R_D, and w=0 throughout the
    interior.  Both inputs are smooth and finite, so the result is robust at all
    grid resolutions while preserving the unsoftened curve where the mode lives.
    """
    r = grid.r
    Omega_unsoft = np.sqrt(np.clip(Omega2_unsoftened, 0.0, None))

    # softened reference curve (smooth, finite kappa^2 to the edge)
    psi_soft = _disk_self_gravity_potential(model, grid, eta)
    Omega2_soft = (
        model.G * model.Mstar / r**3
        + (1.0 / r) * model.dh0_dr(r)
        + ops.ddr(psi_soft) / r
    )
    Omega_soft = np.sqrt(np.clip(Omega2_soft, 0.0, None))

    # smoothstep blend weight w(r): 0 (unsoftened) -> 1 (softened) over the edge
    x = (r - 0.88 * model.Rd) / (0.94 * model.Rd - 0.88 * model.Rd)
    x = np.clip(x, 0.0, 1.0)
    w = x * x * (3.0 - 2.0 * x)
    Omega = (1.0 - w) * Omega_unsoft + w * Omega_soft
    Omega[0] = Omega_soft[0]            # inner endpoint: nearly-Keplerian softened
    return Omega


def build_equilibrium(
    model: DiskModel,
    grid: LogGrid,
    ops: DiffOperators,
    eta: float = 0.1,
    unsoftened_rotation: bool = False,
) -> Equilibrium:
    """Assemble all equilibrium fields and their radial derivatives.

    ``unsoftened_rotation`` selects ARS89's faithful Appendix-A rotation curve
    (unsoftened disk self-gravity with cubic handling of the rho=r singularity);
    set it False to use the simpler eta-softened integral everywhere.
    """
    r = grid.r
    G, Mstar, p, q = model.G, model.Mstar, model.p, model.q

    sigma0 = model.sigma0(r)
    dsigma0 = model.dsigma0_dr(r)
    a0sq = model.a0_sq(r)

    # --- rotation curve Omega^2 = Kepler + pressure + disk self-gravity ----
    Omega2_kep = G * Mstar / r**3                       # (A2)
    Omega2_press = (1.0 / r) * model.dh0_dr(r)           # (A4) = -p a0^2/r^2
    if unsoftened_rotation:
        psi_disk = _disk_self_gravity_potential_unsoftened(model, grid)  # (A3)
    else:
        psi_disk = _disk_self_gravity_potential(model, grid, eta)        # (A3+A5)
    Omega2_disk = ops.ddr(psi_disk) / r                  # (1/r) dpsi_disk/dr

    Omega2 = Omega2_kep + Omega2_press + Omega2_disk
    if unsoftened_rotation:
        Omega = _smooth_outer_rotation(model, grid, ops, eta, Omega2)
    else:
        Omega = np.sqrt(Omega2)

    # --- epicyclic frequency kappa^2 = (1/r^3) d/dr[(r^2 Omega)^2] (5b) -----
    L = (r**2 * Omega) ** 2                              # = r^4 Omega^2
    kappa2 = ops.ddr(L) / r**3
    # The very outermost points (r -> R_D) can have kappa^2 dip below zero where
    # the rotation curve turns over at the disk edge; floor it to keep kappa real
    # for the outer boundary condition (these points are beyond the mode peak).
    floor = 1e-6 * np.median(kappa2[kappa2 > 0])
    kappa2 = np.where(kappa2 > floor, kappa2, floor)
    kappa = np.sqrt(kappa2)

    # radial derivatives of Omega and kappa (2nd-order FD on the grid)
    dOmega = ops.ddr(Omega)
    dkappa = ops.ddr(kappa)

    return Equilibrium(
        model=model,
        grid=grid,
        ops=ops,
        sigma0=sigma0,
        dsigma0=dsigma0,
        a0sq=a0sq,
        Omega=Omega,
        dOmega=dOmega,
        kappa=kappa,
        dkappa=dkappa,
        Omega2=Omega2,
        kappa2=kappa2,
    )
