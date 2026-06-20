"""Equilibrium disk model and unit system (ARS89 Section III).

All quantities are in a unit system with

    G = M_* = R_D = 1,

so that the natural frequency unit is

    Omega_D = (G M_* / R_D^3)^{1/2} = 1,

which is exactly the unit in which ARS89 quote their eigenvalues (e.g. the
Fig. 3 caption value omega = (4.26, -0.232) Omega_D).  Eigenvalues returned by
the solver are therefore already in units of Omega_D.

Profiles (ARS89 eqs 23a, 23b, 24):

    sigma_0(r) = sigma_*  (R_*/r)^p                       (23a)
    sigma_*    = (2 - p) M_D / (2 pi R_*^2 [(R_D/R_*)^{2-p} - 1])   (23b)
    a_0(r)     = a_{0*} (R_*/r)^{q/2}            (sound speed, index q/2)
    Q_*        = Omega_* a_{0*} / (pi G sigma_*),  Omega_* = (G M_*/R_*^3)^{1/2}  (24)
"""

from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np


@dataclass(frozen=True)
class DiskModel:
    """Parameters of the equilibrium disk plus derived normalisation constants.

    Parameters
    ----------
    p : float
        Surface-density power-law index, sigma_0 ~ r^{-p}.
    q : float
        Temperature power-law index; sound speed a_0 ~ r^{-q/2}.
    Rd_over_Rstar : float
        Ratio of outer disk radius to stellar/inner radius, R_D / R_*.
    Md_over_Mstar : float
        Ratio of disk mass to stellar mass, M_D / M_*.
    Qstar : float
        Toomre-like stability parameter at the inner radius (eq. 24).
    m : int
        Azimuthal wavenumber (this code targets m = 1).
    G, Mstar, Rd : float
        Unit-system anchors; defaults give Omega_D = 1.
    """

    p: float = 1.5
    q: float = 0.5
    Rd_over_Rstar: float = 1.0e4
    Md_over_Mstar: float = 1.0
    Qstar: float = 10.0
    m: int = 1

    G: float = 1.0
    Mstar: float = 1.0
    Rd: float = 1.0

    # Derived (filled in __post_init__ via object.__setattr__ because frozen).
    Rstar: float = field(init=False)
    Md: float = field(init=False)
    Omega_D: float = field(init=False)
    Omega_star: float = field(init=False)
    sigma_star: float = field(init=False)
    a0_star: float = field(init=False)

    def __post_init__(self) -> None:
        Rstar = self.Rd / self.Rd_over_Rstar
        Md = self.Md_over_Mstar * self.Mstar
        Omega_D = np.sqrt(self.G * self.Mstar / self.Rd**3)
        Omega_star = np.sqrt(self.G * self.Mstar / Rstar**3)

        p = self.p
        # eq (23b): normalisation of the surface density from the total disk mass.
        sigma_star = (
            (2.0 - p)
            * Md
            / (2.0 * np.pi * Rstar**2 * (self.Rd_over_Rstar ** (2.0 - p) - 1.0))
        )
        # eq (24) solved for the inner sound speed.
        a0_star = self.Qstar * np.pi * self.G * sigma_star / Omega_star

        for name, value in dict(
            Rstar=Rstar,
            Md=Md,
            Omega_D=Omega_D,
            Omega_star=Omega_star,
            sigma_star=sigma_star,
            a0_star=a0_star,
        ).items():
            object.__setattr__(self, name, float(value))

    # --- equilibrium radial profiles (vectorised over r) --------------------

    def sigma0(self, r: np.ndarray) -> np.ndarray:
        """Unperturbed surface density sigma_0(r), eq (23a)."""
        return self.sigma_star * (self.Rstar / r) ** self.p

    def dsigma0_dr(self, r: np.ndarray) -> np.ndarray:
        """Analytic d sigma_0 / dr."""
        return -self.p * self.sigma0(r) / r

    def a0(self, r: np.ndarray) -> np.ndarray:
        """Sound speed a_0(r) ~ r^{-q/2}."""
        return self.a0_star * (self.Rstar / r) ** (self.q / 2.0)

    def a0_sq(self, r: np.ndarray) -> np.ndarray:
        """Squared sound speed a_0^2(r)."""
        return self.a0(r) ** 2

    def h0(self, r: np.ndarray) -> np.ndarray:
        """Unperturbed enthalpy h_0(r) from dh = a^2 dsigma/sigma (eq 4).

        With a_0^2 ~ r^{-q} and sigma_0 ~ r^{-p}, dh_0/dr = (a_0^2/sigma_0)
        dsigma_0/dr = -p a_0^2 / r, which integrates to h_0 = (p/q) a_0^2 for
        q != 0 (an additive constant is irrelevant to the dynamics).
        """
        if self.q == 0:
            # isothermal: h_0 = a_0^2 ln(sigma_0) + const; only dh_0/dr matters.
            return self.a0_sq(r) * np.log(self.sigma0(r))
        return (self.p / self.q) * self.a0_sq(r)

    def dh0_dr(self, r: np.ndarray) -> np.ndarray:
        """d h_0 / dr = (a_0^2 / sigma_0) dsigma_0/dr = -p a_0^2 / r (eq A4)."""
        return self.a0_sq(r) * self.dsigma0_dr(r) / self.sigma0(r)

    @property
    def Mtot(self) -> float:
        """M_* + M_D, the mass that sets the indirect-term coupling (eq 18b)."""
        return self.Mstar + self.Md
