"""Equilibrium "standard reference" disk of LKA98 (§3, Appendix C).

Heavy polytropic disk: Gaussian surface density (eq 14), polytropic EOS with
gamma_p = gamma = 2, central star m* = 0.6, disk mass m_D = 0.4, G = 1.

    Sigma0(r) = S0 exp[-(r-R0)^2 / w^2]                                   (14)
    cs^2(r)   = K gamma_p Sigma0^{gamma_p-1}                              (12)
    h0(r)     = [gamma_p/(gamma_p-1)] K Sigma0^{gamma_p-1}                (5)
    eta(r)^2  = 0.01 [(r-Rin)/(RD-Rin)]^6 + 1e-4   (softening, C6)        (C6)

Parameters (§3): S0=0.372, R0=0.45, w^2=0.05, Rin=0.25, RD=1.0, K=0.25.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class DiskModel:
    S0: float = 0.372
    R0: float = 0.45
    w2: float = 0.05
    Rin: float = 0.25
    RD: float = 1.0
    mstar: float = 0.6
    mdisk: float = 0.4
    G: float = 1.0
    K: float = 0.25
    gamma_p: float = 2.0
    m: int = 2

    # --- equilibrium profiles (vectorised over r) -------------------------

    def Sigma0(self, r: np.ndarray) -> np.ndarray:
        """Unperturbed surface density (eq 14)."""
        return self.S0 * np.exp(-((r - self.R0) ** 2) / self.w2)

    def dSigma0(self, r: np.ndarray) -> np.ndarray:
        """d Sigma0 / dr."""
        return self.Sigma0(r) * (-2.0 * (r - self.R0) / self.w2)

    def cs2(self, r: np.ndarray) -> np.ndarray:
        """Squared sound speed cs^2 = K gamma_p Sigma0^{gamma_p-1} (eq 12)."""
        return self.K * self.gamma_p * self.Sigma0(r) ** (self.gamma_p - 1.0)

    def h0(self, r: np.ndarray) -> np.ndarray:
        """Equilibrium enthalpy (eq 5)."""
        return (self.gamma_p / (self.gamma_p - 1.0)) * self.K * self.Sigma0(r) ** (
            self.gamma_p - 1.0
        )

    def dh0(self, r: np.ndarray) -> np.ndarray:
        """d h0 / dr = gamma_p K Sigma0^{gamma_p-2} dSigma0 (cf. C3)."""
        return (
            self.gamma_p
            * self.K
            * self.Sigma0(r) ** (self.gamma_p - 2.0)
            * self.dSigma0(r)
        )

    def eta2(self, r: np.ndarray) -> np.ndarray:
        """Softening eta(r)^2 (eq C6)."""
        return 0.01 * ((r - self.Rin) / (self.RD - self.Rin)) ** 6 + 1.0e-4

    def dPsistar_dr(self, r: np.ndarray) -> np.ndarray:
        """Stellar contribution to the rotation curve: dPsi*/dr = G m*/r^2."""
        return self.G * self.mstar / r**2

    def disk_mass(self, n: int = 20000) -> float:
        """Check normalisation: m_D = int 2 pi r Sigma0 dr over [Rin, RD]."""
        r = np.linspace(self.Rin, self.RD, n)
        return float(np.trapezoid(2.0 * np.pi * r * self.Sigma0(r), r))
