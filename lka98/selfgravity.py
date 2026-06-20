"""Softened axisymmetric self-gravity of the equilibrium disk (LKA98 App. C).

The equilibrium disk potential (eq 4, axisymmetric part) is

    Psi0(r) = -G int_Rin^RD Sigma0(r') r' dr' int_0^2pi dphi
                  / sqrt(r^2 + r'^2 - 2 r r' cos phi + g^2),   g^2 = eta(r)^2 r^2.

The phi-integral has the closed form (derived here, independent of the paper's
elliptic form C7/C8):

    int_0^2pi dphi / sqrt(A - B cos phi) = (4/sqrt(A+B)) K(k),  k^2 = 2B/(A+B),

with A = r^2 + r'^2 + eta^2 r^2, B = 2 r r', giving

    Iphi(r,r') = 4 / sqrt((r+r')^2 + eta^2 r^2)  *  K(xi^2),
    xi^2 = 4 r r' / ((r+r')^2 + eta^2 r^2).                       (correct modulus)

This is the dimensionally-consistent modulus.  The paper's eq (C8) prints
    xi^2 = 4 r r' / (r + r'^2 + eta^2(r)),
which is dimensionally inconsistent and missing the 2 r r' cross-term; we test
both forms numerically against the stated equilibrium fact Q_min = 1.27 at
r = 0.504 (see selfgravity_check).

scipy.special.ellipk takes the parameter m = k^2.
"""

from __future__ import annotations

import numpy as np
from scipy.special import ellipk

from .model import DiskModel


def phi_integral(r, rp, eta2, modulus: str = "correct"):
    """int_0^2pi dphi / sqrt(r^2 + r'^2 - 2 r r' cos phi + eta^2 r^2).

    modulus="correct":  xi^2 = 4 r r' / ((r+r')^2 + eta^2 r^2)   (derived here)
    modulus="paper_c8": xi^2 = 4 r r' / (r + r'^2 + eta^2)       (literal eq C8)
    The prefactor 4/sqrt(denominator) uses the same denominator as the modulus.
    """
    r = np.asarray(r, dtype=float)
    rp = np.asarray(rp, dtype=float)
    if modulus == "correct":
        denom = (r + rp) ** 2 + eta2 * r**2
    elif modulus == "paper_c8":
        denom = r + rp**2 + eta2
    else:
        raise ValueError(modulus)
    xi2 = 4.0 * r * rp / denom
    return (4.0 / np.sqrt(denom)) * ellipk(xi2)


def disk_potential(model: DiskModel, r_grid: np.ndarray, modulus: str = "correct") -> np.ndarray:
    """Equilibrium disk self-gravity potential Psi0(r) on the grid (eq 4).

    Softening eta is evaluated at the field point r (the paper's g^2(r)).
    Trapezoidal quadrature in r'.
    """
    r = r_grid
    rp = r_grid
    Sig = model.Sigma0(rp)
    eta2 = model.eta2(r)  # field-point softening, shape (N,)

    # (N_field, N_source) matrix of phi-integrals
    R = r[:, None]
    RP = rp[None, :]
    E2 = eta2[:, None]
    Iphi = phi_integral(R, RP, E2, modulus=modulus)  # (N, N)

    integrand = Sig[None, :] * rp[None, :] * Iphi  # Sigma0(r') r' Iphi
    # trapezoidal weights in r'
    w = np.gradient(rp)
    return -model.G * (integrand * w[None, :]).sum(axis=1)


def rotation_curve(model: DiskModel, r_grid: np.ndarray, modulus: str = "correct"):
    """Return (Omega, kappa, Omega2, kappa2) on the grid from radial balance (C1).

    Omega^2(r) = (1/r) d/dr( Psi_disk + Psi* + h0 )
    kappa^2(r) = (1/r^3) d/dr[(r^2 Omega)^2].
    """
    r = r_grid
    Psi_disk = disk_potential(model, r, modulus=modulus)
    dPsi_disk = np.gradient(Psi_disk, r)
    dPsi_star = model.dPsistar_dr(r)
    dh0 = model.dh0(r)

    Omega2 = (dPsi_disk + dPsi_star + dh0) / r
    Omega = np.sqrt(np.clip(Omega2, 0.0, None))

    L = (r**2 * Omega) ** 2
    kappa2 = np.gradient(L, r) / r**3
    kappa = np.sqrt(np.clip(kappa2, 0.0, None))
    return Omega, kappa, Omega2, kappa2


def toomre_Q(model: DiskModel, r_grid: np.ndarray, modulus: str = "correct"):
    """Toomre Q(r) = cs kappa / (pi G Sigma0)."""
    _, kappa, _, _ = rotation_curve(model, r_grid, modulus=modulus)
    cs = np.sqrt(model.cs2(r_grid))
    return cs * kappa / (np.pi * model.G * model.Sigma0(r_grid))
