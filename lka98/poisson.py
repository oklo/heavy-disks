"""Perturbed self-gravity for the linear modes (softened m-harmonic, LKA98).

The perturbed potential from a surface-density perturbation sigma1(r) e^{-i m phi}
is the m-th Fourier harmonic of the softened Green's function (eqs 4, C5):

    Psi1(r) = -2 G int_Rin^RD sigma1(r') r' I_m(r,r') dr',
    I_m(r,r') = int_0^pi cos(m phi) dphi
                / sqrt(r^2 + r'^2 - 2 r r' cos phi + eta(r)^2 r^2)
              = (1/r) J_m(x = r'/r, eps = eta(r)^2),
    J_m(x,eps) = int_0^pi cos(m phi) dphi / sqrt(1 + x^2 + eps - 2 x cos phi).

J_m is evaluated in closed form via complete elliptic integrals K, E (scipy uses
the parameter convention ellip*(k^2)).  All closed forms are unit-tested against
direct quadrature in tests/.
"""

from __future__ import annotations

import numpy as np
from scipy.special import ellipk, ellipe

from .model import DiskModel


def softened_Jm(x, eps, m: int):
    r"""J_m(x,eps) = \int_0^pi cos(m phi) dphi / sqrt(1 + x^2 + eps - 2 x cos phi).

    Substituting 1 + x^2 + eps - 2x cos phi = B (1 - k^2 cos^2(phi/2)),
    B = (1+x)^2 + eps, k^2 = 4x/B, gives
        J_m = (2 (-1)^m / sqrt(B)) * C_m(k),
    with C_m = int_0^{pi/2} cos(2 m psi) / sqrt(1 - k^2 sin^2 psi) dpsi expressed
    through K(k), E(k) via the moments M_{2n} = int sin^{2n}psi/sqrt(...) dpsi:
        M0 = K, M2 = (K-E)/k^2, M4 = [(2+k^2)K - 2(1+k^2)E]/(3 k^4).
    """
    x = np.asarray(x, dtype=float)
    B = (1.0 + x) ** 2 + eps
    k2 = 4.0 * x / B
    K = ellipk(k2)
    E = ellipe(k2)
    if m == 0:
        C = K
    elif m == 1:
        C = K - 2.0 * (K - E) / k2
    elif m == 2:
        M4 = ((2.0 + k2) * K - 2.0 * (1.0 + k2) * E) / (3.0 * k2**2)
        C = K - 8.0 * (K - E) / k2 + 8.0 * M4
    else:
        # m>=3 would need higher moments M6, ...; not required for the m=2
        # reference model. (A naive M6 closed form was found to be wrong vs
        # quadrature, so it is intentionally not provided.)
        raise NotImplementedError(f"softened_Jm only implemented for m=0,1,2 (got {m})")
    return (2.0 * ((-1.0) ** m) / np.sqrt(B)) * C


def softened_Jm_quad(x, eps, m: int, nq: int = 400):
    """J_m(x,eps) by direct Gauss-Legendre quadrature on [0,pi] (any m>=0).

    Softening eps>0 removes the phi=0 (x=1) singularity, so plain quadrature is
    accurate.  Validated against the closed forms for m=0,1,2.
    """
    x = np.asarray(x, dtype=float)[..., None]
    epsa = np.asarray(eps, dtype=float)
    epsa = epsa[..., None] if epsa.ndim else epsa
    nodes, w = np.polynomial.legendre.leggauss(nq)        # on [-1,1]
    phi = 0.5 * np.pi * (nodes + 1.0)                      # map to [0,pi]
    wphi = 0.5 * np.pi * w
    denom = np.sqrt(np.maximum(1.0 + x**2 + epsa - 2.0 * x * np.cos(phi), 1e-300))
    return np.sum(np.cos(m * phi) * wphi / denom, axis=-1)


def build_poisson_matrix(model: DiskModel, r: np.ndarray, m: int,
                         quad: bool = False) -> np.ndarray:
    """P such that Psi1 = P @ sigma1 on the grid (softened m-harmonic).

    Psi1(r_i) = -2 G sum_j sigma1_j r_j (1/r_i) J_m(r_j/r_i, eta(r_i)^2) w_j.
    Trapezoidal weights w_j in r'.  quad=True uses Gauss-Legendre J_m (any m).
    """
    N = r.size
    eps = model.eta2(r)               # field-point softening eta(r)^2, shape (N,)
    X = r[None, :] / r[:, None]       # x_ij = r_j / r_i
    Jfun = softened_Jm_quad if quad else softened_Jm
    Jm = Jfun(X, eps[:, None], m)   # (N, N)
    w = np.gradient(r)
    P = -2.0 * model.G * (1.0 / r[:, None]) * r[None, :] * Jm * w[None, :]
    return P
