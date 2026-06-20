r"""Razor-thin self-gravity kernel and the Poisson / indirect matrices.

Perturbed self-gravity (ARS89 eqs 7d, 7e):

    psi_1(r) = -2 pi G \int K_m(r, rho) sigma_1(rho) d rho,
    K_m(x)   = (x/pi) \int_0^pi cos(m alpha) (1 + x^2 - 2 x cos alpha)^{-1/2} d alpha,
    x = rho / r.

Closed forms (complete elliptic integrals K, E with modulus k, k^2 = 4x/(1+x)^2):

    I_0(x) = \int_0^pi (...)^{-1/2} d alpha               = (2/(1+x)) K(k)
    I_1(x) = \int_0^pi cos alpha (...)^{-1/2} d alpha      = (1+x^2)/(x(1+x)) K(k) - (1+x)/x E(k)
    K_m(x) = (x/pi) I_m(x)

so  K_0(x) = (x/pi)(2/(1+x)) K(k)
    K_1(x) = (1/pi)[ (1+x^2)/(1+x) K(k) - (1+x) E(k) ].

K_m(x) has an integrable logarithmic singularity at x = 1 (rho = r).  Following
ARS89 Appendix B, the cell straddling x = 1 is integrated analytically (here by
adaptive quadrature through the flagged singular point), while the smooth factor
x^{-p} S is held at its on-node value across that narrow cell.

In the reduced variables (ARS89 eqs 25, 27) the Poisson integral becomes the
matrix J with V_i = J_ij S_j:

    V(r) = - \int_{x_*}^{x_D} K_m(x) x^{-p} S(rho = x r) dx.

The indirect (SLING) operator (eq B5a) reuses the same trapezoidal structure on
the smooth integrand x^{2-p}; the omega^2 r^3 prefactor is applied in eigensolve.

scipy's ``ellipk``/``ellipe`` take the parameter ``mm = k^2``.
"""

from __future__ import annotations

import numpy as np
from scipy.special import ellipk, ellipe
from scipy.integrate import quad

from .discretize import LogGrid


# --- kernel closed forms (vectorised over x) -------------------------------

def _k_squared(x: np.ndarray) -> np.ndarray:
    return 4.0 * x / (1.0 + x) ** 2


def K0(x: np.ndarray) -> np.ndarray:
    """Axisymmetric (m=0) kernel K_0(x) = (x/pi)(2/(1+x)) K(k)."""
    x = np.asarray(x, dtype=float)
    mm = _k_squared(x)
    return (x / np.pi) * (2.0 / (1.0 + x)) * ellipk(mm)


def K1(x: np.ndarray) -> np.ndarray:
    """m=1 kernel K_1(x) = (1/pi)[ (1+x^2)/(1+x) K(k) - (1+x) E(k) ]."""
    x = np.asarray(x, dtype=float)
    mm = _k_squared(x)
    return (1.0 / np.pi) * (
        (1.0 + x**2) / (1.0 + x) * ellipk(mm) - (1.0 + x) * ellipe(mm)
    )


def Km(x: np.ndarray, m: int) -> np.ndarray:
    """Kernel K_m for m in {0, 1}."""
    if m == 0:
        return K0(x)
    if m == 1:
        return K1(x)
    raise NotImplementedError(f"closed-form kernel only implemented for m=0,1 (got {m})")


def axisym_unsoftened_integral(x: np.ndarray) -> np.ndarray:
    """\\int_0^{2 pi} d phi / sqrt(1 + x^2 - 2 x cos phi) = (4/(1+x)) K(k).

    The unsoftened axisymmetric kernel for the equilibrium disk self-gravity
    (ARS89 A3); log-divergent at x = 1, handled by cubic product integration.
    """
    x = np.asarray(x, dtype=float)
    mm = _k_squared(x)
    return (4.0 / (1.0 + x)) * ellipk(mm)


def axisym_softened_integral(x: np.ndarray, eta: float) -> np.ndarray:
    r"""\int_0^{2 pi} d phi / sqrt(1 + x^2 - 2 x cos phi + eta^2)  (ARS89 eq A5).

    Used only for the *equilibrium* rotation-curve self-gravity (Appendix A),
    where ARS89 soften with eta = 0.1.  Closed form:

        = 2 * (2 / sqrt((1+x)^2 + eta^2)) K(k_eta),   k_eta^2 = 4x/((1+x)^2+eta^2).
    """
    x = np.asarray(x, dtype=float)
    denom = (1.0 + x) ** 2 + eta**2
    mm = 4.0 * x / denom
    return 2.0 * (2.0 / np.sqrt(denom)) * ellipk(mm)


# --- cubic-spline-over-singularity weights (ARS89 Appendix B) ---------------

def _lagrange(nodes, j, s):
    """Cubic Lagrange basis polynomial l_j(s) for the given nodes."""
    num = 1.0
    den = 1.0
    for k, nk in enumerate(nodes):
        if k == j:
            continue
        num = num * (s - nk)
        den = den * (nodes[j] - nk)
    return num / den


def cubic_singular_band_weights(g, log_f: float):
    r"""Cubic product-integration weights through the x=1 log singularity.

    ARS89 (Appendix B) integrate the two cells straddling x=1 (the intervals
    [1/f, 1] and [1, f]) by representing the smooth factor as a cubic spline and
    integrating the singular kernel analytically.  Working in s = log x (uniform
    spacing ``h = log_f``), the smooth field is interpolated by a local cubic and

        \int_cell g(s) field(s) ds = sum_d w[d] field_{i+d},
        w[d] = \int_cell g(s) l_d(s) ds,

    with ``g`` carrying the (log-singular) kernel times Jacobian.  The two cells
    use the cubic stencils {-2,-1,0,1} and {-1,0,1,2}; the returned weights are
    summed over both cells, are translation-invariant in i, and are computed once.
    """
    weights = {d: 0.0 for d in (-2, -1, 0, 1, 2)}
    h = log_f
    left_nodes = [-2 * h, -h, 0.0, h]
    left_off = [-2, -1, 0, 1]
    right_nodes = [-h, 0.0, h, 2 * h]
    right_off = [-1, 0, 1, 2]
    for nodes, offs, lo, hi in ((left_nodes, left_off, -h, 0.0),
                                (right_nodes, right_off, 0.0, h)):
        for j, off in enumerate(offs):
            val, _ = quad(lambda s: g(s) * _lagrange(nodes, j, s), lo, hi,
                          points=[0.0], limit=200)
            weights[off] += val
    return weights


# --- Poisson / indirect matrix builders ------------------------------------
#
# Trapezoidal rule in x on the geometric grid: the weight for node j is
# (x_{j+1} - x_{j-1})/2 = x_ij (f - 1/f)/2, with half-weights at the global ends.

def build_poisson_matrix(grid: LogGrid, p: float, m: int) -> np.ndarray:
    r"""Assemble J such that V_i = sum_j J_ij S_j  (ARS89 eqs 27, B2).

    V(r) = - \int K_m(x) x^{-p} S dx, trapezoidal in x on the geometric grid,
    with the x=1 cell replaced by the analytic self-term.
    """
    N = grid.N
    log_f = grid.log_f
    f = np.exp(log_f)

    # x_ij = r_j / r_i = f^{(j-i)}; build the offset matrix once.
    j_minus_i = np.arange(N)[None, :] - np.arange(N)[:, None]
    X = f ** j_minus_i  # (N, N) ratios r_j/r_i

    # trapezoidal weight in x: w_ij = x_ij (f - 1/f)/2
    half_band = 0.5 * (f - 1.0 / f)
    W = X * half_band
    W[:, 0] *= 0.5
    W[:, -1] *= 0.5

    # far-field: trapezoidal with the smooth kernel.  Zero the diagonal kernel
    # (x=1) so the singular cells contribute nothing yet; they are added back
    # below by the cubic product-integration correction.
    Kvals = Km(X, m)
    np.fill_diagonal(Kvals, 0.0)
    J = -(Kvals * X ** (-p) * W)

    # cubic-spline-over-singularity correction on the two cells around x=1.
    # Smooth factor interpolated by the cubic is S; fold x^{-p} and the Jacobian
    # e^s into the kernel weight g(s) = K_m(e^s) e^{(1-p)s}.
    g = lambda s: float(Km(np.exp(s), m)) * np.exp((1.0 - p) * s)
    Wcub = cubic_singular_band_weights(g, log_f)
    # trapezoidal singular-cell parts already present in J at d=+-1 (to subtract)
    Wtrap = {
        -1: float(Km(1.0 / f, m)) * f**p * (1.0 - 1.0 / f) / 2.0,
        +1: float(Km(f, m)) * f ** (-p) * (f - 1.0) / 2.0,
    }
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
        J[i, i + d] += val

    return J


def build_indirect_operator(grid: LogGrid, p: float) -> np.ndarray:
    r"""Assemble T such that \int_{x_*}^{x_D} x^{2-p} S dx = sum_j T_ij S_j (eq B5a).

    Smooth integrand -> plain trapezoidal rule in x.  The omega^2 r_i^3 /
    (G (M_*+M_D)) / 2 prefactor is applied per row in eigensolve.assemble.
    """
    N = grid.N
    log_f = grid.log_f
    f = np.exp(log_f)

    j_minus_i = np.arange(N)[None, :] - np.arange(N)[:, None]
    X = f ** j_minus_i
    half_band = 0.5 * (f - 1.0 / f)
    W = X * half_band
    W[:, 0] *= 0.5
    W[:, -1] *= 0.5

    return X ** (2.0 - p) * W
