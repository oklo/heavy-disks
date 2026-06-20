"""Linear global m=2 spiral eigenmode: matrix method (LKA98 §4).

We discretise the reduced first-order system (eqs 26, 27) for the state
(sigma1, u1) on a radial grid, with the perturbed potential Psi1 = P sigma1 a
dense nonlocal coupling (softened m-harmonic, poisson.py), and reflective
boundary conditions u1 = 0 at r = Rin and r = RD.  This gives a 2N x 2N complex
matrix M(omega), nonlinear in omega (through nu, A), whose singularity locates
the eigenfrequency omega = m Omega_p - i gamma1 (eq 22).

Reduced system (sigma1' , u1' written with A = i kappa nu, nu = (omega-mOmega)/kappa,
W = (cs^2/Sigma0) sigma1 + Psi1):
    sigma1' = a11 sigma1 + a12 u1 + b1 Psi1 + c1 Psi1'
    u1'     = a21 sigma1 + a22 u1 + b2 Psi1
with coefficients (derived in EQUATIONS.md):
    a11 = 2 Omega m/(nu kappa r) - (Sigma0/cs^2) d/dr(cs^2/Sigma0)
    a12 = -A Sigma0 (nu^2-1)/(cs^2 nu^2)
    b1  = 2 Omega Sigma0 m/(cs^2 nu kappa r)
    c1  = -Sigma0/cs^2
    a21 = i m^2 cs^2/(nu r^2 kappa Sigma0) - A/Sigma0
    a22 = i m A/(2 nu^2 Omega r) - (1/(Sigma0 r)) d/dr(r Sigma0)
    b2  = i m^2/(nu r^2 kappa)

The eigenfrequency is found by Newton iteration on d/domega log det M (the same
robust det-root method used for ARS89, since the growing mode is a det zero, not
a smallest-singular-value minimum).
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.linalg import lu_factor, lu_solve, svd

from .model import DiskModel
from .selfgravity import rotation_curve
from .poisson import build_poisson_matrix


def derivative_matrix(r: np.ndarray) -> np.ndarray:
    """First-derivative matrix d/dr on a (possibly non-uniform) grid.

    Centred differences in the interior, second-order one-sided at the ends.
    """
    N = r.size
    D = np.zeros((N, N))
    for i in range(1, N - 1):
        h1 = r[i] - r[i - 1]
        h2 = r[i + 1] - r[i]
        D[i, i - 1] = -h2 / (h1 * (h1 + h2))
        D[i, i] = (h2 - h1) / (h1 * h2)
        D[i, i + 1] = h1 / (h2 * (h1 + h2))
    # one-sided second-order at the ends
    h = r[1] - r[0]
    D[0, 0:3] = np.array([-3.0, 4.0, -1.0]) / (2 * h)
    h = r[-1] - r[-2]
    D[-1, -3:] = np.array([1.0, -4.0, 3.0]) / (2 * h)
    return D


@dataclass
class LinearProblem:
    model: DiskModel
    r: np.ndarray
    D: np.ndarray            # d/dr matrix
    P: np.ndarray            # Poisson matrix (Psi1 = P sigma1)
    Omega: np.ndarray
    kappa: np.ndarray
    cs2: np.ndarray
    Sigma0: np.ndarray
    dcs2_over_Sigma0: np.ndarray   # d/dr(cs^2/Sigma0)
    drSigma0: np.ndarray           # d/dr(r Sigma0)

    @classmethod
    def build(cls, model: DiskModel, N: int = 600, modulus: str = "correct") -> "LinearProblem":
        r = np.linspace(model.Rin, model.RD, N)
        D = derivative_matrix(r)
        P = build_poisson_matrix(model, r, model.m)
        Omega, kappa, _, _ = rotation_curve(model, r, modulus=modulus)
        cs2 = model.cs2(r)
        Sigma0 = model.Sigma0(r)
        dcs2_over_Sigma0 = D @ (cs2 / Sigma0)
        drSigma0 = D @ (r * Sigma0)
        return cls(model, r, D, P, Omega, kappa, cs2, Sigma0, dcs2_over_Sigma0, drSigma0)

    def assemble(self, omega: complex) -> np.ndarray:
        m = self.model.m
        r = self.r
        Om, kap = self.Omega, self.kappa
        cs2, Sig = self.cs2, self.Sigma0
        N = r.size

        nu = (omega - m * Om) / kap
        A = 1j * kap * nu                                   # = i(omega - m Omega)
        gamma_sigma = (Sig / cs2) * self.dcs2_over_Sigma0

        a11 = 2.0 * Om * m / (nu * kap * r) - gamma_sigma
        a12 = -A * Sig * (nu**2 - 1.0) / (cs2 * nu**2)
        b1 = 2.0 * Om * Sig * m / (cs2 * nu * kap * r)
        c1 = -Sig / cs2
        a21 = 1j * m**2 * cs2 / (nu * r**2 * kap * Sig) - A / Sig
        a22 = 1j * m * A / (2.0 * nu**2 * Om * r) - self.drSigma0 / (Sig * r)
        b2 = 1j * m**2 / (nu * r**2 * kap)

        D, P = self.D, self.P
        DP = D @ P

        def dg(v):
            return v[:, None] * np.eye(N)  # diag

        # sigma-block: D sigma - a11 sigma - a12 u - b1 (P sigma) - c1 (DP sigma) = 0
        M11 = D - dg(a11) - dg(b1) @ P - dg(c1) @ DP
        M12 = -dg(a12)
        # u-block: D u - a21 sigma - a22 u - b2 (P sigma) = 0
        M21 = -dg(a21) - dg(b2) @ P
        M22 = D - dg(a22)

        M = np.block([[M11, M12], [M21, M22]]).astype(complex)

        # reflective BCs: replace u-ODE rows at the two edges with u=0
        for i in (0, N - 1):
            row = N + i
            M[row, :] = 0.0
            M[row, N + i] = 1.0
        return M

    # --- det-root solver (Newton on d/domega log det M) -------------------

    def dlogdet(self, omega: complex, h: float = 1e-6) -> complex:
        M0 = self.assemble(omega)
        Mp = self.assemble(omega + h)
        Mm = self.assemble(omega - h)
        dM = (Mp - Mm) / (2.0 * h)
        lu = lu_factor(M0)
        return np.trace(lu_solve(lu, dM))

    def smallest_singular_value(self, omega: complex) -> float:
        M = self.assemble(omega)
        scale = np.linalg.norm(M, axis=1)
        scale[scale == 0] = 1.0
        return svd(M / scale[:, None], compute_uv=False)[-1]

    def eigenvector(self, omega: complex):
        M = self.assemble(omega)
        scale = np.linalg.norm(M, axis=1)
        scale[scale == 0] = 1.0
        _, _, Vh = svd(M / scale[:, None])
        x = Vh[-1].conj()
        N = self.r.size
        return x[:N], x[N:]      # sigma1, u1


# eigenvalue convention omega = m Omega_p - i gamma1 (eq 22); standard disk m=2.
LKA98_MATRIX_EIGENVALUE = 2 * 1.70 - 0.804j


def contour_root(problem: LinearProblem, center: complex, hr: float, hi: float,
                 n: int = 160) -> complex:
    """Locate the single eigenfrequency inside a rectangular contour.

    With one analytic zero of det M(omega) inside (winding number 1) and no poles
    (resonances are off-axis for a complex, growing mode), the argument principle
    gives the root directly,
        omega_root = (1/2pi i) oint omega d(log det) / (1/2pi i) oint d(log det),
    which avoids the Newton drift toward the dense neutral real-axis modes.
    """
    re = [center.real - hr, center.real + hr]
    im = [center.imag - hi, center.imag + hi]
    corners = [re[0] + 1j * im[0], re[1] + 1j * im[0],
               re[1] + 1j * im[1], re[0] + 1j * im[1]]
    num = 0.0 + 0j
    den = 0.0 + 0j
    for k in range(4):
        a, b = corners[k], corners[(k + 1) % 4]
        t = (np.arange(n) + 0.5) / n
        w = a + (b - a) * t
        dw = (b - a) / n
        g = np.array([problem.dlogdet(wi) for wi in w])  # d/domega log det
        num += np.sum(w * g) * dw
        den += np.sum(g) * dw
    return num / den


def solve_mode(problem: LinearProblem, omega_guess: complex = LKA98_MATRIX_EIGENVALUE,
               tol: float = 1e-9, max_iter: int = 60, max_step: float = 0.5,
               localize: bool = True, box=(0.35, 0.35)):
    """Find the eigenfrequency: contour-localize the growing root, then Newton-refine."""
    if localize:
        omega = contour_root(problem, complex(omega_guess), box[0], box[1])
    else:
        omega = complex(omega_guess)
    step = np.inf
    for it in range(1, max_iter + 1):
        g = problem.dlogdet(omega)
        if g == 0 or not np.isfinite(g):
            break
        step = -1.0 / g
        if abs(step) > max_step:
            step *= max_step / abs(step)
        omega = omega + step
        if abs(step) < tol:
            break
    sigma1, u1 = problem.eigenvector(omega)
    m = problem.model.m
    info = {
        "nit": it,
        "last_step": abs(step),
        "residual": problem.smallest_singular_value(omega),
        "Omega_p": omega.real / m,
        "gamma1": -omega.imag,
    }
    return omega, sigma1, u1, info
