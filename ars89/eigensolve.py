"""Assemble the discrete eigenproblem and solve for the complex eigenfrequency.

The reduced equation of motion in matrix form (ARS89 eq B6) is, for interior
rows i (the ODE is *not* imposed at the endpoints; rows 0 and N-1 carry the
boundary conditions, eqs 14a & 16):

    M(omega) S = 0,

with (D1, D2 the log-derivative matrices, J the Poisson matrix, T the indirect
operator, all NxN; A, B, C, nu the omega-dependent coefficient fields):

    M = Bracket_V @ J                              # acts on V = J S
      + diag(1/(Sigma r)) @ Bracket_S              # direct S operator
      + diag(-kappa^2 (1-nu^2) r / (2 pi G sigma_0))   # C h_1 term
      + diag((B r^2 + A r) * (1/2) omega^2 r^3/(G Mtot)) @ T   # indirect (SLING)

    Bracket_V = D2 + diag(A r + 2(1-p) - 1) D1 + diag(B r^2 + A r (1-p) + p(p-1))
    Bracket_S = D2 + diag(A r - 2q - 1)     D1 + diag(B r^2 - A r q   + q(q+1))
    Sigma     = 2 pi G sigma_0 / a_0^2      (eq 30; |k| ~ Sigma)

The omega-dependence is non-linear (nu, nu^2 in A,B,C and omega^2 in the indirect
term), so this is a *non-linear* eigenvalue problem.  We locate the complex omega
that makes M(omega) singular by minimising its smallest singular value, seeded at
the ARS89 Fig. 3 value omega = 4.26 - 0.232 i (in units of Omega_D).  The null
right-singular vector at the solution is the eigenfunction S(r).
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.linalg import svd, lu_factor, lu_solve

from .model import DiskModel
from .discretize import LogGrid, DiffOperators, make_log_grid, make_operators
from .kernel import build_poisson_matrix, build_indirect_operator
from .physics import Equilibrium, build_equilibrium


@dataclass
class EigenProblem:
    """Holds the omega-independent ingredients and assembles M(omega)."""

    model: DiskModel
    grid: LogGrid
    ops: DiffOperators
    eq: Equilibrium
    J: np.ndarray          # Poisson matrix (eq 27/B2)
    T: np.ndarray          # indirect operator (eq B5a)
    indirect_scale: float = 1.0   # diagnostic knob on the SLING indirect term
    selfgrav_scale: float = 1.0   # diagnostic knob on the self-gravity (J) term

    @classmethod
    def build(cls, model: DiskModel, N: int = 1200, eta: float = 0.05) -> "EigenProblem":
        grid = make_log_grid(model.Rstar, model.Rd, N)
        ops = make_operators(grid)
        eq = build_equilibrium(model, grid, ops, eta=eta)
        J = build_poisson_matrix(grid, model.p, model.m)
        T = build_indirect_operator(grid, model.p)
        return cls(model=model, grid=grid, ops=ops, eq=eq, J=J, T=T)

    # --- matrix assembly ---------------------------------------------------

    def _enthalpy_plus_potential_operator(self) -> np.ndarray:
        """Matrix Tmat with (psi_1 + h_1)_i = (Tmat S)_i.

        psi_1 = 2 pi G sigma_0 r V = 2 pi G sigma_0 r (J S)  (eq 25b),
        h_1   = a_0^2 S.
        """
        m = self.model
        r = self.grid.r
        pref = 2.0 * np.pi * m.G * self.eq.sigma0 * r           # 2 pi G sigma_0 r
        return pref[:, None] * self.J + np.diag(self.eq.a0sq)

    def assemble(self, omega: complex) -> np.ndarray:
        """Assemble the full complex matrix M(omega) including BC rows."""
        model = self.model
        r = self.grid.r
        p, q, mm = model.p, model.q, model.m
        D1, D2 = self.ops.D1, self.ops.D2
        eq = self.eq

        A, B, C = eq.coefficients(omega)
        Omega, kappa = eq.Omega, eq.kappa

        # --- V-bracket acting on V = J S ---
        c1_D1 = A * r + 2.0 * (1.0 - p) - 1.0
        c1_D0 = B * r**2 + A * r * (1.0 - p) + p * (p - 1.0)
        Bracket_V = D2 + c1_D1[:, None] * D1 + np.diag(c1_D0)
        term_V = Bracket_V @ (self.selfgrav_scale * self.J)

        # --- direct S-bracket, scaled by 1/(Sigma r) ---
        Sigma = 2.0 * np.pi * model.G * eq.sigma0 / eq.a0sq        # eq 30
        c2_D1 = A * r - 2.0 * q - 1.0
        c2_D0 = B * r**2 - A * r * q + q * (q + 1.0)
        Bracket_S = D2 + c2_D1[:, None] * D1 + np.diag(c2_D0)
        term_S = (1.0 / (Sigma * r))[:, None] * Bracket_S

        # --- C h_1 term: -kappa^2 (1-nu^2) r / (2 pi G sigma_0) ---
        nu = eq.nu(omega)
        c_C = -kappa**2 * (1.0 - nu**2) * r / (2.0 * np.pi * model.G * eq.sigma0)
        term_C = np.diag(c_C)

        # --- indirect (SLING) term: (B r^2 + A r) (1/2) omega^2 r^3/(G Mtot) T ---
        c_ind = (B * r**2 + A * r) * 0.5 * omega**2 * r**3 / (model.G * model.Mtot)
        term_ind = self.indirect_scale * c_ind[:, None] * self.T

        M = term_V + term_S + term_C + term_ind

        # --- boundary condition rows (replace rows 0 and N-1) ---
        Tmat = self._enthalpy_plus_potential_operator()
        N = self.grid.N

        # inner BC (eq 16): (d/dr - 2 m Omega/((omega - m Omega) r))(psi_1+h_1)=0 at R_*
        i = 0
        ddr_T_row = (D1 @ Tmat)[i] / r[i]
        coef_in = 2.0 * mm * Omega[i] / ((omega - mm * Omega[i]) * r[i])
        M[i, :] = ddr_T_row - coef_in * Tmat[i, :]

        # outer BC (eq 14a,14b): sigma_1 = (1/D)(dsigma_0/dr)(d/dr - 2 m Omega/((omega-mOmega)r))(psi_1+h_1)
        j = N - 1
        Dscr = kappa[j] ** 2 - (omega - mm * Omega[j]) ** 2          # eq 14b
        ddr_T_rowN = (D1 @ Tmat)[j] / r[j]
        coef_out = 2.0 * mm * Omega[j] / ((omega - mm * Omega[j]) * r[j])
        rhs_row = (eq.dsigma0[j] / Dscr) * (ddr_T_rowN - coef_out * Tmat[j, :])
        sigma1_row = np.zeros(N, dtype=complex)
        sigma1_row[j] = eq.sigma0[j]
        M[j, :] = sigma1_row - rhs_row

        return M

    # --- residual used by the root finder ---------------------------------

    @staticmethod
    def _row_equilibrate(M: np.ndarray) -> np.ndarray:
        """Divide each row by its norm.

        Left-multiplying by a diagonal does not change the null space, so the
        eigenfrequency (where M is singular) is unchanged, but the smallest
        singular value becomes a clean, well-scaled singularity indicator
        instead of being dominated by the disparate row magnitudes of the
        boundary-condition vs ODE rows across four decades in radius.
        """
        scale = np.linalg.norm(M, axis=1)
        scale[scale == 0] = 1.0
        return M / scale[:, None]

    def smallest_singular_value(self, omega: complex) -> float:
        """sigma_min of the row-equilibrated M(omega); zero at an eigenfrequency."""
        M = self._row_equilibrate(self.assemble(omega))
        s = svd(M, compute_uv=False)
        return s[-1]

    def null_vector(self, omega: complex) -> np.ndarray:
        """Right singular vector for the smallest singular value -> S(r)."""
        M = self._row_equilibrate(self.assemble(omega))
        _, _, Vh = svd(M)
        return Vh[-1].conj()

    def dlogdet(self, omega: complex, h: float = 1e-6) -> complex:
        """d/domega log det M = tr(M^{-1} dM/domega).

        Equals sum_k 1/(omega - omega_k) over the eigenfrequencies omega_k, so a
        Newton step -1/dlogdet drives omega to the nearest det-zero (the actual
        eigenfrequency), following the analytic structure rather than the
        smallest-singular-value landscape (whose minima sit on the neutral
        real-axis modes, not on the shielded growing root).
        """
        M0 = self.assemble(omega)
        Mp = self.assemble(omega + h)
        Mm = self.assemble(omega - h)
        dM = (Mp - Mm) / (2.0 * h)
        lu = lu_factor(M0)
        # tr(M0^{-1} dM) without forming the inverse: solve columns, take trace.
        X = lu_solve(lu, dM)
        return np.trace(X)


# ARS89 Fig. 3 eigenvalue (lowest-order m=1 mode, canonical model), units Omega_D.
ARS89_FIG3_EIGENVALUE = 4.26 - 0.232j


def solve_mode(
    problem: EigenProblem,
    omega_guess: complex = ARS89_FIG3_EIGENVALUE,
    tol: float = 1e-8,
    max_iter: int = 50,
    max_step: float = 0.5,
):
    """Find the complex eigenfrequency near ``omega_guess`` by Newton iteration.

    Uses Newton's method on det M(omega) via its logarithmic derivative,

        omega <- omega - 1 / (d/domega log det M),

    which converges to the nearest det-zero (the eigenfrequency) following the
    analytic structure.  This is essential here: the growing SLING mode is a
    genuine det-zero (confirmed by an argument-principle winding number of 1
    around the ARS89 value) but is *not* a smallest-singular-value minimum,
    because the neutral real-axis cavity modes dominate that landscape.

    Returns
    -------
    omega : complex
        Eigenfrequency in units of Omega_D (= 1 in the default unit system).
    S : np.ndarray
        Eigenfunction S(r) = sigma_1/sigma_0 (complex), normalised to unit max.
    info : dict
        Solver diagnostics (iterations, final |step|, residual sigma_min).
    """
    omega = complex(omega_guess)
    step = np.inf
    for it in range(1, max_iter + 1):
        g = problem.dlogdet(omega)
        if g == 0 or not np.isfinite(g):
            break
        step = -1.0 / g
        # damp overly large steps for global robustness
        if abs(step) > max_step:
            step *= max_step / abs(step)
        omega = omega + step
        if abs(step) < tol:
            break

    S = problem.null_vector(omega)
    S = S / S[np.argmax(np.abs(S))]
    info = {
        "nit": it,
        "last_step": abs(step),
        "residual": problem.smallest_singular_value(omega),
        "converged": abs(step) < tol,
    }
    return omega, S, info
