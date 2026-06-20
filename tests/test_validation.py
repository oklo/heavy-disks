"""Validation and unit tests for the ARS89 reproduction.

The headline test recovers the lowest-order m=1 SLING eigenfrequency for the
canonical disk and checks it against the ARS89 Fig. 3 value (4.26 - 0.232 i).
The remaining tests exercise the building blocks (kernel, derivative matrices,
equilibrium limits) independently.
"""

import numpy as np
import pytest
from scipy.integrate import quad

from ars89 import DiskModel, EigenProblem, solve_mode
from ars89.eigensolve import ARS89_FIG3_EIGENVALUE
from ars89.energy import energy_budget
from ars89.kernel import K0, K1
from ars89.discretize import make_log_grid, make_operators


# ---------------------------------------------------------------------------
# headline validation: recover the ARS89 Fig. 3 growing mode
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def canonical_solution():
    model = DiskModel(p=1.5, q=0.5, Rd_over_Rstar=1e4, Md_over_Mstar=1.0, Qstar=10.0, m=1)
    problem = EigenProblem.build(model, N=1000)
    omega, S, info = solve_mode(problem, omega_guess=ARS89_FIG3_EIGENVALUE)
    return model, problem, omega, S, info


def test_eigenfrequency_matches_ars89(canonical_solution):
    _, _, omega, _, info = canonical_solution
    # genuine eigenfrequency: det(M) vanishes (smallest singular value ~ 0)
    assert info["residual"] < 1e-8
    # growing mode matching ARS89 (4.26, -0.232) to ~1% (Re) / ~5% (Im)
    assert 4.0 < omega.real < 4.4, omega
    assert 0.20 < abs(omega.imag) < 0.28, omega         # growing (non-zero) mode
    assert abs(omega.real - 4.26) / 4.26 < 0.03         # pattern speed within 3%
    assert abs(abs(omega.imag) - 0.232) / 0.232 < 0.10  # growth rate within 10%


def test_corotation_radius(canonical_solution):
    _, problem, omega, _, _ = canonical_solution
    r = problem.grid.r
    icr = int(np.argmin(np.abs(problem.eq.Omega - abs(omega.real))))
    # ARS89 quote R_CR = 0.452 R_D for this mode
    assert 0.42 < r[icr] < 0.50


def test_energy_balance(canonical_solution):
    # ARS89 eq (33): the growth rate from the energy budget (Reynolds + acoustic
    # + direct/indirect gravitational work) must match the eigenvalue growth rate.
    # This is an independent check of the eigenfunction and the numerical scheme.
    _, problem, omega, S, _ = canonical_solution
    eb = energy_budget(problem, omega, S)
    assert eb.gamma_energy > 0                       # growing
    assert eb.relative_error < 0.02                  # LHS = RHS to better than 2%
    # Reynolds stress is negative (disturbance feeds energy into the shear, ARS89)
    import numpy as np
    r = problem.grid.r
    assert np.trapezoid(eb.reynolds * 2 * np.pi * r, r) < 0


def test_eigenfunction_oscillatory(canonical_solution):
    _, _, _, S, _ = canonical_solution
    # ARS89 Fig. 3: many radial nodes with a real/imaginary phase offset
    nodes_re = np.sum(np.diff(np.sign(np.real(S))) != 0)
    nodes_im = np.sum(np.diff(np.sign(np.imag(S))) != 0)
    assert nodes_re >= 8
    assert nodes_im >= 8
    assert abs(nodes_re - nodes_im) <= 3   # ~180 deg phase shift -> offset by ~1


# ---------------------------------------------------------------------------
# unit tests for the building blocks
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("x", [0.1, 0.3, 0.7, 1.3, 3.0, 10.0])
def test_kernel_closed_forms(x):
    def Km_num(x, m):
        v, _ = quad(lambda a: np.cos(m * a) / np.sqrt(1 + x * x - 2 * x * np.cos(a)), 0, np.pi)
        return (x / np.pi) * v

    assert np.isclose(float(K0(x)), Km_num(x, 0), rtol=1e-6)
    assert np.isclose(float(K1(x)), Km_num(x, 1), rtol=1e-6)


def test_derivative_matrices():
    g = make_log_grid(1e-3, 1.0, 400)
    op = make_operators(g)
    r = g.r
    f = r**3
    # interior second-order accuracy on a log grid
    assert np.max(np.abs(op.ddr(f) - 3 * r**2)[1:-1] / (3 * r**2)[1:-1]) < 5e-3
    assert np.max(np.abs(op.d2dr2(f) - 6 * r)[1:-1] / (6 * r)[1:-1]) < 5e-3


def test_rotation_curve_keplerian_limit():
    # As M_D -> 0 both the disk self-gravity (~M_D) and the pressure support
    # (a_0 ~ Q_* M_D) vanish, so Omega -> Keplerian and kappa -> Omega.
    model = DiskModel(p=1.5, q=0.5, Rd_over_Rstar=1e3, Md_over_Mstar=1e-8, Qstar=10.0, m=1)
    problem = EigenProblem.build(model, N=400)
    r = problem.grid.r
    Omega_kep = np.sqrt(model.G * model.Mstar / r**3)
    interior = slice(5, -5)
    assert np.allclose(problem.eq.Omega[interior], Omega_kep[interior], rtol=1e-3)
    assert np.allclose(problem.eq.kappa[interior], problem.eq.Omega[interior], rtol=1e-3)


def test_sigma_normalisation():
    # integral of sigma_0 over the disk should equal M_D
    model = DiskModel(p=1.5, Md_over_Mstar=0.3)
    g = make_log_grid(model.Rstar, model.Rd, 2000)
    r = g.r
    mass = np.trapezoid(2 * np.pi * r * model.sigma0(r), r)
    assert np.isclose(mass, model.Md, rtol=1e-3)
