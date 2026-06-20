"""Radial grid and finite-difference operators (ARS89 Appendix B).

The disk spans four decades in radius (R_D/R_* up to 1e4), so a logarithmic
grid is essential.  ARS89 use a geometric grid (eq B1):

    r_{j+1} = f r_j,   r_1 = R_*,   r_N = R_D,   log f = log(R_D/R_*)/(N-1),

i.e. equally spaced in s = log r.  Derivatives are taken with respect to s and
converted back with (eqs B3a, B3b):

    d/dr      = (1/r)   d/ds
    d^2/dr^2  = (1/r^2) (d^2/ds^2 - d/ds)

The difference matrices D1 = d/ds and D2 = d^2/ds^2 (eqs B4a, B4b) use centred
stencils in the interior and one-sided stencils on the first/last rows.  Those
endpoint rows are never used for the ODE itself; the boundary conditions replace
rows 1 and N of the assembled matrix (see eigensolve.assemble).
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class LogGrid:
    """Geometric radial grid, uniform in s = log r (ARS89 eq B1)."""

    r: np.ndarray          # radii, shape (N,)
    s: np.ndarray          # s = log r, shape (N,)
    log_f: float           # uniform spacing in s
    Rstar: float
    Rd: float

    @property
    def N(self) -> int:
        return self.r.size


def make_log_grid(Rstar: float, Rd: float, N: int) -> LogGrid:
    """Build the geometric grid r_j = R_* f^{j-1}, j = 1..N (eq B1)."""
    s = np.linspace(np.log(Rstar), np.log(Rd), N)
    log_f = (np.log(Rd) - np.log(Rstar)) / (N - 1)
    r = np.exp(s)
    return LogGrid(r=r, s=s, log_f=float(log_f), Rstar=Rstar, Rd=Rd)


def d_ds_matrix(N: int, log_f: float) -> np.ndarray:
    """First log-derivative matrix D1 = d/d(log r), eq (B4a).

    Interior rows: centred (-1, 0, 1)/(2 log f).
    First row:  (-3, 4, -1)/(2 log f);  last row: (1, -4, 3)/(2 log f).
    """
    D1 = np.zeros((N, N))
    h = log_f
    # interior centred difference
    idx = np.arange(1, N - 1)
    D1[idx, idx - 1] = -1.0 / (2 * h)
    D1[idx, idx + 1] = +1.0 / (2 * h)
    # one-sided endpoints (second-order)
    D1[0, 0:3] = np.array([-3.0, 4.0, -1.0]) / (2 * h)
    D1[-1, -3:] = np.array([1.0, -4.0, 3.0]) / (2 * h)
    return D1


def d2_ds2_matrix(N: int, log_f: float) -> np.ndarray:
    """Second log-derivative matrix D2 = d^2/d(log r)^2, eq (B4b).

    Interior rows: centred (1, -2, 1)/(log f)^2.
    First row:  (2, -5, 4, -1)/(log f)^2;  last row: (-1, 4, -5, 2)/(log f)^2.
    """
    D2 = np.zeros((N, N))
    h2 = log_f**2
    idx = np.arange(1, N - 1)
    D2[idx, idx - 1] = 1.0 / h2
    D2[idx, idx] = -2.0 / h2
    D2[idx, idx + 1] = 1.0 / h2
    D2[0, 0:4] = np.array([2.0, -5.0, 4.0, -1.0]) / h2
    D2[-1, -4:] = np.array([-1.0, 4.0, -5.0, 2.0]) / h2
    return D2


@dataclass(frozen=True)
class DiffOperators:
    """Finite-difference operators on a LogGrid.

    D1, D2 act in s = log r.  The physical operators are recovered as
    d/dr = diag(1/r) D1 and d^2/dr^2 = diag(1/r^2)(D2 - D1).
    """

    grid: LogGrid
    D1: np.ndarray
    D2: np.ndarray

    def ddr(self, f: np.ndarray) -> np.ndarray:
        """Physical first derivative df/dr of a field sampled on the grid."""
        return (self.D1 @ f) / self.grid.r

    def d2dr2(self, f: np.ndarray) -> np.ndarray:
        """Physical second derivative d^2 f/dr^2 of a grid field."""
        return (self.D2 @ f - self.D1 @ f) / self.grid.r**2


def make_operators(grid: LogGrid) -> DiffOperators:
    D1 = d_ds_matrix(grid.N, grid.log_f)
    D2 = d2_ds2_matrix(grid.N, grid.log_f)
    return DiffOperators(grid=grid, D1=D1, D2=D2)
