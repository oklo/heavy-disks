"""Reproduction of the Adams, Ruden & Shu (1989) m=1 SLING eigenvalue code.

Adams, F. C., Ruden, S. P., & Shu, F. H. 1989, ApJ, 347, 959,
"Eccentric Gravitational Instabilities in Nearly Keplerian Disks"
(bibcode 1989ApJ...347..959A).

The package solves the linear normal-mode eigenvalue problem for the m=1
("eccentric") gravitational instability of a razor-thin, nearly-Keplerian,
self-gravitating disk around a central star, including the indirect (stellar
reflex / SLING) potential term.  See README.md for the equation map.
"""

from .model import DiskModel
from .eigensolve import EigenProblem, solve_mode

__all__ = ["DiskModel", "EigenProblem", "solve_mode"]
