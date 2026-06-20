"""Reproduction of Laughlin, Korchagin & Adams (1998), ApJ 504, 945,
"The Dynamics of Heavy Gaseous Disks".

Linear global spiral-mode analysis of a heavy, polytropic, self-gravitating
disk (the m=2-unstable "standard reference model"), plus the second/third-order
weakly-nonlinear theory and Athena++ hydro comparison.  See EQUATIONS.md for the
transcribed/verified equation set and validation targets.
"""

from .model import DiskModel

__all__ = ["DiskModel"]
