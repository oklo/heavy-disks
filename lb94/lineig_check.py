#!/usr/bin/env python3
"""Verification pass for the LB94 linear eigenproblem:

(1) Argument-principle root COUNT: the number of zeros of det M(omega) inside a rectangle
    is the winding number (1/2 pi i) oint dlogdet domega. This is exhaustive -- if the
    sigma_min hunt missed a fast-growing root, the count exposes it. Run on the full
    window and on the "deep" sub-window (gamma_code > 0.3, i.e. > 0.5/T_unit) where the
    SPH-measured growth (~2.5 code) would live.
(2) Regression at higher N: ars89 vs generalized eigenvalue on the canonical power law.
(3) Sensitivity of the found modes: N, domain boundaries, thickness (softened kernel).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np

from lb94.lineig import (FamilyModel, GeneralEigenProblem, solve_mode, CONV, R_U_AU)
from ars89.model import DiskModel
from ars89.eigensolve import EigenProblem


def winding(problem, re0, re1, im0, im1, n_edge=70):
    """(1/2 pi i) oint dlogdet domega around the rectangle; ~integer = #roots inside."""
    corners = [re0 + 1j * im0, re1 + 1j * im0, re1 + 1j * im1, re0 + 1j * im1]
    total = 0.0 + 0.0j
    for k in range(4):
        a, b = corners[k], corners[(k + 1) % 4]
        ts = np.linspace(0.0, 1.0, n_edge, endpoint=False) + 0.5 / n_edge
        dz = (b - a) / n_edge
        for t in ts:
            total += problem.dlogdet(a + t * (b - a)) * dz
    return total / (2j * np.pi)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"

    if which in ("all", "wind"):
        fam = FamilyModel(m=1)
        pr = GeneralEigenProblem.build_family(fam, N=600, eta_eq=0.1)
        print("=== argument-principle root counts, m=1 (N=600) ===")
        for (re0, re1, im0, im1, tag) in [
                (0.15, 4.0, -3.0, -0.30, "deep:  gamma_code > 0.3 (SPH-like rates)"),
                (0.15, 4.0, -0.30, -0.04, "shallow: 0.04 < gamma_code < 0.3")]:
            w = winding(pr, re0, re1, im0, im1)
            print(f"  [{re0},{re1}] x [{im0},{im1}]i : winding = {w.real:+.3f}{w.imag:+.3f}i   ({tag})")

    if which in ("all", "reg"):
        print("\n=== regression at N=800 (canonical ARS89 power law) ===")
        dm = DiskModel(p=1.5, q=0.5, Rd_over_Rstar=1e4, Md_over_Mstar=1.0, Qstar=10.0, m=1)
        orig = EigenProblem.build(dm, N=800)
        w0, _, i0 = solve_mode(orig, 4.26 - 0.232j)

        class PL(FamilyModel):
            def __init__(s):
                s.p, s.q, s.m = dm.p, dm.q, 1
                s.Rin, s.Rout = 0.0, np.inf
                s.Sigma0 = dm.sigma_star * dm.Rstar ** dm.p
                s.cs0 = dm.a0_star * dm.Rstar ** (dm.q / 2.0)
                s.G, s.Mstar = dm.G, dm.Mstar
                s.Rstar, s.Rd = dm.Rstar, dm.Rd
                s.Md, s.Mtot = dm.Md, dm.Mtot

        gen = GeneralEigenProblem.build_family(PL(), N=800)
        w1, _, i1 = solve_mode(gen, 4.26 - 0.232j)
        print(f"  ars89:       {w0:.5f}  (res {i0['residual']:.1e})")
        print(f"  generalized: {w1:.5f}  (res {i1['residual']:.1e})")
        print(f"  |diff| = {abs(w1-w0):.2e}   ARS89 Fig3: 4.26-0.232j")

    if which in ("all", "sens"):
        print("\n=== sensitivity of the leading m=1 and m=2 modes ===")
        seeds = {1: 1.8988 - 0.0215j, 2: 6.7919 - 0.4793j}
        base = {}
        for m in (1, 2):
            for tag, kw in [("N=700 base", dict(N=700)),
                            ("N=1000", dict(N=1000)),
                            ("rin 6 AU", dict(N=700, rmin=0.06)),
                            ("rout 450 AU", dict(N=700, rmax=4.5)),
                            ("thick eta_p=0.25", dict(N=700, eta_pert=0.25))]:
                mk = {k: v for k, v in kw.items() if k in ("rmin", "rmax")}
                fk = {k: v for k, v in kw.items() if k in ("N", "eta_pert")}
                fam = FamilyModel(m=m, **mk)
                pr = GeneralEigenProblem.build_family(fam, eta_eq=0.1, **fk)
                w, _, info = solve_mode(pr, seeds[m], max_iter=60)
                ok = info["converged"] and info["residual"] < 1e-6
                print(f"  m={m} {tag:18s}: omega = {w.real:7.4f}{w.imag:+8.4f}i  "
                      f"gamma = {-w.imag*CONV:5.2f}/Tu  {'' if ok else '(NOT converged)'}")
