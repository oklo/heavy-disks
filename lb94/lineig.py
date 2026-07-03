#!/usr/bin/env python3
"""Linear global modes (m=1, m=2) of the pinned LB94 basic-state family, reusing the
ars89 eigenvalue machinery generalized from pure power-law Sigma to the tapered family

    Sigma(r) = Sigma0 (r/R0)^(-p) exp[-(Rin/r) - (r/Rout)^2],   c_s ~ r^(-q/2), q = 1/2.

Generalizations vs ars89 (which hardcodes sigma0 ~ r^-p):
  * V-bracket constants (1-p) and p(p-1) -> local fields LP1 = 1 + dln(sigma0)/dlnr and
    LP2 = LP1^2 - LP1 + dLP1/dlnr  (from r P'/P and r^2 P''/P with P = 2 pi G sigma0 r).
  * Poisson matrix weight x^-p -> sigma0(r_j)/sigma0(r_i); singular-band cubic weights
    interpolated over the LOCAL slope p_loc(r_i) = -dln(sigma0)/dlnr.
  * Indirect operator weight x^(2-p) -> x^2 sigma0(r_j)/sigma0(r_i).
  * S-bracket unchanged (c_s is a pure q=1/2 power law); BC rows unchanged (numeric fields).

With tapers off (Rin=0, Rout=inf) everything must reduce EXACTLY to ars89 -- the
regression test builds both matrices on the same grid and diffs them.

Units: G = M_* = 1, R_u = 100 AU  =>  t_u = sqrt(R_u^3/G M_*);  conversions to the SPH
T_unit (474 yr) printed with the results. ARS89 convention omega = m Omega_p - i gamma
(growing modes have Im omega < 0).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.special import ellipk, ellipe
from scipy.integrate import quad

from ars89.discretize import make_log_grid, make_operators
from ars89.kernel import Km, cubic_singular_band_weights, build_poisson_matrix, \
    build_indirect_operator
from ars89.physics import build_equilibrium
from ars89.eigensolve import EigenProblem, solve_mode
from ars89.model import DiskModel

# ---------------- physical constants / unit system ----------------
AUcm, Msun_g, G_cgs, yr_s = 1.496e13, 1.989e33, 6.674e-8, 3.156e7
MSTAR_MSUN = 0.340
R_U_AU = 100.0
SIG_U = MSTAR_MSUN * Msun_g / (R_U_AU * AUcm) ** 2          # surface-density unit, g/cm^2
V_U = np.sqrt(G_cgs * MSTAR_MSUN * Msun_g / (R_U_AU * AUcm))  # cm/s
T_U_YR = (R_U_AU * AUcm / V_U) / yr_s                        # code time unit in years
SPH_TU_YR = 474.0
CONV = SPH_TU_YR / T_U_YR                                    # omega_code -> omega per SPH T_unit


class FamilyModel:
    """LB94 basic-state family; duck-typed to what ars89.build_equilibrium expects."""

    def __init__(self, Sigma0_cgs=48.4, p=2.51, Rin_AU=54.2, Rout_AU=267.0,
                 cs0_kms=0.534, q=0.5, m=1, rmin=0.04, rmax=3.5):
        self.Sigma0 = Sigma0_cgs / SIG_U
        self.p, self.q, self.m = p, q, m
        self.Rin = Rin_AU / R_U_AU
        self.Rout = Rout_AU / R_U_AU
        self.cs0 = cs0_kms * 1e5 / V_U
        self.G, self.Mstar = 1.0, 1.0
        self.Rstar, self.Rd = rmin, rmax                      # grid anchors (deep in tapers)
        # disk mass by quadrature (wide domain)
        rq = np.linspace(1e-3, 10.0, 200000)
        self.Md = float(np.trapezoid(2 * np.pi * rq * self.sigma0(rq), rq))
        self.Mtot = self.Mstar + self.Md

    # --- profiles ---
    def sigma0(self, r):
        r = np.asarray(r, dtype=float)
        out = self.Sigma0 * r ** (-self.p) * np.exp(-self._rin_over(r) - (r / self.Rout) ** 2)
        return out

    def _rin_over(self, r):
        return self.Rin / r if self.Rin > 0 else 0.0 * r

    def dln_sigma_dlnr(self, r):
        r = np.asarray(r, dtype=float)
        return -self.p + self.Rin / r - 2.0 * (r / self.Rout) ** 2

    def d_dlnsig_dlnr(self, r):
        r = np.asarray(r, dtype=float)
        return -self.Rin / r - 4.0 * (r / self.Rout) ** 2

    def dsigma0_dr(self, r):
        return self.sigma0(r) * self.dln_sigma_dlnr(r) / np.asarray(r, dtype=float)

    def a0(self, r):
        return self.cs0 * np.asarray(r, dtype=float) ** (-self.q / 2.0)

    def a0_sq(self, r):
        return self.a0(r) ** 2

    def dh0_dr(self, r):
        return self.a0_sq(r) * self.dln_sigma_dlnr(r) / np.asarray(r, dtype=float)


# ---------------- generalized kernel matrices ----------------

def _Km_scalar(x, m):
    """K_m(x) for any m by quadrature (closed form for m<=1)."""
    if m <= 1:
        return float(Km(np.array([x]), m)[0])
    f = lambda a: np.cos(m * a) / np.sqrt(1.0 + x * x - 2.0 * x * np.cos(a))
    val, _ = quad(f, 0.0, np.pi, limit=200)
    return (x / np.pi) * val


def _Km_offsets(xs, m, eta=0.0):
    """K_m on an array of x values; vectorized alpha-trapezoid for m>1 or eta>0."""
    if eta == 0.0 and m <= 1:
        return Km(xs, m)
    a = np.linspace(0.0, np.pi, 16384)
    ca = np.cos(a)
    out = np.empty_like(xs)
    for i, x in enumerate(xs):
        f = np.cos(m * a) / np.sqrt(1.0 + x * x - 2.0 * x * ca + eta * eta)
        out[i] = (x / np.pi) * np.trapezoid(f, a)
    return out


def build_J_general(grid, model, m, eta_pert=0.0):
    """Poisson matrix with sigma0-ratio weighting (reduces to ars89 for power laws)."""
    N, log_f = grid.N, grid.log_f
    f = np.exp(log_f)
    j_minus_i = np.arange(N)[None, :] - np.arange(N)[:, None]
    X = f ** j_minus_i
    W = X * 0.5 * (f - 1.0 / f)
    W[:, 0] *= 0.5
    W[:, -1] *= 0.5
    sig = model.sigma0(grid.r)
    Srat = sig[None, :] / sig[:, None]                       # sigma0(r_j)/sigma0(r_i)

    # kernel on the 2N-1 unique offsets
    offs = f ** np.arange(-(N - 1), N, dtype=float)
    Koff = _Km_offsets(offs, m, eta_pert)
    Kmat = Koff[j_minus_i + (N - 1)]
    if eta_pert > 0.0:                                        # softened: smooth, plain trapezoid
        return -(Kmat * Srat * W)

    np.fill_diagonal(Kmat, 0.0)
    J = -(Kmat * Srat * W)

    # singular-band cubic weights, interpolated over the local slope p_loc(r_i)
    p_loc = -model.dln_sigma_dlnr(grid.r)
    pgrid = np.arange(np.floor(p_loc.min()) - 1, np.ceil(p_loc.max()) + 1.5, 0.5)
    Wc = {d: [] for d in (-2, -1, 0, 1, 2)}
    Wt1 = {-1: [], +1: []}
    for pg in pgrid:
        g = lambda s: _Km_scalar(np.exp(s), m) * np.exp((1.0 - pg) * s)
        w = cubic_singular_band_weights(g, log_f)
        for d in Wc:
            Wc[d].append(w[d])
        Wt1[-1].append(_Km_scalar(1.0 / f, m) * f ** pg * (1.0 - 1.0 / f) / 2.0)
        Wt1[+1].append(_Km_scalar(f, m) * f ** (-pg) * (f - 1.0) / 2.0)
    idx = np.arange(N)
    for d in (-2, -1, 0, 1, 2):
        val = np.interp(p_loc, pgrid, np.array(Wc[d]))
        corr = -val
        if d in (-1, +1):
            corr = corr + np.interp(p_loc, pgrid, np.array(Wt1[d]))
        i = idx[max(0, -d):N - max(0, d)]
        J[i, i + d] += corr[i]
    return J


def build_T_general(grid, model):
    """Indirect operator with sigma0-ratio weighting: T_ij = x^2 (sig_j/sig_i) w_ij."""
    N, log_f = grid.N, grid.log_f
    f = np.exp(log_f)
    j_minus_i = np.arange(N)[None, :] - np.arange(N)[:, None]
    X = f ** j_minus_i
    W = X * 0.5 * (f - 1.0 / f)
    W[:, 0] *= 0.5
    W[:, -1] *= 0.5
    sig = model.sigma0(grid.r)
    return X ** 2 * (sig[None, :] / sig[:, None]) * W


# ---------------- generalized eigenproblem ----------------

class GeneralEigenProblem(EigenProblem):
    """EigenProblem with the power-law V-bracket constants promoted to local fields."""

    @classmethod
    def build_family(cls, model, N=800, eta_eq=0.1, eta_pert=0.0):
        grid = make_log_grid(model.Rstar, model.Rd, N)
        ops = make_operators(grid)
        eq = build_equilibrium(model, grid, ops, eta=eta_eq)
        J = build_J_general(grid, model, model.m, eta_pert)
        T = build_T_general(grid, model)
        pr = cls(model=model, grid=grid, ops=ops, eq=eq, J=J, T=T)
        pr._LP1 = 1.0 + model.dln_sigma_dlnr(grid.r)
        pr._LP2 = pr._LP1 ** 2 - pr._LP1 + model.d_dlnsig_dlnr(grid.r)
        return pr

    def assemble(self, omega):
        model, r = self.model, self.grid.r
        q, mm = model.q, model.m
        D1, D2 = self.ops.D1, self.ops.D2
        eq = self.eq
        A, B, C = eq.coefficients(omega)
        Omega, kappa = eq.Omega, eq.kappa
        LP1, LP2 = self._LP1, self._LP2

        c1_D1 = A * r + 2.0 * LP1 - 1.0
        c1_D0 = B * r ** 2 + A * r * LP1 + LP2
        Bracket_V = D2 + c1_D1[:, None] * D1 + np.diag(c1_D0)
        term_V = Bracket_V @ (self.selfgrav_scale * self.J)

        Sig = 2.0 * np.pi * model.G * eq.sigma0 / eq.a0sq
        c2_D1 = A * r - 2.0 * q - 1.0
        c2_D0 = B * r ** 2 - A * r * q + q * (q + 1.0)
        Bracket_S = D2 + c2_D1[:, None] * D1 + np.diag(c2_D0)
        term_S = (1.0 / (Sig * r))[:, None] * Bracket_S

        nu = eq.nu(omega)
        term_C = np.diag(-kappa ** 2 * (1.0 - nu ** 2) * r / (2.0 * np.pi * model.G * eq.sigma0))

        c_ind = (B * r ** 2 + A * r) * 0.5 * omega ** 2 * r ** 3 / (model.G * model.Mtot)
        term_ind = (self.indirect_scale if mm == 1 else 0.0) * c_ind[:, None] * self.T

        M = term_V + term_S + term_C + term_ind

        Tmat = self._enthalpy_plus_potential_operator()
        N = self.grid.N
        i = 0
        M[i, :] = (D1 @ Tmat)[i] / r[i] - (2.0 * mm * Omega[i] / ((omega - mm * Omega[i]) * r[i])) * Tmat[i, :]
        j = N - 1
        Dscr = kappa[j] ** 2 - (omega - mm * Omega[j]) ** 2
        rhs = (eq.dsigma0[j] / Dscr) * ((D1 @ Tmat)[j] / r[j]
                                        - (2.0 * mm * Omega[j] / ((omega - mm * Omega[j]) * r[j])) * Tmat[j, :])
        row = np.zeros(N, dtype=complex)
        row[j] = eq.sigma0[j]
        M[j, :] = row - rhs
        return M


# ---------------- regression test vs ars89 ----------------

def regression_test():
    """With tapers off + pure power law, the generalized matrix must equal ars89's."""
    print("=== regression: generalized machinery vs ars89 (canonical power law) ===")
    dm = DiskModel(p=1.5, q=0.5, Rd_over_Rstar=1.0e4, Md_over_Mstar=1.0, Qstar=10.0, m=1)
    N = 400
    orig = EigenProblem.build(dm, N=N)

    class PL(FamilyModel):
        def __init__(s):
            s.p, s.q, s.m = dm.p, dm.q, 1
            s.Rin, s.Rout = 0.0, np.inf
            s.Sigma0 = dm.sigma_star * dm.Rstar ** dm.p       # sigma0 = Sigma0 r^-p
            s.cs0 = dm.a0_star * dm.Rstar ** (dm.q / 2.0)
            s.G, s.Mstar = dm.G, dm.Mstar
            s.Rstar, s.Rd = dm.Rstar, dm.Rd
            s.Md, s.Mtot = dm.Md, dm.Mtot

    pl = PL()
    gen = GeneralEigenProblem.build_family(pl, N=N, eta_eq=0.05)   # match ars89 build default
    for w in (4.26 - 0.232j, 2.0 - 0.5j):
        M0, M1 = orig.assemble(w), gen.assemble(w)
        rel = np.linalg.norm(M1 - M0) / np.linalg.norm(M0)
        print(f"  |M_gen - M_ars|/|M_ars| at omega={w}:  {rel:.2e}")
    w0, S0, _ = solve_mode(orig, 4.26 - 0.232j)
    w1, S1, _ = solve_mode(gen, 4.26 - 0.232j)
    print(f"  eigenvalue: ars89 {w0:.4f}   generalized {w1:.4f}   (ARS89 Fig3: 4.26-0.232j)")
    return abs(w1 - w0) < 1e-6


# ---------------- mode hunt for the LB94 family ----------------

def hunt(problem, re_range, im_range, n_re=25, n_im=14, verbose=True):
    """sigma_min scan + Newton polish from local minima; returns deduped growing roots."""
    res = np.linspace(*re_range, n_re)
    ims = np.linspace(*im_range, n_im)
    S = np.zeros((n_im, n_re))
    for a, wi in enumerate(ims):
        for b, wr in enumerate(res):
            S[a, b] = problem.smallest_singular_value(wr + 1j * wi)
    seeds = []
    for a in range(n_im):
        for b in range(n_re):
            lo_a, hi_a = max(a - 1, 0), min(a + 2, n_im)
            lo_b, hi_b = max(b - 1, 0), min(b + 2, n_re)
            if S[a, b] == S[lo_a:hi_a, lo_b:hi_b].min():
                seeds.append(res[b] + 1j * ims[a])
    roots = []
    for sd in seeds:
        w, Svec, info = solve_mode(problem, sd, max_iter=60)
        if not info["converged"] or info["residual"] > 1e-6:
            continue
        if w.imag > -1e-3:                                    # growing modes only
            continue
        if not (re_range[0] - 0.5 < w.real < re_range[1] + 0.5):
            continue
        if any(abs(w - w2) < 1e-3 for w2, _ in roots):
            continue
        roots.append((w, Svec))
    roots.sort(key=lambda t: t[0].imag)                       # fastest-growing first
    if verbose:
        print(f"  scan {n_re}x{n_im}, {len(seeds)} seeds -> {len(roots)} distinct growing roots")
    return roots, (res, ims, S)


def corotation(problem, w):
    Om, r = problem.eq.Omega, problem.grid.r
    Omp = w.real / problem.model.m
    inside = (r > 0.15) & (r < 3.0)
    return float(np.interp(-Omp, -Om[inside], r[inside]))     # Omega falls with r


def report(problem, roots, label):
    print(f"\n=== {label}: growing modes ===")
    print(f"  {'omega (code)':>22} {'gamma/Tu':>9} {'Omp/Tu':>8} {'R_cor(AU)':>10}")
    for w, _ in roots[:6]:
        g_tu = -w.imag * CONV
        omp_tu = (w.real / problem.model.m) * CONV
        rc = corotation(problem, w) * R_U_AU
        print(f"  {w.real:10.4f} {w.imag:+10.4f}i {g_tu:9.2f} {omp_tu:8.2f} {rc:10.0f}")


if __name__ == "__main__":
    print(f"units: t_u = {T_U_YR:.1f} yr; SPH T_unit = {SPH_TU_YR} yr -> conv x{CONV:.4f}")
    ok = regression_test()
    print(f"  regression {'PASS' if ok else 'FAIL'}\n")

    fam = FamilyModel(m=1)
    print(f"LB94 family: M_d = {fam.Md:.3f} M_* (target 1.10);  domain {fam.Rstar*100:.0f}-{fam.Rd*100:.0f} AU")
    pr1 = GeneralEigenProblem.build_family(fam, N=700, eta_eq=0.1)
    eq = pr1.eq
    Q = eq.a0sq ** 0.5 * eq.kappa / (np.pi * eq.sigma0)
    body = (pr1.grid.r > 0.3) & (pr1.grid.r < 2.2)
    print(f"  equilibrium check: Q_min = {Q[body].min():.2f} at "
          f"{pr1.grid.r[body][np.argmin(Q[body])]*100:.0f} AU (fit_basic_state: 1.35 at 56 AU)")

    # m=1 hunt: Omega_p up to Omega(20 AU); growth up to ~ 2x the SPH-measured rate
    roots1, scan1 = hunt(pr1, (0.1, 4.0), (-3.0, -0.05))
    report(pr1, roots1, "m=1 (indirect ON, razor-thin)")

    fam2 = FamilyModel(m=2)
    pr2 = GeneralEigenProblem.build_family(fam2, N=700, eta_eq=0.1)
    roots2, scan2 = hunt(pr2, (0.2, 8.0), (-3.0, -0.05))
    report(pr2, roots2, "m=2 (razor-thin)")

    np.savez("lb94/lineig_scan.npz",
             res1=scan1[0], ims1=scan1[1], S1=scan1[2],
             res2=scan2[0], ims2=scan2[1], S2=scan2[2],
             roots1=np.array([w for w, _ in roots1]),
             roots2=np.array([w for w, _ in roots2]),
             r=pr1.grid.r,
             modes1=np.array([S for _, S in roots1[:4]]) if roots1 else np.zeros((0,)),
             modes2=np.array([S for _, S in roots2[:4]]) if roots2 else np.zeros((0,)))
    print("\nwrote lb94/lineig_scan.npz")
