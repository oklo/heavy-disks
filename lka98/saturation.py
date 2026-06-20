"""Nonlinear saturation of the m=2 spiral: the full three-wave (m=0,2,4) system.

This is the integration LKA98 could not carry out in 1998 ("the full set of nine
equations is unfortunately quite susceptible to numerical instabilities when
treated as an initial value problem").  We integrate the lab-frame 3-harmonic
truncation of the 2D hydro,

    d z_m / dt = L_m z_m + N_m(z),     m = 0, 2, 4,   z_m = [sigma_m; u_m; v_m],

with an IMEX scheme: the stiff *linear* operator L_m (fast m*Omega rotation,
growth, and the nonlocal self-gravity Psi_m = P_m sigma_m) is treated
IMPLICITLY (cures the explicit instability), and the quadratic coupling N_m is
evaluated pseudo-spectrally (reconstruct on a phi-grid, form the nonlinear terms,
project back onto m=0,2,4) and stepped explicitly.

Linear operator (field ~ e^{-i m phi}, so d/dphi -> -i m):
    d sigma/dt = i m Omega sigma - (1/r) d_r(r Sigma0 u) + (i m Sigma0/r) v
    d u/dt     = i m Omega u + 2 Omega v - d_r[(cs^2/Sigma0) sigma + Psi]
    d v/dt     = i m Omega v - (kappa^2/2Omega) u + (i m/r)[(cs^2/Sigma0) sigma + Psi]
Its eigenvalues are s = i omega (omega the linear eigenfrequencies); the unstable
m=2 mode has Re(s) = gamma1 > 0, Im(s) = 2 Omega_p.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import lu_factor, lu_solve, eig

from .model import DiskModel
from .selfgravity import rotation_curve
from .poisson import build_poisson_matrix
from .eigensolve import LinearProblem, derivative_matrix, solve_mode


class ThreeWave:
    def __init__(self, model: DiskModel, N: int = 400, modulus: str = "correct",
                 nu4: float = 0.0, nu2: float = 3e-4, ms=(0, 2, 4)):
        self.model = model
        self.N = N
        self.ms = tuple(sorted(ms))                       # retained azimuthal harmonics
        self.nu4 = nu4                                    # 4th-order hyperdiffusion
        self.nu2 = nu2                                    # 2nd-order velocity viscosity
        # LOG radial grid (Appendix B): resolves the inner edge where the mode lives
        r = np.geomspace(model.Rin, model.RD, N)
        self.r = r
        self.D = derivative_matrix(r)
        # dissipative diffusion in u=ln r (uniform spacing -> clean [1,-2,1] Laplacian).
        # d^2/du^2 is negative semi-definite; this is the natural smoother on a log grid.
        du = np.log(model.RD / model.Rin) / (N - 1)
        Lap = np.zeros((N, N))
        for i in range(1, N - 1):
            Lap[i, i-1] = Lap[i, i+1] = 1.0 / du**2
            Lap[i, i] = -2.0 / du**2
        Lap[0, 0] = -2.0 / du**2; Lap[0, 1] = 2.0 / du**2     # reflecting (Neumann)
        Lap[-1, -1] = -2.0 / du**2; Lap[-1, -2] = 2.0 / du**2
        self.Lap = Lap
        self.D4 = Lap @ Lap                                   # hyperdiffusion (dissipative)
        Om, kap, _, _ = rotation_curve(model, r, modulus=modulus)
        self.Om, self.kap = Om, kap
        self.cs2 = model.cs2(r)
        self.Sig0 = model.Sigma0(r)
        self.W = cs2_over = self.cs2 / self.Sig0          # cs^2/Sigma0
        # Poisson matrices per harmonic (m>2 by quadrature)
        self.P = {m: build_poisson_matrix(model, r, m, quad=(m > 2)) for m in self.ms}
        self.L = {m: self._build_L(m) for m in self.ms}

    # ---- linear operator -------------------------------------------------
    def _build_L(self, m: int) -> np.ndarray:
        N, r, D = self.N, self.r, self.D
        Om, kap, cs2, Sig0 = self.Om, self.kap, self.cs2, self.Sig0
        dg = lambda v: np.diag(v).astype(complex)
        I = 1j
        Wd = dg(cs2 / Sig0)
        P = self.P[m]
        Z = np.zeros((N, N), complex)
        # sigma row
        Lss = dg(I * m * Om)
        Lsu = -dg(1.0 / r) @ D @ dg(r * Sig0)
        Lsv = dg(I * m * Sig0 / r)
        # u row
        Lus = -D @ (Wd + P)
        Luu = dg(I * m * Om)
        Luv = dg(2.0 * Om)
        # v row
        Lvs = dg(I * m / r) @ (Wd + P)
        Lvu = -dg(kap**2 / (2.0 * Om))
        Lvv = dg(I * m * Om)
        L = np.block([[Lss, Lsu, Lsv],
                      [Lus, Luu, Luv],
                      [Lvs, Lvu, Lvv]])
        # implicit dissipation (treated with the stiff linear part):
        #  - nu2 D^2 viscosity on the velocity equations damps the high-radial-
        #    wavenumber SECONDARY unstable modes (the cascade sink the truncation lacks),
        #    leaving the low-wavenumber dominant mode essentially untouched;
        #  - optional nu4 D^4 hyperdiffusion on all fields for grid-scale cleanup.
        H4 = -self.nu4 * self.D4
        H2 = self.nu2 * self.Lap
        for b in (0, 1, 2):
            L[b*N:(b+1)*N, b*N:(b+1)*N] += H4
        for b in (1, 2):                                  # u, v blocks only
            L[b*N:(b+1)*N, b*N:(b+1)*N] += H2
        # reflecting BC: u_m = 0 at both edges -> freeze those u-rows
        for i in (0, N - 1):
            L[N + i, :] = 0.0
        return L

    # ---- pseudo-spectral nonlinear coupling ------------------------------
    def nonlinear(self, z, Nphi: int = 64):
        """N_m for the retained harmonics from the quadratic terms of (9,10,11)."""
        N, r, D = self.N, self.r, self.D
        rr = r[:, None]
        def split(zm):
            return zm[:N], zm[N:2 * N], zm[2 * N:]
        s = {}; u = {}; v = {}
        for mi, m in enumerate(self.ms):
            s[m], u[m], v[m] = split(z[mi])
        phi = 2.0 * np.pi * np.arange(Nphi) / Nphi
        def recon(am):  # am: dict m-> amplitude(N,)
            f = np.zeros((N, Nphi), complex)
            dphi_f = np.zeros((N, Nphi), complex)
            dr_f = np.zeros((N, Nphi), complex)
            for m in self.ms:
                e = np.exp(-1j * m * phi)[None, :]
                f += am[m][:, None] * e
                dphi_f += (-1j * m) * am[m][:, None] * e
                dr_f += (D @ am[m])[:, None] * e
                if m != 0:                       # +c.c. (m<0 harmonics)
                    ec = np.exp(1j * m * phi)[None, :]
                    f += np.conj(am[m])[:, None] * ec
                    dphi_f += (1j * m) * np.conj(am[m])[:, None] * ec
                    dr_f += np.conj(D @ am[m])[:, None] * ec
            return f.real, dphi_f.real, dr_f.real
        S, dpS, drS = recon(s)
        U, dpU, drU = recon(u)
        V, dpV, drV = recon(v)
        NLs = -(1.0 / rr) * (D @ (rr * S * U)) - (1.0 / rr) * (S * dpV + V * dpS)
        NLu = V * V / rr - U * drU - (V / rr) * dpU
        NLv = -U * V / rr - U * drV - (V / rr) * dpV
        out = {}
        for name, NLf in (("s", NLs), ("u", NLu), ("v", NLv)):
            fh = np.fft.ifft(NLf, axis=1)
            out[name] = {m: fh[:, m] for m in self.ms}
        Nm = {}
        for mi, m in enumerate(self.ms):
            Nv = np.concatenate([out["s"][m], out["u"][m], out["v"][m]])
            Nv[N + 0] = 0.0; Nv[N + (N - 1)] = 0.0       # BC: no nl forcing on u edges
            Nm[mi] = Nv
        return Nm

    # ---- IMEX time stepping ----------------------------------------------
    def integrate(self, z0, dt, nsteps, Nphi=64, record_every=20):
        lu = {m: lu_factor(np.eye(3 * self.N) - dt * self.L[m]) for m in self.ms}
        z = [zc.astype(complex).copy() for zc in z0]
        hist = {"t": []}
        for m in self.ms:
            hist[f"a{m}"] = []
        for n in range(nsteps):
            Nm = self.nonlinear(z, Nphi)
            for mi, m in enumerate(self.ms):
                z[mi] = lu_solve(lu[m], z[mi] + dt * Nm[mi])
            if not all(np.isfinite(zc).all() for zc in z):
                print(f"  [blowup at t={n*dt:.3f}]")
                break
            if n % record_every == 0:
                hist["t"].append(n * dt)
                for mi, m in enumerate(self.ms):
                    sm = z[mi][:self.N]
                    hist[f"a{m}"].append(np.sqrt(np.trapezoid(np.abs(sm)**2, self.r)))
        for k in hist:
            hist[k] = np.array(hist[k])
        return z, hist

    # ---- Stuart-Landau reduction (eigenmode/Galerkin projection) ---------
    def landau(self):
        """Weakly-nonlinear amplitude equation  dA/dt = s2 A - beta |A|^2 A.

        Project onto the dominant m=2 eigenmode phi2 (eigenvalue s2 = gamma1 + i 2 Om_p);
        slave the second-order m=0 and m=4 responses (eqs 36-38, 42-44),
            zeta0 = (2 gamma1 I - L0)^-1 N0(phi2, phi2*),
            zeta4 = (2 s2     I - L4)^-1 N4(phi2, phi2),
        then form the cubic m=2 back-reaction (eqs 39-41) G = N2(zeta0,phi2)+N2(zeta4,phi2*)
        and project on the adjoint mode:  beta = - <psi2,G>/<psi2,phi2>.
        Saturation |A|^2 = gamma1/Re(beta); pattern-speed shift dOm_p = -Im(beta)|A|^2/2.
        """
        N = self.N
        i0, i2, i4 = (self.ms.index(m) for m in (0, 2, 4))
        # physical target eigenvalue s = i*omega from the validated linear solver
        lp = LinearProblem.build(self.model, N=N)
        omega, *_ = solve_mode(lp)
        target = 1j * omega
        w, vl, vr = eig(self.L[2], left=True, right=True)
        i = np.argmin(np.abs(w - target))                # select the physical dominant mode
        s2, phi2, psi2 = w[i], vr[:, i], vl[:, i]
        snorm = np.sqrt(np.trapezoid(np.abs(phi2[:N])**2, self.r))
        phi2 = phi2 / snorm                              # sigma-L2-norm(phi2) = 1
        gamma = s2.real
        # second-order responses (quadratic forcing from m=2 self-interaction only)
        z = [np.zeros(3 * N, complex) for _ in self.ms]
        z[i2] = phi2
        Nq = self.nonlinear(z)
        I3 = np.eye(3 * N)
        zeta0 = np.linalg.solve(2.0 * gamma * I3 - self.L[0], Nq[i0])
        zeta4 = np.linalg.solve(2.0 * s2 * I3 - self.L[4], Nq[i4])
        # cubic m=2 back-reaction: m=0 x m=2 and m=4 x m=2*
        zc = [np.zeros(3 * N, complex) for _ in self.ms]
        zc[i0], zc[i2], zc[i4] = zeta0, phi2, zeta4
        G = self.nonlinear(zc)[i2]
        beta = -np.vdot(psi2, G) / np.vdot(psi2, phi2)
        return dict(s2=s2, gamma=gamma, Omega_p=s2.imag / 2.0, beta=beta,
                    phi2=phi2, zeta0=zeta0, zeta4=zeta4,
                    A2sat=gamma / beta.real if beta.real > 0 else np.nan)

    # ---- initial condition from the linear m=2 eigenmode -----------------
    def eigenmode_ic(self, amp=1e-3):
        lp = LinearProblem.build(self.model, N=self.N)
        omega, s1, u1, _ = solve_mode(lp)
        m = self.model.m
        Wd = self.cs2 / self.Sig0
        W1 = Wd * s1 + self.P[2] @ s1
        A = 1j * (omega - m * self.Om)
        v1 = (-(self.kap**2 / (2.0 * self.Om)) * u1 + (1j * m / self.r) * W1) / A
        u1 = u1.copy(); u1[0] = u1[-1] = 0.0
        z2 = np.concatenate([s1, u1, v1])
        z2 = z2 / np.sqrt(np.trapezoid(np.abs(s1)**2, self.r)) * amp   # normalize sigma-norm to amp
        z0 = [np.zeros(3 * self.N, complex) for _ in self.ms]
        z0[self.ms.index(2)] = z2
        return z0, omega


if __name__ == "__main__":
    tw = ThreeWave(DiskModel(), N=300)
    # 1) validate L_2 spectrum: leading eigenvalue should be i*omega = gamma1 + i 2 Omega_p
    w = eig(tw.L[2], right=False)
    lead = w[np.argmax(w.real)]
    lp = LinearProblem.build(tw.model, N=300)
    omega, *_ = solve_mode(lp)
    print(f"L2 leading eigenvalue s = {lead:.4f}")
    print(f"   i*omega (from linear)  = {1j*omega:.4f}   (gamma1={-omega.imag:.3f}, 2Om_p={omega.real:.3f})")
