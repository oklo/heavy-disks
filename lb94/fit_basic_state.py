#!/usr/bin/env python3
"""Pin down the analytic basic-state family against the 40 kyr viscous disk:

    Sigma(R) = Sigma0 (R/R0)^(-p) exp[ -(R_in/R)^n - (R/R_out)^m ],   R0 = 100 AU
    c_s(R)   = c_s0 (R/R0)^(-q/2),  q = 1/2 (irradiated T ~ s^-1/2)

Fit (Sigma0, p, R_in, n, R_out, m) to the grid Sigma(R); fit c_s0 with q fixed; then check
the family REPRODUCES (not just fits): the rotation curve v_phi(R) computed self-consistently
from star + razor-thin-disk gravity (softened ring quadrature) + midplane pressure gradient,
and Q(R) from the family's own kappa. Neither v_phi nor Q is fitted -- they are predictions.
"""
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

kB, mH, mu, G, AU, Msun = 1.380649e-16, 1.6726e-24, 2.34, 6.674e-8, 1.496e13, 1.989e33
R0 = 100 * AU

# ---------------- load the 2D viscous disk ----------------
fn = "lb94/ybl_grid_a0.010_40kyr.dat"
with open(fn) as f:
    hdr = f.readline()
par = dict(t.split("=") for t in hdr.replace("#", "").split() if "=" in t)
NR, NZ, dZ, Mc = int(par["NR"]), int(par["NZ"]), float(par["dZ_cm"]), float(par["Mc_g"])
d = np.loadtxt(fn)
Rg = d[:, 0].reshape(NR, NZ)[:, 0]
Sig_d = d[:, 2].reshape(NR, NZ).sum(axis=1) * dZ * 2.0
vphi_d = d[:, 3].reshape(NR, NZ)[:, 0]
T_d = d[:, 4].reshape(NR, NZ)[:, 0]
cs_d = np.sqrt(kB * np.maximum(T_d, 1.0) / (mu * mH))

# fit window: disk body -- inside 18 AU the profile is the sink cliff (an artifact whose
# role the SPH heavy particle takes over), outside 280 AU the box edge
fitm = (Rg / AU > 18) & (Rg / AU < 280) & (Sig_d > 0.3)

# ---------------- Sigma family fit: integer taper indices, fit the 4 scales ----------------
def make_ln_sigma(n, m):
    def ln_sigma(theta, R):
        lS0, p, Rin, Rout = theta
        return lS0 - p * np.log(R / R0) - (Rin * AU / R) ** n - (R / (Rout * AU)) ** m
    return ln_sigma

fits = {}
print("taper-index scan (fit window 18-280 AU, %d points):" % fitm.sum())
print("   n  m   Sigma0     p    R_in  R_out   rms(lnSig)")
for n_ in (1, 2, 3, 4):
    for m_ in (1, 2, 3, 4):
        f = make_ln_sigma(n_, m_)
        s = least_squares(lambda th: f(th, Rg[fitm]) - np.log(Sig_d[fitm]),
                          [np.log(25.0), 2.5, 30.0, 220.0],
                          bounds=([0, 0.3, 5, 100], [8, 5, 90, 390]))
        r = np.sqrt(np.mean(s.fun ** 2))
        fits[(n_, m_)] = (r, s)
        print(f"   {n_}  {m_}  {np.exp(s.x[0]):7.1f} {s.x[1]:5.2f} {s.x[2]:6.1f} {s.x[3]:6.0f}"
              f"   {r:.3f}", "  <-- at bound" if (s.x[2] > 89 or s.x[3] > 385) else "")
# the rms floor (~0.12) is the disk's real ring/bump structure, and all 16 members sit on it:
# the data cannot select the taper indices. Canonical choice (n=1, m=2): within Delta-rms
# 0.003 of the formal best, its p matches the directly measured 60-180 AU slope (2.57),
# R_out matches where the disk actually ends, and the outer taper is Gaussian. The (1,1)
# member "wins" by bending p to 1.66 and faking the outer slope with a 106 AU exponential.
n, m = 1, 2
rms, sol = fits[(n, m)]
lS0, p, Rin, Rout = sol.x
Sig0 = np.exp(lS0)
print(f"\nchosen family: n = {n}, m = {m}  (tie-break: interpretability; see comment)")
print(f"  Sigma0 = {Sig0:.1f} g/cm^2 (at 100 AU),  p = {p:.2f},  R_in = {Rin:.1f} AU,  "
      f"R_out = {Rout:.0f} AU,  rms = {rms:.3f} ({100*(np.exp(rms)-1):.0f}%)")

ln_sigma_best = make_ln_sigma(n, m)
def Sigma_f(R):
    return np.exp(ln_sigma_best(sol.x, R))

# c_s normalization with q = 1/2 fixed
csm = (Rg / AU > 25) & (Rg / AU < 220) & (T_d > 10)
cs0 = np.exp(np.mean(np.log(cs_d[csm]) + 0.25 * np.log(Rg[csm] / R0)))
print(f"  c_s0   = {cs0/1e5:.3f} km/s at 100 AU (q = 1/2 fixed)  [T(100 AU) = {mu*mH*cs0**2/kB:.0f} K]")

def cs_f(R):
    return cs0 * (R / R0) ** -0.25

# disk masses
Rq = np.linspace(2 * AU, 400 * AU, 4000)
Md_fam = np.trapezoid(2 * np.pi * Rq * Sigma_f(Rq), Rq)
Md_grid = np.sum(2 * np.pi * Rg[fitm] * Sig_d[fitm]) * (Rg[1] - Rg[0])
print(f"  M_disk(family, 2-400 AU) = {Md_fam/Msun:.3f} Msun   [grid 15-280 AU: {Md_grid/Msun:.3f}]")
print(f"  M_star (grid core)       = {Mc/Msun:.3f} Msun   ->  M_d/M_* = {Md_fam/Mc:.2f}")

# ---------------- self-consistent rotation curve of the family ----------------
# razor-thin axisymmetric disk, softened ring quadrature:
#   g_R(R) = G int a da Sigma(a) int dphi (a cos(phi) - R) / (R^2+a^2-2aR cos(phi)+eps^2)^{3/2}
def disk_gR(Rev, eps_cm):
    a = np.linspace(2 * AU, 400 * AU, 1200)
    da = a[1] - a[0]
    ph = np.linspace(0, np.pi, 720)                       # symmetric: integrate half, double
    dph = ph[1] - ph[0]
    Sa = Sigma_f(a)
    g = np.zeros_like(Rev)
    for i, R in enumerate(Rev):
        den = (R * R + a[:, None] ** 2 - 2 * a[:, None] * R * np.cos(ph[None, :]) + eps_cm ** 2) ** 1.5
        num = a[:, None] * np.cos(ph[None, :]) - R
        g[i] = 2 * G * np.sum(Sa[:, None] * a[:, None] * num / den) * da * dph
    return g

Rev = np.linspace(15 * AU, 300 * AU, 120)
for eps_au, tag in [(1.0, "thin (eps=1 AU)"), (None, "thick (eps=H/2)")]:
    eps = eps_au * AU if eps_au else None
    if eps is None:
        # thickness-mimicking softening: eps = H/2 = cs/(2 Om_kep) -- evaluate iteratively once
        Omk = np.sqrt(G * Mc / Rev ** 3)
        eps_arr = 0.5 * cs_f(Rev) / Omk
        g_disk = np.array([disk_gR(np.array([R]), e)[0] for R, e in zip(Rev, eps_arr)])
    else:
        g_disk = disk_gR(Rev, eps)
    v2 = G * Mc / Rev + np.maximum(-g_disk * Rev, 0) * 0 + (-g_disk) * Rev   # star + disk
    # midplane pressure correction: P_mid ~ Sigma Omega cs (rho_mid ~ Sigma/(sqrt(2pi) H))
    Om_tmp = np.sqrt(np.maximum(v2, 1e-10)) / Rev
    lnP = np.log(Sigma_f(Rev) * Om_tmp * cs_f(Rev))
    dlnP = np.gradient(lnP, np.log(Rev))
    v2p = v2 + cs_f(Rev) ** 2 * dlnP
    vphi_fam = np.sqrt(np.maximum(v2p, 0))
    if eps_au:
        vfam_thin = vphi_fam.copy()
    else:
        vfam_thick = vphi_fam.copy()

# family Q from its own rotation curve (thick version, closer to the 3D disk)
Om_f = vfam_thick / Rev
kap2_f = np.gradient(Rev ** 4 * Om_f ** 2, Rev) / Rev ** 3
Q_f = cs_f(Rev) * np.sqrt(np.maximum(kap2_f, 0)) / (np.pi * G * Sigma_f(Rev))

# grid Q for comparison (real kappa)
Om_d = np.where(Rg > 0, vphi_d / Rg, 0)
kap2_d = np.where(Rg > 0, np.gradient(Rg ** 4 * Om_d ** 2, Rg) / Rg ** 3, 0)
Q_d = np.where(Sig_d > 0, cs_d * np.sqrt(np.maximum(kap2_d, 0)) / (np.pi * G * Sig_d), np.inf)

# residual summaries over the disk body (30-220 AU)
body = (Rev / AU > 30) & (Rev / AU < 220)
vd_i = np.interp(Rev, Rg, vphi_d)
Qd_i = np.interp(Rev, Rg, Q_d)
print("\nPredictions vs grid (30-220 AU):")
print(f"  v_phi: rms (thin)  = {100*np.sqrt(np.mean((vfam_thin[body]/vd_i[body]-1)**2)):.1f}%")
print(f"  v_phi: rms (thick) = {100*np.sqrt(np.mean((vfam_thick[body]/vd_i[body]-1)**2)):.1f}%")
print(f"  Q:     rms (thick) = {100*np.sqrt(np.mean((Q_f[body]/Qd_i[body]-1)**2)):.1f}%")
print(f"  Q_min(family) = {Q_f[body].min():.2f} at {Rev[body][np.argmin(Q_f[body])]/AU:.0f} AU"
      f"   [grid: {Qd_i[body].min():.2f} at {Rev[body][np.argmin(Qd_i[body])]/AU:.0f} AU]")

# ---------------- figure ----------------
fig, axs = plt.subplots(2, 2, figsize=(11.5, 8.5))
ax = axs[0, 0]
ax.loglog(Rg[fitm] / AU, Sig_d[fitm], "k.", ms=4, label="2D viscous disk (40 kyr)")
ax.loglog(Rev / AU, Sigma_f(Rev), "r-", lw=2,
          label=f"family: $\\Sigma_0$={Sig0:.0f}, p={p:.2f},\n$R_{{in}}$={Rin:.0f} AU (n={n:.1f}), "
                f"$R_{{out}}$={Rout:.0f} AU (m={m:.1f})")
ax.set_xlabel("R (AU)"); ax.set_ylabel(r"$\Sigma$ (g cm$^{-2}$)")
ax.set_title("surface density: fit"); ax.legend(fontsize=8, frameon=False); ax.grid(alpha=0.3, which="both")

ax = axs[0, 1]
ax.plot(Rg[fitm] / AU, vphi_d[fitm] / 1e5, "k.", ms=4, label="grid $v_\\phi$ (midplane)")
ax.plot(Rev / AU, vfam_thin / 1e5, "b--", lw=1.5, label="family, thin disk")
ax.plot(Rev / AU, vfam_thick / 1e5, "r-", lw=2, label="family, $\\epsilon$=H/2 (thick)")
ax.plot(Rev / AU, np.sqrt(G * Mc / Rev) / 1e5, "g:", lw=1.5, label="star alone (Kepler)")
ax.set_xlabel("R (AU)"); ax.set_ylabel(r"$v_\phi$ (km/s)")
ax.set_title("rotation curve: PREDICTED (not fitted)"); ax.legend(fontsize=8, frameon=False); ax.grid(alpha=0.3)

ax = axs[1, 0]
ax.plot(Rg[fitm] / AU, np.minimum(Q_d[fitm], 12), "k.", ms=4, label="grid Q")
ax.plot(Rev / AU, Q_f, "r-", lw=2, label="family Q")
ax.axhline(1, color="gray", ls=":", lw=1)
ax.set_xlabel("R (AU)"); ax.set_ylabel("Q"); ax.set_ylim(0, 6)
ax.set_title("Toomre Q: PREDICTED (not fitted)"); ax.legend(fontsize=8, frameon=False); ax.grid(alpha=0.3)

ax = axs[1, 1]
ax.loglog(Rg[csm] / AU, cs_d[csm] / 1e5, "k.", ms=4, label="grid $c_s$")
ax.loglog(Rev / AU, cs_f(Rev) / 1e5, "r-", lw=2, label=f"$c_s$ = {cs0/1e5:.2f} (R/100AU)$^{{-1/4}}$ km/s")
ax.set_xlabel("R (AU)"); ax.set_ylabel(r"$c_s$ (km/s)")
ax.set_title("sound speed: q = 1/2"); ax.legend(fontsize=8, frameon=False); ax.grid(alpha=0.3, which="both")
fig.tight_layout(); fig.savefig("lb94/fig_basic_state.png", dpi=130)
print("\nwrote lb94/fig_basic_state.png")
