#!/usr/bin/env python3
"""Is the SPH m=1 a genuine linear instability or a start-up/noise transient?

Tests, using the existing resolution-audit data:
  A) growth-rate fits per N in a common amplitude window + the noise floor vs the Poisson
     expectation E|c1| = sqrt(pi/4N). If the mode is a real eigenmode seeded by particle
     noise, gamma is N-independent and the onset time shifts by ln(sqrt(N/N0))/gamma.
  B) pattern speed from the 25k snapshot sequence: a global eigenmode rotates rigidly
     (one Omega_p across radii, corotation inside the disk); a material/transient feature
     winds at the local Omega(R).
  C) start-up transient check: Sigma(R) at t = 0, 0.1, 0.2, 0.3 (m=0 rearrangement).
  D) swing amplification: Toomre X_m(R) = kappa^2 R/(2 pi G Sigma m) for m=1,2 from the
     viscous 2D grid (swing is strong only for X <~ 3).
"""
import glob
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

kB, mH, mu, G, AU = 1.380649e-16, 1.6726e-24, 2.34, 6.674e-8, 1.496e13

# ---------------- A) growth vs N ----------------
files = sorted(glob.glob("lb94/sph/modes_N*.dat"))
data = {}
for fn in files:
    N = int(re.search(r"_N(\d+)", fn).group(1))
    data[N] = np.loadtxt(fn)

def onset(t, c, th):
    for i in range(1, len(c)):
        if c[i] >= th > c[i - 1]:
            f = (np.log(th) - np.log(c[i - 1])) / (np.log(c[i]) - np.log(c[i - 1]))
            return t[i - 1] + f * (t[i] - t[i - 1])
    return np.nan

print("=== A) m=1 growth vs N ===")
print(f"{'N':>7} {'noise(t<0.3)':>13} {'Poisson sqrt(pi/4N)':>20} {'gamma':>7} {'t(|c1|=0.2)':>12}")
gam = {}
for N, d in sorted(data.items()):
    t, c1 = d[:, 0], d[:, 1]
    noise = c1[t <= 0.31].mean()
    m = (c1 > 0.03) & (c1 < 0.55)
    g = np.polyfit(t[m], np.log(c1[m]), 1)[0] if m.sum() >= 2 else np.nan
    gam[N] = g
    print(f"{N:7d} {noise:13.4f} {np.sqrt(np.pi/(4*N)):20.4f} {g:7.2f} {onset(t, c1, 0.2):12.2f}")

Ns = sorted(data)
gref = np.nanmean([gam[N] for N in Ns])
print(f"\nmean gamma = {gref:.2f} /T_unit;  onset-delay test vs seed ~ 1/sqrt(N):")
t0 = {N: onset(data[N][:, 0], data[N][:, 1], 0.2) for N in Ns}
for N in Ns[1:]:
    pred = np.log(np.sqrt(N / Ns[0])) / gref
    print(f"  N {Ns[0]:>5} -> {N:>5}:  predicted delay +{pred:.2f}   observed +{t0[N]-t0[Ns[0]]:.2f}")

# ---------------- B) pattern speed from snapshots ----------------
print("\n=== B) pattern speed (25k snapshot sequence) ===")
snaps = sorted(glob.glob("lb94/sph/snap_*.dat"))
annuli = [(0.3, 0.7), (0.7, 1.2), (1.2, 1.8)]
ts, phs, amps = [], [], []
sigma_prof = {}
for fn in snaps:
    with open(fn) as f:
        tt = float(f.readline().split("t=")[1].split()[0])
    d = np.loadtxt(fn)
    x, y = d[:, 0], d[:, 1]
    R, phi = np.hypot(x, y), np.arctan2(y, x)
    row_p, row_a = [], []
    for (r0, r1) in annuli:
        msk = (R > r0) & (R < r1)
        c = np.exp(1j * phi[msk]).sum() / max(msk.sum(), 1)
        row_p.append(np.angle(c)); row_a.append(np.abs(c))
    ts.append(tt); phs.append(row_p); amps.append(row_a)
    if tt < 0.35:                                   # early Sigma(R) for the transient check
        h, edges = np.histogram(R, bins=np.linspace(0.1, 2.3, 34))
        rc = 0.5 * (edges[1:] + edges[:-1])
        sigma_prof[round(tt, 1)] = h / (2 * np.pi * rc * np.diff(edges))
ts = np.array(ts); phs = np.unwrap(np.array(phs), axis=0); amps = np.array(amps)

# Omega(R) from the IC velocities (rotation changes little over the run)
ic = np.loadtxt("lb94/sph/disk_ic.dat")
hv = ic[:, 8] == 1
xi, yi, vxi, vyi = ic[~hv, 0], ic[~hv, 1], ic[~hv, 3], ic[~hv, 4]
Ri = np.hypot(xi, yi)
Omi = (xi * vyi - yi * vxi) / np.maximum(Ri, 1e-9) ** 2
rb = np.linspace(0.2, 2.2, 21)
Om_prof = np.array([np.median(Omi[(Ri > rb[k]) & (Ri < rb[k + 1])]) for k in range(20)])
rbc = 0.5 * (rb[1:] + rb[:-1])

grow = ts >= 0.4
print(f"{'annulus(AU)':>14} {'Omega_p(rad/Tu)':>16} {'local Omega range':>18} {'corotation(AU)':>15}")
for k, (r0, r1) in enumerate(annuli):
    Omp = np.polyfit(ts[grow], phs[grow, k], 1)[0]
    Om_in = np.interp(r1, rbc, Om_prof); Om_out = np.interp(r0, rbc, Om_prof)
    Rcor = np.interp(-Omp if Omp < 0 else Omp, Om_prof[::-1], rbc[::-1])  # Om falls with R
    print(f"  {r0*100:4.0f}-{r1*100:4.0f}     {Omp:16.2f} {Om_out:9.2f}-{Om_in:.2f} {Rcor*100:15.0f}")

# ---------------- C) start-up Sigma(R) ----------------
print("\n=== C) early Sigma(R) evolution (start-up transient check) ===")
keys = sorted(sigma_prof)
if len(keys) >= 2:
    s0 = sigma_prof[keys[0]]
    for kk in keys[1:]:
        rel = np.abs(sigma_prof[kk] / np.maximum(s0, 1) - 1)
        msk = s0 > 20
        print(f"  t={kk:.1f}: max |dSigma/Sigma| (resolved bins) = {rel[msk].max():.3f}, "
              f"median = {np.median(rel[msk]):.3f}")

# ---------------- D) swing X(R) ----------------
print("\n=== D) Toomre swing parameter X_m = kappa^2 R / (2 pi G Sigma m), viscous IC ===")
fn = "lb94/ybl_grid_a0.010_40kyr.dat"
with open(fn) as f:
    hdr = f.readline()
pp = dict(tok.split("=") for tok in hdr.replace("#", "").split() if "=" in tok)
NR, NZ, dZ = int(pp["NR"]), int(pp["NZ"]), float(pp["dZ_cm"])
d = np.loadtxt(fn)
Rg = d[:, 0].reshape(NR, NZ)[:, 0]
rho = d[:, 2].reshape(NR, NZ)
vphi = d[:, 3].reshape(NR, NZ)[:, 0]
Tg = d[:, 4].reshape(NR, NZ)[:, 0]
Sig = d[:, 2].reshape(NR, NZ).sum(axis=1) * dZ * 2
Om = np.where(Rg > 0, vphi / Rg, 0)
kap2 = np.where(Rg > 0, np.gradient(Rg ** 4 * Om ** 2, Rg) / Rg ** 3, 0)
cs = np.sqrt(kB * np.maximum(Tg, 1) / (mu * mH))
Q = np.where(Sig > 0, cs * np.sqrt(np.maximum(kap2, 0)) / (np.pi * G * Sig), np.inf)
X1 = np.where(Sig > 0, kap2 * Rg / (2 * np.pi * G * Sig), np.inf)
disk = (Rg / AU > 25) & (Rg / AU < 220) & (Sig > 1)
print(f"{'R(AU)':>7} {'Q':>6} {'X(m=1)':>7} {'X(m=2)':>7}")
for i in np.where(disk)[0][::4]:
    print(f"{Rg[i]/AU:7.0f} {Q[i]:6.2f} {X1[i]:7.1f} {X1[i]/2:7.1f}")
print("(swing amplification strong only for X <~ 3 at Q ~ 1-2; Toomre 1981)")

# ---------------- figure ----------------
fig, axs = plt.subplots(2, 2, figsize=(12, 9))
cmap = plt.cm.viridis
# (a) raw + time-shifted growth curves
ax = axs[0, 0]
for k, N in enumerate(Ns):
    t, c1 = data[N][:, 0], data[N][:, 1]
    col = cmap(k / (len(Ns) - 1))
    sh = np.log(np.sqrt(N / Ns[0])) / gref
    ax.semilogy(t, np.maximum(c1, 1e-4), "o-", color=col, ms=3, lw=1, alpha=0.35)
    ax.semilogy(t - sh, np.maximum(c1, 1e-4), "s--", color=col, ms=4, lw=2,
                label=f"N={N} (shifted {-sh:+.2f})")
ax.set_xlabel("t (raw: faint;  shifted by $-\\ln\\sqrt{N/N_0}/\\gamma$: bold)")
ax.set_ylabel("$|c_1|$"); ax.legend(fontsize=8, frameon=False)
ax.set_title("common-$\\gamma$, noise-seed collapse test"); ax.grid(alpha=0.3, which="both")
# (b) phase vs t per annulus
ax = axs[0, 1]
for k, (r0, r1) in enumerate(annuli):
    ax.plot(ts, phs[:, k], "o-", label=f"{r0*100:.0f}-{r1*100:.0f} AU")
ax.axvspan(0.4, 1.0, color="gold", alpha=0.15)
ax.set_xlabel("t (T$_{unit}$)"); ax.set_ylabel("m=1 phase (rad, unwrapped)")
ax.set_title("pattern coherence across annuli (25k run)"); ax.legend(frameon=False); ax.grid(alpha=0.3)
# (c) early Sigma(R)
ax = axs[1, 0]
for kk in keys:
    ax.semilogy(rc * 100, np.maximum(sigma_prof[kk], 0.5), label=f"t={kk:.1f}")
ax.set_xlabel("R (AU)"); ax.set_ylabel("particle surface density (arb)")
ax.set_title("start-up m=0 rearrangement"); ax.legend(frameon=False); ax.grid(alpha=0.3)
# (d) Q and X profiles
ax = axs[1, 1]
ax.plot(Rg[disk] / AU, Q[disk], "k-", lw=2, label="Q")
ax.plot(Rg[disk] / AU, X1[disk], "r-", lw=2, label="X (m=1)")
ax.plot(Rg[disk] / AU, X1[disk] / 2, "b--", lw=2, label="X (m=2)")
ax.axhline(3, color="gray", ls=":", lw=1); ax.text(180, 3.15, "swing cutoff", color="gray", fontsize=8)
ax.set_xlabel("R (AU)"); ax.set_ylabel("Q,  X$_m$"); ax.set_ylim(0, 12)
ax.set_title("swing-amplification parameters (viscous IC)"); ax.legend(frameon=False); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig("lb94/sph/fig_m1_analysis.png", dpi=130)
print("\nwrote lb94/sph/fig_m1_analysis.png")
