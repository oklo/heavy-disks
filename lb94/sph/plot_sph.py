#!/usr/bin/env python3
"""LB94 3D SPH results: azimuthal mode growth |c_m|(t) (cf LB94 Fig. 6) and the face-on
surface-density / overdensity map showing the spiral (cf LB94 Fig. 5)."""
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

AU_per_Lu = 100.0

# ---- Fig: mode growth ----
m = np.loadtxt("lb94/sph/modes.dat")
t = m[:, 0]
fig, ax = plt.subplots(figsize=(6.4, 4.6))
for k, (lab, c) in enumerate([("m=1", "C3"), ("m=2", "C0"), ("m=3", "C2"), ("m=4", "C1")], start=1):
    ax.semilogy(t, np.maximum(m[:, k], 1e-4), "o-", color=c, lw=2, ms=4, label=lab)
ax.set_xlabel("time  (T$_{unit}$ = 474 yr)")
ax.set_ylabel(r"$|c_m|$  (azimuthal mode amplitude)")
ax.set_title("LB94 Fig. 6 — spiral mode growth (this work)")
ax.legend(frameon=False); ax.grid(alpha=0.3, which="both")
fig.tight_layout(); fig.savefig("lb94/sph/fig_modes.png", dpi=130)
print("wrote lb94/sph/fig_modes.png")

# ---- Fig: face-on overdensity (initial vs final snapshot) ----
snaps = sorted(glob.glob("lb94/sph/snap_*.dat"))
if snaps:
    def load(fn):
        d = np.loadtxt(fn)
        with open(fn) as f:
            tt = float(f.readline().split("t=")[1].split()[0])
        return d[:, 0] * AU_per_Lu, d[:, 1] * AU_per_Lu, tt
    fig, axs = plt.subplots(1, 2, figsize=(11, 5.2))
    for ax, fn in zip(axs, [snaps[0], snaps[-1]]):
        x, y, tt = load(fn)
        H, xe, ye = np.histogram2d(x, y, bins=180, range=[[-250, 250], [-250, 250]])
        ax.imshow(np.maximum(H.T, 0.3), origin="lower", extent=[-250, 250, -250, 250],
                  cmap="inferno", norm=LogNorm(vmin=0.5, vmax=max(H.max(), 2)))
        ax.set_title(f"t = {tt:.2f} T$_{{unit}}$"); ax.set_xlabel("x (AU)"); ax.set_ylabel("y (AU)")
    fig.suptitle("LB94 Fig. 5 — face-on surface density (this work): m=1 spiral", y=0.99)
    fig.tight_layout(); fig.savefig("lb94/sph/fig_spiral.png", dpi=130)
    print("wrote lb94/sph/fig_spiral.png")
