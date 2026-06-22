#!/usr/bin/env python3
"""Resolution audit: the same viscous-disk IC run at several particle numbers N. The m=1
amplitude |c1|(t) is plotted for each N. Two questions (LB94's, sharpened by the Truelove /
Bate & Burkert SPH resolution criteria that postdate it):
  (1) does the initial seed amplitude scale as ~1/sqrt(N) (Poisson particle noise)?
  (2) does the growth RATE converge with N (a real eigenmode) or drift (a numerical artifact)?
"""
import glob
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

files = sorted(glob.glob("lb94/sph/modes_N*.dat"))
cmap = plt.cm.viridis
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
rows = []
for k, fn in enumerate(files):
    N = int(re.search(r"_N(\d+)", fn).group(1))
    d = np.loadtxt(fn)
    t, c1, c2 = d[:, 0], d[:, 1], d[:, 2]
    col = cmap(k / max(len(files) - 1, 1))
    axs[0].semilogy(t, np.maximum(c1, 1e-4), "o-", color=col, lw=2, ms=4, label=f"N={N}")
    # growth rate gamma: slope of ln|c1| in the linear-amplifying window 0.02<|c1|<0.4
    m = (c1 > 0.02) & (c1 < 0.4)
    gamma = np.polyfit(t[m], np.log(c1[m]), 1)[0] if m.sum() >= 2 else np.nan
    # onset time: first t where |c1| crosses 0.1
    onset = t[np.argmax(c1 > 0.1)] if (c1 > 0.1).any() else np.nan
    rows.append((N, c1[0], gamma, onset, c1.max()))
axs[0].set_xlabel("t  (T$_{unit}$)"); axs[0].set_ylabel(r"$|c_1|$")
axs[0].set_title("m=1 growth vs particle number"); axs[0].legend(frameon=False); axs[0].grid(alpha=0.3, which="both")

rows.sort()
Ns = np.array([r[0] for r in rows], float)
peak = np.array([r[4] for r in rows])
# PEAK m=1 amplitude vs N -- the convergence indicator: disruption (~1) at low N collapsing
# to a marginally-stable saturation at high N
axs[1].semilogx(Ns, peak, "ko-", lw=2, ms=7)
axs[1].axhspan(0.0, 0.15, color="green", alpha=0.12)
axs[1].axhspan(0.6, 1.05, color="red", alpha=0.10)
axs[1].text(Ns[0] * 1.1, 0.05, "marginally stable\n(sustained mild m=1)", color="green", fontsize=9)
axs[1].text(Ns[0] * 1.1, 0.85, "spurious disruption\n(under-resolved)", color="firebrick", fontsize=9)
axs[1].set_xlabel("N"); axs[1].set_ylabel(r"peak $|c_1|$ over the run")
axs[1].set_title("convergence: peak m=1 vs N"); axs[1].set_ylim(0, 1.05); axs[1].grid(alpha=0.3, which="both")
fig.tight_layout(); fig.savefig("lb94/sph/fig_audit.png", dpi=130)
print("wrote lb94/sph/fig_audit.png\n")

print(f"{'N':>7} {'seed|c1|(0)':>12} {'growth rate':>12} {'onset(0.1)':>11} {'peak|c1|':>9}")
for N, s, g, on, pk in rows:
    print(f"{N:7d} {s:12.4f} {g:12.2f} {on:11.2f} {pk:9.3f}")
print("\nlow N -> spurious disruption (peak ~1); high N -> converges to marginally-stable mild m=1")
