#!/usr/bin/env python3
"""Reproduce LB94 Fig. 1 (central-object mass vs time) and Fig. 2 (disk T, Sigma, j
profiles), comparing the viscous (alpha=0.01) and inviscid (alpha=0) collapses."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

G, Msun, AU = 6.67430e-8, 1.989e33, 1.496e13
C_INV, C_VIS = "#1f77b4", "#d62728"

# ---------- Fig. 1: core mass vs time ----------
inv = np.loadtxt("lb94/ybl_mcore_a0.000.dat")
vis = np.loadtxt("lb94/ybl_mcore_a0.010.dat")
fig, ax = plt.subplots(figsize=(6.2, 4.6))
ax.plot(inv[:, 0] / 1e3, inv[:, 1], "-", color=C_INV, lw=2, label=r"inviscid ($\alpha=0$)")
ax.plot(vis[:, 0] / 1e3, vis[:, 1], "--", color=C_VIS, lw=2, label=r"viscous ($\alpha=0.01$)")
ax.axhline(0.68, color="gray", ls=":", lw=1)
ax.text(1, 0.695, "LB94: 0.68 M$_\\odot$", color="gray", fontsize=9, ha="left", va="bottom")
ax.set_xlabel("time  (10$^3$ yr)")
ax.set_ylabel(r"central object mass  $M_c$  (M$_\odot$)")
ax.set_ylim(-0.02, 0.74)
ax.set_title("LB94 Fig. 1 — central object growth")
ax.legend(loc="center right", frameon=False)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("lb94/fig1_mcore.png", dpi=130)
print("wrote lb94/fig1_mcore.png")

# ---------- Fig. 2: disk profiles ----------
def load(fn):
    d = np.loadtxt(fn)
    R, Sig, T, j = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
    m = (Sig > 0.05) & (R < 300) & (R > 30)        # resolved disk, outside the sink
    return R[m], Sig[m], T[m], j[m]

Ri, Si, Ti, ji = load("lb94/ybl_prof_a0.000_40kyr.dat")
Rv, Sv, Tv, jv = load("lb94/ybl_prof_a0.010_40kyr.dat")

fig, axs = plt.subplots(3, 1, figsize=(6.2, 9.2), sharex=True)
axs[0].plot(Ri, Ti, "-", color=C_INV, lw=2, label=r"inviscid ($\alpha=0$)")
axs[0].plot(Rv, Tv, "--", color=C_VIS, lw=2, label=r"viscous ($\alpha=0.01$)")
axs[0].set_ylabel("midplane T  (K)")
axs[0].legend(loc="upper right", frameon=False)
axs[0].set_title("LB94 Fig. 2 — disk structure at 60 kyr")

axs[1].semilogy(Ri, Si, "-", color=C_INV, lw=2)
axs[1].semilogy(Rv, Sv, "--", color=C_VIS, lw=2)
axs[1].set_ylabel(r"$\Sigma$  (g cm$^{-2}$)")

axs[2].plot(Ri, ji / 1e20, "-", color=C_INV, lw=2)
axs[2].plot(Rv, jv / 1e20, "--", color=C_VIS, lw=2)
# Keplerian reference about the final core mass
Mc = inv[-1, 1] * Msun
Rk = np.linspace(Ri.min(), Ri.max(), 100)
jk = np.sqrt(G * Mc * Rk * AU) / 1e20
axs[2].plot(Rk, jk, ":", color="gray", lw=1.2, label=r"Keplerian ($M_c$)")
axs[2].set_ylabel(r"specific $j$  (10$^{20}$ cm$^2$ s$^{-1}$)")
axs[2].set_xlabel("R  (AU)")
axs[2].legend(loc="upper left", frameon=False)

for a in axs:
    a.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("lb94/fig2_profiles.png", dpi=130)
print("wrote lb94/fig2_profiles.png")
