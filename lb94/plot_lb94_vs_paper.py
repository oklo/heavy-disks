#!/usr/bin/env python3
"""Overplot this work on LB94's Fig. 1 and Fig. 2.

LB94 curves are digitized by eye from the scanned page images (lb94/paper/pg-04.png) and
are therefore approximate (~few % in read-off). Fig. 1: core mass vs time. Fig. 2: disk
T, log Sigma, log j vs log10 R(cm) -- LB94 plot this at 20,000 yr, when their disk is
already formed; our disk forms ~30-40 kyr later (see Fig. 1), so for a structural
comparison we show this work at its disk-formed epoch (60 kyr) as well as at 20 kyr."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

G, Msun, AU = 6.67430e-8, 1.989e33, 1.496e13
C_INV, C_VIS, C_LB = "#1f77b4", "#d62728", "#444444"

# ---- LB94 Fig. 1 (digitized): core mass (Msun) vs time (yr) ----
lb1_inv_t = np.array([3,5,7,9,11,13,15,17,19,21,25,30,40,50,60,70]) * 1e3
lb1_inv_M = np.array([0,.04,.10,.17,.24,.31,.36,.39,.41,.415,.42,.42,.42,.42,.42,.42])
lb1_vis_t = np.array([3,5,7,9,11,13,15,17,19,21,25,30,35,40,50,60,70,75]) * 1e3
lb1_vis_M = np.array([0,.04,.11,.18,.25,.32,.37,.40,.43,.46,.49,.52,.535,.55,.575,.59,.60,.605])

# ---- LB94 Fig. 2 (digitized) at 20 kyr; x = log10 R(cm) ----
lb2_T_x = [14.2,14.35,14.5,14.65,14.8,14.95,15.1,15.25,15.4,15.55]
lb2_T_y = [205,175,150,128,105,88,70,55,43,33]
lb2_S_x = [14.2,14.4,14.55,14.7,14.85,15.0,15.15,15.3,15.45,15.55]
lb2_S_y = [1.65,1.8,1.9,1.9,1.85,1.75,1.55,1.2,0.4,-0.8]
lb2_j_x = [14.2,14.4,14.6,14.8,15.0,15.2,15.4,15.55]
lb2_j_y = [20.0,20.2,20.35,20.5,20.6,20.72,20.8,20.82]

# ================= Fig. 1 =================
inv = np.loadtxt("lb94/ybl_mcore_a0.000.dat")
vis = np.loadtxt("lb94/ybl_mcore_a0.010.dat")
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(lb1_inv_t/1e3, lb1_inv_M, "o", ms=4, color=C_LB, mfc="white", label="LB94 inviscid (digitized)")
ax.plot(lb1_vis_t/1e3, lb1_vis_M, "^", ms=4, color=C_LB, label="LB94 viscous (digitized)")
ax.plot(inv[:,0]/1e3, inv[:,1], "-",  color=C_INV, lw=2, label=r"this work, inviscid ($\alpha=0$)")
ax.plot(vis[:,0]/1e3, vis[:,1], "--", color=C_VIS, lw=2, label=r"this work, viscous ($\alpha=0.01$)")
ax.set_xlabel("T (yr)  /  10$^3$"); ax.set_ylabel(r"core mass  $M_c$  (M$_\odot$)")
ax.set_title("LB94 Fig. 1 vs this work — central object growth")
ax.set_xlim(0, 75); ax.set_ylim(0, 0.75)
ax.legend(loc="lower right", frameon=False, fontsize=9); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig("lb94/cmp_fig1.png", dpi=130)
print("wrote lb94/cmp_fig1.png")

# ================= Fig. 2 =================
def load(fn):
    d = np.loadtxt(fn)
    R_AU, logR, Sig, T, j = d[:,0], d[:,1], d[:,2], d[:,3], d[:,4]
    m = (Sig > 0.05) & (logR > 14.15) & (logR < 15.65) & (R_AU > 12)
    return logR[m], Sig[m], T[m], j[m]

# this work with irradiated cooling + Truelove-limited sink, at 40 kyr (both runs stable and
# disk-forming; LB94's Fig. 2 is at 20 kyr, when their disk is comparably developed).
xi, Si, Ti, ji = load("lb94/ybl_prof_a0.000_40kyr.dat")
xv, Sv, Tv, jv = load("lb94/ybl_prof_a0.010_40kyr.dat")

fig, axs = plt.subplots(3, 1, figsize=(7, 10), sharex=True)
# Temperature
axs[0].plot(lb2_T_x, lb2_T_y, "s-", ms=4, color=C_LB, label="LB94 (digitized, 20 kyr)")
axs[0].plot(xi, Ti, "-",  color=C_INV, lw=2, label="this work inviscid (40 kyr)")
axs[0].plot(xv, Tv, "--", color=C_VIS, lw=2, label="this work viscous (40 kyr)")
axs[0].set_ylabel("midplane T (K)"); axs[0].set_ylim(0, 260)
axs[0].legend(loc="upper right", frameon=False, fontsize=8)
axs[0].set_title("LB94 Fig. 2 vs this work (irradiated cooling + Truelove sink) — disk structure")
# log Sigma
axs[1].plot(lb2_S_x, lb2_S_y, "s-", ms=4, color=C_LB)
axs[1].plot(xi, np.log10(Si), "-",  color=C_INV, lw=2)
axs[1].plot(xv, np.log10(Sv), "--", color=C_VIS, lw=2)
axs[1].set_ylabel(r"log$_{10}\,\Sigma$  (g cm$^{-2}$)"); axs[1].set_ylim(-1.5, 4)
# log j
axs[2].plot(lb2_j_x, lb2_j_y, "s-", ms=4, color=C_LB)
axs[2].plot(xi, np.log10(ji), "-",  color=C_INV, lw=2)
axs[2].plot(xv, np.log10(jv), "--", color=C_VIS, lw=2)
axs[2].set_ylabel(r"log$_{10}\,j$  (cm$^2$ s$^{-1}$)"); axs[2].set_ylim(17, 21.3)
axs[2].set_xlabel(r"log$_{10}$ R (cm)"); axs[2].set_xlim(14.15, 15.65)
for a in axs: a.grid(alpha=0.3)
fig.tight_layout(); fig.savefig("lb94/cmp_fig2.png", dpi=130)
print("wrote lb94/cmp_fig2.png")
