#!/usr/bin/env python3
"""Figure for the LB94 linear eigenproblem: sigma_min landscapes (m=1, m=2) with the
located roots, and the leading eigenfunctions with corotation marked."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CONV, R_U_AU = 1.7369, 100.0
d = np.load("lb94/lineig_scan.npz")

fig, axs = plt.subplots(2, 2, figsize=(12, 8.5))
for k, m in enumerate((1, 2)):
    ax = axs[0, k]
    res, ims, S = d[f"res{m}"], d[f"ims{m}"], d[f"S{m}"]
    im = ax.pcolormesh(res, ims, np.log10(S), shading="auto", cmap="viridis")
    plt.colorbar(im, ax=ax, label=r"log$_{10}\,\sigma_{\min}$")
    roots = d[f"roots{m}"]
    if roots.size:
        ax.plot(roots.real, roots.imag, "r*", ms=12, mec="white")
    ax.set_xlabel(r"Re $\omega$  (code)"); ax.set_ylabel(r"Im $\omega$  (code)")
    ax.set_title(f"m={m}: $\\sigma_{{\\min}}$ landscape (stars = roots)")

r = d["r"] * R_U_AU
for k, m in enumerate((1, 2)):
    ax = axs[1, k]
    modes, roots = d[f"modes{m}"], d[f"roots{m}"]
    if modes.size:
        for j in range(min(2, modes.shape[0])):
            w = roots[j]
            lab = (f"$\\omega$={w.real:.2f}{w.imag:+.3f}i:  "
                   f"$\\gamma$={-w.imag*CONV:.2f}/T$_u$")
            ax.semilogx(r, np.abs(modes[j]) / np.abs(modes[j]).max(), lw=2, label=lab)
    ax.set_xlabel("R (AU)"); ax.set_ylabel(r"$|\sigma_1/\sigma_0|$ (normalized)")
    ax.set_xlim(8, 350); ax.legend(fontsize=8, frameon=False)
    ax.set_title(f"m={m} leading eigenfunctions")
    ax.grid(alpha=0.3, which="both")
fig.tight_layout(); fig.savefig("lb94/fig_lineig.png", dpi=130)
print("wrote lb94/fig_lineig.png")
