"""Visualize the 2D collapse: R-Z density field + midplane disk profile."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.loadtxt("lb94/collapse_rho.dat")
R = np.unique(d[:, 0]); Z = np.unique(d[:, 1])
NR, NZ = R.size, Z.size
rho = d[:, 2].reshape(NR, NZ)

mid = np.loadtxt("lb94/collapse_midplane.dat")   # R, rho, vphi, vkep

fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))

# density field, mirrored about the midplane for a full R-Z view
Zfull = np.concatenate([-Z[::-1], Z])
rhofull = np.concatenate([rho[:, ::-1], rho], axis=1)
im = ax[0].pcolormesh(R, Zfull, np.log10(rhofull.T + 1e-6), shading="auto", cmap="magma")
ax[0].set_xlabel("R"); ax[0].set_ylabel("Z"); ax[0].set_aspect("equal")
ax[0].set_title("log10 density — collapsed core + disk")
plt.colorbar(im, ax=ax[0], shrink=0.85)

ax[1].plot(mid[:, 0], mid[:, 1], "k-", lw=1.5, label="midplane density")
ax[1].set_xlabel("R"); ax[1].set_ylabel("rho (midplane)", color="k")
ax[1].set_yscale("log"); ax[1].set_xlim(0, 0.8)
a2 = ax[1].twinx()
m = mid[:, 1] > 1e-2
a2.plot(mid[m, 0], mid[m, 2], "b-", lw=1.3, label="v_phi (gas)")
a2.plot(mid[m, 0], mid[m, 3], "r--", lw=1.0, label="v_Kepler(core)")
a2.set_ylabel("v_phi", color="b"); a2.set_ylim(0, None)
ax[1].set_title("midplane: disk density + rotational support")
a2.legend(loc="upper right", fontsize=8)
plt.tight_layout(); plt.savefig("lb94/collapse.png", dpi=110)
print("disk peak rho =", mid[:, 1].max(), "at R =", mid[mid[:,1].argmax(), 0])
print("wrote lb94/collapse.png")
