"""Stuart-Landau saturation of the m=2 spiral (eigenmode/Galerkin reduction).

Computes the amplitude equation  dA/dt = s2 A - beta |A|^2 A  from the LKA98
three-wave equations, integrates it, and reports the saturation amplitude,
the m=2 density contrast (comparable to the hydro |C2|/|C0|), and the
pattern-speed shift.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lka98.model import DiskModel
from lka98.saturation import ThreeWave

tw = ThreeWave(DiskModel(), N=300, nu2=3e-4, ms=(0, 2, 4))
L = tw.landau()
s2, beta, gamma, Omp = L["s2"], L["beta"], L["gamma"], L["Omega_p"]
A2sat = gamma / beta.real
Asat = np.sqrt(A2sat)
dOmp = -beta.imag * A2sat / 2.0

print(f"dominant m=2 mode:  s2 = {s2.real:.3f}{s2.imag:+.3f}j   (gamma1={gamma:.3f}, Omega_p={Omp:.3f})")
print(f"Landau coefficient: beta = {beta.real:.2f}{beta.imag:+.2f}j   (Re>0 -> SUPERCRITICAL, saturates)")
print(f"saturation:         |A|_sat = {Asat:.4f}")
print(f"pattern speed:      Omega_p {Omp:.3f} -> {Omp+dOmp:.3f}  ({100*dOmp/Omp:+.1f}%)")

# physical m=2 density contrast at saturation: |C2(r)|/|C0(r)| = |A||phi2_sig|/Sigma0
phi2_sig = L["phi2"][:tw.N]
contrast = Asat * np.abs(phi2_sig) / tw.Sig0
print(f"m=2 density contrast |C2|/|C0|: peak {contrast.max():.3f} at r={tw.r[np.argmax(contrast)]:.3f}, "
      f"mean {np.mean(contrast):.3f}")

# integrate the Stuart-Landau ODE from a small seed
dt, T = 0.01, 25.0
t = np.arange(0, T, dt)
A = np.zeros(len(t), complex); A[0] = 1e-3
for n in range(len(t) - 1):
    k = s2 * A[n] - beta * np.abs(A[n])**2 * A[n]
    A[n+1] = A[n] + dt * k            # (phase rotates fast; |A| is what matters)
absA = np.abs(A)

fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].semilogy(t, absA, lw=1.5)
ax[0].axhline(Asat, color="r", ls="--", lw=0.8, label=f"|A|_sat={Asat:.3f}")
ax[0].semilogy(t, 1e-3*np.exp(gamma*t), "k:", lw=0.8, label=f"linear e^{{{gamma:.2f}t}}")
ax[0].set_xlabel("t"); ax[0].set_ylabel("|A| (m=2 amplitude)")
ax[0].set_ylim(1e-4, 1); ax[0].legend(); ax[0].set_title("Stuart-Landau saturation")
ax[1].plot(tw.r, contrast, lw=1.5)
ax[1].set_xlabel("r"); ax[1].set_ylabel("|C2|/|C0| at saturation")
ax[1].set_title("saturated m=2 density contrast")
plt.tight_layout(); plt.savefig("lka98/saturation.png", dpi=110)
print("wrote lka98/saturation.png")
