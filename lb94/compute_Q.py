#!/usr/bin/env python3
"""Toomre Q(R) of the 2D collapse disk, computed properly: Q = c_s kappa / (pi G Sigma),
with the epicyclic frequency kappa from the actual rotation curve (kappa^2 = R^-3 d(R^4 Om^2)/dR)
rather than the kappa=Omega Keplerian approximation. Compares the inviscid and viscous disks
to LB94's stated minimum Q ~ 1.3 (at ~80 AU)."""
import sys
import numpy as np

kB, mH, mu, G, AU = 1.380649e-16, 1.6726e-24, 2.34, 6.674e-8, 1.496e13

def load(fn):
    with open(fn) as f:
        hdr = f.readline()
    p = dict(tok.split("=") for tok in hdr.replace("#", "").split() if "=" in tok)
    NR, NZ, dZ = int(p["NR"]), int(p["NZ"]), float(p["dZ_cm"])
    d = np.loadtxt(fn)
    R = d[:, 0].reshape(NR, NZ)
    rho = d[:, 2].reshape(NR, NZ)
    vphi = d[:, 3].reshape(NR, NZ)
    T = d[:, 4].reshape(NR, NZ)
    Rmid = R[:, 0]
    Sigma = rho.sum(axis=1) * dZ * 2.0                 # int rho dZ over both sides
    cs = np.sqrt(kB * np.maximum(T[:, 0], 1.0) / (mu * mH))   # isothermal sound speed (midplane)
    Om = np.where(Rmid > 0, vphi[:, 0] / Rmid, 0.0)
    # kappa^2 = (1/R^3) d(R^4 Om^2)/dR  via central differences
    f = Rmid**4 * Om**2
    dfdR = np.gradient(f, Rmid)
    kap2 = np.where(Rmid > 0, dfdR / Rmid**3, 0.0)
    kap = np.sqrt(np.maximum(kap2, 0.0))
    return Rmid, Sigma, cs, Om, kap, float(p["alpha"])

def report(fn):
    Rmid, Sigma, cs, Om, kap, alpha = load(fn)
    Q_kap = np.where(Sigma > 0, cs * kap / (np.pi * G * Sigma), np.inf)
    Q_om  = np.where(Sigma > 0, cs * Om  / (np.pi * G * Sigma), np.inf)   # kappa=Omega approx
    disk = (Rmid / AU > 20) & (Rmid / AU < 250) & (Sigma > 0.5)
    print(f"\n{fn}  (alpha={alpha})")
    print("  R(AU)   Sigma     cs(km/s)  kap/Om   Q(real kappa)  Q(kappa=Om)")
    for i in range(len(Rmid)):
        if not disk[i]:
            continue
        if i % 3 == 0:
            print(f"  {Rmid[i]/AU:5.0f}  {Sigma[i]:7.2f}   {cs[i]/1e5:6.2f}   {kap[i]/max(Om[i],1e-30):5.2f}"
                  f"     {Q_kap[i]:6.2f}        {Q_om[i]:6.2f}")
    qd = Q_kap[disk]
    print(f"  --> min Q(real kappa) = {qd.min():.2f} at R = {Rmid[disk][qd.argmin()]/AU:.0f} AU"
          f"   (kappa=Om gives min {Q_om[disk].min():.2f})   [LB94: min Q ~ 1.3 at 80 AU]")

if __name__ == "__main__":
    for fn in sys.argv[1:]:
        report(fn)
