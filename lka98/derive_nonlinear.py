"""Independent symbolic derivation of the LKA98 weakly-nonlinear forcing terms.

Substitute the three-mode ansatz (eq 35) into the nonlinear governing equations
(9,10,11) for gamma_p = 2 (no pressure nonlinearity) and collect the m=0, m=2(=m1),
m=4(=2 m1) harmonics.  This independently reproduces the right-hand sides of the
paper's eqs (29-34) [2nd order], (36-38) [3rd-order m=0], (39-41) [m=2], (42-44) [m=4].

Convention: E = exp(i chi), chi = m (Omega_p t - phi), so d/dphi = -i m E d/dE.
Each field f = f0 + (1/2)(f1 E + f1c/E) + (1/2)(f2 E^2 + f2c/E^2), with f1c = conj(f1).
The amplitude equation for mode n has (1/2)(linear op on f_n) = [E^n coeff of NL] for
n = 1, 2 (so forcing = 2 * coeff), and (linear op on f_0) = [E^0 coeff] for n = 0.
"""

import sympy as sp

r, m = sp.symbols('r m', positive=True)
E = sp.Symbol('E')
I = sp.I

def F(name):
    return sp.Function(name)(r)

# amplitudes (m=0; m=2 and its conjugate; m=4 and its conjugate)
s0, u0, v0 = F('s0'), F('u0'), F('v0')
s1, u1, v1 = F('s1'), F('u1'), F('v1')
s1c, u1c, v1c = F('s1c'), F('u1c'), F('v1c')
s2, u2, v2 = F('s2'), F('u2'), F('v2')
s2c, u2c, v2c = F('s2c'), F('u2c'), F('v2c')

def field(f0, f1, f1c, f2, f2c):
    return f0 + sp.Rational(1, 2)*(f1*E + f1c/E) + sp.Rational(1, 2)*(f2*E**2 + f2c/E**2)

sig = field(s0, s1, s1c, s2, s2c)
uu  = field(u0, u1, u1c, u2, u2c)
vv  = field(v0, v1, v1c, v2, v2c)

def dphi(expr):
    return -I*m*E*sp.diff(expr, E)
def dr(expr):
    return sp.diff(expr, r)

# nonlinear RHS terms (gamma_p = 2; enthalpy linear, so no pressure nonlinearity)
NLu = vv**2/r - uu*dr(uu) - (vv/r)*dphi(uu)              # radial (eq 9)
NLv = -uu*vv/r - uu*dr(vv) - (vv/r)*dphi(vv)             # azimuthal (eq 10)
NLs = -(1/r)*dr(r*sig*uu) - (1/r)*dphi(sig*vv)           # continuity (eq 11)

def harmonic(expr, n):
    """Coefficient of E^n in the Laurent expansion (E^k auto-combine on expand)."""
    ex = sp.expand(expr)
    return sp.simplify(ex.coeff(E, n))

def forcing(expr, n):
    """Paper's RHS forcing for azimuthal mode n."""
    c = harmonic(expr, n)
    return sp.simplify(c if n == 0 else 2*c)

print("="*70)
print("m=0 forcing  (compare to NL1,NL2,NL3 = eqs 30,32,34 / 36-38)")
print("="*70)
print("continuity NL3 =", forcing(NLs, 0))
print("radial     NL1 =", forcing(NLu, 0))
print("azimuthal  NL2 =", forcing(NLv, 0))

print("\n" + "="*70)
print("m=4 (=2 m1) forcing  (compare to eqs 42,43,44)")
print("="*70)
print("continuity (42) =", forcing(NLs, 2))
print("radial     (43) =", forcing(NLu, 2))
print("azimuthal  (44) =", forcing(NLv, 2))

print("\n" + "="*70)
print("m=2 (=m1) forcing  (compare to eqs 39,40,41 RHS cubic couplings)")
print("="*70)
print("continuity (39) =", forcing(NLs, 1))
print("radial     (40) =", forcing(NLu, 1))
print("azimuthal  (41) =", forcing(NLv, 1))
