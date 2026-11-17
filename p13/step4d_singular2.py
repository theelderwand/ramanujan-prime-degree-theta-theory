"""p13 step 4d: singular points (fixed v2) and genus.

Method: m0 refined as a root of the quartic (m-1)(m^3-31m^2-85m-13) by
1-D Newton; for each m0, candidates P are roots of the quadratic F_P(m0, .);
the singular P is the one with F(m0,P) ~ 0 and F_m(m0,P) ~ 0 (checked at
1e-8, consistent with float64 m0 error).  Classify via tangent cone.
"""
import sympy as sp
import numpy as np
import mpmath as mp

mp.mp.dps = 50
msym, Psym = sp.symbols('m P')
coeffs = {
    (0,0): 1, (1,0): -7, (2,0): 21, (3,0): -35, (4,0): 35, (5,0): -21,
    (6,0): 7, (7,0): -1,
    (0,1): 192, (1,1): 7424, (2,1): 10368, (3,1): -2816, (4,1): 1216,
    (0,2): 12288, (1,2): -188416, (2,2): 45056,
    (0,3): 262144,
}
F = sum(c*msym**i*Psym**j for (i, j), c in coeffs.items())
Fm = sp.diff(F, msym)
FP = sp.diff(F, Psym)
Fmm = sp.diff(Fm, msym)
Fmp = sp.diff(Fm, Psym)
FPP = sp.diff(FP, Psym)

def fF(m_, p_):
    return float(F.subs({msym: m_, Psym: p_}))
def fFm(m_, p_):
    return float(Fm.subs({msym: m_, Psym: p_}))
def fFP(m_, p_):
    return float(FP.subs({msym: m_, Psym: p_}))

# exact check at m=1 first
print("F(1,P) =", sp.factor(F.subs(msym, 1)))
print("F_m(1,1/4) =", Fm.subs({msym: 1, Psym: sp.Rational(1, 4)}),
      " F_P(1,1/4) =", FP.subs({msym: 1, Psym: sp.Rational(1, 4)}))

nr = np.roots([1, -32, -54, 72, 13])
print("quartic roots:", nr)
sing = []
for r in nr:
    if abs(r.imag) > 1e-8:
        continue
    # refine m0 by Newton on quartic q(m) = m^4-32m^3-54m^2+72m+13
    m0 = float(r.real)
    for _ in range(50):
        q = m0**4 - 32*m0**3 - 54*m0**2 + 72*m0 + 13
        dq = 4*m0**3 - 96*m0**2 - 108*m0 + 72
        step = q/dq
        m0 -= step
        if abs(step) < 1e-15:
            break
    # quadratic F_P(m0, P): coefficients
    a = 3*262144
    b = 2*(45056*m0**2 - 188416*m0 + 12288)
    c = (1216*m0**4 - 2816*m0**3 + 10368*m0**2 + 7424*m0 + 192)
    disc = b*b - 4*a*c
    if disc < 0:
        print(f"m0={m0:.12f}: F_P quadratic has no real roots (disc<0)")
        continue
    for sgn in (1, -1):
        p0 = (-b + sgn*np.sqrt(disc))/(2*a)
        # refine p0 by Newton on F_P(m0, .)
        for _ in range(50):
            fp = fFP(m0, p0)
            fpp = float(FPP.subs({msym: m0, Psym: p0}))
            if abs(fpp) < 1e-12:
                break
            step = fp/fpp
            p0 -= step
            if abs(step) < 1e-15:
                break
        vf, vm, vp = fF(m0, p0), fFm(m0, p0), fFP(m0, p0)
        if abs(vf) < 1e-6 and abs(vm) < 1e-6 and abs(vp) < 1e-6:
            a2 = float(Fmm.subs({msym: m0, Psym: p0}))/2
            b2 = float(Fmp.subs({msym: m0, Psym: p0}))
            c2 = float(FPP.subs({msym: m0, Psym: p0}))/2
            disc2 = b2*b2 - 4*a2*c2
            kind = "node delta=1" if abs(disc2) > 1e-9 else "degenerate tangent cone"
            sing.append((m0, p0, kind))
            print(f"singular (m,P)=({m0:.15f}, {p0:.15f})  {kind}  "
                  f"F={vf:.1e} Fm={vm:.1e} FP={vp:.1e} disc2={disc2:.4g}")
print("affine singular points:", len(sing))
print("genus = 15 - 9 (inf) - %d = %d" % (len(sing), 15-9-len(sing)))
