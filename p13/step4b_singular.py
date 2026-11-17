"""p13 step 4b: singular points and genus (robust numerics).

From step 4: gcd(res_P(F,F_m), res_P(F,F_P)) = const*(m-1)^2*(m^3-31m^2-85m-13)^2.
Square-free part: (m-1)*(m^3-31m^2-85m-13).  Find its roots with numpy,
refine, and for each m0 solve F(m0,P)=0, keeping the points where
F_m = F_P = 0.  Determine each singularity type via the tangent cone
(quadratic part of F at the point): two distinct tangents => node, delta=1;
one tangent (double) => cusp-like, delta computed from multiplicity/branches.
Also: discriminant disc_P(F) for the Riemann-Hurwitz check on m: C -> P^1.
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
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

# square-free singular locus in m: (m-1)*(m^3-31m^2-85m-13)
sqf = (msym-1)*(msym**3 - 31*msym**2 - 85*msym - 13)
print("square-free gcd:", sp.expand(sqf))
nr = np.roots([1, -32, -54, -13])  # (m-1)(m^3-31m^2-85m-13) expanded: m^4-32m^3-54m^2+72m+13
print("numpy roots of sqf:", nr)
# refine each real root with mpmath findroot on the cubic factor + m=1
def cubic(x):
    return x**3 - 31*x**2 - 85*x - 13
m_roots = [mp.mpf(1)]
for r in nr:
    if abs(r.imag) < 1e-8:
        m_roots.append(mp.findroot(cubic, mp.mpf(float(r.real))))
print("refined m-roots:", [mp.nstr(x, 15) for x in m_roots])

sing = []
for m0 in m_roots:
    # F(m0, P) as polynomial in P; find its roots numerically
    Fm0 = sp.Poly(F.subs(msym, float(m0)), Psym)
    coeffsP = [complex(c) for c in Fm0.all_coeffs()]
    proots = np.roots(coeffsP)
    for pr in proots:
        if abs(pr.imag) < 1e-6:
            p0 = mp.mpf(float(pr.real))
            # refine (m0, p0) jointly by Newton on (F, Fm) or (F, FP)
            def fvec(X):
                mm, pp = X[0], X[1]
                return (float(F.subs({msym: mm, Psym: pp})),
                        float(Fm.subs({msym: mm, Psym: pp})))
            # Newton with mpmath
            m_, p_ = m0, p0
            for _ in range(30):
                f0 = F.subs({msym: m_, Psym: p_})
                f1 = Fm.subs({msym: m_, Psym: p_})
                f2 = FP.subs({msym: m_, Psym: p_})
                # Jacobian of (F, Fm): [[Fm, FP],[Fmm, Fmp]]
                a = Fm.subs({msym: m_, Psym: p_}); b = FP.subs({msym: m_, Psym: p_})
                c = Fmm.subs({msym: m_, Psym: p_}); d = Fmp.subs({msym: m_, Psym: p_})
                det = a*d - b*c
                if abs(det) < mp.mpf(10)**(-40):
                    break
                dm = (-f0*d + f1*b)/det
                dp = (f0*c - f1*a)/det
                m_ += dm; p_ += dp
                if abs(dm) < mp.mpf(10)**(-40) and abs(dp) < mp.mpf(10)**(-40):
                    break
            vm = Fm.subs({msym: m_, Psym: p_})
            vp = FP.subs({msym: m_, Psym: p_})
            vf = F.subs({msym: m_, Psym: p_})
            if abs(vm) < mp.mpf(10)**(-30) and abs(vp) < mp.mpf(10)**(-30):
                # tangent cone: quadratic part coefficients
                a2 = Fmm.subs({msym: m_, Psym: p_})/2
                b2 = Fmp.subs({msym: m_, Psym: p_})
                c2 = FPP.subs({msym: m_, Psym: p_})/2
                disc2 = b2**2 - 4*a2*c2
                kind = "node (delta=1)" if disc2 != 0 else "cusp-like (tangent cone degenerate)"
                sing.append((m_, p_, kind))
                print(f"singular (m,P)=({mp.nstr(m_,20)}, {mp.nstr(p_,20)})  {kind}  disc2={mp.nstr(disc2,6)}")
print("affine singular points:", len(sing))

# discriminant of F as cubic in P
print("computing disc_P(F)...")
disc = sp.discriminant(F, Psym)
print("disc_P(F) degree in m:", sp.Poly(disc, msym).degree())
print("disc_P(F) factorisation:")
print("  ", sp.factor(disc))
