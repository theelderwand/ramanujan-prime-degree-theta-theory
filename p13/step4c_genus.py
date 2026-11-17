"""NOTE (2026-11-14 re-verification): SUPERSEDED. This script finds 0 affine
singular points and reports genus 6, which is wrong: it misses all four nodes.
The correct genus analysis is in step4f_node_check.py (4th node) plus the
independent 50-digit node verification in paper/REVERIFICATION_LOG.md.
Correct: 4 affine nodes, genus = 15 - 9 - 4 = 2. Kept for audit trail.
"""
"""p13 step 4c: singular points (fixed) and genus.

Singular m-coordinates are roots of the square-free quartic
(m-1)(m^3-31m^2-85m-13) = m^4-32m^3-54m^2+72m+13.
For each m0, a singular point has P a DOUBLE root of F(m0, .), hence a root
of the quadratic F_P(m0, .) that also kills F and F_m.  Classify via the
tangent cone (quadratic part): distinct tangents => node, delta=1.
Genus: plane model degree 7, arith genus 15; delta_inf = 9 (unibranch (4,7));
g = 15 - 9 - sum(delta_affine).  Cross-check by Riemann-Hurwitz on m: C->P^1
(degree 3) using disc_P(F), and by the degree-2 map X0(52) -> C.
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

# quartic coefficients (correct): m^4 - 32 m^3 - 54 m^2 + 72 m + 13
nr = np.roots([1, -32, -54, 72, 13])
print("quartic roots (numpy):", nr)
m_roots = []
for r in nr:
    if abs(r.imag) < 1e-8:
        m_roots.append(mp.mpf(float(r.real)))
    else:
        m_roots.append(mp.mpc(float(r.real), float(r.imag)))
print("m-roots:", [mp.nstr(x, 12) for x in m_roots])

sing = []
for m0 in m_roots:
    # quadratic F_P(m0, P) = 0
    q = sp.Poly(FP.subs(msym, complex(m0)), Psym)
    qc = [complex(c) for c in q.all_coeffs()]
    proots = np.roots(qc)
    for pr in proots:
        p0 = mp.mpc(float(pr.real), float(pr.imag))
        # Newton refine on (F, FP) as a system in (m, P) near (m0, p0)
        mm, pp = m0, p0
        for _ in range(40):
            f0 = complex(F.subs({msym: complex(mm), Psym: complex(pp)}))
            f1 = complex(FP.subs({msym: complex(mm), Psym: complex(pp)}))
            a = complex(Fm.subs({msym: complex(mm), Psym: complex(pp)}))
            b = complex(FP.subs({msym: complex(mm), Psym: complex(pp)}))
            # Jacobian of (F, FP) wrt (m, P): [[Fm, FP],[Fmp, FPP]]
            c_ = complex(Fmp.subs({msym: complex(mm), Psym: complex(pp)}))
            d_ = complex(FPP.subs({msym: complex(mm), Psym: complex(pp)}))
            det = a*d_ - b*c_
            if abs(det) < 1e-30:
                break
            dm = (-f0*d_ + f1*b)/det
            dp = (f0*c_ - f1*a)/det
            mm = mm + dm; pp = pp + dp
            if abs(dm) < 1e-35 and abs(dp) < 1e-35:
                break
        vm = complex(Fm.subs({msym: complex(mm), Psym: complex(pp)}))
        vp = complex(FP.subs({msym: complex(mm), Psym: complex(pp)}))
        vf = complex(F.subs({msym: complex(mm), Psym: complex(pp)}))
        if abs(vm) < 1e-20 and abs(vp) < 1e-20 and abs(vf) < 1e-20:
            a2 = complex(Fmm.subs({msym: complex(mm), Psym: complex(pp)}))/2
            b2 = complex(Fmp.subs({msym: complex(mm), Psym: complex(pp)}))
            c2 = complex(FPP.subs({msym: complex(mm), Psym: complex(pp)}))/2
            disc2 = b2**2 - 4*a2*c2
            kind = "node delta=1" if abs(disc2) > 1e-12 else "degenerate tangent cone"
            sing.append((mm, pp, kind))
            print(f"singular (m,P)=({mp.nstr(mm,18)}, {mp.nstr(pp,18)})  {kind}")
print("affine singular points:", len(sing))

# m=0 fiber: F(0,P) = (1+64P)^3, triple root P=-1/64 -> totally ramified, not singular
print("F(0,P) =", sp.factor(F.subs(msym, 0)))
print("F_m(0,-1/64) =", Fm.subs({msym: 0, Psym: sp.Rational(-1, 64)}), "(nonzero => smooth, e=3)")
# quartic factor of disc: simple ramification points
q4 = 27*msym**4 - 212*msym**3 + 1194*msym**2 - 2756*msym + 4563
print("quartic disc factor roots (numpy):", np.roots([27, -212, 1194, -2756, 4563]))
print("genus: 15 - 9 (inf) - %d (affine nodes) = %d" % (len(sing), 15-9-len(sing)))
