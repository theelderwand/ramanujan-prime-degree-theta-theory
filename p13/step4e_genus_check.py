"""NOTE (2026-11-14 re-verification): this script's CONCLUSION IS WRONG.
It evaluates F_m at (m0, P0=-56.45), the double root of F_P, which is NOT on
the curve, and concludes the fibre over m0=33.545 is a smooth ramification
point (e=3, R contribution 2), giving g=3, R=10, R'=0. The actual point on
the curve is (33.545, -93.17), where F=F_m=F_P=0 with nonzero tangent-cone
discriminant: it is a NODE (see step4f_node_check.py, and the independent
50-digit check in paper/REVERIFICATION_LOG.md). Correct: 4 nodes, genus 2,
R=8, R'=4. Do not use this script's genus conclusion. Kept for audit trail.
"""
"""p13 step 4e: resolve the m0=33.545 fiber and confirm genus = 3.

Hypothesis: over m0 = 33.54542973 (a root of the cubic factor), F(m0,P) has a
TRIPLE root (smooth point, e=3), not a node.  Check: F_P(m0,P) quadratic has
discriminant ~ 0; F(m0,P0) ~ 0; F_m(m0,P0) != 0.
Also re-verify the three nodes' tangent cones and recompute R for
Riemann-Hurwitz on m: C -> P^1 (degree 3).
"""
import sympy as sp
import numpy as np

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

def fF(m_, p_):
    return float(F.subs({msym: m_, Psym: p_}))
def fFm(m_, p_):
    return float(Fm.subs({msym: m_, Psym: p_}))
def fFP(m_, p_):
    return float(FP.subs({msym: m_, Psym: p_}))

# refine m0 = 33.54542973 as root of cubic m^3-31m^2-85m-13
m0 = 33.54542973
for _ in range(60):
    q = m0**3 - 31*m0**2 - 85*m0 - 13
    dq = 3*m0**2 - 62*m0 - 85
    s = q/dq
    m0 -= s
    if abs(s) < 1e-16:
        break
print("m0 =", repr(m0))
# F_P(m0,P) quadratic coefficients: a P^2 + b P + c
a = 3*262144
b = 2*(45056*m0**2 - 188416*m0 + 12288)
c = (1216*m0**4 - 2816*m0**3 + 10368*m0**2 + 7424*m0 + 192)
print("F_P(m0,P) quadratic disc =", b*b - 4*a*c, " (expect ~0 for triple root)")
p0 = -b/(2*a)
print("P0 (double root of F_P) =", p0)
print("F(m0,P0) =", fF(m0, p0), " F_m(m0,P0) =", fFm(m0, p0),
      " F_P(m0,P0) =", fFP(m0, p0))
# F(m0,P) should be ~ 262144*(P-P0)^3: check roots
roots = np.roots([262144,
                  (45056*m0**2 - 188416*m0 + 12288),
                  (1216*m0**4 - 2816*m0**3 + 10368*m0**2 + 7424*m0 + 192),
                  (-(m0-1)**7)])
print("roots of F(m0,P):", roots)
print("=> fiber over m0=33.545: single point, e=3, smooth (F_m != 0). Contributes 2 to R.")
print()
print("Riemann-Hurwitz on m (deg 3):")
print("  m=0: e=3 -> 2 ; m=1 node -> 0 ; cubic roots -2.383, -0.163 nodes -> 0 ;")
print("  m0=33.545 e=3 -> 2 ; quartic 4 simple ramification -> 4 ; m=inf e=3 -> 2")
print("  R = 2+0+0+2+4+2 = 10 ; 2g-2 = -6+10 = 4 => g = 3")
print("Plane formula: 15 - 9 (inf, unibranch (4,7)) - 3 (nodes) = 3. CONSISTENT.")
print("X0(52) check: g=5, degree-2 map, 2*5-2 = 2*(2*3-2) + R' => R' = 0: involution fixed-point-free.")
