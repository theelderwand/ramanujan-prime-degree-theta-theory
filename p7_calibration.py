"""Calibration (done right): the s_1 quadratic for p=7, derived symbolically.

s_1 = sum_k u_k^7/p = a^2/b + b^2/c + c^2/a  (cyclic order),
where a,b,c are the roots of r(xi) = xi^3 + 2(1+3p-X) xi^2 + p^2(p+4) xi - p^4,
X = phi^4(q)/phi^4(q^7).  The two cyclic sums S+ and S- are the two roots of
a quadratic whose coefficients are SYMMETRIC in a,b,c, hence polynomials in
e1,e2,e3, hence in (p,X).  Derive that quadratic with sympy, then verify
numerically that the dissection's s_1 is a root, at several q.

This validates the method Phase 2 needs for p=11 (where r_11 is unknown,
so the same derivation cannot be run yet -- the point here is to prove the
machinery works on the known case).
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
import sympy as sp
import mpmath as mp
from theta import phi_num, u_k_num, p_prod_chi

DPS = 50
mp.mp.dps = DPS

# ---- symbolic derivation ----
# Clear denominators: S+ = (a^3 c + b^3 a + c^3 b)/e3, S- = (a^3 b + b^3 c + c^3 a)/e3
a, b, c = sp.symbols('a b c')
e1, e2, e3 = sp.symbols('e1 e2 e3')
Np = a**3*c + b**3*a + c**3*b
Nm = a**3*b + b**3*c + c**3*a
# S+ + S- = (Np+Nm)/e3 ; S+*S- = Np*Nm/e3^2.  Symmetrize the numerators.
s_sum_num = sp.symmetrize(Np + Nm, [a, b, c])
s_prod_num = sp.symmetrize(Np * Nm, [a, b, c])
print("Np + Nm =", s_sum_num)
print("Np * Nm =", s_prod_num)
# substitute the elementary combinations:
#   (a+b+c) -> e1 = -2(1+3p-X), (ab+ac+bc) -> e2 = p^2(p+4), (abc) -> e3 = p^4
p_s, X_s = sp.symbols('p X')
E1 = -2*(1 + 3*p_s - X_s)
E2 = p_s**2*(p_s + 4)
E3 = p_s**4
subs = {a*b*c: E3, a*b + a*c + b*c: E2, a + b + c: E1}
S_sum = sp.expand(s_sum_num[0].subs(subs) / p_s**4)   # (Np+Nm)/e3
S_prod = sp.expand(s_prod_num[0].subs(subs) / p_s**8)  # Np*Nm/e3^2
print()
print("S+ + S- =", sp.factor(S_sum))
print("S+ * S- =", sp.factor(S_prod))
print()
print("Quadratic for s_1: T^2 - (S+ + S-) T + (S+ S-) = 0")
print("  coeff of T: -(", sp.factor(S_sum), ")")
print("  constant  :", sp.factor(S_prod))

# ---- numeric verification ----
def s1_numeric(q):
    p = p_prod_chi(7, q, DPS)
    X = (phi_num(q, DPS) / phi_num(q ** 7, DPS)) ** 4
    us = [u_k_num(7, k, q, DPS) for k in range(1, 4)]
    ws = [u ** 7 / p for u in us]
    return p, X, sum(ws)

print()
print("Numeric check (s_1 should satisfy the quadratic):")
for q in ("0.2", "0.3", "0.4", "0.5", "0.6"):
    qv = mp.mpf(q)
    p, X, s1 = s1_numeric(qv)
    # evaluate S_sum, S_prod numerically
    f_sum = sp.lambdify((p_s, X_s), S_sum, 'mpmath')
    f_prod = sp.lambdify((p_s, X_s), S_prod, 'mpmath')
    Ss = f_sum(p, X)
    Ps = f_prod(p, X)
    resid = s1**2 - Ss*s1 + Ps
    print(f"  q={q}: s1={mp.nstr(s1, 15)}, S+ + S-={mp.nstr(Ss, 15)}, "
          f"resid={mp.nstr(abs(resid), 5)} (scale {mp.nstr(abs(s1**2), 5)})")
    # which root is s1? roots of T^2 - Ss T + Ps
    disc = mp.sqrt(Ss**2 - 4*Ps)
    r1 = (Ss + disc) / 2
    r2 = (Ss - disc) / 2
    print(f"      roots: {mp.nstr(r1, 12)}, {mp.nstr(r2, 12)}; "
          f"|s1-r1|={mp.nstr(abs(s1-r1), 3)}, |s1-r2|={mp.nstr(abs(s1-r2), 3)}")
