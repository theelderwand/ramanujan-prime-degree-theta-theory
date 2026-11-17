"""p13 step 4f: settle the m0=33.545 fiber at 80-digit precision.

m0 refined as a root of the cubic m^3-31m^2-85m-13 to 80 digits; P0 as the
root of the quadratic F_P(m0,.) near -93.17, refined by Newton at 80 digits.
Then F(m0,P0) and F_m(m0,P0): if both ~0 (to ~1e-60), (m0,P0) is a node
(delta=1); if F_m is clearly nonzero, it is smooth ramification.
"""
import sympy as sp
import mpmath as mp

mp.mp.dps = 80
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

def ev(expr, m_, p_):
    return expr.subs({msym: m_, Psym: p_})

# refine m0: root of cubic near 33.54542973
m0 = mp.mpf('33.54542973')
for _ in range(100):
    q = m0**3 - 31*m0**2 - 85*m0 - 13
    dq = 3*m0**2 - 62*m0 - 85
    s = q/dq
    m0 -= s
    if abs(s) < mp.mpf(10)**(-70):
        break
print("m0 =", mp.nstr(m0, 30))
print("cubic(m0) =", mp.nstr(m0**3 - 31*m0**2 - 85*m0 - 13, 5))

# F_P(m0,P) quadratic: a P^2 + b P + c (mpmath)
a = mp.mpf(3*262144)
b = 2*(45056*m0**2 - 188416*m0 + 12288)
c = (1216*m0**4 - 2816*m0**3 + 10368*m0**2 + 7424*m0 + 192)
disc = b*b - 4*a*c
print("disc of F_P quadratic =", mp.nstr(disc, 10))
p0 = (-b - mp.sqrt(disc))/(2*a)
print("P0 candidate =", mp.nstr(p0, 30))
# Newton refine p0 on F_P(m0, .) = 0
for _ in range(100):
    fp = ev(FP, m0, p0)
    fpp = ev(FPP, m0, p0)
    s = fp/fpp
    p0 -= s
    if abs(s) < mp.mpf(10)**(-70):
        break
print("P0 refined =", mp.nstr(p0, 30))
print("F(m0,P0)   =", mp.nstr(ev(F, m0, p0), 5))
print("F_m(m0,P0) =", mp.nstr(ev(Fm, m0, p0), 5))
print("F_P(m0,P0) =", mp.nstr(ev(FP, m0, p0), 5))
# tangent cone discriminant at (m0,P0)
a2 = ev(Fmm, m0, p0)/2
b2 = ev(Fmp, m0, p0)
c2 = ev(FPP, m0, p0)/2
print("tangent-cone disc =", mp.nstr(b2**2 - 4*a2*c2, 5))
print("=> node (delta=1) iff F and F_m both ~0 and disc != 0")
