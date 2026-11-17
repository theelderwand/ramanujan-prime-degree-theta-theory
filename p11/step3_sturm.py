"""p11 step 3: Sturm-bound PROOF of the relation G(Y,P) = 0.

Setup (all verified in step 0 / step 2):
  Y = m^2 = phi^4(q)/phi^4(q^11),  P = p/32,  p = p_prod.
  G(Y,P) = sum_{(i,j) in triangle, i+j<=5} c_ij Y^i P^j   (21 terms, exact ints)
  with G(Y,P) = 0 as a q-series to order 600 and numerically to 1e-50.

Sturm argument on Gamma0(44) (index mu = 72):
  K = eta(2z)^2 eta(22z)^2  (weight 2, trivial character).
  K^25 has weight 50 and ord_{1/22}(K^25) = 50.
  Y^i P^j has its ONLY pole at cusp 1/22, of order 10(i+j) <= 50.
  Hence every term K^25 Y^i P^j is a HOLOMORPHIC modular form of weight 50,
  trivial character, on Gamma0(44) (Newman + D_eff square + Ligozat >= 0,
  all checked symbolically below), and
  K^25 * G(Y,P) is a weight-50 form whose q-expansion vanishes to
  q^300 = the Sturm bound ceil(50*72/12) = 300  =>  it is identically 0
  => G(Y,P) = 0 (K is not identically zero).
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
sys.path.insert(0, "/home/hatch/workspace/modular-forms/scripts")
from fractions import Fraction
from math import gcd
from theta import (phi_list, chi_list, series_pow, series_div, series_mul)
from eta import eta_quotient_qexp
from sturm import index_gamma0
import sympy as sp

N = 400
# ---- exact coefficients from step 2 ----
coeffs = {
    (0,0): 1, (1,0): -5, (2,0): 10, (3,0): -10, (4,0): 5, (5,0): -1,
    (0,1): 160, (1,1): 10624, (2,1): 19392, (3,1): 2432, (4,1): 160,
    (0,2): 10240, (1,2): -796672, (2,2): -55296, (3,2): -10240,
    (0,3): 327680, (1,3): 6914048, (2,3): 327680,
    (0,4): 5242880, (1,4): -5242880,
    (0,5): 33554432,
}
mons = sorted(coeffs.keys())

# ---- eta edicts ----
m_ed = {1: -4, 2: 10, 4: -4, 11: 4, 22: -10, 44: 4}     # m
Y_ed = {d: 2*e for d, e in m_ed.items()}                # Y = m^2
P_ed = {1: -1, 2: 2, 4: -1, 11: 11, 22: -22, 44: 11}    # P = p/32
K_ed = {2: 2, 22: 2}                                    # K
NLEV = 44

def add_ed(*eds):
    out = {}
    for e in eds:
        for d, x in e.items():
            out[d] = out.get(d, 0) + x
    return {d: x for d, x in out.items() if x}

def scal_ed(ed, k):
    return {d: k*x for d, x in ed.items()}

def is_square_frac(fr):
    fr = Fraction(fr)
    if fr < 0:
        return False
    return sp.integer_nthroot(fr.numerator, 2)[1] and sp.integer_nthroot(fr.denominator, 2)[1]

def cusp_ord(ed, c):
    if c == 0:
        return sum(Fraction(d*e, 24) for d, e in ed.items())
    fac = Fraction(NLEV, 24*gcd(c*c, NLEV))
    return fac * sum(Fraction(e*gcd(c, d)**2, d) for d, e in ed.items())

def divisors(n):
    return [d for d in range(1, n+1) if n % d == 0]

print("index [SL2:Gamma0(44)] =", index_gamma0(44))
print("weight 50 Sturm bound = ceil(50*72/12) =", 50*72//12)

# ---- symbolic checks for every term K^25 Y^i P^j ----
all_ok = True
for (i, j) in mons:
    ed = add_ed(scal_ed(K_ed, 25), scal_ed(Y_ed, i), scal_ed(P_ed, j))
    wt = Fraction(sum(ed.values()), 2)
    s1 = sum(d*e for d, e in ed.items())
    s2 = sum((NLEV//d)*e for d, e in ed.items())
    D = Fraction(1)
    for d, e in ed.items():
        D *= Fraction(d)**e
    ords = {c: cusp_ord(ed, c) for c in [0]+divisors(NLEV)}
    ok = (wt == 50 and s1 % 24 == 0 and s2 % 24 == 0 and is_square_frac(D)
          and all(o >= 0 for o in ords.values()))
    if not ok:
        all_ok = False
        print(f"  term Y^{i}P^{j}: FAIL wt={wt} s1%24={s1%24} s2%24={s2%24} "
              f"D={D} minord={min(ords.values())}")
print("all 21 terms weight 50, Newman ok, trivial character, holomorphic at "
      "all cusps:", all_ok)
# pole orders of Y^i P^j at cusp 1/22 (the only pole)
print("ord_{1/22}(Y^i P^j) = -10(i+j); ord_{1/22}(K^25) =",
      cusp_ord(scal_ed(K_ed, 25), 22))

# ---- spot-check: K^25 series via eta.py equals theta-route series ----
# K core series (shift 2): eta quotient {2:2,22:2}
_, CK = eta_quotient_qexp({2: 2, 22: 2}, N)
_, CK25 = eta_quotient_qexp({2: 50, 22: 50}, N)
print("K^25 via eta.py shift =", sum(Fraction(d*e, 24) for d, e in {2: 50, 22: 50}.items()))
print("K^25 = (K)^25 as series:", series_pow(CK, 25, N) == CK25)

# ---- q-expansion of S = K^25 * G(Y,P) to order N ----
phi = phi_list(N)
phi2 = series_pow(phi, 2, N)
phi11 = [0]*(N+1)
for i, c in enumerate(phi_list(N//11)):
    phi11[11*i] = c
m = series_div(phi2, series_pow(phi11, 2, N), N)
Y = series_mul(m, m, N)
chi = chi_list(N)
chi11 = [0]*(N+1)
for i, c in enumerate(chi_list(N//11)):
    chi11[11*i] = c
ratio = series_div(chi, series_pow(chi11, 11, N), N)   # P = q^5 * ratio
Ypow = [[1]+[0]*N]
for i in range(1, 6):
    Ypow.append(series_mul(Ypow[-1], Y, N))
Rpow = [[1]+[0]*N]
for j in range(1, 6):
    Rpow.append(series_mul(Rpow[-1], ratio, N))
CK25 = series_pow(CK, 25, N)
total = [0]*(N+1)
for (i, j), c in coeffs.items():
    core = series_mul(series_mul(CK25, Ypow[i], N), Rpow[j], N)
    shift = 50 + 5*j   # K^25 shift 50, P^j shift 5j
    for r in range(N+1-shift):
        total[shift+r] += c*core[r]
B = 300
bad = [(r, total[r]) for r in range(B+1) if total[r] != 0]
print(f"Sturm check: coefficients q^0..q^{B} of K^25*G(Y,P) all zero: {not bad}"
      + ("" if not bad else f" FIRST BAD {bad[0]}"))
nz = [(r, total[r]) for r in range(N+1) if total[r] != 0]
print(f"first nonzero coefficient within q^{N}: {nz[0] if nz else None}")
print("=> K^25*G(Y,P) = 0 as a weight-50 form on Gamma0(44) by Sturm, so G(Y,P)=0.")

# ---- structural form: G in terms of U = 1+p = 1+32P, V = Y-1 ----
Ysym, Psym = sp.symbols('Y P')
G = sum(c*Ysym**i*Psym**j for (i, j), c in coeffs.items())
U, V = sp.symbols('U V')
Guv = sp.expand(G.subs({Ysym: V+1, Psym: (U-1)/32}))
print("G in (U=1+p, V=m^2-1):")
print("  ", Guv)
print("  factor(Guv) =", sp.factor(Guv))
# remainder after removing U^5 - V^5
R = sp.expand(G - (1+32*Psym)**5 + (Ysym-1)**5)
print("mixed remainder R = G - (1+32P)^5 + (Y-1)^5 =", sp.factor(R))
print("  R / (P*Y) =", sp.factor(R/(Psym*Ysym)))
# cusp values
print("G(1, 0) =", G.subs({Ysym: 1, Psym: 0}), " (cusp inf: Y=1, P=0)")
print("G(11, 1) =", G.subs({Ysym: 11, Psym: 1}), " (cusp 0: Y=11, P=1 i.e. p=32)")
msym, psym = sp.symbols('m p')
F = sp.expand(G.subs({Ysym: msym**2, Psym: psym/32}))
print("F(m,p) =", F)
print("F(1,0) =", F.subs({msym: 1, psym: 0}), " F(11,32) =", F.subs({msym: 11, psym: 32}))
