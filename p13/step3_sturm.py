"""p13 step 3: Sturm-bound PROOF of the relation F(m, P) = 0.

Setup (all verified in step 0 / step 2):
  m = phi^2(q)/phi^2(q^13),  P = p_prod/64,  p_prod = 64 q^7 chi/chi(q^13)^13.
  F(m,P) = sum_{(i,j) in support} c_ij m^i P^j   (17 terms, exact ints)
  with F(m,P) = 0 as a q-series to order 600 and numerically to 1e-50.
  Pure parts: F(m,0) = -(m-1)^7,  F(0,P) = (1+64P)^3.

Sturm argument on Gamma0(52) (index mu = 84):
  K = eta(2z)^2 eta(26z)^2  (weight 2, trivial character, ord 7/3 at every cusp).
  K^18 has weight 36 and ord_{1/26}(K^18) = 42.
  m^i P^j has its ONLY pole at cusp 1/26, of order 6i+14j <= 42 (max over
  support: (7,0) and (0,3) both give 42).
  Hence every term K^18 m^i P^j is a HOLOMORPHIC modular form of weight 36,
  trivial character, on Gamma0(52) (Newman + D_eff square + Ligozat >= 0,
  all checked symbolically below), and
  K^18 * F(m,P) is a weight-36 form whose q-expansion vanishes to
  q^252 = the Sturm bound ceil(36*84/12) = 252  =>  it is identically 0
  => F(m,P) = 0 (K is not identically zero).
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

P_PRIME = 13
N = 400
# ---- exact coefficients from step 2 ----
coeffs = {
    (0,0): 1, (1,0): -7, (2,0): 21, (3,0): -35, (4,0): 35, (5,0): -21,
    (6,0): 7, (7,0): -1,
    (0,1): 192, (1,1): 7424, (2,1): 10368, (3,1): -2816, (4,1): 1216,
    (0,2): 12288, (1,2): -188416, (2,2): 45056,
    (0,3): 262144,
}
mons = sorted(coeffs.keys())

# ---- eta edicts ----
m_ed = {1: -4, 2: 10, 4: -4, 13: 4, 26: -10, 52: 4}    # m
P_ed = {1: -1, 2: 2, 4: -1, 13: 13, 26: -26, 52: 13}  # P = p_prod/64
K_ed = {2: 2, 26: 2}                                   # K
NLEV = 52
R_POW = 18   # K^18: ord 42 at 1/26, weight 36, Newman ok (18*56 = 1008 = 42*24)

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

print("index [SL2:Gamma0(52)] =", index_gamma0(52))
WT = 2*R_POW
print(f"weight {WT} Sturm bound = ceil({WT}*84/12) =", (WT*84 + 11)//12)

# ---- symbolic checks for every term K^18 m^i P^j ----
all_ok = True
for (i, j) in mons:
    ed = add_ed(scal_ed(K_ed, R_POW), scal_ed(m_ed, i), scal_ed(P_ed, j))
    wt = Fraction(sum(ed.values()), 2)
    s1 = sum(d*e for d, e in ed.items())
    s2 = sum((NLEV//d)*e for d, e in ed.items())
    D = Fraction(1)
    for d, e in ed.items():
        D *= Fraction(d)**e
    ords = {c: cusp_ord(ed, c) for c in [0]+divisors(NLEV)}
    ok = (wt == WT and s1 % 24 == 0 and s2 % 24 == 0 and is_square_frac(D)
          and all(o >= 0 for o in ords.values()))
    if not ok:
        all_ok = False
        print(f"  term m^{i}P^{j}: FAIL wt={wt} s1%24={s1%24} s2%24={s2%24} "
              f"D={D} minord={min(ords.values())}")
print(f"all {len(mons)} terms weight {WT}, Newman ok, trivial character, holomorphic at "
      "all cusps:", all_ok)
# pole orders of m^i P^j at cusp 1/26 (the only pole)
print("max pole order of m^i P^j at 1/26 over support:",
      max(6*i+14*j for (i, j) in mons),
      "; ord_{1/26}(K^18) =", cusp_ord(scal_ed(K_ed, R_POW), 26))
print("ord of K at every cusp:",
      {c: cusp_ord(K_ed, c) for c in [0]+divisors(NLEV)})

# ---- spot-check: K^18 series via eta.py equals theta-route series ----
_, CK = eta_quotient_qexp({2: 2, 26: 2}, N)
_, CK18 = eta_quotient_qexp({2: 2*R_POW, 26: 2*R_POW}, N)
print("K^18 via eta.py shift =", sum(Fraction(d*e, 24) for d, e in {2: 2*R_POW, 26: 2*R_POW}.items()))
print("K^18 = (K)^18 as series:", series_pow(CK, R_POW, N) == CK18)

# ---- q-expansion of S = K^18 * F(m,P) to order N ----
phi = phi_list(N)
phi2 = series_pow(phi, 2, N)
phip = [0]*(N+1)
for i, c in enumerate(phi_list(N//P_PRIME)):
    phip[P_PRIME*i] = c
m = series_div(phi2, series_pow(phip, 2, N), N)
chi = chi_list(N)
chip = [0]*(N+1)
for i, c in enumerate(chi_list(N//P_PRIME)):
    chip[P_PRIME*i] = c
ratio = series_div(chi, series_pow(chip, P_PRIME, N), N)   # P = q^7 * ratio
SHIFT_P = (P_PRIME*P_PRIME - 1)//24
mpow = [[1]+[0]*N]
for i in range(1, 8):
    mpow.append(series_mul(mpow[-1], m, N))
Rpow = [[1]+[0]*N]
for j in range(1, 4):
    Rpow.append(series_mul(Rpow[-1], ratio, N))
CK18 = series_pow(CK, R_POW, N)
SHIFT_K = R_POW * 7 // 3   # 42; K^18 shift = 18 * 7/3
total = [0]*(N+1)
for (i, j), c in coeffs.items():
    core = series_mul(series_mul(CK18, mpow[i], N), Rpow[j], N)
    shift = SHIFT_K + SHIFT_P*j   # K^18 shift 42, P^j shift 7j
    for r in range(N+1-shift):
        total[shift+r] += c*core[r]
B = (WT*84 + 11)//12
bad = [(r, total[r]) for r in range(B+1) if total[r] != 0]
print(f"Sturm check: coefficients q^0..q^{B} of K^18*F(m,P) all zero: {not bad}"
      + ("" if not bad else f" FIRST BAD {bad[0]}"))
nz = [(r, total[r]) for r in range(N+1) if total[r] != 0]
print(f"first nonzero coefficient within q^{N}: {nz[0] if nz else None}")
print("=> K^18*F(m,P) = 0 as a weight-36 form on Gamma0(52) by Sturm, so F(m,P)=0.")

# ---- structural form: F = (1+64P)^3 - (m-1)^7 + P*m*S(m,P) ----
msym, Psym = sp.symbols('m P')
F = sum(c*msym**i*Psym**j for (i, j), c in coeffs.items())
R = sp.expand(F - (1+64*Psym)**3 + (msym-1)**7)
print("mixed remainder R = F - (1+64P)^3 + (m-1)^7 =", sp.factor(R))
print("  R / (P*m) =", sp.factor(R/(Psym*msym)))
print("F(1, 0) =", F.subs({msym: 1, Psym: 0}), " (cusp inf: m=1, P=0)")
print("F(13, 1) =", F.subs({msym: 13, Psym: 1}), " (cusp 0: m=13, P=1 i.e. p_prod=64)")
