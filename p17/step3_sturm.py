"""p17 step 3: Sturm-bound PROOF of the relation F(m, P) = 0.

Setup (all verified in step 0 / step 2):
  m = phi^2(q)/phi^2(q^17),  P = p_prod/256,  p_prod = 256 q^12 chi/chi(q^17)^17.
  F(m,P) = sum_{(i,j): i+3j<=12} c_ij m^i P^j   (35 terms, exact ints)
  with F(m,P) = 0 as a q-series to order 600 and numerically to 1e-50.
  Pure parts: F(m,0) = (m-1)^12,  F(0,P) = (256P-1)^4.

Sturm argument on Gamma0(68) (index mu = 108):
  K = eta(2z)^2 eta(34z)^2  (weight 2, trivial character, ord 3 at every cusp).
  K^32 has weight 64 and ord_{1/34}(K^32) = 96.
  m^i P^j has its ONLY pole at cusp 1/34, of order 8i+24j = 8(i+3j) <= 96
  (the support triangle i+3j<=12 is exactly the constant-pole-order triangle).
  Hence every term K^32 m^i P^j is a HOLOMORPHIC modular form of weight 64,
  trivial character, on Gamma0(68) (Newman + D_eff square + Ligozat >= 0,
  all checked symbolically below), and
  K^32 * F(m,P) is a weight-64 form whose q-expansion vanishes to
  q^576 = the Sturm bound ceil(64*108/12) = 576  =>  it is identically 0
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

P_PRIME = 17
N = 800
# ---- exact coefficients from step 2 ----
coeffs = {
    (0,0): 1, (1,0): -12, (2,0): 66, (3,0): -220, (4,0): 495, (5,0): -792,
    (6,0): 924, (7,0): -792, (8,0): 495, (9,0): -220, (10,0): 66, (11,0): -12,
    (12,0): 1,
    (0,1): -1024, (1,1): -391168, (2,1): -5329920, (3,1): -10408960,
    (4,1): 5120, (5,1): 19456, (6,1): -744448, (7,1): 99328, (8,1): -24576,
    (9,1): -1024,
    (0,2): 393216, (1,2): -195100672, (2,2): 919994368, (3,2): -192806912,
    (4,2): 107347968, (5,2): -19464192, (6,2): 393216,
    (0,3): -67108864, (1,3): -6643777536, (2,3): -3690987520, (3,3): -67108864,
    (0,4): 4294967296,
}
mons = sorted(coeffs.keys())

# ---- eta edicts ----
m_ed = {1: -4, 2: 10, 4: -4, 17: 4, 34: -10, 68: 4}    # m
P_ed = {1: -1, 2: 2, 4: -1, 17: 17, 34: -34, 68: 17}  # P = p_prod/256
K_ed = {2: 2, 34: 2}                                   # K
NLEV = 68
R_POW = 32   # K^32: ord 96 at 1/34, weight 64, Newman ok

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

print("index [SL2:Gamma0(68)] =", index_gamma0(68))
WT = 2*R_POW
print(f"weight {WT} Sturm bound = ceil({WT}*108/12) =", (WT*108 + 11)//12)

# ---- symbolic checks for every term K^32 m^i P^j ----
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
print("max pole order of m^i P^j at 1/34 over support:",
      max(8*i+24*j for (i, j) in mons),
      "; ord_{1/34}(K^32) =", cusp_ord(scal_ed(K_ed, R_POW), 34))
print("ord of K at every cusp:",
      {c: cusp_ord(K_ed, c) for c in [0]+divisors(NLEV)})

# ---- spot-check: K^32 series via eta.py equals theta-route series ----
_, CK = eta_quotient_qexp({2: 2, 34: 2}, N)
_, CK32 = eta_quotient_qexp({2: 2*R_POW, 34: 2*R_POW}, N)
print("K^32 via eta.py shift =", sum(Fraction(d*e, 24) for d, e in {2: 2*R_POW, 34: 2*R_POW}.items()))
print("K^32 = (K)^32 as series:", series_pow(CK, R_POW, N) == CK32)

# ---- q-expansion of S = K^32 * F(m,P) to order N ----
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
ratio = series_div(chi, series_pow(chip, P_PRIME, N), N)   # P = q^12 * ratio
SHIFT_P = (P_PRIME*P_PRIME - 1)//24
mpow = [[1]+[0]*N]
for i in range(1, 13):
    mpow.append(series_mul(mpow[-1], m, N))
Rpow = [[1]+[0]*N]
for j in range(1, 5):
    Rpow.append(series_mul(Rpow[-1], ratio, N))
CK32 = series_pow(CK, R_POW, N)
SHIFT_K = R_POW * 3   # 96; K^32 shift = 32 * 3
total = [0]*(N+1)
for (i, j), c in coeffs.items():
    core = series_mul(series_mul(CK32, mpow[i], N), Rpow[j], N)
    shift = SHIFT_K + SHIFT_P*j   # K^32 shift 96, P^j shift 12j
    for r in range(N+1-shift):
        total[shift+r] += c*core[r]
B = (WT*108 + 11)//12
bad = [(r, total[r]) for r in range(B+1) if total[r] != 0]
print(f"Sturm check: coefficients q^0..q^{B} of K^32*F(m,P) all zero: {not bad}"
      + ("" if not bad else f" FIRST BAD {bad[0]}"))
nz = [(r, total[r]) for r in range(N+1) if total[r] != 0]
print(f"first nonzero coefficient within q^{N}: {nz[0] if nz else None}")
print("=> K^32*F(m,P) = 0 as a weight-64 form on Gamma0(68) by Sturm, so F(m,P)=0.")

# ---- structural form: F = (m-1)^12 + (256P-1)^4 - 1 + P*m*S(m,P) ----
msym, Psym = sp.symbols('m P')
F = sum(c*msym**i*Psym**j for (i, j), c in coeffs.items())
R = sp.expand(F - (msym-1)**12 - (256*Psym-1)**4 + 1)
print("mixed remainder R = F - (m-1)^12 - (256P-1)^4 + 1 =", sp.factor(R))
print("  R / (P*m) =", sp.factor(R/(Psym*msym)))
print("F(1, 0) =", F.subs({msym: 1, Psym: 0}), " (cusp inf: m=1, P=0)")
print("F(17, 1) =", F.subs({msym: 17, Psym: 1}), " (cusp 0: m=17, P=1 i.e. p_prod=256)")
