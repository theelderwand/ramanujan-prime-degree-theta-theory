"""p11 step 0: verify eta-quotient identities for phi and chi against theta.py,
compute exact series for m and p, and Ligozat cusp orders at level 44.

Conventions: q is the theta q (theta.py's q).  Claimed identities:
  phi(q) = eta(2z)^5 / (eta(z)^2 eta(4z)^2)      [e4, edict {2:5,1:-2,4:-2}]
  chi(q) = q^{1/24} eta(2z)^2 / (eta(z) eta(4z)) [e4, edict {2:2,1:-1,4:-1}, shift -1/24]
with q = e^{2 pi i z} the eta.py variable.
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
sys.path.insert(0, "/home/hatch/workspace/modular-forms/scripts")
from fractions import Fraction
from math import gcd
from theta import phi_list, chi_list, series_pow, series_div, series_mul, series_scale
from eta import eta_quotient_qexp
import mpmath as mp
from theta import phi_num, chi_num, p_prod_chi

N = 200

# --- 1. phi identity ---
shift, coeffs = eta_quotient_qexp({2: 5, 1: -2, 4: -2}, N)
print("phi eta quotient: shift =", shift)
print("  matches phi_list:", coeffs == phi_list(N))

# --- 2. chi identity: q^{1/24} * eta quotient has shift 0 ---
shift2, coeffs2 = eta_quotient_qexp({2: 2, 1: -1, 4: -1}, N)
print("chi eta quotient: shift =", shift2, "(+1/24 from prefactor -> 0)")
print("  matches chi_list:", coeffs2 == chi_list(N))

# --- 3. m and p series (exact integers) ---
N2 = 600
phi = phi_list(N2)
phi2 = series_pow(phi, 2, N2)
phi11 = [0]*(N2+1)
for i, c in enumerate(phi_list(N2//11)):
    if 11*i <= N2:
        phi11[11*i] = c
phi11_2 = series_pow(phi11, 2, N2)
m = series_div(phi2, phi11_2, N2)
print("m series first 12:", m[:12])
print("m[0] =", m[0], " m[1] =", m[1])

chi = chi_list(N2)
chi11 = [0]*(N2+1)
for i, c in enumerate(chi_list(N2//11)):
    if 11*i <= N2:
        chi11[11*i] = c
chi11_11 = series_pow(chi11, 11, N2)
# p = 32 q^5 chi / chi(q^11)^11
ratio = series_div(chi, chi11_11, N2)
p = [0]*(N2+1)
for i, c in enumerate(ratio):
    if i+5 <= N2:
        p[i+5] = 32*c
print("p series first 12:", p[:12])
print("p[5] =", p[5], "(expect 32)")

# cross-check p against eta quotient: p = 32 * eta-quotient {1:-1,2:2,4:-1,11:11,22:-22,44:11}, shift 5
shiftp, coeffp = eta_quotient_qexp({1: -1, 2: 2, 4: -1, 11: 11, 22: -22, 44: 11}, N2)
print("p eta quotient shift =", shiftp, "(expect 5)")
pet = [0]*(N2+1)
for i, c in enumerate(coeffp):
    if i+5 <= N2:
        pet[i+5] = 32*c
print("p via eta quotient matches theta.py p:", pet == p)

# cross-check m against eta quotient {1:-4,2:10,4:-4,11:4,22:-10,44:4}, shift 0
shiftm, coeffm = eta_quotient_qexp({1: -4, 2: 10, 4: -4, 11: 4, 22: -10, 44: 4}, N2)
print("m eta quotient shift =", shiftm, "(expect 0)")
print("m via eta quotient matches theta.py m:", coeffm == m)

# --- 4. Ligozat cusp orders at level 44 ---
def cusp_ord(edict, Nlev, c):
    if c == 0:
        return sum(Fraction(d*e, 24) for d, e in edict.items())
    fac = Fraction(Nlev, 24*gcd(c*c, Nlev))
    return fac * sum(Fraction(e*gcd(c, d)**2, d) for d, e in edict.items())

def divisors(n):
    return [d for d in range(1, n+1) if n % d == 0]

m_ed = {1: -4, 2: 10, 4: -4, 11: 4, 22: -10, 44: 4}
p_ed = {1: -1, 2: 2, 4: -1, 11: 11, 22: -22, 44: 11}
for name, ed in (("m", m_ed), ("p", p_ed)):
    print(f"--- Ligozat orders for {name} (level 44) ---")
    tot_pole = Fraction(0)
    for c in [0] + divisors(44):
        o = cusp_ord(ed, 44, c)
        lab = "inf" if c == 0 else f"1/{c}"
        print(f"  cusp {lab}: ord = {o}")
        if o < 0:
            tot_pole += -o
    print(f"  total pole order (sum of -ord over negative cusps) = {tot_pole}")
    # Newman checks
    s1 = sum(d*e for d, e in ed.items())
    s2 = sum((44//d)*e for d, e in ed.items())
    D = Fraction(1)
    for d, e in ed.items():
        D *= Fraction(d)**e
    print(f"  Newman: sum d*e = {s1} (mod24 {s1%24}), sum (44/d)*e = {s2} (mod24 {s2%24}), D_eff = {D}")

# --- 5. cusp-0 asymptotics (numerics): t = e^{-pi u}, u -> 0 ---
print("--- cusp 0 limits: t = e^{-pi u} ---")
mp.mp.dps = 30
for u in ("0.5", "0.3", "0.2", "0.1", "0.05"):
    u = mp.mpf(u)
    t = mp.e**(-mp.pi*u)
    mm = (phi_num(t, 30)/phi_num(t**11, 30))**2
    pp = p_prod_chi(11, t, 30)
    v = mp.e**(-mp.pi/u)   # local parameter at cusp 0 (up to 11th root)
    print(f"  u={u}: m={mp.nstr(mm,12)} (m-11={mp.nstr(mm-11,6)}), "
          f"p={mp.nstr(pp,12)} (p-32={mp.nstr(pp-32,6)}), "
          f"log|m-11|/log v = {mp.nstr(mp.log(abs(mm-11))/mp.log(v),6)}, "
          f"log|p-32|/log v = {mp.nstr(mp.log(abs(pp-32))/mp.log(v),6)}")
