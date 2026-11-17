"""p13 step 0: verify eta-quotient identities for phi and chi, compute exact
series for m and P, and Ligozat cusp orders at level 52.

Setup (generalising p11):
  m = phi^2(q)/phi^2(q^13), eta quotient {1:-4,2:10,4:-4,13:4,26:-10,52:4}, shift 0
  p_prod = 64 q^7 chi(q)/chi(q^13)^13   [general formula: 2^6 q^7 chi/chi^13]
  P = p_prod/64 = q^7 * chi/chi(q^13)^13, eta quotient {1:-1,2:2,4:-1,13:13,26:-26,52:13}, shift 7
  Y = m^2
Level 52, index [SL2:Gamma0(52)] = 84.
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
sys.path.insert(0, "/home/hatch/workspace/modular-forms/scripts")
from fractions import Fraction
from math import gcd
from theta import phi_list, chi_list, series_pow, series_div, series_mul
from eta import eta_quotient_qexp
from sturm import index_gamma0
import mpmath as mp
from theta import phi_num, chi_num, p_prod_chi

P_PRIME = 13
NLEV = 4 * P_PRIME
print("index [SL2:Gamma0(%d)] =" % NLEV, index_gamma0(NLEV))

N = 200
# --- 1. phi identity ---
shift, coeffs = eta_quotient_qexp({2: 5, 1: -2, 4: -2}, N)
print("phi eta quotient: shift =", shift)
print("  matches phi_list:", coeffs == phi_list(N))

# --- 2. chi identity ---
shift2, coeffs2 = eta_quotient_qexp({2: 2, 1: -1, 4: -1}, N)
print("chi eta quotient: shift =", shift2, "(+1/24 from prefactor -> 0)")
print("  matches chi_list:", coeffs2 == chi_list(N))

# --- 3. m and P series (exact integers) ---
N2 = 600
phi = phi_list(N2)
phi2 = series_pow(phi, 2, N2)
phip = [0]*(N2+1)
for i, c in enumerate(phi_list(N2//P_PRIME)):
    if P_PRIME*i <= N2:
        phip[P_PRIME*i] = c
phip_2 = series_pow(phip, 2, N2)
m = series_div(phi2, phip_2, N2)
print("m series first 12:", m[:12])
print("m[0] =", m[0], " m[1] =", m[1])

chi = chi_list(N2)
chip = [0]*(N2+1)
for i, c in enumerate(chi_list(N2//P_PRIME)):
    if P_PRIME*i <= N2:
        chip[P_PRIME*i] = c
chip_p = series_pow(chip, P_PRIME, N2)
# P = q^7 * chi / chi(q^13)^13
ratio = series_div(chi, chip_p, N2)
SHIFT_P = (P_PRIME*P_PRIME - 1)//24
print("SHIFT_P =", SHIFT_P, "(expect 7)")
P = [0]*(N2+1)
for i, c in enumerate(ratio):
    if i+SHIFT_P <= N2:
        P[i+SHIFT_P] = c
print("P series first 12:", P[:12])
print("P[%d] =" % SHIFT_P, P[SHIFT_P], "(expect 1)")

# cross-check P against eta quotient
p_ed = {1: -1, 2: 2, 4: -1, 13: 13, 26: -26, 52: 13}
shiftp, coeffp = eta_quotient_qexp(p_ed, N2)
print("P eta quotient shift =", shiftp, "(expect 7)")
pet = [0]*(N2+1)
for i, c in enumerate(coeffp):
    if i+SHIFT_P <= N2:
        pet[i+SHIFT_P] = c
print("P via eta quotient matches theta.py P:", pet == P)

# cross-check m against eta quotient
m_ed = {1: -4, 2: 10, 4: -4, 13: 4, 26: -10, 52: 4}
shiftm, coeffm = eta_quotient_qexp(m_ed, N2)
print("m eta quotient shift =", shiftm, "(expect 0)")
print("m via eta quotient matches theta.py m:", coeffm == m)

# --- 4. Ligozat cusp orders at level 52 ---
def cusp_ord(edict, Nlev, c):
    if c == 0:
        return sum(Fraction(d*e, 24) for d, e in edict.items())
    fac = Fraction(Nlev, 24*gcd(c*c, Nlev))
    return fac * sum(Fraction(e*gcd(c, d)**2, d) for d, e in edict.items())

def divisors(n):
    return [d for d in range(1, n+1) if n % d == 0]

for name, ed in (("m", m_ed), ("P", p_ed)):
    print(f"--- Ligozat orders for {name} (level 52) ---")
    tot_pole = Fraction(0)
    for c in [0] + divisors(NLEV):
        o = cusp_ord(ed, NLEV, c)
        lab = "inf" if c == 0 else f"1/{c}"
        print(f"  cusp {lab}: ord = {o}")
        if o < 0:
            tot_pole += -o
    print(f"  total pole order (degree) = {tot_pole}")
    s1 = sum(d*e for d, e in ed.items())
    s2 = sum((NLEV//d)*e for d, e in ed.items())
    D = Fraction(1)
    for d, e in ed.items():
        D *= Fraction(d)**e
    print(f"  Newman: sum d*e = {s1} (mod24 {s1%24}), sum (52/d)*e = {s2} (mod24 {s2%24}), D_eff = {D}")

# --- 5. cusp-0 asymptotics (numerics) ---
print("--- cusp 0 limits: t = e^{-pi u} ---")
mp.mp.dps = 30
for u in ("0.5", "0.3", "0.2", "0.1", "0.05"):
    u = mp.mpf(u)
    t = mp.e**(-mp.pi*u)
    mm = (phi_num(t, 30)/phi_num(t**13, 30))**2
    pp = p_prod_chi(13, t, 30)
    v = mp.e**(-mp.pi/u)
    print(f"  u={u}: m={mp.nstr(mm,12)} (m-13={mp.nstr(mm-13,6)}), "
          f"p={mp.nstr(pp,12)} (p-64={mp.nstr(pp-64,6)}), "
          f"log|m-13|/log v = {mp.nstr(mp.log(abs(mm-13))/mp.log(v),6)}, "
          f"log|p-64|/log v = {mp.nstr(mp.log(abs(pp-64))/mp.log(v),6)}")
