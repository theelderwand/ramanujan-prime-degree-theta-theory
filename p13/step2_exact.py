"""p13 step 2: recover the EXACT integer relation.

From step 1 (mod 999983): the relation first appears at box (7,4), nullspace
dim 2 (the relation plus its multiple by P).  The first null vector is
supported on 17 monomials:
  m^0..7 P^0,  m^0..4 P^1,  m^0..2 P^2,  m^0 P^3
with P^3 coefficient 1 (mod prime).  Here we solve the exact rational
nullspace with sympy on exact integer series over that support, then verify
exactly to order 600 and numerically at several q to 50 digits.

Note: unlike p=11, odd powers of m occur, so the relation is genuinely in
(m, P), not in (Y=m^2, P).
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
from theta import (phi_list, chi_list, series_pow, series_div, series_mul)
import sympy as sp
import mpmath as mp
from theta import phi_num, p_prod_chi

P_PRIME = 13
N = 600
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
ratio = series_div(chi, series_pow(chip, P_PRIME, N), N)
SHIFT_P = (P_PRIME*P_PRIME - 1)//24
P = [0]*(N+1)                    # P = p_prod/64
for i, c in enumerate(ratio):
    if i+SHIFT_P <= N:
        P[i+SHIFT_P] = c

# support from step 1 (mod-prime null vector)
mons = ([(i, 0) for i in range(8)] + [(i, 1) for i in range(5)]
        + [(i, 2) for i in range(3)] + [(0, 3)])
print("support monomials (i=m-exp, j=P-exp):", mons, " count:", len(mons))
mpow = [[1]+[0]*N]
for i in range(1, 8):
    mpow.append(series_mul(mpow[-1], m, N))
Ppow = [[1]+[0]*N]
for j in range(1, 4):
    Ppow.append(series_mul(Ppow[-1], P, N))
mon_series = {}
for (i, j) in mons:
    mon_series[(i, j)] = series_mul(mpow[i], Ppow[j], N)

# exact nullspace: use coefficients t^0..t^R-1, R = 90
R = 90
A = sp.Matrix([[mon_series[(i, j)][r] for (i, j) in mons] for r in range(R)])
ns = A.nullspace()
print("nullspace dimension (exact, R=90):", len(ns))
assert len(ns) == 1
v = ns[0]
from math import gcd as igcd
from functools import reduce
dens = [sp.denom(x) for x in v]
L = reduce(lambda a, b: a*b//igcd(a, b), dens, 1)
ints = [int(x*L) for x in v]
g = reduce(igcd, [abs(x) for x in ints if x != 0])
ints = [x//g for x in ints]
idx03 = mons.index((0, 3))
if ints[idx03] < 0:
    ints = [-x for x in ints]
print("primitive integer coefficients (P^3 coeff = %d):" % ints[idx03])
for (i, j), c in zip(mons, ints):
    if c:
        print(f"  m^{i} P^{j}: {c}")

# exact verification to order N: sum c_ij m^i P^j == 0
tot = [0]*(N+1)
for (i, j), c in zip(mons, ints):
    if c:
        s = mon_series[(i, j)]
        for r in range(N+1):
            tot[r] += c*s[r]
nz = [(r, c) for r, c in enumerate(tot) if c != 0]
print("exact check to order", N, ": all zero =", not nz,
      ("first nonzero at %s" % (nz[0],) if nz else ""))

# numeric verification at several q, 50 digits
print("numeric verification (50 digits):")
mp.mp.dps = 50
for qstr in ("0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.85"):
    q = mp.mpf(qstr)
    mm = (phi_num(q, 50)/phi_num(q**13, 50))**2
    pp = p_prod_chi(13, q, 50)/64
    tot = mp.mpf(0)
    scale = mp.mpf(0)
    for (i, j), c in zip(mons, ints):
        t = c * mm**i * pp**j
        tot += t
        scale = max(scale, abs(t))
    print(f"  q={qstr}: |F| = {mp.nstr(abs(tot), 5)}, max term = {mp.nstr(scale, 5)}, "
          f"rel = {mp.nstr(abs(tot)/scale, 5)}")
# flagship point q1 = e^{-pi/sqrt(13)}: m = sqrt(13), p_prod = 2^{3} * G_1/G_169^13... use chi form
q1 = mp.e**(-mp.pi/mp.sqrt(13))
mm = (phi_num(q1, 50)/phi_num(q1**13, 50))**2
pp = p_prod_chi(13, q1, 50)/64
tot = sum(c * mm**i * pp**j for (i, j), c in zip(mons, ints))
print(f"  flagship q1=e^(-pi/sqrt13): m={mp.nstr(mm,20)} (sqrt13={mp.nstr(mp.sqrt(13),20)}), "
      f"P={mp.nstr(pp,20)}, |F| = {mp.nstr(abs(tot), 5)}")

# print the polynomial nicely in (m, P)
msym, Psym = sp.symbols('m P')
F = sum(c * msym**i * Psym**j for (i, j), c in zip(mons, ints) if c)
print("F(m,P) =", sp.expand(F))
print("factor F(m,P):", sp.factor(F))
# pure parts
print("F(m,0) =", sp.factor(F.subs(Psym, 0)))
print("F(0,P) =", sp.factor(F.subs(msym, 0)))
# Y = m^2 substitution: does it factor / simplify?
Ysym = sp.symbols('Y')
print("F in terms of Y=m^2 where possible: even-m part + odd-m part")
Fe = sum(c * Ysym**(i//2) * Psym**j for (i, j), c in zip(mons, ints) if c and i % 2 == 0)
Fo = sum(c * Ysym**(i//2) * Psym**j for (i, j), c in zip(mons, ints) if c and i % 2 == 1)
print("  even part (as poly in Y,P):", sp.expand(Fe))
print("  odd part / m (as poly in Y,P):", sp.expand(Fo))
