"""p11 step 2: recover the EXACT integer relation.

From step 1 (mod 999983): the relation lives on the triangle
  Y^i P^j,  i + j <= 5,   Y = m^2,  P = p/32,
with P^5 coefficient 1 (mod prime).  Here we solve the exact rational
nullspace with sympy on exact integer series, then verify exactly to
order 600 and numerically at several q to 50 digits.
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
from theta import (phi_list, chi_list, series_pow, series_div, series_mul)
import sympy as sp
import mpmath as mp
from theta import phi_num, p_prod_chi

N = 600
phi = phi_list(N)
phi2 = series_pow(phi, 2, N)
phi11 = [0]*(N+1)
for i, c in enumerate(phi_list(N//11)):
    phi11[11*i] = c
m = series_div(phi2, series_pow(phi11, 2, N), N)
Y = series_mul(m, m, N)          # Y = m^2
chi = chi_list(N)
chi11 = [0]*(N+1)
for i, c in enumerate(chi_list(N//11)):
    chi11[11*i] = c
ratio = series_div(chi, series_pow(chi11, 11, N), N)
P = [0]*(N+1)                    # P = p/32
for i, c in enumerate(ratio):
    if i+5 <= N:
        P[i+5] = c

# triangle monomials Y^i P^j, i+j<=5
mons = [(i, j) for j in range(6) for i in range(6-j)]
print("monomials (i=Y-exp, j=P-exp):", mons)
Ypow = [[1]+[0]*N]
for i in range(1, 6):
    Ypow.append(series_mul(Ypow[-1], Y, N))
Ppow = [[1]+[0]*N]
for j in range(1, 6):
    Ppow.append(series_mul(Ppow[-1], P, N))
mon_series = {}
for (i, j) in mons:
    mon_series[(i, j)] = series_mul(Ypow[i], Ppow[j], N)

# exact nullspace: use coefficients t^0..t^R-1, R = 90
R = 90
A = sp.Matrix([[mon_series[(i, j)][r] for (i, j) in mons] for r in range(R)])
ns = A.nullspace()
print("nullspace dimension (exact, R=90):", len(ns))
assert len(ns) == 1
v = ns[0]
# scale to primitive integer vector with P^5 coeff positive
from math import gcd as igcd
from functools import reduce
dens = [sp.denom(x) for x in v]
L = reduce(lambda a, b: a*b//igcd(a, b), dens, 1)
ints = [int(x*L) for x in v]
g = reduce(igcd, [abs(x) for x in ints if x != 0])
ints = [x//g for x in ints]
idx55 = mons.index((0, 5))
if ints[idx55] < 0:
    ints = [-x for x in ints]
print("primitive integer coefficients (P^5 coeff = %d):" % ints[idx55])
for (i, j), c in zip(mons, ints):
    if c:
        print(f"  Y^{i} P^{j}: {c}")

# exact verification to order N: sum c_ij Y^i P^j == 0
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
    mm = (phi_num(q, 50)/phi_num(q**11, 50))**2
    pp = p_prod_chi(11, q, 50)
    YY = mm**2
    PP = pp/32
    tot = mp.mpf(0)
    scale = mp.mpf(0)
    for (i, j), c in zip(mons, ints):
        t = c * YY**i * PP**j
        tot += t
        scale = max(scale, abs(t))
    print(f"  q={qstr}: |F| = {mp.nstr(abs(tot), 5)}, max term = {mp.nstr(scale, 5)}, "
          f"rel = {mp.nstr(abs(tot)/scale, 5)}")
# flagship point
q1 = mp.e**(-mp.pi/mp.sqrt(11))
mm = (phi_num(q1, 50)/phi_num(q1**11, 50))**2
pp = p_prod_chi(11, q1, 50)
YY = mm**2; PP = pp/32
tot = sum(c * YY**i * PP**j for (i, j), c in zip(mons, ints))
print(f"  flagship q1=e^(-pi/sqrt11): m={mp.nstr(mm,20)} (sqrt11={mp.nstr(mp.sqrt(11),20)}), "
      f"p={mp.nstr(pp,20)}, |F| = {mp.nstr(abs(tot), 5)}")

# print the polynomial nicely in (Y, P) and in (m, p)
Ysym, Psym = sp.symbols('Y P')
G = sum(c * Ysym**i * Psym**j for (i, j), c in zip(mons, ints) if c)
print("G(Y,P) =", sp.expand(G))
print("factor G(Y,P):", sp.factor(G))
msym, psym = sp.symbols('m p')
F = sp.expand(G.subs({Ysym: msym**2, Psym: psym/32}))
print("F(m,p) = G(m^2, p/32), expanded, #terms =", len(F.args))
print("F(m,p) =", F)
