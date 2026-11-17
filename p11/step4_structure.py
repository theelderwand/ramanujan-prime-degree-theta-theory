"""p11 step 4: structure of the curve G(Y,P)=0.
  (a) irreducibility over Q (factor_list).
  (b) invariance under the Fricke-type involution (Y,P) -> (121/Y, 1-P):
      check Y^5 * G(121/Y, 1-P) == +/- G(Y,P) symbolically.
  (c) numeric action of Atkin-Lehner involutions w4, w11, w44 on (Y,P):
      tau = i*u/2 (q = e^{-pi u}); w_Q sends u -> u' with
        w4:  u' = 1/u,   w11: u' = 4/(11u),   w44: u' = 1/(11u).
  (d) genus of the plane quintic via singular points (resultants).
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
import sympy as sp
import mpmath as mp
from theta import phi_num, p_prod_chi

Y, P = sp.symbols('Y P')
coeffs = {
    (0,0): 1, (1,0): -5, (2,0): 10, (3,0): -10, (4,0): 5, (5,0): -1,
    (0,1): 160, (1,1): 10624, (2,1): 19392, (3,1): 2432, (4,1): 160,
    (0,2): 10240, (1,2): -796672, (2,2): -55296, (3,2): -10240,
    (0,3): 327680, (1,3): 6914048, (2,3): 327680,
    (0,4): 5242880, (1,4): -5242880,
    (0,5): 33554432,
}
G = sum(c*Y**i*P**j for (i, j), c in coeffs.items())
print("(a) factor_list(G):", sp.factor_list(G))
print("    G(121, 1) =", G.subs({Y: 121, P: 1}), " (cusp 0: Y=m^2=121, P=p/32=1)")

# (b) involution
Ginv = sp.expand(Y**5 * G.subs({Y: 121/Y, P: 1-P}))
print("(b) Y^5*G(121/Y, 1-P) - G(Y,P) =", sp.expand(Ginv - G))
print("    Y^5*G(121/Y, 1-P) + G(Y,P) =", sp.expand(Ginv + G))

# (c) numeric Atkin-Lehner action
mp.mp.dps = 40
def YP_of_u(u):
    q = mp.e**(-mp.pi*u)
    m = (phi_num(q, 40)/phi_num(q**11, 40))**2
    p = p_prod_chi(11, q, 40)
    return m**2, p/32
for u in ("0.7", "1.3"):
    u = mp.mpf(u)
    Y0, P0 = YP_of_u(u)
    print(f"(c) u={u}: Y={mp.nstr(Y0,15)}, P={mp.nstr(P0,15)}")
    for name, up in (("w4", 1/u), ("w11", 4/(11*u)), ("w44", 1/(11*u))):
        Y1, P1 = YP_of_u(up)
        print(f"    {name} (u'={mp.nstr(up,8)}): Y'={mp.nstr(Y1,15)}, P'={mp.nstr(P1,15)}"
              f" | Y'/Y={mp.nstr(Y1/Y0,6)}, 121/Y={mp.nstr(121/Y0,10)},"
              f" 1-P={mp.nstr(1-P0,10)}, P'/P={mp.nstr(P1/P0,6) if P0 else 'n/a'}")

# (d) singular points: resultants
print("(d) computing resultants (may take a while)...")
GY = sp.diff(G, Y)
GP = sp.diff(G, P)
R1 = sp.resultant(G, GP, P)   # in Y; vanishes at Y-coords where G, G_P share a P-root
R2 = sp.resultant(G, GY, P)
g = sp.gcd(R1, R2)
print("    deg R1 =", sp.Poly(R1, Y).degree(), " deg R2 =", sp.Poly(R2, Y).degree(),
      " deg gcd =", sp.Poly(g, Y).degree())
print("    gcd factorisation over Q:", sp.factor_list(g))
# numeric roots of gcd
roots = sp.nroots(g, n=60, maxsteps=300)
print("    numeric Y-roots of gcd:", [complex(r) for r in roots])
for r in roots:
    yv = complex(r)
    if abs(yv.imag) > 1e-6:
        continue
    yv = mp.mpf(str(yv.real))
    # solve G(yv, P) = 0 for P numerically: use mp.polyroots on coeffs in P
    cp = [mp.mpf(coeffs.get((i, j), 0))*yv**i for j in range(6) for i in range(6)]
    # polynomial in P: sum_j (sum_i c_ij y^i) P^j
    poly = [sum(mp.mpf(coeffs.get((i, j), 0))*yv**i for i in range(6)) for j in range(6)]
    try:
        proots = mp.polyroots(poly[::-1], maxsteps=200, extraprec=50)
    except Exception as e:
        print(f"    Y={mp.nstr(yv,12)}: polyroots failed: {e}")
        continue
    for pr in proots:
        pv = complex(pr)
        if abs(pv.imag) > 1e-6:
            continue
        pv = mp.mpf(str(pv.real))
        # evaluate G, G_Y, G_P at (yv, pv) with high precision via sympy subs (exact-ish)
        vals = [float(sp.N(G.subs({Y: sp.nsimplify(yv, rational=False),
                                    P: sp.nsimplify(pv, rational=False)}), 20)) for _ in [0]]
        # simpler: evaluate with mpmath lambdified
        fG = sp.lambdify((Y, P), G, 'mpmath')
        fGY = sp.lambdify((Y, P), GY, 'mpmath')
        fGP = sp.lambdify((Y, P), GP, 'mpmath')
        print(f"    candidate singular point Y={mp.nstr(yv,14)}, P={mp.nstr(pv,14)}:"
              f" G={mp.nstr(fG(yv,pv),3)}, G_Y={mp.nstr(fGY(yv,pv),3)},"
              f" G_P={mp.nstr(fGP(yv,pv),3)}")
