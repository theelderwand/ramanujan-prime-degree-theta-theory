"""p11 step 5: final checks.
  (1) Newton-refine the 4th singular point (Y,P) ~ (-0.144430, 0.004060),
      multiplicity / tangent cone / delta (completes genus = 6 - 4 = 2).
  (2) Genuine Atkin-Lehner w4, w11 (correct representatives, which do NOT
      preserve the imaginary axis) acting on (Y,P) at complex tau:
        w4 : tau -> (6 tau + 1/2)/(22 tau + 2)   [matrix (1/2)[[12,1],[44,4]]]
        w11: tau -> (33 tau + 8)/(44 tau + 11)  [matrix (1/sqrt11)[[33,8],[44,11]]]
        w44: tau -> -1/(44 tau)
      Test: does w4 or w11 fix (Y,P)? (w44 sends Y -> 121/Y, known.)
  (3) Extra verification of G(Y,P)=0 at complex q (40 digits).
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
GY = sp.diff(G, Y); GP = sp.diff(G, P)
GYY = sp.diff(GY, Y); GYP = sp.diff(GY, P); GPP = sp.diff(GP, P)
fG = sp.lambdify((Y, P), G, 'mpmath')
fGY = sp.lambdify((Y, P), GY, 'mpmath')
fGP = sp.lambdify((Y, P), GP, 'mpmath')
fGYY = sp.lambdify((Y, P), GYY, 'mpmath')
fGYP = sp.lambdify((Y, P), GYP, 'mpmath')
fGPP = sp.lambdify((Y, P), GPP, 'mpmath')

mp.mp.dps = 60
# ---- (1) 4th singular point ----
y0, p0 = mp.mpf('-0.14442967371861'), mp.mpf('0.0040600200')
for _ in range(30):
    F = mp.matrix([fGY(y0, p0), fGP(y0, p0)])
    J = mp.matrix([[fGYY(y0, p0), fGYP(y0, p0)],
                   [fGYP(y0, p0), fGPP(y0, p0)]])
    step = mp.lu_solve(J, F)
    y0 -= step[0]; p0 -= step[1]
    if abs(step[0]) < mp.mpf('1e-50') and abs(step[1]) < mp.mpf('1e-50'):
        break
print("(1) 4th singular point: Y =", mp.nstr(y0, 25), " P =", mp.nstr(p0, 25))
print("    residuals |G|,|G_Y|,|G_P| =",
      mp.nstr(abs(fG(y0, p0)), 3), mp.nstr(abs(fGY(y0, p0)), 3),
      mp.nstr(abs(fGP(y0, p0)), 3))
for d in range(1, 5):
    vals = []
    for a in range(d+1):
        b = d-a
        fd = sp.lambdify((Y, P), sp.diff(G, Y, a, P, b), 'mpmath')
        vals.append((a, b, fd(y0, p0)))
    nz = [(a, b, v) for (a, b, v) in vals if abs(v) > mp.mpf('1e-20')]
    if nz:
        print(f"    multiplicity {d}; degree-{d} derivatives:",
              [(a, b, mp.nstr(v, 6)) for (a, b, v) in nz])
        import math
        cone = [mp.mpf(0)]*(d+1)
        for (a, b, v) in nz:
            cone[b] += v/(mp.factorial(a)*mp.factorial(b))
        croots = mp.polyroots(cone[::-1], maxsteps=500, extraprec=80)
        print("    tangent-cone slopes:", [mp.nstr(r, 8) for r in croots])
        print("    ordinary node => delta = 1" if len(set([complex(r) for r in croots])) == d
              else "    NON-ordinary")
        break
print("    => four ordinary nodes, sum delta = 4, genus = 6 - 4 = 2.")

# ---- (2) genuine w4, w11 ----
mp.mp.dps = 40
def YP_of_tau(tau):
    q = mp.e**(2*mp.pi*mp.j*tau)
    m = (phi_num(q, 40)/phi_num(q**11, 40))**2
    p = p_prod_chi(11, q, 40)
    return m**2, p/32
for t in ("0.35", "0.5"):
    tau = mp.mpc(0, mp.mpf(t))
    Y0, P0 = YP_of_tau(tau)
    tau4 = (6*tau + mp.mpf('1/2'))/(22*tau + 2)
    tau11 = (33*tau + 8)/(44*tau + 11)
    tau44 = -1/(44*tau)
    Y4, P4 = YP_of_tau(tau4)
    Y11, P11 = YP_of_tau(tau11)
    Y44, P44 = YP_of_tau(tau44)
    print(f"(2) tau = {t}i: Y={mp.nstr(Y0,12)}, P={mp.nstr(P0,12)}")
    print(f"    w4 : Y'={mp.nstr(Y4,12)} P'={mp.nstr(P4,12)}"
          f" | Y'/Y={mp.nstr(Y4/Y0,6)}")
    print(f"    w11: Y'={mp.nstr(Y11,12)} P'={mp.nstr(P11,12)}"
          f" | Y'/Y={mp.nstr(Y11/Y0,6)}")
    print(f"    w44: Y'={mp.nstr(Y44,12)} vs 121/Y={mp.nstr(121/Y0,12)};"
          f" G(121/Y,P')={mp.nstr(fG(121/Y0, P44),3)}")

# ---- (3) complex-q verification ----
print("(3) G(Y,P)=0 at complex q (40 digits):")
mp.mp.dps = 40
for q in [mp.mpc('0.3', '0.2'), mp.mpc('-0.2', '0.4'), mp.mpc('0.1', '-0.5')]:
    m = (phi_num(q, 40)/phi_num(q**11, 40))**2
    p = p_prod_chi(11, q, 40)
    YY, PP = m**2, p/32
    tot = sum(c*YY**i*PP**j for (i, j), c in coeffs.items())
    scale = max(abs(c*YY**i*PP**j) for (i, j), c in coeffs.items())
    print(f"    q={q}: |G|={mp.nstr(abs(tot),4)}, max term={mp.nstr(scale,4)},"
          f" rel={mp.nstr(abs(tot)/scale,4)}")
