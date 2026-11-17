"""NOTE (2026-11-14 re-verification): SUPERSEDED. This script finds only 3
nodes and reports genus 3; it misses the 4th node at Y=-0.1444, P=0.00406
(its P-coordinate is a double root, found by step5_final.py). Correct: 4
ordinary nodes, genus = 6 - 4 = 2. Kept for audit trail.
"""
"""p11 step 4c: refine singular points by Newton on (G_Y, G_P) = 0,
then multiplicity / tangent-cone / delta, and genus = 6 - sum(delta).
Thresholds are relative: near a singular point the degree-1 derivatives
are ~1e6-1e7 at generic nearby points, so |G_Y|,|G_P| < 1e-12 signals one.
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
import numpy as np
import sympy as sp
import mpmath as mp

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
fGY = sp.lambdify((Y, P), GY, 'mpmath')
fGP = sp.lambdify((Y, P), GP, 'mpmath')
fGYY = sp.lambdify((Y, P), GYY, 'mpmath')
fGYP = sp.lambdify((Y, P), GYP, 'mpmath')
fGPP = sp.lambdify((Y, P), GPP, 'mpmath')
fG = sp.lambdify((Y, P), G, 'mpmath')

mp.mp.dps = 60
starts = [(-2.9398849906605, 0.26961640660751),
          (-1.9622462430645, -0.011516281126006),
          (0.29758971402802, 0.041980800954304)]
total_delta = 0
for (y00, p00) in starts:
    y0, p0 = mp.mpf(str(y00)), mp.mpf(str(p00))
    # Newton on (GY, GP) = 0
    for _ in range(30):
        F = mp.matrix([fGY(y0, p0), fGP(y0, p0)])
        J = mp.matrix([[fGYY(y0, p0), fGYP(y0, p0)],
                       [fGYP(y0, p0), fGPP(y0, p0)]])
        try:
            step = mp.lu_solve(J, F)
        except Exception as e:
            print("  Newton solve failed:", e)
            break
        y0 -= step[0]; p0 -= step[1]
        if abs(step[0]) < mp.mpf('1e-50') and abs(step[1]) < mp.mpf('1e-50'):
            break
    print(f"singular point: Y={mp.nstr(y0, 30)}, P={mp.nstr(p0, 30)}")
    print(f"   residuals: |G|={mp.nstr(abs(fG(y0,p0)),3)} |G_Y|={mp.nstr(abs(fGY(y0,p0)),3)}"
          f" |G_P|={mp.nstr(abs(fGP(y0,p0)),3)}")
    # multiplicity: smallest degree d with a nonzero derivative
    mult = None
    for d in range(1, 6):
        vals = []
        for a in range(d+1):
            b = d-a
            fd = sp.lambdify((Y, P), sp.diff(G, Y, a, P, b), 'mpmath')
            vals.append((a, b, fd(y0, p0)))
        nz = [(a, b, v) for (a, b, v) in vals if abs(v) > mp.mpf('1e-20')]
        if nz:
            mult = d
            print(f"   multiplicity {d}; degree-{d} derivative values:")
            for (a, b, v) in nz:
                print(f"      d^{d}/dY^{a}dP^{b} = {mp.nstr(v, 8)}")
            # tangent cone: homogeneous degree-d part, as polynomial in s = dP/dY
            # term x^a z^b -> coefficient of s^b is the derivative / (a! b!)
            import math
            cone = [mp.mpf(0)]*(d+1)
            for (a, b, v) in nz:
                cone[b] += v/(mp.factorial(a)*mp.factorial(b))
            # roots of sum_b cone[b] s^b
            croots = mp.polyroots(cone[::-1], maxsteps=500, extraprec=80)
            print(f"   tangent-cone slopes s=dP/dY: {[mp.nstr(r, 8) for r in croots]}")
            # distinct slopes (cluster within 1e-6)
            slopes = sorted([complex(r) for r in croots], key=lambda z: (z.real, z.imag))
            distinct = 0
            for s in slopes:
                if all(abs(s-t) > 1e-6 for t in slopes[:slopes.index(s)]):
                    distinct += 1
            print(f"   distinct tangents: {distinct} of {d}")
            if distinct == d:
                delta = d*(d-1)//2
                print(f"   ORDINARY => delta = {delta}")
                total_delta += delta
            else:
                print("   NON-ordinary: delta not computed by this formula")
            break
    if mult is None:
        print("   multiplicity > 5?? (unexpected)")
print("sum delta =", total_delta, "=> genus = 6 -", total_delta, "=", 6-total_delta)
