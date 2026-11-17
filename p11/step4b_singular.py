"""NOTE (2026-11-14 re-verification): SUPERSEDED. This script reports
"singular points found: 0" and genus 6; its root-finding misses the nodes.
The correct genus analysis is step5_final.py (4 ordinary nodes, genus 2),
confirmed by the discriminant factorisation in paper/REVERIFICATION_LOG.md.
The w44 check at the end of this script is still valid. Kept for audit trail.
"""
"""p11 step 4b: numeric singular-point analysis of G(Y,P)=0, and w44 check.
  - Y-roots of the resultant-gcd quartic via numpy, refined by mpmath.
  - For each, P-roots of G(Y,.) and singular-point test (G_Y = G_P = 0).
  - Multiplicity + tangent-cone analysis at each singular point (mpmath
    derivatives), delta for ordinary points; genus = 6 - sum(delta).
  - Check w44 preserves the curve: G(121/Y, P') = 0 with P' = P(w44 tau).
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
import numpy as np
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
fG = sp.lambdify((Y, P), G, 'mpmath')
fGY = sp.lambdify((Y, P), GY, 'mpmath')
fGP = sp.lambdify((Y, P), GP, 'mpmath')

mp.mp.dps = 50
# quartic from step 4: 59049 Y^4 + 280422 Y^3 + 293768 Y^2 - 64614 Y - 14641
qcoef = np.array([59049, 280422, 293768, -64614, -14641], dtype=float)
yroots = np.roots(qcoef)
print("quartic Y-roots (numpy):", yroots)
sing = []
for yr in yroots:
    if abs(yr.imag) > 1e-6:
        print(f"  Y={yr}: complex root, skipping affine-real analysis "
              f"(singular point would be non-real)")
        continue
    y0 = mp.mpf(str(yr.real))
    # refine Y-root with mpmath Newton on the quartic
    f = lambda t: 59049*t**4 + 280422*t**3 + 293768*t**2 - 64614*t - 14641
    y0 = mp.findroot(f, y0)
    # P-roots of G(y0, P): coefficients in P
    poly = [sum(mp.mpf(coeffs.get((i, j), 0))*y0**i for i in range(6)) for j in range(6)]
    proots = mp.polyroots(poly[::-1], maxsteps=500, extraprec=80)
    for pr in proots:
        if abs(mp.im(pr)) > mp.mpf('1e-30'):
            continue
        p0 = mp.re(pr)
        g, gy, gp = fG(y0, p0), fGY(y0, p0), fGP(y0, p0)
        print(f"  Y={mp.nstr(y0,14)} P={mp.nstr(p0,14)}: |G|={mp.nstr(abs(g),3)}"
              f" |G_Y|={mp.nstr(abs(gy),3)} |G_P|={mp.nstr(abs(gp),3)}")
        if abs(g) < mp.mpf('1e-30') and abs(gy) < mp.mpf('1e-30') and abs(gp) < mp.mpf('1e-30'):
            sing.append((y0, p0))
print("singular points found:", len(sing))

# multiplicity + tangent cone at each singular point
def taylor_deg(y0, p0, maxdeg=6):
    """Lowest total degree with a nonzero partial derivative (mpmath)."""
    for d in range(1, maxdeg+1):
        terms = []
        for a in range(d+1):
            b = d-a
            der = sp.diff(G, Y, a, P, b)
            fd = sp.lambdify((Y, P), der, 'mpmath')
            v = fd(y0, p0)
            terms.append((a, b, v))
        nz = [(a, b, v) for (a, b, v) in terms if abs(v) > mp.mpf('1e-25')]
        if nz:
            return d, terms
    return None, None
total_delta = mp.mpf(0)
for (y0, p0) in sing:
    d, terms = taylor_deg(y0, p0)
    print(f"singular point Y={mp.nstr(y0,12)}, P={mp.nstr(p0,12)}: multiplicity {d}")
    # tangent cone: homogeneous part of degree d as polynomial in t = z/x
    # (coeffs of x^a z^b, a+b=d) -> univariate in s=z/x: sum c_a s^{d-a}
    if d is not None:
        cone = [0]*(d+1)
        for (a, b, v) in terms:
            cone[d-a] += float(v)  # coefficient of s^b
        troots = np.roots(np.array(cone[::-1], dtype=float)) if d > 0 else []
        print(f"   tangent-cone roots (s = dP/dY): {troots}")
        distinct = len(set(np.round(troots, 6)))
        print(f"   distinct tangents: {distinct} of {d}")
        if distinct == d:
            delta = d*(d-1)//2
            print(f"   ordinary => delta = {delta}")
            total_delta += delta
        else:
            print("   NON-ordinary tangent cone: delta needs blow-up analysis")
print("sum of delta (ordinary points only) =", total_delta)
print("=> genus = 6 - sum(delta) =", 6 - total_delta,
      "(if all singular points are ordinary and none were missed)")

# w44 check: G(121/Y, P') = 0 ?
print("w44 check:")
def YP_of_u(u):
    q = mp.e**(-mp.pi*u)
    m = (phi_num(q, 40)/phi_num(q**11, 40))**2
    p = p_prod_chi(11, q, 40)
    return m**2, p/32
for u in ("0.7", "1.3", "0.9"):
    u = mp.mpf(u)
    Y0, P0 = YP_of_u(u)
    Y1, P1 = YP_of_u(1/(11*u))
    print(f"  u={u}: Y'={mp.nstr(Y1,12)} vs 121/Y={mp.nstr(121/Y0,12)};"
          f" G(121/Y, P') = {mp.nstr(fG(121/Y0, P1), 3)};"
          f" G(Y,P) = {mp.nstr(fG(Y0, P0), 3)}")
    print(f"     P'={mp.nstr(P1,12)}, 1-P={mp.nstr(1-P0,12)}")
