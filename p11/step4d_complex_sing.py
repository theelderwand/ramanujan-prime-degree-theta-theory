"""p11 step 4d (v2): check ALL P-roots (incl. complex) at each quartic Y-root
for singular points, using numpy (float64 is ample: singular points have
|G_Y|,|G_P| ~ 0 vs ~1e6 at smooth points)."""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
import numpy as np

coeffs = {
    (0,0): 1, (1,0): -5, (2,0): 10, (3,0): -10, (4,0): 5, (5,0): -1,
    (0,1): 160, (1,1): 10624, (2,1): 19392, (3,1): 2432, (4,1): 160,
    (0,2): 10240, (1,2): -796672, (2,2): -55296, (3,2): -10240,
    (0,3): 327680, (1,3): 6914048, (2,3): 327680,
    (0,4): 5242880, (1,4): -5242880,
    (0,5): 33554432,
}

def ev(Yv, Pv, dY=0, dP=0):
    """Evaluate partial derivative of G at (Yv, Pv) by explicit sums."""
    s = 0j
    for (i, j), c in coeffs.items():
        if i < dY or j < dP:
            continue
        import math
        f = math.factorial(i)//math.factorial(i-dY) if dY else 1
        g = math.factorial(j)//math.factorial(j-dP) if dP else 1
        s += c*f*g*(Yv**(i-dY))*(Pv**(j-dP))
    return s

qcoef = [59049, 280422, 293768, -64614, -14641]
for yr in np.roots(np.array(qcoef, dtype=float)):
    y0 = yr.real
    poly = [sum(coeffs.get((i, j), 0)*y0**i for i in range(6)) for j in range(6)]
    proots = np.roots(np.array(poly[::-1], dtype=float))
    print(f"Y = {y0:.12f}:")
    for pr in proots:
        g, gy, gp = ev(y0, pr), ev(y0, pr, 1, 0), ev(y0, pr, 0, 1)
        flag = "  <-- SINGULAR" if (abs(gy) < 1e-6 and abs(gp) < 1e-6) else ""
        print(f"    P = {pr.real:+.12f}{pr.imag:+.12f}j: |G|={abs(g):.2e}"
              f" |G_Y|={abs(gy):.2e} |G_P|={abs(gp):.2e}{flag}")
