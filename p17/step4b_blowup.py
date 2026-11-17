"""p17 step 4b: delta-invariant at infinity via iterated blow-ups.

The point at infinity [0:1:0] of F(m,P)=0, in the chart P=1 with
coordinates (u, z) = (m/P, 1/P), has multiplicity 8 (tangent cone z^8).
We resolve it by iterated blow-ups, summing m_Q(m_Q-1)/2 over all
infinitely near points Q (including the original point).  Exact sympy
arithmetic; exceptional points located via univariate solving (roots may
be algebraic -- tracked numerically at 50 digits for the recursion, with
multiplicities computed exactly from the translated polynomial).

Also: affine singular points (from the degree-12 gcd factor) are located
the same way (their m-coordinates are its roots); nodes contribute delta=1.
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
import sympy as sp
import numpy as np
import mpmath as mp

mp.mp.dps = 50
x, y = sp.symbols('x y')
msym, Psym = sp.symbols('m P')
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
F = sum(c*msym**i*Psym**j for (i, j), c in coeffs.items())

def multiplicity(poly, pt=(0, 0)):
    """Lowest total degree of poly translated so pt is the origin."""
    p = sp.expand(poly.subs({x: x+pt[0], y: y+pt[1]}))
    p = sp.Poly(p, x, y)
    if p is None or p.is_zero:
        return None
    return min(sum(m) for m in p.monoms())

def blowup_strict(poly, chart):
    """Strict transform of poly (singular at origin, multiplicity m) under
    blow-up.  chart 'A': x = x1*y (new coords (x1, y)); chart 'B': y = y1*x
    (new coords (x, y1)).  Returns the strict transform polynomial in (x, y)
    (reusing symbols: in chart A the coords are (x1, y) -> rename to (x, y);
    in chart B (x, y1) -> (x, y))."""
    m = multiplicity(poly)
    if chart == 'A':
        sub = poly.subs({x: x*y})
        p = sp.Poly(sp.expand(sub), x, y)
        # divide by y^m
        q, r = sp.div(sp.Poly(p.as_expr(), y), y**m)
        assert r.is_zero
        return sp.expand(q.as_expr())
    else:
        sub = poly.subs({y: y*x})
        p = sp.Poly(sp.expand(sub), x, y)
        q, r = sp.div(sp.Poly(p.as_expr(), x), x**m)
        assert r.is_zero
        return sp.expand(q.as_expr())

def exceptional_points(poly, chart):
    """Points of the strict transform on the exceptional divisor.
    Chart A: y=0, solve poly(x, 0)=0 for x.  Chart B: x=0, solve poly(0, y)=0.
    Returns list of exact sympy roots (Rational or RootOf)."""
    if chart == 'A':
        uni = sp.Poly(poly.subs({y: 0}), x)
    else:
        uni = sp.Poly(poly.subs({x: 0}), y)
    if uni is None or uni.is_zero:
        return []
    if uni.degree() < 1:
        return []
    roots = sp.roots(uni.as_expr())  # dict root -> multiplicity; distinct roots
    return list(roots.keys())

def resolve(poly, pt=(0, 0), depth=0, max_depth=12, label="root"):
    """Recursively blow up; returns (total delta contribution, tree log).
    Each DISTINCT infinitely near point is visited once: roots are
    deduplicated (sp.roots as a dict), and a point found in both charts
    (chart-A coord q <-> chart-B coord 1/q) is counted once."""
    p = sp.expand(poly.subs({x: x+pt[0], y: y+pt[1]}))
    m = multiplicity(p)
    if m is None:
        return 0, []
    delta = m*(m-1)//2
    log = [(label, depth, pt, m, delta)]
    if m <= 1 or depth >= max_depth:
        return delta, log
    total = delta
    seen = set()  # canonical representatives of points already visited
    for chart in ('A', 'B'):
        st = blowup_strict(p, chart)
        pts = exceptional_points(st, chart)
        for q in pts:
            if not q.is_real:
                log.append((f"{label}/{chart} complex pt {q}", depth+1, q, None, None))
                continue
            # canonical representative: chart-A coordinate of the point
            canon = q if chart == 'A' else (1/q if q != 0 else sp.oo)
            if canon in seen:
                continue
            seen.add(canon)
            ptq = (q, 0) if chart == 'A' else (0, q)
            mq = multiplicity(st, ptq)
            if mq is not None and mq >= 1:
                d, l = resolve(st, ptq, depth+1, max_depth,
                               f"{label}/{chart}@{sp.N(q, 8)}")
                total += d
                log.extend(l)
    return total, log

# --- point at infinity: chart P=1, (u, z) = (m/P, 1/P) ---
# F_tilde(u, z) = z^12 * F(u/z, 1/z)
u, z = sp.symbols('u z')
Ftilde = 0
for (i, j), c in coeffs.items():
    Ftilde += c * u**i * z**(12-i-j)
Ftilde = sp.expand(Ftilde)
print("Ftilde lowest-degree terms check: mult at origin =",
      min(sum(m) for m in sp.Poly(Ftilde, u, z).monoms()))
# rename (u, z) -> (x, y) for the resolver
Fxy = sp.expand(Ftilde.subs({u: x, z: y}))
delta_inf, log = resolve(Fxy, (0, 0), label="inf")
print("delta at infinity =", delta_inf)
for entry in log:
    print("  ", entry)
