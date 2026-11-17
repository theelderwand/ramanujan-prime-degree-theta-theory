"""p17 step 4: structure of the curve F(m,P) = 0.

- Irreducibility, pure parts: F(m,0) = (m-1)^12, F(0,P) = (256P-1)^4;
  F = (m-1)^12 + (256P-1)^4 - 1 - 1024*P*m*S(m,P), S explicit (step 3).
- Singular points: resultants res_P(F, F_m), res_P(F, F_P); gcd gives
  m-coordinates of affine singular points.
- Discriminant disc_P(F) for the Riemann-Hurwitz check on m: C -> P^1 (deg 4).
- Point at infinity: Newton-polygon analysis (edge (0,0)-(12,4), edge
  polynomial (w^3-256)^4); delta computed via genus cross-checks.
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
import sympy as sp

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
print("irreducible over Q:", sp.factor(F) == F)

Fm = sp.diff(F, msym)
FP = sp.diff(F, Psym)
print("computing resultants...")
res_m = sp.resultant(F, Fm, Psym)
res_P = sp.resultant(F, FP, Psym)
print("deg res_P(F, F_m) in m:", sp.Poly(res_m, msym).degree())
print("deg res_P(F, F_P) in m:", sp.Poly(res_P, msym).degree())
g = sp.gcd(sp.Poly(res_m, msym), sp.Poly(res_P, msym))
print("gcd degree:", g.degree())
print("gcd factorisation:", sp.factor(g.as_expr()))
print("computing disc_P(F)...")
disc = sp.discriminant(F, Psym)
print("disc_P(F) degree in m:", sp.Poly(disc, msym).degree())
print("disc_P(F) factorisation:")
print("  ", sp.factor(disc))
