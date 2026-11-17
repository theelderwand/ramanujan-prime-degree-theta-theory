"""p13 step 4: structure of the curve F(m,P) = 0.

- Irreducibility (sympy factor, already seen in step 2).
- Pure parts: F(m,0) = -(m-1)^7, F(0,P) = (1+64P)^3;
  F = (1+64P)^3 - (m-1)^7 - 1 + P*m*S(m,P), S explicit.
- Singular points: resultants res_P(F, F_m), res_P(F, F_P); gcd gives
  m-coordinates of affine singular points.
- Discriminant disc_P(F) and Riemann-Hurwitz for the degree-3 map m: C -> P^1.
- Point at infinity: Puiseux analysis (single branch (s^4, s^7), delta = 9).
- Genus via plane-curve formula, cross-checked by Riemann-Hurwitz.
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
import sympy as sp
import mpmath as mp

msym, Psym = sp.symbols('m P')
coeffs = {
    (0,0): 1, (1,0): -7, (2,0): 21, (3,0): -35, (4,0): 35, (5,0): -21,
    (6,0): 7, (7,0): -1,
    (0,1): 192, (1,1): 7424, (2,1): 10368, (3,1): -2816, (4,1): 1216,
    (0,2): 12288, (1,2): -188416, (2,2): 45056,
    (0,3): 262144,
}
F = sum(c*msym**i*Psym**j for (i, j), c in coeffs.items())
print("irreducible over Q:", sp.factor(F) == F)

# structural form
R = sp.expand(F - (1+64*Psym)**3 + (msym-1)**7 + 1)
print("R+1 check: F - (1+64P)^3 + (m-1)^7 + 1 =", sp.factor(R))
S = sp.factor(R/(Psym*msym))
print("S(m,P) = (F - (1+64P)^3 + (m-1)^7 + 1)/(P*m) =", S)

# partials and resultants
Fm = sp.diff(F, msym)
FP = sp.diff(F, Psym)
print("computing resultants (may take a while)...")
res_m = sp.resultant(F, Fm, Psym)
res_P = sp.resultant(F, FP, Psym)
print("deg res_P(F, F_m) in m:", sp.Poly(res_m, msym).degree())
print("deg res_P(F, F_P) in m:", sp.Poly(res_P, msym).degree())
g = sp.gcd(sp.Poly(res_m, msym), sp.Poly(res_P, msym))
print("gcd degree:", g.degree())
print("gcd factorisation:", sp.factor(g.as_expr()))
# singular m-coordinates (roots of gcd)
roots_m = sp.nroots(g.as_expr(), n=60, maxsteps=300)
print("singular m-coordinates (numeric, 60 digits -> shown 12):")
sing = []
for r in roots_m:
    if abs(sp.im(r)) < mp.mpf(10)**(-30):
        m0 = float(sp.re(r))
        # solve F(m0, P) = 0 for P
        Fp = sp.Poly(F.subs(msym, m0), Psym)
        proots = sp.nroots(Fp.as_expr(), n=50, maxsteps=200)
        for pr in proots:
            if abs(sp.im(pr)) < 1e-25:
                p0 = float(sp.re(pr))
                # check singularity: Fm and FP vanish
                vm = float(Fm.subs({msym: m0, Psym: p0}))
                vp = float(FP.subs({msym: m0, Psym: p0}))
                sing.append((m0, p0, vm, vp))
                print(f"  (m,P) = ({m0:.12f}, {p0:.12f})  Fm={vm:.2e} FP={vp:.2e}")
print("affine singular points found:", len(sing))

# discriminant of F as cubic in P
disc = sp.discriminant(F, Psym)
print("disc_P(F) degree in m:", sp.Poly(disc, msym).degree())
print("disc_P(F) factorisation:", sp.factor(disc))
