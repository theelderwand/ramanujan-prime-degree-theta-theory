"""p17 step 4c: affine singular points of F(m,P)=0 and the genus.

Singular m-coordinates are the roots of the degree-12 factor
  S(m) = 15625 m^12 + 57350 m^11 + 435636 m^10 + 4017474 m^9
         - 7136763 m^8 + 58172316 m^7 - 2075760 m^6 + 39725924 m^5
         + 142214563 m^4 + 63554366 m^3 + 7265596 m^2 - 57222 m - 4913
(the gcd of res_P(F, F_m) and res_P(F, F_P), up to a constant).
For each root m0, the singular P0 is the common root of F(m0, P) and
F_P(m0, P); the tangent-cone discriminant F_mm*F_PP - F_mP^2 decides
node (delta=1) vs. worse singularity.

Genus: g = 55 - delta_inf - N  (plane degree 12, delta_inf = 40 from
step4b).  Cross-check via Riemann-Hurwitz on m (deg 4): g = 3 + T where
T = 12 - N is the number of triple-ramification points.
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
import sympy as sp
import mpmath as mp

mp.mp.dps = 80
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
Fm = sp.diff(F, msym); FP = sp.diff(F, Psym)
Fmm = sp.diff(Fm, msym); FPP = sp.diff(FP, Psym); FmP = sp.diff(Fm, Psym)

S = (15625*msym**12 + 57350*msym**11 + 435636*msym**10 + 4017474*msym**9
     - 7136763*msym**8 + 58172316*msym**7 - 2075760*msym**6 + 39725924*msym**5
     + 142214563*msym**4 + 63554366*msym**3 + 7265596*msym**2 - 57222*msym - 4913)
scoeffs = [mp.mpf(int(c)) for c in sp.Poly(S, msym).all_coeffs()]
mroots = mp.polyroots(scoeffs, maxsteps=200, extraprec=50)
print("number of roots of S:", len(mroots))

def ev(expr, m0, P0):
    return expr.subs({msym: m0, Psym: P0})

nodes = 0
triple = 0
for k, m0 in enumerate(mroots):
    # univariate in P at m0: P0 is a root of F_P(m0, .) that also zeroes F
    gP = sp.Poly(FP.subs({msym: m0}), Psym)
    gc = [mp.mpc(complex(c)) for c in gP.all_coeffs()]
    try:
        gr = mp.polyroots(gc, maxsteps=500, extraprec=80, error=False)
    except Exception as e:
        print(f"root {k}: polyroots failed on F_P: {e}")
        continue
    best = None
    for b in gr:
        Fv = complex(ev(F, m0, b))
        if best is None or abs(Fv) < best[0]:
            best = (abs(Fv), b)
    _, P0 = best
    d = abs(complex(ev(F, m0, P0)))
    # verify the point is singular
    Fv = complex(ev(F, m0, P0)); Fmv = complex(ev(Fm, m0, P0)); FPv = complex(ev(FP, m0, P0))
    # tangent-cone discriminant
    D = complex(ev(Fmm, m0, P0))*complex(ev(FPP, m0, P0)) - complex(ev(FmP, m0, P0))**2
    is_node = abs(D) > mp.mpf('1e-20')
    if is_node:
        nodes += 1
        kind = "NODE delta=1"
    else:
        triple += 1
        kind = "NOT a node (tangent cone degenerate) -- check"
    print(f"root {k}: m0={mp.nstr(m0, 12)} P0={mp.nstr(P0, 12)} "
          f"|F|={abs(Fv):.2e} |Fm|={abs(Fmv):.2e} |FP|={abs(FPv):.2e} "
          f"|D|={abs(D):.3e} {kind}")

print(f"\nnodes N = {nodes}, non-nodes T = {triple}")
print(f"genus g = 55 - 40 - {nodes} = {55-40-nodes}")
print(f"RH cross-check: g = 3 + T = {3+triple}")
