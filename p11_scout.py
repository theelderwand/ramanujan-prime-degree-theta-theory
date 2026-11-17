"""Phase-2 scouting for p = 11.

Questions (report §3.2, §7 Phase 2 Steps 1-2):
  (a) At the flagship point q1 = e^{-pi/sqrt 11}: compute u_k (k=1..5),
      w_k = u_k^{11}/p_prod. Are the w_k Galois conjugates in the real
      cyclotomic field Q(zeta_11)^+ = Q(2 cos(2pi/11))? Test: minimal
      polynomial of each w_k via PSLQ (degree should be 5, same polynomial
      for all five if they are conjugates).
  (b) At generic q: are the power sums s_j = sum_k w_k^j rational functions
      of p = p_prod and m = phi^2(q)/phi^2(q^11) of low degree? Fit with a
      small monomial basis over Q from many data points.
  (c) Values needed later: m at the flagship point; G_121 minimal polynomial
      guess (PSLQ at 80 digits).

All of this is EXPLORATORY (scouting, not proof). Results are reported
honestly, including failures.
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
import mpmath as mp
from theta import phi_num, u_k_num, p_prod_num, p_prod_chi, G_num

DPS = 50
mp.mp.dps = DPS


def minpoly_pslq(x, deg=5, dps=50):
    """Try to find the minimal polynomial of x (degree <= deg) via PSLQ.
    Returns list of integer coeffs [a0..ad] with sum a_i x^i = 0, or None."""
    mp.mp.dps = dps
    powers = [x ** i for i in range(deg + 1)]
    rel = mp.pslq(powers, maxcoeff=10 ** 6, maxsteps=10000)
    return rel


def main():
    print("=== (a) flagship point q1 = e^{-pi/sqrt 11} ===")
    q1 = mp.e ** (-mp.pi / mp.sqrt(11))
    p = p_prod_chi(11, q1, DPS)
    print(f"  p = {mp.nstr(p, 30)}  (2/G_11^4 = {mp.nstr(2 / G_num(11, DPS) ** 4, 30)})")
    m = (phi_num(q1, DPS) / phi_num(q1 ** 11, DPS)) ** 2
    print(f"  m = phi^2(q)/phi^2(q^11) = {mp.nstr(m, 30)}")
    us = [u_k_num(11, k, q1, DPS) for k in range(1, 6)]
    print(f"  u_k = {[mp.nstr(u, 20) for u in us]}")
    print(f"  decreasing u1>...>u5: {all(us[i] > us[i+1] for i in range(4))}")
    ws = [u ** 11 / p for u in us]
    print(f"  w_k = u_k^11/p = {[mp.nstr(w, 20) for w in ws]}")
    print(f"  prod w_k = {mp.nstr(mp.fprod(ws), 20)}  (p^6 = {mp.nstr(p**6, 20)})")
    print()
    print("  minimal polynomials of w_k (PSLQ, degree <= 5, 50 digits):")
    polys = []
    for k, w in enumerate(ws, 1):
        rel = minpoly_pslq(w, 5, DPS)
        polys.append(rel)
        print(f"    w_{k}: {rel}")
    if all(r is not None for r in polys):
        # normalise sign and compare
        def norm(r):
            r = list(r)
            if r[-1] < 0:
                r = [-x for x in r]
            return r
        same = all(norm(r) == norm(polys[0]) for r in polys)
        print(f"  all five w_k share one minimal polynomial: {same}")
        if same:
            print(f"  common minpoly coeffs (const..x^5): {norm(polys[0])}")
    print()
    print("  sanity: minimal polynomial of 2cos(2pi/11) (should be degree 5):")
    c = 2 * mp.cos(2 * mp.pi / 11)
    print(f"    2cos(2pi/11) = {mp.nstr(c, 25)}; PSLQ: {minpoly_pslq(c, 5, DPS)}")
    print(f"    known: x^5 + x^4 - 4x^3 - 3x^2 + 3x + 1 = 0")

    print()
    print("=== (b) power sums s_j = sum w_k^j at generic q ===")
    import sympy as sp
    # gather data at several q values
    qs = [mp.mpf(x) for x in ("0.2", "0.3", "0.4", "0.5", "0.6")]
    data = []
    for q in qs:
        p_q = p_prod_chi(11, q, DPS)
        m_q = (phi_num(q, DPS) / phi_num(q ** 11, DPS)) ** 2
        us_q = [u_k_num(11, k, q, DPS) for k in range(1, 6)]
        ws_q = [u ** 11 / p_q for u in us_q]
        ss = [sum(w ** j for w in ws_q) for j in range(1, 6)]
        data.append((q, p_q, m_q, ss))
        print(f"  q={q}: p={mp.nstr(p_q, 12)}, m={mp.nstr(m_q, 12)}")
        print(f"     s_1..s_5 = {[mp.nstr(s, 12) for s in ss]}")
    print()
    print("  Fit attempt: s_j as a polynomial in (p, m) with rational coeffs.")
    print("  Monomial basis: 1, p, p^2, p^3, p^4, m, m^2, p*m, p^2*m, p*m^2.")
    print("  (5 data points, 10 unknowns -> underdetermined; use exact rational")
    print("   solve on the first 10 points would be needed. We have only 5 q values,")
    print("   so instead test the SIMPLEST hypotheses first.)")
    # Hypothesis tests with the 5 points: is s_1 a polynomial in p alone? in m alone?
    # Use sympy rational interpolation: fit s_1 = a0 + a1 p + a2 p^2 + a3 p^3 + a4 p^4
    # (5 points, 5 unknowns) then check residual is ~0 (it will be exact by construction
    # for 5 points; the real test is whether the coeffs are small rationals).
    for j in range(5):
        xs = [float(d[1]) for d in data]   # p values
        ys = [float(d[3][j]) for d in data]  # s_j
        # Vandermonde solve in rationals via sympy with high-precision inputs
        # use mpmath values converted to sympy Rational via nsimplify on 50-digit data
        import sympy as sp
        A = sp.Matrix([[sp.nsimplify(mp.mpf(d[1]) ** a, rational=False) for a in range(5)]
                       for d in data])
        b = sp.Matrix([sp.nsimplify(mp.mpf(d[3][j]), rational=False) for d in data])
        try:
            sol = A.LUsolve(b)
            print(f"  s_{j+1} = sum_{{a=0..4}} c_a p^a with c = {[str(sp.nsimplify(x, [sp.sqrt(2)])) for x in sol]}")
        except Exception as e:
            print(f"  s_{j+1}: fit failed: {e}")
    print()
    print("  NOTE: a 5-point quartic fit is vacuous (interpolation). The honest test")
    print("  needs >10 points and a check that the fitted coefficients are SMALL")
    print("  rationals stable across point sets. See PHASE1_REPORT.md.")

    print()
    print("=== (c) G_121 minimal polynomial guess (80 digits) ===")
    mp.mp.dps = 80
    g121 = G_num(121, 80)
    print(f"  G_121 = {mp.nstr(g121, 60)}")
    for deg in (2, 3, 4, 5, 6, 8):
        rel = minpoly_pslq(g121, deg, 80)
        print(f"  degree <={deg}: PSLQ relation: {rel}")
    # also G_121^? : try G_121^4, G_121^8 etc. (G_n often has nicer powers)
    for pw in (2, 4):
        rel = minpoly_pslq(g121 ** pw, 4, 80)
        print(f"  G_121^{pw}: degree<=4 PSLQ: {rel}")


if __name__ == "__main__":
    main()
