"""Phase-2 scouting, part 2: power-sum relations via PSLQ.

Key idea (calibrated on p=7 first): let w_k = u_k^p / p_prod and
s_j = sum_k w_k^j.  The w_k depend on the cyclic order of the auxiliary
roots, so s_j need not be a polynomial in (p, m) -- but it should satisfy
a LOW-DEGREE polynomial relation over Q(p, m) (for p=7, s_1 satisfies a
quadratic over Q(p, X), X = phi^4(q)/phi^4(q^7), because the two cyclic
orders give the two roots of that quadratic).

Method: at ONE generic q, run PSLQ on the vector of monomials
{p^a m^b s_1^c : a+b <= D, c <= 2}.  A small-coefficient relation found
there is then VERIFIED at several other q values.  This separates a real
relation from interpolation noise.

Also: minimal polynomial of m = phi^2(q)/phi^2(q^11) at the flagship point.
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
import mpmath as mp
from theta import phi_num, u_k_num, p_prod_chi

DPS = 60
mp.mp.dps = DPS


def s1_and_pm(n, q):
    p = p_prod_chi(n, q, DPS)
    m = (phi_num(q, DPS) / phi_num(q ** n, DPS)) ** 2
    us = [u_k_num(n, k, q, DPS) for k in range(1, (n - 1) // 2 + 1)]
    ws = [u ** n / p for u in us]
    s1 = sum(ws)
    return p, m, s1, ws


def pslq_relation(n, q, D, cdeg=2):
    """Find integer relation among {p^a m^b s1^c} at one q. Returns
    (monomials, coeffs) or (None, None)."""
    p, m, s1, _ = s1_and_pm(n, q)
    mons = []
    labels = []
    for c in range(cdeg + 1):
        for a in range(D + 1):
            for b in range(D + 1 - a):
                mons.append(p ** a * m ** b * s1 ** c)
                labels.append(f"p^{a} m^{b} s1^{c}")
    rel = mp.pslq(mons, maxcoeff=10 ** 8, maxsteps=20000)
    return labels, rel, (p, m, s1)


def verify_relation(n, labels, rel, qs):
    """Check the relation sum rel_i * mon_i = 0 at other q values."""
    ok = True
    for q in qs:
        q = mp.mpf(q)
        p, m, s1, _ = s1_and_pm(n, q)
        vals = {}
        for lab in labels:
            # parse "p^a m^b s1^c"
            parts = lab.split()
            v = mp.mpf(1)
            for part in parts:
                base, ex = part.split("^")
                v *= {"p": p, "m": m, "s1": s1}[base] ** int(ex)
            vals[lab] = v
        tot = sum(c * vals[lab] for c, lab in zip(rel, labels))
        scale = max(abs(vals[lab]) for lab in labels)
        print(f"    q={q}: |sum| = {mp.nstr(abs(tot), 5)}  (scale {mp.nstr(scale, 5)})"
              f"  rel err = {mp.nstr(abs(tot) / scale, 5)}")
        if abs(tot) / scale > mp.mpf(10) ** (-30):
            ok = False
    return ok


def main():
    print("=== calibration: p=7, s_1 over Q(p, X) ===")
    print("  (X = phi^4(q)/phi^4(q^7); here m as computed IS X^2, so we use X.)")
    # for p=7 use X = phi^4(q)/phi^4(q^7) directly
    def s1_X7(q):
        p = p_prod_chi(7, q, DPS)
        X = (phi_num(q, DPS) / phi_num(q ** 7, DPS)) ** 4
        us = [u_k_num(7, k, q, DPS) for k in range(1, 4)]
        ws = [u ** 7 / p for u in us]
        return p, X, sum(ws)
    q0 = mp.mpf("0.4")
    p0, X0, s0 = s1_X7(q0)
    print(f"  at q=0.4: p={mp.nstr(p0, 15)}, X={mp.nstr(X0, 15)}, s1={mp.nstr(s0, 15)}")
    # NOTE: the first attempt (X-degree <= 3) found only the modular equation
    # itself, factored as s1 * Q(p,X) = 0.  To force a genuine s1-relation we
    # reduce X^2 via the modular equation X^2 = (2+5p)X - (1-p)^3, so X only
    # appears to power 0 or 1, and Q cannot hide in the monomial set.
    mons = []
    labels = []
    for c in range(3):
        for a in range(5):
            for b in range(2):
                mons.append(p0 ** a * X0 ** b * s0 ** c)
                labels.append(f"p^{a} X^{b} s1^{c}")
    rel = mp.pslq(mons, maxcoeff=10 ** 8, maxsteps=20000)
    print(f"  PSLQ relation (X-degree <= 1 via modular eq, {len(mons)} monomials):")
    if rel is None:
        print("    None found")
    else:
        terms = [(int(c), l) for c, l in zip(rel, labels) if c != 0]
        for c, l in terms:
            print(f"    {c:+d} * {l}")
        # genuineness check: does the relation involve s1^2, or is it
        # s1 * (multiple of the modular equation)?  If only s1^1 terms occur,
        # test whether the coefficient polynomial is a multiple of Q.
        import sympy as sp
        p_s, X_s = sp.symbols('p X')
        Q = X_s**2 - (2 + 5*p_s)*X_s + (1 - p_s)**3
        B = sum(int(c) * p_s**int(l.split()[0].split('^')[1])
                * X_s**int(l.split()[1].split('^')[1])
                for c, l in terms if l.endswith('s1^1'))
        A = sum(int(c) * p_s**int(l.split()[0].split('^')[1])
                * X_s**int(l.split()[1].split('^')[1])
                for c, l in terms if l.endswith('s1^2'))
        print(f"  s1^2-coefficient polynomial A(p,X) = {sp.expand(A)}")
        print(f"  s1^1-coefficient polynomial B(p,X) = {sp.expand(B)}")
        print(f"  A is zero: {sp.expand(A) == 0}; "
              f"B divisible by modular eq Q: {sp.rem(sp.Poly(B, X_s), sp.Poly(Q, X_s)).as_expr() == 0 if B != 0 else 'n/a (B=0)'}")
        print("  verifying at other q:")
        ok = True
        for q in ("0.2", "0.3", "0.5", "0.6"):
            p, X, s = s1_X7(mp.mpf(q))
            tot = mp.mpf(0)
            scale = mp.mpf(0)
            for c, lab in zip(rel, labels):
                if c == 0:
                    continue
                parts = lab.split()
                v = mp.mpf(1)
                for part in parts:
                    base, ex = part.split("^")
                    v *= {"p": p, "X": X, "s1": s}[base] ** int(ex)
                tot += c * v
                scale = max(scale, abs(v))
            print(f"    q={q}: rel err = {mp.nstr(abs(tot) / scale, 5)}")
            if abs(tot) / scale > mp.mpf(10) ** (-30):
                ok = False
        print(f"  calibration relation VERIFIED at all q: {ok}")

    print()
    print("=== p=11, s_1 over Q(p, m) ===")
    for D in (2, 3):
        labels, rel, _ = pslq_relation(11, mp.mpf("0.4"), D)
        print(f"  D={D} ({len(labels)} monomials), q=0.4:")
        if rel is None:
            print("    No relation found.")
            continue
        terms = [(int(c), l) for c, l in zip(rel, labels) if c != 0]
        print(f"    {len(terms)} nonzero terms:")
        for c, l in terms:
            print(f"      {c:+d} * {l}")
        print("    verifying at other q:")
        ok = verify_relation(11, labels, rel, ["0.2", "0.3", "0.5", "0.6"])
        print(f"    relation VERIFIED at all q: {ok}")

    print()
    print("=== minimal polynomial of m at the flagship point q1 = e^{-pi/sqrt 11} ===")
    q1 = mp.e ** (-mp.pi / mp.sqrt(11))
    m1 = (phi_num(q1, DPS) / phi_num(q1 ** 11, DPS)) ** 2
    print(f"  m = {mp.nstr(m1, 40)}")
    for deg in range(1, 9):
        rel = mp.pslq([m1 ** i for i in range(deg + 1)], maxcoeff=10 ** 8,
                      maxsteps=20000)
        print(f"  degree <={deg}: {rel}")


if __name__ == "__main__":
    main()
