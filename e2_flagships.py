"""E2. The three flagship evaluations, end-to-end from the machine.

Flagships (as stated in Berndt-Rebak 2025 Cor. 3.4, Cor. 5.4; Rebak 2023 Thm 4.1):
  (3) phi(e^{-3 pi sqrt 3}) = 3^{-3/4} phi(e^{-pi sqrt 3}) {1 + 2^{1/3}}
  (5) phi(e^{-5 pi sqrt 5}) = 5^{-3/4} phi(e^{-pi sqrt 5})
          {1 + (2(1+tan pi/5)/(1-sin pi/5))^{1/5} + (2(1-tan pi/5)/(1+sin pi/5))^{1/5}}
  (7) phi(e^{-7 pi sqrt 7}) = 7^{-3/4} phi(e^{-pi sqrt 7})
          {1 + (cos pi/7 / (2 cos^2 2pi/7))^{2/7} + (cos 2pi/7 / (2 cos^2 3pi/7))^{2/7}
               + (cos 3pi/7 / (2 cos^2 pi/7))^{2/7}}

Machine reproduction: at q1 = e^{-pi/sqrt(p)} the dissection gives
  phi(q1^{1/p})/phi(q1^p) = 1 + sum u_k,  q1^{1/p} = e^{-pi/(p sqrt p)},
  q1^p = e^{-pi sqrt p},
and the theta transformation phi(e^{-pi/(p sqrt p)}) = p^{3/4} phi(e^{-p pi sqrt p})
turns this into the flagship. We verify, at q1:
  - p_prod from the chi-quotient vs the G-formula,
  - the multiplier ratio phi^4(q1)/phi^4(q1^p) (septic: = 7),
  - the auxiliary polynomial, its roots, and the cyclic root ordering
    matching u_k = (monomial in roots)^{1/p} to the published trig terms.
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
import mpmath as mp
from theta import phi_num, f_num, chi_num, u_k_num, p_prod_num, p_prod_chi, G_num

DPS = 50
mp.mp.dps = DPS


def err(a, b):
    return abs(a - b) / max(abs(b), mp.mpf(1))


def flagship_cubic():
    print("--- cubic flagship: phi(e^{-3 pi sqrt 3}) ---")
    lhs = phi_num(mp.e ** (-3 * mp.pi * mp.sqrt(3)), DPS)
    rhs = (mp.mpf(3) ** mp.mpf('-3/4') * phi_num(mp.e ** (-mp.pi * mp.sqrt(3)), DPS)
           * (1 + mp.mpf(2) ** mp.mpf('1/3')))
    print(f"  LHS = {mp.nstr(lhs, 30)}")
    print(f"  RHS = {mp.nstr(rhs, 30)}")
    print(f"  rel err = {mp.nstr(err(lhs, rhs), 5)}")


def flagship_quintic():
    print("--- quintic flagship: phi(e^{-5 pi sqrt 5}) ---")
    t1 = (2 * (1 + mp.tan(mp.pi / 5)) / (1 - mp.sin(mp.pi / 5))) ** mp.mpf('1/5')
    t2 = (2 * (1 - mp.tan(mp.pi / 5)) / (1 + mp.sin(mp.pi / 5))) ** mp.mpf('1/5')
    lhs = phi_num(mp.e ** (-5 * mp.pi * mp.sqrt(5)), DPS)
    rhs = (mp.mpf(5) ** mp.mpf('-3/4') * phi_num(mp.e ** (-mp.pi * mp.sqrt(5)), DPS)
           * (1 + t1 + t2))
    print(f"  LHS = {mp.nstr(lhs, 30)}")
    print(f"  RHS = {mp.nstr(rhs, 30)}")
    print(f"  rel err = {mp.nstr(err(lhs, rhs), 5)}")
    return t1, t2


def flagship_septic():
    print("--- septic flagship: phi(e^{-7 pi sqrt 7}) ---")
    c1, c2, c3 = mp.cos(mp.pi / 7), mp.cos(2 * mp.pi / 7), mp.cos(3 * mp.pi / 7)
    s1 = (c1 / (2 * c2 ** 2)) ** mp.mpf('2/7')
    s2 = (c2 / (2 * c3 ** 2)) ** mp.mpf('2/7')
    s3 = (c3 / (2 * c1 ** 2)) ** mp.mpf('2/7')
    lhs = phi_num(mp.e ** (-7 * mp.pi * mp.sqrt(7)), DPS)
    rhs = (mp.mpf(7) ** mp.mpf('-3/4') * phi_num(mp.e ** (-mp.pi * mp.sqrt(7)), DPS)
           * (1 + s1 + s2 + s3))
    print(f"  LHS = {mp.nstr(lhs, 30)}")
    print(f"  RHS = {mp.nstr(rhs, 30)}")
    print(f"  rel err = {mp.nstr(err(lhs, rhs), 5)}")
    return s1, s2, s3


def machine_cubic():
    print("--- machine: cubic at q1 = e^{-pi/sqrt 3} ---")
    q1 = mp.e ** (-mp.pi / mp.sqrt(3))
    p_chi = 2 * q1 ** mp.mpf('1/3') * chi_num(q1, DPS) / chi_num(q1 ** 3, DPS) ** 3
    G3 = mp.mpf(2) ** mp.mpf('1/12')
    p_G = mp.sqrt(2) / G3 ** 2   # sqrt(2) G_{1/3} / G_3^3, G_{1/3}=G_3
    print(f"  p via chi-quotient = {mp.nstr(p_chi, 30)}")
    print(f"  p via G-formula    = {mp.nstr(p_G, 30)}")
    print(f"  2^(1/3)            = {mp.nstr(mp.mpf(2) ** mp.mpf('1/3'), 30)}")
    print(f"  rel err chi vs 2^(1/3) = {mp.nstr(err(p_chi, mp.mpf(2) ** mp.mpf('1/3')), 5)}")
    u = u_k_num(3, 1, q1, DPS)
    print(f"  u from dissection  = {mp.nstr(u, 30)}  (should equal p)")
    # multiplier relation phi^4(q)/phi^4(q^3) = 1 + p^3
    m4 = (phi_num(q1, DPS) / phi_num(q1 ** 3, DPS)) ** 4
    print(f"  phi^4(q)/phi^4(q^3) = {mp.nstr(m4, 30)}")
    print(f"  1 + p^3             = {mp.nstr(1 + p_chi ** 3, 30)}")
    print(f"  rel err = {mp.nstr(err(m4, 1 + p_chi ** 3), 5)}")


def machine_quintic(t1, t2):
    print("--- machine: quintic at q1 = e^{-pi/sqrt 5} ---")
    q1 = mp.e ** (-mp.pi / mp.sqrt(5))
    p_chi = p_prod_chi(5, q1, DPS)
    G5 = ((1 + mp.sqrt(5)) / 2) ** mp.mpf('1/4')
    p_G = 2 / G5 ** 4   # 2 G_{1/5} / G_5^5 = 2/G_5^4
    p_exp = mp.sqrt(5) - 1
    print(f"  p via chi-quotient = {mp.nstr(p_chi, 30)}")
    print(f"  p via G-formula    = {mp.nstr(p_G, 30)}")
    print(f"  sqrt(5)-1          = {mp.nstr(p_exp, 30)}")
    print(f"  rel err chi vs sqrt5-1 = {mp.nstr(err(p_chi, p_exp), 5)}")
    # multiplier relation phi^2(q)/phi^2(q^5) = 1 + p
    m2 = (phi_num(q1, DPS) / phi_num(q1 ** 5, DPS)) ** 2
    print(f"  phi^2(q)/phi^2(q^5) = {mp.nstr(m2, 30)}")
    print(f"  1 + p               = {mp.nstr(1 + p_chi, 30)}")
    print(f"  rel err = {mp.nstr(err(m2, 1 + p_chi), 5)}")
    # auxiliary quadratic xi^2 - ((p-1)^2+7) xi + p^3 = 0
    p = p_chi
    a = -((p - 1) ** 2 + 7)
    b = p ** 3
    disc = mp.sqrt(a * a - 4 * b)
    alpha = (-a + disc) / 2
    beta = (-a - disc) / 2
    print(f"  aux quad roots: alpha = {mp.nstr(alpha, 25)}, beta = {mp.nstr(beta, 25)}")
    u = u_k_num(5, 1, q1, DPS)
    v = u_k_num(5, 2, q1, DPS)
    print(f"  u,v from dissection = {mp.nstr(u, 25)}, {mp.nstr(v, 25)}")
    u_pred = (alpha * p) ** mp.mpf('1/5')
    v_pred = (beta * p) ** mp.mpf('1/5')
    print(f"  (alpha p)^(1/5)     = {mp.nstr(u_pred, 25)}")
    print(f"  (beta  p)^(1/5)     = {mp.nstr(v_pred, 25)}")
    print(f"  published t1        = {mp.nstr(t1, 25)}")
    print(f"  published t2        = {mp.nstr(t2, 25)}")
    print(f"  rel err u vs t1     = {mp.nstr(err(u, t1), 5)}")
    print(f"  rel err v vs t2     = {mp.nstr(err(v, t2), 5)}")


def machine_septic(s1, s2, s3):
    print("--- machine: septic at q1 = e^{-pi/sqrt 7} ---")
    q1 = mp.e ** (-mp.pi / mp.sqrt(7))
    p_chi = p_prod_chi(7, q1, DPS)
    G7 = mp.mpf(2) ** mp.mpf('1/4')
    p_G = 2 * mp.sqrt(2) / G7 ** 6   # 2 sqrt2 G_{1/7} / G_7^7 = 2 sqrt2 / G_7^6
    print(f"  p via chi-quotient = {mp.nstr(p_chi, 30)}")
    print(f"  p via G-formula    = {mp.nstr(p_G, 30)}  (should be 1)")
    print(f"  rel err vs 1       = {mp.nstr(err(p_chi, mp.mpf(1)), 5)}")
    # multiplier relation: X = phi^4(q)/phi^4(q^7) satisfies X^2-(2+5p)X+(1-p)^3=0
    X = (phi_num(q1, DPS) / phi_num(q1 ** 7, DPS)) ** 4
    print(f"  X = phi^4(q)/phi^4(q^7) = {mp.nstr(X, 30)}  (should be 7)")
    print(f"  rel err vs 7 = {mp.nstr(err(X, mp.mpf(7)), 5)}")
    resid = X ** 2 - (2 + 5 * p_chi) * X + (1 - p_chi) ** 3
    print(f"  modular-eq residual = {mp.nstr(resid, 5)}")
    # auxiliary cubic r(xi) = xi^3 + 2 xi^2 (1+3p-X) + xi p^2 (p+4) - p^4
    p = p_chi
    c2 = 2 * (1 + 3 * p - X)
    c1 = p ** 2 * (p + 4)
    c0 = -p ** 4
    print(f"  r(xi) = xi^3 + ({mp.nstr(c2, 10)}) xi^2 + ({mp.nstr(c1, 10)}) xi + ({mp.nstr(c0, 10)})")
    print(f"  (expected xi^3 - 6 xi^2 + 5 xi - 1)")
    # roots numerically (mpmath polyroots)
    roots = mp.polyroots([1, c2, c1, c0], maxsteps=200, extraprec=50)
    roots = sorted([mp.re(r) for r in roots], reverse=True)
    print(f"  roots of r: {[mp.nstr(r, 25) for r in roots]}")
    c1t, c2t, c3t = mp.cos(mp.pi / 7), mp.cos(2 * mp.pi / 7), mp.cos(3 * mp.pi / 7)
    pred_roots = sorted([1 / (2 * c) ** 2 for c in (c1t, c2t, c3t)], reverse=True)
    print(f"  1/(2cos k pi/7)^2: {[mp.nstr(r, 25) for r in pred_roots]}")
    for r, pr in zip(roots, pred_roots):
        print(f"    root rel err vs trig = {mp.nstr(err(r, pr), 5)}")
    # u,v,w from dissection
    u = u_k_num(7, 1, q1, DPS)
    v = u_k_num(7, 2, q1, DPS)
    w = u_k_num(7, 3, q1, DPS)
    print(f"  u,v,w from dissection = {mp.nstr(u, 25)}, {mp.nstr(v, 25)}, {mp.nstr(w, 25)}")
    print(f"  published s1,s2,s3   = {mp.nstr(s1, 25)}, {mp.nstr(s2, 25)}, {mp.nstr(s3, 25)}")
    print(f"  NOTE: published order is (v,u,w) in dissection order: u=s2, v=s1, w=s3")
    print(f"  rel err u vs s2 = {mp.nstr(err(u, s2), 5)}")
    print(f"  rel err v vs s1 = {mp.nstr(err(v, s1), 5)}")
    print(f"  rel err w vs s3 = {mp.nstr(err(w, s3), 5)}")
    # cyclic root ordering: u = (alpha^2 p / beta)^{1/7} etc.
    # roots sorted desc: r1 > r2 > r3 correspond to 1/(2cos k pi/7)^2 for k=3,1,2?
    # cos pi/7 > cos 2pi/7 > cos 3pi/7 so 1/(2cos)^2 is smallest for k=1, largest for k=3.
    # Try all cyclic assignments of (alpha,beta,gamma) = permutation of roots
    import itertools
    best = None
    for perm in itertools.permutations(roots):
        a, b, g = perm
        up = (a ** 2 * p / b) ** mp.mpf('1/7')
        vp = (b ** 2 * p / g) ** mp.mpf('1/7')
        wp = (g ** 2 * p / a) ** mp.mpf('1/7')
        e = err(up, u) + err(vp, v) + err(wp, w)
        if best is None or e < best[0]:
            best = (e, perm, up, vp, wp)
    e, perm, up, vp, wp = best
    print(f"  best cyclic assignment total rel err = {mp.nstr(e, 5)}")
    print(f"  (alpha,beta,gamma) = {[mp.nstr(r, 12) for r in perm]}")
    print(f"    = 1/(2cos k pi/7)^2 for k = "
          f"{[int(round(7 * mp.acos(1 / (2 * mp.sqrt(r))) / mp.pi)) for r in perm]}")
    print(f"  (alpha^2 p/beta)^(1/7) = {mp.nstr(up, 25)} vs u = {mp.nstr(u, 25)}")
    # Lemma 3.8: sum of the 7th powers of the terms = 41, i.e.
    # sum (cos k pi/7 / (2 cos^2 ...))^2 = 41. (s_i = (ratio)^{2/7}, so s_i^7 = ratio^2.)
    ssum = s1 ** 7 + s2 ** 7 + s3 ** 7
    print(f"  sum s_i^7 = {mp.nstr(ssum, 30)}  (Lemma 3.8: 41)")
    print(f"  rel err vs 41 = {mp.nstr(err(ssum, mp.mpf(41)), 5)}")


def main():
    flagship_cubic()
    print()
    t1, t2 = flagship_quintic()
    print()
    s1, s2, s3 = flagship_septic()
    print()
    machine_cubic()
    print()
    machine_quintic(t1, t2)
    print()
    machine_septic(s1, s2, s3)


if __name__ == "__main__":
    main()
