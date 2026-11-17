"""E4. Modular equations via the modern (Sturm-bound) route.

(A) Sturm-bound proofs of the multiplier relations of degrees 3, 5, 7
    (Berndt, Ramanujan's Notebooks Part III; Berndt-Rebak 2025 Thm 2.1(iii),
    Thm 4.2(iii); Son 2001 / Rebak 2023 Entry 1.1(iii)):
      p=3: phi^4(q) = phi^4(q^3) + p^3 phi^4(q^3),   p = 2 q^{1/3} chi(q)/chi^3(q^3)
      p=5: phi^2(q) = phi^2(q^5) + p phi^2(q^5),     p = 4 q chi(q)/chi^5(q^5)
      p=7: X^2-(2+5p)X+(1-p)^3 = 0, X = phi^4(q)/phi^4(q^7),
           p = 8 q^2 chi(q)/chi^7(q^7); cleared: phi^8(q) - (2+5p)phi^4(q)phi^4(q^7)
           + (1-p)^3 phi^8(q^7) = 0.
    Each term is an eta quotient (phi(q) = eta(2z)^5/(eta(z)^2 eta(4z)^2),
    chi(q) = q^{1/24} eta(2z)^2/(eta(z) eta(4z))).  The prover checks, for every
    term: Newman congruences, character (D_eff = (-1)^k prod d^{e_d}; two terms
    share a character iff the ratio of D_eff's is a square in Q), holomorphicity
    at every cusp (Ligozat: ord_{1/c} = (N/(24 gcd(c^2,N))) sum e_d gcd(c,d)^2/d),
    then compares q-expansions to the Sturm bound B = ceil(k [SL2:Gamma0(N)]/12).

(B) Numerical verification (50 digits) of Schroeter's degree-11 modular equation
    (Berndt Part III, Entry 7(i)):
      (ab)^{1/4} + ((1-a)(1-b))^{1/4} + 2 (16 a b (1-a)(1-b))^{1/12} = 1,
    a = lambda(tau), b = lambda(11 tau), lambda = 16 eta(tau/2)^8 eta(2tau)^16/eta(tau)^24.
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/modular-forms/scripts")
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
from fractions import Fraction
from math import gcd
import mpmath as mp
from eta import eta_quotient_qexp
from sturm import index_gamma0


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def is_square_frac(fr):
    fr = Fraction(fr)
    if fr < 0:
        return False
    n, d = fr.numerator, fr.denominator
    rn = int(n ** 0.5)
    rd = int(d ** 0.5)
    # exact integer sqrt check
    import math
    rn = math.isqrt(n)
    rd = math.isqrt(d)
    return rn * rn == n and rd * rd == d


def D_eff(edict, k):
    """(-1)^k prod d^{e_d} as a Fraction."""
    D = Fraction(1)
    for d, e in edict.items():
        D *= Fraction(d) ** e
    if k % 2 == 1:
        D = -D
    return D


def newman_ok(edict, N):
    s1 = sum(d * e for d, e in edict.items())
    s2 = sum((N // d) * e for d, e in edict.items())
    return s1 % 24 == 0 and s2 % 24 == 0, s1, s2


def cusp_ord(edict, N, c):
    """Order at cusp 1/c (c | N); c=0 means infinity (use shift)."""
    if c == 0:
        return sum(Fraction(d * e, 24) for d, e in edict.items())
    fac = Fraction(N, 24 * gcd(c * c, N))
    s = sum(Fraction(e * gcd(c, d) ** 2, d) for d, e in edict.items())
    return fac * s


def weight(edict):
    return Fraction(sum(edict.values()), 2)


def check_term(edict, N, label):
    k = weight(edict)
    ok_n, s1, s2 = newman_ok(edict, N)
    D = D_eff(edict, int(k) if k.denominator == 1 else 0)
    ords = {c: cusp_ord(edict, N, c) for c in divisors(N)}
    ords[0] = cusp_ord(edict, N, 0)
    holo = all(o >= 0 for o in ords.values())
    print(f"  {label}: edict={dict(sorted(edict.items()))}")
    print(f"    weight={k}, Newman(sum d*e={s1}, sum (N/d)*e={s2}): {ok_n}, "
          f"D_eff={D}, cusp orders>={min(ords.values())} (all>=0: {holo})")
    return {"k": k, "D": D, "holo": holo, "newman": ok_n, "ords": ords}


def series_with_shift(edict, Nterms):
    """Return dict {exponent (Fraction): coeff} for the eta quotient."""
    shift, coeffs = eta_quotient_qexp(edict, Nterms)
    return {shift + i: c for i, c in enumerate(coeffs)}


def sturm_prove(terms, N, label, check_to=None):
    """terms: list of (scalar, edict). Prove sum scalar*etaquotient == 0."""
    print(f"--- {label} (level N={N}) ---")
    infos = []
    for i, (s, e) in enumerate(terms):
        infos.append(check_term(e, N, f"term {i} (x{s})"))
    ks = {str(i["k"]) for i in infos}
    assert len(ks) == 1, f"weights differ: {ks}"
    k = infos[0]["k"]
    assert all(i["newman"] and i["holo"] for i in infos), "modularity check failed"
    D0 = infos[0]["D"]
    for i in infos[1:]:
        assert is_square_frac(D0 / i["D"]), f"characters differ: {D0} vs {i['D']}"
    print(f"    all terms weight {k}, same character, holomorphic on Gamma0({N})")
    mu = index_gamma0(N)
    B = int(k * mu / 12)
    if (k * mu / 12) != B:
        B += 1
    print(f"    index mu={mu}, Sturm bound B = ceil({k}*{mu}/12) = {B}")
    nterms = (check_to or B) + 5
    # sum the series with shifts
    total = {}
    for s, e in terms:
        ser = series_with_shift(e, nterms)
        for exp, c in ser.items():
            total[exp] = total.get(exp, 0) + s * c
    bad = [(exp, c) for exp, c in sorted(total.items()) if exp <= B and c != 0]
    print(f"    coefficients 0..{B} of the sum all zero: {not bad}"
          + ("" if not bad else f"  FIRST BAD: {bad[0]}"))
    # also report how far it stays zero (only exponents every term covers: <= nterms)
    mx = None
    for exp in sorted(total):
        if exp <= nterms and total[exp] != 0:
            mx = exp
            break
    print(f"    first nonzero coefficient (if any, within computed range) at q^{mx}")
    return not bad


def eta_prod(tau, dps=50):
    """eta(tau) for tau in the upper half-plane, via the product."""
    mp.mp.dps = dps
    q = mp.e ** (2 * mp.pi * mp.j * tau)
    p = mp.mpc(1)
    n = 1
    while True:
        t = q ** n
        p *= (1 - t)
        if abs(t) < mp.mpf(10) ** (-(dps + 10)):
            break
        n += 1
    return mp.e ** (mp.pi * mp.j * tau / 12) * p


def lam(tau, dps=50):
    """Elliptic lambda function: (theta2/theta3)^4
    = 16 eta(tau/2)^8 eta(2tau)^16 / eta(tau)^24."""
    mp.mp.dps = dps
    return 16 * eta_prod(tau / 2, dps) ** 8 * eta_prod(2 * tau, dps) ** 16 / eta_prod(tau, dps) ** 24


def main():
    mp.mp.dps = 50
    # ---- p = 3 ----
    phi4_q = {1: -8, 2: 20, 4: -8}
    phi4_q3 = {3: -8, 6: 20, 12: -8}
    p3_phi4_q3 = {1: -3, 2: 6, 3: 1, 4: -3, 6: 2, 12: 1}  # p^3 phi^4(q^3)/8
    sturm_prove([(1, phi4_q), (-1, phi4_q3), (-8, p3_phi4_q3)], 12,
                "p=3: phi^4(q) - phi^4(q^3) - p^3 phi^4(q^3) = 0", check_to=60)
    print()
    # ---- p = 5 ----
    phi2_q = {1: -4, 2: 10, 4: -4}
    phi2_q5 = {5: -4, 10: 10, 20: -4}
    p_phi2_q5 = {1: -1, 2: 2, 4: -1, 5: 1, 20: 1}  # p phi^2(q^5)/4
    sturm_prove([(1, phi2_q), (-1, phi2_q5), (-4, p_phi2_q5)], 20,
                "p=5: phi^2(q) - phi^2(q^5) - p phi^2(q^5) = 0", check_to=60)
    print()
    # ---- p = 7 ----
    Ep = {1: -1, 2: 2, 4: -1, 7: 7, 14: -14, 28: 7}  # p/8
    phi8_q = {1: -16, 2: 40, 4: -16}
    phi4_q_phi4_q7 = {1: -8, 2: 20, 4: -8, 7: -8, 14: 20, 28: -8}
    phi8_q7 = {7: -16, 14: 40, 28: -16}

    def add_edicts(*es):
        out = {}
        for e in es:
            for d, x in e.items():
                out[d] = out.get(d, 0) + x
        return {d: x for d, x in out.items() if x != 0}

    Ep2 = {d: 2 * x for d, x in Ep.items()}
    Ep3 = {d: 3 * x for d, x in Ep.items()}
    # terms: phi^8(q) - 2 phi^4(q)phi^4(q^7) - 40 Ep*phi^4phi^4 + phi^8(q^7)
    #        - 24 Ep*phi^8(q^7) + 192 Ep^2*phi^8(q^7) - 512 Ep^3*phi^8(q^7)
    terms7 = [
        (1, phi8_q),
        (-2, phi4_q_phi4_q7),
        (-40, add_edicts(Ep, phi4_q_phi4_q7)),
        (1, phi8_q7),
        (-24, add_edicts(Ep, phi8_q7)),
        (192, add_edicts(Ep2, phi8_q7)),
        (-512, add_edicts(Ep3, phi8_q7)),
    ]
    sturm_prove(terms7, 28,
                "p=7: phi^8(q)-(2+5p)phi^4(q)phi^4(q^7)+(1-p)^3 phi^8(q^7) = 0",
                check_to=60)
    print()
    # ---- degree 11: Schroeter, numerical ----
    print("--- degree 11 (Schroeter, Entry 7(i)): numerical, 50 digits ---")
    for y in ("1", "0.5", "2", "0.3"):
        tau = mp.mpc(0, mp.mpf(y))
        a = lam(tau)
        b = lam(11 * tau)
        lhs = (a * b) ** mp.mpf('1/4') + ((1 - a) * (1 - b)) ** mp.mpf('1/4') \
            + 2 * (16 * a * b * (1 - a) * (1 - b)) ** mp.mpf('1/12')
        print(f"  tau = {y}i: a = {mp.nstr(mp.re(a), 12)}, b = {mp.nstr(mp.re(b), 12)}")
        print(f"    LHS = {mp.nstr(lhs, 30)}  (should be 1); |LHS-1| = {mp.nstr(abs(lhs - 1), 5)}")


if __name__ == "__main__":
    main()
