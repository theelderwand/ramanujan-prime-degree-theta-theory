"""E3. Class invariants G_n.

Verify the published closed forms (Berndt-Rebak 2025, Cors. 3.13/3.18/5.18/5.19;
Watson for G_243; Ramanujan for G_3, G_5, G_7, G_11, G_25, G_27) against the
definition G_n = 2^{-1/4} q^{-1/24} chi(q), q = e^{-pi sqrt(n)}, computed
numerically to 50 digits.

Published closed forms used:
  G_3   = 2^{1/12}
  G_5   = ((1+sqrt5)/2)^{1/4}
  G_7   = 2^{1/4}
  G_11  = ((3 sqrt33+17)^{1/3} - (3 sqrt33-17)^{1/3} + 2)/(3*2^{1/4})
  G_25  = (sqrt5+1)/2
  G_27  = 2^{1/12} (2^{1/3}-1)^{-1/3}
  G_99^3  = 2^{1/4}/3 * ((7822+1362 sqrt33)^{1/3} + (7957+1383 sqrt33)^{1/3}
                         + 2 sqrt33 + 13)
  G_729^3 = sqrt2/2 * R^{1/3} * (3(sqrt3+1)/(sqrt3+1-R^{1/3}) - 2),
            R = ((2(sqrt3+1))^{1/3}+1)/((2(sqrt3-1))^{1/3}-1)
  G_125 = (1/2)((11+5 sqrt5)/2)^{1/4}
          * ((1/sqrt5){1+(2(1+tan pi/5)/(1-sin pi/5))^{1/5}
                        +(2(1-tan pi/5)/(1+sin pi/5))^{1/5}}^2 - 1)
  G_625 = (11+5 sqrt5)/4
          * ((tan^2 pi/5 / sqrt5){1+(8(3 cos pi/5+sin pi/5))^{1/5}
                        +(8(3 cos pi/5-sin pi/5))^{1/5}}^2 - 1)
  G_243 = 2^{1/12}(2^{1/3}+2^{2/3}+3^{1/3})/(9-2*3^{4/3})^{1/3}

Also: G_n as a function of n is checked for the reciprocity G_{1/n} = G_n at
n = 1/3, 1/5, 1/7 (used in the flagship machine).
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
import mpmath as mp
from theta import G_num, chi_num

DPS = 50
mp.mp.dps = DPS


def rel(a, b):
    return abs(a - b) / abs(b)


def cbrt(x):
    return mp.sign(x) * abs(x) ** mp.mpf('1/3')


def check(name, n, closed):
    g = G_num(n, DPS)
    print(f"  {name}: G_num = {mp.nstr(g, 25)}")
    print(f"  {'':{len(name)}}  closed = {mp.nstr(closed, 25)}")
    print(f"  {'':{len(name)}}  rel err = {mp.nstr(rel(g, closed), 5)}")
    return rel(g, closed)


def main():
    print("=== E3: class invariants vs definition, 50 digits ===")
    s5 = mp.sqrt(5)
    s33 = mp.sqrt(33)
    s3 = mp.sqrt(3)
    check("G_3", 3, mp.mpf(2) ** mp.mpf('1/12'))
    check("G_5", 5, ((1 + s5) / 2) ** mp.mpf('1/4'))
    check("G_7", 7, mp.mpf(2) ** mp.mpf('1/4'))
    check("G_11", 11, (cbrt(3 * s33 + 17) - cbrt(3 * s33 - 17) + 2) / (3 * mp.mpf(2) ** mp.mpf('1/4')))
    check("G_25", 25, (s5 + 1) / 2)
    check("G_27", 27, mp.mpf(2) ** mp.mpf('1/12') * (mp.mpf(2) ** mp.mpf('1/3') - 1) ** mp.mpf('-1/3'))
    # G_99 (cube)
    g99_cubed = (mp.mpf(2) ** mp.mpf('1/4') / 3) * (
        cbrt(7822 + 1362 * s33) + cbrt(7957 + 1383 * s33) + 2 * s33 + 13)
    g99 = g99_cubed ** mp.mpf('1/3')
    check("G_99", 99, g99)
    # G_729 (cube)
    R = ((2 * (s3 + 1)) ** mp.mpf('1/3') + 1) / ((2 * (s3 - 1)) ** mp.mpf('1/3') - 1)
    g729_cubed = (mp.sqrt(2) / 2) * R ** mp.mpf('1/3') * (
        3 * (s3 + 1) / (s3 + 1 - R ** mp.mpf('1/3')) - 2)
    g729 = g729_cubed ** mp.mpf('1/3')
    check("G_729", 729, g729)
    # G_125
    t1 = (2 * (1 + mp.tan(mp.pi / 5)) / (1 - mp.sin(mp.pi / 5))) ** mp.mpf('1/5')
    t2 = (2 * (1 - mp.tan(mp.pi / 5)) / (1 + mp.sin(mp.pi / 5))) ** mp.mpf('1/5')
    g125 = mp.mpf('1/2') * ((11 + 5 * s5) / 2) ** mp.mpf('1/4') * ((1 + t1 + t2) ** 2 / s5 - 1)
    check("G_125", 125, g125)
    # G_625
    u1 = (8 * (3 * mp.cos(mp.pi / 5) + mp.sin(mp.pi / 5))) ** mp.mpf('1/5')
    u2 = (8 * (3 * mp.cos(mp.pi / 5) - mp.sin(mp.pi / 5))) ** mp.mpf('1/5')
    g625 = ((11 + 5 * s5) / 4) * ((mp.tan(mp.pi / 5) ** 2 / s5) * (1 + u1 + u2) ** 2 - 1)
    check("G_625", 625, g625)
    # G_243 (Watson)
    g243 = (mp.mpf(2) ** mp.mpf('1/12')
            * (mp.mpf(2) ** mp.mpf('1/3') + mp.mpf(2) ** mp.mpf('2/3') + mp.mpf(3) ** mp.mpf('1/3'))
            / (9 - 2 * mp.mpf(3) ** mp.mpf('4/3')) ** mp.mpf('1/3'))
    check("G_243", 243, g243)

    print()
    print("=== reciprocity G_{1/n} = G_n ===")
    for n in (3, 5, 7, 11):
        g1 = G_num(mp.mpf(1) / n, DPS)
        g2 = G_num(n, DPS)
        print(f"  n={n}: rel err = {mp.nstr(rel(g1, g2), 5)}")

    print()
    print("=== G_n needed for the p=11 flagship (scouting) ===")
    # G-formula (report §2.3, verified analytically in E1):
    #   p_prod = 2^{(p-1)/4} G_n / G_{p^2 n}^p  for q = e^{-pi sqrt n}.
    # At the flagship point q1 = e^{-pi/sqrt 11} (n = 1/11):
    #   p = 2^{5/2} G_{1/11} / G_{11}^{11} = 4 sqrt(2) / G_11^{10}.
    g11 = G_num(11, DPS)
    p11 = 4 * mp.sqrt(2) / g11 ** 10
    print(f"  G_11 = {mp.nstr(g11, 25)}")
    print(f"  p_11(flagship q1) = 4 sqrt(2)/G_11^10 = {mp.nstr(p11, 25)}")
    print(f"  (cross-check vs chi-quotient in p11_scout.py: 0.416439434275864095678569298496)")
    g121 = G_num(121, DPS)
    print(f"  G_121 (numeric, no closed form published) = {mp.nstr(g121, 25)}")
    # minimal polynomial guess for G_121 via sympy nsimplify on high-precision value
    import sympy as sp
    mp.mp.dps = 80
    g121_hi = G_num(121, 80)
    print(f"  G_121 to 60 digits: {mp.nstr(g121_hi, 60)}")


if __name__ == "__main__":
    main()
