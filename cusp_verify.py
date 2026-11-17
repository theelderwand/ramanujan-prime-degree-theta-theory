"""Fresh re-verification of cusp values (CONJECTURE.md §2.6).

Claim: m(1/p) = (-1/p) and P(1/2) = (2/p)/2^{(p-1)/2},
for p = 11, 13, 17, 19, 23.
Direct numerical limit: tau = 1/p + i*t (t -> 0) for m; tau = 1/2 + i*t for P.
High precision (200 digits) to survive the cancellation; Richardson check in t.
"""
from mpmath import mp, mpf, mpc, exp, pi, fabs
mp.dps = 200

def eta(tau, terms=300):
    q = exp(2 * pi * 1j * tau)
    s = mpc(0)
    for k in range(-terms, terms + 1):
        e = k * (3 * k - 1) // 2
        s += (-1) ** k * q ** e
    return exp(pi * 1j * tau / 12) * s

def m_val(p, tau):
    return (eta(tau) ** -4 * eta(2 * tau) ** 10 * eta(4 * tau) ** -4
            * eta(p * tau) ** 4 * eta(2 * p * tau) ** -10 * eta(4 * p * tau) ** 4)

def P_val(p, tau):
    return (eta(tau) ** -1 * eta(2 * tau) ** 2 * eta(4 * tau) ** -1
            * eta(p * tau) ** p * eta(2 * p * tau) ** (-2 * p) * eta(4 * p * tau) ** p)

print("Cusp values (200 digits)")
print("== m(1/p) should be (-1/p) = -1 for p=3 mod 4, +1 for p=1 mod 4 ==")
for p in [11, 13, 17, 19, 23]:
    target = mpf(1) if p % 4 == 1 else mpf(-1)
    vals = []
    for t in [mpf('0.1'), mpf('0.05'), mpf('0.02')]:
        tau = mpc(mpf(1) / p, t)
        vals.append(m_val(p, tau))
    # Richardson extrapolation in t (assume error ~ c1*t + c2*t^2): use t=0.05,0.02
    v1, v2 = vals[1], vals[2]
    t1, t2 = mpf('0.05'), mpf('0.02')
    # linear extrapolation to t=0 from the two smallest t
    extr = (v2 * t1 - v1 * t2) / (t1 - t2)
    print(f"p={p}: t=0.1 -> {mp.nstr(vals[0],10)} | t=0.05 -> {mp.nstr(vals[1],10)} | "
          f"t=0.02 -> {mp.nstr(vals[2],10)} | extrap -> {mp.nstr(extr,10)} | "
          f"target {mp.nstr(target,3)} | err {mp.nstr(fabs(extr-target),3)}")

print("== P(1/2) should be (2/p)/2^{(p-1)/2} ==")
for p in [11, 13, 17, 19, 23]:
    # (2/p) = +1 iff p = +-1 mod 8, -1 iff p = +-3 mod 8
    leg2 = 1 if p % 8 in (1, 7) else -1
    target = mpf(leg2) / mpf(2) ** ((p - 1) // 2)
    vals = []
    for t in [mpf('0.1'), mpf('0.05'), mpf('0.02')]:
        tau = mpc(mpf('0.5'), t)
        vals.append(P_val(p, tau))
    v1, v2 = vals[1], vals[2]
    t1, t2 = mpf('0.05'), mpf('0.02')
    extr = (v2 * t1 - v1 * t2) / (t1 - t2)
    print(f"p={p}: t=0.1 -> {mp.nstr(vals[0],10)} | t=0.05 -> {mp.nstr(vals[1],10)} | "
          f"t=0.02 -> {mp.nstr(vals[2],10)} | extrap -> {mp.nstr(extr,10)} | "
          f"target {mp.nstr(target,10)} | err {mp.nstr(fabs(extr-target),3)}")
