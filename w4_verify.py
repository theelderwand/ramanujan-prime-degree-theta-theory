"""Fresh re-verification of the Atkin-Lehner w4 mechanism (CONJECTURE.md §2.2).

Claim: w4(m) = (-1/p)*m and w4(P) = P for p >= 5
       (p=3: w4(m) = -m, w4(P) = zeta_3 * P).
w4 is represented by W4 = [[4x, y], [4pz, 4w]] with 4xw - pyz = 1, det = 4.
  p = 3 mod 4: y = z = 1, x = (p+1)/4, w = 1.
  p = 1 mod 4: y = 3, z = 1, x = (3p+1)/4, w = 1.
eta via pentagonal series at 40 digits. Branch: exp(pi*i*tau/12), never q**(1/24).
"""
from mpmath import mp, mpf, mpc, exp, pi, fabs
mp.dps = 60

def eta(tau, terms=400):
    """Dedekind eta via pentagonal theorem: eta = e^{pi i tau/12} * sum_k (-1)^k q^{k(3k-1)/2}."""
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

def w4_matrix(p):
    if p % 4 == 3:
        x, y, z, w = (p + 1) // 4, 1, 1, 1
    else:
        x, y, z, w = (3 * p + 1) // 4, 3, 1, 1
    assert 4 * x * w - p * y * z == 1, (p, x, y, z, w)
    assert (4 * x) * (4 * w) - y * (4 * p * z) == 4
    return (4 * x, y, 4 * p * z, 4 * w)

def w4_tau(p, tau):
    a, b, c, d = w4_matrix(p)
    return (a * tau + b) / (c * tau + d)

taus = [mpc(0, 1), mpc(0, 2), mpc('0.3', '1.5')]
print("w4 verification (40 digits)")
print("p | tau | m(w4t)/m(t) vs (-1/p) | P(w4t)/P(t) vs 1")
all_ok = True
for p in [3, 5, 7, 11, 13, 17]:
    leg = 1 if p % 4 == 1 else -1  # (-1/p) = +1 iff p = 1 mod 4
    for tau in taus:
        t2 = w4_tau(p, tau)
        rm = m_val(p, t2) / m_val(p, tau)
        rP = P_val(p, t2) / P_val(p, tau)
        if p == 3:
            target_P = exp(2 * pi * 1j / 3)
            err_m = fabs(rm - leg)
            err_P = fabs(rP - target_P)
            ok = err_m < mpf('1e-30') and err_P < mpf('1e-30')
            print(f"p={p} tau={tau}: m-ratio err={mp.nstr(err_m,3)} P-ratio err vs zeta3={mp.nstr(err_P,3)} {'OK' if ok else 'FAIL'}")
        else:
            err_m = fabs(rm - leg)
            err_P = fabs(rP - 1)
            ok = err_m < mpf('1e-30') and err_P < mpf('1e-30')
            print(f"p={p} tau={tau}: m-ratio err={mp.nstr(err_m,3)} P-ratio err={mp.nstr(err_P,3)} {'OK' if ok else 'FAIL'}")
        all_ok = all_ok and ok
print("ALL OK" if all_ok else "SOME FAILED")
