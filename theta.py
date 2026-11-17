"""Core theta-function machinery for Phase 1 of the Ramanujan prime-degree
theta theory programme.

Exact integer q-series (lists of int coefficients) plus high-precision
mpmath numerics.

Conventions (following Berndt-Rebak 2025, Reb\'ak 2023):
  f(a,b)   = sum_{m=-inf}^{inf} a^{m(m+1)/2} b^{m(m-1)/2}   (Jacobi triple product)
  phi(q)   = f(q,q) = sum_{m} q^{m^2}
  chi(q)   = (-q;q^2)_inf = prod_{m>=0} (1+q^{2m+1})
  G_n      = 2^{-1/4} q^{-1/24} chi(q),  q = e^{-pi sqrt(n)}

The universal dissection (Entry 31, Ch. 16; Kim-Rebak Lemma 2.3), for odd n,
in the variable t = q^{1/n}:
  phi(t)/phi(t^{n^2}) = 1 + sum_{k=1}^{(n-1)/2} 2 t^{k^2} f(t^{n(n+2k)}, t^{n(n-2k)})/phi(t^{n^2})
with u_k(t) = 2 t^{k^2} f(t^{n(n+2k)}, t^{n(n-2k)})/phi(t^{n^2}).
"""
import mpmath as mp

# ---------------------------------------------------------------- exact series
# A series is a list c[0..N] of ints (coeff of t^i).

def _trim(c):
    while c and c[-1] == 0:
        c.pop()
    return c


def phi_list(N):
    """phi(t) = sum_m t^{m^2} to order N."""
    c = [0] * (N + 1)
    m = 0
    while m * m <= N:
        c[m * m] += 1 if m == 0 else 2
        m += 1
    return c


def phi_pow_list(pow_n, N):
    """phi(t^{pow_n}) = sum_m t^{pow_n * m^2} to order N."""
    c = [0] * (N + 1)
    m = 0
    while pow_n * m * m <= N:
        c[pow_n * m * m] += 1 if m == 0 else 2
        m += 1
    return c


def f_list(A, B, N):
    """f(t^A, t^B) = sum_m t^{A m(m+1)/2 + B m(m-1)/2} to order N."""
    c = [0] * (N + 1)
    # exponent e(m) = (A+B)/2 m^2 + (A-B)/2 m ; grows like m^2, symmetric-ish
    m = 0
    while True:
        e = A * m * (m + 1) // 2 + B * m * (m - 1) // 2
        if e > N and m > 0:
            # for m>0 exponents increase; but negative m also contributes
            break
        if 0 <= e <= N:
            c[e] += 1
        m += 1
        if m > 10 * N:
            break
    m = -1
    while True:
        e = A * m * (m + 1) // 2 + B * m * (m - 1) // 2
        if e > N and m < 0:
            break
        if 0 <= e <= N:
            c[e] += 1
        m -= 1
        if m < -10 * N:
            break
    return c


def chi_list(N):
    """chi(t) = prod_{m>=0} (1 + t^{2m+1}) to order N."""
    c = [1] + [0] * N
    m = 0
    while 2 * m + 1 <= N:
        d = 2 * m + 1
        # multiply by (1 + t^d)
        for i in range(N - d, -1, -1):
            c[i + d] += c[i]
        m += 1
    return c


def series_mul(a, b, N):
    c = [0] * (N + 1)
    for i in range(N + 1):
        if a[i]:
            ai = a[i]
            for j in range(N + 1 - i):
                if b[j]:
                    c[i + j] += ai * b[j]
    return c


def series_pow(a, e, N):
    """a^e for integer e (a[0] must be 1 for e<0)."""
    if e == 0:
        return [1] + [0] * N
    if e < 0:
        a = series_inv(a, N)
        e = -e
    r = [1] + [0] * N
    base = a[:]
    while e:
        if e & 1:
            r = series_mul(r, base, N)
        base = series_mul(base, base, N)
        e >>= 1
    return r


def series_inv(a, N):
    """1/a for a[0] == 1."""
    assert a[0] == 1
    inv = [0] * (N + 1)
    inv[0] = 1
    for n in range(1, N + 1):
        inv[n] = -sum(a[k] * inv[n - k] for k in range(1, n + 1))
    return inv


def series_div(a, b, N):
    """a/b with b[0] == 1."""
    return series_mul(a, series_inv(b, N), N)


def series_shift(c, s, N):
    """t^s * c(t), truncated to order N."""
    if s > N:
        return [0] * (N + 1)
    return [0] * s + c[: N + 1 - s]


def series_scale(c, k, N):
    return [k * x for x in c[: N + 1]]


def series_add(*cs):
    N = len(cs[0]) - 1
    out = [0] * (N + 1)
    for c in cs:
        for i in range(N + 1):
            out[i] += c[i]
    return out


def series_sub(a, b):
    return [x - y for x, y in zip(a, b)]


def is_zero(c):
    return all(x == 0 for x in c)


def first_nonzero(c):
    for i, x in enumerate(c):
        if x:
            return i, x
    return None


# ------------------------------------------------------------- mpmath numerics
def phi_num(q, dps=50):
    """phi(q) = sum_m q^{m^2}, high precision. q an mpmath number, |q|<1."""
    mp.mp.dps = dps
    s = mp.mpf(1)
    m = 1
    while True:
        t = 2 * q ** (m * m)
        s += t
        if abs(t) < mp.mpf(10) ** (-(dps + 10)):
            break
        m += 1
        if m > 10000:
            break
    return s


def f_num(a, b, dps=50):
    """f(a,b) = sum_m a^{m(m+1)/2} b^{m(m-1)/2}, high precision."""
    mp.mp.dps = dps
    s = mp.mpf(1)
    # m > 0 and m < 0
    for sgn in (1, -1):
        m = sgn
        while True:
            e1 = m * (m + 1) // 2
            e2 = m * (m - 1) // 2
            t = a ** e1 * b ** e2
            s += t
            if abs(t) < mp.mpf(10) ** (-(dps + 10)):
                break
            m += sgn
            if abs(m) > 10000:
                break
    return s


def chi_num(q, dps=50):
    """chi(q) = prod_{m>=0} (1+q^{2m+1}), high precision."""
    mp.mp.dps = dps
    p = mp.mpf(1)
    m = 0
    while True:
        t = q ** (2 * m + 1)
        p *= (1 + t)
        if abs(t) < mp.mpf(10) ** (-(dps + 10)):
            break
        m += 1
        if m > 100000:
            break
    return p


def G_num(n, dps=50):
    """Class invariant G_n = 2^{-1/4} q^{-1/24} chi(q), q = e^{-pi sqrt(n)}."""
    mp.mp.dps = dps
    q = mp.e ** (-mp.pi * mp.sqrt(n))
    return mp.mpf(2) ** mp.mpf('-1/4') * q ** mp.mpf('-1/24') * chi_num(q, dps)


def u_k_num(n, k, q, dps=50):
    """u_k(q) = 2 q^{k^2/n} f(q^{n+2k}, q^{n-2k}) / phi(q^n)."""
    mp.mp.dps = dps
    return (2 * q ** (mp.mpf(k * k) / n)
            * f_num(q ** (n + 2 * k), q ** (n - 2 * k), dps)
            / phi_num(q ** n, dps))


def p_prod_num(n, q, dps=50):
    """p_prod(q) = prod_k u_k(q), k=1..(n-1)/2."""
    mp.mp.dps = dps
    p = mp.mpf(1)
    for k in range(1, (n - 1) // 2 + 1):
        p *= u_k_num(n, k, q, dps)
    return p


def p_prod_chi(n, q, dps=50):
    """Closed form p_prod = 2^{(n-1)/2} q^{(n^2-1)/24} chi(q)/chi(q^n)^n."""
    mp.mp.dps = dps
    return (mp.mpf(2) ** ((n - 1) // 2)
            * q ** (mp.mpf(n * n - 1) / 24)
            * chi_num(q, dps) / chi_num(q ** n, dps) ** n)
