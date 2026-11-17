"""E1. Universal dissection + p_prod formula.

(a) Exact q-series check in t = q^{1/n}:
      phi(t)/phi(t^{n^2}) = 1 + sum_k u_k(t),
      prod_k u_k(t) = 2^{(n-1)/2} t^{n(n^2-1)/24} chi(t^n)/chi(t^{n^2})^n,
    for n = 3, 5, 7, 9, 11, 13, 17.
(b) mpmath numeric check of the dissection and the p_prod chi-formula at
    several q values, 50 digits.
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
import mpmath as mp
from theta import (phi_list, phi_pow_list, f_list, chi_list, series_mul,
                   series_pow, series_div, series_shift, series_scale,
                   series_add, series_sub, is_zero, first_nonzero,
                   phi_num, f_num, chi_num, u_k_num, p_prod_num, p_prod_chi)

DPS = 50
N = 120  # series order in t


def dissection_series(n, N):
    """Returns (LHS - RHS) as a series in t, to order N. Should be zero."""
    phi_t = phi_list(N)
    phi_tn2 = phi_pow_list(n * n, N)
    inv = series_div([1] + [0] * N, phi_tn2, N)  # 1/phi(t^{n^2})
    LHS = series_mul(phi_t, inv, N)
    terms = []
    for k in range(1, (n - 1) // 2 + 1):
        fk = f_list(n * (n + 2 * k), n * (n - 2 * k), N)
        uk = series_mul(series_shift(fk, k * k, N), inv, N)
        uk = series_scale(uk, 2, N)
        terms.append(uk)
    RHS = series_add([1] + [0] * N, *terms)
    return series_sub(LHS, RHS), terms


def pprod_series(n, N):
    """Returns (prod u_k - chi-formula) as a series in t. Should be zero."""
    phi_tn2 = phi_pow_list(n * n, N)
    inv = series_div([1] + [0] * N, phi_tn2, N)
    prod = [1] + [0] * N
    for k in range(1, (n - 1) // 2 + 1):
        fk = f_list(n * (n + 2 * k), n * (n - 2 * k), N)
        uk = series_mul(series_shift(fk, k * k, N), inv, N)
        uk = series_scale(uk, 2, N)
        prod = series_mul(prod, uk, N)
    # chi formula: 2^{(n-1)/2} t^{n(n^2-1)/24} chi(t^n) / chi(t^{n^2})^n
    chi_tn = chi_list(N)
    # chi(t^n): substitute t -> t^n
    chi_tn_sub = [0] * (N + 1)
    for i, c in enumerate(chi_tn):
        if n * i <= N:
            chi_tn_sub[n * i] = c
    chi_tn2 = chi_list(N)
    chi_tn2_sub = [0] * (N + 1)
    for i, c in enumerate(chi_tn2):
        if n * n * i <= N:
            chi_tn2_sub[n * n * i] = c
    ratio = series_div(chi_tn_sub, series_pow(chi_tn2_sub, n, N), N)
    exp = n * (n * n - 1) // 24
    assert n * (n * n - 1) % 24 == 0
    rhs = series_scale(series_shift(ratio, exp, N), 2 ** ((n - 1) // 2), N)
    return series_sub(prod, rhs)


def main():
    print("=== E1(a): exact series check, order t^%d ===" % N)
    for n in (3, 5, 7, 9, 11, 13, 17):
        d, _ = dissection_series(n, N)
        ok = is_zero(d)
        fnz = first_nonzero(d)
        print(f"n={n:2d} dissection LHS-RHS == 0: {ok}"
              + ("" if ok else f"  FIRST NONZERO at t^{fnz[0]}: {fnz[1]}"))
        dp = pprod_series(n, N)
        okp = is_zero(dp)
        fnzp = first_nonzero(dp)
        print(f"n={n:2d} p_prod  prod-chi  == 0: {okp}"
              + ("" if okp else f"  FIRST NONZERO at t^{fnzp[0]}: {fnzp[1]}"))

    print()
    print("=== E1(b): mpmath numeric check, 50 digits ===")
    mp.mp.dps = DPS
    for n in (3, 5, 7, 11, 13, 17):
        for qval in ("0.3", "0.5", "0.7"):
            q = mp.mpf(qval)
            lhs = phi_num(q ** (mp.mpf(1) / n), DPS) / phi_num(q ** n, DPS)
            rhs = mp.mpf(1) + sum(u_k_num(n, k, q, DPS)
                                  for k in range(1, (n - 1) // 2 + 1))
            err = abs(lhs - rhs)
            pp = p_prod_num(n, q, DPS)
            pc = p_prod_chi(n, q, DPS)
            errp = abs(pp - pc) / abs(pc)
            print(f"n={n:2d} q={qval}: |dissection err|={mp.nstr(err, 5)}"
                  f"  |p_prod rel err|={mp.nstr(errp, 5)}")


if __name__ == "__main__":
    main()
