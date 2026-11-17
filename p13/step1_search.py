"""p13 step 1: search for the polynomial relation F(m, P) = 0 by exact-series
linear algebra, detected modulo primes.

Degree bounds from step 0 (Ligozat, level 52):
  deg(m) = 6  (sole pole: order 6 at cusp 1/26)
  deg(P) = 14 (sole pole: order 14 at cusp 1/26)
Hence for the minimal relation F(m,P)=0: deg_m F <= deg(P) = 14,
deg_P F <= deg(m) = 6.  Search boxes (I,J) with I<=14, J<=6, in increasing
total size.  Variables: m and P = p_prod/64 (integer series, starts q^7).

Exactness: convolutions done in int64 with residues mod prime ~1e6.
Each product of two residues < 1e12; sums of <=401 terms < 4e14 < 2^63,
so int64 np.convolve is exact.  (p11 used float64; int64 is strictly safer.)
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
import numpy as np
from theta import (phi_list, chi_list, series_pow, series_div)

P_PRIME = 13
N = 400  # series order

def series_mod(c, prime):
    return np.array([x % prime for x in c[:N+1]], dtype=np.int64)

def conv_mod(a, b, prime):
    c = np.convolve(a, b)[:N+1]
    return np.mod(c, prime)

def build_monomials(mser, pser, I, J, prime):
    mpow = [np.ones(N+1, dtype=np.int64)]
    for i in range(1, I+1):
        mpow.append(conv_mod(mpow[-1], mser, prime))
    ppow = [np.ones(N+1, dtype=np.int64)]
    for j in range(1, J+1):
        ppow.append(conv_mod(ppow[-1], pser, prime))
    labels, cols = [], []
    for j in range(J+1):
        for i in range(I+1):
            cols.append(conv_mod(mpow[i], ppow[j], prime))
            labels.append((i, j))
    return labels, cols

def nullspace_mod(cols, nrows, prime):
    """Nullspace of A (nrows x ncols) over F_prime, Gauss-Jordan.
    Returns list of null vectors (int lists mod prime), or []."""
    A = [[int(cols[j][r]) % prime for j in range(len(cols))] for r in range(nrows)]
    ncols = len(cols)
    pivots = []
    r = 0
    for c in range(ncols):
        piv = None
        for i in range(r, nrows):
            if A[i][c] % prime != 0:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], -1, prime)
        A[r] = [(x * inv) % prime for x in A[r]]
        for i in range(nrows):
            if i != r and A[i][c] % prime != 0:
                f = A[i][c]
                A[i] = [(A[i][k] - f * A[r][k]) % prime for k in range(ncols)]
        pivots.append(c)
        r += 1
        if r == nrows:
            break
    free = [c for c in range(ncols) if c not in pivots]
    vecs = []
    for f in free:
        v = [0]*ncols
        v[f] = 1
        for i, pc in enumerate(pivots):
            v[pc] = (-A[i][f]) % prime
        vecs.append(v)
    return vecs

def main():
    phi = phi_list(N)
    phi2 = series_pow(phi, 2, N)
    phip = [0]*(N+1)
    for i, c in enumerate(phi_list(N//P_PRIME)):
        phip[P_PRIME*i] = c
    m = series_div(phi2, series_pow(phip, 2, N), N)
    chi = chi_list(N)
    chip = [0]*(N+1)
    for i, c in enumerate(chi_list(N//P_PRIME)):
        chip[P_PRIME*i] = c
    ratio = series_div(chi, series_pow(chip, P_PRIME, N), N)
    SHIFT_P = (P_PRIME*P_PRIME - 1)//24
    P = [0]*(N+1)
    for i, c in enumerate(ratio):
        if i+SHIFT_P <= N:
            P[i+SHIFT_P] = c
    print("m[:6] =", m[:6], " P[:10] =", P[:10])
    print("max |m| coeff:", max(abs(x) for x in m), " max |P| coeff:",
          max(abs(x) for x in P))

    primes = [999983, 999979]
    boxes = [(1,1),(2,1),(2,2),(3,2),(3,3),(4,2),(4,3),(5,3),(5,4),
             (6,3),(6,4),(7,4),(8,4),(8,5),(9,5),(10,5),(10,6),(12,6),(14,6)]
    for (I, J) in boxes:
        ncols = (I+1)*(J+1)
        nrows = ncols + 40
        found = None
        for prime in primes:
            mser = series_mod(m, prime)
            pser = series_mod(P, prime)
            labels, cols = build_monomials(mser, pser, I, J, prime)
            vecs = nullspace_mod(cols, nrows, prime)
            if not vecs:
                print(f"box I={I} J={J} ({ncols} monomials): no relation mod {prime}")
                found = False
                break
            found = vecs
        if found:
            print(f"box I={I} J={J} ({ncols} monomials): NULLSPACE dim {len(found)} mod {primes[0]}")
            v = found[0]
            nz = [(labels[k], v[k]) for k in range(len(labels)) if v[k] != 0]
            print("  first null vector (mod %d), nonzero terms:" % primes[0])
            for (i, j), c in nz:
                print(f"    m^{i} P^{j}: {c}")
            break
    print("done")

if __name__ == "__main__":
    main()
