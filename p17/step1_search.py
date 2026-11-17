"""p17 step 1: search for the polynomial relation F(m, P) = 0 by exact-series
linear algebra, detected modulo primes.

Degree bounds from step 0 (Ligozat, level 68):
  deg(m) = 8  (sole pole: order 8 at cusp 1/34)
  deg(P) = 24 (sole pole: order 24 at cusp 1/34)
Hence for the minimal relation F(m,P)=0: deg_m F <= deg(P) = 24,
deg_P F <= deg(m) = 8.  Search boxes (I,J) in increasing total size.
Variables: m and P = p_prod/256 (integer series, starts q^12).

Exactness: convolutions in int64 with residues mod prime ~1e6
(products < 1e12, sums of <=401 terms < 4e14 < 2^63: exact).
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
import numpy as np
from theta import (phi_list, chi_list, series_pow, series_div)

P_PRIME = 17
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
    print("m[:6] =", m[:6], " P[:14] =", P[:14])
    print("max |m| coeff:", max(abs(x) for x in m), " max |P| coeff:",
          max(abs(x) for x in P))

    primes = [999983, 999979]
    boxes = [(4,2),(6,2),(6,3),(8,3),(8,4),(10,4),(12,4),(12,5),
             (14,5),(16,5),(16,6),(18,6),(20,6),(20,7),(24,7),(24,8)]
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
                print(f"box I={I} J={J} ({ncols} monomials): no relation mod {prime}",
                      flush=True)
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
