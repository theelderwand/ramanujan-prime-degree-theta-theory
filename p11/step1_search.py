"""p11 step 1: search for the polynomial relation F(m, p) = 0 by exact-series
linear algebra, detected modulo primes (numpy float64 convolution is exact
for prime ~1e6: partial sums < 6e14 < 2^53).

Degree bounds from step 0 (Ligozat, level 44):
  deg(m) = 5  (sole pole: order 5 at cusp 1/22)
  deg(p) = 10 (sole pole: order 10 at cusp 1/22)
Hence for the minimal relation F(m,p)=0: deg_m F <= deg(p) = 10,
deg_p F <= deg(m) = 5.  Search boxes (I,J) with I<=10, J<=5, in increasing
total size.  Variables: m and P = p/32 (integer series, starts q^5) to keep
coefficients small; the relation is then converted back to (m, p).
"""
import sys
sys.path.insert(0, "/home/hatch/workspace/theta-phase1")
import numpy as np
from theta import (phi_list, chi_list, series_pow, series_div, series_scale)

N = 400  # series order

def series_mod(c, prime):
    return np.array([x % prime for x in c[:N+1]], dtype=np.float64)

def conv_mod(a, b, prime):
    c = np.convolve(a, b)[:N+1]
    return np.mod(np.round(c), prime)

def build_monomials(mser, pser, I, J, prime):
    """Return (labels, matrix rows=coeffs t^0..t^{R-1}) for monomials m^i P^j."""
    mpow = [np.ones(N+1)]
    for i in range(1, I+1):
        mpow.append(conv_mod(mpow[-1], mser, prime))
    ppow = [np.ones(N+1)]
    for j in range(1, J+1):
        ppow.append(conv_mod(ppow[-1], pser, prime))
    labels, cols = [], []
    for j in range(J+1):
        for i in range(I+1):
            cols.append(conv_mod(mpow[i], ppow[j], prime))
            labels.append((i, j))
    return labels, cols

def nullspace_mod(cols, nrows, prime):
    """Nullspace of the matrix A (nrows x ncols) over F_prime, via Gauss
    elimination on A^T? We solve A x = 0 with A rows=coeffs, cols=monomials.
    Returns list of null vectors (as int lists mod prime), or []."""
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
    phi11 = [0]*(N+1)
    for i, c in enumerate(phi_list(N//11)):
        phi11[11*i] = c
    m = series_div(phi2, series_pow(phi11, 2, N), N)
    chi = chi_list(N)
    chi11 = [0]*(N+1)
    for i, c in enumerate(chi_list(N//11)):
        chi11[11*i] = c
    ratio = series_div(chi, series_pow(chi11, 11, N), N)
    P = [0]*(N+1)  # P = p/32 = q^5 * ratio
    for i, c in enumerate(ratio):
        if i+5 <= N:
            P[i+5] = c
    print("m[:6] =", m[:6], " P[:8] =", P[:8])

    primes = [999983, 999979]
    boxes = [(1,1),(2,1),(2,2),(3,2),(3,3),(4,2),(4,3),(5,3),(5,4),
             (6,3),(6,4),(7,4),(8,4),(8,5),(9,5),(10,5)]
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
