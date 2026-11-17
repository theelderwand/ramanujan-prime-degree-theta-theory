"""p17: summary of results (for the final report).

THEOREM (degree-17 modular equation).  Let m = phi^2(q)/phi^2(q^17),
p_prod = 256 q^12 chi(q)/chi(q^17)^17, P = p_prod/256.  Then F(m,P) = 0,
where (support: 35 monomials, triangular pattern i + 3j <= 12)
  F(m,P) = (m-1)^12 + (256P-1)^4 - 1
           - 1024*P*m*(65536*P^2*m^2 + 3604480*P^2*m + 6488064*P^2
             - 384*P*m^5 + 19008*P*m^4 - 104832*P*m^3 + 188288*P*m^2
             - 898432*P*m + 190528*P
             + m^8 + 24*m^7 - 97*m^6 + 727*m^5 - 19*m^4 - 5*m^3
             + 10165*m^2 + 5205*m + 382)
i.e. explicitly (coefficients by (i, j) = (m-power, P-power)):
  P^0: m^12 - 12 m^11 + 66 m^10 - 220 m^9 + 495 m^8 - 792 m^7 + 924 m^6
       - 792 m^5 + 495 m^4 - 220 m^3 + 66 m^2 - 12 m + 1   [= (m-1)^12]
  P^1: -1024 - 391168 m - 5329920 m^2 - 10408960 m^3 + 5120 m^4 + 19456 m^5
       - 744448 m^6 + 99328 m^7 - 24576 m^8 - 1024 m^9
  P^2: 393216 - 195100672 m + 919994368 m^2 - 192806912 m^3
       + 107347968 m^4 - 19464192 m^5 + 393216 m^6
  P^3: -67108864 - 6643777536 m - 3690987520 m^2 - 67108864 m^3
  P^4: 4294967296
Proved by Sturm bound on Gamma0(68): K = eta(2z)^2 eta(34z)^2 (ord 3 at
every cusp), K^32 weight 64, bound 576, all 35 terms checked (Newman,
trivial character, holomorphic at all cusps), coefficients q^0..q^576 of
K^32 F vanish (no nonzero coefficient to q^800).  F is irreducible over Q.
Pure parts: F(m,0) = (m-1)^12, F(0,P) = (256P-1)^4.  Cusp values:
F(1,0) = 0, F(17,1) = 0.  Flagship: m = sqrt(17),
P = p_prod/256 = 0.00017769354832319060602 (rel. err ~1e-51, 7 points).
Curve F=0: plane curve of degree 12, arithmetic genus 55; one point at
infinity, unibranch, multiplicity sequence 8, 4, 4 (delta=40); twelve
affine nodes (delta=1 each); genus = 55 - 40 - 12 = 3 (confirmed by
Riemann-Hurwitz on m, R=10, and by the degree-2 map X0(68) -> C,
g(X0(68))=7, R'=4).  Discriminant:
  disc_P(F) = -1267650600228229401496703205376 * m^3
    * (16 m^6 - 49 m^5 + 372 m^4 + 1338 m^3 + 6324 m^2 - 14161 m + 78608)
    * S(m)^2
where S(m) = 15625 m^12 + 57350 m^11 + 435636 m^10 + 4017474 m^9
  - 7136763 m^8 + 58172316 m^7 - 2075760 m^6 + 39725924 m^5
  + 142214563 m^4 + 63554366 m^3 + 7265596 m^2 - 57222 m - 4913
is the degree-12 polynomial whose roots are the node m-coordinates.
The fibre over m=0 is the single point (0, 1/256), smooth, with e=4.
Unlike p=13, no odd-power obstruction: the relation is in (m, P) with the
triangular support i + 3j <= 12.
"""
