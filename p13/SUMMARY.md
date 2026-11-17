"""p13: summary of results (for the final report).

THEOREM (degree-13 modular equation).  Let m = phi^2(q)/phi^2(q^13),
p_prod = 64 q^7 chi(q)/chi(q^13)^13, P = p_prod/64.  Then F(m,P) = 0, where
  F(m,P) = -(m-1)^7 + (1+64P)^3 - 1 + 64*P*m*(704*P*m - 2944*P
             + 19*m^3 - 44*m^2 + 162*m + 116)
i.e. explicitly:
  F(m,P) = 262144 P^3 + 45056 P^2 m^2 - 188416 P^2 m + 12288 P^2
         + 1216 P m^4 - 2816 P m^3 + 10368 P m^2 + 7424 P m + 192 P
         - m^7 + 7 m^6 - 21 m^5 + 35 m^4 - 35 m^3 + 21 m^2 - 7 m + 1.
Proved by Sturm bound on Gamma0(52): K = eta(2z)^2 eta(26z)^2, K^18 weight 36,
bound 252, all 17 terms checked, coefficients q^0..q^252 of K^18 F vanish
(no nonzero coefficient to q^400).  F is irreducible over Q.
Pure parts: F(m,0) = -(m-1)^7, F(0,P) = (1+64P)^3.
Curve F=0: plane septic, arithmetic genus 15; one point at infinity,
unibranch (4,7), delta=9; four affine nodes
  (1, 1/4), (33.54542973059973, -93.17083713060966),
  (-2.382790973655804, -0.086883662251872), (-0.162638756943925, 0.007720792861530),
each delta=1; genus = 15 - 9 - 4 = 2 (confirmed by Riemann-Hurwitz on m,
R=10, and by the degree-2 map X0(52) -> C, g(X0(52))=5, R'=4).
Note: unlike p=11, odd powers of m occur; the relation is in (m,P), not (Y,P).
"""
