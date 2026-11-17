# Degree-13 and Degree-17 Modular Equations in Ramanujan's Prime-Degree Theta Function Theory

**Status: both complete — explicit polynomials found and proved by Sturm bound.**
Date: 2026-11-14. Method: identical to the p=11 proof (mod-prime linear algebra
on q-expansions to find the relation; Sturm bound on Γ₀(4p) to prove it).
Scripts: `p13/`, `p17/` (mirroring `p11/`). Paper updated:
`paper/paper.tex` (Theorems 5.2, 5.3).

## Common setup

For prime p, with m = φ²(q)/φ²(q^p) and p_prod = 2^{(p−1)/2} q^{(p²−1)/24}
χ(q)/χ^p(q^p) (Theorem 2.1 of the paper — proved for all odd n), write
P = p_prod / 2^{(p−1)/2}. Both m and P are eta-quotients on Γ₀(4p) with
trivial character (D_eff(m) = 1/p², D_eff(P) = 1, both squares) and their
only pole at the cusp 1/(2p). The relation F(m,P) = 0 is found by exact
linear algebra on q-expansions (mod-prime nullspace, then exact rational
nullspace) and proved by multiplying by a power of K = η(2z)²η(2pz)² and
checking the Sturm bound.

---

## 1. Degree 13

- Level Γ₀(52), index [SL₂(ℤ):Γ₀(52)] = 84.
- m = η(z)^{−4}η(2z)^{10}η(4z)^{−4}η(13z)^4η(26z)^{−10}η(52z)^4,
  deg(m) = 6 (only pole 1/26, order −6), D_eff = 1/169.
- P = p_prod/64 = q⁷χ(q)/χ¹³(q¹³), eta quotient
  {1:−1, 2:2, 4:−1, 13:13, 26:−26, 52:13}, shift 7,
  deg(P) = 14 (only pole 1/26, order −14), D_eff = 1.
- Ligozat cusp orders: m: ∞ 0, 0 0, 1/2 6, 1/4 0, 1/13 0, 1/26 −6;
  P: ∞ 7, 0 0, 1/2 0, 1/4 0, 1/13 7, 1/26 −14.
  (Γ₀(52) has six cusps; 1/52 is Γ₀(52)-equivalent to ∞.)

### Theorem (degree-13 modular equation)

F(m,P) = 0, where

    F(m,P) = 262144·P³ + 45056·P²m² − 188416·P²m + 12288·P²
           + 1216·Pm⁴ − 2816·Pm³ + 10368·Pm² + 7424·Pm + 192·P
           − m⁷ + 7m⁶ − 21m⁵ + 35m⁴ − 35m³ + 21m² − 7m + 1.

Structure: F(m,0) = −(m−1)⁷, F(0,P) = (1+64P)³, and

    F = (1+64P)³ − (m−1)⁷ − 1
        + 64·P·m·(704·Pm − 2944·P + 19m³ − 44m² + 162m + 116).

F is irreducible over ℚ. Unlike p=11, odd powers of m occur — the relation
is in (m,P), not (Y,P) with Y = m². Support: 17 monomials
(m⁰…⁷P⁰, m⁰…⁴P¹, m⁰…²P², m⁰P³); the relation first appears at the monomial
box (7,4), where the nullspace is 2-dimensional (the relation and its
multiple by P).

### Proof (Sturm bound)

K = η(2z)²η(26z)² has ord 7/3 at every cusp of Γ₀(52). K¹⁸ has weight 36
and ord_{1/26}(K¹⁸) = 42, equal to the maximum pole order of m^i P^j over
the support (42). All 17 terms K¹⁸·m^i·P^j are holomorphic weight-36 forms
with trivial character on Γ₀(52) (Newman, character, and Ligozat checks all
pass). Sturm bound ⌈36·84/12⌉ = 252: the coefficients q⁰…q²⁵² of K¹⁸·F are
all zero (in fact no nonzero coefficient to q⁴⁰⁰). Hence K¹⁸·F ≡ 0, so F = 0.

### Verification

- Exact: F(m,P) = 0 as a formal q-series to q⁶⁰⁰.
- Numerical (50 digits): relative error ~10⁻⁵¹ at 7 points; flagship
  m = √13, P ≈ 0.0034695471649933081995.
- Cusp values: F(1,0) = 0 (cusp ∞: m=1, P=0), F(13,1) = 0 (cusp 0: m=13,
  P=1, i.e. p_prod = 64).

### The curve F(m,P) = 0 has genus 2

Plane septic, arithmetic genus 15. One point at infinity, unibranch with
Puiseux pair (4,7), δ = 9. Four affine nodes (δ = 1 each):
(1, 1/4), (33.54542973059973, −93.17083713060966),
(−2.382790973655804, −0.086883662251872),
(−0.162638756943925, 0.007720792861530).
Genus = 15 − 9 − 4 = 2, confirmed by Riemann–Hurwitz on m (R = 8: 2 at m = 0,
4 at the four simple ramification points, 2 at m = ∞; the nodes contribute 0)
and by the degree-2 map X₀(52) → C (g(X₀(52)) = 5, R′ = 4).
Discriminant:
disc_P(F) = −68719476736·m²·(m−1)²·(m³−31m²−85m−13)²·(27m⁴−212m³+1194m²−2756m+4563).
The fibre over m = 0 is the single point (0, −1/64), smooth, with e = 3
(F_m(0,−1/64) = −169 = −13² ≠ 0). [Correction 2026-11-14: an earlier draft
gave this fibre as (0, 1/4); F(0,1/4) = 4913 ≠ 0. The value 1/4 is the
P-coordinate of the node at m = 1.]

Scripts: `p13/step0_explore.py` (eta quotients, cusp orders),
`p13/step1_search.py` (mod-prime nullspace), `p13/step2_exact.py` (exact and
numerical verification), `p13/step3_sturm.py` (Sturm proof),
`p13/step4_structure.py`, `p13/step4b_singular.py`–`p13/step4f_node_check.py`
(genus analysis), `p13/SUMMARY.md`.

---

## 2. Degree 17

- Level Γ₀(68), index [SL₂(ℤ):Γ₀(68)] = 108.
- m = η(z)^{−4}η(2z)^{10}η(4z)^{−4}η(17z)^4η(34z)^{−10}η(68z)^4,
  deg(m) = 8 (only pole 1/34, order −8), D_eff = 1/289.
- P = p_prod/256 = q¹²χ(q)/χ¹⁷(q¹⁷), eta quotient
  {1:−1, 2:2, 4:−1, 17:17, 34:−34, 68:17}, shift 12,
  deg(P) = 24 (only pole 1/34, order −24), D_eff = 1.
- Ligozat cusp orders: m: ∞ 0, 0 0, 1/2 8, 1/4 0, 1/17 0, 1/34 −8;
  P: ∞ 12, 0 0, 1/2 0, 1/4 0, 1/17 12, 1/34 −24.
  (Γ₀(68) has six cusps; 1/68 is Γ₀(68)-equivalent to ∞.)

### Theorem (degree-17 modular equation)

F(m,P) = 0, where (35 monomials; triangular support i + 3j ≤ 12)

    P⁰: m¹² − 12m¹¹ + 66m¹⁰ − 220m⁹ + 495m⁸ − 792m⁷ + 924m⁶
        − 792m⁵ + 495m⁴ − 220m³ + 66m² − 12m + 1          [= (m−1)¹²]
    P¹: −1024 − 391168m − 5329920m² − 10408960m³ + 5120m⁴ + 19456m⁵
        − 744448m⁶ + 99328m⁷ − 24576m⁸ − 1024m⁹
    P²: 393216 − 195100672m + 919994368m² − 192806912m³
        + 107347968m⁴ − 19464192m⁵ + 393216m⁶
    P³: −67108864 − 6643777536m − 3690987520m² − 67108864m³
    P⁴: 4294967296

Structure: F(m,0) = (m−1)¹², F(0,P) = (256P−1)⁴, and

    F = (m−1)¹² + (256P−1)⁴ − 1 − 1024·P·m·S(m,P),

    S(m,P) = 65536·P²m² + 3604480·P²m + 6488064·P²
           − 384·Pm⁵ + 19008·Pm⁴ − 104832·Pm³ + 188288·Pm²
           − 898432·Pm + 190528·P
           + m⁸ + 24m⁷ − 97m⁶ + 727m⁵ − 19m⁴ − 5m³
           + 10165m² + 5205m + 382.

F is irreducible over ℚ. The relation first appears at the monomial box
(12,4), where the nullspace is 1-dimensional.

### Proof (Sturm bound)

K = η(2z)²η(34z)² has ord 3 at every cusp of Γ₀(68). K³² has weight 64
and ord_{1/34}(K³²) = 96, equal to the maximum pole order of m^i P^j over
the support (96, attained uniformly on the triangle boundary i + 3j = 12).
All 35 terms K³²·m^i·P^j are holomorphic weight-64 forms with trivial
character on Γ₀(68) (Newman, character, and Ligozat checks all pass).
Sturm bound ⌈64·108/12⌉ = 576: the coefficients q⁰…q⁵⁷⁶ of K³²·F are all
zero (in fact no nonzero coefficient to q⁸⁰⁰). Hence K³²·F ≡ 0, so F = 0.

### Verification

- Exact: F(m,P) = 0 as a formal q-series to q⁶⁰⁰.
- Numerical (50 digits): relative error ~10⁻⁵¹ at 7 points; flagship
  m = √17, P ≈ 0.00017769354832319060602.
- Cusp values: F(1,0) = 0 (cusp ∞: m=1, P=0), F(17,1) = 0 (cusp 0: m=17,
  P=1, i.e. p_prod = 256).

### The curve F(m,P) = 0 has genus 3

Plane curve of degree 12, arithmetic genus 55. One point at infinity,
unibranch, with multiplicity sequence 8, 4, 4 under iterated blow-up,
δ = 28 + 6 + 6 = 40. Twelve affine nodes (δ = 1 each); their
m-coordinates are the 12 roots of

    S(m) = 15625m¹² + 57350m¹¹ + 435636m¹⁰ + 4017474m⁹ − 7136763m⁸
         + 58172316m⁷ − 2075760m⁶ + 39725924m⁵ + 142214563m⁴
         + 63554366m³ + 7265596m² − 57222m − 4913

(each verified a node by a nonzero tangent-cone discriminant
F_mm·F_PP − F_mP² at 80 digits). Genus = 55 − 40 − 12 = 3, confirmed by
Riemann–Hurwitz on m (R = 12: 3 at m = 0, 6 at the six simple ramification
points, 3 at m = ∞; the nodes contribute 0) and by the degree-2 map X₀(68) → C
(g(X₀(68)) = 7, R′ = 4). Discriminant:
disc_P(F) = −1267650600228229401496703205376·m³·(16m⁶−49m⁵+372m⁴+1338m³
  +6324m²−14161m+78608)·S(m)².
The fibre over m = 0 is the single point (0, 1/256), smooth, with e = 4
(F_m(0,1/256) = −4913 = −17³ ≠ 0).

Scripts: `p17/step0_explore.py`, `p17/step1_search.py`, `p17/step2_exact.py`,
`p17/step3_sturm.py`, `p17/step4_structure.py` (resultants, discriminant),
`p17/step4b_blowup.py` (δ at infinity by iterated blow-up),
`p17/step4c_affine.py` (affine nodes), `p17/SUMMARY.md`.

### Note on Ramanujan's degree-17 equations (Entry 12)

Ramanujan already recorded degree-17 modular equations (Notebook 2,
Entry 12) — but in the classical (α, β) variables, not in the theta
variables (m, p_prod). The relation above is the first in Ramanujan's
theta variables, and it is new. (17 is in Ramanujan's "quadratics only"
list: (17−1)/2 = 8 = 2³.)

---

## 3. The emerging general pattern — and why "any prime p" is not yet proved

| p | level | index | relation | support | weight | Sturm bound | curve genus |
|---|-------|-------|----------|---------|--------|-------------|-------------|
| 11 | 44 | 72 | G(Y,P), deg 5 | Y^iP^j, i+j ≤ 5 (21 mon.) | 50 | 300 | 2 |
| 13 | 52 | 84 | F(m,P), deg 7 | 17 mon. | 36 | 252 | 2 |
| 17 | 68 | 108 | F(m,P), deg 12 | i+3j ≤ 12 (35 mon.) | 64 | 576 | 3 |

Common features: pure parts are perfect powers with exponents equal to the
degrees (deg_m F and deg_P F); the mixed remainder is divisible by P·m;
the map X₀(4p) → C has degree 2 in all three cases (index
[ℂ(X₀(4p)):ℚ(m,P)] = 2); the genus is 2, 2, 3.

But a uniform theorem for all primes p is **not** established:
- The relation is in (m,P) for p=13, 17 but in (Y,P) = (m²,P) for p=11
  (only even powers of m occur at p=11); there is no proved rule for
  when odd powers of m appear.
- The support pattern differs: triangular i+j ≤ 5 (p=11), an irregular
  17-monomial support (p=13), triangular i+3j ≤ 12 (p=17).
- The genus sequence 2, 2, 3 has no proved formula.
- The genus-1 obstruction on X₀(p) for p ≥ 11 (the reason the p=11
  relation could not be a rational formula in (m, p_prod) on X₀(11))
  persists: the functions live on X₀(4p), and each prime still needs its
  own computation — the Sturm bound grows (300, 252, 576) and the linear
  algebra grows with the support, so the method is effective prime by
  prime but is not a uniform proof.
- The auxiliary polynomials r₁₃, r₁₇ (degrees 6 and 8) and the root
  formulas remain open, exactly as r₁₁ is open. The modular equation is
  Layer 3 of the uniform template; Layers 4–5 are untouched for p=13, 17.

Next primes by the same method: p=19 (level 76) and p=23 (level 92).
The Sturm bounds will be larger (weight and index both grow), but nothing
in the p=13/17 computations suggests a hard obstruction — the bottleneck
is the size of the exact q-series linear algebra and the Sturm check,
both of which scaled comfortably (seconds) at p=17.
