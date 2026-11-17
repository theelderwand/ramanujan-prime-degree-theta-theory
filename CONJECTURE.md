# A Conjecture for the Degree‑p Modular Equation in Ramanujan's Prime‑Degree Theta Theory

**Status:** conjecture, formulated 2026‑11‑14 from the six proved cases p = 3, 5, 7, 11, 13, 17.
**Source cases:** `paper/paper.tex` (Theorems 5.1–5.3 for p = 11, 13, 17; §4 for p = 3, 5, 7),
`PHASE1_REPORT.md`, `P13_P17_REPORT.md`. New computations for this document:
`/tmp/w4_test.py` (Atkin–Lehner w₄ action), `/tmp/cusp_test.py` (cusp values, support counts).

**How to read the labels.** Every claim carries one of:

- **[Theorem]** — proved in general (all odd primes p), proof sketched or cited.
- **[Proved for p ≤ 17]** — proved for the six cases (by Sturm bound / exact computation).
- **[Verified numerically]** — checked numerically to high precision for p ≤ 17 (and, where
  stated, for p = 19, 23); no proof written yet.
- **[Pattern from data]** — holds in all cases where it can be checked, but no proof and no
  mechanism is known; it is a pattern, not a theorem.
- **[Pure guess]** — an extrapolation with thin or no evidence; stated precisely so it can be
  falsified.

Nothing labelled [Pattern] or [Pure guess] is presented as a theorem. The document's job is to
say exactly what the data forces, what it suggests, and what it does not decide.

---

## 1. The data

Let p be an odd prime. Define the two eta‑quotients on Γ₀(4p)

- m = φ(q)²/φ(qᵖ)² = η(z)⁻⁴ η(2z)¹⁰ η(4z)⁻⁴ η(pz)⁴ η(2pz)⁻¹⁰ η(4pz)⁴,
- P = p_prod / 2^((p−1)/2), where p_prod = 2^((p−1)/2) q^((p²−1)/24) χ(q)/χ(qᵖ)ᵖ,
  so P = η(z)⁻¹ η(2z)² η(4z)⁻¹ η(pz)ᵖ η(2pz)^{−2p} η(4pz)ᵖ
  (an eta‑quotient of q‑order (p²−1)/24 at ∞; the formula p_prod = 2^((p−1)/2)·(that
  eta‑quotient) is proved for all odd p in Phase 1, E1, via the Jacobi triple product).

The proved relations F(X, P) = 0, with X = m or Y = m²:

| p | X | relation | e_X = deg_X F | e_P = deg_P F | monomials | support | genus of C | g(X₀(4p)) | level / index |
|---|---|----------|---------------|---------------|-----------|---------|-----------|-----------|---------------|
| 3 | Y = m² | Y − 1 − p³ = 0 (p = 2P) | 1 | 1 | 3 | i+j ≤ 1 | 0 | 0 | Γ₀(12) |
| 5 | m | m − 1 − p = 0 (p = 4P) | 1 | 1 | 3 | i+j ≤ 1 | 0 | 1 | Γ₀(20) |
| 7 | Y = m² | Y² − (2+5p)Y + (1−p)³ = 0 (p = 8P) | 2 | 3 | 7 | 3i+2j ≤ 6 | 0 | 2 | Γ₀(28) |
| 11 | Y = m² | G(Y,P) = 0 | 5 | 5 | 21 | i+j ≤ 5 | 2 | 4 | Γ₀(44) / 72 |
| 13 | m | F(m,P) = 0 | 7 | 3 | 17 | 3i+7j ≤ 21 | 2 | 5 | Γ₀(52) / 84 |
| 17 | m | F(m,P) = 0 | 12 | 4 | 35 | i+3j ≤ 12 | 3 | 7 | Γ₀(68) / 108 |

Here e_X = deg_X F is the degree in X, e_P = deg_P F the degree in P. "Support" is the set
of (i, j) with a nonzero coefficient of XⁱPʲ. In every case **all** monomials of the stated
triangle occur (the counts above are exactly the triangle counts).

Pure parts (normalised so F(0, 0) = 1):

| p | F(X, 0) | F(0, P) |
|---|---------|---------|
| 3 | 1 − Y | 1 + p³ = 1 + 8P³ |
| 5 | 1 − m | 1 + p = 1 + 4P |
| 7 | (1 − Y)² | (1 − p)³ = (1 − 8P)³ |
| 11 | (1 − Y)⁵ | (1 + 32P)⁵ |
| 13 | (1 − m)⁷ | (1 + 64P)³ |
| 17 | (1 − m)¹² | (1 − 256P)⁴ |

---

## 2. The conjectural template for general prime p

### 2.1 The functions and their degrees — [Theorem]

For every odd prime p, m and P are modular functions on Γ₀(4p) (Newman + character checks),
with cusp orders (Ligozat, computed in general — see §2.1 proof sketch) at the six cusps
∞, 0, 1/2, 1/4, 1/p, 1/(2p):

| cusp | ord(m) | ord(P) |
|------|--------|--------|
| ∞ | 0 | (p²−1)/24 |
| 0 | 0 | 0 |
| 1/2 | (p−1)/2 | 0 |
| 1/4 | 0 | 0 |
| 1/p | 0 | (p²−1)/24 |
| 1/(2p) | −(p−1)/2 | −(p²−1)/12 |

Consequences: the only pole of m is at 1/(2p), of order (p−1)/2, so **deg(m) = (p−1)/2**;
the only pole of P is at 1/(2p), of order (p²−1)/12, so **deg(P) = (p²−1)/12**; P's only
zeros are at ∞ and 1/p, each of order (p²−1)/24. The genus of X₀(4p) is **(p−3)/2**
(index μ = 6(p+1), 6 cusps, no elliptic points since 4 | 4p and 2 | 4p).

*Proof sketch.* The Ligozat order at a cusp 1/c is (N/(24·gcd(c², N)))·Σ_e e·gcd(c, d)²/d
for N = 4p. Each entry above is a short exact computation in the edicts (done by hand for
general p in the working notes; the p = 11 instance is in `p11/step0_explore.py`). The genus
formula is the standard μ/12 − c/2 + 1 with μ = 6(p+1), c = 6, e₂ = e₃ = 0. ∎

### 2.2 The Atkin–Lehner involution w₄ — [Verified numerically; proof outstanding]

Let w₄ be the Atkin–Lehner involution of X₀(4p) attached to Q = 4, represented by
W₄ = [[4x, y], [4pz, 4w]] with 4xw − pyz = 1. Then

> **w₄(m) = (−1/p)·m,  w₄(P) = P**   (for p ≥ 5; for p = 3, w₄(m) = −m but w₄(P) = ζ₃·P).

*Evidence.* Computed numerically (eta via pentagonal series, 40 digits) at τ = i, 2i,
0.3+1.5i for every p ∈ {3, 5, 7, 11, 13, 17}: the ratio m(w₄τ)/m(τ) equals (−1/p) and
P(w₄τ)/P(τ) equals 1, to within 10⁻³⁰ (p = 3: P-ratio = e^{2πi/3}). A proof should go
through the eta‑multiplier (Dedekind sums) for W₄; it has not been written. This is the
single most important missing proof in the template, because (2.3)–(2.5) all follow from it.

*Remark (p = 3 degeneracy).* For p = 3 the normalised P is not w₄‑invariant; only P³ is.
That is why the p = 3 relation is in (Y, p³), not (Y, P), and why the formulae
e_X = (p²−1)/24 and deg(P) = (p²−1)/12 fail at p = 3 ((p²−1)/24 = 1/3 is not an integer).
p = 3 is degenerate in every formula that follows; it is excluded from the general
statements and listed separately.

### 2.3 Which variables: (m, P) or (Y, P) — [Follows from 2.2]

Since w₄ fixes P and sends m to (−1/p)·m:

- if **p ≡ 1 (mod 4)**, then (−1/p) = +1, w₄ fixes m, and the relation is in **(m, P)**;
- if **p ≡ 3 (mod 4)**, then (−1/p) = −1, w₄ sends m ↦ −m, so only **Y = m²** is
  w₄‑invariant, and the relation is in **(Y, P)**.

This fits all six cases (p = 3, 7, 11 → Y; p = 5, 13, 17 → m). The rule is p mod 4, not
p mod 8: the mechanism is the Legendre symbol (−1/p), which depends only on p mod 4.
[Verified numerically for p ≤ 17 via 2.2.]

### 2.4 The degree‑2 map — [Proved for p ≤ 17; conjectured in general]

> The map X₀(4p) → C, where C is the plane curve F(X, P) = 0, has degree 2, and
> C ≅ X₀(4p)/w₄.

*Why for p ≤ 17.* [Q(X₀(4p)) : Q(P)] = deg(P) = (p²−1)/12, and [Q(X, P) : Q(P)] = e_X
(the degree of F in X), so the map has degree ((p²−1)/12)/e_X = 2 in every proved case
(using e_X = (p²−1)/24, §2.5). Since Q(X, P) ⊆ Fix(w₄) and both are index‑2 subfields,
Q(X, P) = Fix(w₄). [Conjectured for general p: the identity e_X = (p²−1)/24 is the only
input, and it is a theorem (§2.5), so the degree‑2 statement is one step from proved —
it needs only the existence of a relation F with deg_X F = (p²−1)/24, which is what the
Sturm computation produces.]

### 2.5 The pure‑part exponents — [Theorem, given 2.2 and 2.4]

> **e_X = (p²−1)/24**  (p ≥ 5), and
> **e_P = (p−1)/4 if p ≡ 1 (mod 4) [X = m];  e_P = (p−1)/2 if p ≡ 3 (mod 4) [X = Y].**

*Proof.* deg_C(P) = deg_{X₀(4p)}(P)/2 = (p²−1)/24, and deg_C(P) = [Q(X,P) : Q(P)] = e_X
(the degree of the minimal polynomial of X over Q(P)). Similarly deg_C(X) =
deg_{X₀}(X)/2 = e_P: for X = m, deg_{X₀}(m)/2 = (p−1)/4; for X = Y = m²,
deg_{X₀}(Y)/2 = (p−1)/2. ∎

Note the p = 13 "anomaly" (e_X = 7 > e_P = 3) is not an anomaly at all: e_X = (p²−1)/24
grows like p²/24 while e_P grows like p/4 or p/2, so e_X > e_P for all p ≥ 11, and the
crossover is exactly at p = 11 (e_X = e_P = 5). The exponents are **not** (p−1)/2 or
(p+1)/2 in any simple way; the correct pair is ((p²−1)/24, (p−1)/4 or (p−1)/2).

### 2.6 The pure parts — [Pattern from data, with a proved mechanism for the X‑part]

Normalised by F(0, 0) = 1, the template is

> **F(X, 0) = (1 − X)^{e_X},   F(0, P) = (1 − (2/p)·2^{(p−1)/2}·P)^{e_P},   F(0, 0) = 1.**

i.e. the full structure is

> **F(X, P) = (1 − X)^{e_X} + (1 − (2/p)·2^{(p−1)/2}·P)^{e_P} − 1 + P·X·S(X, P)**

for some polynomial S. Check against the table in §1: p = 11 gives (1−Y)⁵ + (1+32P)⁵ − 1
((2/11) = −1, so −(2/11) = +1); p = 13 gives (1−m)⁷ + (1+64P)³ − 1 ((2/13) = −1);
p = 17 gives (1−m)¹² + (1−256P)⁴ − 1 ((2/17) = +1, so −(2/17) = −1). All six cases fit,
including the genus‑0 cases (p = 7: (1−Y)² + (1−8P)³ − 1; p = 5: (1−m) + (1+4P) − 1;
p = 3: (1−Y) + (1+8P³) − 1 in the variable p³).

*Why the X‑part is a perfect power* [Theorem, given 2.2]: P's zeros on C are the single
point (X, P) = (1, 0) — the image of the two cusps ∞ and 1/p, which w₄ swaps, and at both
of which X = 1. (X(∞) = 1 since m(∞) = 1; X(1/p) = 1 because m(1/p) = (−1/p), verified
numerically below, and X is m or m² according as (−1/p) = +1 or −1.) Since deg_C(P) = e_X,
the fibre P = 0 is the single point X = 1 with multiplicity e_X, i.e. F(X, 0) = (1−X)^{e_X}
up to the normalisation F(0,0) = 1. Similarly the P‑part is a perfect power because X's only
zero on X₀(4p) is the cusp 1/2 (order (p−1)/2), which is fixed by w₄ (hence ramified in the
degree‑2 map), so X has a single zero on C, of order e_P, at the P‑value P(1/2).

*The sign (2/p)* [Pattern from data — 6/6, plus 2 new primes verified]: the P‑value at the
cusp 1/2 is P(1/2) = (2/p)/2^{(p−1)/2}, so the pure P‑part is (1 − (2/p)·2^{(p−1)/2}P)^{e_P}.
Verified numerically (cusp‑limit computation, `/tmp/cusp_test.py`) for p = 5, 7, 11, 13, 17
**and for p = 19, 23** (ratio to prediction 0.9999993, the deviation being the finite‑t
error, identical across all p). The sign is the Legendre symbol (2/p) = +1 iff p ≡ ±1
(mod 8), −1 iff p ≡ ±3 (mod 8). No proof is known; the mechanism should be the
eta‑quotient cusp‑value formula at 1/2. The cusp values m(1/p) = (−1/p) were likewise
verified numerically for p = 11, 13, 17, 19, 23 (to ~10⁻¹³; the computation needs 400
digits because of catastrophic cancellation in the pentagonal series at tiny Im τ).

### 2.7 The support — [Pattern from data — 6/6]

> **The support of F is the full triangle {(i, j) : i·e_P + j·e_X ≤ e_X·e_P},**
> i.e. the triangle with vertices (0,0), (e_X, 0), (0, e_P). Every monomial of the
> triangle occurs with nonzero coefficient.

Check: p = 11: i·5 + j·5 ≤ 25, i.e. i+j ≤ 5, 21 monomials ✓. p = 13: 3i+7j ≤ 21, 17 ✓.
p = 17: 4i+12j ≤ 48, i.e. i+3j ≤ 12, 35 ✓. p = 7: 3i+2j ≤ 6, 7 ✓. p = 5, 3: i+j ≤ 1, 3 ✓.
The number of monomials is N(e_X, e_P) = #{(i,j) : i·e_P + j·e_X ≤ e_X·e_P}.

No proof is known. It is a genericity statement (the relation uses the whole triangle), and
it is the most "accidental"‑looking part of the template — but 6/6, including three
genus‑0 cases where the polynomials are completely classical, is real evidence. Note the
triangle is exactly the Newton polygon forced by the pure parts plus the condition that
F(0,0) = 1: the pure parts pin the two legs (e_X, 0) and (0, e_P), and the hypotenuse
i·e_P + j·e_X = e_X·e_P is the line joining them. What is *not* forced is that every
interior monomial occurs — that is the conjectural part.

### 2.8 The mixed remainder — [Theorem, given 2.6]

> The mixed part F(X,P) − [(1−X)^{e_X} + (1 − (2/p)·2^{(p−1)/2}P)^{e_P} − 1] is
> divisible by P·X.

*Proof.* It vanishes when X = 0 (since F(0,P) is the pure P‑part) and when P = 0 (since
F(X,0) is the pure X‑part), hence is divisible by X and by P, hence by X·P. ∎
The nontrivial content is therefore entirely in the perfect‑power pure parts (§2.6), not
in the divisibility, which is trivial algebra once the pure parts are known.

### 2.9 Irreducibility — [Theorem, given 2.4–2.5]

F is irreducible over ℚ: deg_X F = e_X = [Q(X,P) : Q(P)], so F is the minimal polynomial
of X over Q(P). (Verified directly for p = 11, 13, 17.)

### 2.10 The genus — [Pure guess, and the weakest part of the conjecture]

Let R′ = #Fix(w₄) be the number of fixed points of w₄ on X₀(4p). By Riemann–Hurwitz for
the degree‑2 map X₀(4p) → C,

> **g(C) = (p − 1 − R′)/4**,   where g(X₀(4p)) = (p−3)/2.

Two of the fixed points are cusps: w₄ swaps ∞ ↔ 1/p and 0 ↔ 1/4, and fixes 1/2 and
1/(2p) (checked by direct computation of the cusp permutation). So R′ = 2 + N_CM(p),
where N_CM(p) is the number of fixed **CM points** of w₄. The data:

| p | 3 | 5 | 7 | 11 | 13 | 17 |
|---|---|---|---|----|----|----|
| R′ | 2 | 4 | 6 | 2 | 4 | 4 |
| N_CM(p) | 0 | 2 | 4 | 0 | 2 | 2 |
| g(C) | 0 | 0 | 0 | 2 | 2 | 3 |

**No formula for N_CM(p) was found.** It is not a function of p mod 4 (p = 7 and p = 11 are
both ≡ 3 mod 4 but have N_CM = 4 and 0), nor of p mod 8 (p = 3 and p = 11 are both
≡ 3 mod 8 with N_CM = 0; p = 5 and p = 13 are both ≡ 5 mod 8 with N_CM = 2; p = 7 is
≡ 7 mod 8 with N_CM = 4; p = 17 is ≡ 1 mod 8 with N_CM = 2 — a suggestive but
single‑data‑point‑per‑class pattern). The count should come from the CM theory of w₄
(class numbers of orders in ℚ(i) and the number of w₄‑admissible level structures), but
the derivation was not completed. Integrality forces R′ ≡ p−1 (mod 4), so for p ≡ 3
(mod 4), R′ ∈ {2, 6, 10, …} and g(C) ∈ {(p−3)/4, (p−7)/4, …}; for p ≡ 1 (mod 4),
R′ ∈ {4, 8, …} and g(C) ∈ {(p−5)/4, (p−9)/4, …}.

**Primary prediction (minimal R′):** R′ = 2 if p ≡ 3 (mod 4), R′ = 4 if p ≡ 1 (mod 4),
i.e. g(C) = (p−3)/4 for p ≡ 3 (mod 4) and g(C) = (p−5)/4 for p ≡ 1 (mod 4). This fits
p = 11, 13, 17 but **fails p = 7** (which has R′ = 6, g = 0 instead of the predicted
g = 1). The genus prediction is therefore the least reliable part of this document; treat
the values in §3 as the minimal‑R′ guess, with the true genus possibly smaller by
(R′ − R′_min)/4.

---

## 3. Predictions for p = 19 and p = 23

| | p = 19 | p = 23 |
|---|---|---|
| p mod 4 | 3 | 3 |
| variables X | Y = m² | Y = m² |
| e_X = (p²−1)/24 | 15 | 22 |
| e_P = (p−1)/2 | 9 | 11 |
| pure X‑part | (1 − Y)¹⁵ | (1 − Y)²² |
| (2/p) | −1 (19 ≡ 3 mod 8) | +1 (23 ≡ 7 mod 8) |
| pure P‑part | (1 + 512P)⁹ | (1 − 2048P)¹¹ |
| support | 3i + 5j ≤ 45 (i·9 + j·15 ≤ 135) | i + 2j ≤ 22 (i·11 + j·22 ≤ 242) |
| monomials | **82** | **144** |
| level / index | Γ₀(76) / 120 | Γ₀(92) / 144 |
| g(X₀(4p)) = (p−3)/2 | 8 | 10 |
| genus g(C), minimal‑R′ guess | **4** (R′ = 2) | **5** (R′ = 2) |
| genus g(C), if R′ = 6 | 3 | 4 |
| cusp values | Y(∞) = Y(1/19) = 1, Y(0) = 361, P(0) = 1, P(1/2) = −1/512 | Y(∞) = Y(1/23) = 1, Y(0) = 529, P(0) = 1, P(1/2) = +1/2048 |

Full predicted template (p = 19):

> F₁₉(Y, P) = (1 − Y)¹⁵ + (1 + 512P)⁹ − 1 + P·Y·S₁₉(Y, P),
> support 3i + 5j ≤ 45, 82 monomials, irreducible over ℚ, genus 4 (if R′ = 2).

Full predicted template (p = 23):

> F₂₃(Y, P) = (1 − Y)²² + (1 − 2048P)¹¹ − 1 + P·Y·S₂₃(Y, P),
> support i + 2j ≤ 22, 144 monomials, irreducible over ℚ, genus 5 (if R′ = 2).

Already verified for p = 19, 23 (before any polynomial is computed): the cusp values
m(1/p) = (−1/p) = −1 and P(1/2) = (2/p)/2^{(p−1)/2} (p = 19: −1/512; p = 23: +1/2048),
which are exactly the inputs the template needs for the pure parts.

---

## 4. What would falsify the conjecture

Each part has its own falsifier; the conjecture is designed to be killed piece by piece,
not defended as a package.

1. **w₄ action (§2.2).** If for some prime p ≥ 5, m(w₄τ)/m(τ) is not the constant (−1/p),
   or P(w₄τ)/P(τ) is not the constant 1 (computed numerically at a generic τ), then the
   (m,P)/(Y,P) rule and the degree‑2 identification both collapse. *Test: the same
   `/tmp/w4_test.py` computation at p = 19, 23 (takes seconds).*
2. **Variables (§2.3).** If the p = 19 or p = 23 relation needs odd powers of m (i.e. is
   not a polynomial in Y = m² and P), the p‑mod‑4 rule is wrong.
3. **Exponents (§2.5).** If e_X ≠ (p²−1)/24 or e_P ≠ (p−1)/2 (for p ≡ 3 mod 4) at
   p = 19 or 23, the degree computation is wrong — this would be a genuine surprise,
   since §2.5 is a theorem given §2.2 and §2.4.
4. **Pure parts (§2.6).** If F(X,0) is not the perfect power (1−X)^{e_X}, or F(0,P) is
   not (1 − (2/p)2^{(p−1)/2}P)^{e_P} — in particular if the P‑part sign is not −(2/p) —
   then the cusp‑value picture is wrong. *The sign is already verified at the cusp‑value
   level for p = 19, 23, so a failure here would mean the pure part is not a perfect
   power at all, i.e. P has a zero on C away from X = 1.*
5. **Support (§2.7).** If any monomial of the triangle i·e_P + j·e_X ≤ e_X·e_P is absent
   (coefficient zero), or any monomial outside the triangle occurs, the support rule is
   wrong. *This is the most likely part to fail: it is pure genericity with no mechanism.*
6. **Genus (§2.10).** If g(C) at p = 19 is not 4, or at p = 23 is not 5, the minimal‑R′
   guess is wrong (but the framework g(C) = (p−1−R′)/4 survives — only the value of R′
   changes). A genus *larger* than the minimal‑R′ prediction would falsify the framework
   itself (R′ ≥ 2 always, since the two cusps 1/2 and 1/(2p) are fixed).
7. **Irreducibility (§2.9).** If F₁₉ or F₂₃ factors over ℚ, the degree argument is wrong.

---

## 5. Where the pattern is thin, and what could break

- **Only three hard data points.** The genus‑0 cases (p = 3, 5, 7) constrain the exponents
  and the support, but they say nothing about the genus, and their curves are rational, so
  the "template" is really extrapolated from p = 11, 13, 17 alone. Any feature that first
  appears at p = 19 or 23 (a new kind of singularity, a non‑triangular support, a second
  w₄‑orbit of cusps) is invisible to us now.
- **The genus is genuinely unpredicted.** N_CM(p) has no formula; p = 7 already breaks the
  naive pattern. The genus predictions for p = 19, 23 are the minimal‑R′ values and could
  be too large by 1 or more. Computing N_CM(19) and N_CM(23) from the CM fixed‑point theory
  of w₄ is the single most valuable theoretical step before running the p = 19, 23
  computations — it would turn the genus from a guess into a prediction.
- **The support triangle has no mechanism.** It holds 6/6, but there is no a priori reason
  the relation could not miss an interior monomial (e.g. a coefficient vanishing by an
  accident of the linear algebra). The p = 19 case (82 monomials) is the first real test:
  a single zero coefficient inside the triangle kills §2.7.
- **The w₄ action is unverified as a proof.** It is numerically certain (10⁻³⁰ across
  three τ and six primes), but the eta‑multiplier proof has not been written. If the proof
  attempt fails, it would mean the numerical identity is subtler than it looks (e.g. the
  multiplier depends on the representative in a way the numerics did not detect — unlikely,
  but the proof is the check).
- **p = 3 warns that degeneracies exist.** The p = 3 case breaks the exponent formulae,
  the P‑invariance, and the degree computation, all at once, because (p²−1)/24 is not an
  integer. There could be further degeneracies at larger p that the six cases do not reveal
  (e.g. if (p−1)/4 or (p²−1)/24 interacts with the level structure in an unforeseen way —
  no such interaction is visible, but the p = 3 precedent says to watch for it).
- **No uniform theorem.** The template is a per‑prime statement. The equation size grows
  (monomials ~ e_X·e_P/2 ~ p³/192), the genus grows (~p/4), and there is no uniform bound
  or closed form for S(X, P). A proof "for any prime p" in the sense of a single uniform
  argument is not in sight; what is in sight is a uniform *method* (Sturm bound on
  Γ₀(4p)) plus this template as the predicted answer for each prime.

---

## 6. How to test it (the p = 19 and p = 23 computations)

Same method as p = 11, 13, 17 (`p11/`, `p13/`, `p17/` scripts):

1. **Step 0.** Eta‑quotient edicts for p = 19, 23 (same shape as §1); Ligozat cusp orders
   (should reproduce the table in §2.1); Newman + character checks; cusp values
   m(1/p) = −1, P(1/2) = −1/512 (p = 19) / +1/2048 (p = 23) — already verified here.
2. **Step 1–2.** Linear algebra over the predicted triangle (82 monomials for p = 19,
   144 for p = 23) to find the exact relation; check the pure parts, the support, and
   irreducibility against §3.
3. **Step 3.** Sturm proof on Γ₀(76) (index 120) / Γ₀(92) (index 144), with
   K = η(2z)²η(2pz)² raised to a power clearing the maximal pole order over the support
   (following the p = 11/13/17 pattern: K²⁵ weight 50 bound 300; K¹⁸ weight 36 bound 252;
   K³² weight 64 bound 576; bound = μ·k/12).
4. **Step 4.** Genus from the plane curve (singularities / delta invariants) and comparison
   with §3; the genus result determines R′ = p − 1 − 4·g(C) and hence N_CM(p), the first
   new data on the genus problem.

Estimated effort: the same as the p = 17 computation, scaled by the monomial count
(144 vs 35) and the Sturm bound; no new ideas are needed, only compute time.

---

## Appendix: corrections to the paper / reports found while writing this

1. **p = 13 report, fibre over m = 0.** `P13_P17_REPORT.md` states the fibre over m = 0
   is the point (0, 1/4). That point is not on the curve: F(0, 1/4) = 4913 ≠ 0. The
   correct fibre is the single point **(0, −1/64)**, since F(0, P) = (1+64P)³, and it is
   smooth: F_m(0, −1/64) = −169 = −13². (The value −1/64 is P(1/2), the image of the
   cusp 1/2, as in §2.6.)
2. **p = 11, notation.** The paper writes "p_prod = 32q⁵·η…" which reads as if q⁵ were a
   separate factor; in fact the eta‑quotient itself has q‑order 5, i.e. p_prod = 32 ×
   (eta‑quotient with shift 5), and P = p_prod/32 is exactly that eta‑quotient
   (`p11/step0_explore.py`: `p_ed = {1:−1, 2:2, 4:−1, 11:11, 22:−22, 44:11}`, shift 5).
   The mathematics is unaffected; only the notation is misleading.
3. **p = 11 cusp value.** The paper's "G(11,1) = 76026720" label (in the summary, since
   corrected) is a slip: the cusp‑0 point is (Y, P) = (121, 1) since m(0) = 11, Y = m² = 121.
   G(121, 1) = 0, as stated in the paper.
4. **p = 13 cusp table.** The table lists both "∞" and "1/52" with the same orders; 1/52 is
   Γ₀(52)‑equivalent to ∞, so the "1/52" row duplicates the "∞" row (the cusp 0 is the
   c = 1 row). Harmless redundancy, not an error in the proof (the Sturm argument uses the
   Ligozat orders, which are correct).
