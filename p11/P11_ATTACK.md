# Phase 2, Step 1 — The degree-11 modular equation F(m, p_prod) = 0

**Date:** 2026-11-14. **Tools:** Python 3, sympy 1.12, mpmath 1.2.1, numpy 1.26.4 (no SageMath / PARI / Mathematica).
**Scripts:** `~/workspace/theta-phase1/p11/` (`step0_explore.py`, `step1_search.py`, `step2_exact.py`,
`step3_sturm.py`, `step4_structure.py`, `step4b_singular.py`, `step4c_genus.py`, `step4d_complex_sing.py`,
`step5_final.py`). Builds on `~/workspace/theta-phase1/theta.py`, `e4_modulareq.py`, and
`~/workspace/modular-forms/scripts/{eta,sturm}.py`.

## Headline result — the relation is FOUND and PROVED by Sturm bound

With `m = φ²(q)/φ²(q¹¹)`, `p = p_prod = 32q⁵χ(q)/χ¹¹(q¹¹)`, `Y = m²`, `P = p/32`,

```
G(Y, P) = 0,  where
G(Y,P) = −Y⁵ + 5Y⁴ − 10Y³ + 10Y² − 5Y + 1
       + P·(160Y⁴ + 2432Y³ + 19392Y² + 10624Y + 160)
       + P²·(−10240Y³ − 55296Y² − 796672Y + 10240)
       + P³·(327680Y² + 6914048Y + 327680)
       + P⁴·(−5242880Y + 5242880)
       + 33554432·P⁵.
```

Equivalently, in the task's variables `(m, p)`:

```
F(m,p) = −m¹⁰ + 5m⁸p + 5m⁸ − 10m⁶p² + 76m⁶p − 10m⁶
       + 10m⁴p³ − 54m⁴p² + 606m⁴p + 10m⁴
       − 5m²p⁴ + 211m²p³ − 778m²p² + 332m²p − 5m²
       + p⁵ + 5p⁴ + 10p³ + 10p² + 5p + 1  =  0.
```

Structural form (exact): `G(Y,P) = (1−Y)⁵ + (1+32P)⁵ − 1 + P·Y·S(Y,P)` with

```
S(Y,P) = −32·(163840P³ − 10240P²Y − 216064P² + 320PY² + 1728PY + 24896P
              − 5Y³ − 76Y² − 606Y − 332),
```

i.e. in `(m,p)`: `F(m,p) = (1+p)⁵ − (m²−1)⁵ + (p/32)·m²·S(m², p/32)`.
The pure parts are perfect fifth powers: the `P`-only part is `(1+32P)⁵ = (1+p)⁵`
(coefficients `1, 160, 10240, 327680, 5242880, 33554432 = C(5,k)·2^{5k}`) and the
`Y`-only part is `(1−Y)⁵`. Only **even** powers of `m` occur: the relation is really
in `Y = m² = φ⁴(q)/φ⁴(q¹¹)`, the exact analogue of the p = 7 variable
`X = φ⁴(q)/φ⁴(q⁷)`. Total degree 5 in `(Y,P)` jointly (`Y^i P^j` with `i+j ≤ 5` only);
`G` is irreducible over ℚ (sympy `factor_list`).

## 1. Setup: where m and p live

Verified exactly (series to `q^200`, `step0_explore.py`):

- `φ(q) = η(2z)⁵/(η(z)²η(4z)²)`, `χ(q) = q^{1/24}η(2z)²/(η(z)η(4z))`, with `q = e^{2πiz}`.
- `m` is the eta quotient `{1:−4, 2:10, 4:−4, 11:4, 22:−10, 44:4}`, shift 0.
- `p` is the eta quotient `32·{1:−1, 2:2, 4:−1, 11:11, 22:−22, 44:11}`, shift 5
  (i.e. `p = 32q⁵·(integer series)`; `p/32 = q⁵·(1 + q + q³ + …)`).

Both satisfy the Newman congruences at level **44** with trivial character
(`D_eff(m) = 1/121`, `D_eff(p) = 1`, both squares), so `m, p ∈ ℂ(X₀(44))` — level
`4·11`, exactly as in the p = 7 case (level 28), **not** level 11. The genus-1
obstruction of the Phase-1 report is therefore not quite the right picture for
`(m, p)`: see §4.

Ligozat cusp orders on `Γ₀(44)` (cusps `∞`, `0`, `1/2`, `1/4`, `1/11`, `1/22`;
note `1/44 ∼ ∞` and `1/1 ∼ 0`, so the table's seven entries cover six cusps):

| cusp | ord(m) | ord(p) |
|---|---|---|
| ∞ | 0 (value 1) | 5 (zero) |
| 0 | 0 (value 11) | 0 (value 32) |
| 1/2 | 5 (zero) | 0 |
| 1/4 | 0 | 0 |
| 1/11 | 0 | 5 (zero) |
| 1/22 | −5 (pole) | −10 (pole) |

Hence `deg(m) = 5`, `deg(p) = 10`, and both functions have their **only** pole at
the cusp `1/22`. Degree bounds for the minimal relation: `deg_m F ≤ deg(p) = 10`,
`deg_p F ≤ deg(m) = 5`. Cusp values check out: `F(1, 0) = 0`, `F(11, 32) = 0`
(equivalently `G(1,0) = G(121,1) = 0`).

## 2. How the relation was found

1. **Mod-prime detection** (`step1_search.py`): exact integer series for `m` and
   `P = p/32` to `q^400`; monomials `m^i P^j` built by numpy convolution mod
   `p = 999983` (float64 is exact here: partial sums `< 6·10¹⁴ < 2⁵³`); nullspace
   mod `p` over boxes `(I,J)` increasing. No relation for any box below
   `(10,5)`; at `(10,5)` (66 monomials, 106 rows) the nullspace is 1-dimensional,
   supported only on even powers of `m` with the triangular pattern
   `Y^i P^j`, `i+j ≤ 5`, and `P⁵`-coefficient 1 mod `p`.
2. **Exact recovery** (`step2_exact.py`): sympy rational nullspace on the
   21 triangle monomials with exact integer series (90 rows) gives a unique
   primitive integer vector — the coefficients quoted above (largest
   `33554432 = 2²⁵`; all coefficients are small).
3. **Verification** (three independent routes):
   - **Exact series:** `G(Y,P) = 0` as a formal power series to `q^600`
     (exact integer arithmetic, `theta.py` series machinery).
   - **Numerics, 50 digits:** relative error `≤ 7×10⁻⁵⁰` at `q = 0.2, …, 0.85`
     (8 points); at the flagship `q₁ = e^{−π/√11}` (where `m = √11`,
     `p = 4√2/G₁₁¹⁰ = 0.41643943427586409568…`) `|F| = 1.6×10⁻⁴⁵`
     (relative `~10⁻⁵⁰`); at three complex `q` (40 digits) relative error
     `~10⁻⁴⁰`.
   - **Sturm-bound proof** (§3) — this is a proof, not just evidence.

Approaches 2 (Fricke quotient) and 3 (Weber variable) from the task were **not
needed**: the direct attack succeeded on the first attempt at the maximal
degree box. No PSLQ/LLL was needed either — exact linear algebra sufficed,
because the coefficients are small (contra the "large coefficients" worry in
the attack plan; the Bailey–Borwein–Borwein remark about three-digit
coefficients is consistent with what was found).

## 3. Sturm-bound proof

On `Γ₀(44)`, index `μ = 72`. Let `K = η(2z)²η(22z)²` (weight 2, trivial
character). Then `K²⁵` has weight 50 and `ord_{1/22}(K²⁵) = 50`, while
`Y^i P^j` has its only pole at `1/22`, of order `10(i+j) ≤ 50`. Hence every
term `K²⁵·Y^i·P^j` (21 terms) is a **holomorphic modular form of weight 50,
trivial character, on `Γ₀(44)`** — each term's Newman congruences, `D_eff`
square, and Ligozat orders `≥ 0` at every cusp were checked symbolically
(`step3_sturm.py`, all 21 pass). The Sturm bound is `⌈50·72/12⌉ = 300`.
The `q`-expansion of `K²⁵·G(Y,P)` (exact integer arithmetic: `theta.py` series
for `Y, P`, `eta.py` for `K`, cross-checked `K²⁵` against its eta quotient)
has **all coefficients `q⁰…q³⁰⁰` equal to zero** (in fact no nonzero coefficient
to `q⁴⁰⁰`). Therefore `K²⁵·G(Y,P) ≡ 0` as a modular form, and since
`K ≢ 0`, **`G(Y,P) = 0` identically**. This is a complete proof in the same
style as the E4 proofs of the p = 3, 5, 7 relations.

## 4. Structure of the curve G(Y,P) = 0 (genus 2, not 1)

- **Irreducible** over ℚ; plane quintic, degree 5 in each variable.
- **Point at infinity:** the degree-5 homogeneous part is `−(Y−32P)⁵`, so the
  projective curve meets `Z = 0` only at `[32:1:0]`; that point is **smooth**
  (in the chart `P = 1` the linear term is `z·254803968 ≠ 0`) — the line at
  infinity is tangent with contact order 5.
- **Singular points (affine):** exactly four, all **ordinary double points**
  (nodes, δ = 1 each; tangent slopes are complex-conjugate pairs, so they are
  isolated real points). Refined by Newton on `(G_Y, G_P) = 0` to 25–30 digits
  (residuals `≤ 10⁻⁵⁵`):
  - `(Y, P) ≈ (−2.93988499066053333581914968158, 0.26961640660750642547444000427)`
  - `(Y, P) ≈ (−1.96224624306452032731430997765, −0.0115162811260059299781515215795)`
  - `(Y, P) ≈ (0.297589714028020910873560077037, 0.0419808009543044510091461004548)`
  - `(Y, P) ≈ (−0.1444296737186051078224049, 0.004060020066252666663289697)`
  Their `Y`-coordinates are the four roots of
  `59049Y⁴ + 280422Y³ + 293768Y² − 64614Y − 14641`
  (`59049 = 3¹⁰`, `14641 = 11⁴`), obtained as `gcd(res_P(G,G_P), res_P(G,G_Y))`
  (the gcd is that quartic squared). The fourth point was initially missed: its
  `P`-coordinate is a double root of `G(Y,·)`, which `mp.polyroots` splits into
  a spurious complex pair; numpy roots plus the resultant-gcd argument locate it.
- **Genus:** `(5−1)(5−2)/2 − Σδ = 6 − 4 = 2`. So the `(m,p)`-relation is a
  **genus-2** plane curve, *not* genus 1. Correction to the Phase-1/Phase-2
  scouting expectation ("plane model of an elliptic curve"): that was stated
  for `X₀(11)`; the functions `m, p` actually live on `X₀(44)`, and their image
  curve has genus 2. Consistency check (Riemann–Hurwitz): the map
  `X₀(44) → C` has degree 2 (tower `[ℂ(X₀(44)):ℚ(P)] = 10`,
  `[ℚ(Y,P):ℚ(P)] = 5`), `g(X₀(44)) = 4`, `g(C) = 2`, so
  `2·4 − 2 = 2(2·2 − 2) + R` gives `R = 2` ramification points — consistent.
  (An earlier genus-3 count, from only three nodes, violated Riemann–Hurwitz
  and is what forced the discovery of the fourth node.)
- **The degree-2 map is the `w₄` quotient:** the genuine Atkin–Lehner
  involution `w₄` of `X₀(44)` (representative `(1/2)[[12,1],[44,4]]`,
  `τ ↦ (6τ+1/2)/(22τ+2)`) **fixes both `Y` and `P`** (verified numerically to
  40 digits at two `τ`), so `ℚ(Y,P) = ℂ(X₀(44))^{w₄}` and `C ≅ X₀(44)/w₄`.
  The field extension is `ℚ(m,P)/ℚ(Y,P)` with Galois involution `m ↦ −m`.
- **`w₁₁` / `w₄₄` action:** `w₄₄` (Fricke, `τ ↦ −1/(44τ)`) sends `Y ↦ 121/Y`
  (verified to 15 digits at three points; `w₁₁` acts identically on `(Y,P)`
  since it differs from `w₄₄` by `w₄`). It preserves the curve:
  `G(121/Y, P') = 0` to `~10⁻³⁰` absolute (`~10⁻⁴⁸` relative) where
  `P' = P(w₄₄τ)`; `P'` is *not* `1−P` — it is determined implicitly by `G`
  (no simple closed form found). Note: the naive maps `τ ↦ −1/(4τ)` and
  `τ ↦ −1/(11τ)` do **not** normalise `Γ₀(44)` and are not involutions of
  `X₀(44)`; an early numeric test using them gave a spurious near-identity
  `Y(1/u) ≈ u²Y(u)` (agreement only to `~10⁻¹⁰`, an artefact of
  `φ(e^{−11πu}) ≈ 1` in that range) — discarded, correct representatives used.

## 5. What remains (for Step 2 and beyond)

1. **Auxiliary quintic `r₁₁`** (Step 2): with `G(Y,P) = 0` in hand, the
   coefficients of `r₁₁(ξ)` are to be found as functions of `(Y,P)` — test
   (H1) vs (H2) (§1.5 of the generalisation report) using the power sums of
   `u_k¹¹/p` and cyclic sums, now with the exact relation available to reduce
   powers of `Y`. The `w₄`-quotient picture says the natural coefficient field
   is `ℚ(Y,P)`, the function field of the genus-2 curve `C`.
2. **Root formula / ordering** (Step 3): the exponent pair `(a,b)` in
   `u_k = (α_k^a p/α_{k+1}^b)^{1/11}` and the cyclic ordering, by the p = 7
   method (monotonicity + a cyclic-sum identity).
3. **Flagship / `G₁₂₁`** (Step 4): at `q₁`, `Y = 11`, `P = p/32` with
   `p = 4√2/G₁₁¹⁰`; `r₁₁` at that point should be a quintic over `ℚ(G₁₁)`.
   `G₁₂₁`'s minimal polynomial is still open (Phase 1).
4. **Explicit `w₁₁`-action on `P`:** find a closed form for `P*` in
   `(Y,P)` (it lies in `ℚ(Y,P)`); relatedly, a hyperelliptic model of `C`
   (genus 2) would make the arithmetic of the curve more transparent.
5. **Fricke-quotient / Weber routes:** not needed for the relation itself;
   the `w₄`-quotient above is the precise sense in which a Fricke-type
   quotient enters. The degree-11 Schröter proof (Γ(2) infrastructure) and
   degree-17 entries remain open from Phase 1.

## 6. Files

| file | contents |
|---|---|
| `p11/step0_explore.py` | eta-quotient identities (exact), `m`/`p` series, Ligozat orders, cusp-0 limits |
| `p11/step1_search.py` | mod-prime linear-algebra search over monomial boxes; finds the triangle at (10,5) |
| `p11/step2_exact.py` | exact integer coefficients (sympy nullspace), exact check to `q⁶⁰⁰`, 50-digit numerics |
| `p11/step3_sturm.py` | **Sturm proof** (weight 50, bound 300) + structural forms |
| `p11/step4_structure.py` | irreducibility, `w₄₄` numeric action, resultant-gcd quartic |
| `p11/step4b_singular.py` | numeric singular-point scan (found 3; 4th hidden by a double root) |
| `p11/step4c_genus.py` | Newton refinement of 3 nodes, multiplicity/tangent cones |
| `p11/step4d_complex_sing.py` | complex-`P` scan locating the 4th node |
| `p11/step5_final.py` | 4th-node refinement, genuine `w₄`/`w₁₁` action, complex-`q` verification |
| `p11/P11_ATTACK.md` | this report |
