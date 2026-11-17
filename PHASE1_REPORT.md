# Phase 1 Report — Ramanujan Prime-Degree Theta Theory: Computational Reproduction

**Date:** 2026-11-14. **Tools:** Python 3, sympy 1.12, mpmath 1.2.1 (NO SageMath / PARI / Mathematica on this machine). **Scripts:** `~/workspace/theta-phase1/` (`theta.py`, `e1_dissection.py`, `e2_flagships.py`, `e3_classinvariants.py`, `e4_modulareq.py`, `p11_scout.py`, `p11_scout2.py`, `p7_calibration.py`). Sources: the deep-research report `../ramanujan-prime-degree-theta-theory-20261114-0146/report.md`, Berndt–Rebák 2025 (arXiv:2312.15513), Rebák 2023 (arXiv:2112.11884), Berndt *Notebooks* Part III.

## Summary

| Item | Status |
|---|---|
| E1 dissection (exact series + 50-digit numerics), n = 3,5,7,9,11,13,17 | ✅ verified |
| E1 p_prod χ-formula, general odd n (analytic derivation + checks) | ✅ verified, **new**: proved for all odd n, checked at 11,13,17 |
| E2 three flagships (cubic, quintic, septic) | ✅ verified to 50 digits |
| E2 machine reproduction (p from G, aux polys, root ordering, Lemma 3.8) | ✅ verified |
| E3 class invariants G_3…G_729 (11 values incl. new G_99, G_729, G_125, G_625) | ✅ verified to 50 digits |
| E4 multiplier relations p=3,5,7 | ✅ **proved by Sturm bound** (Newman + character + Ligozat + coefficient check) |
| E4 degree-11 Schröter equation | ⚠️ verified numerically (50 digits, 4 points); full Sturm proof NOT done — needs Γ(2) infrastructure |
| Phase-2 scouting (p=11) | ⚠️ mixed: flagship data computed; PSLQ power-sum fishing **failed** (spurious relations); calibration on p=7 **succeeded**; m = √11 at flagship |

---

## E1. Universal dissection and p_prod

### What was done
Implemented the dissection in the variable t = q^{1/n} (so all exponents are integers):
φ(t)/φ(t^{n²}) = 1 + Σ_{k=1}^{(n−1)/2} 2t^{k²} f(t^{n(n+2k)}, t^{n(n−2k)})/φ(t^{n²}),
with exact integer-coefficient series (φ, f via direct sparse sums; χ via product; series division by inversion), plus mpmath numerics at 50 digits.

### Results
- **Exact series, order t¹²⁰:** LHS − RHS = 0 for n = 3, 5, 7, 9, 11, 13, 17. (Exact integer arithmetic — this is a proof to that order.)
- **Numerics, 50 digits:** |dissection error| ≤ 2×10⁻⁴⁹ at q = 0.3, 0.5, 0.7 for n = 3, 5, 7, 11, 13, 17.
- **p_prod:** Π_k u_k = 2^{(n−1)/2} q^{(n²−1)/24} χ(q)/χⁿ(qⁿ) verified exactly (series to t¹²⁰) and numerically (rel err ≤ 6×10⁻⁵⁰) for the same n.

### New: general p_prod formula (analytic derivation)
The report stated the χ-formula only for p = 3, 5, 7. It holds for **every odd n**, by the Jacobi triple product:
f(q^{n+2k}, q^{n−2k}) = (−q^{n+2k};q^{2n})∞(−q^{n−2k};q^{2n})∞(q^{2n};q^{2n})∞.
As k runs over 1..(n−1)/2, the exponents n±2k run over all odd residues mod 2n except n, so
Π_k (−q^{n+2k};q^{2n})∞(−q^{n−2k};q^{2n})∞ = χ(q)/χ(qⁿ),
and φ(qⁿ)^{(n−1)/2} = χ(qⁿ)^{n−1}(q^{2n};q^{2n})∞^{(n−1)/2} (triple product with a=b=qⁿ), cancelling the (q^{2n};q^{2n})∞ factors. With Σk² = n(n²−1)/24:
**p_prod = 2^{(n−1)/2} q^{(n²−1)/24} χ(q)/χⁿ(qⁿ) for all odd n.**
Corollary (G-form): substituting χ(q) = 2^{1/4}q^{1/24}G_n gives
**p_prod = 2^{(n−1)/4} G_n / G_{n²}^{n}** for q = e^{−π√n},
which reproduces c₃=√2, c₅=2, c₇=2√2 and gives c₁₁ = 2^{5/2} = 4√2. (This corrects an error in the first draft of E3, which had used the p=5 shape 2/G₁₁⁴ for p=11; the correct flagship value is p = 4√2/G₁₁¹⁰ = 0.416439434275864095678569298496, now agreeing with the χ-quotient to 50 digits.)

---

## E2. The three flagship evaluations (machine reproduction)

All three flagships verified as stated, 50 digits (relative errors 0 to 3×10⁻⁵¹):
- Cubic: φ(e^{−3π√3}) = 3^{−3/4}φ(e^{−π√3}){1+2^{1/3}}.
- Quintic: φ(e^{−5π√5}) = 5^{−3/4}φ(e^{−π√5}){1+(2(1+tan π/5)/(1−sin π/5))^{1/5}+(2(1−tan π/5)/(1+sin π/5))^{1/5}}.
- Septic: φ(e^{−7π√7}) = 7^{−3/4}φ(e^{−π√7}){1+(cos π/7/(2cos²2π/7))^{2/7}+(cos 2π/7/(2cos²3π/7))^{2/7}+(cos 3π/7/(2cos²π/7))^{2/7}}.

Machine reproduction at q₁ = e^{−π/√p} (the point where the dissection + theta transformation give the flagship):
- **Cubic:** p via χ-quotient = p via G-formula (√2/G₃²) = 2^{1/3}, rel err 2×10⁻⁵¹; u from dissection = p; φ⁴(q)/φ⁴(q³) = 3 = 1+p³, rel err 1×10⁻⁵⁰.
- **Quintic:** p = √5−1 (both routes), rel err 4×10⁻⁵¹; φ²(q)/φ²(q⁵) = √5 = 1+p; auxiliary quadratic ξ²−((p−1)²+7)ξ+p³ = 0 has roots α ≈ 6.7770610073, β ≈ 0.2786670827; (αp)^{1/5} and (βp)^{1/5} reproduce the two published tan/sin terms exactly (rel err ≤ 4×10⁻⁵¹).
- **Septic:** p = 1 (both routes, exact); X = φ⁴(q)/φ⁴(q⁷) = 7 (rel err 5×10⁻⁵¹); modular-equation residual −2×10⁻⁴⁹; r(ξ) = ξ³−6ξ²+5ξ−1 exactly; its roots match 1/(2cos kπ/7)² (rel err ≤ 2×10⁻⁵⁰); the cyclic assignment (α,β,γ) = roots for k = 3,2,1 gives u = (α²p/β)^{1/7} etc. matching the dissection's u,v,w (total rel err 2×10⁻⁵⁰); **Lemma 3.8 confirmed: Σ s_i⁷ = 41.0** (rel err 1×10⁻⁵⁰), where s_i are the published trig terms.
- **Ordering note (worth recording):** the published septic terms appear in the order (v,u,w) relative to the dissection's u,v,w — i.e. s₁ = v, s₂ = u, s₃ = w. The sum is unaffected; the cyclic assignment above is stated in dissection order. Also verified: u > v > w at q₁, as the theory requires.

---

## E3. Class invariants

Verified against the definition G_n = 2^{−1/4}q^{−1/24}χ(q), q = e^{−π√n}, 50 digits:

| n | closed form source | rel err |
|---|---|---|
| 3 | 2^{1/12} | 3×10⁻⁵¹ |
| 5 | ((1+√5)/2)^{1/4} | 2×10⁻⁵¹ |
| 7 | 2^{1/4} | 0 |
| 11 | [(3√33+17)^{1/3}−(3√33−17)^{1/3}+2]/(3·2^{1/4}) | 0 |
| 25 | (√5+1)/2 | 7×10⁻⁵¹ |
| 27 | 2^{1/12}(2^{1/3}−1)^{−1/3} | 3×10⁻⁵¹ |
| 99 | G₉₉³ = 2^{1/4}/3·[(7822+1362√33)^{1/3}+(7957+1383√33)^{1/3}+2√33+13] | 2×10⁻⁵¹ |
| 729 | G₇₂₉³ = (√2/2)·R^{1/3}·(3(√3+1)/(√3+1−R^{1/3})−2), R = ((2(√3+1))^{1/3}+1)/((2(√3−1))^{1/3}−1) | 2×10⁻⁴⁹ |
| 125 | Cor. 5.18 (tan/sin π/5 fifth roots) | 6×10⁻⁵¹ |
| 625 | Cor. 5.19 | 6×10⁻⁵¹ |
| 243 | Watson: 2^{1/12}(2^{1/3}+2^{2/3}+3^{1/3})/(9−2·3^{4/3})^{1/3} | 2×10⁻⁵¹ |

Reciprocity G_{1/n} = G_n confirmed (n = 3,5,7,11; rel err ≤ 7×10⁻⁵¹). G₁₂₁ computed numerically (no closed form published): 3.54880113413074478010727440546246876515786749179993826540053 (60 digits). PSLQ found no minimal polynomial for G₁₂₁ (or G₁₂₁², G₁₂₁⁴) of degree ≤ 8 at 80 digits — expected: its degree is large (class number of Q(√−121) is 5, but G₁₂₁ lives in a larger ring class field); CM-theory computation of its minimal polynomial is a Phase-2 item, not attempted here.

---

## E4. Modular equations via Sturm bound

### Proved (full Sturm argument)
A prover was written (`e4_modulareq.py`) that, for each term, checks: Newman congruences Σd·e_d ≡ 0 and Σ(N/d)e_d ≡ 0 (mod 24); character via D_eff = (−1)^k∏d^{e_d} (two terms share a character iff the ratio of D_eff's is a square in Q — this correctly gives the χ_{−4} character for the weight-1 quintic terms); holomorphicity at every cusp via Ligozat's formula ord_{1/c} = (N/(24·gcd(c²,N)))·Σ e_d·gcd(c,d)²/d (plus the shift at ∞); then compares q-expansions to the Sturm bound B = ⌈k·[SL₂:Γ₀(N)]/12⌉. Using φ(q) = η(2z)⁵/(η(z)²η(4z)²) and χ(q) = q^{1/24}η(2z)²/(η(z)η(4z)) (both verified against their series first):

| identity | level | weight | Sturm bound | result |
|---|---|---|---|---|
| p=3: φ⁴(q) − φ⁴(q³) − p³φ⁴(q³) = 0 | Γ₀(12) | 2 | 4 | ✅ proved; sum is 0 to q⁶⁰ |
| p=5: φ²(q) − φ²(q⁵) − pφ²(q⁵) = 0 | Γ₀(20) | 1 | 3 | ✅ proved; sum is 0 to q⁶⁰ |
| p=7: φ⁸(q) − (2+5p)φ⁴(q)φ⁴(q⁷) + (1−p)³φ⁸(q⁷) = 0 | Γ₀(28) | 4 | 16 | ✅ proved (7 terms); sum is 0 to q⁶⁰ |

All Newman/character/cusp checks pass for every term (characters: trivial for p=3,7; χ_{−4} for all three p=5 terms). These are genuine proofs of Berndt's multiplier relations, by the modern route — the calibration the programme needs.

### Degree 11: numerical only (honest gap)
Schröter's equation (Berndt Part III, Entry 7(i)),
(αβ)^{1/4} + ((1−α)(1−β))^{1/4} + 2(16αβ(1−α)(1−β))^{1/12} = 1,
α = λ(τ), β = λ(11τ), λ = 16η(τ/2)⁸η(2τ)¹⁶/η(τ)²⁴,
was verified numerically to 50 digits at τ = i, 0.5i, 2i, 0.3i (|LHS−1| ≤ 2×10⁻⁴⁷). A full Sturm-bound proof was **not** done: the equation lives on Γ(2)∩Γ₀(11) (λ is a modular function for Γ(2)), involves 4th/12th roots, and clearing them to a polynomial in (α,β) needs a resultant elimination the present eta machinery (integer-d η-quotients on Γ₀(N) only) does not support. Building that infrastructure — half-integer eta products, the Γ(2) valence/Sturm bound, and the resultant elimination — is a Phase-2 infrastructure item. The degree-17 Entry 12 equations were not attempted for the same reason (and the OCR of Berndt Part III is too garbled to extract them reliably).

---

## Phase-2 scouting (p = 11)

### Flagship data (computed, 50 digits)
- q₁ = e^{−π/√11}: p = 0.416439434275864095678569298496 (χ-quotient = 4√2/G₁₁¹⁰); m = φ²(q₁)/φ²(q₁¹¹) = **√11** exactly (PSLQ degree 2: [11, 0, −1]; this is expected from the theta transformation, a consistency check).
- u_k = 1.8352467839, 1.4190121186, 0.9294325201, 0.5336551405, 0.3223982758 (strictly decreasing ✓); w_k = u_k¹¹/p = 1910.23738646, 112.79641607, 1.07360532339, 0.00240057800, 0.00000939221; Πw_k = p⁶ ✓ (0.005215681788).
- Naive PSLQ minimal-polynomial test on the w_k (degree ≤ 5, 50 digits): **no relation found** — but this test is inconclusive as designed: the w_k span 8 orders of magnitude (ill-conditioned) and, more importantly, at the flagship point the auxiliary quintic's coefficients lie in Q(p, m) with p = 4√2/G₁₁¹⁰ and m = √11, so the w_k live in a compositum, not purely in Q(ζ₁₁)⁺. The correct test is the power-sum relation below.

### Power-sum relations: honest negative + validated method
- **p=11:** PSLQ fishing for a quadratic relation A(p,m)s₁² + B(p,m)s₁ + C(p,m) = 0 with A,B,C of total degree ≤ 2 (18 monomials) and ≤ 3 (30 monomials) found candidate relations at q=0.4 — but both **failed verification catastrophically** at other q (relative errors 6–322). They are spurious (PSLQ artefacts at 60 digits with 18–30 monomials). Conclusion: s₁ for p=11 satisfies NO low-degree quadratic relation over Q(p,m) with small coefficients. This is consistent with the genus-1 obstruction: X₀(11) has genus 1, so the function field is not Q(p,m) — a third generator is probably needed. Finding it is Phase-2 Step 2's real work.
- **p=7 calibration (method validated):** the same PSLQ fishing also produced a spurious relation — so the method was fixed by deriving the relation symbolically instead. With r(ξ) known, sympy gives the exact quadratic for s₁ = Σα²/β (the two roots are the two cyclic orders): T² − (S⁺+S⁻)T + S⁺S⁻ = 0 with S⁺+S⁻ and S⁺S⁻ explicit polynomials in (p,X) (degree up to X⁵ in the constant term). Numerically, the dissection's s₁ satisfies this quadratic at q = 0.2…0.6 with residuals ≤ 1×10⁻⁴⁵, and is consistently the smaller root. The machinery works — it just needs r₁₁ as input, which is exactly what Phase 2 must produce.

### What remains for Phase 2 (in priority order)
1. **The modular equation (Step 1):** find the explicit polynomial relation between m = φ²(q)/φ²(q¹¹) and p_prod. Strategy per the report: write both as eta-quotients on Γ₀(11), compute q-expansions, find the relation by linear algebra/Gröbner elimination over the function field. Note the genus-1 warning: the relation may need a third generator (a coordinate on X₀(11)); the failed PSLQ test above is evidence for this.
2. **The auxiliary quintic (Step 2):** compute power sums of u_k¹¹/p as modular functions; determine whether they lie in Q(p,m) or need the third generator; identify the Galois group (test cyclicity C₅). The p=7 calibration script (`p7_calibration.py`) is the template: once r₁₁'s coefficients are known as functions of (p,m[,t]), the same symbolic derivation gives the resolvent for s₁.
3. **Root formula and ordering (Step 3):** find the monomial pattern u_k = (⋯)^{1/11} by matching q-expansions; the ordering will need the p=7 method (monotonicity + a symmetric identity analogous to Lemma 3.8's "=41").
4. **Flagship + G₁₃₃₁ (Step 4):** at q₁, p = 4√2/G₁₁¹⁰ and m = √11 are now known exactly; the remaining input is r₁₁ at that point. G₁₂₁'s minimal polynomial (CM theory) is also still open.
5. **Infrastructure:** extend the Sturm prover to Γ(2) (half-integer eta products) to finish the degree-11 Schröter proof; extract Entry 12 (degree 17) cleanly from a better source than the garbled OCR.

### Errors caught and fixed during Phase 1
- E3's first draft used p = 2/G₁₁⁴ for the p=11 flagship (the p=5 shape). Correct is p = 4√2/G₁₁¹⁰ (from the general formula p_prod = 2^{(p−1)/4}G_n/G_{p²n}^p). Fixed; both routes now agree to 50 digits.
- The p=7 "Lemma 3.8" check first used Σs_i² (= 4.27); the lemma is Σs_i⁷ = 41 (the s_i are (ratio)^{2/7}, so s_i⁷ = ratio²). Fixed; Σs_i⁷ = 41.0 confirmed.
- The E4 "first nonzero coefficient" line initially reported truncation artefacts past the computed range; fixed to report only within the computed range (no nonzero coefficient to q⁶⁰ in any of the three proofs).
- The p=11 PSLQ power-sum relations (D=2, D=3) are spurious — reported as such, not as results.
