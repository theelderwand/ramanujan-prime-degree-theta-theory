# Deep re-verification log — 2026-11-14 (paper v4)

Every check below was re-run fresh in this session (2026-11-14 evening).
Nothing was taken on trust from earlier reports. All scripts were run with
`~/workspace/.venv-theta/bin/python` unless noted.

## Part 1: re-verification results

### 1. p=11 Sturm proof — PASS
- `p11/step3_sturm.py`: Γ₀(44), index 72; weight 50; Sturm bound
  ⌈50·72/12⌉ = 300; all 21 terms weight 50, Newman ok, trivial character,
  holomorphic at all cusps; coefficients q⁰…q³⁰⁰ of K²⁵·G all zero; no
  nonzero coefficient to q⁴⁰⁰. ⇒ K²⁵·G ≡ 0, so G(Y,P) = 0.
- `p11/step2_exact.py`: nullspace dimension 1 (exact, R=90); exact check to
  q⁶⁰⁰ all zero; numerics at 8 points, rel err ≤ 6.5×10⁻⁵⁰.
- Cusp values: G(1,0) = 0, G(121,1) = 0, F(1,0) = 0, F(11,32) = 0.
  (The script prints "G(11,1) = 76026720" — a labelling slip; the cusp-0
  point is (Y,P) = (121,1) since m(0) = 11, Y = m². Not a math error.)

### 2. p=13 Sturm proof — PASS
- `p13/step3_sturm.py`: Γ₀(52), index 84; K¹⁸ weight 36; Sturm bound 252;
  all 17 terms pass Newman/character/cusp checks; coefficients q⁰…q²⁵² of
  K¹⁸·F all zero; no nonzero coefficient to q⁴⁰⁰. ⇒ F(m,P) = 0.
- `p13/step2_exact.py`: exact check to q⁶⁰⁰ all zero; numerics 50 digits.
- Cusp values: F(1,0) = 0, F(13,1) = 0.

### 3. p=17 Sturm proof — PASS
- `p17/step3_sturm.py`: Γ₀(68), index 108; K³² weight 64; Sturm bound 576;
  all 35 terms pass; coefficients q⁰…q⁵⁷⁶ of K³²·F all zero; no nonzero
  coefficient to q⁸⁰⁰. ⇒ F(m,P) = 0.
- `p17/step2_exact.py`: exact check to q⁶⁰⁰ all zero; numerics 50 digits
  (flagship m = √17, P ≈ 0.00017769354832319060602).
- Cusp values: F(1,0) = 0, F(17,1) = 0.

### 4. General p_prod formula (Phase 1, E1) — PASS
- `e1_dissection.py`: dissection LHS−RHS == 0 and p_prod χ-formula == 0
  exactly (series to t¹²⁰) for n = 3, 5, 7, 9, 11, 13, 17; numerics at
  50 digits (rel err ≤ 1.7×10⁻⁴⁹). The analytic proof via the Jacobi triple
  product is in `paper/paper.tex` §2 (Theorem thm:pprod).

### 5. p=3,5,7 multiplier relations (E4) — PASS
- `e4_modulareq.py`: p=3: Γ₀(12), weight 2, Sturm bound 4, coeffs 0…4 zero.
  p=5: Γ₀(20), weight 1, bound 3, coeffs 0…3 zero. p=7: Γ₀(28), weight 4,
  bound 16, coeffs 0…16 zero. Degree-11 Schröter relation is numerical
  only (50 digits, |LHS−1| ≤ 2×10⁻⁴⁷), as stated in the paper.

### 6. Flagships and class invariants (E2, E3) — PASS
- `e2_flagships.py`: three flagships reproduced to 50 digits
  (rel err ≤ 4.8×10⁻⁵¹); Lemma 3.8 Σs_i⁷ = 41 confirmed (rel err 1.5×10⁻⁵⁰);
  cyclic root ordering (α,β,γ) = k = 3,2,1 confirmed.
- `e3_classinvariants.py`: 11 class invariants verified to ~50 digits
  (incl. new G_99, G_729, G_125, G_625); reciprocity G_{1/n} = G_n for
  n = 3,5,7,11; G₁₂₁ = 3.54880113413074478010727440546246876515786749179993826540053
  (60 digits).

### 7. Conjecture claims — PASS (with one method note)
- **w₄ mechanism** (`w4_verify.py`, fresh script, 60 digits, eta via
  pentagonal series): w₄(m) = (−1/p)·m and w₄(P) = P for p = 5,7,11,13,17
  (p=3: m→−m, P→ζ₃P), 3 test points each, 18/18 checks pass, worst error
  2.7×10⁻⁴⁶. NOTE: a first attempt at 40 digits / 60 pentagonal terms FAILED
  at larger p (errors up to 10⁻⁶) — a series-truncation artifact, not a
  mathematical failure: Im(w₄τ) can be as small as ~4×10⁻⁴, where 60 terms
  leave an O(1) tail. Resolved with 400 terms / 60 digits. The eta-multiplier
  proof is still unwritten.
- **Pure-part exponents**: e_X = (p²−1)/24, e_P = (p−1)/4 (p≡1 mod 4) or
  (p−1)/2 (p≡3 mod 4) — verified against all six cases (5: 1,1; 7: 2,3;
  11: 5,5; 13: 7,3; 17: 12,4). p=19: (15,9); p=23: (22,11).
- **Support triangles**: full triangle {(i,j): i·e_P + j·e_X ≤ e_X·e_P},
  all monomials present — p=11 (21), p=13 (17), p=17 (35), verified by
  exact coefficient extraction.
- **Pure parts**: G11(Y,0) = −(Y−1)⁵, G11(0,P) = (1+32P)⁵;
  F13(m,0) = −(m−1)⁷, F13(0,P) = (1+64P)³;
  F17(m,0) = (m−1)¹², F17(0,P) = (256P−1)⁴. All exact.
- **P(1/2) = (2/p)/2^{(p−1)/2}** (`cusp_verify.py`, 60 digits): at
  τ = 1/2 + 0.005i, rel err 1.7×10⁻²⁴ (p=11), 1.3×10⁻²⁰ (p=13),
  1.5×10⁻¹⁵ (p=17), 8.3×10⁻¹⁴ (p=19), 3.1×10⁻¹¹ (p=23). PASS.
- **m(1/p) = (−1/p)**: direct numerical limit is INFEASIBLE with the
  pentagonal series (the local parameter at cusp 1/p is e^{−2π/(4p²t)},
  still O(1) at t = 0.1 and needing 400+ digits at t ≪ 0.001). Verified
  DEDUCTIVELY instead: the pure parts (explicit above) plus P(1/p) = 0
  (Ligozat, ord_{1/p}(P) = (p²−1)/24 > 0) force X(1/p) = 1; for p≡1 mod 4
  this gives m(1/p) = 1 = (−1/p); for p≡3 mod 4 the sign follows from
  w₄ swapping ∞ ↔ 1/p (checked by representative: W₄ = [[4,1],[4p,12]]
  sends ∞ ↦ 1/p) together with w₄(m) = (−1/p)·m and m(∞) = 1.
  The paper states this deduction, not a direct numerical limit.

### 8. Genus computations — PASS (with two buggy scripts found, see §10)
- **p=11, genus 2**: `p11/step5_final.py` finds the 4th node
  (Y,P) = (−0.1444296737186051, 0.0040600200662527), ordinary, δ=1;
  four ordinary nodes, genus = 6−4 = 2. Independently confirmed by the
  discriminant: disc_P(G) = 1267650600228229401496703205376·Y⁴·(3125Y⁴+…+45753125)
  ·(59049Y⁴+280422Y³+293768Y²−64614Y−14641)² — the squared quartic's roots
  (−2.9399, −1.9622, −0.1444, 0.2976) are the node Y-coordinates.
  Riemann–Hurwitz on Y (deg 5): R = 4 (Y=0, e=5) + 4 (simple ramif) +
  4 (Y=∞, e=5) = 12 ⇒ g = 2. ✓
- **p=13, genus 2**: four affine nodes verified independently at 50 digits
  (m,P) = (1, 1/4), (33.5454297306, −93.1708371306),
  (−2.38279097366, −0.0868836622519), (−0.162638756944, 0.00772079286153);
  each has |F|,|F_m|,|F_P| ≤ 10⁻³⁰ and nonzero tangent-cone discriminant.
  Discriminant: disc_P(F) = −68719476736·m²·(m−1)²·(m³−31m²−85m−13)²·
  (27m⁴−212m³+1194m²−2756m+4563); the squared factors mark the nodes, the
  separable quartic marks 4 simple ramification points.
  Riemann–Hurwitz on m (deg 3): R = 2 (m=0, e=3) + 4 (simple ramif) +
  2 (m=∞, e=3) = 8 ⇒ g = 2. ✓
- **p=17, genus 3**: `p17/step4b_blowup.py`: δ at infinity = 40
  (multiplicity sequence 8,4,4); `p17/step4c_affine.py`: 12 affine nodes,
  each verified by nonzero tangent-cone discriminant; genus = 55−40−12 = 3.
  Independently confirmed by the discriminant:
  disc_P(F) = −1267650600228229401496703205376·m³·(sextic)·S(m)², where the
  sextic is separable (6 simple ramification points) and S(m) is the separable
  degree-12 node polynomial.
  Riemann–Hurwitz on m (deg 4): R = 3 (m=0, e=4) + 6 (simple ramif) +
  3 (m=∞, e=4) = 12 ⇒ g = 3. ✓

### 9. Irreducibility — PASS
- G₁₁, F₁₃, F₁₇ all irreducible over ℚ (sympy `Poly.factor_list()` gives a
  single factor of multiplicity 1 in each case). NOTE: a naive
  `factor(F) == F` test returns False for G₁₁ because `factor()` reorders
  terms; the correct test is `factor(F) − F == 0` or `factor_list()`.

## Part 2: corrections made (all verified above)

| # | Location | Was | Now (correct) | Evidence |
|---|----------|-----|---------------|----------|
| 1 | paper §6, p=13; P13_P17_REPORT.md §1 | fibre over m=0 is (0, 1/4) | (0, −1/64), F_m(0,−1/64) = −169 = −13² | F(0,1/4) = 4913 ≠ 0; F(0,−1/64) = 0. The value 1/4 is the P-coordinate of the node at m=1. |
| 2 | paper §6, p=13; report §1 | Riemann–Hurwitz on m, R = 10 | R = 8 (2 at m=0, 4 simple ramif, 2 at m=∞) | R=10 would give g=3, contradicting genus 2. Source: `p13/step4e_genus_check.py` bug (see §10). |
| 3 | paper §6, p=17; report §2 | Riemann–Hurwitz on m, R = 10 | R = 12 (3 at m=0, 6 simple ramif, 3 at m=∞) | R=10 would give g=2, contradicting genus 3. Discriminant sextic is separable (6 points). |
| 4 | paper §5–6 (three places) | `\pprod = 32q^5·η…` (and 64q⁷, 256q¹²) | `32·η…` with a note that the eta-quotient itself has q-order 5 (resp. 7, 12) at ∞ | Σ e·d/24 = 5, 7, 12 for the three edicts. The χ-form `\pprod = 32q⁵χ/χ¹¹` is unaffected (χ/χ¹¹ has constant term 1). |
| 5 | paper §6 (p=13, p=17); report §1–2 | cusp lists include 1/52 (resp. 1/68) as a separate cusp | removed; 1/52 ~ ∞ in Γ₀(52), 1/68 ~ ∞ in Γ₀(68) | Γ₀(52) and Γ₀(68) each have 6 cusps, not 7. |

## Part 3: paper v4 changes

- **New §7, "A conjectural template for general prime p"** (Conjecture 7.1):
  w₄ mechanism [Verified], variable rule [Theorem given w₄], exponents
  [Theorem given w₄ + degree-2], pure parts [Pattern with proved X-mechanism],
  support triangle [Pattern 6/6], degree-2 map [Proved p≤17], genus [Guess —
  no formula for N_CM(p)], predictions for p=19 (82 monomials, genus 4) and
  p=23 (144 monomials, genus 5), per-part falsification conditions.
- **Title** updated to include the conjecture; **abstract** updated;
  **MSC 2020** (11F27, 11F11, 33D10, 11Y40) and **keywords** added.
- **Conclusion**: uniform-template open problem now cites §7; added open
  problems for the w₄ eta-multiplier proof and for N_CM(p) from CM theory.
- **Vasuki–Yathirajsharma reference verified** against the web: Ramanujan J.
  61 (2023), 1145–1150, DOI 10.1007/s11139-022-00630-z — the bibliography
  entry was already correct.
- **Build**: pdflatex was missing from this machine; installed via apt
  (texlive-latex-base/recommended, fonts-recommended, lmodern). Final build:
  **16 pages, 0 overfull boxes, 0 undefined refs/cites**. (Two overfull boxes
  appeared on the first build of v4 — a 13.6pt overfull hbox in the new §7
  and a 1.8pt overfull vbox on p.1 from the longer abstract — both fixed by
  rewording/shortening. The 26 hyperref "Token not allowed in a PDF string"
  warnings are the standard math-in-section-title bookmark warnings and are
  harmless.)

## §10. Scripts found buggy or superseded (annotated in place, kept for audit)

- `p13/step4e_genus_check.py` — **WRONG CONCLUSION (g=3, R=10, R′=0)**.
  It evaluates F_m at (m0, P0=−56.45), the double root of F_P, which is NOT
  on the curve, and concludes the fibre over m0=33.545 is a smooth
  ramification point. The actual point on the curve is (33.545, −93.17),
  where F=F_m=F_P=0 with nonzero tangent-cone discriminant: a node
  (verified by `step4f_node_check.py` and independently here at 50 digits).
  Correct: 4 nodes, genus 2, R=8, R′=4. The report's old "R=10" came from
  this script.
- `p13/step4c_genus.py` — **SUPERSEDED**: finds 0 affine nodes, genus 6.
- `p11/step4b_singular.py` — **SUPERSEDED**: finds 0 singular points, genus 6
  (its w44 check at the end is still valid).
- `p11/step4c_genus.py` — **SUPERSEDED**: finds only 3 nodes, genus 3
  (misses the 4th node at Y=−0.1444; `step5_final.py` is authoritative).
- `p11/step3_sturm.py` output — labelling slip only: prints "G(11,1) =
  76026720" for cusp 0; the cusp value is Y=121, G(121,1)=0. Not a math error.

## Open items (unchanged, now in the paper)

- Auxiliary polynomials r₁₁, r₁₃, r₁₇ and root formulas; G₁₂₁'s minimal
  polynomial; w₁₁-action on P; hyperelliptic models.
- The w₄ eta-multiplier proof (most important missing proof).
- N_CM(p) from CM fixed-point theory (would turn the genus into a prediction).
- p=19 and p=23 computations (template predicts everything checkable).
- Pre-submission: consider asking Berndt/Rebák for a read (not done; needs
  the user's word). Author affiliation/address block still to be confirmed
  for submission.
