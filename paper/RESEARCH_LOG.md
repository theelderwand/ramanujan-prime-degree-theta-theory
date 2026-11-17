# Research Log — Ramanujan Prime-Degree Theta Function Paper

## Project
Paper: "A general product formula and computational verification of Ramanujan's
prime-degree theta function identities" (target: The Ramanujan Journal;
backup: Hardy–Ramanujan Journal).
Author: Ashwin Sathya Raghunathan. Computational work and drafting done with AI
assistance (Muse); all mathematical claims verified against cited sources and by
recomputation. AI assistance is disclosed in the paper's Acknowledgments.

## Timeline
- 2026-11-08: Churn exploration in ~/workspace/modular-forms/ (554 eta-quotient
  identities); notebook hunt launched.
- 2026-11-13: Notebook analysis — no unexamined groundbreaking identity; one open
  claim (Lost Notebook p. 336, second identity) + Candidate 2 (prime-degree theta
  theory). User chose Candidate 2 for deep research.
- 2026-11-14: Deep research on prime-degree theory (report in
  ~/workspace/research_notes/ramanujan-prime-degree-theta-theory-20261114-0146/).
- 2026-11-14: Phase 1 computational reproduction (~/workspace/theta-phase1/,
  PHASE1_REPORT.md). Key new result: p_prod chi-formula proved for ALL odd n.
- 2026-11-14: Generalizations research (report in
  ~/workspace/research_notes/ramanujan-theta-generalizations-20261114-1856/).
  Key finding: X0(11) genus 1 is the obstruction.
- 2026-11-14: Paper drafted (paper.tex, 9 pages, amsart). p=11 attack (Phase 2
  Step 1) launched in parallel as a subagent.

## Decisions
1. **Paper shape**: NOT the full p=11 solution (thesis-scale, may fail). Instead a
   paper whose central new theorem is the general product formula (Theorem 2.1),
   plus computational verification, Sturm-bound proofs, and the p=11 frontier as
   a documented open problem with a concrete attack plan.
2. **Author**: Ashwin Sathya Raghunathan (sole). AI assistance disclosed in
   Acknowledgments.
3. **Target journal**: The Ramanujan Journal (primary). Hardy–Ramanujan Journal
   as backup (Berndt–Rebák 2021 survey appeared there).
4. **LaTeX**: amsart class, pdflatex. Installed texlive-latex-base +
   texlive-latex-recommended via apt (2026-11-14).
5. **Honesty rules applied**: every theorem stated is verified; the degree-11
   Schröter equation is numerical-only and labelled as such; the PSLQ negative
   search is reported as evidence, not a theorem; H1/H2 is stated as open;
   G_121 has no closed form (numerical value only).

## What the paper proves (new)
- Theorem 2.1: p_prod = 2^{(n-1)/2} q^{(n^2-1)/24} chi(q)/chi^n(q^n) for ALL odd n
  (analytic proof via Jacobi triple product). Corollary: G-form
  p_prod = 2^{(n-1)/4} G_n/G_{n^2}^n. Hence c_11 = 4√2 (was a prediction, now a
  theorem).
- Sturm-bound proofs of the p=3,5,7 multiplier relations (new prover:
  Newman + character via D_eff + Ligozat cusp orders + coefficient check to q^60).
- G_121 numerically to 60 digits (new value; no closed form known).

## What the paper does NOT prove
- No p=11 theta identity (auxiliary quintic r_11 unknown; H1 vs H2 unresolved).
- No (m, p_prod) relation on Gamma0(11) (genus-1 obstruction; third generator
  likely needed).
- No closed form for G_121 or G_1331.
- Degree-11 Schröter equation: numerical only (needs Gamma(2) infrastructure).
- Degree-17 Entry 12 equations: not attempted (OCR garbled).

## Dead ends / errors caught
- E3 first draft used p = 2/G_11^4 for the p=11 flagship (p=5 shape). Correct is
  4√2/G_11^10. Fixed; both routes agree to 50 digits.
- Lemma 3.8 is Σs_i^7 = 41, not Σs_i^2. Fixed.
- p=11 PSLQ power-sum relations (D=2, D=3) are spurious — reported as such, not
  as results. Consistent with genus-1 obstruction.
- Published septic terms appear in order (v,u,w) relative to dissection's u,v,w.
  Noted in paper §3.2.
- "Roots lie in Q(zeta_p)^+" is too strong (false at p=5 special point). Paper
  states the corrected cyclic-over-coefficient-field version (§5.3).

## p=11 attack — SOLVED (2026-11-14, evening)
- **The relation F(m, p_prod) = 0 was FOUND and PROVED by Sturm bound.**
  G(Y,P) = 0 with Y = m², P = p/32: a plane quintic, irreducible over Q,
  total degree 5, largest coefficient 33554432 = 2^25. Pure parts are perfect
  fifth powers: (1−Y)^5 and (1+32P)^5 = (1+p)^5.
- **Proof**: m, p are eta-quotients on Gamma0(44) (level 4·11, NOT level 11);
  K = eta(2z)^2 eta(22z)^2, K^25 weight 50, Sturm bound 300; all 21 terms
  K^25·Y^i·P^j holomorphic weight-50 forms; coefficients q^0…q^300 all zero
  (checked to q^400). Exact series check to q^600; numerics ≤7×10⁻⁵⁰.
- **Two corrections to prior expectation**:
  (a) Genus is 2, NOT 1: the plane quintic has exactly four ordinary nodes
      (delta=1 each; Y-coords are roots of 59049Y^4+280422Y^3+293768Y^2−64614Y−14641),
      so genus = 6−4 = 2. The fourth node was initially missed (its P-coordinate
      is a double root); a genus-3 count violated Riemann–Hurwitz and forced
      its discovery.
  (b) The degree-2 map is the w4 quotient: C ≅ X0(44)/w4; the extension is
      Q(m,P)/Q(Y,P), m ↦ −m. w11/w44 sends Y ↦ 121/Y and preserves C; no simple
      closed form for its action on P was found.
- **Paper revised** (v2, 11 pages): title changed to "The degree-eleven modular
  equation in Ramanujan's prime-degree theta function theory"; new §5 with
  Theorem 5.1 (G(Y,P)=0), structural form, Sturm proof, genus-2 analysis;
  §6 (frontier) rewritten — genus-1 obstruction corrected, attack plan updated
  (Steps 0 and 2 done; H1/H2 now Step 1); open problems updated (F(m,p)=0
  removed; flagship/G_1331, w11-action, hyperelliptic model added).
- Scripts + report: ~/workspace/theta-phase1/p11/ (P11_ATTACK.md, step0–step5).
- One labelling slip in step3_sturm.py output: it prints "G(11,1) = 76026720"
  for cusp 0, but the cusp value is m=11 so Y=m²=121; G(121,1)=0 verified
  directly. Not a mathematical error.

## Before submission — checklist
- [x] Re-run all Phase 1 scripts end-to-end and confirm the numbers in the paper
      match the script output exactly (esp. Table in §3.3, flagship values).
      DONE 2026-11-14 (deep re-verification; see paper/REVERIFICATION_LOG.md).
      All Sturm proofs, flagships, class invariants re-confirmed.
- [x] Decide p=11 attack outcome (see above); revise §5 if it succeeds.
      DONE — Theorem 5.1 proved; p=13, p=17 added (Theorems 5.2, 5.3).
- [x] Proofread for typos, reference completeness (check DOI for Vasuki–
      Yathirajsharma IJNT paper — currently cited without volume/year in the
      research report; verify before citing or drop it).
      DONE 2026-11-14 — verified against the web: Ramanujan J. 61 (2023),
      1145–1150, DOI 10.1007/s11139-022-00630-z. Bibliography entry correct.
- [x] Check The Ramanujan Journal formatting requirements (they accept amsart;
      confirm reference style and whether they want MSC codes + keywords).
      The Ramanujan Journal accepts amsart; MSC codes + keywords are standard
      for Springer journals and are now included.
- [x] Add MSC 2020 codes (11F27, 11F11, 33D10, 11Y40) and keywords.
      DONE 2026-11-14 (paper v4).
- [ ] Confirm author affiliation/address block for submission.
- [ ] Consider asking Berndt or Rebák for a read (they are the experts; a
      pre-submission check would catch errors and is standard courtesy given the
      paper builds directly on their work). NOT DONE — needs the user's word
      before contacting anyone.

## 2026-11-14 — deep re-verification + paper v4 (16 pages)
- Every script re-run fresh; see paper/REVERIFICATION_LOG.md for the full log.
- **Errors found and fixed** (all in genus cross-checks / labels, none in the
  theorems themselves): p=13 fibre over m=0 is (0,−1/64), not (0,1/4);
  Riemann–Hurwitz R for p=13 is 8, not 10; R for p=17 is 12, not 10;
  eta-quotient displays no longer write q⁵/q⁷/q¹² as a separate factor;
  1/52 and 1/68 removed as duplicate cusp rows (Γ₀(52), Γ₀(68) have 6 cusps).
  Root cause of the R errors: p13/step4e_genus_check.py evaluated F_m at a
  point off the curve and concluded the m0=33.545 fibre was smooth
  ramification; it is in fact the 4th node (annotated in place, kept for audit).
- **New §7**: conjectural template for general prime p (Conjecture 7.1), with
  w₄ mechanism, exponents, pure parts, support triangle, genus problem, and
  p=19/23 predictions. Title/abstract/conclusion updated; MSC + keywords added.
- Build: pdflatex was missing on this machine; installed via apt
  (texlive-latex-base/recommended, fonts-recommended, lmodern). Final: 16 pages,
  0 overfull boxes, 0 undefined refs/cites.

## 2026-11-14 — p=13 and p=17 added (Theorems 5.2, 5.3)
- p=13: F(m,P) = 0 (plane septic, genus 2), proved on Gamma0(52), K^18 weight 36,
  Sturm bound 252. p=17: F(m,P) = 0 (degree 12, genus 3), proved on Gamma0(68),
  K^32 weight 64, Sturm bound 576. Both verified exactly to q^600 and numerically
  to 50 digits (~1e-51). Scripts: ../p13/, ../p17/; report: ../P13_P17_REPORT.md.
- Paper v3: 14 pages, 0 overfull boxes, 0 undefined refs/cites. Title/abstract/
  intro/conclusion/open problems updated. New section 6 (§6, "The degree-thirteen
  and degree-seventeen modular equations") inserted before the p=11 frontier.
- Note: Ramanujan's Entry 12 degree-17 equations are in (alpha,beta) variables;
  Theorem 5.3 is the first in the theta variables (m, p_prod). The Gamma(2) open
  problem was reworded accordingly.
- A uniform theorem for all primes p is NOT claimed: the (m,P) vs (Y,P) split,
  the support patterns, and the genus sequence 2,2,3 have no proved rule.
  Next cases by the same method: p=19 (level 76), p=23 (level 92).
