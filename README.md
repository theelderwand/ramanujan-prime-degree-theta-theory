# The Degree-11, 13, and 17 Modular Equations, and a Conjecture for General Prime Degree, in Ramanujan's Prime-Degree Theta Function Theory

By Ashwin Sathya Raghunathan.

In his second notebook, Ramanujan recorded that the theta functions
φ(q^{1/n}), ψ(q^{1/n}), f(−q^{1/n}) for prime n can be expressed as a sum of
½(n−1) nth roots plus φ(q^n), and that the auxiliary equations needed are
quadratics only for degrees 5, 17, 257, cubics for 7, 13, 19, 37, 73, 97,
109, 163, 193, and quintics for 11, 41, 101. The cubic, quintic, and septic
cases have been proved in recent work of Berndt–Rebák and Rebák; no case for
any prime p ≥ 11 was known.

This repository proves the first three missing cases and states a precise
conjecture for general prime degree.

## Results

- **Theorem 5.1 (degree 11).** An explicit modular equation G(Y, P) = 0 with
  Y = φ⁴(q)/φ⁴(q¹¹) and P = p_prod/32, a plane quintic of genus 2, proved by
  the Sturm-bound method on Γ₀(44).
- **Theorem 5.2 (degree 13).** An explicit modular equation F(m, P) = 0 with
  m = φ²(q)/φ²(q¹³) and P = p_prod/64, genus 2, proved by the Sturm-bound
  method on Γ₀(52).
- **Theorem 5.3 (degree 17).** An explicit modular equation F(m, P) = 0 with
  m = φ²(q)/φ²(q¹⁷) and P = p_prod/256, genus 3, proved by the Sturm-bound
  method on Γ₀(68).
- **General product formula (Theorem 2.x).** For every odd n,
  p_prod = 2^{(n−1)/2} q^{(n²−1)/24} χ(q)/χⁿ(qⁿ), proved analytically from
  the Jacobi triple product.
- **Conjecture 7.1.** A uniform template for the degree-p modular equation
  for general prime p: the Atkin–Lehner w₄ mechanism, the pure-part
  exponents e_X = (p²−1)/24 and e_P = (p−1)/4 or (p−1)/2, the full
  triangular support, and explicit predictions for p = 19 and p = 23
  (82 and 144 monomials; genera 4 and 5). Each part is labelled by its
  status — theorem, verified, pattern, or guess — with per-part falsifiers.

The paper also verifies the dissection, the three flagship evaluations of
degree 7, and eleven class invariants (including G₁₂₁ to 60 digits)
computationally to 50 digits, and gives Sturm-bound proofs of the multiplier
relations of degrees 3, 5, and 7.

## Contents

- `paper/paper.tex` — the paper (v4, 16 pages). Build with `paper/build.sh`.
- `paper/REVERIFICATION_LOG.md` — the deep re-verification log: every check
  re-run fresh, exact numbers, and all corrections made.
- `CONJECTURE.md` — the full conjecture document with the w₄ mechanism,
  pure-part exponents, support, genus problem, and p = 19/23 predictions.
- `PHASE1_REPORT.md`, `P13_P17_REPORT.md` — computational reports.
- `theta.py`, `e1_dissection.py`, `e2_flagships.py`, `e3_classinvariants.py`,
  `e4_modulareq.py` — Phase 1 scripts.
- `p11/`, `p13/`, `p17/` — per-prime scripts: search, exact verification,
  Sturm proof, structure and genus analysis.
- `w4_verify.py`, `cusp_verify.py` — fresh verification of the w₄ action
  and the cusp values used by the conjecture.

## Reproducing the results

All scripts run with plain Python (`sympy`, `mpmath`); no SageMath, PARI, or
Mathematica is needed. For example, the degree-11 Sturm proof:

```
python p11/step3_sturm.py
```

and the degree-13 and degree-17 proofs:

```
python p13/step3_sturm.py
python p17/step3_sturm.py
```

## Status

The three modular equations and the general product formula are proved.
The auxiliary quintics r₁₁, r₁₃, r₁₇ and the root formulas remain open, as
does the genus formula (the number N_CM(p) of fixed CM points of w₄ has no
known formula). The paper is intended for submission to *The Ramanujan
Journal*.
