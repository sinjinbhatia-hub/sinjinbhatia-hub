# Review: ST-1 Inference foundations

Reviewer: adversarial pass (independent). Module folder: curriculum/03-statistical-modeling/01-inference-foundations/.
Status line is at the bottom. Items are numbered in order found; each has file, location, severity (HIGH/MED/LOW), what is wrong, fix, and resolution.

## Phase 1 findings

### Verification performed
- Ran `MPLBACKEND=Agg timeout 300 python3 project.py`: exit 0 in 14 s, 6 figures written. Imports: numpy, scipy, matplotlib, sympy, pathlib only.
- Recomputed every number in notes §3.3–§3.12, Examples 4.1–4.7 and all solutions (A1–A6, B1–B8, C1–C3) in /tmp/claude-0/st1_check.py and st1_check2.py (numpy/scipy/sympy). All agree: Ex 4.1 (0.2222/0.19/0.1818; n=5 σ²=4: 8.0 vs 5.76; n=100 spread 2.0%), Ex 4.2 (−238.1, 4.7619, 0.0648, (0.173, 0.427), 0.0789/0.0614/0.0825), Ex 4.3 (0.07, 0.08, 114.29, 816.33, 0.00875, 0.001225), Ex 4.4 (0.00833 vs 0.0333, ratio (n+2)/3, factor 12), Ex 4.5 (4.4190, 2.5994, 0.3928, 0.1378, 3.7865, 2.2273, 6.90, 12, 0.575), Ex 4.6 (3.2, 0.4, 6.25, (2.416, 3.984), 0.04076, 0.01630), Ex 4.7/B5 (0.03910, 2.883e-4, 2.658e-4, 1.0844), §3.11/B6 (0.3077, 0.003942 vs 0.0042, roots 0.1447/0.8553), §3.12 (D = 0.4189, KL quadrature agrees), A2 (1/3, 0.1217), A3 (1.75, 2.2913, 2.1), B2 (0.03333, 0.01515, 0.00833, 0.00826, c* = 12/11), B3 (0.1111, 0.1543, 0.1250, 0.1, 0.8), B4 (−0.8473, 0.3086), C2 (B(3/2,3/2) = π/8, I = 1/2, π²/4), C3 (ψ'(2) = 0.6449, 1.5505, 6.900).
- SymPy: (3.12) Fisher matrix and inverse; Exponential/Poisson/Geometric/Rayleigh information; Gamma Fisher matrix and [I⁻¹]_αα = α/(αψ'(α)−1); MoM asymptotic variance of the Gamma shape = 2α(α+1) by the delta method on (m̂₁, m̂₂); Uniform max moments and Var((n+1)/n·max) = θ²/(n(n+2)); B1 c* = 1/(n+1), MSE 2/(n+1); B2 c* = (n+2)/(n+1), MSE 1/(n+1)²; H(Exp(λ)) = 1 − log λ; KL gap; shrinkage identity (3.25); B5 variance via pgf; biased bound (3.15) below actual Var(σ̂²) by 2(n−1)σ⁴/n³; Normal exponential-family A(η).
- Step-by-step read of every derivation in §3 (Theorems 3.3–3.28) and every solution: signs, factors, index ranges, hypotheses all check. Exercise/solution numbering A1–A6, B1–B8, C1–C3 matches; no "similarly"/"left to the reader"/"it can be shown".
- Cross-references checked against the actual files: F5 Theorems 3.4, 3.8, 3.9, 3.12, 3.17, 3.19, 3.22, 3.26, 3.28, 3.29, 3.30, 3.32, 3.34, 3.45, 3.46, Definition 3.33, tags (3.14), (3.21), (3.24), (3.29), (3.30), §3.4, §3.11, §3.13; LA-14 Theorems 3.2, 3.6, 3.13, tags (3.6), (3.17), (3.38)–(3.39), §3.1 (Riesz), §3.9 (uses p for dimension, so "p = 1" is right); LA-8 Theorems 3.4, 3.19; LA-5 Gram–Schmidt; LA-12 covariant; F4 Poisson process; ROADMAP/briefs for ST-2 (Wilks, bootstrap), ST-3 (sandwich), ST-4, ST-5, ST-6 (variational), ML-1, ML-2 (natural gradient named in brief), ML-6 (EM), ML-10 (ELBO), SM-4 (ln Z). All exist as named, except item 1 below.
- Template: §0–§10 headings present; conventions box; 14 Try-it-first / 13 What-just-happened; 7 worked examples; A6/B8/C3; §8 specific; §10 matches brief's five sources; §9 has 11 questions.
- History: Gauss 1809, Laplace 1810, Pearson 1894 (Weldon's 1000 Naples crabs, nonic), Fisher 1912/1922/1925, Rao 1945, Cramér 1946, Fréchet 1943, Darmois 1945, Blackwell 1947, Wald 1949, Neyman–Scott 1948, Lehmann–Scheffé 1950, Efron–Hinkley 1978, Huber 1967, White 1982, Minka's Gamma initializer: all consistent with my knowledge. One soft spot: item 2.

### Findings

1. **notes.md §3.4, Gamma MoM bullet — MED (cross-reference).** "Mean α/β, variance α/β² (F5 (3.9))": F5's (3.9) is the Gamma *MGF* (λ/(λ−t))^a, not the mean and variance. The object cited does yield them, but only after differentiating. Fix: cite it as "differentiate the MGF, F5 (3.9)".

2. **notes.md §1, first paragraph — LOW (history).** "least squares, the method behind his recovery of the asteroid Ceres": whether Gauss used least squares in the 1801 Ceres orbit is contested (Gauss claimed to have used the method since 1795; the surviving Ceres computation does not show it, Stigler 1981). Cannot confirm from my own knowledge. Fix: soften to "which he said he had been using in orbit computations since 1795" and keep the Theoria motus date, which is certain.

3. **notes.md §3.7, Theorem 3.17 — LOW (math, incomplete proof).** The statement is an "if and only if"; the derivation proves only "attained ⇒ (3.16)" (steps 1–2) and the exponential-family consequence (step 3). The converse is one line (Var(s_n) = (nI₁)²Var(θ̂) = nI₁) and is absent. Also "with probability 1" is in step 1 but not in the statement. Fix: add the converse line and the qualifier.

4. **notes.md §3.10, Theorem 3.25 step 4 — LOW (math, hypothesis not stated).** "A PSD Hessian everywhere makes A convex" needs the domain to be convex (a PSD Hessian on a non-convex set does not give the supporting-hyperplane inequality). The natural parameter space {η : A(η) < ∞} is convex (Hölder), but the notes never say so. Fix: add the parenthetical.

5. **notes.md §3.6, regularity paragraph — LOW (consistency).** "Every family in this module satisfies (R) except Uniform(0, θ)": the exercises' Uniform(−θ, θ) (A3) and Uniform(θ, θ+1) (A6c) also fail (R1). Fix: reword to "every family in these notes", and name the two exercise families.

6. **notes.md §3.12, misspecification example — LOW (self-containedness).** Cites F5 §3.13 for H(Exp(λ)) = 1 − log λ; F5 gives it as h = 1 + log μ for Exponential(1/μ), one substitution away. Fix: add the one-line integral so the value is derived here.

7. **notes.md §3.9, Theorem 3.22 — LOW (template).** Has a Try-it-first but no What-just-happened (the only key derivation without one). Fix: add one sentence.

8. **project.py Part 6 — LOW (code, weakened check).** Example 4.2's exact Binomial 0.0789, normal 0.0614 and continuity-corrected 0.0825 are printed, not asserted; these are exact (non-Monte-Carlo) numbers quoted in the notes. Fix: assert each to 4 decimals.

9. **project.py Part 7 — LOW (code, weakened check).** Example 4.6 quotes "covers the truth 94.5% of the time"; coverage is printed but not asserted. Fix: assert coverage lies in (0.93, 0.965) at n = 20 and n = 200.

10. **project.py Part 4 — LOW (code, unchecked claim).** Notes §3.4 and Example 4.5 state the MoM shape estimator's asymptotic variance 2α(α+1) "(delta method on (m̂₁, m̂₂))" and the Gamma Fisher matrix inverse; neither is derived in the notes nor checked in code. Fix: add a SymPy cell that derives both (delta method with the Gamma moments; symbolic inverse of the 2×2 Fisher matrix) and asserts them.

11. **notes.md length — LOW (style, not fixed).** 9,452 words against the STYLE guideline of 3,500–7,000 for a module of this weight. The extra length is correct material tied to the brief's long topic list; trimming would mean rewriting correct sections, which this review does not do.

No HIGH-severity findings: every derivation step, every quoted number, and every assert in project.py checked out.
