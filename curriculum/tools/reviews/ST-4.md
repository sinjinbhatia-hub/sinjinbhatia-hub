# Review: ST-4 Generalized linear models

Reviewer: adversarial pass (independent). Module folder: curriculum/03-statistical-modeling/04-generalized-linear-models/.
Status line is at the bottom. Items are numbered in order found; each has file, location, severity (HIGH/MED/LOW), what is wrong, fix, and resolution.

## Phase 1 findings

### Verification performed
- Ran `MPLBACKEND=Agg timeout 300 python3 project.py`: exit 0 in 28 s, 6 figures written. Imports: pathlib, matplotlib, numpy, sympy, scipy only.
- Recomputed every number in notes §1–§4 and every solution independently (own IRLS loop, own BFGS maximization, closed forms, scipy.stats; not the project's `glm_fit`) in the scratch script /tmp/claude-0/st4_check.py. All agree: §1 (0.493 + 0.197x, 16.0 %, 0.031/0.235, skew −2.26, p = 1.9e−14; 14.5 %, 0.58/6.79), §2 (0 → 3 → 2.1991 → 1.6427 → 1.4165 → 1.3867 → 1.3863, errors 4.5e−4 then 1e−7), §3.4 (logit trace 1.13e−1, 4.66e−3, 8.45e−6, 2.78e−11, 4.7e−16; separation trace at iterations 1, 2, 3, 5, 10, 25 and the slope gain 2.00 per pass from iteration 9; se 1.64e5, z 0.0003; ridge (−3.9221, 1.1206) in 5 passes, first pass (−2.9302, 0.8372); C1 Fisher/Newton error sequences), §3.5 (234.16/182.44/169.27/168.25; 51.7/13.2/1.02; p 6.4e−13/2.8e−4/0.31; rejection 0.046, KS 0.83), §3.6 (D̄ = 5.99 and 4.6 % for grouped data; mean 70.5 sd 11.3 for Bernoulli; 61.09), §3.7 (1.596, 1.814, 1.701; 0.998/1.701/1.705; AUC 0.8229; HL 9.60 (p 0.29) vs 2639; softmax 0.36764/0.862), §3.8 (1.127 ± 0.079; 651.7 vs 399.9; coverage 0.743/0.927/0.957; φ̂ 2.52, r̂ 2.04; se 0.0460/0.0633/0.0666, AIC 871.2/823.1; ZIP 0.307, 39.9 %/22.3 %; Gamma (0.978, 0.521) vs (0.793, 0.517), −0.176), Ex 4.1 (−0.8473, 0.3567, 0.21, 0.021, 1.3083), Ex 4.2 (XᵀX, det 105, score (0, 3.5), z, (−0.4, 0.8) both routes, (−0.5549, 1.1098), MLE (−0.6070, 1.2140), se (1.2181, 0.9126), D 4.9560, 12 log 2 = 8.3178, LRT 3.3618 p 0.067, Wald 1.330 p 0.18), Ex 4.3 (3.34 [2.20, 5.08]; RR table 1.9935/1.8284/1.3364), Ex 4.4, Ex 4.5 (38/15, μ̂, the four terms, D 0.2400 p 0.97; 19.008, X² 20.53, p 2.7e−4), Ex 4.6 (32/243 = 0.13169; 12), Ex 4.7 (16.6578, 10.3157, 3.2118, 3.7311, all four r_D), Ex 4.8 (7/9), A1–A5, B1–B7 (B2 vertices and trapezoids; B3 both passes; B5 2.0/3.125; B6 61.086 and 0.115; B7 0.9332/0.9276), B8 deviance, C1(a), C2(a)(b)(d) (−0.1758, 0.8388, ν̂ = 2.873 from D/(2n) = 0.1840), C3(b) (1.4, 2.24).
- SymPy checks: (3.11) (observed weight for general V, g) ✓; every row of the table (b′ = μ, b″ = V) ✓; Gamma and Inverse Gaussian log-densities in the (3.1) form ✓; NB b′ and b″ ✓; (3.25) ✓; closed-form unit deviances Normal/Poisson/Gamma/IG from (3.19) ✓; softmax b derivatives ✓; C2(a) and C2(d) derivatives ✓.
- Step-by-step read of every derivation in §3 (Theorems 3.1–3.8, the Newton→IRLS algebra, the ridge step, the Hauck–Donner mechanism, the NB exponential-family rewrite, the max-ent Lagrangian) and every solution: signs, factors of 2 and π, transposes/inverses, index ranges and hypotheses check. Exercise/solution numbering A1–A5, B1–B8, C1–C3 matches; no "similarly"/"left to the reader"/"it can be shown"; no banned filler words (the grep hits are the "What just happened" headers). Equation tags 3.1–3.31 unique.
- Cross-references checked against the actual files: ST-1 §3.10, (3.21), (3.22), (3.23), Theorems 3.26, 3.27 ✓; ST-2 §3.6, Theorem 3.9, (3.16) ✓; ST-3 hat matrix, Gauss–Markov, F test, leave-one-out/Cook, WLS, AIC ✓; LA-14 (3.21) softmax gradient, (3.23) Newton ✓; LA-8 Theorem 3.16 ✓; LA-5 ✓; F5 Theorems 3.26, 3.28, 3.32, §3.4 (Poisson + Poisson), §3.5 (Gamma–Poisson) ✓; F4 (3.29) ✓; ST-5 folder exists; ST-6, ST-9, ML-3, ML-4, ML-7, SM-4, SM-7 are in ROADMAP with the named topics. Two soft spots: items 5 and 6.
- Template: §0–§10 with the template headings; conventions box; 6 Try-it-first / 7 What-just-happened pairs around Theorems 3.1, 3.2, 3.3, 3.5, 3.7, 3.8; 8 worked examples (brief asks ≥7), two from other fields (Geiger counter, crack detection); A1–A5, B1–B8, C1–C3 (16); B8 and C1 extend/use project.py; §8 specific; §9 has 12 questions; §10 matches the brief's five sources. Brief topics and project items all present. No REVIEW FOCUS section in the brief.
- History: Bliss 1935 with Fisher's appendix (probit by iterative weighting), Berkson 1944 (logit), Nelder and Wedderburn 1972, Jørgensen 1987 (exponential dispersion models/Tweedie), Albert and Anderson 1984, Hauck and Donner 1977, Firth 1993, Fahrmeir and Kaufmann 1985, Pregibon 1981, Hosmer and Lemeshow 1980, Wedderburn 1974: all consistent with my knowledge. No findings.

### Findings

1. **notes.md §3.4 (ridge paragraph) and solutions.md B3(c) — MED (wrong sign/meaning of a stated number).** Notes: "boundary at −3.922/1.121 = −3.5, midway between x = 3 and x = 4"; solutions: "−2.9302/0.8372 = −3.5: the boundary sits midway between x = 3 and x = 4". The ratio β₀/β₁ is −3.5, but the decision boundary is at x = −β₀/β₁ = +3.5; "−3.5" is not between 3 and 4. Fix: write x = −β₀/β₁ = 3.922/1.121 = 3.5 (and 2.9302/0.8372 = 3.5).

2. **project.py (Parts 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, exercise checks, Your turn, B8) — MED (asserts weaker than what the notes claim).** Dozens of numbers quoted in notes.md and solutions.md are only printed, or asserted qualitatively: Part 1's 0.493/0.197, 16.0 %, 0.031/0.235, −2.26, 14.5 %, 0.58/6.79 (asserted as > 0.1, ratio > 5, < −1.5, > 0.05, ratio > 10); the logit IRLS trace (only the quadratic-rate inequality); Ex 4.2's β⁽²⁾, MLE, se, deviance 4.9560, p-values, 6 passes; Part 5's 1.2057 and [0.7871, 1.6243], c = 1.701 (asserted 1.69–1.71), 0.998/1.701/1.705 (asserted 1.55–1.85); AUC 0.8229 and HL 9.60/0.29/2639 (asserted p > 0.05, p < 1e−4); the separation trace values, se 1.6e5, z 0.0003, ridge (−3.922, 1.121) in 5 passes; softmax 0.36764/0.862 (asserted within 2e−4 of L-BFGS, > 0.85); offsets 1.127 ± 0.079, 651.7/399.9 (asserted CI covers 1, gap > 100); Ex 4.5's D 0.2400, p 0.97, 19.008, 20.53, 2.7e−4 and Ex 4.7's 10.3157/3.2118/3.7311; Ex 4.6's 0.13169, 4.006, 12.014; coverage 0.743/0.927/0.957, φ̂ 2.52, r̂ 2.04, se 0.0460/0.0633/0.0666, AIC 871.2/823.1 (asserted < 0.8, > 0.9, > 0.9, 1.2 < r̂ < 3.5); the analysis-of-deviance table, rejection 0.046, KS 0.83, 70.5/11.3, 5.99/0.046, 61.086/0.115; A2, A3, A4, A5 errors, B4, B5, B7 values; C2's (0.978, 0.521), 0.334, (0.793, 0.517), −0.1758, 0.8388, 2.873, 2.991, 184.009, 0.8357; C3's (0.831, 0.513), 0.307, 0.399/0.223, −2507.3/−2845.8; C1's (1.9143, 2.6611), the two error sequences, 6 passes, se (0.3857, 0.3730)/(0.3791, 0.3670); B8's (0.4312, 0.3919), 0.222. Fix: assert each quoted value to its quoted precision (the script is seeded, so these are deterministic).

3. **notes.md §3.4 "The ridge fix" — LOW (notation).** "The gradient is u − λβ" with the intercept unpenalized; it is u − λI₀β (the next sentence uses I₀ correctly). Fix: write λI₀β.

4. **notes.md §3.7 first paragraph — LOW (false for negative coefficients).** "∂p/∂x_j = β_j p(1−p) ≤ β_j/4" fails when β_j < 0. Fix: |∂p/∂x_j| ≤ |β_j|/4.

5. **notes.md §3.4 "Convergence" — LOW (cross-reference to content that is not there).** "local convergence is quadratic (LA-14 §3.5)": LA-14 §3.5 states Newton's step (3.23) and says "ML-2 proves the convergence rates". Fix: cite LA-14 (3.23) for the method and say the quadratic rate is proved in ML-2.

6. **notes.md §3.7 "Loss, physics, boundary" — LOW (imprecise cross-reference).** "minimizes the KL divergence ... ST-1 §3.12": the object is ST-1 Theorem 3.28 (the MLE targets the KL projection), which sits in §3.12. Fix: name the theorem.

7. **exercises.md conventions line — LOW (wrong equation number).** "w_i and z_i from (3.10)": w_i is defined in (3.6); (3.10) defines z_i. Fix: "(3.6) and (3.10)".

8. **project.py markdown headers and one comment — LOW (wrong notes sections).** Part 10 says "(notes 3.6, Example 4.6)" (overdispersion is §3.8); Part 11 says "(notes 3.8)" (analysis of deviance is §3.5, the Bernoulli failure §3.6); Part 12 says "(notes 3.8)" (residuals/hat matrix are §3.6). Part 2's comment lists the Inverse Gaussian canonical link as "mu^{-2}"; it is −1/(2μ²). Fix the four strings.

9. **solutions.md A5(b) — LOW (incomplete description of a number list).** "errors 0.81, 0.26, 0.030, 4.5×10⁻⁴, 10⁻⁷" are the errors after passes 2–6; the list silently skips pass 1, whose error (1.61) is larger than the starting error (1.39) because the first step overshoots. Fix: say which passes, and name the overshoot (it is the "undamped Newton can overshoot" remark of §3.4).

10. **project.py Part 6 `roc_curve` — LOW (portability).** `np.trapezoid` exists only in NumPy ≥ 2.0; older installs have `np.trapz`. Fix: one-line fallback.

11. **notes.md / project.py length — LOW.** notes.md is ≈8,250 words (guide 3,500–7,000) and project.py 665 lines (guide 150–400). Driven by the brief's topic list and 12-hour weight; previous module reviews took the same position. Not fixed by design.

## Phase 2 resolutions

Before each change the problem was re-derived or recomputed (scratch script /tmp/claude-0/st4_check.py; the ridge boundary x = +3.5 confirmed by an independent BFGS maximization of the penalized log-likelihood).

1. FIXED: notes §3.4 now reads "boundary at x = −β₀/β₁ = 3.922/1.121 = 3.5"; solutions B3(c) now reads "the boundary x = −β₀/β₁ = 2.9302/0.8372 = 3.5 sits midway between x = 3 and x = 4". project.py already asserted β₀/β₁ = −3.5 and the ridge values; it now also asserts (−3.9221, 1.1206) in 5 passes and the first pass (−2.9302, 0.8372).
2. FIXED: 45 tight asserts added across Parts 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, the exercise checks, Your turn (C1, C2, C3) and the B8 reference, each at the precision quoted in notes.md/solutions.md (0.6 × the last quoted digit, or rtol 1–6 % for values quoted to 2–3 significant figures). Three loose asserts were replaced by tight ones (c = 1.701; probit/logit 0.998/1.701/1.705). `se_obs` was pulled out of a print into a variable so C1(c)'s observed-information standard errors can be asserted.
3. FIXED: "the gradient is u − λI₀β ... (u − λI₀β); add and subtract λI₀β".
4. FIXED: "∂p/∂x_j = β_j p(1−p), so |∂p/∂x_j| ≤ |β_j|/4".
5. FIXED: "Newton's method (LA-14 (3.23)) ... local convergence is quadratic (the rate is proved in ML-2; the trace below shows it)".
6. FIXED: "ST-1 Theorem 3.28".
7. FIXED: "w_i from (3.6) and z_i from (3.10)".
8. FIXED: Part 10 → "(notes 3.8, Example 4.6)", Part 11 → "(notes 3.5-3.6)", Part 12 → "(notes 3.6)", comment → "-1/(2 mu^2)".
9. FIXED: A5(b) now names the pass-1 overshoot (1.39 → 1.61, the undamped-Newton remark of §3.4) and says the listed errors are those after passes 2 to 6; project.py asserts both.
10. FIXED: `roc_curve` falls back to `np.trapz` when `np.trapezoid` is absent.
11. NOT FIXED: notes ≈8,250 words and project.py 710 lines (after the added asserts) are forced by the brief's topic and project lists for a 12-hour module; cutting would remove brief-required content.

Re-ran `MPLBACKEND=Agg timeout 300 python3 project.py` after all fixes: exit 0 in 28 s, every printed value byte-identical to the pre-fix run, all asserts pass. Edits: 5 in notes.md, 2 in solutions.md, 1 in exercises.md, 53 lines changed in project.py.

STATUS: review complete; 11 findings (0 HIGH, 2 MED, 9 LOW); 10 fixed, 1 not fixed (item 11, notes/project length, by design); project.py exits 0.
