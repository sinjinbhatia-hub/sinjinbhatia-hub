# Review: ST-2 Intervals and tests

Reviewer: adversarial pass (independent). Module folder: curriculum/03-statistical-modeling/02-intervals-and-tests/.
Status line is at the bottom. Items are numbered in order found; each has file, location, severity (HIGH/MED/LOW), what is wrong, fix, and resolution.

## Phase 1 findings

### Verification performed
- Ran `MPLBACKEND=Agg timeout 300 python3 project.py`: exit 0 in 10 s, 6 figures written. Imports: itertools, pathlib, matplotlib, numpy, sympy, scipy only.
- Recomputed every number in notes §2–§4 and every solution independently (own enumeration/formulas, not the project's functions) in the scratchpad script st2_check.py (numpy/scipy/sympy). All agree: §2 (t₄ = 2.776, 0.95²⁰ = 0.3585), §3.2 (replication-mean 0.834; Wilson quadratic coefficients and discriminant by SymPy), §3.4 ((3.11) = 0.636 at n = 50; permutation 16/924 = 0.0173, pooled t p = 0.0113), §3.6 (Poisson trinity 14.650/11.130/17.067 at x̄ = 92/30), §3.8 (z = 3 → p = 0.0027; replication at p = 0.05 → 0.50), §3.9 ((3.21) by SymPy; winner's curse s = 0.2236, c = 0.4383, a = 1.0655, power 0.143, E = 0.5528 = 2.76×; ∫ₐ^∞ zφ = φ(a) by SymPy), §3.10 (1.35×10⁻³, 0.126, 2.87×10⁻⁷, 0.64/0.23), §3.11 (B₀₁ formula by SymPy; 0.1936 … 43.94; n* = 510.7; p = 0.0124), Ex 4.1 (8.660, −2.121, −3.674; 18), Ex 4.2 (98.06, 0.01564, 0.0017378, 0.04169, 0.01318, 2.2622, 0.0298, [9.7762, 9.8358], 19.023/2.700, [0.0287, 0.0761]), Ex 4.3 (Wald [−0.0315, 0.2315], Wilson [0.0279, 0.3010], 0.1921/0.1645/0.1366; exact coverage 0.876/0.957 at p = 0.10, 0.639/0.925 at p = 0.05 with the stated x-sets {1,2,3,4} and {0,1,2}; Wald(5,20) = [0.060, 0.440]; grid means 0.858/0.954 at n = 20), Ex 4.4 (ideal bootstrap sd 0.7977; atoms 0.207/0.183/0.183; s/√n = 1.181), Ex 4.5 (0.8224, k = 3.632, 0.6388), Ex 4.6 (1.9723, 3.890, 1.4322, 0.1659, 3.592, 0.0800/0.0580), Ex 4.7 (E = 24/26, X² = 5.769, G = 5.826, p = 0.0163/0.0158, z = 2.402, CI [0.050, 0.430], OR 2.667), Ex 4.8 (0.8416, 5.603, 31.40, z-power 0.807, t-power 0.782 at 32, 0.795 at 33, 0.808 at 34), Ex 4.9 (2/3/5 with own Bonferroni/Holm/BH), A1 (2.1314, [49.742, 54.858], [49.850, 54.750], 27.488/6.262, [12.573, 55.19], [3.546, 7.429]), A2 ([0, 0.2775], 0.13877, 0.3487), A3 (0.0547, 0.3828, 0.2335/0.1211/0.0282, 0.0107, 0.1493), A4 (1/1/4; BH thresholds), A6 (1.2816, 116.75, 233.5), B2 (123.0, 0.06504, 6.908/28.845, [0.02808, 0.11726], [8.53, 35.61], 0.037/0.052), B3 (26.296, 131.48, 0.0772, 13.148, 0.6619), B4 (−0.43016, 1.7175, 2.3112, 2.3112, 0.9777, p 0.190/0.128/0.323, exact 0.154; I₁(θ) = θ⁻² by SymPy), B5 (27 triples: 243, 1107 + 1458 = 2565, 14/9, 1.2472 vs 1.5275, √(2/3)), B7 (1.6449/1.7507/3.0902, 0.8038, 3.496, 99.5, 0.222; g(u) by SymPy), C1 (the v-integral reproduces the t density by SymPy, k = 1 Cauchy), C4 (as §3.11).
- Step-by-step read of every derivation in §3 (Theorems 3.1–3.13, Wilson's quadratic, the bootstrap-max atom, the basic interval, Wilks one-parameter and k-parameter, the trinity algebra, the Pearson expansion, the ANOVA decomposition, (3.21), (3.22), the BH proof) and every solution (B1, B2(a), B3, B4, B5(a), B6, B7(a)(b), C1, C2, C4): signs, factors of 2 and π, index ranges, hypotheses all check. Theorem 3.8 step 1 checked by SymPy. Exercise/solution numbering A1–A6, B1–B8, C1–C4 matches; no "similarly"/"left to the reader"/"it can be shown"; no banned filler words.
- Cross-references checked against the actual files: LA-1 Example 4.5 (projection onto a line) ✓; LA-5 Theorems 3.2 (Pythagoras), 3.12 (orthogonal projections), 3.19 (Gram–Schmidt) ✓; LA-8 Theorem 3.11 (spectral) ✓; F4 Theorems 3.6 (Bayes), 3.14 (variance of a sum), 3.15(b) (universality of the Uniform), 3.19(c) (Γ(½)) ✓; F5 Theorems 3.6 (change of variables), 3.8 (sum of Exponentials), 3.22, 3.23 (MVN), 3.35 (CLT), §3.5 (χ² = Gamma(k/2, ½)) ✓; ST-1 I₁, regularity (R), Theorem 3.19 asymptotic normality, Uniform MLE X₍ₙ₎ ✓; ST-3/ST-4 exist; ST-5, ST-7, ML-1 exist in ROADMAP with the named topics (Bayes factors; unit roots; validation). Three problems: items 11, 12, 13.
- Template: §0–§10 headings present; conventions box; 8 Try-it-first / 8 What-just-happened; 9 worked examples; A1–A6, B1–B8, C1–C4 (18); B8 and C3 extend project.py; §8 specific; §9 has 12 questions; §10 matches the brief's six sources. Brief topics all present (sampling distributions by projection, pivots, Wald/Wilson with enumeration, Welch, duality, bootstrap se/intervals/max failure/jackknife/permutation, NP lemma, UMP for Normal mean, Wilks, trinity, χ² GOF/independence, ANOVA by projection, p-values Uniform, misreadings list, power/sample size, winner's curse, Bonferroni/Holm/BH with proof, look-elsewhere, p-hacking, pre-registration, Bayes factor/Jeffreys–Lindley, reporting checklist). Project covers every item in the brief's list.
- History: Gosset 1908 (Macdonell's 3,000 criminals, 750 samples of 4), Fisher's 1912 geometric proof, Fisher 1925 (P = 0.05 "convenient"), Fisher 1956 quotation, Neyman–Pearson 1933, Neyman 1937, Efron 1979, Ioannidis 2005, Open Science Collaboration 2015 (97%/36%), Wasserstein–Lazar 2016, Wilson 1927, Welch 1947/Satterthwaite, Brown–Cai–DasGupta 2001, Efron 1987 (BCa), Quenouille 1949/Tukey 1958, Fisher 1935, Wilks 1938, Pearson 1900, Greenland et al. 2016 (25 misinterpretations), Gelman–Carlin 2014, Holm 1979, Benjamini–Hochberg 1995, Benjamini–Yekutieli 2001, Lyons 2013, 2012 Higgs global significances, Simmons–Nelson–Simonsohn 2011 (60.7%), Gelman–Loken 2013, Jeffreys 1939, Lindley 1957: all consistent with my knowledge. One soft spot: item 17.

### Findings

1. **project.py Part 1 and notes.md §3.1 (What just happened), §7 — MED (weak assert; missing hypothesis).** Notes claim corr(X̄, S²) ≈ 0.68 for Exponential data; the script asserts only `r_exp > 0.3`. The value also depends on n (exact: Cov = μ₃/n, Var S² = μ₄/n − σ⁴(n−3)/(n(n−1)), giving 0.686 at n = 5, 0.698 at n = 10, → 1/√2), and the notes do not say n = 5. Fix: assert |r_exp − 0.686| < 0.02 and write "at n = 5" in the notes.

2. **project.py Part 2 and notes.md Example 4.2 — MED (weak/missing asserts).** Exact (deterministic) coverages quoted in Example 4.3 — Wald 0.639 and Wilson 0.925 at p = 0.05, 0.876/0.957 at p = 0.10, grid means 0.858/0.954 (n = 20) and Wald 0.931 (n = 100) — are asserted only as `< 0.70`, `> 0.90` or not at all. The "z with S" coverage 0.912 at n = 10 (Example 4.2) is printed, not asserted; its exact value P(|t₉| ≤ 1.960) = 0.918 is one line and should be quoted. Fix: tight asserts (atol 5×10⁻⁴ on the exact numbers), assert the MC value within 0.01 of 0.918, add the exact value to Example 4.2.

3. **project.py Part 3 — MED (weak asserts).** Notes quote percentile/basic coverage 0.91/0.89 for the mean, basic 0.82 for the median, 0.000/0.877 for the maximum and 0.948 for the exact pivot; the asserts are `> 0.87`, `> 0.78`, `< 0.90`, `> 0.92`. Fix: assert within 0.015 of the quoted values (seeded, deterministic).

4. **project.py Part 4 — LOW (weak assert).** Permutation p = 0.0173 (= 16/924) and pooled-t p = 0.0113 are only asserted to agree within 0.01. Fix: assert the exact values.

5. **project.py Part 6 — MED (missing asserts).** Example 4.6's t = 1.972, λ = 0.166, −2 log λ = 3.592, p = 0.080 (t₉) and 0.058 (χ²₁), and §3.6's trinity 14.65/11.13/17.07 on x̄ = 3.067, are printed but not asserted. Fix: add asserts.

6. **project.py Part 7 — LOW (missing asserts).** Example 4.7's G = 5.826, p = 0.0163/0.0158, CI [0.050, 0.430] and odds ratio 2.67 are printed, not asserted. Fix: add asserts.

7. **project.py Part 5 — LOW (missing asserts).** Example 4.8's t-test power 0.782 at n = 32 and smallest n = 34 (power 0.808) are computed but not asserted. Fix: add asserts.

8. **project.py Part 8 — LOW (weak assert).** Winner's curse: notes quote 2.7× (simulation) and type-S share 1.8%; the assert allows 2.5–3.3 and does not test the type-S share. Fix: assert |ratio − 2.71| < 0.05 and |share − 0.018| < 0.005.

9. **project.py Part 9 — MED (weak asserts).** §3.10 quotes BH FDR 0.041 (theory 0.040), power 0.72 vs Holm 0.40, BH FWER 0.47; the asserts are qualitative (FWER ≤ 0.065, FDR ≤ 0.06, ordering). Fix: assert each quoted value within 0.02 and Holm/Bonferroni FWER within 0.02 of their printed values.

10. **project.py Part 10 and reference cells — LOW (printed, not asserted; weak asserts).** Solution numbers A1, A2, A3, A6, B2, B3, B4 p-values, B7, C4 are printed only; B8 coverages 0.931/0.927/0.952 asserted as `> 0.88`; C3's 0.198/0.027 asserted as 0.17–0.22 and `< 0.05`. Fix: add asserts to the quoted precision.

11. **solutions.md C4 — MED (wrong cross-reference).** "independent Normals add, F4 §3.14": F4 §3.14 is moment generating functions. The Normal + Normal result is F5 §3.4 ("Normal + Normal", proved by MGF). Fix: cite F5 §3.4.

12. **solutions.md B8 — LOW (cross-reference to content that does not exist).** "ST-3 also covers the residual bootstrap": ST-3's notes do not mention a residual bootstrap. Fix: describe it in one clause and cite Efron and Tibshirani ch. 9 instead of ST-3.

13. **notes.md §3.10, Theorem 3.12 proof — LOW (imprecise cross-reference).** "the union bound (F4 Theorem 3.4)": F4 Theorem 3.4 has only the two-event formula P(A∪B) = P(A)+P(B)−P(A∩B) (item 4); the m₀-event bound needs it iterated. Fix: "(F4 Theorem 3.4(4), P(A∪B) ≤ P(A)+P(B), applied m₀−1 times)".

14. **notes.md §3.4, jackknife paragraph — LOW (math).** "the median, which takes only two values when one point is deleted": true for even n; for odd n the leave-one-out median takes three values (delete below, at, or above the middle). Checked numerically. Fix: "two or three values".

15. **notes.md §3.3, last line — LOW (convention inconsistency).** "is p ≥ 0.05?": the module rejects when p ≤ α (Theorem 3.11: "reject when p ≤ α has size exactly α"), so the non-rejection region is p > 0.05. Fix: "is p > 0.05?".

16. **notes.md §7 — LOW (consistency).** "in nine parts": project.py has ten numbered parts (Part 10 prints the exercise numbers) plus the reference cells. Fix: "ten parts, the tenth printing the exercise numbers".

17. **notes.md §1, last history paragraph — LOW (history, soften).** "the first time it had taken a position on a specific matter of statistical practice": the ASA described the 2016 statement as its first explicit recommendations on a foundational matter of practice; I can confirm the ASA said so, not the stronger unqualified claim. Fix: attribute the claim to the ASA ("by its own account, the first time …").

18. **notes.md and project.py length — LOW (style).** notes.md is ≈9,350 words against STYLE's 3,500–7,000; project.py is 499 lines against 150–400. The brief's topic list (twelve subsections, nine examples, fifteen project items) forces this and the 12-hour weight is the top of the range. NOT FIXED: by design; cutting would remove brief-required content.

## Phase 2 resolutions

Before each fix the problem was re-verified (exact Exponential correlation 2√((n−1)/(8n−6)) = 0.686 at n = 5 by formula and by a 2×10⁵-sample simulation; P(|t₉| ≤ 1.960) = 0.918; leave-one-out median value counts 3/2 for n = 5/6 numerically; F4 §3.14 and ST-3 notes re-read; F4 Theorem 3.4 re-read). Every quoted value the new asserts test was first reproduced independently in st2_check.py (see Phase 1).

1. FIXED: project.py Part 1 computes the exact correlation `r_exp_th = 2*sqrt((n-1)/(8n-6))`, asserts it equals 0.686 to 10⁻³ and the simulated value within 0.02 of it; notes §3.1 now says "at n = 5", gives the exact value and formula source (Cov(X̄,S²) = μ₃/n) and the n→∞ limit 1/√2; §7 says "at n = 5".
2. FIXED: Part 2 asserts the exact coverages [0.876, 0.639] (Wald) and [0.957, 0.925] (Wilson) to 10⁻³, the grid means 0.858/0.954 (n = 20) and 0.931/0.950 (n = 100) to 10⁻³, and the n = 10 "z with S" Monte Carlo coverage within 0.01 of both 0.912 and the exact P(|t₉| ≤ 1.96) = 0.918; Example 4.2 now quotes the exact 0.918 with the simulated 0.912 in parentheses.
3. FIXED: Part 3 asserts percentile 0.914 / basic 0.892 (mean) and basic 0.824 (median) within 0.015, all coverages > 0.80, and the maximum's basic 0.877 and exact-pivot 0.948 within 0.015.
4. FIXED: Part 4 asserts p_perm = 16/924 exactly and the pooled-t p = 0.0113 to 5×10⁻⁵.
5. FIXED: Part 6 asserts t = 1.972, λ = 0.166, −2 log λ = 3.592, p = 0.080 and 0.058 (Example 4.6), and x̄ = 92/30 with the trinity [14.65, 11.13, 17.07] to 5×10⁻³.
6. FIXED: Part 7 asserts G = 5.826, both p-values, z = 2.402, the CI [0.050, 0.430] and the odds ratio 2.667.
7. FIXED: Part 5 asserts n_t = 34, pow_t(32) = 0.782 and pow_t(34) = 0.808 to 10⁻³.
8. FIXED: Part 8 asserts the simulated ratio within 0.05 of 2.71, the formula value 0.553σ to 5×10⁻⁴ (2.76×), and the type-S share within 0.005 of 0.018.
9. FIXED: Part 9 asserts BH FDR 0.041 (±0.005), BH power 0.72, Holm power 0.40 and BH FWER 0.47 (±0.02), keeping the control asserts.
10. FIXED: Part 10 now asserts A1 (t₁₅, both intervals, χ² quantiles, σ interval), A2 (Wilson [0, 0.2775], 0.3487), A3 (four probabilities), A6 (116.75), B2 (λ̂, quantiles, both intervals), B3 (26.296, 131.48, 0.0772, 0.6619), B4 (the three statistics and four p-values), B7 (0.804, 3.50, 99.5, 0.22), C4 (five Bayes factors to 0.1%, n* = 510.7, p = 0.0124); the reference cells assert B8's 0.931/0.927/0.952 within 0.01 and C3's 0.198/0.027 within 0.01/0.005.
11. FIXED: solutions C4 cites F5 §3.4.
12. FIXED: solutions B8 describes the residual bootstrap in one clause and cites Efron and Tibshirani ch. 9; the claim about ST-3 is gone.
13. FIXED: Theorem 3.12's proof cites F4 Theorem 3.4(4) with the two-event inequality written out and "applied m₀−1 times".
14. FIXED: "two values (three when n is odd)".
15. FIXED: "is p > 0.05?".
16. FIXED: "in ten parts (the tenth prints and asserts the numbers behind the exercise solutions)".
17. FIXED: "by its own account the first time it had issued explicit recommendations on a specific matter of statistical practice".
18. NOT FIXED: length (notes ≈9,400 words, project.py 522 lines after the added asserts) is forced by the brief's topic list and the module's 12-hour weight; cutting would remove brief-required content.

Re-ran `MPLBACKEND=Agg timeout 300 python3 project.py` after all fixes: exit 0 (about 10 s), 32 edits applied (8 notes, 2 solutions, 22 project), all asserts pass, every printed value identical to the pre-fix run.

STATUS: review complete; 18 findings (0 HIGH, 6 MED, 12 LOW); 17 fixed, 1 not fixed (item 18, notes/project length, by design); project.py exits 0.
