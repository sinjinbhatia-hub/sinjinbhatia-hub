# Review: ST-5 Bayesian inference I

Reviewer: independent adversarial pass (2026-10-31). Module: `curriculum/03-statistical-modeling/05-bayesian-inference-1/`.
Method: every derivation in notes §3.1–§3.13 and every solution re-derived by hand; every quoted number recomputed
with scipy in `/tmp/claude-0/st5_check.py`; RNG-dependent quantities extracted from an actual run of `project.py`
(`/tmp/claude-0/st5_extract.py`, numpy 2.5.3 / scipy 1.18.1); cross-references grepped in ST-1..ST-4, F4, F5, LA-5, LA-14.
Baseline: `MPLBACKEND=Agg python3 project.py` exits 0 in 16.6 s, imports numpy/scipy/matplotlib/sympy only.

## Findings

### 1. [HIGH] notes.md Example 4.8, two-parameter paragraph — numbers contradict project.py
Quoted: MAP (−0.742, 1.266), Laplace sd (0.313, 0.416), grid mean (−0.777, 1.367), sd (0.322, 0.434), log evidence −35.622 vs −35.597.
Script (Part 9) actually prints: MAP (0.037, 1.670), Laplace sd (0.307, 0.464), grid mean (0.046, 1.780), sd (0.314, 0.481), log evidence −35.847 vs −35.822.
Cause: the simulated dataset depends on the RNG stream, which has shifted since the prose was written. Fix: replace with the script's values; pin them with asserts (see 6).
STATUS: FIXED: Example 4.8 two-parameter numbers replaced with the run's values (MAP (0.037, 1.670), sd (0.307, 0.464), grid (0.046, 1.780), sd (0.314, 0.481), log evidence −35.847 vs −35.822); pinned by asserts in Part 9.

### 2. [HIGH] notes.md Example 4.9 table — every entry contradicts project.py
Quoted log evidence row −111.10, −70.83, −14.86, −16.18, …; script prints −130.06, −82.02, −16.64, −17.80, −18.12, −18.64, −19.00, −19.44 (all five rows differ).
The sentence after the table ("fit gains 0.13, Occam factor loses 1.45, BF₂₃ = e^{1.32} = 3.7") is wrong for the actual data: fit gains 0.43, Occam loses 1.59, BF₂₃ = e^{1.16} = 3.2.
Qualitative claims (evidence and BIC both peak at degree 2; max log-lik monotone; BIC penalizes more heavily) remain true and are asserted. Fix: regenerate the table and sentence from the run; pin.
STATUS: FIXED: table regenerated from the run (all five rows); sentence now 'fit gains 0.43, Occam loses 1.59, BF₂₃ = e^{1.16} = 3.2'; whole table pinned by an assert in Part 10.

### 3. [HIGH] notes.md Example 4.11 — numbers contradict project.py
Quoted sample skewness 1.703 and Normal-data p_B = 0.134; script prints 2.712 and 0.665 (p_B for the lognormal data is 0 of 4000 in both). Fix: update; pin.
STATUS: FIXED: skewness 2.712, Normal-data p_B = 0.665; pinned in Part 12.

### 4. [HIGH] solutions.md C3 — all numbers contradict project.py Part 14
Quoted MLE (−0.754, 1.294), ℓ̂ = −31.242, −BIC/2 = −35.336, grid −35.597, Laplace −35.622, dropped terms −4.160/1.838/2.067 (sum −0.255), "−35.591 within 0.006 of grid, 0.031 of Laplace".
Script: MLE (0.041, 1.711), ℓ̂ = −31.533, −BIC/2 = −35.628, grid −35.822, Laplace −35.847, dropped −4.198/1.838/2.176 (sum −0.184), −35.812 within 0.010 of grid and 0.035 of Laplace; BIC gap 0.19 (factor ≈1.2), not 0.26. Same RNG cause as 1. Fix: update; pin.
STATUS: FIXED: C3 now reads MLE (0.041, 1.711), ℓ̂ = −31.533, −BIC/2 = −35.628, grid −35.822, Laplace −35.847, gap 0.19 (factor ≈1.2), dropped −4.198/1.838/2.176 (sum −0.184), −35.812 within 0.010 of grid and 0.035 of Laplace; pinned in Part 14.

### 5. [MEDIUM] solutions.md B7(b) — Monte Carlo value stale
"Monte Carlo with 400,000 paired draws gives 0.9496"; the script's Part 14 gives 0.9501 (quadrature 0.9502 is correct). Fix: update to 0.9501; pin loosely.
STATUS: FIXED: 0.9501; pinned (together with P(best) = (0.022, 0.659, 0.319)) in Part 14.

### 6. [MEDIUM] project.py — quoted RNG-dependent numbers printed but never asserted
Parts 9, 10, 12 and 14 print the numbers quoted in Examples 4.8/4.9/4.11 and solution C3/B7 but only check loose qualitative tolerances, so the drift in items 1–5 went undetected. Also deterministic Example 4.1/4.3/4.7/4.8(one-parameter) values are printed, not asserted. Fix: add asserts pinning each quoted value at its quoted precision (|actual − quoted| ≤ half a unit in the last quoted digit), so any future stream change fails loudly instead of silently desynchronizing the prose.
STATUS: FIXED: asserts added pinning every quoted number in Examples 4.1, 4.3, 4.7, 4.8 (both cases), 4.9, 4.11 and solutions B7, C3 at half a unit in the last quoted digit. Script still exits 0 (16 s).

### 7. [LOW] notes.md Example 4.7 — rounding
Predictive variance at x* = 8 is 2.84648, quoted as "2.847"; should be 2.846 (sd 1.69 and "three times larger" are fine).
STATUS: FIXED: 2.846.

### 8. [LOW] notes.md §3.9 lasso example — implicit σ² = 1
"orthonormal design with z = 0.3 and prior scale b = 2 (λ = ½)" uses λ = σ²/b, which equals ½ only with σ² = 1 (the script's `nlp` indeed uses unit noise variance). Fix: state σ² = 1.
STATUS: FIXED: 'with σ² = 1' added.

### 9. [LOW] notes.md §3.12 Jeffreys–Lindley — ambiguous comparator
"sits in a region where H₁'s density is a hundred times thinner": the factor 100 is H₀'s peak density (39.9) over H₁'s density at x̄ (0.40); at x̄ itself the ratio is 13.5. Fix: add "than H₀'s peak".
STATUS: FIXED: 'than H₀'s peak' added.

### 10. [LOW, not fixed] notes.md length
9108 words against the STYLE guideline of 3500–7000. A 12-hour module with a long brief; trimming would mean rewriting correct material, which this review does not do.
STATUS: NOT FIXED: out of scope (would require rewriting correct sections)

## Verified correct (no action)
- Math: Theorems 3.1–3.14 and Lemma 3.4 checked step by step (Beta–Binomial moments and Eve's-law predictive variance; precision-weighted Normal update and its F5 (3.19) form; NIG completion of the square incl. the constant check; the t scale-mixture lemma and both t marginals; reference-prior limit α₀ → −½ reproducing ST-2's t and χ² intervals; Poisson–Gamma / Negative Binomial; exponential-family conjugate update, the boundary-term argument for E[∇A] = χ/ν, the table's Jacobians, and the NIG tie κ₀ = ν, α₀ = (ν+3)/2 with b = (χ₂ − χ₁²/ν)/2; the logit Jacobian (3.21); Jeffreys invariance and all four table rows incl. 1/(2π cosh(η/2)); HPD-shortest proof; calibration theorem; Bernstein–von Mises heuristic; loss minimizers; ridge λ = σ²/τ² against ST-3 (3.38) and lasso λ = σ²/b against ST-3 (3.43) in their actual normalizations; BLR posterior, predictive, leverage form, NIG regression b_n, g-prior shrinkage, Laplace approximation; evidence (3.38) incl. the (2π)^{d/2} bookkeeping, (3.39), BF (3.40), the maximizer (3.41), the e^{z²/2} bound, BIC derivation).
- Numbers: Examples 4.1–4.7, 4.10, 4.12 and Solutions A1–A5, B1–B8, C1, C2 recomputed independently; all agree to the quoted precision (A2 predictive interval 16.4307/29.0550, B2 half-width 0.09890, C2 u* = 49.44 etc.).
- Cross-references resolve to the named objects: ST-1 (3.11)–(3.13), (3.21)–(3.22), Thms 3.19/3.27/3.28, §3.10, prior written π(θ); ST-2 (3.4), (3.23) B₀₁, Example 4.2 (the ten g values), §3.6 Wald; ST-3 (3.16)–(3.17), (3.38)–(3.39), (3.43), Thm 3.5(b), Example 4.1; ST-4 (3.6), (3.13) IRLS with λI₀, §3.3 concavity, Thms 3.4/3.8 (3.31); F4 Thms 3.5–3.7; F5 Thms 3.6, 3.12, 3.22, 3.24 (3.19), 3.26, 3.28, 3.32, Def 3.44; LA-14 Thm 3.13; LA-5 normal equations. ST-6/ML-4/ML-10/SM-4/ST-8/ST-9 are forward references to modules not yet written.
- History: Bayes d. 1761 / Price 1763; Laplace 1774, male births, Saturn's mass; Fisher 1920s, Neyman 1930s; Jeffreys 1939; Metropolis 1953; Gelfand & Smith 1990; Edwards–Lindman–Savage 1963; Berger & Sellke 1987; Kass & Wasserman (1995); Park & Casella Bayesian lasso; Haldane prior; Zellner g-prior; Bernardo reference priors; van der Vaart ch. 10; MacKay ch. 28; Laplace on the sunrise. All confirmable.
- Template: §0–§10 with template headings; conventions box consistent with usage (variance-parametrized Normal, shape–rate Gamma, IG mean b/(a−1), t location–scale, BF₁₀); 6 Try-it-first / 9 What-just-happened; 12 worked examples; exercises 5 A / 8 B / 3 C, numbering matches solutions, all solutions complete (no "similarly"/"left to the reader"); §10 matches the brief's sources; code requirement met (B7, C2, C3).
- No REVIEW FOCUS section in the brief (not a programming module).

## Status

DONE: 10 findings (4 high, 2 medium, 4 low); 9 fixed, 1 not fixed (notes length, out of scope). `MPLBACKEND=Agg python3 project.py` exits 0 in about 16 s with the new pinning asserts.
