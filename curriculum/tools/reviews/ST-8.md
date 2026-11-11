# Review: ST-8 Time series II (03-statistical-modeling/08-time-series-2)

Reviewer: independent adversarial pass (2026-11-11). Phase 1 = findings; Phase 2 = fixes marked FIXED / NOT FIXED.
Format: `[severity] file:location — what is wrong — fix`. Severity: high (wrong math/number/assert), medium (incomplete solution, inconsistency, misleading or unconfirmable claim, assert weaker than the prose), low (pedagogy/style).
Scratch recomputation: `/tmp/claude-0/st8/recompute.py` (every number in §3–§4 and the solutions), `chk2.py` (exact values of the seeded numbers), `c3c.py` (solution C3(c) code block), run log `run1.txt`.

## Method

- Baseline run `MPLBACKEND=Agg timeout 300 python3 project.py`: exit 0, 9.7 s, "All checks passed". Imports: math, time, numpy, scipy, matplotlib, pathlib only.
- Recomputed independently (numpy/scipy/sympy): AR(2) denominator D(c) and (3.11) symbolically; Examples 4.1–4.8; chi-square quantiles and interval factors; Example 4.6 and B5 by a hand-rolled predict/update (plus Joseph and information forms); the Fibonacci gains of Example 4.7; Riccati roots for A5; GARCH/ARCH kurtosis, the ARMA(1,1)→GARCH rho(1) identity and E(αz²+β)⁴ symbolically; A1–A6, B2, B7–B9, C1 (two-sensor Riccati iteration → gain matrix, column ratio 16); Fisher's exact null formula against a simulation; the 1/Δ² and 1/Δ⁶ side-lobe decay exponents (1.98, 5.95); H_500 = 6.79. All stated numbers agree with my recomputation and with today's script output.
- Read every derivation in §3.1–§3.17 and every solution line by line (signs, 2π factors, transposes, conjugates, index ranges). One wrong inequality found (item 1); the rest hold.
- Cross-references checked in the referenced notes: ST-7 Example 4.2 (φ=(1,−0.5)), (3.11) companion form, (3.15) quasi-period, (3.22) ARMA(1,1) ACF, Theorem 3.6 triangle, Theorem 3.12 SES≡ARIMA(0,1,1), Ljung–Box (3.27); F5 Theorems 3.22/3.24 and the Schur complement; LA-10 Schur complement/Cholesky; LA-13 cancellation; LA-14 Jacobian; ST-5 §3.3 (precisions add) and the Beta–Binomial shrinkage weight; ST-6 §3.1 importance sampling; ST-3 sandwich; ST-1 observed information and Theorem 3.20; CM-1 §3.8 (ω_max² = ω0² − γ²/2); F3 §3.1 convention table, §3.9, §3.16 Hann window. Dirichlet kernel is in F2, not F3 (item 5). ML-4 and ML-6 have no notes yet (forward links prescribed by the brief; not actionable).
- Solution C3(c)'s code block extracted and run against the project's `bootstrap_pf`: runs (4 s), profile maximum at φ = 0.94 on the first data set, flat within one log-likelihood unit over 0.03–0.05; consistent with the prose.
- Template: §0–§10 present with the template headings; conventions box; 8 worked examples; A/B/C = 6/9/3; §8 specific and covers every item in the brief; §10 lists the brief's six sources. Brief has no REVIEW FOCUS section.

## Findings

(Phase 1 list; status appended in Phase 2.)

1. **[high] notes.md §3.3 AR(2) step 5, its What-just-happened, and §9 Q4.** "Since (1+r²)/(2r) ≥ 1, cos 2πω\* ≥ cos θ, so ω\* ≤ θ/2π" and "damping pushes its peak toward lower frequency". Multiplying by a factor ≥ 1 increases cos θ only when cos θ ≥ 0. For φ1 < 0 (θ > π/2) the inequality reverses: φ=(−1,−0.5) has ω\* = 0.3850 > θ/2π = 0.375 (recomputed; also forced by the φ→−φ mirror symmetry the notes state for AR(1)). Correct statement: |cos 2πω\*| ≥ |cos θ|, the peak is pushed away from ω = ¼ toward the nearer end of [0, ½]; below θ/2π iff φ1 > 0. Fix the step, the summary, Q4, and assert both cases in project.py.

2. **[medium] project.py — asserts weaker than the prose, quoted seeded numbers not pinned.** Quoted in notes/exercises/solutions but only printed or loosely bounded: leakage ratios ("invisible"/"three orders of magnitude": assert was `> 50`); Fisher p ≈ 1e-7 (assert `< 1e-3`) and the peak location (never checked, though the brief asks for "periodogram peak vs Fisher's test"); transfer-function error 1.2 % and resonance 0.078 (assert `< 0.1`, `|w−0.08| < 0.003`); tracker RMSE 2.85/1.76/0.91 (only ordering asserted); periodogram variances 1.05/1.11/1.03; coverage 94.8 %; ML estimates 0.58±0.07, 0.90±0.08; "local linear trend's three variances within their standard errors" (assert was 3 se); ARMA log-likelihood −276.614139; EKF 0.054 vs 1.48 vs raw 0.208 ("worse than the raw sensor" not asserted); GARCH QML 0.097/0.856, sample kurtosis 3.58, lag-1 ACF of squares 0.145 (exercise B9), 0.191 over 40 paths; Example 4.3's 5-point MA 0.04; Example 4.6's K and P_{1|1}; B5's numbers; C1(c) gain matrix; C3 −190.51 (sd 0.44) vs −190.35; A-series numbers. Also unasserted claims: Joseph form = (3.23), precisions add, P_{t|n} ⪯ P_{t|t}, smoothed sd peaks mid-gap (§7), H_500 ≈ 6.8 (§8). Fix: add asserts pinning each to today's output (no extra draws from the shared RNG, so nothing shifts).

3. **[medium] notes.md §3.13.** "The connection to ST-7's forecasting folklore is exact … adding the seasonal gives additive Holt–Winters." Exact for the local level (SES) and the local linear trend (Holt). Not exact for the seasonal case: additive Holt–Winters is the steady state of a single-source-of-error seasonal model (equivalently a restricted seasonal ARIMA, which is all ST-7 §3.13 claims); the multi-disturbance trend + dummy-seasonal model here gives close but not identical forecasts. Fix: restrict "exact" to the first two and say what is true of the third.

4. **[medium] notes.md §1 history.** (a) "he used it on earthquake records and, later, sunspots": Schuster's earthquake-periodicity paper (1897) predates the 1898 periodogram paper; the sunspot periodogram is 1906. (b) "many were noise, and nobody could yet say which": Schuster himself derived the chance (exponential) law for a single trial frequency; what was missing was the test for the largest of many ordinates (Fisher 1929, already in §3.6). (c) "Schmidt's group at NASA Ames adapted it … for the Apollo midcourse navigation computer": Schmidt's group developed the extended filter for the Apollo midcourse navigation problem; the onboard implementation was MIT Instrumentation Laboratory's. (d) "runs in every GPS receiver": unverifiable universal. Fix: reword each to what I can confirm.

5. **[medium] notes.md §3.4, §3.6 cross-reference.** "F3's Dirichlet kernel": the Dirichlet kernel is F2 Definition 3.11 (§3.8); F3 has its continuous limit, the sinc (§3.16). Fix the citation in both places.

6. **[medium] notes.md §10 (and §3.1 citation).** Book-content claims I cannot confirm: Shumway & Stoffer ch. 4 "multitapers" and ch. 6 "the diffuse start"; Brockwell & Davis *Introduction* "Herglotz's theorem and the periodogram's asymptotics proved carefully" (the careful proofs are in *Theory and Methods* ch. 4 and 10, which §3.5 already cites); state space is ch. 8 only in the 2nd edition of the *Introduction*. Fix: soften to what the chapters certainly contain; point the Herglotz proof at *Theory and Methods* ch. 4.

7. **[medium] notes.md conventions box.** Says ω, α, β are GARCH parameters "in §3.15–3.17 only", but they are also used that way in §2, Example 4.8, §8, §9 Q12 and exercises A6, B8, B9. The box also omits other reuses the notes rely on: γ_t and ω_t for the seasonal state and its shock (§3.8, §3.13), ν for slope / degrees of freedom / Matérn smoothness, p for dim y / scaled variance / position, B for backshift and bandwidth, θ for MA coefficient / root angle / parameter vector. Fix: correct the scope and list the overloads.

8. **[medium] notes.md Theorem 3.8 and §3.9–3.10: unstated hypotheses.** The RTS gain needs P_{t+1|t}⁻¹, which does not exist for the exactly observed AR companion form of Example 4.5 (P_{t+1|t} = diag(σ², 0) there); the update needs F_t invertible; "precisions add" needs P_{t|t−1} and H invertible. Fix: state each condition at the point of use.

9. **[low] notes.md equation numbering.** Tags jump (3.24) → (3.26); (3.25) does not exist. Fix without renumbering: tag the closed-form steady state as (3.25).

10. **[low] notes.md Try-it-first / What-just-happened coverage.** Theorem 3.7 (innovations) and Theorem 3.9 (likelihood) have no Try-it-first; §3.4's derivations have neither. Fix: add one-line prompts.

11. **[low] notes.md Theorem 3.2 step 5 and Theorem 3.5 Try-it-first.** Step 5 uses E x_t = 0 without saying so; the orthogonality list omits the cos·cos and sin·sin sums for j ≠ k that step 3 of the proof uses. Fix: add both.

12. **[low] notes.md §3.16 step 2.** "ρ(h) = (α+β)ρ(h−1)" holds for h ≥ 2 only. Fix: say so.

13. **[low] notes.md §5 ST-6 bullet.** "forward filtering, backward sampling, which is (3.26) with draws in place of means": FFBS draws from the conditional in step 3 of Theorem 3.8's proof, not from (3.26). Fix the pointer.

14. **[low] notes.md §8 / ST-7 §3.14 mismatch.** ST-7's preview writes the spectral density in the radian convention with 1/2π and says the complex-root AR(2) has "a peak at ω of (3.15)"; here the convention is cycles with no 2π, and the peak is not at θ (item 1). ST-7 is outside this module's folder; add a conversion line to the §8 bullet.

15. **[low] project.py.** `ekf` docstring cites (3.29) for the filter (it is (3.28); (3.29) is the pendulum map). "Your turn" 3 says "resample when the effective sample size drops below 1000 … reusing `bootstrap_pf`", but `bootstrap_pf` resamples every step. Fix both.

16. **[low] notes.md Example 4.7.** "simple exponential smoothing with α=0.618": α is a GARCH parameter in this module's conventions; say "smoothing constant".

17. **[low] notes.md §3.17.** "GARCH is the special case in which next period's variance is a deterministic function of the past": GARCH is not nested in the log-AR(1) SV model just written. Fix wording.

18. **[low] lengths.** notes.md ≈ 9,200 words and project.py 744 lines exceed STYLE's 3,500–7,000 words / 150–400 lines. The brief's topic and project lists require it; runtime is 10 s.

## Status

Phase 1 complete; Phase 2 in progress.
