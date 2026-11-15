# Review: ST-8 Time series II (03-statistical-modeling/08-time-series-2)

Two passes. Pass 1 (2026-11-11) completed Phase 1 and was interrupted before any fix reached the module (the four files were still as written on 2026-11-09). Pass 2 (2026-11-15) re-verified every pass-1 item independently, added items 19–23, and did Phase 2.
Format: `[severity] file:location — what is wrong — fix`. Severity: high (wrong math/number/assert), medium (incomplete solution, inconsistency, misleading or unconfirmable claim, assert weaker than the prose), low (pedagogy/style).
Scratch recomputation: pass 1 in `/tmp/claude-0/st8/` (`recompute.py`, `chk2.py`, `c3c.py`); pass 2 in `/tmp/claude-0/st8/v2/` (`mycheck.py` symbolic and numeric re-derivations, `exact.py` and `more.py` exact values of every seeded number, `leak.py` and `leak2.py` the leakage claim, `c3c.py` the solution C3(c) code block).

## Method

- Baseline run `MPLBACKEND=Agg timeout 300 python3 project.py`: exit 0, 12 s, "All checks passed". Imports: math, time, numpy, scipy, matplotlib, pathlib only.
- Pass 2 recomputed independently (numpy/scipy/sympy and by hand): the AR(2) denominator D(c), (3.11) and its root form symbolically; the peak against θ/2π for five parameter pairs; Examples 4.1–4.8 by hand; chi-square quantiles; Example 4.6 and B5 (also in Joseph and information form); the Fibonacci gains of Example 4.7; A1–A6, B2, B7–B9; the ARMA(1,1)→GARCH ρ(1) identity and E(αz²+β)², E(αz²+β)⁴ symbolically; the C1(c) gain matrix by iterating the Riccati recursion (column ratio 16); Fisher's exact null formula against 200,000 simulated maxima; side-lobe decay exponents (rectangular −1.98, Hann −5.95 to −5.99); H_500 = 6.793; the bilinear-warped resonance 0.0782. Every stated number agrees with today's script output; exact values of the seeded numbers were read out of the script's namespace (e.g. periodogram variances 1.0505, 1.1053, 1.0257; coverage 0.9477; Fisher p = 1.036e-7).
- Every derivation in §3.1–§3.17 and every solution read line by line (signs, 2π factors, transposes, conjugates, index ranges, hypotheses). One wrong inequality (item 1) and one claim contradicted by the project's own figure (item 19); the rest hold.
- Cross-references checked in the referenced notes: ST-7 Example 4.2, (3.11) companion form, (3.15) quasi-period, (3.22) ARMA(1,1) ACF, Theorem 3.6 triangle, Theorem 3.12 SES≡ARIMA(0,1,1), Ljung–Box (3.27); F5 Theorems 3.22/3.24 and the Schur complement; LA-10 Cholesky/Schur complement; LA-13 cancellation; LA-14 Jacobian; ST-5 "precisions add" and the shrinkage weight; ST-6 §3.1 importance sampling; ST-3 sandwich; ST-1 observed information and Theorem 3.20; CM-1 resonance (ω_max² = ω0² − γ²/2); F3 §3.1 convention table, §3.9, §3.11, §3.16 Hann window. The Dirichlet kernel is F2 Definition 3.11, not F3 (item 5). ST-7 cites ST-8 (3.3), (3.8), (3.11) by number, so no equation here may be renumbered. ML-4 and ML-6 have no notes yet (forward links prescribed by the brief).
- Solution C3(c)'s code block extracted and run against the project's `bootstrap_pf`: runs in 5 s, profile maximum at φ = 0.94, within one log-likelihood unit over 0.94–0.95; consistent with the prose.
- Template: §0–§10 present with the template headings; conventions box; 8 worked examples; A/B/C = 6/9/3; §8 specific and covers every item in the brief; §10 lists the brief's six sources. The brief has no REVIEW FOCUS section.

## Findings

1. **[high] notes.md §3.3 AR(2) step 5, its What-just-happened, and §9 Q4.** "Since (1+r²)/(2r) ≥ 1, cos 2πω\* ≥ cos θ, so ω\* ≤ θ/2π" and "damping pushes its peak toward lower frequency". Multiplying by a factor ≥ 1 increases cos θ only when cos θ ≥ 0. For φ1 < 0 (θ > π/2) the inequality reverses: φ=(−1,−0.5) has ω\* = 0.3850 > θ/2π = 0.375 (pass 2 recomputed: argmax of the spectrum on a 10⁶-point grid agrees). Correct statement: |cos 2πω\*| ≥ |cos θ|, the peak is pushed away from ω = ¼ toward the nearer end of [0, ½]; below θ/2π iff φ1 > 0. Fix the step, the summary, Q4, and assert both cases in project.py.
   - **Status:** FIXED: step 5 now proves |cos 2πω*| ≥ |cos θ| and states both cases (below θ/2π for φ1 > 0, above for φ1 < 0, equal for φ1 = 0, endpoint maximum when the right side leaves (−1,1)); the What-just-happened and §9 Q4 follow; project.py Part 1 asserts both cases (ω* = 0.1150 < 0.125 and 0.3850 > 0.375) and that each peak is further from ¼ than θ/2π.

2. **[medium] project.py — asserts weaker than the prose, quoted seeded numbers not pinned.** Quoted in notes/exercises/solutions but only printed or loosely bounded: Fisher p ≈ 1e-7 (assert `< 1e-3`) and the peak location (never checked, though the brief asks for "periodogram peak vs Fisher's test"); transfer-function error 1.2 % and resonance 0.078 (assert `< 0.1`, `|w−0.08| < 0.003`); tracker RMSE 2.85/1.76/0.91 (only ordering asserted); periodogram variances 1.05/1.11/1.03; coverage 94.8 %; ML estimates 0.58±0.07, 0.90±0.08; "local linear trend's three variances within their standard errors" (assert was 3 se); ARMA log-likelihood −276.614139; EKF 0.054 vs 1.48 vs raw 0.208 ("worse than the raw sensor" not asserted); GARCH QML 0.097/0.856, sample kurtosis 3.58, lag-1 ACF of squares 0.145 (exercise B9), 0.191 over 40 paths; Example 4.3's 5-point MA 0.04; Example 4.6's K and P_{1|1}; Example 4.7's 2/3 and 5/8; B5's numbers; C1(c) gain matrix; C3 −190.51 (sd 0.44) vs −190.35; A-series numbers. Also unasserted claims: Joseph form = (3.23), precisions add, P_{t|n} ⪯ P_{t|t}, smoothed sd peaks mid-gap (§7), H_500 ≈ 6.8 (§8), side lobes ∝ Δ⁻² and Δ⁻⁶ (§3.6). Fix: add asserts pinning each to today's output (no extra draws from the shared RNG, so nothing shifts).
   - **Status:** FIXED: project.py now pins every quoted number to today's output (tolerances at the quoted precision) and asserts each listed claim: Example 4.1 budget, 4.2, 4.3 (5-point MA = 0.04; ∫4sin² = 2), 4.4 (variances, quantiles, coverage 0.948), 4.5 (F = 2.4, 4/3), 4.6 (K, P_{1|1}, Joseph form, precisions add), 4.7 (2/3, 5/8, 13/21, 0.618), 4.8; Fisher peak j = 37 and 5e-8 < p < 2e-7; transfer error 0.012 ± 0.001 and resonance 0.0782; tracker RMSE; P_{t|t} − P_{t|n} ⪰ 0; smoothed sd largest mid-gap; ML estimates and SEs; trend variances within one SE; −276.614139; EKF 0.054 / 1.48 / raw 0.208 and r_l > raw; GARCH 0.145, 3.58, 0.097, 0.856, forecast vs Monte Carlo at every horizon; H_500; side-lobe exponents; all A/B/C solution numbers (A1–A6, B2, B5, B7–B9, C1 incl. K = P_{t|t}ZᵀH⁻¹, C2 incl. posterior variance, C3). No new draws from the shared RNG; every previously printed line is unchanged except the two reworked in items 19 and 22.

3. **[medium] notes.md §3.13.** "The connection to ST-7's forecasting folklore is exact … adding the seasonal gives additive Holt–Winters." Exact for the local level (SES) and the local linear trend (Holt). Not exact for the seasonal case: additive Holt–Winters is the steady state of a single-source-of-error seasonal model (equivalently a restricted seasonal ARIMA, which is all ST-7 §3.13 claims); the multi-disturbance trend + dummy-seasonal model here gives close but not identical forecasts. Fix: restrict "exact" to the first two and say what is true of the third.
   - **Status:** FIXED: "exact" restricted to the local level and local linear trend; the seasonal case described as Holt–Winters-type recursions, close but not identical, with the single-noise-source model named as the exact one.

4. **[medium] notes.md §1 history.** (a) "he used it on earthquake records and, later, sunspots": Schuster's earthquake-periodicity paper (1897) predates the 1898 periodogram paper; the sunspot periodogram is 1906. (b) "many were noise, and nobody could yet say which": Schuster himself derived the chance (exponential) law for a single trial frequency; what was missing was a test for the largest of many ordinates (Fisher 1929, already in §3.6). (c) "Schmidt's group at NASA Ames adapted it … for the Apollo midcourse navigation computer": Schmidt's group developed the extended filter for the Apollo midcourse navigation problem; the onboard implementation was MIT Instrumentation Laboratory's. (d) "runs in every GPS receiver": unverifiable universal. Fix: reword each to what can be confirmed.
   - **Status:** FIXED: (a) "He had already hunted for periods in earthquake records, and in 1906 he turned the periodogram on sunspots"; (b) Schuster's single-ordinate chance law acknowledged, the missing piece named as the test for the largest ordinate (Fisher 1929); (c) "to the Apollo midcourse navigation problem, and a version of it ran in the Apollo guidance computer"; (d) "Variants of the same recursion now run in GPS receivers, phone motion sensors and autopilots."

5. **[medium] notes.md §3.4, §3.6 cross-reference.** "F3's Dirichlet kernel": the Dirichlet kernel is F2 Definition 3.11 (§3.8); F3 has its continuous limit, the sinc (§3.16). Fix the citation in both places.
   - **Status:** FIXED in §3.4 ("F2's Dirichlet kernel (Definition 3.11) rescaled, the discrete counterpart of F3's sinc"). §3.6 names the kernel by (3.12) without attributing it to F3, so it needed no change.

6. **[medium] notes.md §10 (and §3.1 citation).** Book-content claims that cannot be confirmed: Shumway & Stoffer ch. 4 "multitapers" and ch. 6 "the diffuse start"; Brockwell & Davis *Introduction* "Herglotz's theorem and the periodogram's asymptotics proved carefully" (the careful proofs are in *Theory and Methods* ch. 4 and 10, which §3.5 already cites); state space is ch. 8 only in the 2nd edition of the *Introduction*. Fix: soften to what the chapters certainly contain; point the Herglotz proof at *Theory and Methods* ch. 4.
   - **Status:** FIXED: §3.1 cites *Theory and Methods* ch. 4; §10 entries reworded to chapter contents that are certain (S&S: tapering, coherence, ML and EM; B&D *Introduction* 2nd ed. numbering with the proofs pointed at *Theory and Methods* ch. 4 and 10; D&K chapter-by-chapter; Tsay "a section on stochastic volatility").

7. **[medium] notes.md conventions box.** Says ω, α, β are GARCH parameters "in §3.15–3.17 only", but they are also used that way in §2, Example 4.8, §8, §9 Q12 and exercises A6, B8, B9. The box also omits other reuses the notes rely on: γ_t and ω_t for the seasonal state and its shock (§3.8, §3.13), ν for slope / degrees of freedom / Matérn smoothness, p for dim y / scaled variance / position, B for backshift and bandwidth, θ for MA coefficient / root angle / parameter vector. Fix: correct the scope and list the overloads.
   - **Status:** FIXED: the box now gives the true scope of the GARCH letters and lists the other reused letters (γ_t, ω_t, ν, p, B, θ, f, h).

8. **[medium] notes.md Theorem 3.8 and §3.9–3.10: unstated hypotheses.** The RTS gain needs P_{t+1|t}⁻¹, which does not exist for the exactly observed AR companion form of Example 4.5 (pass 2 confirmed: P_{t+1|t} = diag(σ², 0) there, determinant 0); the update needs F_t invertible; "precisions add" needs P_{t|t−1} and H invertible. Fix: state each condition at the point of use.
   - **Status:** FIXED: F_t invertible stated at step 6 of §3.9; "when P_{t|t−1} and H are invertible" added to the precision form; Theorem 3.8 now assumes P_{t+1|t} invertible, with Example 4.5 named as the degenerate case (asserted in project.py: P_{t+1|t} = diag(1, 0)).

9. **[low] notes.md equation numbering.** Tags jump (3.24) → (3.26); (3.25) does not exist. Fix without renumbering: tag the closed-form steady state as (3.25).
   - **Status:** FIXED: the closed-form steady state p = ½(q + √(q²+4q)), K = p/(p+1) is now (3.25); tags run 3.1–3.36 with no gap and nothing renumbered.

10. **[low] notes.md Try-it-first / What-just-happened coverage.** Theorem 3.7 (innovations) and Theorem 3.9 (likelihood) have no Try-it-first; §3.4's derivations have neither. Fix: add one-line prompts.
   - **Status:** FIXED: Try-it-first added to §3.4, Theorem 3.7 and Theorem 3.9; What-just-happened added to §3.4 and after Theorem 3.7 (17 of each now).

11. **[low] notes.md Theorem 3.2 step 5 and Theorem 3.5 Try-it-first.** Step 5 uses E x_t = 0 without saying so; the orthogonality list omits the cos·cos and sin·sin sums for j ≠ k that step 3 of the proof uses. Fix: add both.
   - **Status:** FIXED: step 5 says to take E x_t = 0 and why that loses nothing; the Try-it-first lists all three orthogonality sums and cites F3 Theorem 3.19.

12. **[low] notes.md §3.16 step 2.** "ρ(h) = (α+β)ρ(h−1)" holds for h ≥ 2 only. Fix: say so.
   - **Status:** FIXED: "for h ≥ 2" added.

13. **[low] notes.md §5 ST-6 bullet.** "forward filtering, backward sampling, which is (3.26) with draws in place of means": FFBS draws from the conditional in step 3 of Theorem 3.8's proof, not from (3.26). Fix the pointer.
   - **Status:** FIXED: the ST-6 bullet now describes FFBS as draws from the conditional in step 3 of Theorem 3.8's proof.

14. **[low] notes.md §8 / ST-7 §3.14 mismatch (pass 1).** Pass 1 saw ST-7's preview in the radian convention with the AR(2) peak "at ω of (3.15)". Pass 2: ST-7 has since been corrected by its own review (its lines 109 and 406 now use cycles per sample, no 2π, and say the peak sits slightly below ω/2π for its Example 4.2, which is the φ1 > 0 case). Nothing left to do here.
   - **Status:** NOT FIXED: nothing to fix (resolved in ST-7 by its own review).

15. **[low] project.py.** `ekf` docstring cites (3.29) for the filter (it is (3.28); (3.29) is the pendulum map). "Your turn" 3 says "resample when the effective sample size drops below 1000 … reusing `bootstrap_pf`", but `bootstrap_pf` resamples every step. Fix both.
   - **Status:** FIXED: docstring cites (3.28) for the filter and (3.29) for the pendulum map; "Your turn" 3 now says `bootstrap_pf` resamples every step and gives ESS-triggered resampling as an extension.

16. **[low] notes.md Example 4.7.** "simple exponential smoothing with α=0.618": α is a GARCH parameter in this module's conventions; say "smoothing constant".
   - **Status:** FIXED: "smoothing constant 0.618".

17. **[low] notes.md §3.17.** "GARCH is the special case in which next period's variance is a deterministic function of the past": GARCH is not nested in the log-AR(1) SV model just written. Fix wording.
   - **Status:** FIXED: "GARCH differs in having no separate volatility shock: next period's variance is a deterministic function of past returns, so nothing is latent".

18. **[low] lengths.** notes.md ≈ 9,200 words and project.py 744 lines exceed STYLE's 3,500–7,000 words / 150–400 lines. The brief's topic and project lists require it; runtime is about 12 s.
   - **Status:** NOT FIXED: the brief's topic and project lists set the length; project.py grew to about 850 lines with the pinning asserts and still runs in about 11 s.

19. **[medium] notes.md §3.6, project.py Part 4 (pass 2).** "a sinusoid of amplitude 0.01 … is invisible without a taper". The project's own data and `taper.png` contradict this: the untapered ordinate at ω = 0.3 is 6.2 times the median leakage around it (1.27e-2 against 2.0e-3), a visible one-bin blip; the script's 1.9× compares it with the higher leakage at ω = 0.25. With the same noise and amplitude 0.003 the untapered line is 1.7 times its surroundings (lost) and the Hann-tapered one 966 times. Fix: measure the line against its own surroundings, add the amplitude-0.003 case (same noise vector, no new RNG draws), assert both, and reword the sentence to what the numbers show.
   - **Status:** FIXED: Part 4 measures the line against the median of its own surroundings for amplitudes 0.01 and 0.003 with the same noise vector (RNG stream unchanged): raw 6.2× / Hann 11044×, and raw 1.7× / Hann 966×, all asserted; §3.6 reworded to those numbers; the side-lobe exponents (−1.95, −5.98) are measured, asserted and quoted.

20. **[low] solutions.md C1(a) (pass 2).** Uses a_{1|1} = P_{1|1}ZᵀH⁻¹y without saying where it comes from: it needs K = P_{1|1}ZᵀH⁻¹ (proved only in part (b)) and I − KZ = P_{1|1}P_{1|0}⁻¹ → 0 to drop the prior mean. Fix: add the two lines.
   - **Status:** FIXED: the two missing lines added (a_{1|1} = (I − K_1Z)a_{1|0} + K_1y, K_1 = P_{1|1}ZᵀH⁻¹ from part (b), I − K_1Z = P_{1|1}P_{1|0}⁻¹ = 0).

21. **[low] notes.md §3.12 (pass 2).** "Durbin & Koopman ch. 5 treat the exact diffuse likelihood": ch. 5 is the exact initialization of the filter and smoother; the diffuse log-likelihood itself is in ch. 7. Fix: cite ch. 5 for the exact diffuse initialization.
   - **Status:** FIXED: "Durbin & Koopman ch. 5 give the exact treatment of a diffuse start".

22. **[low] project.py Part 7 / notes.md §3.11 (pass 2).** Raw-sensor RMSE is computed over all 300 steps (2.85) but filtered and smoothed RMSE over t ≥ 20. Fix: same window for all three (raw becomes 2.86), update the sentence in §3.11, pin all three.
   - **Status:** FIXED: all three RMSEs over t ≥ 20 (2.86, 1.76, 0.91), asserted; §3.11 updated.

23. **[low] outside this folder, not editable here.** ST-7 §3.14 says of the AR(2) peak "damping pulls it to a slightly lower frequency (ST-8 (3.11))" without the φ1 > 0 qualifier. True for its own example; the general statement is item 1's.
   - **Status:** NOT FIXED: ST-7 is outside this module's folder; its sentence is true for the example it cites.

24. **[low] notes.md §3.8 ARMA bullet, solutions.md A4 (pass 2).** The state dimension was called r = max(p, q+1) two lines after Definition 3.6 used r for the number of state shocks (R is k×r), and Q, H were not given for this form.
   - **Status:** FIXED: the bullet and solution A4 use k = max(p, q+1), the bullet states Q = σ², H = 0, and the conventions box lists the two meanings of r and the AR-order meaning of p.

25. **[low] project.py Part 1 comment, §7 runtime, pinned numbers and the reader (pass 2).** The Part 1 comment cited (3.4) for the Fourier series of the ACF (it is (3.3)); §7 said "about 10 seconds" (11–12 s here); and once numbers are pinned, a reader who draws from the shared `rng` in the middle of the file would trip later asserts without knowing why.
   - **Status:** FIXED: comment cites (3.3); §7 says about 12 seconds; §7 and the "Your turn" header say where to put new code and why.

## Verification after the fixes

- `MPLBACKEND=Agg timeout 300 python3 project.py`: exit 0 in 11–12 s, 872 lines, imports unchanged (math, time, numpy, scipy, matplotlib, pathlib). Output compared line by line with the baseline run: every previously printed number is identical; the only changed lines are the leakage report (item 19), the tracker RMSE line (item 22) and the Fisher line (now also prints the peak index), plus new lines for the AR(2) pair, H_500, the side-lobe exponents, the null rejection rate and the smoothed sd across the gap.
- notes.md equation tags run 3.1–3.36 with no gap; ST-7's citations of ST-8 (3.3), (3.8), (3.11) still point at the same equations. 17 Try-it-first and 17 What-just-happened blocks.
- Solution C3(c)'s code block re-run after the fixes is unaffected (it uses only `rng` after everything else, and `bootstrap_pf`).

## Status

Complete: 25 items (1 high, 8 medium, 16 low); 22 fixed, 3 not fixed by design (14 nothing left to fix, 18 length set by the brief, 23 lives in ST-7's folder); project.py exits 0.
