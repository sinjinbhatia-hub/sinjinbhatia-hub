# Review: ST-7 (03-statistical-modeling/07-time-series-1)

Reviewer: independent adversarial pass. Phase 1 = findings; Phase 2 = fixes marked FIXED / NOT FIXED.
Format: `[severity] file:location — what is wrong — fix`. Severity: high (wrong math/number/assert), medium (incomplete solution, inconsistency, misleading claim), low (pedagogy/style).
(The file existed with this header and an empty Findings section; nothing to re-verify from the interrupted pass.)

Scratch work: `/tmp/claude-0/st7/` (`recompute.py` = independent recomputation of every §4 and solutions number; `run0.txt` = pre-fix script log).

## Method (so far)

- `MPLBACKEND=Agg timeout 300 python3 project.py`: exit 0, 17.6 s. Imports: time, pathlib, numpy, scipy, matplotlib, sympy only.
- Recomputed Examples 4.1–4.8 and solutions A1–A6, B1–B8, C1, C3 independently (numpy/scipy/sympy, not the project's functions): every stated value matches.
- Cross-references checked in the referenced notes: ST-3 Theorem 3.4 (FWL), (3.36) (AIC/BIC); F5 Theorems 3.17, 3.19, 3.21, 3.24, 3.29, §3.12; LA-7 (3.8); LA-8 Theorem 3.16; LA-10 Theorem 3.5; ST-6 (3.18), Example 4.5; CM-1 §3.8 (notation γ, ω₁ matches); ST-8 conventions and (3.11). SM-7, ML-8 not written yet (roadmap entries match the claims).

## Findings

1. **[medium] notes.md §3.3 (■ after Theorem 3.4) and §3.14.** The spectral density is quoted as $\frac{\sigma^2}{2\pi}|\psi(e^{-i\omega})|^2$ and $f(\omega)=\sum_h\gamma(h)e^{-i\omega h}/2\pi$ and attributed to ST-8. ST-8's conventions box and (3.8)–(3.9) use $f(\omega)=\sum_h\gamma(h)e^{-2\pi i\omega h}$ with $\omega$ in cycles per sample and *no* $2\pi$ ("no factors of 2π anywhere"). The forward reference names an object in a convention the target module does not use. Fix: quote ST-8's form in both places.

2. **[medium] notes.md §3.14.** "a complex-root AR(2) into a peak at $\omega$ of (3.15)". ST-8 (3.11) and its ■ show the spectral peak is *not* at the root argument: $\cos2\pi\omega^*=\cos\omega\cdot\frac{1+r^2}{2r}$, so the peak is at a lower frequency (for Example 4.2: spectral period 8.69 vs quasi-period 8). Fix: "a peak near (slightly below) the frequency $\omega/2\pi$ of (3.15)".

3. **[medium] notes.md §5 (ST-6 bullet).** "effective sample size (ST-6 (3.18)) is $n/\sum_h\rho(h)$": the sum is ambiguous (every other $\sum_h$ in the module starts at 1, which would give the wrong formula). ST-6 (3.18) is $S/(1+2\sum_{k\ge1}\rho_k)$. Fix: write $n/(1+2\sum_{h\ge1}\rho(h))$.

4. **[medium] project.py — simulated numbers quoted in notes/solutions are not pinned.** Asserts are much weaker than the prose: random-walk ACF 0.996/0.983/0.965/0.924 (only "> band"); Part 3 means/sds 0.978, 0.988, 0.0622/0.0608, 0.004 (only atol 0.03 / rtol 15%); coverage "between 0.92 and 0.96" (assert allows 0.915–0.985); 0.638/0.930 (assert `< 0.85`, `±0.03`); Ljung–Box 0.052 and "rejects every time"; DF quantiles −1.94/−2.89; power 0.059/0.082/0.129/0.34/0.873; n(φ̂−1) mean −1.772, 5% −7.833; spurious 0.669/0.121/0.052; EG −3.39, 1.9993, 0.149, rejection 1.000; B7 0.657/0.752/0.830/0.876, 0.049; B1 eigenvalues; C3. All values match today's output (run0.txt) but nothing stops drift. Fix: add a `pin()` helper and pin each quoted value to its quoted precision.

5. **[low] project.py Part 0 — worked-example values printed but not asserted.** Ex 4.5 "0.479"; Ex 4.6 interval endpoints, $\sqrt P=1.57465$, half-width 3.08631, ratio 0.945; Ex 4.7 the three terms and total −3.62496; Ex 4.2 ρ(5..8) (covered only through the closed form). Fix: assert them.

6. **[medium] notes.md §3.8 (■ after (3.25)).** "The project checks (3.25) against the 200-dimensional Toeplitz MVN evaluated by Cholesky." The n = 200 check in Part 3 is for the AR(2) (the general `ar_loglik_exact`), and (3.25) itself (the AR(1) formula) is only checked at n = 3 (Example 4.7). Fix: add a literal implementation of (3.25) to the project and assert it against the 200-dimensional MVN for an AR(1) (separate seeded generator so the main stream is unchanged); also assert $\det\Gamma_n=\sigma^{2n}/(1-\phi^2)$ and the tridiagonal precision claimed in the same paragraph.

7. **[low] notes.md §3.5, Theorem 3.6.** Has ▶ Try-it-first but no ■ What-just-happened (template: after each derivation). Fix: add one.

8. **[low] notes.md §3.7 step 1 of ARMA(1,1).** "$\psi_j=\phi\psi_{j-1}=(\phi+\theta)\phi^{j-1}$ for $j\ge1$": the recursion $\psi_j=\phi\psi_{j-1}$ holds only for $j\ge2$ (at $j=1$, $\psi_1=\theta+\phi\psi_0$). The closed form is right. Fix: state the range of each.

9. **[low] notes.md §3.3 (■ after Theorem 3.4).** "Wold's theorem says every stationary, purely nondeterministic process has this form": Wold gives $\sum\psi_j^2<\infty$, not the absolute summability of Definition 3.3, and white (uncorrelated) innovations. Fix: say so in a clause.

10. **[low] notes.md conventions box vs usage.** $\alpha$ is the DF regression constant in (3.33), the adjustment speed in (3.36) and the smoothing constant in §3.13; $\gamma$ is the autocovariance but also a regression coefficient in (3.36)/C2 and CM-1's damping rate in Example 4.2; $\omega$, $r$ (modulus), $\varrho$, $\tau$ are not in the box. Fix: add a sentence to the box listing the overloaded letters; add a parenthetical at (3.36).

11. **[low] notes.md §3.9 truncated predictor.** "$|z_{\min}|^{-n}$": $z_{\min}$ is never defined. Fix: define it (the root of $\theta(z)$ closest to the unit circle).

12. **[low] notes.md §3.4 MA($q$).** (3.10) uses $\theta_0$ without defining it. Fix: "$\theta_0=1$".

13. **[low] solutions.md B6(a).** Equality case of Cauchy–Schwarz: only the case $(x_2,\dots,x_n)=c\,(x_1,\dots,x_{n-1})$ is treated; linear dependence also allows $(x_1,\dots,x_{n-1})=\mathbf 0$ with the other vector arbitrary. (That case is also killed by $x_n=0$, so the conclusion stands, but the proof as written skips it.) Fix: add the case.

14. **[low] solutions.md A3.** "decaying roughly like $1-h/n$ times a constant" is ambiguous ($c(1-h/n)$ or $1-ch/n$?). The notes' own numbers ($0.996,0.983,0.965,0.924$ at $n=2000$) are $1-ch/n$ with $c\approx7$. Fix: write $1-c\,h/n$.

15. **[low] notes.md Example 4.2 (physics reading) and §2, §5.** "AR(2) with complex roots is the sampled noisy underdamped oscillator": true for the AR part (eigenvalues). An oscillator driven by continuous-time white noise and sampled every Δ is an ARMA(2,1) with that AR polynomial; it is a pure AR(2) only if the kicks arrive at the sampling instants. Fix: one clause in Example 4.2.

16. **[low] project.py "Your turn" 2.** "Simulate Yule's model directly: an AR(2) with complex roots, modulus 0.9 and period 11". Those are not Yule's fitted values (his 1927 fit to the sunspot deviations is about $x_t=1.34x_{t-1}-0.65x_{t-2}$, modulus ≈ 0.81, period ≈ 10.6 years). Fix: "a Yule-style sunspot model".

17. **[low] project.py Part 0, Ex 4.3 line.** Two trailing comments run together on the wrong line. Fix: tidy.

18. **[low] notes.md §8 (random-walk ACF).** "Rule: if the ACF has not dropped below about 0.5 by lag 10, difference before identifying anything" is stated as a rule; a stationary AR(1) with φ = 0.95 has ρ(10) = 0.60, and §3.10 says to prefer a test to eyeballing. Fix: call it a rule of thumb that triggers a unit-root test.

