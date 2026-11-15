# Review: PY-7 Matplotlib and scientific communication

Module: curriculum/00-python/07-matplotlib-and-communication/ (notes.md, exercises.md, solutions.md, project.py, plot_helpers.py)
Brief: curriculum/research/briefs/PY-7.md. Reviewer: independent adversarial pass, 2026-10-26.
Method: ran project.py headless with -W default and timing; extracted and ran every ```python block of notes.md and solutions.md in one session each (isolated copies of the folder) and diffed stated ```text output against real output; recomputed every quoted number in a scratch script; checked cross-references against the referenced modules' files.

Format: one item per finding. Severity: high (wrong output/number/claim, failing or weakened assert, missing required section) / medium (misleading or inconsistent) / low (pedagogy, polish). Status: OPEN / FIXED / NOT FIXED: reason.

## Findings

### 1. [medium] Quiver "page direction" claim is wrong unless the panel is square — notes.md Example 4.7 (text at "One subtlety in step 1", code `du, dv = (c / (2 * L) ...)` / `angles="uv"`), project.py Part 3 (comment + quiver call), solutions.md B3 (same construction), exercises.md B3 ("page-normalised")
Claim: normalising (theta', omega') by the axis ranges gives "the components in axes fractions", and drawing them with `angles="uv"` makes arrows tangent to the curves on the page. Measured (scratch, same figure): the panel is 1.36 times wider than tall, so the arrows are off the true tangent by up to 8.77 deg (median 7.2 deg); solution B3's panel (1.32:1) is off by up to 7.9 deg. Axes fractions equal page direction only for a square panel. Fix: unit-normalise (u, v) in data units and use `angles="xy"` (direction computed through transData; measured error 2e-3 deg); reword the explanation; add the measurement to the notes block (real output) and an assert in project.py; update B3 code/text and the exercise wording. STATUS: FIXED: notes Example 4.7 text, code and stated output (block now prints the measurement: panel 1.33 times wider than tall, axes-fraction arrows off by up to 8.1 deg; real output pasted); project.py Part 3 uses angles="xy" with unit vectors, prints W/H = 1.40 / 9.5 deg and asserts 5 < max discrepancy <= the analytic bound atan(1/sqrt(k)) - atan(sqrt(k)) for k = H/W; solutions B3 code and exercises B3 wording updated; B3 stated output unchanged and re-verified.

### 2. [medium] Cross-reference attributes the Euler–Maclaurin bracket forms to PY-5 §3.10 — notes.md Example 4.2 first paragraph; project.py Part 2 markdown comment
PY-5 §3.10 (fact 2, its eq. (3.2)) derives T - I ~ (b-a)h^2/12 f''(xi) and S - I ~ (b-a)h^4/180 f^(4)(xi) by a Taylor argument; it never mentions Euler–Maclaurin and never writes [f'(b) - f'(a)]. The bracket forms follow by one step (the per-panel sum is a Riemann sum for the integral of f''). The math is right; the attribution is loose. Fix: state the step. STATUS: FIXED: notes Example 4.2 now cites PY-5 §3.10 eq. (3.2) and states the Riemann-sum step to the bracket forms; project.py Part 2 comment matches.

### 3. [medium] Wrong number — solutions.md B2 first paragraph: "h ≈ sqrt(2 eps |f|/|f''|) ≈ 1.5e-8"
With |f| = |f''| this is sqrt(2 eps) = 2.107e-8 (sympy argmin of h/2 F2 + eps F/h confirms sqrt(2 eps F/F2)); 1.5e-8 is sqrt(eps). Fix the number. STATUS: FIXED: solutions B2 now reads sqrt(2 eps) ≈ 2.1e-8 (within a factor of 2 of the grid point 1e-8).

### 4. [low] Wrong digit — solutions.md B4(a) prose: "P(|t_4| < 1.96) = 0.879"
Exact value 2 t.cdf(1.96, 4) - 1 = 0.87845 -> 0.878, which is also what the block prints. Fix prose. STATUS: FIXED: 0.878.

### 5. [medium] `MaxNLocator(5)` described as "at most five ticks" — exercises.md B6; project.py "Your turn" item 1
MaxNLocator(5) bounds the number of *intervals* at 5, so up to six ticks (checked: tick_values(0, 1) returns 6 ticks); the solution correctly asserts <= 6, so exercise and solution disagree. Fix the wording in both places. STATUS: FIXED: exercises B6 and project.py "Your turn" now say at most five tick intervals (MaxNLocator(5), so at most six ticks), consistent with the solution's assert.

### 6. [medium] "project.py exercises every function in it" (notes §3.7, Build 5) is false — project.py never calls `ph.label_end` (only a "Your turn" comment names it) and never tests the ValueError path that §3.4 and §8 rely on ("plot_convergence raises instead")
Fix: in Part 2, call label_end on the slope-guide check figure and assert the annotation's text, colour and anchor; call plot_convergence with a zero error and assert it raises ValueError. STATUS: FIXED: project.py Part 2 calls ph.label_end and asserts the annotation's text, colour and anchor point; calls plot_convergence with a zero error and asserts the ValueError and its message.

### 7. [medium] Asserts weaker than the notes' claims — project.py Part 2 (`abs(c_t - 1) < 1e-4 and abs(c_s - 1) < 2e-3` vs notes Example 4.2 "within 2e-5 and 1.3e-4"; recomputed 1.785e-5 and 1.275e-4), Part 6 (pixel test `< 1e-2` vs notes "within 0.004, one of its 256 steps"; real 3.8e-3 / 1.8e-3), Part 4 (only `tv_kde < tv_gauss` vs notes "0.132 ... 0.019")
Fix: tighten to the claimed values (2e-5, 1.3e-4; 5e-3; ranges around 0.132 and 0.019). STATUS: FIXED: Part 2 asserts 2e-5 / 1.3e-4 (measured 1.785e-5 / 1.275e-4); Part 6 asserts < 5e-3; Part 4 asserts 0.12 < TV(Gaussian) < 0.14 and 0.015 < TV(KDE) < 0.025. Script re-run: exit 0, 21.9 s, no warnings.

### 8. [low] notes.md §6 summary names "a legend decision" among A1–A6; no such exercise exists (A6 is norms and colour-map choice). Fix wording. STATUS: FIXED: "norms and colour-map choice".

### 9. [low] Signature mismatch — exercises.md C2 asks for `label_lines(ax, lines, min_gap=4)`; solutions.md defines `label_lines(ax, lines, names, gap=4.0, fontsize=9)`. Fix: exercise gains `names`, solution renames `gap` -> `min_gap`. STATUS: FIXED: exercise signature label_lines(ax, lines, names, min_gap=4); solution renamed gap -> min_gap (C2 stated output re-verified).

### 10. [low] Colorbar unit for the phase-plane speed is inconsistent and wrong in both places — notes.md Example 4.7 "phase speed (rad per $L/g$)" vs project.py Part 3 "(rad/(L/g)^{1/2})"; hypot(theta', omega') mixes rad/sqrt(L/g) and rad/(L/g). Fix both to "phase-plane speed ... (time in units of sqrt(L/g))", matching the nondimensionalisation stated in the text. STATUS: FIXED: both labels now "phase-plane speed sqrt(theta'^2 + omega'^2) (time in units of sqrt(L/g))".

### 11. [low] Project figures missing labels the notes carry — project.py Part 7 axis labels lack "(dimensionless)" (notes Build 6 has it); Part 11 `size_check` figure has no axis labels (brief: every project figure labelled). Fix locally. STATUS: FIXED: Part 7 labels carry "(dimensionless)"; size_check has x/y labels (tight=False, so the asserted 450x600 pixel size is unchanged).

### 12. [low, no action] solutions.md C3 stated timing column differs from a re-run by a few hundredths of a second (file sizes identical to the byte); the solution header already states timings are machine-dependent. NOT FIXED: not an error.

### 13. [info, no action] notes.md is 9,455 words by `wc -w` (code blocks counted), above the 6,000–8,000 target but under the 9,500 ceiling; nothing is padding. NOT FIXED: within the ceiling.

## Checked and correct (so the next reviewer need not redo them)
- Versions in the conventions box match the environment (Python 3.13.16, NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.2).
- project.py: exit 0 in 21.9 s wall, no warnings under `-W default`; imports only numpy/scipy/matplotlib + stdlib; no `plt.show`, no pyplot state-machine drawing in project.py or plot_helpers.py; jet appears only in the comparison (Part 8, B5); figures saved to figures/ (gitignored).
- All 24 stated outputs in notes.md and 16/17 in solutions.md reproduce exactly when every block is run in one session (the 17th is item 12); no warnings.
- Math recomputed: Euler–Maclaurin leading terms -h^2/6 and h^4/90 (sympy); the residual (T-I)/(-h^2/6) - 1 = h^2/60 = 1.785e-5 at n=96 matches the next E–M term exactly; (3.7) L* constants (841/108, 4/29, (6/29)^3); (4.2) 16p(1-p) = pi(4-pi), sqrt = 1.642; (3.9); lie factor 0.75/0.02184 = 34.3 and z = 2.12; pendulum energies 0.1224/1.8011/2.6450 and separatrix omega = +-2cos(theta/2); damped-oscillator closed forms satisfy x(0)=1, x'(0)=0; angles 14.04/75.96 deg; A3(b) slope -0.0869; A6 TwoSlopeNorm 0.75; B7 mass loss 7.8e-6 at t=8; C1 median |r| ~ 0.6745/sqrt(199) = 0.048; B5 turbo max L* 90.9 at 0.49, one ridge at the centre, ends 12.0/24.5 (so it fails the diverging test only on the ends, as stated); jet max L* 95.9 at 0.64.
- Template: sections 0–10 with the template headings; conventions box; Try-it-first/What-just-happened pairs around every Build; 7 worked examples (brief asks >= 6); A1–A6, B1–B7, C1–C3 (16 problems; template 12–18); B7 extends project.py; §8 specific; §10 lists exactly the brief's six sources; exercises heading style matches PY-5/PY-6.
- Cross-references exist: PY-5 §3.10 and integrators.py (trapezoid/simpson, orders 2 and 4; project cross-checks them at all ten n); PY-6 §3.10 random-walk ensemble (10^4 walkers, 10^3 steps, MSD slope 1); PY-2 csv/json; PY-9 `df.plot(ax=ax)`; F3/F4/F5, ST-1/ST-2, LA-13, CM-1/CM-5, EM-1/EM-2 folders exist.
- Attributions: Tufte's lie factor, Crameri–Shephard–Heron (Nat. Commun. 11, 2020), Okabe–Ito palette (hex values match the canonical eight), Scott's rule n^(-1/5), CIE 1976 L*, Heun/RK4, Van der Pol amplitude 2.0086 for mu=1, pool-adjacent-violators; nothing invented. Honesty section is framed as guidance ("rules of thumb ... not measured laws"; the 2 L* threshold called a convention).
- API claims verified by the runs: `sharex` joins axes; `origin`/`extent` pixel rule (3.4) and the half-pixel failure (9 of 63); bare imshow extent; LogNorm/Normalize values; log-axis non-positive points map ~ -1.1e5 px with no warning; constrained layout + `fig.legend(loc="outside ...")`; PNG size = figsize x dpi.

## Verification after fixes
- `MPLBACKEND=Agg python3 -W default project.py`: exit 0, 21.9 s wall, no warnings.
- Every ```python block of notes.md (25) and solutions.md (18) re-run in one session each: all stated outputs reproduce exactly except the C3 timing column (item 12).

STATUS: review complete; 13 findings (0 high, 6 medium, 7 low/info); 11 fixed, 2 need no action; project.py passes.
