# EM-1 review (Vector calculus toolkit)

Reviewer: independent adversarial pass. Module folder: `02-physics/b-electromagnetism/01-vector-calculus-toolkit/`. Brief: `research/briefs/EM-1.md`.
Format: one item per finding — file, location, severity (HIGH / MED / LOW), what is wrong, fix. Status marked in Phase 2.

## Phase 1 — complete

(Phase 2 status is marked on each finding below.)

### Checked so far
- Read brief, STYLE, TEMPLATE, notes.md (all sections), exercises.md, solutions.md, project.py in full.
- Ran `MPLBACKEND=Agg timeout 300 python3 project.py` on the untouched module: exit 0 in 9.3 s; imports are numpy/scipy/matplotlib/sympy/itertools/pathlib only.
- Hand-checked every step of notes §3.1–§3.13 (Thm 3.1 proof, BAC-CAB indices, §3.2 transformation laws, Thm 3.3/3.4 face and edge pairings and signs, (P4) and v×(∇×A) index derivations, Thm 3.5 lemma, Thm 3.6, sphere/disk da cross products, Thm 3.8–3.11, spherical scale factors, r̂×θ̂=φ̂, Thm 3.12–3.13 face areas, (3.29), the §3.10 table incl. spherical curl, Thm 3.14 steps 1–6 incl. the sign of da on r=ε, Helmholtz steps 1–4 incl. the surface-term order, (3.35) and its checks, Thm 3.16–3.17, gauge paragraph).
- Hand-checked Examples 4.1–4.7 (all derivatives, face table, curl of (z−y,x,xy), hemisphere integral, erf potential algebra incl. 4π(2π)^{-3/2}=√(2/π), −4π(1+2σ²)^{-3/2} values, U for Ex 4.6, Fourier coefficients in Ex 4.7).
- Hand-checked solutions A1–C4 every line (A2 cross products, A3, A4, A5(c) Gaussian mean 1+3σ², A6 helix integral and the ẑ=∇×A trick, B2(c) weak-curl argument, B3, B4 vector area, B5, B6(d), B7, B8 discrete product-rule algebra, C1(b), C2 bounds, C3 steps 3–4 signs, C4 index gymnastics and the reflection check).
- Cross-reference and numeric recomputation: below.

### Cross-reference and numeric checks (all pass)
- LA-12 §3.17 (transformation law), §3.19 (Levi-Civita), (3.23) cross product, Appendix A.1 (rotation group), Exercise B6 (ε–δ identity): all exist as cited. LA-6 §3.11 (Jacobian): exists. F3 Definition 3.12, Theorem 3.14, Theorem 3.17, (3.8) derivative rule with ξ (k = 2πξ), §3.10 (3.17) distributional derivative: all exist as cited. CM-1 Theorem 3.3, §3.5, Exercise B2 (the vortex): exist. CM-3 Example 4.8 uses the v×(∇×A) identity and cites EM-1. EM-3 §3.1 has the mean-value property. EM-6/7, QM-5/6, ML-10 are not yet written (forward references consistent with the roadmap).
- Recomputed independently in /tmp/claude-0/em1/recompute.py (sympy/scipy): Example 4.2 face table and total 3; Example 4.3 curl (x, 1−y, 2), loop 8π, hemisphere 8π; Example 4.5 −6.8403, −11.1963, −12.1986, −12.4727, 4π(2π)^{-3/2} = √(2/π), ∇²u_σ = −4πg_σ symbolically; Example 4.6 curl 0, U(1,2,3) = −17; Example 4.7 Û_k = 1/4; A1 (2,2,3), 8, 3, (−y,−z,−x); A2 13, (−4,2,4), (−3,−4,6); A3 √154 = 12.4097, 10/√3 = 5.7735; A5(c) −8π, −26.6407, −25.2270; A6 1/2 − 2π = −5.7832, 1/2, curl (0,0,−2); B3 πR⁴, π(R⁴+R²), −πR², R=2 values; B4 edges −1/2 each, flux −3/2; B7 W, W′, ∇·W = −2yz; C1(b) 1.840302 and the a = 10, 100, 1000 values, 2π/3; C3 ∇²u = μ²u, Yukawa pairing −4π; spherical-curl table entries vs Cartesian at a random point.
- Equation tags (3.1)–(3.36) sequential, every referenced tag exists; Theorems 3.1–3.17 sequential and all referenced ones exist; A1–A6, B1–B8, C1–C4 match between exercises.md and solutions.md; no "similarly"/"left to the reader" in solutions.md.
- History (§1, C3 Yukawa 1935, Debye screening, §10 sources incl. Marsden–Tromba ch. 8, Griffiths App. B, Purcell–Morin App. F): every date and attribution confirmed from own knowledge; nothing to soften.
- Template: §0–§10 with template headings; conventions box; 15 Try-it-first and 15 What-just-happened; 7 worked examples (4.6 is a CM-1 force); A/B/C with 6/8/4; B8 and C1(c) extend project.py; §8 matches the brief's honest list; §10 matches the brief's four sources.
- Probed actual assert margins (/tmp/claude-0/em1/probe.py): off-centre enclosing sphere error 8.7e-14, non-enclosing −1e-16 and −2.7e-15; Helmholtz residuals ≤ 3.3e-15 with field scale 0.42; Theorem 3.4 limit errors 2.4e-2, 1.9e-4, 1.9e-6 (ratios 126, 100); grid Laplacian −12.56820 with 99.86 % within 4h.

### Findings

1. **MED — notation collision φ (test function) vs φ (azimuth).** notes.md §3.11: Theorem 3.14 paragraph (line 339), Try-it-first (341), derivation steps 1–4 (344–349) and Example 4.5(b) (473–475); exercises.md C3(b),(c); solutions.md B2(c) steps 1–3 and C3(b),(c); project.py Part 7 docstring/prints. The conventions box fixes φ as the azimuth, and in notes step 3 and solutions B2(c) steps 2–3 one integral contains both φ (test function) and dφ (azimuth); B2(b),(d) also use φ for the azimuth/potential arctan(y/x). Fix: rename the test function to ψ in all these places (ψ is free there; Green's identities' (φ, ψ) are generic functions and are untouched). **→ FIXED (renamed the test function to ψ in notes §3.11 and Example 4.5(b), exercises C3(b),(c), solutions B2(c) and C3(b),(c), project.py Part 7 docstring and print; Green's-identity (φ, ψ) and C2(b)'s φ = ψ = U untouched; notes line 347 now reads ψ … dθ dφ … ψ dΩ)**
2. **LOW — project.py Part 5 assert weaker than the notes.** Notes §3.4 (line 136) say the error in (3.8) "falls by about 100 each time" as the parallelogram shrinks by 10; the assert only checks `errs[0] > errs[1] > errs[2] and errs[2] < 1e-5`. Observed ratios 126 and 100. Fix: assert both successive ratios lie in [50, 200]. **→ FIXED (asserts both successive error ratios in (50, 200); observed 126 and 100)**
3. **LOW — project.py Part 6 assert weaker than Example 4.4.** Notes say the non-enclosing spheres give "0 to eight decimals" and the off-centre enclosing sphere gives 4π; the assert uses `abs(val − want) < 1e-6`. Actual errors ≤ 1e-13. Fix: tolerance 5e-9. **→ FIXED (tolerance 5e-9; observed errors 8.7e-14, 1e-16, 2.7e-15)**
4. **LOW — project.py Part 8 assert vs Example 4.7's "checks every property to 10⁻¹⁴".** Assert is `< 1e-10 * scale` (≈ 4e-11); actual residuals ≤ 3.3e-15. Fix: assert absolute residual < 1e-13 and say "to better than 10⁻¹³" in the notes so the claim is exactly what is asserted. **→ FIXED (assert < 1e-13 absolute, unused `scale` removed; notes Example 4.7 now says "to better than 10⁻¹³"; observed ≤ 3.3e-15)**
5. **LOW — project.py Part 7(c) asserts weaker than Example 4.5's numbers.** Notes quote "−12.568, 99.9 % of it within four grid spacings"; asserts allow |total + 4π| < 1e-2 (does not pin −12.568) and a fraction as low as 98 %. Observed −12.56820 and 0.9986. Fix: assert |total − (−12.568)| < 5e-4 and near/total > 0.998. **→ FIXED (asserts |total + 12.568| < 5e-4 and near/total > 0.998; observed −12.56820, 0.9986)**
6. **LOW — notes.md §3.1 (line 61) and §3.2 item 3 (line 85): "(LA-6)" cites the module without an object.** The triple-product determinant is LA-6 Theorem 3.15(a); the identity ε_lmn R_il R_jm R_kn = (det R) ε_ijk is the Leibniz formula (LA-6 Theorem 3.9) applied to the matrix whose rows are rows i, j, k of R, and Exercise C4(a) proves it here. Fix: cite those. **→ FIXED (now "LA-6 Theorem 3.15(a)" and "the Leibniz formula, LA-6 Theorem 3.9, applied to the matrix whose rows are rows i, j, k of R; Exercise C4(a) proves it")**
7. **LOW — notes.md Example 4.5(b) (line 471): "erf(r/√2σ)" is ambiguous about the grouping** (the intended argument is r/(√2 σ), as E′ and project.py confirm). Fix: write erf(r/(√2σ)) in the two places. **→ FIXED (erf(r/(√2 σ)) in both places)**
8. **LOW — length.** project.py is 686 lines (STYLE: 150–400) and notes.md ≈ 9100 words (STYLE: 3,500–7,000). Every part of the project and every topic of the notes is demanded by the brief (ten project items, ≥ 6 worked examples, Helmholtz, curvilinear, delta, potentials). Not fixed: cutting would remove brief-required content; no local fix exists. **→ NOT FIXED: every part is required by the brief; no local fix. Runtime is 9 s, well under the 60 s limit, so the length costs nothing operationally.**

No HIGH findings. No mathematical error, wrong number, missing section, incomplete solution, convention inconsistency (other than item 1) or dubious attribution found.

## Phase 2 — complete

- All edits confined to the module folder and this file. `MPLBACKEND=Agg timeout 300 python3 project.py` re-run after the fixes: exit 0, ~9 s, every tightened assert passes with the margins noted above.
- Re-verified before changing: each of items 1–7 was confirmed real (notation collision visible on notes line 347 and solutions B2(c) steps 2–3; assert margins probed numerically; LA-6 theorem numbers checked in LA-6 notes.md).

**Status: 8 findings (0 HIGH, 1 MED, 7 LOW); 7 FIXED, 1 NOT FIXED (length, brief-mandated); script exits 0; no mathematical, numerical, historical or completeness errors found.**
