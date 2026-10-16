# LA-8 review (Inner product spaces and the spectral theorem)

Reviewer: independent adversarial pass. Module folder: `01-linear-algebra/08-inner-product-spaces-spectral-theorem/`. Brief: none.
Format: one item per finding — file, location, severity (HIGH / MED / LOW), what is wrong, fix. Status marked in Phase 2.

## Phase 1 — in progress

(Checks and findings are appended below as they are made. If this section still says "in progress", the reviewer was interrupted; re-verify each item before trusting it.)

### Checked so far
- Read all four files in full; `MPLBACKEND=Agg python3 project.py` exits 0 on the unmodified module (~10 s).
- Hand-checked every derivation step in notes §3.1–§3.12 (Thms 3.4–3.22, Weyl corollary, remark after 3.13) and every step in solutions A1–C4. Checked by hand: Legendre norms 8/45, 8/175, 2/5, 2/7; Ex 4.2 characteristic polynomial, Gram–Schmidt, A^{-1} = I − J/4, √A = I + J/3; Ex 4.3 eigenpairs and tip checks; Ex 4.4 projectors, n·σ, |n,±⟩, p± = 3/4, 1/4; Ex 4.5 50 ± 30√2, tan 22.5° = √2 − 1, σq₁ check, Mohr formulas; Ex 4.6 eigenvalues, sum of squares, pivots, Cholesky entries; Prop 3.22 eigenvector verification incl. |n,−⟩; A1–A6 numbers; B3 Hessians; B7(b) a'_pq formula, a'_pp − a_pp = t a_pq, tan 2θ, small root t; C1 interlacing check 4 ≥ 3 ≥ 1 ≥ 1 ≥ 1; C2 polarization signs in both conventions; C3 2λ² − 6λ + 3, (3 ± √3)/2, M-orthogonality.
- Numerical recomputation and cross-reference checks: see below.

### Findings

(Second reviewer, continuing after the first pass was interrupted before any finding was recorded. Every item in "Checked so far" above was re-verified; all of it holds. Additional verification done in this pass:)

- Numerical recomputation in `/tmp/claude-0/la8_check.py` (numpy/scipy/sympy): Ex 4.1 (u_2, u_3, norms 8/45, 8/175, P_2, P_3, 2/5, 2/7), Ex 4.2 (char. poly (μ−1)²(μ+2), eigenvalues 4,1,1, Q orthogonal, A^{-1} = I − J/4, √A = I + J/3), Ex 4.3 (λ = 6, 1; angles 26.57°, −63.43°/116.57°; axis-tip checks; ratio √6 ≈ 2.45), Ex 4.4 (n, n·σ, |n,±⟩, p± = 3/4, 1/4, σ_y projectors), Ex 4.5 (50 ± 30√2 = 92.43/7.57, r = 42.43, θ_p = 22.5°, tan 22.5° = √2 − 1, σq₁ = (85.4, 35.4)), Ex 4.6 (eigenvalues 2 ∓ √2, 2; sampled-sine eigenvectors for k = 1,2,3; pivots 2, 3/2, 4/3; det 4; Cholesky entries), (3.15) and both Bloch eigenvectors symbolically for general θ, φ; A1 (2, √6, √11, 0.2462), A2 (1/4, 1/3, 1/5, √15/4, 14.48°, √(1/15) = 0.258), A3 (√A, A^{-1}, e^A), A4 (all four eigenpairs, −3 ± 2√2, semi-axes 2.41/0.41, asymptotes −2 ± √3), A5 (±i, N² = 0, rank 1, σ_y eigenvalues ±1), A6 (3, 4, 2, 3.8), B3 (gradient, Hessian, critical points, h(±1,0) = ∓2), B6(c) against expm, B7(b) a'_pq and a'_pp formulas symbolically, C1 (4,1,1 vs 3,1), C3 (2λ² − 6λ + 3, (3 ± √3)/2, v₁, v₂, M-orthogonality 0, Euclidean dot 1/2, scipy.linalg.eigh order). All agree with the text.
- Cross-references checked against the actual files: LA-1 Def 3.9, Thm 3.8, Thm 3.11, Ex 4.5 (t* = u·v/‖v‖²), §3.6, Ex B6 (parallelogram), C1 (weighted), C3 (Frobenius; (c) asks for the operator norm); LA-2 Thm 3.15 (one-sided inverse suffices); LA-3 Thm 3.35 and its consequence (ii) (dim(U∩W) ≥ dim U + dim W − n, explicitly flagged there as what Courant–Fischer uses); LA-5 (S^⊥)^⊥ = S, P = QQ^⊤, Gram–Schmidt §3.8, (3.7)/(3.8); LA-6 Thm 3.11 and the complex-matrices remark (det A^† = conj det A, "LA-8 uses this"); LA-7 trace = eigenvalue sum, e^{At} power series, rotation eigenvalues e^{±iθ}; LA-9 §3.11 polar decomposition, Ex C2(c) (singular-value Weyl); LA-10 Cholesky/LDL^⊤/Schur/normal; LA-12 Thm 3.12 (Sylvester), Appendix A (SO(3), SU(2)); LA-13 Weyl, Rayleigh-quotient critical points, shifted QR; LA-14 Thm 3.14 (second-derivative test), ML-2 bullet on Hessian eigenvalues and convexity; F2 (inner product with conjugate first slot, self-adjoint d²/dx², beads preview with the tridiagonal K); F5 (Σ PSD, spectral theorem for the covariance ellipse, symmetric square root from LA-8); CM-6 (Ka = ω²Ma, tridiagonal K in derivation 3, continuum limit); CM-7 (inertia tensor). QM-1/4/5/6, QC-1, ML-2/4/6, ST-3 are not written yet; their roadmap rows match what the notes say about them (QM-5: spherical harmonics, Pauli matrices, Bloch sphere; QM-6: hydrogen; QM-4: uncertainty principle derived; ML-2: convexity; ML-4: kernels; ML-6: PCA; ST-3: regression, ridge; QC-1: Bloch sphere, single-qubit gates).
- CODE: `MPLBACKEND=Agg timeout 300 python3 project.py` exits 0 in ~10 s; imports numpy, scipy, matplotlib, sympy only; seeded; headless; six figures saved. Every assert compared with the §7 claim it is meant to test; none is weaker than the notes' claim except item 7 below (a general formula asserted only at one angle).
- TEMPLATE: §0–§10 with template headings; conventions box (and the same conventions restated at the top of exercises.md and solutions.md); ▶ Try it first / ■ What just happened around Thms 3.4, 3.6, 3.8, 3.9–3.10, 3.11 (both parts), 3.13, 3.14, 3.15, 3.16, 3.17, 3.18, 3.19, 3.20, 3.22; six worked examples including one mechanics example (4.5, stress tensor) and one QM example (4.4); exercises A1–A6, B1–B8, C1–C4 (18 total; B7(d) and B2/B5 extend project.py); solutions mirror the numbering and every part of every exercise is answered; no "similarly"/"left to the reader"/"it can be shown"; §8 is specific; §9 has 12 questions; §10 matches the roadmap row (Axler ch. 6–7, Strang §6.4–6.5) plus 18.06 lectures 25–28 and Trefethen & Bau 24/27/30, all correctly described.
- HISTORY: Fourier 1807 memoir, examiners Lagrange/Laplace/Monge/Lacroix, 1822 book; Euler 1777 orthogonality trick; Cauchy 1829 (secular inequalities); Hermite 1855; Hilbert 1904–1910 and the word "spectrum", Wirtinger 1897 (Hill's equation); Schmidt 1907/1908; Heisenberg 1925, Born–Jordan and Born–Heisenberg–Jordan 1925–26, Dirac/Jordan transformation theory 1926–27, von Neumann 1927 and 1932, Born 1926; Jordan–von Neumann 1935; Jacobi 1846; Forsythe–Henrici 1960. All consistent with my knowledge; nothing to flag.

### Findings

1. **notes.md, Example 4.6 (line 562) — MED — vague cross-link.** "F2 and the PDE modules meet it again." There is no module code for "the PDE modules" (the roadmap has none), and STYLE lists "cross-links that name a field but not a module code" as a send-back item. The matrix does appear in F2 (the beads-on-a-string preview of CM-6, F2 §3) and in CM-6 (derivation 3: the tridiagonal stiffness matrix of N masses on a string, and its continuum limit to the wave equation). Fix: name F2 and CM-6 with those objects.

2. **notes.md, §5 "LA-14 / ML-2" bullet — LOW — misattribution.** "Theorem 3.19: a function is convex iff its Hessian is PSD everywhere." Theorem 3.19 is the second-derivative test at one critical point; the global statement is ML-2's (roadmap: "convexity"), and LA-14 restates the local test as its Theorem 3.14. Fix: say Theorem 3.19 reappears as LA-14 Theorem 3.14 and that ML-2 upgrades it to the global convexity criterion.

3. **notes.md, Example 4.5 last paragraph and §7 item 9 — LOW — misdescription of the code.** "The project draws it from the eigen decomposition." Part 9 computes (σ_n, τ_n) by sweeping θ over all planes with the original σ and only marks the eigenvalues; the eigen decomposition is used for the principal stresses and θ_p, not to draw the circle. Fix: describe what the code does.

4. **notes.md, §7 item 7 — LOW — code mismatch.** "random SPD matrices B^⊤B": the code uses B^⊤B + 0.05 I (to keep the matrices away from singular). Fix: say so.

5. **notes.md, §3.2 (line 165) — LOW — imprecise cross-link.** "(LA-5 §3)" for Gram–Schmidt; LA-5's Gram–Schmidt section is §3.8. Fix: "LA-5 §3.8".

6. **solutions.md, B6(b) — LOW — gap in the argument.** "Two complementary orthogonal projections whose ranges lie in the two eigenspaces must be the spectral projectors" is stated without the reason. Fix: add it — the ranges sum to C² because P₊ + P₋ = I, the eigenspaces E₊, E₋ meet only in 0 (Theorem 3.10), so range(P±) = E± by a dimension count, and an orthogonal projection is determined by its range.

7. **project.py, Part 8 — LOW — assert weaker than the notes' claim.** Example 4.4 states the general Born probabilities p± = ½(1 ± cos θ) for |0⟩ measured along n; the code asserts only p₊ = 3/4 at θ = π/3. Fix: assert the general formula inside the random (θ, φ) loop.

8. **notes.md / project.py — LOW — length over STYLE targets.** notes.md is ~14,500 words (STYLE: 3,500–7,000); project.py is 830 lines (STYLE: 150–400). The content is correct and the brief forbids rewriting correct sections. Not fixed.

No HIGH findings. No wrong numbers, no broken derivation steps, no failing asserts, no missing sections, no incomplete solutions, no convention inconsistencies, no unconfirmable history.

## Phase 1 — complete

## Phase 2 — status

1. FIXED — notes.md Example 4.6 now reads: "You meet it again in F2, as the stiffness matrix of the beads-on-a-string preview of CM-6, and in CM-6 itself, where it is $K$ for $N$ equal masses on a string and its continuum limit is the wave equation." (Both objects verified in F2 §3 and CM-6 derivation 3.)
2. FIXED — notes.md §5 bullet: "Theorem 3.19 reappears as LA-14 Theorem 3.14 (the second-derivative test in matrix-calculus form), and ML-2 upgrades it from one critical point to every point: a $C^2$ function is convex iff its Hessian is PSD everywhere; ..." (rest of the bullet unchanged).
3. FIXED — notes.md Example 4.5 closing line and §7 item 9 now describe what Part 9 does: sweep θ over all planes, assert the circle equation, assert the two ends are the eigenvalues at θ_p = 22.5°.
4. FIXED — notes.md §7 item 7: "random SPD matrices $B^\top B + 0.05I$ (SPD by Theorem 3.15 (iii), shifted to stay away from singular), and on random symmetric matrices that are usually indefinite".
5. FIXED — notes.md §3.2: "(LA-5 §3.8)".
6. FIXED — solutions.md B6(b): added the dimension-count argument (range(P±) ⊆ E±, ranges sum to C² because P₊ + P₋ = I, eigenspaces meet only in 0 by Theorem 3.10, so range(P±) = E±; a projection is determined by its range; H = P₊ − P₋ is (3.10) with equal eigenvalues grouped).
7. FIXED — project.py Part 8: inside the 50-sample random (θ, φ) loop, `assert np.isclose(abs(inner(vp, e1))**2, 0.5*(1 + np.cos(th)))` and the matching p₋ assert. Passes.
8. NOT FIXED: length over STYLE targets (notes ~14,500 words, project 830 lines); fixing would mean rewriting correct sections, which the brief forbids.

Re-run after fixes: `MPLBACKEND=Agg timeout 300 python3 project.py` → exit 0 (~10 s), "All checks passed."

STATUS: 8 findings (0 HIGH, 1 MED, 7 LOW); 7 fixed, 1 deliberately left (length); script passes.
