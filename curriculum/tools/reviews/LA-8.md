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
