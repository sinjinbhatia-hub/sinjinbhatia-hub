# LA-9 review (SVD and low-rank approximation)

Reviewer: independent adversarial pass. Module folder: `01-linear-algebra/09-svd/`. Brief: none.
Format: one item per finding — file, location, severity (HIGH / MED / LOW), what is wrong, fix. Status marked in Phase 2.

## What was checked (Phase 1)

- MATH: notes §3.1–§3.13 every derivation step (Thms 3.2–3.13), Examples 4.1–4.5, and all 17 solutions read line by line. Every numerical value recomputed independently in `/tmp/claude-0/la9_check*.py` (numpy/scipy): Ex 4.1 (eigenvalues 45, 5; u, v; A_1; NumPy's sign pattern), Ex 4.2, Ex 4.3 (uncentered eigenvalues 36±√1000, top eigenvector 35.78°, mean 33.69°), Ex 4.4 (σ = 2.00005000125 / 4.999875e-5, κ = 40002.00, δx, ratios, log10 κ = 4.60), Ex 4.5 (φ, angles 58.28/31.72/26.57°, R, S, V_mech, E eigenvalues ±), B6(c) (303.02, 1.98, 84.3°), B7(c) (κ2 = φ², κ1 = κ∞ = 4), B7(d) shear-100 (κ ≈ 1.0002e4), B8 (k_gd = 24; predicted error 8.806; optimum 24–26 and 8.73–8.78 over seeds; noisy count > 4/√3·τ√N = 24 on the project's draw; RMT algebra x+1/x, (x²−1)/x², 2+3/x², threshold √3 re-derived), C1(c) code (asserts pass for 40 seeds; "seeds 6 and 31" claim confirmed exactly), C3(b) Hilbert table (all 9 rows match to the digit; slope 1.49; 4 log10(1+√2) = 1.53), B3(a) code runs, B3(b) "typical output" magnitudes confirmed over 3 seeds.
- CONSISTENCY: exercise/solution numbering A1–A6, B1–B8, C1–C3 match; every part of every exercise has a complete solution; no "similarly"/"left to the reader"/"it can be shown". Conventions box matches usage (S, λ double duty flagged; A = QS vs F = RU). Cross-references verified against the actual files: LA-8 Thm 3.11/3.12/3.13 (step 1 outer product)/3.16/3.17, (3.11), normal matrices §3.10; LA-3 Ex C3(b), (3.8) dim formula; LA-1 Ex C3; LA-5 (3.7)/(3.8), Gram–Schmidt; LA-7 Thm 3.4 (trace); LA-12 Thm 3.16; LA-13 Golub–Kahan/randomized SVD, Ex B(c) distance to singularity; LA-14 (3.25) and its pseudoinverse form, Ex B2 ridge; ST-1 Thm 3.6; F3 Thm 3.19; F5; CM-6; CM-7 (3.8); EM-3 Green's functions. ST-3, ML-6, ML-9, QM-7, QC-2 do not exist yet; roadmap rows match the descriptions used.
- CODE: `MPLBACKEND=Agg python3 project.py` exits 0 in ~4 s; imports numpy/scipy/matplotlib/sympy only; seeded; headless; figures saved. Asserts compared against the §7 claims one by one.
- TEMPLATE: §0–§10 present with template headings; conventions box; ▶/■ around all 11 theorems; 5 worked examples (one mechanics); 6 A / 8 B / 3 C; §8 specific; §10 matches roadmap row (Strang ch. 7, LALFD I.8–I.9, 18.065 lectures 4–7) plus T&B and Stewart.
- HISTORY: Beltrami 1873 / Jordan 1874 / Sylvester 1889 "canonical multipliers" / Schmidt 1907 / Weyl 1912 / Stewart SIAM Review 35(4) 1993 / Eckart–Young Psychometrika 1936 / Spearman 1904 / Mirsky 1960 / Golub–Kahan 1965 / Golub–Reinsch 1970 / Moore 1920 / Penrose 1955 / Pearson 1901 / Hotelling 1933 / Deerwester 1990 / Marchenko–Pastur 1967 / Vershynin 2018 ch. 4 / Gavish–Donoho IEEE TIT 60 (2014) / Kabsch — all consistent with my knowledge. One date I cannot confirm (item 2).

## Findings

1. **notes.md §3.11, Derivation step 4 (line 496) — MED — wrong cross-reference claim.** Says "LA-8 constructed √M = Σ√λ_i P_i in (3.11) *but did not prove uniqueness*, so here is the proof." LA-8 does prove it: the Remark after Theorem 3.13, step 3 ("Let S be symmetric with S² = A and all eigenvalues nonnegative … S = √A") and the remark's closing line "a positive semidefinite matrix has exactly one positive semidefinite square root. LA-9's polar decomposition relies on this." The proof given here is correct and self-contained; only the attribution is false. Fix: say LA-8 proved it and present this as a second proof (the eigenspace form is what the polar decomposition uses).

2. **notes.md §1 (line 15) — LOW — unconfirmable date.** "The name 'singular value' settled in during the 1930s." I cannot confirm the decade from my own knowledge (the term entered from the integral-equations literature; the exact timing is Stewart's territory). Fix: soften to a statement without a decade.

3. **notes.md §5, ST-3 bullet (line 684) — LOW — over-specific claim about an unwritten module.** "ST-3 lists LA-9 as a prerequisite for exactly these four facts." ST-3 has not been written; the roadmap lists LA-9 as a prerequisite of ST-3 but says nothing about "exactly four facts". Fix: drop "exactly these four facts".

4. **project.py Part 3 (line 377) — LOW — assert weaker than the notes' claim.** `assert k_gd == 24 and abs(k_opt - k_gd) <= 1` allows the measured denoising optimum to be 23–25, while notes §7 item 3 says the optimum k = 24 "is asserted" and solutions B8 says "Part 3 asserts the 24 and the 8.8". The seeded run gives k_opt = 24 with den_err(24) = 8.7547 vs 8.7818 (k = 23) and 8.7583 (k = 25): a margin of 3.6e-3, far above platform-level floating-point differences. Fix: assert `k_opt == 24`.

5. **notes.md / project.py — LOW — length over STYLE targets (not fixed).** notes.md is ~15,200 words (STYLE: 3,500–7,000 "depending on module weight"; this is the 12-hour capstone); project.py is 605 lines (STYLE: 150–400). The content is correct and the instructions say not to rewrite correct sections, so left as is.

No HIGH findings. No wrong numbers, no broken derivation steps, no failing asserts, no missing sections, no incomplete solutions were found.

## Phase 2 — status

1. FIXED — notes.md §3.11 step 4 now reads: LA-8 "proved uniqueness in the remark after its Theorem 3.13, by comparing eigenbases. Here is a second proof, which shows directly that any PSD square root acts as √λ on each eigenspace of M; that is the form the polar decomposition uses." Proof body unchanged.
2. FIXED — notes.md §1: "The name 'singular value' came later, from the integral-equations literature." (decade removed).
3. FIXED — notes.md §5 ST-3 bullet: "ST-3 lists LA-9 as a prerequisite for these facts."
4. FIXED — project.py Part 3: `assert k_gd == 24 and k_opt == 24`. Re-run: exit 0, prints "best rank 24 (clean σ_i > √3·τ√N: 24), error 8.755 vs 12.808".
5. NOT FIXED: length over STYLE targets (notes ~15,200 words, project 605 lines); fixing would mean rewriting correct sections, which the review brief forbids.

Re-run after fixes: `MPLBACKEND=Agg python3 project.py` → exit 0 (~4 s), "All checks passed."

STATUS: 5 findings (0 HIGH, 1 MED, 4 LOW); 4 fixed, 1 deliberately left (length); script passes.
