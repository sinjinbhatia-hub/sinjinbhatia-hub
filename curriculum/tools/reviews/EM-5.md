# Review: EM-5 Magnetostatics

Reviewer pass (adversarial, independent). Module: `curriculum/02-physics/b-electromagnetism/05-magnetostatics/`.
Brief: `curriculum/research/briefs/EM-5.md`. Scratch checks: `/tmp/claude-0/em5_check.py`, `/tmp/claude-0/em5_run.txt`.

## What was checked

- **A. Math.** Every derivation in notes §3.1–§3.17 re-derived step by step (cycloid, Biot–Savart straight wire / (3.12) general segment / loop / solenoid, Helmholtz f'' and (3.15)–(3.16), the div/curl of B including the surface term and the delta, Ampère applications, Theorem 3.5 steps 1–5, Theorem 3.6, the solenoid A and its slope jump, the spinning shell steps 1–6 with the Legendre orthogonality, Theorem 3.7 including the loop orientation n̂×t̂, the vector-area identity (3.31) and its companion sign (checked on an explicit circle), the dipole curl in spherical components and the delta coefficient 4π − 4π/3 → 2μ0/3, the ε–δ contraction in (3.35), diamagnetic Δv and χ_m, Curie χ_m, Theorem 3.8 steps 1–4 and the curl theorem, (3.40)–(3.44), (3.45)–(3.47)). f'' and f''' of (R²+u²)^(−3/2) checked with SymPy. Every number in Examples 4.1–4.7 and in every solution recomputed with scipy.constants (`/tmp/claude-0/em5_check.py`); all agree to the quoted digits.
- **B. Consistency.** Exercise/solution numbering A1–A5, B1–B8, C1–C3 matches; every part of every exercise has a worked solution with intermediate values; no "similarly"/"left to the reader". Conventions box matches usage (n̂ outward for K_b, above-minus-below jumps, Coulomb gauge, m = Ia). Cross-references verified against the referenced files: EM-1 (P4)–(P6), (3.1), (3.11), (3.13), (3.18)–(3.20), (3.32)–(3.33), (3.36), Theorem 3.17, §3.2 pseudovector, §3.5, §3.11; EM-2 (3.3), §3.7, Example 4.2; EM-3 (3.27), (3.32)–(3.34), Example 4.5; EM-4 (3.9), (3.12), Theorem 3.3; CM-3 (4.17); LA-12 field tensor pointer. SM-3/SM-6/QM-4/QC-1 do not exist yet (forward references named in the brief).
- **C. Code.** `MPLBACKEND=Agg python3 project.py` exits 0 in 12.7 s; imports only numpy/scipy/matplotlib/sympy/time/pathlib. Asserts test the notes' claims at the stated strength ((3.2), |v| const to 1e−12, (3.13), (3.14), (3.15)/(3.16) symbolically, Ampère circulation to 1e−7, div/curl B on grids, (3.28) wall jump, curl A = B, ∮A·dl = Φ, (3.27), (3.40), (3.44), H⊥ jump = M, (3.46) vs 2D FD within 4 %, B3's d = √3R symbolically, B8 force/torque vs (3.34)/(3.35) with O(side²) convergence, C1's volume integral, C2 root, C3 bias-sign prediction). Figures produced; bar-magnet H panel shows H opposing M inside as §3.15 claims.
- **D. Template.** §0–§10 with template headings; conventions box; 8 Try-it-first / 6 What-just-happened; 7 worked examples; A/B/C with 5/8/3 exercises; B8 extends project.py; §8 covers every item in the brief's §8 list; §10 matches the brief.
- **E. History.** Ørsted 1820, Ampère (Sept 1820), Biot–Savart 1/s, Faraday 1831, Maxwell 1861, electrotonic state → A, Heaviside/Hertz, Ehrenberg–Siday 1949, Aharonov–Bohm 1959, Lawrence, Bohr 1911 / van Leeuwen 1919, Weiss 1907, Hopkinson, Curie temperatures (1043/1388/627 K), SI redefinition 20 May 2019, μ0 to better than 1e−9: all consistent with my knowledge. One dating problem, item 1.
- **F.** Brief has no REVIEW FOCUS section (physics module).

## Findings

1. **notes.md §1, line 13 — MEDIUM (history).** "That last guess is the bound current J_b = ∇×M of §3.13, made quantitative a century and a half later." The equivalence of a magnetized body to Ampèrian currents (J_b = ∇×M, K_b = M×n̂) is 19th-century (Thomson/Maxwell; Maxwell's *Treatise*, 1873), not ~1970. Unconfirmable/incorrect dating. Fix: drop the date ("which §3.13 makes quantitative"). **FIXED** (notes.md §1).

2. **solutions.md B3(e), line 59 — LOW (number/consistency).** "B_z'(0) = 0.6413 × 0.5027 = 0.3224 T/m (the integrator gives 0.3224 T/m; the formula 0.3223)" states two different values for the formula. Exact: 48√3/7^{5/2} · μ0I/R² = 0.32235 T/m, printed by project.py as 0.3223; the integrator's central difference prints 0.3224. Fix: write the product as 0.3223 and keep the parenthetical. **FIXED** (solutions.md B3(e): "0.3223 T/m (more precisely 0.32235; the integrator's central difference gives 0.3224)").

3. **solutions.md B1(a), line 36 — LOW (sloppy orientation claim).** "its outward normal agrees with S_1's on S_1 and is opposite to S_2's on S_2" is backwards for the configuration of part (b) (hemisphere S_2 bulges outside the disk S_1, so the outward normal agrees with S_2 and opposes S_1). Which surface agrees depends on which lies outside; the result (difference = ±∮ = 0) is unaffected. Fix: say so. **FIXED** (solutions.md B1(a): orientation stated as "whichever surface lies outside", with the sign made explicit as ±∮ = 0).

4. **project.py Part 5 (`E_surf = sig / 8.8541878128e-12`) — LOW (code hygiene).** ε0 hard-coded although `EPS0` is imported from scipy.constants. Fix: use `EPS0`. **FIXED** (project.py Part 5; printed value unchanged, 1.129e5 V/m).

5. **notes.md length — LOW, out of scope.** 9,478 words vs STYLE's 3,500–7,000 guideline. The module is a 12-hour one and every section is correct and required by the brief; trimming would mean rewriting correct sections. NOT FIXED: left as is.

No mathematical errors, wrong numbers, failing or weakened asserts, missing sections, or incomplete solutions found.

## Status

project.py re-run after fixes: exit 0, 12.7 s, printed output identical apart from timings. **STATUS: 5 findings (0 high, 1 medium, 4 low); 4 fixed; 1 not fixed (notes length, out of scope — correct content required by the brief).**
