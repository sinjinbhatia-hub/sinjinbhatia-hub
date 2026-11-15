# Review: EM-6 Electrodynamics (Faraday, Maxwell, energy/momentum, gauge)

Reviewer pass started. Script `project.py` runs clean (exit 0, 11 s) on first run.
Format: ID | file | location | severity | problem | fix | status

## Scope covered
- A. Math: notes §3.1–3.13 every derivation step re-checked (signs, factors, index identities); all numbers in Examples 4.1–4.6 and solutions A1–C3 recomputed independently in /tmp/claude-0/em6_check.py (scipy.constants) — all agree with notes/solutions/project output to the quoted digits.
- B. Consistency: exercises A1–A6, B1–B8, C1–C3 all have complete solutions (no "similarly"/"left to the reader"); conventions box matches usage; every cross-reference checked by grep: EM-1 (3.1),(3.11)–(3.13),(3.19), Thm 3.15/3.16, (P4); EM-2 (3.15)–(3.17), §3.7; EM-4 (3.16), §3.5; EM-5 (3.13),(3.21),(3.22),(3.24),(3.25),(3.26),(3.28),(3.32),(3.42), Thm 3.1/3.2, §3.6/3.7, Ex B1(c), Aharonov–Bohm; CM-3 (4.17), §5; CM-4 §3.10; LA-12 §3.17/3.18 — all exist and name the stated objects. SM-5, QM-4, QC-1, EM-7, EM-8, QM-1 are roadmap codes.
- C. Code: `MPLBACKEND=Agg python3 project.py` exits 0 in 11 s; imports numpy/scipy/matplotlib/sympy only; asserts test the claimed values (Example numbers, M12=M21, Neumann vs elliptic, Poynting I^2R, stress force = Coulomb, gauge invariance, FDTD speed/energy/R,T/flux, group velocity, Courant growth).
- D. Template: §0–§10 with template headings; conventions box; Try-it-first/What-just-happened at every key derivation; 6 worked examples; A/B/C exercises (6/8/3); §8 specific; §10 matches the brief's four sources.
- E. History: Faraday 29 Aug 1831, Henry 1832, Lenz 1834, Neumann 1845, Weber–Kohlrausch 1856, Fizeau 1849, Maxwell 1861–62/1873 and the 1862 quotation, Hertz 1887–88, Poynting 1884, Lorenz 1867, Drude 1900, Dirac 1931, Kerst 1940, Yee 1966, Feynman II §17-2/§17-4/§27-6 — all consistent with my knowledge. One attribution flagged (item 5).
- F. No REVIEW FOCUS section in the brief (not a programming module).

## Findings

1. notes.md | §2 ASCII sketch (lines 23–29) | MEDIUM | The sketch is captioned "induced E circulates clockwise" (correct for B out of page and growing) but the arrows as drawn (top leftward, left side down, bottom rightward, right side up) run counterclockwise. The Lenz intuition picture contradicts its own caption. | Reverse the four arrows so the loop is clockwise. | status: FIXED

2. notes.md | §3.12 Coulomb gauge paragraph (line 298) | LOW | "the instantaneous part of A's source (ε₀∇∂ₜV in (3.33))": the term in (3.33) is c⁻²∇∂ₜV = μ₀ε₀∇∂ₜV; the quoted coefficient is off by μ₀. | Write "the c⁻²∇∂ₜV term of (3.33)". | status: FIXED

3. notes.md | §3.10 What-just-happened (line 262) | LOW | "T_ij da_j is the i-th component of the force on a surface element whose normal is along j" — j is summed, so the sentence conflates T_ij (force per area in direction i on a face with normal j) with the contraction T_ij da_j (force through the element da). | Reword to state both correctly. | status: FIXED

4. project.py | Part 3, Example 4.3(b) labels (M12/M21) | LOW | The notes define Φ₂ = M₂₁I₁ with loop 1 = solenoid, so the flux of the solenoid's B through the small loop is M₂₁ and the sum of loop integrals of the small loop's A is M₁₂; the project names them the other way round (M12 = flux of B_sol). Same number, swapped label. | Swap variable names and print labels to match the notes. | status: FIXED

5. notes.md | §10 Pair-with, Purcell & Morin entry (line 451) | LOW (attribution) | Credits Purcell & Morin ch. 9 with "an electromagnetic wave built from a current sheet switched on at t = 0". I can confirm that construction in Feynman Lectures II ch. 18 §18-4 ("Travelling fields"), not in Purcell ch. 9 (which posits the plane wave and checks it against Maxwell's equations). | Move the current-sheet attribution to the Feynman entry; describe Purcell ch. 7/9 by content I can confirm (the moving rod / moving loop / moving source sequence of ch. 7 and its own reciprocity proof). | status: FIXED

6. project.py | Part 9 markdown cell | LOW | Writes the TE system as "dB_z/dt = -(dE_y/dx - dE_x/dy), eps0 dE_x/dt = dH_z/dy, ..." mixing B_z and H_z while the array evolved is H_z (update uses dt/(mu0 dx)). | Write mu0 dH_z/dt = -(...). | status: FIXED

7. notes.md | Theorem 3.3 proof step 2 (line 118) | LOW | Cites "EM-5 (3.24)" as the thin-wire formula for A; EM-5 (3.24) is the volume form (μ₀/4π)∫J/𝓇 dτ'. | Say "(EM-5 (3.24) with J dτ' → I dl')". | status: FIXED

8. notes.md | whole file | LOW | 8,568 words against the style guide's 3,500–7,000. | Not fixed: this is the 12-hour module and the brief lists every topic present; cutting would drop required derivations. Noted only. | status: NOT FIXED (by design)

Checked and found correct (no action): flux-rule orientation and sign (Thm 3.1), Faraday's disk EMF, (3.9) sign, (3.11), Neumann symmetry, L_sol = 3.158 mH, L_tor = 0.2027 mH, both energy routes, u_B/u_E = 10⁴, displacement-current derivation, polarization current, all four boundary conditions incl. the K_f·(n×ℓ) orientation, wave equation, duality rotation, Dirac condition qΦ_m = nh, Poynting derivation (every sign), stress-tensor derivation (ε-identity contraction, product rule backwards), (3.30)–(3.31), (3.33), gauge theorem and converse, Lorenz-gauge reachability, Lagrangian density (3.38) and its Gauss half, Example 4.3 table percentages (+8%, +51%, +0.5%, +2.0%), 0.99682 finite-solenoid factor, Example 4.6 integrals (1/4d², 1/4d⁴), Yee dispersion relation and Courant argument, discrete-energy conservation (verified by summation by parts, also with spatially varying ε), group-velocity formula, growth factor e^{2η} = 1.4908, C1 betatron condition and numbers, C2 c_a/c_b and α = σZ₀/2, C3 Euler–Lagrange index bookkeeping, A2(c) voltmeter readings ∓ℰ/2.

## Phase 2 result
Fixes applied to notes.md (items 1, 2, 3, 5, 7) and project.py (items 4, 6); exercises.md and solutions.md needed no change. Re-ran `MPLBACKEND=Agg python3 project.py`: exit 0, 9.9 s, all asserts pass, output numbers unchanged (only the Example 4.3(b) labels moved to match the notes' M21/M12 convention).

STATUS: reviewed and fixed — 8 findings (1 medium, 7 low), 7 fixed, 1 not fixed by design (notes length 8.6k words for a 12 h module); project.py exits 0.
