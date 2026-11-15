# Review: EM-8 (Potentials, radiation, relativity)

Reviewer notes. One item per finding: file, location, severity, what is wrong, fix. Status marks added in Phase 2.

Phase 1: complete. Phase 2: complete. project.py run 1 (before edits): exit 0 in 5.4 s; run 2 (after edits): exit 0 in 3.9 s, printed output identical.

## Findings

### A. Math (notes §3–§4, solutions) — checked line by line
Verified by hand and/or /tmp/claude-0/em8_check.py (numpy/scipy): Theorem 3.1 steps 1–7 (signs, the cancelling ρ̇/𝓇² terms); Jefimenko (3.4)–(3.5) incl. the curl remark; L–W uniqueness and Jacobian 1 − 𝓇̂·v/c; (3.9); Theorem 3.4 steps 1–6 (𝓇u = cR, (𝓇·u)² = c²R² − v²R²sin²θ, B = v×E/c²); dipole steps 1–8 ((3.13)–(3.16), ∫sin³ = 4/3, 2π, ½); (3.17), (3.18); Larmor steps 1–3; (3.20) γ⁴ for circular motion; t_f = a₀³/(4r_e²c) = 1.556e-11 s; half-wave R = 73.08 Ω; (3.21) τ = 6.266e-24 s, cτ = (2/3)r_e; Lorentz transformation steps 1–6; interval invariance; Theorem 3.7 steps 1–4; Purcell Theorem 3.8 steps 1–7 (γ-addition identity verified algebraically); Theorem 3.9 steps 1–7 (sheet field directions checked with the right-hand rule, rotation y→z, z→−y); Theorem 3.10 steps 1–3; (3.27), Theorem 3.11 components 1–3; (3.29); Theorem 3.12 μ = 0, 1 and G μ = 0, 1; continuity; index-placement remark; (3.31) two components; −□A^μ; (3.32) μ = 1 and μ = 0; (3.33) expansion and the μ = 0 Euler–Lagrange equation; ∂L/∂(∂_νA_μ) = −F^{νμ}/μ₀. Examples 4.1–4.6: every number recomputed and matches (Ex 4.2: t_r = −8.339 ns, 𝓇 = 2.5 m, V = 8.988 V, E = 5.752 V/m; Ex 4.3: 43.90 W, 87.79 Ω, 5.240e-6 W/m², 0.0628 V/m; Ex 4.4: 9.044e22, 4.669e-8 W, 1.556e-11 s, 1.520e-16 s, 2.05e5, 3.3e-6; Ex 4.5: 2.000e-4 T, 3.2044e-17 N, −1.113e-10 C/m, −200.0 V/m, 5.56e-6, 1.362e4 C/m, 7.34e-4 m/s, 8.2e-15; Ex 4.6: 2.502e-3 T). Solutions A1–C3: every derivation step and number checked (A1 10.01/22.00 μs/6.562 km/655.4 m/0.9987 km/0.2178; A2 2.0658/0.8750/1696.0; A3 5.855; A4 1.759e17/1.766e-19 W/2.579 s/3.6e-11 s; A5; A6; B1 (Griffiths Ex. 10.2) incl. ∂_s algebra; B2 (a)–(d); B3 0.19/2.294/12.07; B4 (3.17), 0.06584, 4.335e-3; B5 ṙ algebra, 1.556e-11, 2.048e5, 3.26e-6; B6 γβ² cancellation; B7 B/γ, 1e3 m/s; B8 (a)–(d); B9 ratio 2, (1 + cos²θ)/2, 1.7729e-3 W; C1 6.266e-24, 1.879e-15, e^{-1}; C2 aberration, 5.739°, 0.0573°; C3 all partials).
No mathematical error or wrong number found.

### B. Consistency
Exercise/solution numbering A1–A6, B1–B9, C1–C3 matches; every solution complete (no "similarly"/"left to the reader"/"it can be shown" — grep clean). Conventions box matches usage (metric, F^{0i} = E_i/c, ∂_νF^{μν}, Lorenz gauge, t_r). Cross-refs resolve: EM-6 (3.37) and §3.12; EM-7 §3.13 blue sky; LA-12 (3.16), (3.18), (3.19), §3.8 (Sylvester), §3.15, §3.16; CM-2 §3.9; F6 = 00-foundations/06; EM-5 sheet field; EM-1 4πδ³.

### C. Code
`MPLBACKEND=Agg timeout 300 python3 project.py` exits 0 in 5.4 s; imports numpy/scipy/matplotlib/sympy + stdlib only; seeded; headless; figures saved (figures/ is gitignored). Asserts test the notes' claims at the stated tolerances (Ex 4.2 to 1e-12, Ex 4.5 to 1e-10, all six tensor components for 20 random boosts, both invariants, Lorenz gauge, ∂_νF = ∂_νG = 0 off source, power equal through r = 0.3–50 λ, exponent 4 ± 2e-3, Liénard total, YT asserts).

### Items

1. [notes.md §3.13, after Theorem 3.7] LOW (history/citation). "Taylor §15.9 gives the collision argument" — I cannot confirm the section number (the four-momentum/collision material is in a later section of Taylor ch. 15; §15.9 is, to my recollection, the invariant scalar product/light cone). Fix: cite "Taylor ch. 15 (the section introducing four-momentum)".
   → FIXED — now cites "Taylor ch. 15, in the section that introduces four-momentum".
2. [notes.md §3.2, after Theorem 3.1] LOW (citation). "Griffiths Problem 10.8" — the Lorenz-condition problem is 10.8 in the 3rd edition and 10.10 in the 4th. Fix: "one of Griffiths' ch. 10 problems".
   → FIXED — now "one of Griffiths' ch. 10 problems".
3. [notes.md conventions box and §8 bullet 5] LOW (convention statement). Says Jackson uses "F^{0i} = −E_i (Gaussian or SI)"; in SI that convention reads −E_i/c, and §3.18 itself says −E_i/c. Fix: "F^{0i} = −E_i (Gaussian; −E_i/c in SI)" in both places.
   → FIXED — conventions box, §3.17 warning and §8 bullet all read "F^{0i} = −E_i (Gaussian; −E_i/c in SI)"; §8 also gives the SI form μ₀J^ν.
4. [notes.md §3.13 interval paragraph] LOW (logic of a remark). Sylvester's law is cited as the reason timelike/spacelike agree in every frame; that follows from the invariance of Δs² just shown. Sylvester's content (LA-12 Thm 3.12) is that no coordinates make η definite. Fix: reword the clause.
   → FIXED — invariance of Δs² now carries the timelike/spacelike claim; Sylvester (LA-12 §3.8, Thm 3.12) cited for "η cannot be made definite".
5. [notes.md Theorem 3.7, Theorem 3.12] LOW (template). Both have a Try-it-first but no "■ What just happened" line. Fix: add one sentence after each.
   → FIXED — one-sentence ■ What just happened added after each.
6. [solutions.md B4(a)] LOW (presentation). The intermediate S expression omits the minus sign of B_θ and of φ̂×θ̂ = −r̂ (they cancel; the result is right). Fix: show both signs explicitly.
   → FIXED — both minus signs shown explicitly; result unchanged.
7. [exercises.md B2(d)] LOW (wording). "the naive ∂/∂t of q/(4πε₀𝓇) with 𝓇 treated as fixed" is self-contradictory (a fixed 𝓇 gives zero); the solution answers the intended question (𝓇 changing at the instantaneous rate −𝓇̂·v). Fix: reword the exercise to match.
   → FIXED — exercise now asks exactly what the solution answers (𝓇 changing at −𝓇̂·v, as if t_r advanced at the rate of t).
8. [project.py Parts 6–7] LOW (checkability). The Example 4.3/4.4 numbers quoted in notes.md (43.90 W, 87.79 Ω, 5.240e-6 W/m², 0.0628 V/m, 73.1 Ω, 9.044e22, 4.669e-8 W, 1.556e-11 s, 1.520e-16 s, 2.05e5, 3.3e-6) and the solutions numbers (A1, A2, A4, B3, B4, C1, C2) are printed but not asserted. Fix: add asserts against the quoted values.
   → FIXED — 13 assert blocks added (Ex 4.3, 4.4, 4.5, half-wave 73.1 Ω, A1–A4, B3, B4, B9, C1, C2, YT1 widths), all at ≤1e-3 relative tolerance against the quoted values; script exits 0 in 3.9 s; printed output unchanged.
9. [project.py] LOW (style). 557 lines vs. the 150–400 guideline; runs in 5 s. Not worth cutting the reference solutions the notes point to.
   → NOT FIXED: length guideline only; cutting the reference solutions the notes and solutions point to would make the module worse.
10. [notes.md] LOW (style). ~9,500 words vs. 3,500–7,000 guideline for a 10-hour module whose brief spans retardation, radiation and relativity.
   → NOT FIXED: length guideline only; the brief mandates this scope and no section is padding.
11. [exercises.md] LOW (template). 9 B problems vs. the 5–8 target (18 total, within 12–18).
   → NOT FIXED: 18 problems total is within the 12–18 target; dropping a B exercise would remove content the notes reference (B2, B5, B6, B8 essential; B9 is the code exercise).

### D. Template — §0–§10 with template headings present; conventions box; Try-it/What-just-happened around every key derivation except item 5; six worked examples (brief asks ≥6); A/B/C; §8 specific; §10 matches the brief's sources.
### E. History — Lorenz 1867, Liénard 1898/Wiechert 1900, Larmor 1897, Rutherford 1911, Michelson–Morley 1887, Lorentz 1895–1904, Poincaré, Einstein 1905 (magnet-and-conductor opening), Jefimenko 1966, J. J. Thomson kink construction, Landau–Lifshitz, Morin's book: all consistent with my knowledge. Only item 1 (a section number) could not be confirmed.
### F. Not a programming module; brief has no REVIEW FOCUS section.

Phase 1 complete: 11 items, 0 high, 0 medium, 11 low.

Phase 2 complete: 8 of 11 items fixed; the 3 unfixed are length/count guidelines, not errors. Final: project.py exits 0 with all asserts (including the new ones) passing.

STATUS: DONE — 11 low-severity findings, 0 high, 0 medium; 8 fixed, 3 left as documented style deviations; script passes.
