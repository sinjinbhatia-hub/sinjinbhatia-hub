# Review: EM-7 Electromagnetic waves

Reviewer: adversarial pass (Phase 1 in progress). Module folder: curriculum/02-physics/b-electromagnetism/07-electromagnetic-waves/

Baseline: `MPLBACKEND=Agg python3 project.py` exits 0 in ~4 s before any edits.

## Findings (one per item: file, location, severity, problem, fix)

### A. Math (notes §3–§4, solutions) — every derivation stepped through, every number recomputed in /tmp scratch with scipy constants

Verified correct (no action): d'Alembert (3.1–3.2); Lemma 3.2; string junction (3.4) incl. k↔v forms; Theorem 3.3 steps 1–3 (sign of Ampère step, k̂×(k̂×E) identity); (3.7)–(3.8); Jones matrices (3.10), QWP = S, HWP(22.5°) = Hadamard, Stokes (3.11), Bloch ↔ Poincaré axis map; (3.12); normal incidence (3.13) incl. β = μ1n2/(μ2n1); Theorem 3.4; Vandermonde lemma and Snell (3.15); p-Fresnel (3.16) (all three cross products, D⊥ = step 4 verified algebraically); s-Fresnel (3.17) and B2's full derivation; R + T = 1 both polarizations; Brewster (3.18) and the ê_T ∥ k̂_R dipole argument; TIR root choice, (3.19), |r| = 1; conductor (3.20)–(3.22), κ formula, phase lag direction, |β̃| ≈ 3×10⁴, R = 0.999914; Lorentz (3.23)–(3.25) incl. the α factor; Sellmeier/Cauchy (3.26) incl. B = B1λ1²/(2n0) = 0.00307; Larmor/Rayleigh 5.86; plasma numbers (8.98 MHz, 1.64e16 rad/s, 115 nm); (3.29)–(3.30) incl. Re(1/a) = τ0²/|a|², chirp sign, n_g(ω0) = 1 − F/γ²; waveguide (3.31)–(3.33) all six component equations; Theorem 3.6; coax (3.34) and Z = 50.08 Ω; all six worked examples (every intermediate value matches to the printed precision); all solution numbers A1–C3 (B4(c) r_s = −0.1000 − 0.9950i, phase −95.74°; B5(c) 39.73 µs; B8 R = 0.0141, 99.6 nm; C1; C2 48.62°/54.62°, max 45.94° at 51.34°; C3 2.992e-16 s, 1.04e12 yr). Cauchy-fit gap explanation checked: UV-only fit gives B = 0.00325 (expansion's next term), IR tail adds the rest to 0.00350.

1. project.py Part 3, `close(Rs, 0.09202, ...)` — LOW. Reference value is 0.09202; exact R_s = 0.30334² = 0.092013, and notes Example 4.4 / solutions B3 say 0.09201. Passes only because rtol = 2e-4. Fix: reference 0.09201.

2. notes §3.6, circular polarization sentence — LOW. "for 1/√2(1, i) at z = 0, E = (E0/√2)(cos ωt, sin ωt)": E0 appears from nowhere since J was given normalized. Fix: write J = (E0/√2)(1, i).

3. solutions.md C2(b) — LOW/MEDIUM (completeness). The 45° angles and the maximum retardance are taken from the project's numerical scan, with only one hand check; a "find" question should be solvable by hand. Fix: add the quadratic in s = sin²θ for Δ = 45° (roots s = 0.5631, 0.6648 → 48.62°, 54.62°) and the analytic maximum tan(Δ_max/2) = (1 − n²)/(2n) at sin²θ = 2n²/(1 + n²) (→ 45.94° at 51.34°); assert both in project.py.

### B. Consistency

Numbering A1–A5, B1–B8, C1–C3 matches between exercises.md and solutions.md; every solution complete (no "similarly"/"left to the reader"; the only grep hits are the template heading and "simply connected"). Conventions box matches usage (e^{−iωt}, Griffiths r_p, interface geometry, TE_mn indexing). Cross-references checked against the referenced notes: EM-6 (3.1), (3.22), (3.23), (3.24), (3.26), (3.31), §3.8, §7, Example 4.5 — all exist and say what EM-7 says; F3 (3.12) Gaussian pair, (3.23) transfer function 1/(k − mω² + icω), "pure tone" eigenfunction remark; F1 §3.8 sign-convention warning box; F2 §1 string; F6 Example 4.3 Fermat → Snell; CM-1 §3.8; CM-6 Example 4.5 ω(κ); EM-1 Theorem 3.16; EM-3 Theorem 3.3; EM-4 (3.5), (3.18); LA-1 composition, LA-5 projections — all exist. QM-1/2/3/9, QC-1, EM-8 are forward references; EM-8 folder exists, QM/QC folders are empty but ROADMAP.md lists tunnelling in QM-3, path integrals in QM-9, Bloch sphere in QC-1, so the objects named exist in the plan.

4. notes §3.15 last paragraph — LOW. "the retarded potentials of EM-6 (3.37)": EM-6 (3.37) is the pair of Lorenz-gauge inhomogeneous wave equations □V = −ρ/ε0, □A = −μ0J; EM-6 itself says EM-8 solves them with retarded potentials. Fix wording: "the inhomogeneous wave equations EM-6 (3.37), solved by retarded potentials".

5. notes §7 — LOW. "runs in about 3 seconds"; measured 4.0 s. Fix: "about 4 seconds".

6. notes §3.11 fibre paragraph — LOW. "light inside the critical angle bounces down the core by TIR" reads as if TIR happens for incidence angles below θ_c. Fix: "light that meets the core wall beyond the critical angle".

### C. Code

`MPLBACKEND=Agg python3 project.py` exits 0 in 4.0 s; imports numpy, scipy.constants, matplotlib, sympy, pathlib, time only. Asserts checked against the notes' claims: Maxwell/Lemma 3.2/string junction (SymPy), Example 4.1 numbers, QWP = S, HWP = Hadamard, rotating-QWP trajectory formula, Malus, three polarizers = 1/8, Fresnel numbers at 45° (rtol 2e-4), Brewster zero (<1e-14), R + T = 1 at 500 angles both ways (<1e-12), all four interface conditions on assembled fields (<1e-12, incl. beyond θ_c), |r| = 1 beyond θ_c, (3.19) = 1/Im k_z, FDTD r, t, R, T (0.5–1%), chirped-pulse widths FDTD vs FFT vs (3.30) (3%), compress-then-spread ordering, v_g (1%), exact spreading for quadratic k (1e-4), Cauchy fit error < 2e-4 and B within 15% (gap explained in notes), ionosphere f_p, skin depth exact vs (3.21) (1e-6) and vs quoted values, TE fields vs Maxwell + E∥ = 0 (SymPy), cutoffs, v_p v_g = c², AR coating and plasma mirror, C1 group index. FFT sign handling (conj of e^{ikz} on rfft bins) verified to be correct for e^{−iωt} physics.

7. project.py Part 8 SymPy check — LOW. Uses m ≥ 1, n ≥ 0 only, so TE_0n is not covered, and B⊥ = 0 on the walls (which §3.15 step 4 says "holds automatically by (3.31)") is not checked. Fix: run the same check for (m ≥ 0, n ≥ 1) and add the four B⊥ wall checks (measured cost < 1 s).

8. project.py — LOW. Numbers that the notes/solutions state but the script only prints: Example 4.5's R = 0.999914 and cB/E = 3.27e4; Example 4.6's λ_g = 3.953 cm, k = 158.96 rad/m, 4.398 ns/m; §3.13's n(589 nm) = 1.4571, copper ω_p = 1.642e16 rad/s and λ_p = 115 nm; B4(c), B5(c), C2, C3. Fix: add `close`/assert lines for each.

9. project.py length 506 lines vs STYLE's 150–400 — LOW, NOT FIXED: the brief requires eight parts plus three reference "Your turn" solutions and Part 9; trimming would remove required checks.

### D. Template

Sections 0–10 present with template headings; conventions box; 14 Try-it-first and 12 What-just-happened; six worked examples (brief's six); A5/B8/C3 split; B8 requires project.py; §8 specific (six items from the brief); §9 has 12 questions; §10 lists exactly the brief's sources. notes.md is ~9,200 `wc` words (LaTeX-inflated) against the 3,500–7,000 guideline — LOW, NOT FIXED (12-hour module; the brief's topic list fills it).

### E. History

Confirmed from own knowledge: d'Alembert 1747; Newton's prism 1666; Malus 1808; Brewster 1815; Fresnel 1821–23 and the elastic-ether account; Cauchy 1836; Hertz at Karlsruhe 1887–88 (spark gap, zinc-sheet standing waves, pitch prism, polarization); Sellmeier 1871; Rayleigh 1871; Lorentz's electron theory 1890s; Sommerfeld–Brillouin 1914; Jackson §7.10 Kramers–Kronig; Feynman I ch. 31; Griffiths §9.3.3/Fig. 9.15; IKAROS 2010; solar constant 1361 W/m²; Cu σ = 5.96e7 S/m, N = 8.47e28 m⁻³; fibre loss 0.2 dB/km at 1.55 µm; Hecht vs IEEE circular-handedness conventions. Nothing to soften.

### F. REVIEW FOCUS

Brief has no REVIEW FOCUS section (physics module). Not applicable.

## Phase 2 — fixes

1. FIXED — project.py Part 3 reference for R_s is now 0.09201.
2. FIXED — notes §3.6 now writes J = (E0/√2)(1, i) before giving the field.
3. FIXED — solutions C2(b) now derives the two 45° angles from the quadratic (1 + τ)s² − (1 + n²)s + n² = 0 (s = 0.56312, 0.66480 → 48.62°, 54.62°) and the maximum tan(Δ_max/2) = (1 − n²)/(2n) at sin²θ = 2n²/(1 + n²) (45.94° at 51.34°); project.py Part 9 asserts the scan reproduces all three to 0.01° / 1e-3°.
4. FIXED — notes §3.15 now says "the inhomogeneous wave equations EM-6 (3.37), by retarded potentials".
5. FIXED — notes §7 now says "about 4 seconds" (measured 3.3–4.0 s).
6. FIXED — notes §3.11 fibre sentence now reads "light that meets the core wall beyond the critical angle".
7. FIXED — project.py Part 8 SymPy check now runs for (m ≥ 1, n ≥ 0) and (m ≥ 0, n ≥ 1) and also asserts B⊥ = 0 on all four walls; printed line updated.
8. FIXED — project.py now asserts: n(589 nm) = 1.4571 (1e-4), Cu ω_p = 1.642e16 (1e-3) and λ_p = 115 nm (3e-3, the quote is 3 s.f.), R_Cu = 0.999914 (1e-6), cB/E = 3.27e4 (2e-3), λ_g = 3.953 cm, k = 158.96 rad/m, 4.398 ns/m, B4(c) r_s, phase and 121.5 nm, B5(c) n = 0.8936 and 39.73 µs, C2(b) angles and maximum, C3 doubling times.
9. NOT FIXED — project.py length (now 519 lines) exceeds the 150–400 guideline; the brief's eight parts plus reference solutions require it.
10. NOT FIXED — notes.md word count above guideline (LaTeX-inflated count for a 12-hour module).

Re-run after fixes: `MPLBACKEND=Agg python3 project.py` exits 0 in 3.3 s; every printed number identical to the pre-fix run except the Part 8 message.

STATUS: review complete — 10 findings (0 high, 1 medium, 9 low), 8 fixed, 2 deliberately left (length guidelines); script passes.
