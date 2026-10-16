# Review: EM-2 Electrostatics

Reviewer pass, 2026-10-16. Module: `02-physics/b-electromagnetism/02-electrostatics/`.
Method: every derivation in notes §3–§4 and solutions.md re-derived by hand; every quoted number recomputed
independently in `/tmp/claude-0/em2_check.py` (numpy/scipy) and compared with the text and with `project.py` output;
`MPLBACKEND=Agg python3 project.py` run (exit 0, 12.4 s, 4 figures); all cross-module references opened and checked.

## Verified correct (no action)

- A. MATH: Gauss derivation (3.4)–(3.5); converse symmetry argument; shell/plane/coax fields; curl-free derivation and
  (3.7)/(3.8); Poisson consistency check; boundary conditions (3.10)–(3.12) and signs; (3.13)–(3.16) with the ½ and the
  surface-term fall-off; self-energy integral and electron numbers (1.41 m_ec² at 1 fm, 1.41e3 at 1e-3 fm); conductor
  theorem and force (3.17); capacitor formulas (3.18a–c), (3.19); §3.12 gravity/heat/Jacobi dictionary.
  Ex 4.1 (E, V, Poisson check, 898.8 V/m, 89.88 V, 134.81 V); Ex 4.2 (asinh table, 17.975 V, 41.389 V); Ex 4.3
  (expansion signs, (4.2), (4.3), 0.99 at r = 5d, 2.2e-3 at 20d); Ex 4.4 (3Q²/20πε₀R both ways, 1/6 inside,
  5.3925e-8 J, a_C = 0.720 MeV, 973/606/367 MeV, 0.9968 pair sum); Ex 4.5 (88.54 pF, 8.854 nC, 0.8854 µC/m², 1e5 V/m,
  4.427e-7 J, 4.427e-4 N, 45 mg, 0.0443 Pa, fringing table); Ex 4.6 (σ_wall formula, ratio 27, 7.96 nC/m², 89.88 V).
  Solutions A1–A6, B1–B8, C1–C3: every step and every number re-derived (A1 86.32 V/m, 26.96 V, 1.0243 m; A2 q/24ε₀;
  A4 5.647e5; A5 all; A6 4.866e-7 J, 243.3 V; B1 84.3 min, 7.91 km/s; B2 709 µF, −4.52e5 C, coax ratios; B3 bracket
  (R+r)−|R−r|; B4 4.4e-5/1.76e-4; B5 kq²/2(1/a−1/b); B6 0.1123 N; B8 345.9 V/m at 0.0707 m, disk table incl. 1−¾(R/z)²;
  C1 1.691/1.409 fm, 386.16 fm = 137.036 r_e; C2 g_pq stability bound, 2686/28554, c₁₁ = 4/π² via sinh π/2π, 2510;
  C3 both proofs).
- B. CONSISTENCY: exercises/solutions numbering A1–A6, B1–B8, C1–C3 match; all solutions complete, no "similarly"/
  "left to the reader"/"it can be shown" anywhere; conventions box matches use (n̂ below→above, V(∞)=0 localized,
  𝓇 = r − r′, 𝒱 for volume). Cross-refs resolve: EM-1 Thms 3.7, 3.11, 3.14 (3.32)/(3.33), 3.15, 3.16; CM-1 Thm 3.3;
  LA-13 Thm 3.19 and (3.24) = cos πh; F6 Exercise C3 (Laplace as Euler–Lagrange); F3 impulse response; F2 sine series;
  CM-5 Kepler. EM-6, QM-6, QM-8 are roadmap forward references (ROADMAP.md lines 112, 137, 139).
- C. CODE: exits 0; imports only numpy/scipy/matplotlib/sympy; seeded; headless; asserts test the notes' stated
  tolerances (1e-7 Gauss, 1e-9 V, 1e-10 energies, 1e-10 wall charge, 1e-3 centre field/σ, Jacobi ratio ≈ 4,
  3/5 kQ²/R, 1/6, 0.99, a_C, 973/606 MeV, 27, hemisphere force, C2 sweeps and c₁₁).
- D. TEMPLATE: §0–§10 with template headings; conventions box; ▶/■ around every key derivation; 6 worked examples
  (one from nuclear physics, B1 from gravity); 6 A / 8 B / 3 C = 17 exercises; three code exercises; §8 specific
  with fixes; §10 lists exactly the brief's four sources.
- E. HISTORY: Priestley 1767, Cavendish ~1773 (2%), Maxwell 1879, Coulomb 1785, Newton–Bentley 1693, Poisson 1812/1813,
  Green 1828/Thomson 1845, Gauss 1835/1867 (as in brief), Faraday 1836 foil room and 1843 ice pail, Williams–Faller–Hill
  1971 (~1e-16), Principia I Props 70–71, Poincaré stresses 1906 — all consistent with my knowledge; nothing to soften.

## Findings

1. notes.md §3.4, line 95 (spherical shell) — LOW — "previewing §3.8" points at Work and Energy; the field jump
   σ/ε₀ is the boundary condition of §3.7. Fix: change to §3.7.
2. notes.md §3.5 — LOW (style) — equation tags appear out of order: (3.8) [V = k∫ρ/𝓇] is defined before (3.7)
   [V(b)−V(a) = −∫E·dl]. Fix: swap the two tags and all 17 references in notes.md (exercises.md and solutions.md
   contain no (3.7)/(3.8) references, verified by grep).
3. project.py lines 291 and 396 — LOW — Part 5 header cites "§3.11" (cavities) for the capacitor, should be §3.10;
   Part 6 header cites "§3.10" (capacitors) for the conductor in a field, should be §3.9. Fix the two comments.
4. notes.md §7 — LOW — the build list is numbered 1–7 while the code and every other reference in notes/exercises/
   solutions use the script's Part 0–8 numbering ("project Part 3", "Part 5", "Part 6", "Part 8"). Fix: label the
   list items Part 0, Part 1, Part 2, Parts 3–4, Part 5, Part 6, Parts 7–8 so a reader can find each cell.

No high- or medium-severity findings: no mathematical error, no wrong number, no weakened assert, no missing section,
no incomplete solution, no convention inconsistency, no dubious attribution.
