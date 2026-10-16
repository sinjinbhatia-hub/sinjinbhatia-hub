# CM-5 review (Central forces and the two-body problem)

Reviewer: independent adversarial pass (second reviewer; the first left only the progress stub below with no findings, so the review was redone from scratch). Module folder: `02-physics/a-classical-mechanics/05-central-forces/`. Brief: `research/briefs/CM-5.md`.
Format: one item per finding — file, location, severity (HIGH / MED / LOW), what is wrong, fix. Status marked in Phase 2.

## Phase 1 progress log

- [x] Read brief, STYLE, TEMPLATE, notes.md (all sections), exercises.md, solutions.md, project.py (all 882 lines).
- [x] Ran `MPLBACKEND=Agg python3 project.py`: exit 0 in 5.3 s on the untouched module; all printed numbers agree with the values quoted in notes.md and solutions.md.
- [x] MATH: notes §3.1–3.13 re-derived step by step (reduction, planarity/area law, radial equation and the Lagrangian-substitution sign trap, Binet chain rule, conic solution and focus–directrix, e(E, ℓ), a = −k/2E, Kepler III, vis-viva, Hohmann (3.19), stability ratio 3 + rf'/f and √(n+3), 1/r³ precession and its near-circular check, GR formula, dσ/dΩ Jacobian, hard sphere, Coulomb asymptote argument, Rutherford, σ(θ > θ₀), LRL magnitude). No errors found.
- [x] MATH: worked examples 4.1–4.6 and every number in solutions A1–C4 recomputed independently with mpmath in `/tmp/claude-0/.../scratchpad/recompute.py` (170 values). All agree to the quoted precision except the labelling slip in A4 below.
- [x] MATH: solutions A1–C4 read line by line (B2 sec-φ algebra, B4 Yukawa log-derivative, B5 oscillator, B6 closest approach, B8 dA/dt, C1 resonant term and second-order coefficient, C2 Kepler's equation, C3 1/r² cross-section and impulse approximation incl. I₁ = 2, I₂ = π/2, C4 L4 gradient). No errors found.
- [x] CONSISTENCY: exercises A1–A5, B1–B8, C1–C4 all have complete solutions; numbering matches; no "similarly"/"left to the reader"; conventions box matches usage (F < 0 attractive, E < 0 bound, V = +k/r for scattering, θ vs φ); every equation reference (3.1)–(3.29) in all three files resolves to a tag in notes.md; cross-references checked against CM-3 Ex 4.7, CM-4 Thm 3.1/3.6/3.7/§3.8/Ex 4.4/Ex 4.6/B2(e)/B6(d)/C2, CM-6 (K a = ω² M a), CM-7 line 303 (same sign trap), CM-8 §3.4 (Routhian) and C2, LA-1 column picture, roadmap entries for QM-1, QM-6, EM-8, SM-4.
- [x] CODE: imports numpy/scipy/matplotlib/sympy only; seeded RNG; headless; figures saved. Asserts audited against the notes' claims (findings 1–3, 12).
- [x] TEMPLATE: §0–§10 with template headings; conventions box; 7 Try-it-first and 13 What-just-happened blocks; 6 worked examples; A/B/C exercises (5/8/4); §8 specific; §10 matches the brief (Taylor 8 & 14, Goldstein 3, MIT 8.223 lectures 7–10).
- [x] HISTORY: Tycho/Kepler 1601/1609/1619 and the eight-minutes quote; Halley–Newton August 1684, De motu, Principia 1687, Props. 1 and 11, Book I Sec. XI; Poincaré c. 1890; Schrödinger 1926; Geiger–Marsden 1909 (1 in 8000, platinum), Rutherford 1911, Geiger–Marsden 1913 (5°–150°); Hohmann 1925; bi-elliptic threshold 11.94; Cotes's spirals; Bertrand 1873; Le Verrier 1859, Newcomb ~43″, ~532″ planetary; Einstein Nov 1915; Lagrange 1772 (L4/L5), Euler 1767 (collinear); Routh/Gascheau bound ½(1 − √(23/27)). All consistent with my knowledge except the two softenings below (findings 8, 9).

## Findings

### 1. [MED] project.py Part 4 (n = −3.5 branch, ~lines 300–312) — instability claim never asserted
Theorem 3.9 / notes §7 / Ex 4.5 figure claim that for n = −3.5 the circular orbit is unstable: a 2% slower launch falls in, a 2% faster one escapes. The code integrates both with a terminal `stop` event at r = 0.03 or r = 12 and plots them, but asserts nothing. If the force law or event were wrong the script would still pass. Verified in scratch that both runs terminate by the event (status 1; r_event = 0.03 for ε = −0.02, 12 for ε = +0.02).
Fix: assert `s.status == 1` and that the event radius is < 0.05 (ε < 0) / > 11.9 (ε > 0).
**FIXED.**

### 2. [LOW] project.py Part 6 (~lines 465–466) — Rutherford Monte Carlo asserts far weaker than the quoted result
Notes Ex 4.6 quotes "median relative error 2.7%" (max 12.4% printed). Asserts allow median < 0.12 and max < 0.35. The RNG is seeded and the integration deterministic, so the thresholds can be tightened without flakiness.
Fix: median < 0.05, max < 0.20.
**FIXED.**

### 3. [LOW] project.py Part 4 (~lines 320–323) — quoted measured ratios not asserted
Notes Ex 4.5 quotes measured ω_r/ω_φ = 1.414319 (n = −1) and 0.948627 (n = −2.1); the asserts only check |ratio − √(n+3)| < 3e-3.
Fix: add asserts pinning the quoted values to 1e-5 (they are deterministic).
**FIXED.**

### 4. [LOW] notes.md §3.7 "What just happened" (line ~272) and project.py Part 2 markdown (~line 196) — "random a"
Both say Part 2 uses random a and e. The code uses `a = np.logspace(log10 0.3, log10 5, 10)` (a fixed grid) and random e only.
Fix: "ten values of a from 0.3 to 5 and random e".
**FIXED.**

### 5. [LOW] notes.md §0 (line 9) — wrong CM-4 theorem number
"a time-independent Lagrangian conserves h (CM-4 Theorem 3.7)": in CM-4, Theorem 3.6 is the energy-function conservation law; Theorem 3.7 is h = T + V for fixed constraints. §3.3 step 3 correctly needs both.
Fix: cite "CM-4 Theorems 3.6–3.7".
**FIXED.**

### 6. [LOW] solutions.md A4 — ½q mislabelled
"ΔT = T(√(1+q) − 1) ≈ ½qT = 4332.59 d × 4.772×10⁻⁴": ½q = 4.773×10⁻⁴; 4.772×10⁻⁴ is √(1+q) − 1 (the script prints the latter). Both give 2.07 d.
Fix: write "ΔT = T(√(1+q) − 1) ≈ ½qT = 4332.59 d × 4.77×10⁻⁴ = 2.07 d".
**FIXED.**

### 7. [LOW] solutions.md A1(d) — sign of Earth's barycentric offset depends on an unstated labelling
"By (3.2) Earth sits at −(m_Moon/M) r from the barycenter": the sign is only right if Earth is body 2 and r points Earth→Moon, which is not stated. The distance (4668 km) is correct.
Fix: state the distance and the side ("(m_Moon/M)·384,400 km from the barycenter, on the side away from the Moon").
**FIXED.**

### 8. [LOW] solutions.md A1(b) — date of deuterium's spectroscopic identification
"identified spectroscopically in 1932": Urey, Brickwedde and Murphy observed the shifted Balmer lines in late 1931 and published in 1932; the discovery is usually dated 1931.
Fix: soften to "in 1931–32".
**FIXED.**

### 9. [LOW] notes.md §3.13 — "no further conserved quantities" overstated
Poincaré (and Bruns 1887 for algebraic integrals) showed there are no further *analytic* integrals beyond the classical ten; stated without the qualifier it reads as a stronger theorem than exists.
Fix: "no further analytic conserved quantities".
**FIXED.**

### 10. [LOW] notes.md §3.6, §3.7, §3.11 (Theorems 3.5, 3.7, 3.14) — no Try-it-first prompt before three derivations the brief names as central
Kepler's first law, Kepler's third law and Rutherford's formula are derived without the template's ▶ Try-it-first prompt (7 of 13 derivations have one).
Fix: add a one-sentence prompt before each.
**FIXED.**

### 11. [LOW] Length — notes.md 8,945 words (style guide: 3,500–7,000); project.py 882 lines (style guide: 150–400)
Both exceed the style guide. The brief's topic list (eleven derivations, six worked examples, ten project parts) is what drives the length, and nothing in the excess is wrong or redundant.
**NOT FIXED: cutting correct, brief-mandated material is not a local fix; left for the author/brief owner.**

### 12. [LOW] project.py Part 9 — headline worked-example numbers printed but not asserted
Ex 4.1 (47.4 s, 30.286 km/s, 365.25635 d), Ex 4.3 (112.8°, 87.5 km/s), Ex 4.5 (e = 1.3033e-3), Ex 4.6 (3.84e-5, σ(>90°) = 1626 fm²), Part 5 (0.076163) and C1 (ratio 1.00252) are printed only. The template allows a printed comparison, but pinning the quoted values makes the script fail if a constant or formula changes.
Fix: add asserts for these values at the quoted precision.
**FIXED.**

## Status
STATUS: 12 findings (0 HIGH, 1 MED, 11 LOW); 11 fixed, 1 not fixed (length, deliberately); project.py re-run exits 0 after the fixes.
