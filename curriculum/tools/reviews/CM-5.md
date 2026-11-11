# CM-5 review (Central forces and the two-body problem)

Reviewer: third independent adversarial pass. Module folder: `02-physics/a-classical-mechanics/05-central-forces/`. Brief: `research/briefs/CM-5.md`.

**State on arrival.** A previous review file existed with 12 findings, 11 marked FIXED. None of those fixes were present in the module files (git: module last touched 2026-10-14, review file 2026-10-16; none of the fix markers — "Theorems 3.6–3.7", "ten values of a", `status == 1`, "1931", "analytic" — occurred in the module). So the earlier reviewer's edits were lost before being written to disk. Every earlier item was re-verified from scratch below before being trusted; the ones that survive keep their numbering (1–12) and new findings follow (13–15). This pass applied the fixes to disk and re-ran the script (see Phase 2 log).

Format: one item per finding — file, location, severity (HIGH / MED / LOW), what is wrong, fix. Status marked in Phase 2.

## Phase 1 progress log

- [x] Read brief, STYLE, TEMPLATE, previous review, notes.md (606 lines), exercises.md (87), solutions.md (212), project.py (882).
- [x] Ran `MPLBACKEND=Agg timeout 300 python3 project.py` on the untouched module: exit 0, 13 s wall (7.3 s reported); all printed numbers agree with notes.md and solutions.md.
- [x] MATH notes §3.1–3.13 re-derived by hand: reduction cross term and μ coefficient; r·L = 0; radial equation and the Lagrangian-substitution sign; Binet chain rule (two steps); (3.11) solution, focus–directrix, completing the square to (3.13); e(E,ℓ) and a = −k/2E; Kepler III incl. e cancellation and μ/k = 1/GM; v_p, vis-viva, both Hohmann burns; (3.20) and Theorem 3.9 steps 1–5; 1/r³ precession (3.22) and its near-circular check (SymPy: 3 + rf'/f = kr/(kr+λ) ≈ 1 − λ/kr); dσ/dΩ Jacobian; hard sphere; Coulomb hyperbola energy, asymptote argument, (3.26) both forms; Rutherford algebra; (3.28); LRL magnitude. No errors found.
- [x] MATH: every number in Examples 4.1–4.6 and solutions A1–C4 recomputed independently with mpmath (scratch `recompute.py`, ~150 values, no project.py code reused). All agree to the quoted precision except the A4 label slip (finding 6).
- [x] MATH: solutions read line by line (B1 tidal expansion, B2 sec-φ algebra, B4 Yukawa log-derivative, B5 oscillator, B6 closest approach, B7, B8 Ȧ and apsidal |A|, C1 resonant term and second-order coefficient 2.5δ, C2 Kepler's equation and averages, C3 cross-section algebra and impulse approximation with I₁ = 2, I₂ = π/2, C4 gradient and Routh bound). Found the wording problems 13–14 only.
- [x] CONSISTENCY: A1–A5, B1–B8, C1–C4 numbered identically in both files, all solutions complete (no "similarly"/"left to the reader"); conventions box matches usage; every (3.x) tag referenced in the three files exists in notes.md. Cross-references checked against the actual files: CM-3 Example 4.7 (spherical pendulum, same sign trap) ✓; CM-4 Theorem 3.1 ✓, Theorems 3.6/3.7 (finding 5), §3.8 LRL ✓, Example 4.4 (Kepler III from scaling) ✓, Example 4.6 (LRL conserved) ✓, Exercise B2 (two bodies, ten charges) ✓, B6(d) ⟨r⟩ = a(1+e²/2) ✓, C2 negative heat capacity ✓; CM-6 Ka = ω²Ma (3.5) ✓; CM-7 line 303 (same trap) ✓; CM-8 §3.4 Routhian and Exercise C2 (Kepler action-angle) ✓; LA-1 column picture ✓; roadmap entries for QM-1, QM-6, SM-4, EM-8 (dipole radiation; module not yet written) ✓.
- [x] CODE: imports numpy/scipy/matplotlib/sympy only; seeded RNG; headless; figures saved; 7.3 s. Assert audit → findings 1, 2, 3, 12. Scratch probe confirmed the n = −3.5 runs end by the terminal event at r = 0.03 (ε = −0.02) and r = 12.0 (ε = +0.02), so finding 1's assert is safe.
- [x] TEMPLATE: §0–§10 with template headings; conventions box; 7 Try-it-first (finding 10) and 13 What-just-happened; 6 worked examples (brief asks ≥ 6); A5/B8/C4 = 17 exercises, B8 extends project.py; §8 has the brief's five honest items; §10 matches the brief.
- [x] HISTORY: Tycho/Kepler 1601, 1609, 1619, eight-minutes quote, second law found first; Halley Aug 1684, De motu, Principia 1687, Props. 1 and 11, Book I Sec. XI; Poincaré c. 1890 (finding 9); Schrödinger 1926; Geiger–Marsden 1909 (platinum, 1 in 8000), Rutherford 1911, Geiger–Marsden 1913 (5°–150°); Hohmann 1925; bi-elliptic threshold 11.94; Cotes's spirals; Bertrand 1873; Le Verrier 1859, Newcomb ≈ 43″, ≈ 532″ planetary; Einstein Nov 1915; Lagrange 1772 (L4/L5), Euler 1767 (collinear); Routh bound ½(1 − √(23/27)); ʻOumuamua 2017 numbers; G uncertainty 2×10⁻⁵. All consistent with my knowledge except finding 8 (deuterium date).

## Findings

### 1. [MED] project.py Part 4, n = −3.5 branch (lines ~368–375) — instability claim never asserted
Theorem 3.9, §7 and the Example 4.5 figure claim that for n = −3.5 a 2% slower launch falls in and a 2% faster one escapes. The code integrates both with a terminal event at r = 0.03 or r = 12 and plots them, but asserts nothing; a wrong force law or a broken event would still pass. Verified in scratch: both runs terminate (status 1) with r_event = 0.03 and 12.0 respectively.
Fix: assert `s.status == 1` and that the event radius is < 0.05 (ε < 0) / > 11.9 (ε > 0).
**FIXED.** Negative test (same branch integrated with the stable n = −2.5) now fails with "no stop event".

### 2. [LOW] project.py Part 6 (lines ~502–503) — Rutherford Monte Carlo asserts far weaker than the quoted result
Example 4.6 quotes "median relative error 2.7%" (printed max 12.4%); the asserts allow median < 0.12 and max < 0.35. Seeded RNG and deterministic integration, so they can be tightened without flakiness.
Fix: median < 0.05, max < 0.20.
**FIXED.** Negative test (Rutherford prefactor scaled by 1.1) now fails the median assert.

### 3. [LOW] project.py Part 4 (lines ~385–386) — the measured ratios quoted in Example 4.5 are not asserted
Notes quote ω_r/ω_φ = 1.414319 (n = −1) and 0.948627 (n = −2.1) as measured; the code only checks |ratio − √(n+3)| < 3e-3, so the quoted six-figure values could drift silently.
Fix: pin them to 1e-5.
**FIXED.** Negative test (pin moved to 0.948700) fails.

### 4. [LOW] notes.md §3.7 What-just-happened (line 271) and project.py Part 2 markdown (lines 219–220) — "random a"
Both say Part 2 uses random a and random e. The code uses `np.logspace(log10 0.3, log10 5, 10)` for a (a fixed grid) and random e only.
Fix: "a from 0.3 to 5 and random e" / "a on a log-spaced grid from 0.3 to 5 and random e".
**FIXED** in both files.

### 5. [LOW] notes.md §0 (line 9) and §3.3 step 3 (line 144) — wrong CM-4 theorem number
"a time-independent Lagrangian conserves h (CM-4 Theorem 3.7)": in CM-4, Theorem 3.6 is the energy-function conservation law; Theorem 3.7 is h = T + V for fixed constraints. Step 3 of Theorem 3.3 needs both and cites only 3.7. (Confirmed by reading CM-4 notes lines 198 and 234.)
Fix: cite "CM-4 Theorems 3.6–3.7" in both places.
**FIXED.**

### 6. [LOW] solutions.md A4 (line 26) — ½q mislabelled
"≈ ½qT = 4332.59 d × 4.772×10⁻⁴": ½q = 4.773×10⁻⁴; 4.772×10⁻⁴ is √(1+q) − 1 (what the script prints). Both give 2.07 d.
Fix: equate T(√(1+q) − 1) to the 4.772×10⁻⁴ product and give ½q = 4.773×10⁻⁴ as the approximation in parentheses.
**FIXED.**

### 7. [LOW] solutions.md A1(d) (line 18) — barycentric offset sign depends on an unstated labelling
"By (3.2) Earth sits at −(m_Moon/M) r from the barycenter": the sign is only right if Earth is body 2 and r points Earth→Moon, which is not stated. The distance (4668 km) is correct.
Fix: state the distance (m_Moon/M)·384,400 km and the side ("away from the Moon") without the sign.
**FIXED.**

### 8. [LOW] solutions.md A1(b) (line 14) — date of deuterium's spectroscopic identification
"identified spectroscopically in 1932": Urey, Brickwedde and Murphy observed the shifted Balmer lines in late 1931 and published in 1932; the discovery is usually dated 1931.
Fix: "1931–32 (Urey, Brickwedde and Murphy)".
**FIXED.**

### 9. [LOW] notes.md §3.13 (line 473) — Poincaré's result overstated
"there are no further conserved quantities and the general motion is chaotic": Poincaré (and Bruns 1887 for algebraic integrals) showed there are no further *analytic* integrals beyond the classical ten, and not every three-body motion is chaotic.
Fix: "no further analytic conserved quantities and that the motion can be chaotic".
**FIXED.**

### 10. [LOW] notes.md §3.6, §3.7, §3.11 (Theorems 3.5, 3.7, 3.14) — no Try-it-first prompt before three derivations the brief names as central
Kepler's first law, Kepler's third law and Rutherford's formula are derived without the template's ▶ Try-it-first prompt (7 of 13 derivations have one).
Fix: add a one-sentence prompt before each.
**FIXED** (now 10 prompts).

### 11. [LOW] Length — notes.md 8,945 words (style guide 3,500–7,000); project.py 882 lines, now 904 (style guide 150–400)
Both exceed the style guide. The brief's topic list (eleven derivations, six worked examples, ten project parts) drives the length; nothing in the excess is wrong or redundant.
**NOT FIXED: cutting correct, brief-mandated material is not a local fix; left for the author/brief owner.**

### 12. [LOW] project.py Parts 5, 8, 9 and C1 — headline worked-example numbers printed but not asserted
Example 4.1 (0.999996997, 449 km, 30.286/29.291 km/s, 365.25635 d, 47.4 s), 4.3 (e = 1.200, 112.8°, 87.5 km/s), 4.4 (Δv 2399.5 + 1457.2 = 3856.7 m/s, 5.290 h), 4.5 (e = 1.3033×10⁻³), 4.6 (227.5 MeV fm, σ(>90°) = 1626 fm², 3.84×10⁻⁵, 518 fm²/sr, 1.03°), §3.10 (0.076163 rad, 42.98″) and C1 (ratio 1.00252) are printed only. The template allows a printed comparison, but pinning them makes the script fail if a constant or formula changes under the quoted text.
Fix: add asserts at the quoted precision.
**FIXED** (asserts added at each print site, tolerances at the quoted precision).

### 13. [LOW] solutions.md B1(c) (line 54) — "external energy" has the wrong sign in words
"The external energy is −[m₁Φ(r₁) + m₂Φ(r₂)]": with Φ the potential per unit mass (so that part (a)'s Φ = −g·x gives +m g·r in L), the potential energy is +[m₁Φ(r₁) + m₂Φ(r₂)]; the minus sign belongs to its appearance in L. The rest of the solution (the −½μ rᵀHr term "in L") is consistent with the latter reading.
Fix: "The external potential energy is m₁Φ(r₁) + m₂Φ(r₂) (it enters L with a minus sign)".
**FIXED.**

### 14. [LOW] solutions.md C1(d) (line 152) — residual attributed entirely to eccentricity
"the remaining difference is the e = 0.2 correction": the near-circular second-order prediction is 1.00250 and the measured ratio 1.00252; integrating the same equation with e → 0 (scratch) gives 1.00251, so roughly half the residual is O(δ³), not eccentricity.
Fix: "the remaining 2×10⁻⁵ is shared between the O(δ³) terms and the e = 0.2 correction".
**FIXED.**

### 15. [LOW] notes.md §2 (line 60) — guess attributed to an argument that cannot give it
"their energy and (by CM-4 Example 4.4's scaling) their period depend only on a": mechanical similarity gives T² ∝ a³ among geometrically similar orbits (fixed e); it cannot show e-independence, and §3.7's What-just-happened says so explicitly ("which CM-4's scaling argument could not show"). The two sentences contradict each other.
Fix: "their energy and, as CM-4 Example 4.4's scaling suggests, their period depend only on a".
**FIXED.**

## Phase 2 log

- Fixes applied with an exact-match replacement script (8 edits notes.md, 5 solutions.md, 12 project.py; each old string matched exactly once).
- `MPLBACKEND=Agg timeout 300 python3 project.py`: exit 0, 4.5 s. Printed output byte-identical to the pre-fix run (only timing lines differ), so no quoted number changed; the quoted numbers are now pinned by asserts.
- Negative tests on scratch copies: breaking the n = −3.5 force law, the pinned Example 4.5 ratio, or the Rutherford prefactor each makes the script exit 1 at the intended assert.
- Only the three module files and this review file were modified by this pass.

## Status
STATUS: 15 findings (0 HIGH, 1 MED, 14 LOW); 14 fixed, 1 not fixed (length, deliberately); project.py exits 0 after the fixes with unchanged output.
