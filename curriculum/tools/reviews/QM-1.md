# Review: QM-1 (01-why-quantum)

Independent adversarial review. Module: `curriculum/02-physics/d-quantum-mechanics/01-why-quantum/`.
Method: every derivation in notes §3–§4 and solutions.md re-derived by hand; every number recomputed in an
independent scratch script (hand-entered formulas, scipy constants); `project.py` run headless (exit 0, 2.2 s);
cross-references opened in LA-7, LA-8, LA-12, CM-4, CM-5, CM-8, EM-5, EM-7, EM-8, QM-2..5 and the SM/QM briefs.

Severity: HIGH = wrong math / wrong number / failing code; MEDIUM = assert weaker than the claim, unpinned
seeded number, inconsistency a reader would trip on, dubious attribution; LOW = rounding, wording, pedagogy.

## What was checked and found correct (no action)

- Derivations: mode count and (3.1), Planck average energy and (3.2) with both limits, Bohr levels (3.4)–(3.5)
  and R = me^4/(8 eps0^2 h^3 c), Compton (3.7) steps 1–5, de Broglie standing waves, Theorems 3.3, 3.4, 3.5, 3.6,
  Proposition 3.7 (both matrix products), the Robertson sketch, Larmor precession sign and direction.
- Solutions: all 17 (A1–A5, B1–B9, C1–C3) present, numbered as in exercises.md, complete, no banned phrases.
  Algebra of B2 (reduced mass), B4 (projector identity), B6, B8, B9, C1, C2 (lune area), C3 re-derived.
- Numbers: every quoted value in notes §3–§4 and solutions.md agrees with an independent recomputation except
  the rounding items listed below.
- Cross-references: LA-8 Theorems 3.4, 3.6, 3.9, 3.10, 3.12, 3.13, 3.20, Proposition 3.22, (3.10), (3.15), (3.16),
  §3.1, §3.3, §3.12; LA-7 §3.9; LA-12 (tensor product); CM-4 virial theorem; CM-5 reduced mass; CM-8 (3.24),
  §3.6, §3.9; EM-5 (3.35); EM-7 §3.6 and Poincaré sphere; EM-8 Example 4.4 (1.556e-11 s), (3.25), §3.8–3.9;
  QM-2 §3.5, QM-4 Theorems 3.6/3.7 all exist and say what QM-1 claims. Other modules cite QM-1 §1.4, §3.3, §3.8,
  §3.11, (3.13), Proposition 3.7, Exercise B6: none of these labels were changed.
- Template: sections 0–10 with template headings, conventions box, six worked examples, A/B/C exercises
  (5/9/3), §8 specific and matching the brief, §10 matches the brief's five sources.

## Findings

### Code (project.py)

1. **MEDIUM. project.py Part 3, `check()`.** For p = 0 or p = 1 the tolerance was `4*0.5 = 2` counts, so
   "repeat measurement gives 100%" and "recombined beam gives 0%" would pass with 9999/10000 or 1/10000. The notes
   claim these are exact. Fix: require exact equality when p is 0 or 1.
2. **MEDIUM. project.py, seeded numbers quoted in prose but not pinned.** Notes Example 4.2 (h = 6.61e-34,
   "within 0.3%", W = 2.29 eV), Example 4.6 (4946, 2423, 0.4899, 0.7505), solutions B5(d) (1874, 3752), solutions
   C2 table (0.6670, 0.4990, 0.3337). All match today's output, but only loose asserts (2% on h, 0.08 eV on W,
   4-sigma binomial) protected them. Fix: add exact pins.
3. **MEDIUM. project.py Parts 1–2, asserts weaker than the notes' numbers.** `54.0 < ratio < 55.5` for a quoted
   54.70; `1.5e-11 < t_f < 1.6e-11` for 1.56e-11; `0.12 < z < 0.135 mm` for 0.127 mm; H-beta to 1e-3 only; Planck
   and Rayleigh–Jeans values, Wien peak 579.6 nm, v_max, f0, threshold, E_R, alpha, Balmer list, series limits,
   Compton lam', photon energies, recoil 576 eV, Davisson–Germer 1.3%, baseball wavelength, B7 omega0/period and
   the remaining solutions numbers were printed but never asserted. Fix: pin each quoted value.
4. **MEDIUM. notes §2 / §3.8 step 4 / Example 4.6 vs project.py Part 3.** The notes say Part 3 simulates the
   *SGx* analyzer with recombined beams (P(down) = 0). The code only recombined the 60-degree analyzer. Fix: add the
   SGx recombination (path-amplitude sum and completeness), without consuming RNG so seeded numbers do not move.
5. **MEDIUM. project.py Part 5 / Part 8, (3.16) never tested with a complex overlap.** Every overlap used was
   real (1, 0.5, 0, cos beta; the `exp(0.4j)` in the B8 tags sits on the component orthogonal to t1 and drops out),
   so the conjugate placement in (3.16) and B8(a)'s fringe shift by delta were untested. Fix: add a tag pair with
   overlap r e^{i delta} and assert P = 2A^2[1 + r cos(2kx + delta)] against the explicit joint state.
6. **LOW. project.py Part 5.** "3000 electrons one at a time" were sampled and plotted but nothing checked that
   the samples show fringes with no tag and none with orthogonal tags. Fix: assert the fraction of hits near the
   dark fringes (about 2.5% vs about 25%).
7. **LOW. project.py Part 4 / Part 6.** P(+n'|n) = (1 + n.n')/2 (notes §3.9, Exercise B4) and Theorem 3.3 for a
   non-diagonal H were not asserted. Fix: add both (deterministic).
8. **LOW. project.py Part 2 header comment.** Says "E_n = -E_R/n^2 ... (3.5)"; in the notes E_n is part of (3.4)
   and (3.5) is the Rydberg formula. Fix the labels. Also Part 2 header says "notes 3.2-3.6" (should be 3.3–3.6).
9. **LOW. project.py B7 vs exercises B7.** Exercise gives g = 2.0023; the code used the CODATA g and printed
   1.7609e11, while the hand computation written in the solution (2.0023 x 9.274e-24 / 1.0546e-34) gives 1.7608e11.
   Fix: code uses g = 2.0023 as the exercise states; solution says 1.7608e11 (period 35.68 ps unchanged).
10. **LOW. project.py Part 5 x-label** reads "units of lambda L / d = 20", which is not what the axis is. Fix.

### Notes (notes.md)

11. **MEDIUM (history). §3.3 What-just-happened:** "correct to the precision of 1913 spectroscopy." False as
    stated: spectroscopy was good to parts per million; Bohr's computed R agreed only within the few-percent
    uncertainty of the 1913 values of e, m, h (Bohr says so in the 1913 paper). Fix the sentence.
12. **MEDIUM (history). §1.3:** "a jump between orbits emits one photon" attributed to Bohr 1913. Bohr did not
    accept light quanta in 1913; his postulate was radiation of frequency f = dE/h. Also "Rutherford's 1911
    scattering experiments": the experiments were Geiger and Marsden's (1909–13); 1911 is Rutherford's analysis.
    Fix both wordings.
13. **LOW (history). §1.2:** the threshold frequency is attributed to Lenard 1902 by sentence structure; Lenard
    showed the maximum energy is independent of intensity and rises with frequency. "h to within about half a
    percent" for Millikan is ambiguous (his quoted uncertainty was about 0.5%; his value 6.57e-34 is 0.9% from
    today's). Fix: detach the threshold from Lenard; say "quoted uncertainty".
14. **MEDIUM (consistency). Header vs §0 vs ROADMAP.** Header lists "SM-5 helpful" as a prerequisite, §0 says
    SM-5 comes later, and the ROADMAP makes QM-3 a prerequisite of SM-5 (circular). What §3.1 uses is the
    Boltzmann factor and equipartition, which are SM-4. Fix the header.
15. **LOW (consistency). §3.6:** "...and why wave mechanics came first" reads as contradicting §1.8 (matrix
    mechanics 1925, wave mechanics 1926). The intended comparison (with spin) is made correctly in §3.12. Fix.
16. **LOW (rounding). Example 4.1.** The displayed constants do not reproduce the displayed results: with
    k_B = 1.381e-23 the quotient is x = 5.754, not 5.7551; 8 pi (1.986e-25)/(3.125e-32) = 1.597e8, not 1.598e8;
    k_B T is 6.903e-20, not 6.904e-20. Final values are right. Fix the displayed constants. Also "the temperature
    of the Sun's surface" for 5000 K (it is 5772 K): say "roughly".
17. **LOW (rounding). Example 4.5.** p shown as 5.40e-24, which gives lambda = 0.1227 nm, not the stated (correct)
    0.1226 nm. Fix: show p = 5.403e-24.
18. **LOW (template). §3.11** has a Try-it-first but no What-just-happened. Add one.
19. **LOW (history, not changed).** §1.7 dates Merli–Missiroli–Pozzi to 1974 (experiment; paper 1976) and both
    it and Tonomura 1989 used an electron biprism rather than two material slits. §3.12 "spin was not even
    recognized as a separate two-dimensional degree of freedom until Pauli ... 1927" is a strong phrasing
    (Uhlenbeck–Goudsmit 1925; Pauli's two-component formalism 1927). Both defensible; see status below.

### Exercises / solutions

20. **LOW. solutions A1.** "about 10^-9 of what is needed" (it is 1.8e-9) and "about 5e8 simultaneous photons"
    (it is 5.56e8, which rounds to 6e8). Fix the two rounded figures.
21. **LOW. solutions A5.** Displayed (1.381e-23)(300)/(1.602e-19) gives 0.02586, not the stated (correct)
    0.02585. Fix the displayed constants.
22. **LOW. exercises B8(b) / solutions B8(b).** The tag phase is called delta, the same symbol as the overlap
    phase in part (a), but here it multiplies the component orthogonal to t1 and has no effect on the overlap or
    the fringes. A reader expects the shift of part (a). Fix: rename to gamma and say why it drops out.
23. **LOW. solutions B7(d).** See item 9 (1.7609e11 vs the displayed arithmetic 1.7608e11).

## Status
