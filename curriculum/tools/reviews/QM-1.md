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

## Second pass (the first reviewer was interrupted after writing the list above)

State found on disk: `project.py` already had the `pin()` helper and the Part 1 pins (Example 4.1 and 4.2
numbers, the seeded fitted h and W); nothing else in the list had been applied, and notes.md, exercises.md and
solutions.md were untouched. The Part 1 pins quote 1.5976e8 and 6.903e-20, which the notes did not yet say
(item 16), so code and prose disagreed until item 16 was applied.

Every item 1–23 was re-verified before acting on it: derivations re-derived (SymPy: Planck mean energy, Bohr
r_n/E_n/v_1, reduced-mass L and kinetic energy, all six Pauli products, Bloch expectation values and the n.sigma
eigenvector, precession phase, B7, B8(a), B9, C3 ratio and its 1 + 3/(2n) expansion; NumPy: brute-force mode count
against k^3 V/(3 pi^2), (3.16) against an explicit joint state with complex tags, the B4 projector identity, the
C2 lune fraction) and every number recomputed in /tmp/claude-0/qm1/check1.py with CODATA constants. All 23 items
are real. None is a wrong final number or a wrong derivation: the module has 0 HIGH findings.

### Additional findings (second pass)

24. **LOW (rounding). solutions C1(a).** The fixed-point iterates from x = 5 are 4.96631, 4.96516, 4.96512; the
    solution lists "4.966, 4.9651", but the second iterate rounds to 4.9652. Fix: list 4.9663, 4.9652, 4.9651.
25. **LOW (consistency). solutions C1(c).** Says project.py checks the Stefan–Boltzmann integral "to 0.2%"; the
    assert is 1e-6 and the actual error is 1.5e-7. Fix the prose.
26. **LOW (style). notes Theorem 3.2.** The statement is the bare formula; the hypotheses (free electron,
    initially at rest, photon scattered through theta) are only in the Try-it-first. Fix: state them.
27. **LOW (math, unstated hypothesis). notes §3.10 What-just-happened.** "The fringe visibility is exactly
    |<t1|t2>|" holds only for equal slit amplitudes; in general it is 2|phi1||phi2||<t1|t2>|/(|phi1|^2+|phi2|^2).
    Fix: add the qualifier.
28. **LOW (template). notes §3.4 and §3.6.** The brief names four key derivations (Bohr levels, Compton shift,
    de Broglie standing waves, the Stern–Gerlach force); the last two have no Try-it-first / What-just-happened.
    Fix: add one-line prompts (no renumbering).
29. **LOW (consistency). exercises.md header vs solutions.md header.** One says the numbers are printed by
    "Part 8", the other by "Parts 2, 6 and 8"; they are printed in Parts 1, 2, 3, 6 and 8. Fix both.
30. **LOW (code). project.py.** R = m e^4/(8 eps0^2 h^3 c) (notes §3.3), the de Broglie standing-wave step
    2 pi r = n h/p => L = n hbar (§3.6), Parseval for the Born probabilities (Postulate 3), and the classical
    Stern–Gerlach band (uniform cos alpha, §3.4) were stated in the notes and never checked. Fix: add asserts.
31. **LOW (history). notes §3.12.** "Spin was not even recognized as a separate two-dimensional degree of
    freedom until ... 1927" contradicts §1.4 (spin proposed 1925). First pass listed it under item 19 as
    "not changed"; it is an internal inconsistency, so fix the wording. Same for §1.7: say the 1974 and 1989
    experiments used an electron biprism in place of slits, and that the Bologna paper appeared in 1976.

Checked again and found correct (no action): section F does not apply (no REVIEW FOCUS in the brief), but the
two code blocks in solutions.md (B8(b), C2(c)) were run and work. Imports are numpy/scipy/matplotlib/sympy plus
stdlib only. No banned phrases. Inbound references from other modules (QM-1 §1.4, §3.3, §3.8, §3.11, (3.13),
Theorem 3.x, Proposition 3.7, Exercise B6, "Part 3") are untouched by the fixes: no label is renumbered.

34. **LOW (length, not an error).** notes.md is about 9,400 words (STYLE: 3,500–7,000) and project.py about 810
    lines (STYLE: 150–400; it was 589 before this review, and the pins and new checks added about 220). exercises.md
    has 9 B problems (TEMPLATE target 5–8; the total, 17, is inside STYLE's 12–18).

Totals after both passes: 34 findings, 0 HIGH, 8 MEDIUM (items 1, 2, 3, 4, 5, 11, 12, 14), 26 LOW
(items 32 and 33 were found while fixing and are listed in the resolution log).

## Resolution log

(Each item is marked here as it is applied.)

- 11 FIXED. §3.3 What-just-happened now says the 1913 agreement was within the several-percent uncertainty of
  e, m, h, and that with today's constants it is right to the reduced-mass correction (5 parts in 10^4).
- 12 FIXED. §1.3: Geiger and Marsden's experiments (from 1909) as analyzed by Rutherford in 1911; Bohr's jump
  "emits radiation of frequency f = dE/h" with a parenthesis that Bohr did not accept light quanta until the
  mid-1920s.
- 13 FIXED. §1.2: the threshold is a separate sentence, not attributed to Lenard; Millikan "quoted uncertainty".
- 14 FIXED. Header says "EM-8 and SM-4 helpful"; §0 says which two SM-4 results §3.1 quotes and that SM-5
  comes after QM-3.
- 15 FIXED. §3.6 no longer says "why wave mechanics came first"; §3.12 says matrix mechanics came a few months
  earlier and why wave mechanics was adopted for calculating.
- 16 FIXED. Example 4.1 now displays hc = 1.9864e-25 J m, k_B T = 0.43087 eV = 6.903e-20 J, x = 1239.84/(500 x
  0.43087) = 5.7551, 8 pi hc/lam^5 = 1.5976e8, RJ with 6.903e-20; each displayed quotient reproduces the stated
  value (rechecked numerically). "Roughly the temperature of the Sun's surface (5772 K)". Matches the Part 1 pins.
- 17 FIXED. Example 4.5 displays p = 5.403e-24 from 9.1094e-31 and 1.6022e-19, and lambda = 6.626e-34/5.403e-24.
- 18 FIXED. §3.11 has a What-just-happened after Proposition 3.7.
- 19/31 FIXED. §1.7: "(Bologna, 1974; published 1976) ... both with an electron biprism playing the role of the
  two slits". §3.12: "Spin, proposed in 1925 (§1.4), did not get its two-component description ... until Pauli
  introduced his matrices in 1927".
- 26 FIXED. Theorem 3.2 states its hypotheses.
- 27 FIXED. §3.10: "for slits of equal strength (|phi1| = |phi2|)" and the role of the overlap's phase.
- 28 FIXED. §3.4 and §3.6 each have a Try-it-first and a What-just-happened; §3.6's two steps are numbered.
  No equation, theorem or section number changed.
- 20 FIXED. solutions A1: "about 2e-9 of what is needed"; "2.3/(4.136e-9) = 5.6e8 photons absorbed at once".
- 21 FIXED. solutions A5 displays k_B T = (1.3806e-23)(300) = 4.142e-21 J = 0.02585 eV.
- 22 FIXED. exercises B8(b) and solutions B8(b): the tag phase is gamma; the exercise asks whether it shifts
  the fringes; the solution explains why it drops out and how to get part (a)'s shift (phase on the whole tag).
- 23 FIXED (solutions side). B7(d) says 1.7608e11 rad/s, matching the displayed arithmetic with g = 2.0023.
- 24 FIXED. solutions C1(a) iterates: 4.9663, 4.9652, 4.9651.
- 25 FIXED. solutions C1(c): "agrees with a T^4 to better than one part in 10^6" (the assert's tolerance).
- 29 FIXED. Both headers now say the numbers are printed and asserted by project.py (Parts 1, 2, 3, 6, 8).
- 32 (new, found by a pin) LOW (rounding). notes Example 4.1: "5.08e-4 J/m^3 per nm" is a double rounding of
  5.07496e5 J/m^4 (5.075e5 -> 5.08e5); to three figures it is 5.07e-4. FIXED in the notes and pinned.
- 1 FIXED. `check()` requires an exact count when the predicted probability is 0 or 1.
- 2 FIXED. Seeded numbers are pinned exactly: Example 4.2 (h = 6.611e-34, W = 2.287, within 0.3%), Example 4.6
  (4946, 2423, 0.4899, 7505 -> 0.7505, sigmas 50 and 35), B5(d) (1874, 3752, sigmas 39 and 48), the C2 Monte
  Carlo column (0.6670, 0.4990, 0.3337). Every new check that needs random numbers uses its own generator, so
  the main stream, and with it every quoted number, is where it was (verified: the run prints the same values
  as before the edits).
- 3 FIXED. Every quoted value in Examples 4.1–4.5, §1.3, §3.3–3.5 and solutions A1–C3 is pinned to half a unit
  in its last quoted digit with `pin()`; the loose range asserts are gone.
- 4 FIXED. `recombined()` builds the analyzer from its two path amplitudes; for SGx and for 60 degrees the code
  asserts completeness, that the output is still |up>, that the two path amplitudes to |down> cancel, that the
  recorded-path probabilities add to 0.5 / 0.375, and that 0 of 10^4 atoms exit z- (exact).
- 5 FIXED. Part 5 tests (3.16) with <t1|t2> = 0.6 e^{0.9 i} against the explicit joint state, asserts that the
  wrong conjugate fails, that P = 2A^2[1 + r cos(2kx + delta)], visibility r, and the fringe shift -delta/2k;
  also the unequal-slit visibility formula of item 27.
- 6 FIXED. The 3000 sampled hits are checked: 2.5% near a dark fringe when amplitudes add, 25% with orthogonal
  tags or classical probabilities (4-sigma binomial against the pattern's own integral).
- 7 FIXED. Part 4 asserts P(+n'|n) = (1 + n.n')/2 = cos^2(gamma/2) and the projector identity on a 20 x 20 grid
  of axis pairs; Part 6 asserts Theorem 3.3 for a non-diagonal H (unitarity, norm, Born sum, the Schrodinger
  equation by finite difference) and that a non-Hermitian generator is not unitary.
- 8 FIXED. Part 2 header: sections 3.3–3.6 and the right equation labels.
- 9/23 FIXED. Code uses g = 2.0023 as the exercise states (checked against CODATA to 1e-5); 1.7608e11 rad/s,
  35.68 ps and 28 GHz are pinned, as is the arithmetic displayed in the solution.
- 10 FIXED. Part 5 x-label.
- 30 FIXED. Asserts added for R = m e^4/(8 eps0^2 h^3 c), the standing-wave condition (SymPy), Parseval and
  (3.12) for the Born probabilities, mutual unbiasedness of the x, y, z bases, the collapse states, the flat
  classical Stern–Gerlach band, |up> saturating the spin uncertainty bound, and the conjugate-slot example of §8.
- 33 (new) LOW. project.py "Your turn" 2 said the He+ 4->3 line "sits at 468.6 nm"; the formula the reader is
  told to use gives 468.65 nm (infinite mass, vacuum), which is what the script prints. 468.6 nm is the line in
  air with the reduced mass. FIXED: comment says both, both pinned.
- notes §7 updated to describe the new checks (Parts 2–6, Part 8 "prints and pins").
- 34 NOT FIXED: length guidelines, not errors. The brief's topic list (six experiments plus the postulates plus
  spin-1/2) and project list are what make the notes and the script long; the script runs in 2 s. Cutting a B
  exercise would remove content the notes point to (B3, B5, B6, B8 are marked essential in §6; B2, B4, B6, B7
  and B8 are cited from §3) or the only practice on Compton kinematics (B1) and on expectation values (B9).

Final check: `MPLBACKEND=Agg timeout 300 python3 project.py` exits 0 in about 2 s, twice in a row with identical
output; 69 `pin()` calls and about 90 asserts. The seeded numbers printed today are the ones the prose quotes
(4946, 2423, 0.4899, 7505, 1874, 3752; h = 6.611e-34, W = 2.287; C2 column 0.6670, 0.4990, 0.3337). Exercise and
solution numbering match (17 = 17); equation tags (3.1)–(3.18) and Theorems 3.1–3.8 are unchanged; math
delimiters balance in all three Markdown files; both code blocks in solutions.md run in the project namespace.

## Status

STATUS: DONE. 34 findings (0 high, 8 medium, 26 low); 33 fixed, 1 not fixed (item 34, length guidelines); project.py exits 0 with all asserts, including the new pins, passing.
