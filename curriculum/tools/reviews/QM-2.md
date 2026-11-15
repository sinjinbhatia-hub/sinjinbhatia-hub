# QM-2 review (independent adversarial pass)

Module: `curriculum/02-physics/d-quantum-mechanics/02-wavefunction-schrodinger/`

**Pass 2 (this file was left by an interrupted reviewer).** The first reviewer wrote items 1–21 and
had started on `project.py` only (Parts 1–9 and the Part 3 markdown); `notes.md`, `exercises.md` and
`solutions.md` were untouched. I re-verified every item below before acting on it (scratch:
`/tmp/claude-0/qm2b/` `nums.py`, `sims.py`, `box.py`, `dtcheck.py`, `quartic.py`), then continued.
Each item carries a status line. Baseline at the start of pass 2: `project.py` exits 0 in about 10 s
(5 s of it is Part 2's SymPy).

What was checked and found clean (no item below): every derivation step in notes §3.1–§3.11
(continuity, momentum operator, Hermiticity, (3.11), Ehrenfest, Theorem 3.7 unit conversion,
Theorem 3.8, group/phase velocity, the Gaussian packet's completing-the-square step, stationary
states, Theorem 3.10, the Madelung split), redone by hand in pass 2; every number in §3.12 and
Examples 4.1–4.6 and in solutions A1–C3 (recomputed with numpy/sympy in `nums.py`: all agree with
the text to the digits quoted, except items 4 and 22); exercise/solution numbering (A1–A5, B1–B8,
C1–C3 on both sides, every solution complete); cross-references F3 (3.8), (3.12), (3.15),
Theorem 3.11 and its centering remark, F3 Exercise B2; F1 §3.8 warning; LA-4 change of basis; LA-8
Theorems 3.12, 3.20; LA-13; CM-8 (3.28), §3.10, Exercise A3; EM-6 continuity; EM-7 (3.29), EM-7 C3;
EM-8 §3.8; QM-1 de Broglie and Davisson–Germer; QM-3 and QM-4 (Robertson is QM-4 Theorem 3.7);
forward references QM-8, QM-9, SM-5, QC-5, QC-9 against their briefs; template sections 0–10,
conventions box, 14 Try-it-first and 12 What-just-happened blocks, six worked examples, 5 A / 8 B /
3 C exercises; §8 and §10 against the brief; no banned words; imports numpy/scipy/sympy/matplotlib
only; history (Bloch's 1976 recollection of Debye, Schrödinger's four 1926 papers, Born's footnote,
Ehrenfest 1927, Kennard 1927, Madelung 1927, Bohm 1952: all consistent with what I know; item 13 is
the one overstatement). QM-3 cites QM-2 (3.4), (3.5), (3.16), (3.28)–(3.30), Theorems 3.4 and 3.10,
Examples 4.2 and 4.5 and Exercise B3 by number, so nothing was renumbered.

## Findings

1. **project.py, Parts 5–9: asserts weaker than the numbers the notes quote.** Severity: medium.
   The notes quote (3.7) "to 10^-13" (assert was `e1 < 1e-4`), quartic drift "1.98" (assert
   `dev > 0.3`), `2 cos t` "to 1.6e-5" (atol 1e-4), energy "5e-6" (1e-4), continuity "5e-6" (1e-4),
   H–J residual "1.5e-8" and "5%" (`< 1e-4`, `> 0.02`), phase velocity "1.005" (`|v-1| < 0.01`),
   step reflection "13.2% vs 13.1%" (only their difference bounded). A later edit could move any of
   these and the script would still pass. Fix: pin each quoted value.
   Status: FIXED. Pass 1 pinned Parts 5–9 (phase velocity against its closed form 1.0048, continuity < 1e-5, 2 cos t to 2e-5, energy to 6e-6, drift 1.981, step 0.1322 and 0.1310, H–J 1e-7 and 0.048 against a closed form); pass 2 re-ran and confirmed each printed value (and that the continuity residual is the O(dt^2) error of the centred time difference: 5.25e-6, 1.31e-6, 3.28e-7 for dt, dt/2, dt/4), tightened the (3.7) check to 2e-12 (see 26), and turned two one-sided bounds into pins: continuity in (4e-6, 7e-6), H–J residual in (5e-9, 5e-8).
2. **project.py, Parts 1, 10, 11: quoted numbers printed but not asserted.** Severity: medium.
   n lambda^3 = 1.7e-7, air and Cu spacings, 579 km, 2.2 d, 1.128 eV, 1.10 um, 0.940 eV, A2's
   150.4 eV, A4's 1.45 and 3.862e-13 m, B6(c)'s 0.56/0.30/0.20, B8(b)'s three errors. Fix: assert.
   Status: FIXED. Part 1 by pass 1. Pass 2 added asserts in Part 10 (sqrt(3) tau, 5.8e5 m/s, 579 km, 1.90e13 s, 1.90e5 s, 2.2 d, 6.02e-20 J, 1.128 eV, 2.728e14 Hz, 3.67e-15 s, 1.10 um, 0.940 eV, 0.4244, 0.1801) and Part 11 (A2: k, v_g, v_p, 2.410e-17 J, 150.4 eV; A4: 3.862e-13 m, m, n, 1.45; A5(b): i hbar/2 on the grid; B4: 14.80, 0.1801, 2.667; B6(b) and (c); B7 from the grid instead of typed-in arithmetic; B8(b) errors and ratios; B8(c) 67.0).
3. **notes §8 and project Part 3 markdown: time-step rule stated with the grid's maxima, which the
   project itself violates.** Severity: medium. Part 4 has (k_max^2/2) dt = 5.2, Part 7 has 11.2,
   Part 8 has 6.6, all above pi, and all accurate; and in Part 7's quartic run the grid's
   V_max dt = 25.9 (pass 2, `dtcheck.py`), where §8 says "if V_max dt/hbar approaches pi ... the
   result is garbage". The rule is about momenta and potential values where the state has weight
   (solutions B8(c) says this correctly). Fix: reword both, and show the failure with a number.
   Status: FIXED. Project markdown by pass 1. Notes §8 rewritten: the rule is about where Psi has weight; it quotes Part 7's own grid maxima (26 and 11, both > pi, energy still good to 6e-5) and a dt = 0.4 run (<H> dt = 3.0, state error 0.79, norm 1 to 1e-12). All five numbers are asserted in Parts 7 and 11. Exercise B8(c) and its solution now use the same run.
4. **exercises B8(c): dx = 0.047 gives p_max = 66.8, the solution says 67.0** (which is for the
   project's dx = 24/512 = 0.046875). Severity: low. Fix: state the grid as `gh`, dx = 24/512.
   Status: FIXED. Exercise B8(c) names the grid (512 points on [-12, 12), dx = 24/512); solution shows 512 pi/24 = 67.0; project asserts dx and 67.0.
5. **notes §3.12: "Localizing an electron below lambda_C requires momentum spread >~ m_e c".**
   Severity: low. With lambda_C = h/m_e c, (3.15) gives only sigma_p > m_e c/4 pi. The statement is
   right for the reduced length hbar/m_e c (sigma_p > m_e c/2), which is what solution A4(a) uses.
   Fix: say hbar/m_e c.
   Status: FIXED. §3.12 now says sigma_x < hbar/m_e c forces sigma_p > m_e c/2 by (3.15), and points to Exercise A4.
6. **notes §3.10 step 4: Sturm–Liouville identification "p = 1, q = V".** Severity: low. In F2's
   form Lf = (1/w)[(pf')' + qf], (3.28) is p = 1, q = -2mV/hbar^2, w = 1, eigenvalue -2mE/hbar^2.
   Fix: write that.
   Status: FIXED. §3.10 step 4 writes (3.28) in F2's form with p = 1, q = -2mV/hbar^2, w = 1, eigenvalue -2mE/hbar^2.
7. **notes §3.11: "This is the derivation CM-8 §3.10 promised."** Severity: low. CM-8 §3.10 carries
   the split out itself (its "capstone"); nothing is promised there. Fix: say CM-8 did it once and
   this section redoes it with A = sqrt(rho) and keeps the term CM-8 dropped.
   Status: FIXED. §3.11 now says CM-8 §3.10 did the split as its capstone and dropped the hbar^2 term; here it is kept.
8. **notes §10: Griffiths "§1.5–1.6 derivation of p-hat".** Severity: low. The derivation is §1.5
   (Momentum); §1.6 is the uncertainty principle. Fix.
   Status: FIXED. §10 says §1.5.
9. **solutions B5(c): "local momentum hbar S_x ∝ 2 hbar C x".** Severity: low. With Psi = R e^{iS/hbar}
   the local momentum is S_x = 2 hbar C x. Fix.
   Status: FIXED. B5(c) now reads S/hbar = C x^2, local momentum S_x = 2 hbar C x, with the signs spelled out.
10. **solutions C1(c): exactness argument incomplete.** Severity: low. "[p^2, x] ∝ p" is neither zero
    nor a number; it drops out because the symmetric splitting cancels the single commutator
    (B8(a)), and only then do the remaining nested commutators reduce to a number. Fix: say so.
    Status: FIXED. C1(c) now says the single commutator is cancelled by the symmetric arrangement (B8(a)), the third-order terms are a number and zero, and everything higher vanishes; adds why the unsymmetric split is not exact.
11. **solutions B8(b): ratio 4.20 left unexplained.** Severity: low. It is exactly what second order
    predicts against a dt = 0.000625 reference: (16^2-1)/(8^2-1) = 4.05, (8^2-1)/(4^2-1) = 4.20. Fix: add.
    Status: FIXED. B8(b) derives 255/63 = 4.05 and 63/15 = 4.20; the project asserts the measured ratios against these to 0.01.
12. **project "Your turn" 2: V = 1e4 wall with the script's usual dt = 0.01 aliases** (V dt/2 = 50 rad;
    pass 1 measured slosh periods 0.33 and 0.41). Severity: low. Pass 2 (`box.py`): with V = 1e4 the
    answer also depends on the grid (0.430 to 0.448 for 1024 to 4096 points, dt from 1e-4 to
    2.5e-5), because dx = 0.003 does not resolve the decay length 1/kappa = 0.007. With V = 1e3 on
    `Grid(-1, 2, 1024)` the period is 0.461 to 0.464 for every dt from 1e-3 to 5e-5 and 512 to 4096
    points, against (1 + 2/kappa)^2 * 4/(3 pi) = 0.463. Fix: use V = 1e3, put the grid, dt and the
    expected number in the prompt, and pin the number with a reference run.
    Status: FIXED. Prompt rewritten for V = 1000 on Grid(-1, 2, 1024), dt = 2e-4, with the reason the old settings fail and the expected 0.46; a reference run now measures 0.462 over 6 periods and asserts |T - 0.46| < 0.005 and agreement with (1 + 2/kappa)^2 * 4/(3 pi) = 0.463 to 1%. Notes §7 updated to match.
13. **notes §1: "He did not get it by guessing. He got it from the analogy".** Severity: low
    (history). Schrödinger's route was messier (a relativistic attempt first; the first paper
    used a variational argument on the Hamilton–Jacobi equation; the optical analogy is the
    second paper's justification). Fix: soften to "his guide was".
    Status: FIXED. §1 now reads "It was not a blind guess. His guide was the analogy ...".
14. **solutions A5(b): "Theorem 3.4(4)"** is step 4 of the derivation, not a part of the theorem.
    Severity: low. Fix wording.
    Status: FIXED. A5(b) cites step 4 of Theorem 3.4's derivation.
15. **project Part 11 comment "energy error ratio"** on the B8 block: it is the state error.
    Severity: low. Fix comment.
    Status: FIXED. Comment now says state error.
16. **notes §3.3: "§3.6 proves it for any Q(p)".** Severity: low. §3.6 shows the operator form
    equals the |Phi|^2 average, given that |Phi|^2 is the momentum density (itself motivated, not
    proved). Fix wording.
    Status: FIXED. §3.3 now says §3.6 shows agreement with the |Phi|^2 average for any power of p.
17. **notes §7 item (8): the 13.2% vs 13.1% gap is not explained.** Severity: low. It is grid
    resolution of the sharp step (dx = 0.12): with 16384 points I get 13.11% (dt = 0.02) and 13.10%
    (dt = 0.005); the plane-wave average is 13.10%. (Pass 2 reproduced 0.13215, 0.13112, 0.13104.)
    Fix: say so.
    Status: FIXED. §7 item (8) explains the gap; Part 8 now also runs the 16384-point grid and asserts 0.1311 within 2e-4 of the plane-wave average (adds 1.8 s).
18. **Theorem 3.8(b),(c) with a nonzero mean (x-hat = i hbar d/dp) is never checked by the
    project.** Severity: low. Fix: add a SymPy check.
    Status: FIXED by pass 1 (Part 2: for a Gaussian with symbolic x_0, p_0, int|Phi|^2 = 1, int p|Phi|^2 = p_0, int Phi*(i hbar d/dp)Phi = x_0). Pass 2 read the code and confirmed it tests (3.17)–(3.19); notes §7 now lists it.
19. **solutions B6(c): "a very narrow packet leaves the classical orbit sooner" is asserted with no
    number.** Severity: low. I ran a = 0.1: the time is 0.33, below a = 0.3's 0.56 (and a = 0.2 gives
    0.48, a = 0.5 gives 0.47), so the claim is true and the optimum sits near a = 0.3. (Pass 2:
    0.325, 0.475, 0.555, 0.465 on the dt = 0.005 grid; a = 0.35 gives 0.565.) Fix: add the a = 0.1
    run to the project and quote it.
    Status: FIXED. a = 0.1 added to Part 11's loop (0.325, asserted below a = 0.3's 0.555); solution B6(c) quotes it against (b)'s 1.83 and gives the times to the step (0.555, 0.300, 0.200).
20. **project: `np.trapezoid` needs NumPy >= 2.0.** Severity: low. Fix: fall back to `np.trapz`.
    Status: FIXED. `trapz = getattr(np, "trapezoid", None) or np.trapz`.
21. **notes §3.4: LA-8 cross-reference.** Severity: low. "Real eigenvalues and orthogonal
    eigenvectors" are LA-8 Theorems 3.9 and 3.10; Theorem 3.12 is the spectral theorem built on them.
    Fix: name all three.
    Status: FIXED. §3.4 names LA-8 Theorems 3.9, 3.10 and 3.12.

### Added in pass 2

22. **notes Example 4.6: "by t = 10 the two disagree by 1.98".** Severity: medium (a quoted
    simulation number that is false as written). 1.98 is the largest gap over 0 <= t <= 10, reached
    at t = 5.57; at t = 10 the gap is 0.43 (`quartic.py`). §3.4 words it correctly ("by up to 1.98
    within t = 10"). Fix: reword the example and pin the time of the maximum.
    Status: FIXED. Example 4.6 now says the gap peaks at 1.98 near t = 5.6 and is 0.43 at t = 10; Part 7 asserts both (5.57 within 0.02, 0.426 within 5e-3).
23. **notes §7: "`project.py` (about 4 s)".** Severity: low. It takes about 10 s since the
    Theorem 3.8 SymPy check went in. Fix: update, and list that check among Part 2's contents.
    Status: FIXED. §7 says about 12 s (5 s SymPy) and lists the Theorem 3.8 check.
24. **solutions header: "Every number below is printed or asserted in Part 11 (or the part named)".**
    Severity: low. B6(b)'s estimates 0.61, 0.26, 0.18 are not computed anywhere in the script, and
    A4(a) and B4(c) do not name the part that checks them. Fix: compute and assert B6(b)'s numbers,
    name the parts.
    Status: FIXED. B6(b)'s three estimates are computed and asserted; A4(a) and B4(c) numbers are asserted in Part 11; header reworded.
25. **project Part 3 markdown: the brief asks for the splitting to be derived in a markdown cell;
    the cell only asserts that the dt^2 terms cancel.** Severity: low. Fix: add the two BCH lines
    (the full argument is solutions B8(a)).
    Status: FIXED. Part 3's markdown cell now carries the two-line BCH derivation and the order of the unsymmetric split.
26. **project Part 7: the quartic run's energy drift (5.9e-5) is bounded only by `< 1e-4` together
    with the harmonic one, and (3.7)'s rounding-level check is `< 1e-11` against a quoted 1e-13.**
    Severity: low. Fix: tighten both, since the reworded §8 quotes the quartic drift.
    Status: FIXED. `e1 < 2e-12`; quartic drift asserted in (5e-5, 7e-5).

27. **solutions B3(b) and B7(c): two loose pointers.** Severity: low. B3(b) said "the number the
    project's step uses" without the formula or the part; B7(c) credited the constant i hbar/2 in
    <xp> to A5(b), which is about real Psi only, when it holds for every state by A5(a).
    Status: FIXED. B3(b) shows R = ((k-q)/(k+q))^2 = 1/9 with q = 1 and names Part 8; B7(c) cites
    A5(a): xp = (xp + px)/2 + i hbar/2.
28. **project.py is about 700 lines against STYLE's 150–400 guideline** (617 before this review).
    Severity: low.
    Status: NOT FIXED: every part is required by the brief's project list or pins a number the notes
    or solutions quote; cutting to 400 lines means removing checks, which is the wrong trade.

Run at the end of pass 2: `MPLBACKEND=Agg python3 project.py` exits 0 in about 12 s, twice, with identical
output apart from the timing line.

Status: review complete. 28 findings (0 high, 4 medium, 24 low); 27 FIXED, 1 NOT FIXED (item 28, project length guideline).

