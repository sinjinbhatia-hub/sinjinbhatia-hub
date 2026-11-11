# QM-2 review (independent adversarial pass)

Module: `curriculum/02-physics/d-quantum-mechanics/02-wavefunction-schrodinger/`
Scratch work: `/tmp/claude-0/qm2/` (check1..check5.py). Baseline: `project.py` exits 0 in 4.3 s.

What was checked and found clean (no item below): every derivation step in notes §3.1–§3.11
(continuity, momentum operator, Hermiticity, (3.11), Ehrenfest, Theorem 3.7 unit conversion,
Theorem 3.8, group/phase velocity, the Gaussian packet's completing-the-square step, stationary
states, Theorem 3.10, the Madelung split); every number in §3.12 and Examples 4.1–4.6 and in
solutions A1–C3 (recomputed with numpy/sympy: all agree with the text to the digits quoted);
exercise/solution numbering (A1–A5, B1–B8, C1–C3 on both sides); cross-references F3 (3.8), (3.12),
(3.15), Theorem 3.11 and its centering remark, F3 Exercise B2; LA-8 Theorems 3.12, 3.20; CM-8 (3.28),
§3.10, Exercise A3; F1 §3.8 warning; EM-7 (3.29), EM-7 C3; QM-3 Theorem 3.2 and (3.23); QM-4
Robertson; template sections 0–10; §8 and §10 against the brief; history (Bloch's 1976 recollection
of Debye, Schrödinger's four 1926 papers, Born's footnote, Ehrenfest 1927, Kennard 1927, Madelung
1927, Bohm 1952: all consistent with what I know).

## Findings

1. **project.py, Parts 5–9: asserts weaker than the numbers the notes quote.** Severity: medium.
   The notes quote (3.7) "to 10^-13" (assert was `e1 < 1e-4`), quartic drift "1.98" (assert
   `dev > 0.3`), `2 cos t` "to 1.6e-5" (atol 1e-4), energy "5e-6" (1e-4), continuity "5e-6" (1e-4),
   H–J residual "1.5e-8" and "5%" (`< 1e-4`, `> 0.02`), phase velocity "1.005" (`|v-1| < 0.01`),
   step reflection "13.2% vs 13.1%" (only their difference bounded). A later edit could move any of
   these and the script would still pass. Fix: pin each quoted value.
2. **project.py, Parts 1, 10, 11: quoted numbers printed but not asserted.** Severity: medium.
   n lambda^3 = 1.7e-7, air and Cu spacings, 579 km, 2.2 d, 1.128 eV, 1.10 um, 0.940 eV, A2's
   150.4 eV, A4's 1.45 and 3.862e-13 m, B6(c)'s 0.56/0.30/0.20, B8(b)'s three errors. Fix: assert.
3. **notes §8 and project Part 3 markdown: time-step rule stated with the grid's k_max, which the
   project itself violates.** Severity: medium. Part 4 has (k_max^2/2) dt = 5.2, Part 7 has 11.2,
   Part 8 has 6.6, all above pi, and all accurate. The rule is about momenta and potential values
   where the state has weight (solutions B8(c) says this correctly). Fix: reword both.
4. **exercises B8(c): dx = 0.047 gives p_max = 66.8, the solution says 67.0** (which is for the
   project's dx = 24/512 = 0.046875). Severity: low. Fix: state the grid as `gh`, dx = 24/512.
5. **notes §3.12: "Localizing an electron below lambda_C requires momentum spread >~ m_e c".**
   Severity: low. With lambda_C = h/m_e c, (3.15) gives only sigma_p > m_e c/4 pi. The statement is
   right for the reduced length hbar/m_e c (sigma_p > m_e c/2), which is what solution A4(a) uses.
   Fix: say hbar/m_e c.
6. **notes §3.10 step 4: Sturm–Liouville identification "p = 1, q = V".** Severity: low. In F2's
   form Lf = (1/w)[(pf')' + qf], (3.28) is p = 1, q = -2mV/hbar^2, w = 1, eigenvalue -2mE/hbar^2.
   Fix: write that.
7. **notes §3.11: "This is the derivation CM-8 §3.10 promised."** Severity: low. CM-8 §3.10 carries
   the split out itself (its "capstone"); nothing is promised there. Fix: say CM-8 did it once and
   this section redoes it with A = sqrt(rho) and keeps the term CM-8 dropped.
8. **notes §10: Griffiths "§1.5–1.6 derivation of p-hat".** Severity: low. The derivation is §1.5
   (Momentum); §1.6 is the uncertainty principle. Fix.
9. **solutions B5(c): "local momentum hbar S_x ∝ 2 hbar C x".** Severity: low. With Psi = R e^{iS/hbar}
   the local momentum is S_x = 2 hbar C x. Fix.
10. **solutions C1(c): exactness argument incomplete.** Severity: low. "[p^2, x] ∝ p" is neither zero
    nor a number; it drops out because the symmetric splitting cancels the single commutator
    (B8(a)), and only then do the remaining nested commutators reduce to a number. Fix: say so.
11. **solutions B8(b): ratio 4.20 left unexplained.** Severity: low. It is exactly what second order
    predicts against a dt = 0.000625 reference: (16^2-1)/(8^2-1) = 4.05, (8^2-1)/(4^2-1) = 4.20. Fix: add.
12. **project "Your turn" 2: V = 1e4 wall with the script's usual dt = 0.01 aliases** (V dt/2 = 50 rad;
    I measured slosh periods 0.33 and 0.41). Severity: low. With dt <= 1e-4 the period is 0.43–0.44
    (ideal 0.4244; effective width L + 2/kappa predicts 0.4365). Fix: put dt and the expected number
    in the prompt.
13. **notes §1: "He did not get it by guessing. He got it from the analogy".** Severity: low
    (history). Schrödinger's route was messier (a relativistic attempt first; the first paper
    used a variational argument on the Hamilton–Jacobi equation; the optical analogy is the
    second paper's justification). Fix: soften to "his guide was".
14. **solutions A5(b): "Theorem 3.4(4)"** is step 4 of the derivation, not a part of the theorem.
    Severity: low. Fix wording.
15. **project Part 11 comment "energy error ratio"** on the B8 block: it is the state error.
    Severity: low. Fix comment.
16. **notes §3.3: "§3.6 proves it for any Q(p)".** Severity: low. §3.6 shows the operator form
    equals the |Phi|^2 average, given that |Phi|^2 is the momentum density (itself motivated, not
    proved). Fix wording.
17. **notes §7 item (8): the 13.2% vs 13.1% gap is not explained.** Severity: low. It is grid
    resolution of the sharp step (dx = 0.12): with 16384 points I get 13.11% (dt = 0.02) and 13.10%
    (dt = 0.005); the plane-wave average is 13.10%. Fix: say so.
18. **Theorem 3.8(b),(c) with a nonzero mean (x-hat = i hbar d/dp) is never checked by the
    project.** Severity: low. Fix: add a SymPy check.
19. **solutions B6(c): "a very narrow packet leaves the classical orbit sooner" is asserted with no
    number.** Severity: low. I ran a = 0.1: the time is 0.33, below a = 0.3's 0.56 (and a = 0.2 gives
    0.48, a = 0.5 gives 0.47), so the claim is true and the optimum sits near a = 0.3. Fix: add the
    a = 0.1 run to the project and quote it.
20. **project: `np.trapezoid` needs NumPy >= 2.0.** Severity: low. Fix: fall back to `np.trapz`.
21. **notes §3.4: LA-8 cross-reference.** Severity: low. "Real eigenvalues and orthogonal
    eigenvectors" are LA-8 Theorems 3.9 and 3.10; Theorem 3.12 is the spectral theorem built on them.
    Fix: name all three.

