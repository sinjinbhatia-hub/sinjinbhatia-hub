# QM-3 review (independent adversarial pass)

Module: `curriculum/02-physics/d-quantum-mechanics/03-one-dimensional-problems/`

Baseline: `MPLBACKEND=Agg python3 project.py` exits 0 in about 5 s. Scratch work:
`/tmp/claude-0/qm3/` (`nums.py` recomputes every quoted number from CODATA constants with
independent root-finding; `sym.py` redoes the algebra in SymPy; `blocks.py` runs the code blocks of
solutions.md).

Each item: file, location, severity, what is wrong, the fix, status.

## Findings

1. **notes.md §3.12 "The bands": band edges quoted from a coarse grid, and "4.936 = π²/2" is false.**
   Severity: medium. The notes quote the first band as 1.954 to 4.936 and the second from 9.463;
   `project.py` reads these off a 12000-point energy grid (spacing 0.005). Root-finding
   cos k + 3 sin k / k = ±1 gives 1.9532, 4.9348 (= π²/2), 9.4580, 19.7392. Three of the four quoted
   values are wrong in the last digit(s), and the project only prints them. Fix: root-find the edges
   in `project.py`, assert them (and that the upper edges are n²π²/2), correct the notes to
   1.953, 4.935, 9.458, 19.74.
   Status: open.

2. **notes.md Example 4.8: "the absolute half-lives move by a factor of several if you change the
   1.2 fm" understates the sensitivity, and "within a factor of about 3" is 3.4 for Po-212.**
   Severity: medium. Recomputed: r0 = 1.1 / 1.2 / 1.3 fm gives U-238 half-lives 2.65e11 / 9.87e9 /
   4.28e8 yr and Po-212 1.36e-6 / 8.85e-8 / 6.80e-9 s: a factor of 13 to 27 per 0.1 fm. The absolute
   agreement is therefore partly a matter of the radius choice; the ratio of the two half-lives is
   the robust prediction (10^24.8, 10^24.5, 10^24.3). Model/measured is 2.2 (U) and 1/3.4 (Po).
   Fix: say so, and assert the sensitivity and the two factors in `project.py`.
   Status: open.

3. **project.py: quoted numbers that are printed but not asserted, or asserted loosely.** Severity:
   medium. P(E3) = 0.001370; the six z_n and E_n/V0 of Example 4.5 and κa = 0.83, 12.6%; double-well
   slope −3.924 (assert is "within 5%"); ammonia 9.87e-5 eV; ln T slope −1.99988 (assert is 1e-3);
   single-barrier T = 1.6% at the double-barrier resonance (assert is `< 0.05`); comb edges (item 1);
   thick-barrier 4.2354e-5 and the ratio 9.44e4; STM 1.087e10, 8.79, 1.022; field emission 13.04,
   2.2e-6, 26.08, 4.7e-12; Gamow γ = 44.4, 16.3, 9.87e9 yr, 8.85e-8 s, 10^24.5, 10^23.7 (assert is
   `22 < log10 < 26`); packet 0.3483 / 0.3482 / 0.3470 (assert is a 1e-3 difference); A1 0.3760 eV,
   1.1281 eV, 1099 nm; B2 partial sums 34, 300, 3000; B5 z0 = 8.100, z1 = 1.3974, −9.7024 eV,
   −9.6240 eV; C1 table; C3 56.1 and 10^24.4. All of these agree with my recomputation today (except
   item 1), but nothing stops them drifting. Fix: pin each one.
   Status: open.

4. **notes.md §3.11, What-just-happened after Theorem 3.2: Wronskian sign.** Severity: low. With the
   theorem's own definition W = ψ1'ψ2 − ψ2'ψ1 and ψ1 = e^{ikx}, ψ2 = e^{−ikx}, W = +2ik, not −2ik
   (solutions B4 has +2ik and +2iq). Fix: 2ik.
   Status: open.

5. **notes.md Theorem 3.1 proof: the case split is wrong at z0 = π/2.** Severity: low. "If z0 ≤ π/2,
   g(z0) = tan z0 > 0" fails at z0 = π/2 (tan undefined). Correct split: z0 < π/2 uses g(z0) =
   tan z0 > 0; z0 ≥ π/2 uses g → +∞ as z → π/2⁻.
   Status: open.

6. **notes.md §3.10: "the odd state has a node in the barrier, where the even state has a small
   bump".** Severity: low. Inside the barrier the even state is ∝ cosh κx: a shallow minimum between
   the two humps, small but nonzero. It is not a bump. Fix the wording.
   Status: open.

7. **notes.md §3.5 What-just-happened: "a theme of QM-8's scattering theory".** Severity: low
   (dangling cross-reference). QM-8's brief is approximation methods (perturbation theory,
   variational, WKB, time-dependent PT); QM-9's brief lists scattering theory among the things the
   course does not cover. Fix: drop the pointer and point to the project's "Your turn" item 3, which
   finds bound states as zeros of M22 at k = iκ.
   Status: open.

8. **notes.md §3.6: "in 3D a spherical well needs a minimum depth (QM-6)".** Severity: low (dangling
   cross-reference). QM-6's notes treat only the infinite spherical well (§3.3, Exercise B4). What
   QM-6 does have is u(0) = 0 (QM-6 (3.5)), which makes the ℓ = 0 states the odd states of this
   section, and those need z0 > π/2. Fix: cite that.
   Status: open.

9. **notes.md §3.10: "QC-9's superconducting flux qubits".** Severity: low (dangling
   cross-reference). QC-9's brief (Hamiltonian simulation, VQE, QAOA, adiabatic computing, NISQ
   landscape) has no flux qubits. Fix: keep the physics, attach QC-9 only to what its brief has.
   Status: open.

10. **notes.md §1: "It was the first quantum calculation of a decay rate."** Severity: low
    (history). Dirac's 1927 radiation theory already gave spontaneous-emission rates. What is safe:
    it was the first application of quantum mechanics to the nucleus. Fix: say that.
    Status: open.

11. **notes.md §1: "Mulliken's 1924 isotope shifts".** Severity: low (history). The boron monoxide
    band-spectrum work appeared in 1924 (a note in Nature) and 1925 (Physical Review); the
    half-integer conclusion is in the 1925 paper. Fix: "1924–25" and name BO.
    Status: open.

12. **notes.md §3.8: Ramsauer–Townsend "electrons of about 1 eV".** Severity: low. The cross-section
    minima sit below 1 eV (roughly 0.3 eV in argon, higher in krypton and xenon). Fix: "below about
    1 eV".
    Status: open.

13. **notes.md §3.8: "the form of Geiger and Nuttall's empirical rule".** Severity: low. Their 1911
    rule related the decay constant to the alpha range; ln t½ linear in Z/√E is the modern form.
    Fix: "the modern form of".
    Status: open.

14. **notes.md §10: Townsend "ch. 6–7 ... the double well and ammonia in more depth".** Severity:
    low (attribution). Townsend's ammonia molecule and maser are in his ch. 4 (§4.5), not ch. 6–7.
    Fix: say what ch. 6–7 do and cite §4.5 for ammonia.
    Status: open.

15. **notes.md Example 4.1: reason given for the vanishing boundary terms.** Severity: low. Both
    boundary terms are said to vanish "because u(1−u) is zero at both ends"; the second one vanishes
    because sin(nπu) is zero at both ends. Fix the sentence.
    Status: open.

16. **solutions.md B1(c): wrong reason for discarding words that start with a₊.** Severity: low. "A
    word whose leftmost letter is a₊ ends in a state with n ≥ 1" is false for a word with two of
    each (it would have to pass through n = −1, so it is zero). Correct reason:
    ⟨ψ0|a₊φ⟩ = ⟨a₋ψ0|φ⟩ = 0.
    Status: open.

17. **solutions.md B2(d): only one of the two discontinuities is named.** Severity: low.
    Ψ(x,0) = √(2/a) on (0, a/2) also jumps at x = 0 (the wall forces ψ(0) = 0). Fix: name both.
    Status: open.

18. **solutions.md B3: broken paragraph.** Severity: low. A dangling "Then" before "The two cross
    terms...", and "Check:" runs on after a display equation on the same line. Fix the layout.
    Status: open.

19. **solutions.md C1(d): "10% off at b = 1 and exact to four figures by b = 3", "Both levels move
    symmetrically".** Severity: low. Ratios asymptotic/exact are 0.910, 0.993, 0.9997: 9% low, 0.7%
    low, 0.03% low (0.004958 vs 0.004959 differ in the fourth figure). The shifts from −½ are 0.115
    down and 0.183 up at b = 1 (not symmetric), 0.00245 and 0.00251 at b = 3. Fix both sentences.
    Status: open.

20. **notes.md §3.1 revivals: odd-n states revive at T_rev/8, not first at T_rev/2.** Severity: low.
    The text ("already revive at T_rev/2") is true but suggests T_rev/2 is the earliest. For odd n,
    n² ≡ 1 (mod 8), so at T_rev/8 every phase is e^{−iπ/4}. Checked numerically (overlap 1.000000 at
    T_rev/8, 0.997 at T_rev/16). Fix: add the fact and assert it.
    Status: open.

## Checked and found clean

(filled in at the end)

Status: phase 1 in progress.
