# Review: EM-4 — Electric fields in matter

Reviewer pass started 2026-10-20. Module folder: `curriculum/02-physics/b-electromagnetism/04-electric-fields-in-matter/`.
Method: every number in notes §4 and solutions.md recomputed with python3 (see scratch run; all 60-odd values agree to the printed precision); every derivation step in §3.1–§3.15 checked by hand; every EM-1/EM-2/EM-3 cross-reference checked against the actual tag/section in that module's notes.md; `project.py` run (exit 0, 16.8 s) and each assert compared with the sentence in the notes it is supposed to test.

Severity: HIGH = wrong central math/number; MEDIUM = wrong statement, wrong reference, or assert that does not test the notes' claim; LOW = wording/pedagogy/local.

## Findings

1. **MEDIUM — notes.md §0 line 9 and §3.4 step 6 (line 111): wrong cross-reference.** Both cite "EM-2 (3.8)" as the potential of a charge distribution. In EM-2, (3.8) is $V(\mathbf b)-V(\mathbf a) = -\int\mathbf E\cdot d\mathbf l$; the potential of a distribution $V = \frac{1}{4\pi\varepsilon_0}\int\rho/\mathscr r\,d\tau'$ is (3.7). Fix: change both to (3.7).

2. **MEDIUM — notes.md §1 line 17: historical overclaim, internally inconsistent.** "Debye (1912) added molecules with permanent dipole moments and explained why water's $\varepsilon_r$ is 80 while glass's is 5." Debye's theory explains qualitatively why polar liquids have large $\varepsilon_r$ and its $1/T$ dependence; it does not give 80 for water — the module's own §3.14 shows Debye + Clausius–Mossotti returns a negative $\varepsilon_r$ for water and that Onsager/Kirkwood are needed. Fix: soften to "explained qualitatively why polar liquids have $\varepsilon_r$ in the tens while nonpolar glass has about 5 (the quantitative story for water is harder: §3.14)".

3. **MEDIUM (history) — notes.md §1 line 17: "with Lorentz's 1880 correction".** Mossotti/Clausius already obtained the relation from a conducting-spheres model; Lorentz (1880) re-derived it from molecular dipoles via the local field (the route §3.14 takes). Calling it a "correction" is not something I can confirm. Fix: reword to "re-derived by Lorentz (1880) from the field one molecule actually feels, the route §3.14 takes".

4. **MEDIUM — notes.md §3.7 line 185 and §3.14 line 323: stated number contradicted by the code.** Notes say the lattice gives $|\mathbf E|\le0.006\,P/\varepsilon_0$ at interior sites; the code gives $0.00637$ (printed as 0.006 at 3 decimals). Also "at the lattice sites" (line 185) — the probe is only the deep-interior sites $r\le R/2$. Fix: write "about $0.006$" and "at interior lattice sites ($r\le R/2$)"; assert `micro < 0.007`.

5. **MEDIUM — project.py: asserts weaker than what the notes claim.**
   - line 207 `err_in < 2e-3` vs notes "uniform to $2\times10^{-10}$" (actual 2.3e-10) → `< 1e-9`.
   - line 213 `err_out < 2e-3` vs notes "dipole to $2\times10^{-14}$" (actual 1.8e-14) → `< 1e-13`.
   - line 243 `micro < 0.03` vs notes "0.006" (actual 0.0064) → `< 0.007`.
   - line 107 `abs(Qb) < 1e-3*Qabs` vs notes "$Q_b = 10^{-17}$ against 1.8" (actual 7.8e-18 relative) → `< 1e-12*Qabs`.
   - line 129 `rtol=2e-3` on the B2 grid dipole moment vs solutions "the grid gives 0.47872" (5 s.f.) → `rtol=1e-4`.
   - line 589 `abs(slope-3) < 0.03*3` vs notes "slope 3.000" (actual 3.0000000000003) → `< 1e-6`.
   - line 596 `rtol=3e-3` on fixed-Q vs fixed-V force vs notes "agree to 0.02%" (actual 1.75e-4) → `rtol=3e-4`.
   - line 612 `rtol=0.02` on the fringing-field integral vs notes "1.5000 against 1.5000" (actual 2.6e-7) → `rtol=1e-5`; `cum(-1)<0.05, cum(1)>0.97` vs notes "0.6%, 121%, 100.0%" → `cum(-1)<0.01`, `1.15<cum(0)<1.3`, `abs(cum(1)-1)<1e-4`.
   - line 690 `eps_ons ≈ 29.2 ± 0.3` vs solutions "29.1" (actual 29.085) → `29.1 ± 0.05`.
   - lattice count 2109 (notes §3.6, §3.7, §7) and B6 grid fraction 0.156 (solutions B6(d), project comment) are printed but never asserted → add asserts.

6. **LOW — notes.md §3.2 line 75: "EM-3 Exercise A6 did the integral".** EM-3 A6 states the Langevin function and asks the reader to evaluate it; it does not derive it. Fix: "EM-3 Exercise A6 evaluated the Langevin function; SM-4 derives it".

7. **LOW — notes.md §3.14 line 330: "for liquid argon it is a 20% effect".** Clausius–Mossotti exceeds the dilute $\chi_e$ by 17% (0.5087 vs 0.4349; solutions C1(a) quote the inverse, 14.5%). Fix: "a 17% effect".

8. **LOW — notes.md §7 line 429: "matches $C'$ to 1%".** Actual −1.05%. Fix: "to about 1%"; assert `< 0.015`.

9. **LOW — project.py Part 3 markdown (lines 285–288)** says "here eps is sampled at face midpoints", but Example 4.3 uses `harmonic_faces` on node values with the interface midway between nodes; `face_eps` is what Parts 4 and 6 use. Fix the comment.

10. **LOW — solutions.md B4(b) line 57: "just by Gauss for $\mathbf E$"** — style guide bans "just" in this sense. Fix: "directly by Gauss's law for $\mathbf E$".

11. **LOW — project.py length 805 lines** vs STYLE "150–400". Sibling EM modules are 570–690 lines and the brief's project list is long; shortening would delete checks. NOT FIXED (course-wide convention; flag for the author).

## Checked and found correct (no action)
- All numerical values in Examples 4.1–4.5 and solutions A1–C3 (recomputed independently; agree to the printed digits).
- Derivations §3.1 (3.2), §3.4 Theorem 3.3 incl. (3.9) sign, Corollary 3.4, §3.5 counting, §3.6 both routes, Lemma 3.6 and Theorem 3.7, §3.8–§3.10 incl. Theorem 3.10, §3.11 both routes and limits, §3.12 incl. (3.24), §3.13 incl. battery bookkeeping and $(\mathbf P\cdot\nabla)\mathbf E = \tfrac12\varepsilon_0\chi_e\nabla E^2$, §3.14 (3.26)–(3.27), §3.15 (3.28).
- Cross-references: EM-1 §3.5/(P2), (3.18), §3.12; EM-2 (3.2), (3.4)–(3.5), (3.9), (3.10)–(3.12), (3.13)–(3.14), (3.16), (3.18a)–(3.19), (4.1), (4.3), Example 4.1, §3.5, §3.7, §3.9, §3.10; EM-3 §3.2, (3.1), (3.15), (3.24), (3.27), (3.32)–(3.34), Examples 4.4–4.5, (4.3), (4.5). LA-12, LA-13 exist. SM-*, QM-*, EM-5/6/7 are forward references to modules not yet written.
- Exercise/solution numbering A1–A5, B1–B8, C1–C3 match; no "similarly"/"left to the reader"; every solution complete.
- Template: sections 0–10 with template headings; conventions box consistent with usage; Try-it-first/What-just-happened around every key derivation; 5 worked examples (4.5 is the cross-field one); §8 specific; §10 matches the brief; 12 check-yourself questions.
- History: Faraday 1837/Whewell/"specific inductive capacity", Mossotti 1850, Clausius 1879, Debye 1912, Onsager 1936, Kirkwood g — consistent with what I know.
- Code: imports numpy/scipy/matplotlib/sympy only; seeded; headless; 16.8 s; exit 0; figures saved.
