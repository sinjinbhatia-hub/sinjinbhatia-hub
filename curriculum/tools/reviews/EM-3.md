# Review: EM-3 Boundary-value problems

Reviewer: independent adversarial pass. Module folder: curriculum/02-physics/b-electromagnetism/03-boundary-value-problems/.
Script status at start: `MPLBACKEND=Agg python3 project.py` exits 0 in 6.7 s.

Format: one item per finding. Severity: HIGH (wrong math/number/assert), MEDIUM (misleading claim, dubious attribution, template gap), LOW (pedagogy/cosmetic).

## Findings

### A. Math (notes §3–§4, solutions)
Every derivation step in §3.1–§3.10 re-checked by hand and the following re-verified with SymPy/NumPy in /tmp/claude-0/em3_check{1,2}.py: (3.3) and the Q(<S) integral; (3.6)–(3.7) and the all-θ cancellation; B2(b) σ(θ) and its integral; B2(c) closed-form neutral-sphere force; (3.14) ratio asymptotics 1 − 2/k, (k+2)a_{k+2}/(k a_k) = 1 + O(k⁻²); Theorem 3.5 step 4 integral; the binomial expansion of the generating function; (4.1) closed form; Example 4.6(b) V₂ = 2kqd²P₂/r³; C3 ⟨g²⟩ = 2/3; Jacobi-vs-SOR factor (≈51 at h = a/40). Every numerical value in Examples 4.1–4.6 and solutions A2, A4–A6, B2(c)(d), B3(b), B5(c), B6(b), B8, C3(c) recomputed independently; all agree with the text and with project.py's output to the printed precision. No mathematical or numerical errors found.

### B. Consistency
- Exercise/solution numbering matches (A1–A6, B1–B8, C1–C3); every solution is complete; no "similarly"/"left to the reader".
- Conventions box matches usage (σ = −ε₀∂V/∂n with n̂ out of the metal, P_ℓ(1) = 1, θ from the symmetry axis, Griffiths moments, Condon–Shortley Y_ℓ^m).
- Cross-references verified against the referenced notes: EM-1 Theorems 3.8, 3.11, (3.22), (3.23), §3.10; EM-2 (3.9)–(3.12), (3.14), (3.16)–(3.17), Theorems 3.4(i)(iv)(v), 3.5, §3.7, §3.8, Examples 4.2, 4.3, 4.6, (4.1), (4.3), project Part 6 (0.989); F2 §3.5, §3.12, Theorems 3.2, 3.5, 3.9, 3.10, 5.2, §5 Sturm–Liouville; LA-8 Theorems 3.10, 3.19, Example 4.1; LA-13 §3.9, Theorem 3.19, (3.24); CM-7 inertia tensor; EM-4 uses EM-3 Example 4.5, (3.15), (3.32)–(3.34) consistently. All exist. Three do not resolve (items 1–3 below).

1. [MEDIUM] FIXED — notes.md §3.9 step 2 (line ~267) and project.py `V_pts` docstring (line ~232): cite "EM-2 (3.8)" for V = k∫ρ/𝓇 dτ'. In EM-2, (3.8) is V(b) − V(a) = −∫E·dl; the potential of a distribution is EM-2 (3.7). FIX: change both to EM-2 (3.7).
2. [LOW] FIXED — notes.md §3.9 quadrupole paragraph (line ~279): "In LA-12's language (§3.18), Q_ij is a symmetric traceless tensor of type (0,2); under rotations its five components mix only among themselves (an irreducible piece)". LA-12 §3.18 works the stress and inertia tensors and the (0,2) law (3.19); it says nothing about tracelessness or irreducible pieces. FIX: cite LA-12 (3.19)/§3.18 only for the (0,2) tensor law and the inertia-tensor analogue; state the irreducible-decomposition fact as a fact (proved via QM-5's angular-momentum algebra), not as LA-12 content.
3. [LOW] FIXED — notes.md Example 4.3 (line ~337) "red-black SOR, LA-13 §3.9" and §5 LA-13 bullet (line ~381); project.py Part 3 markdown "(LA-13)". LA-13 §3.9 covers Jacobi and Gauss–Seidel (Theorem 3.19, (3.24)) but not SOR; red-black SOR was introduced in EM-2 project Part 6. FIX: attribute SOR to EM-2 Part 6 and LA-13 §3.9 to the convergence theory.
4. [LOW] FIXED — notes.md §5 CM-5 bullet (line ~385): names CM-5 but not the object; CM-5 §3.10 is "Precession under small perturbations" (inverse-cube), the relevant place. FIX: cite CM-5 §3.10 explicitly.
5. [LOW, cosmetic, NOT FIXED] Equation tags in §3.6 run (3.15), (3.17), (3.18), (3.19), (3.16), (3.20): Theorem 3.5 carries (3.16) but appears after (3.19). References to (3.16) are consistent everywhere (27 cross-file uses of 3.16–3.19); renumbering is not a local change. Left as is.

### C. Code
- `MPLBACKEND=Agg timeout 300 python3 project.py` exits 0 in 6.7 s; imports are numpy/scipy/matplotlib/sympy/stdlib only; seeded RNG; headless; figures saved. Asserts checked against the notes' claims one by one.
6. [LOW] FIXED — project.py line ~291 `assert np.abs(V_pts(qs_sp, pos_sp, R*u)).max() < 1e-12 * K * q / R` tests |V| < 9×10⁻¹¹ V, while notes §3.3 step 4 and solutions B2(a) claim "|V| < 10⁻¹³ V" (observed 5.9×10⁻¹⁴). Assert weaker than the claim. FIX: assert |V| < 10⁻¹² V absolute and state "below 10⁻¹² V" in notes and solutions (17× margin over the observed value).
7. [LOW] FIXED — project.py Part 7 C3 loop: `assert abs(exact / lead - 1) < 0.1` for r ≥ 1 nm, while solutions C3(c) claims the leading term is "exact to a part in 10³ or better" there (observed 6×10⁻⁴ and 2×10⁻⁵). FIX: tighten to 2×10⁻³.

### D. Template
- Sections 0–10 present with the template headings; conventions box; 10 Try-it-first and 8 What-just-happened; six worked examples (brief asks ≥6); 6 A / 8 B / 3 C exercises; B8 extends project.py; §8 specific; §10 matches the brief's four sources.
8. [LOW] FIXED — notes.md §3.9: the multipole derivation (a key derivation with a ▶ prompt) has no ■ What just happened after its ∎. FIX: add one.
9. [LOW, NOT FIXED] §9 has 14 check-yourself questions (template: 8–12). Harmless; left.
10. [LOW, NOT FIXED] STYLE asks for one worked example "from another field"; all six are electrostatics (the brief prescribes them). Cross-field content lives in Example 4.3's LA-13 comparison, §5, and C3 (SM-6). Not a local fix; left.
11. [LOW] FIXED — notes.md §10 Griffiths line: lists 4 + 8 + 7 = 19 problems and says "those are the fifteen". FIX: wording.

### E. History / attributions
- Confirmed from own knowledge: Laplace 1780s; Legendre 1782 spheroid; Poisson 1812 two spheres; Green 1828 essay; Thomson 1845 (Green's essay and images); Fourier 1807; Earnshaw 1842; Paul/Penning traps; classical collapse ~10⁻¹¹ s; Earth J₂ = 1.08×10⁻³; water 1.85 D, HCl 1.08 D; Griffiths Table 3.3 values; Stein & Shakarchi ch. 3 for completeness.
12. [MEDIUM] FIXED — notes.md §8 (line ~416) "Arfken & Weber and many chemistry texts do not [include the Condon–Shortley phase]" and §10 (line ~441) "Its Y_ℓ^m convention omits the Condon–Shortley phase; be warned." Arfken & Weber define Y_ℓ^m with an explicit (−1)^m (6th ed. eq. 12.153; their Table 12.3 has Y_1^1 = −√(3/8π) sin θ e^{iφ}), i.e. they DO include it. Cannot confirm the claim; believe it is wrong. FIX: remove the Arfken attribution; keep a generic "check the convention" warning.
13. [LOW] FIXED — notes.md §3.9 (line ~284) "Griffiths §3.4.4 and Problem 3.48": problem numbers differ between Griffiths editions (3.42 in the 3rd, 3.48 in the 4th). FIX: soften to the section plus a description of the problem.
14. [LOW] FIXED — notes.md §1 (line ~17) "Riemann later named the first form Dirichlet's problem": Riemann (1857) coined "Dirichlet's principle" for the variational form; the boundary-value problem acquired Dirichlet's name through the same lectures. FIX: soften.
15. [LOW] FIXED — notes.md §3.6 (line ~193) "proof in Arfken & Weber §12.1": the Rodrigues-formula proof is not in §12.1 (generating function) in the 6th edition. FIX: cite ch. 12 without a section number.


## Resolution
- FIXED (12): items 1, 2, 3, 4, 6, 7, 8, 11, 12, 13, 14, 15. Edits: notes.md (13 replacements), solutions.md (1), project.py (4). Each replacement was applied by exact-once string match; `git diff` reviewed.
- NOT FIXED (3, all LOW, by design): item 5 (equation-tag order (3.16) after (3.19): cosmetic, 27 consistent cross-file references, renumbering is not a local change); item 9 (14 check-yourself questions vs template 8–12: harmless); item 10 (no worked example from another field: the brief prescribes six electrostatics examples; adding one is not a local fix).
- Re-run after fixes: `MPLBACKEND=Agg timeout 300 python3 project.py` exits 0 in 6.3 s, all asserts pass with the two tightened thresholds (|V| on the image sphere < 1e-12 V; C3 leading-order ratio within 2e-3 for r >= 1 nm).

STATUS: review complete; 15 findings (0 high, 2 medium, 13 low), 12 fixed, 3 low cosmetic items deliberately left; math and all numbers verified correct; script passes.
