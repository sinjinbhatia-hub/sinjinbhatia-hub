# Review: ST-6 Bayesian computation

Reviewer: independent adversarial pass (2026-11-04). Module: `curriculum/03-statistical-modeling/06-bayesian-computation/`.
Scratch recomputation: `/tmp/claude-0/st6/recompute.py`, `/tmp/claude-0/st6/flat.py`, script run log `/tmp/claude-0/st6/run1.txt`.

## Method

- Ran `MPLBACKEND=Agg timeout 300 python3 project.py`: exit 0, 15.5 s, "All checks passed". Imports: numpy, scipy, matplotlib, sympy, time, pathlib only.
- Recomputed every number in notes §3–§4 and solutions A1–C2 independently (numpy/scipy/sympy): (3.3) integral symbolically; Ex 4.1–4.7; A1–A5, B1–B9, C1–C2. Eight-schools exact posterior recomputed with scipy.quad instead of the project's trapezoid grid: E[τ|y]=3.5977, E[θ_j|y]=(6.21, 4.94, 3.93, 4.76, 3.62, 4.04, 6.30, 4.85), μ̂(10)=3.6625, B_j and E[θ_j|τ=10,y] all match. CAVI iterates (0.14286, 14.2671), (0.37058, 15.0642), (0.39128, 15.1366) → (0.39335, 15.1439) match. B3 leapfrog matrix, det 1, trace 2−w² verified in SymPy. B4(b) posterior of η verified in SymPy. Shadow Hamiltonian (3.15) verified in SymPy. (3.36) ratio (α_n−1)/α_n verified symbolically.
- All stated numbers agree with the script output and with my recomputation. No wrong numbers found.
- Checked every derivation step in §3.1–§3.8 (MH from detailed balance, Gibbs acceptance 1, Albert–Chib conditionals, Theorem 3.3 involution/volume/detailed balance, leapfrog Jacobians, shadow Hamiltonian, n_eff from Var(θ̄), R̂ expectations, shrinkage/(3.23)/(3.24)/p(τ|y), non-centered Jacobian, ELBO, CAVI update, CAVI Normal fixed point, reparametrization gradient). Signs, factors, index ranges correct except item 2 below.
- Cross-references verified in the referenced modules: F5 Theorems 3.22/3.23/3.24/3.28/3.32/3.40/3.41, Definition 3.44, eqs (3.26)/(3.27), §3.12, §3.13 (Gibbs' inequality is Theorem 3.45 there); CM-8 Theorem 3.4 (Liouville), §3.5, §3.11 (symplectic integrators, kick/drift), eqs (3.5), (3.20); ST-5 Theorem 3.3 (NIG), eqs (3.3), (3.10), (3.14) [combined tag 3.14, 3.15], (3.30), §3.8 (Bernstein–von Mises), §3.10 (PPC), Example 4.10 (Jeffreys–Lindley); ST-4 §3.7 (latent-variable probit); ST-3 (ridge, F); LA-8, LA-13 (preconditioning), LA-14 exist. SM-6, SM-7, ML-10, ST-8 have no notes yet (forward links prescribed by the brief; not actionable).
- History (§1, §3): Metropolis et al. 1953 (hard disks, MANIAC), Hastings 1970, Geman & Geman 1984, Duane et al. 1987, Gelfand & Smith 1990, Stan 2012, Rubin 1981, Stein 1956 / James–Stein 1961, Efron–Morris 1975, Albert–Chib 1993, Roberts–Gelman–Gilks 1997 (0.234, 2.38/√d), 0.44 in 1D, Beskos et al. 2013 (0.65), Hoffman–Gelman 2014, Gelman–Rubin 1992, Geyer 1992, Vehtari et al. 2017/2021, Watanabe 2010, Hairer–Lubich–Wanner ch. IX, Peskun — all consistent with my knowledge. One wrong section number in a citation (item 1).
- Template: §0–§10 present with template headings; conventions box; ≥3 worked examples (7); A/B/C exercises (5/9/2 = 16; B slightly above the template's 5–8 target but within STYLE's 12–18 total); §8 specific; §10 matches the brief's sources. Try-it-first/What-just-happened pairs around Theorems 3.1–3.5, Albert–Chib, n_eff, R̂, shrinkage, non-centered, CAVI; two derivations lack a Try-it-first (item 7).
- Brief has no REVIEW FOCUS section (not a programming module).

## Findings

1. **[medium] notes.md §3.4 (scaling paragraph) and §10 (Neal entry).** "Neal (2011, §5.4.4)" / "its scaling (§5.4.4)". Neal's chapter has no §5.4.4; scaling with dimensionality is §4.4 (§5.4 is "Acceptance using windows of states"). Fix: §4.4 in both places. FIXED.

2. **[medium] notes.md §3.6, non-centered derivation step 2.** "The map η_j=(θ_j−μ)/τ has Jacobian τ^{−J}, which cancels the τ^{−J}." Wrong direction: ∂η/∂θ has determinant τ^{−J}; what multiplies the density is |∂θ/∂η| = τ^{J}. (Solution B4(a) does it correctly.) Fix: state that the change of variables multiplies by |det ∂θ/∂η| = τ^J, which cancels the τ^{−J}. FIXED.

3. **[low] notes.md equation numbering.** Tags jump (3.6)→(3.8) and (3.24)→(3.26); (3.7) and (3.25) do not exist, which looks like a deleted equation. Fix without renumbering (exercises/solutions/project cite (3.14), (3.15), (3.18), (3.20), (3.22), (3.34)–(3.36)): tag the Gibbs proposal in Theorem 3.2 as (3.7) and display the exact p(τ|y) in §3.6 step 4 as (3.25); cite (3.25) from the project's `tau_post` docstring. FIXED.

4. **[low] notes.md conventions box vs usage.** Box says "ρ (no subscript) is momentum", but §3.3, §3.4, B2, B3 and project use ρ for the correlation of the bivariate test Gaussian ("the ρ=0.9 Gaussian"). Fix: add the second use to the box. FIXED.

5. **[low] notes.md §3.3 (regression paragraph).** "Project Part 3 matches the exact posterior means to three decimals": script gives 1.0418 vs 1.0412 (they differ in the third decimal). Fix: "to within 10^{-3}". FIXED.

6. **[low] notes.md §3.7 (WAIC paragraph).** "in a regular model with weak priors p_WAIC is close to the parameter count", immediately followed by no pooling p_WAIC = 4.02 for 8 parameters. Not false (asymptotic in observations per parameter) but contradicted by the example on the page: with one observation per θ_j, Var[log p(y_j|θ_j)] = Var(½χ²₁) = ½, so p_WAIC = 4 exactly. Fix: add that clause. FIXED.

7. **[low] notes.md §3.4 (leapfrog) and §3.8 (CAVI for the Normal model).** Numbered derivations without a ▶ Try-it-first prompt (template: before each key derivation; the brief names both as derive-fully items). Fix: add one prompt each. FIXED.

8. **[low] notes.md §3.8 Theorem 3.4 step 3.** "Gibbs' inequality, F5 §3.13" — STYLE asks for the specific object; F5 has it as Theorem 3.45. Fix: cite F5 Theorem 3.45. FIXED.

9. **[low] solutions.md B8.** "Per gradient evaluation the non-centered sampler is roughly 100 times more efficient": 5969/(8000·25) = 0.0299 vs 120/(3200·75) = 0.00050, ratio 60. Fix: "roughly 60 times". FIXED.

10. **[low] solutions.md C2(b).** "log π drops by about 8 between them": log π(4) − log π(0) = 8 − log 2 = 7.31. Fix: "about 7". FIXED.

11. **[low] solutions.md B6(a).** "The chains clearly sit in different places" — "clearly" is banned by STYLE. Fix: drop the word. FIXED.

12. **[low] project.py asserts weaker than the notes' claims.** (a) Example 4.2's μ_n=4.1714, β_n=13.7671 and Example 4.6's μ̂(10)=3.663, school A 11.15, E[τ|y]=3.60, E[θ_j|y] list, Example 4.7's second iterate (0.37058, 15.0642) are printed but not asserted; (b) the scaling assert only requires ratio(d=50) > 2·ratio(d=2) while the notes state 3.8× → 28× and the ordering of RWM and HMC efficiencies. Fix: assert the stated values to the printed precision; require ratio(50) > 4·ratio(2), ratio(10) > ratio(2), and both efficiencies decreasing in d. FIXED.

## Verified, no change

- All 60+ stated numerical values (list above). `d=2` RWM efficiency 0.037–0.042 across acceptance 0.10–0.30 (claim "flat from 0.09 to 0.3" holds). p_B=0.72, LOO/WAIC values, divergence counts, n_eff and R̂ values, B8 table, C2 numbers all reproduce from the seeded script.
- Equation references in exercises/solutions all resolve to the right equations.
- Forward links to SM-6, SM-7, ML-10, ST-8 cannot be verified (modules not written); they are prescribed by the brief and left as is.

## Status

Status: 12 findings (2 medium, 10 low), all FIXED; `project.py` re-run after fixes exits 0 in 14.8 s with output identical to the pre-fix run (log `/tmp/claude-0/st6/run2.txt`).
