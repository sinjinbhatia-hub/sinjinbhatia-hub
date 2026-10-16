# EM-1 review (Vector calculus toolkit)

Reviewer: independent adversarial pass. Module folder: `02-physics/b-electromagnetism/01-vector-calculus-toolkit/`. Brief: `research/briefs/EM-1.md`.
Format: one item per finding — file, location, severity (HIGH / MED / LOW), what is wrong, fix. Status marked in Phase 2.

## Phase 1 — in progress

(If this section still says "in progress", the reviewer was interrupted; re-verify each item before trusting it.)

### Checked so far
- Read brief, STYLE, TEMPLATE, notes.md (all sections), exercises.md, solutions.md, project.py in full.
- Ran `MPLBACKEND=Agg timeout 300 python3 project.py` on the untouched module: exit 0 in 9.3 s; imports are numpy/scipy/matplotlib/sympy/itertools/pathlib only.
- Hand-checked every step of notes §3.1–§3.13 (Thm 3.1 proof, BAC-CAB indices, §3.2 transformation laws, Thm 3.3/3.4 face and edge pairings and signs, (P4) and v×(∇×A) index derivations, Thm 3.5 lemma, Thm 3.6, sphere/disk da cross products, Thm 3.8–3.11, spherical scale factors, r̂×θ̂=φ̂, Thm 3.12–3.13 face areas, (3.29), the §3.10 table incl. spherical curl, Thm 3.14 steps 1–6 incl. the sign of da on r=ε, Helmholtz steps 1–4 incl. the surface-term order, (3.35) and its checks, Thm 3.16–3.17, gauge paragraph).
- Hand-checked Examples 4.1–4.7 (all derivatives, face table, curl of (z−y,x,xy), hemisphere integral, erf potential algebra incl. 4π(2π)^{-3/2}=√(2/π), −4π(1+2σ²)^{-3/2} values, U for Ex 4.6, Fourier coefficients in Ex 4.7).
- Hand-checked solutions A1–C4 every line (A2 cross products, A3, A4, A5(c) Gaussian mean 1+3σ², A6 helix integral and the ẑ=∇×A trick, B2(c) weak-curl argument, B3, B4 vector area, B5, B6(d), B7, B8 discrete product-rule algebra, C1(b), C2 bounds, C3 steps 3–4 signs, C4 index gymnastics and the reflection check).
- Cross-reference and numeric recomputation: below.

### Findings

(appended as found)
