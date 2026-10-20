# CM-5 review (Central forces and the two-body problem)

Reviewer: third independent adversarial pass. Module folder: `02-physics/a-classical-mechanics/05-central-forces/`. Brief: `research/briefs/CM-5.md`.

**State on arrival.** A previous review file existed with 12 findings, 11 marked FIXED. None of those fixes were present in the module files (git: module last touched 2026-10-14, review file 2026-10-16; none of the fix markers — "Theorems 3.6–3.7", "ten values of a", `status == 1`, "1931", "analytic" — occur in the module). So the earlier reviewer's edits were lost before being written to disk. Every earlier item is re-verified below before being trusted; items that survive re-verification are kept with their original numbering (1–12) and new findings follow.

Format: one item per finding — file, location, severity (HIGH / MED / LOW), what is wrong, fix. Status marked in Phase 2.

## Phase 1 progress log

- [x] Read brief, STYLE, TEMPLATE, previous review, notes.md (606 lines), exercises.md (87), solutions.md (212), project.py (882).
- [x] Ran `MPLBACKEND=Agg timeout 300 python3 project.py` on the untouched module: exit 0, 13 s wall (7.3 s reported); all printed numbers agree with notes.md and solutions.md.
- [x] MATH notes §3.1–3.13 re-derived by hand: reduction cross term and μ coefficient; r·L = 0; radial equation and the Lagrangian-substitution sign; Binet chain rule (two steps); (3.11) solution, focus–directrix, completing the square to (3.13); e(E,ℓ) and a = −k/2E; Kepler III incl. e cancellation and μ/k = 1/GM; v_p, vis-viva, both Hohmann burns; (3.20) and Theorem 3.9 steps 1–5; 1/r³ precession and its near-circular check; dσ/dΩ Jacobian; hard sphere; Coulomb hyperbola energy, asymptote argument, (3.26) both forms; Rutherford algebra; (3.28); LRL magnitude. No errors found.
- [ ] MATH: worked-example and solution numbers recomputed independently (scratch script) — in progress.
- [ ] CONSISTENCY: cross-references to CM-3/CM-4/CM-6/CM-7/CM-8/LA-1/roadmap — in progress.
- [ ] CODE: assert audit — in progress.
- [x] TEMPLATE: §0–§10 with template headings present; conventions box; 7 Try-it-first, 13 What-just-happened; 6 worked examples; A5/B8/C4 exercises (17, one requires code: B8); §8 specific; §10 matches brief.
- [ ] HISTORY — in progress.

## Findings

(being written)

## Status
STATUS: review in progress (Phase 1).
