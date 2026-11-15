# Review: ST-9 Causal inference and experimental design (03-statistical-modeling/09-causal-inference)

One pass (2026-11-15), Phase 1 then Phase 2. No earlier review file existed.
Format: `[severity] file:location — what is wrong — fix`. Severity: high (wrong math/number/assert), medium (incomplete solution, inconsistency, wrong cross-reference, misleading or unconfirmable claim, assert weaker than the prose), low (pedagogy/style).
Scratch recomputation in `/tmp/claude-0/st9/` (`check1.py` worked examples and symbolic formulas, `check2.py` exercises and solutions, `check3.py` Monte Carlo claims with independent seeds).

## Method (filled in as the review proceeds)

- Baseline run `MPLBACKEND=Agg timeout 300 python3 project.py`: exit 0, 1.8 s, "All checks passed". Imports: itertools, time, numpy, matplotlib, pathlib, scipy only.

## Findings

1. **[high] solutions.md C3 table and prose vs project.py Part 13.** The solution quotes, for k=1: mean F 3.22, median bias 0.173, mean bias 0.391; for k=5: mean F 3.26, median bias 0.277, mean bias 0.240. Today's script prints k=1: 3.16, 0.077, −1.163; k=5: 3.30, 0.265, 0.233 (the shared RNG stream moved after the prose was written). The derived ratios (0.240/0.947 = 0.25, 1/3.26 = 0.31, "0.277 vs 0.173", "the sample mean of 1000 draws (0.391)") are stale too. None of the C3 numbers is asserted. Fix: update the table and prose to today's output and pin every C3 number with an assert.
   - **Status:** open

2. **[medium] notes.md §3.5 What-just-happened.** "The project's 'famous' sample shows correlation 0.000 in the population": the script prints −0.006 (and −0.683 among the selected; theory for the cut T+L > 1.5 is −0.678). Fix the number; pin both with asserts (the present assert is `abs(r_all) < 0.02 and r_sel < -0.4`).
   - **Status:** open

3. **[medium] notes.md §3.6 five-DAG table, last row.** "{Z3} alone blocks two paths and opens one": Example 4.6 lists four backdoor paths; Z3 is a non-collider on paths 1, 2 and 3 and the collider on path 4, so {Z3} blocks three and opens one. Fix the count.
   - **Status:** open

4. **[medium] notes.md §3.2 blocking: wrong cross-reference.** "The ANOVA identity (ST-2 Example 4.7's decomposition)": ST-2 Example 4.7 is a 2×2 contingency table. The sum-of-squares decomposition is ST-2 (3.19) in §3.7. Fix the citation.
   - **Status:** open

5. **[medium] notes.md §6 vs exercises.md C2.** §6 says "C1–C3: LATE with defiers, the propensity score theorem, and weak instruments". Exercise C2 is the front-door derivation; no exercise is on the propensity score theorem. Fix §6.
   - **Status:** open
