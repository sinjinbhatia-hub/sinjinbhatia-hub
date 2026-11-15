# Review: PY-9 pandas and working with data

Module: curriculum/00-python/09-pandas-and-data/ (notes.md, exercises.md, solutions.md, project.py)
Brief: curriculum/research/briefs/PY-9.md. Reviewer: independent adversarial pass, 2026-10-31.
Environment: Python 3.13.16, pandas 3.0.5, NumPy 2.5.3, Matplotlib 3.11.2, SciPy 1.18.1, no pyarrow (same as the module's stated run).
Method: ran project.py headless with -W default and timing; extracted and ran every ```python block of notes.md (25 blocks) and solutions.md (15 blocks) in one session each (isolated copy of the folder for solutions) and diffed stated ```text output against real output; recomputed quoted numbers in scratch scripts under /tmp/claude-0/py9; checked semantics claims (loc inclusivity, alignment, merge, resample closed/label, rolling, CoW) against pandas 3.0.5 behaviour and the user guide; checked cross-references against the referenced modules' files; applied the brief's REVIEW FOCUS list.

Format: one item per finding. Severity: high (wrong output/number/claim, failing or weakened assert, missing required section) / medium (misleading or inconsistent) / low (pedagogy, polish). Status: OPEN / FIXED / NOT FIXED: reason.

## Findings

### 1. [medium] §8 overstates the `apply` slowdown — notes.md §8 bullet "`apply(lambda row: ...)` is the slow path ... (§3.7: hundreds of times slower)"
§3.7's own stated table gives `apply(axis=1)` 172.8 ms against `where` 4.8 ms (36x); re-run 231.6 / 4.1 ms (56x). Only `iterrows` is in the hundreds (196x stated, 280x re-run). Fix: "tens of times slower, and `iterrows` hundreds". STATUS: FIXED: notes §8 now reads "(§3.7: tens of times slower, and `iterrows` hundreds)".

### 2. [medium] "Both require sorted time" is wrong for `resample` — solutions.md B4 last sentence
`resample` groups by bin and is indifferent to row order (verified: an unsorted 3-point series resamples correctly in pandas 3.0.5); only `merge_asof` requires sorted keys ("left keys must be sorted"). Fix wording. STATUS: FIXED: solutions B4 now says `resample` does not care about row order and `merge_asof` needs both tables sorted on the key.

### 3. [medium] `sort_values` called "stable" without qualification — notes.md §3.7 first paragraph
The multi-key sort is stable (verified on 2000 rows with ties), but the default single-key `sort_values` uses quicksort and is not (verified: ties come out scrambled); pandas docs: only `mergesort`/`stable` are stable and `kind` applies to single-column sorts. A reader will generalise the parenthetical to the single-key case. Fix: say which sort is stable and name `kind="stable"`. STATUS: FIXED: notes §3.7 now says the multi-key sort is stable, the default single-key sort (quicksort) is not, and names `kind="stable"`.

### 4. [low] Wrong reason for sorting before comparing merges — solutions.md B2 explanation: "pandas orders rows by key and then by left position"
pandas 3.0.5 preserves the left frame's row order for `inner` and `left` (verified: keys b,a,a,b stay in that order), the right frame's for `right`, and sorts keys only for `outer`. The sort in `canonical` is still right (order is not part of a join's definition); the stated reason is not. Fix wording. STATUS: FIXED: solutions B2 now gives the real orders (left order for inner/left, right order for right, sorted keys for outer) and says the multiset of rows is the definition.

### 5. [low] Prose contradicts its own table — solutions.md B6: "`itertuples` and the plain loop over arrays ... are 20 to 50 times faster than `iterrows`"
Stated table: 2068/51.4 = 40x for `itertuples`, 2068/17.9 = 116x for the array loop (re-run 37x / 100x). Fix: "roughly 40 and 100 times". STATUS: FIXED: solutions B6 now reads "roughly 40 and 100 times faster than `iterrows`", matching its table.

### 6. [low] Complexity claim — exercises.md C3 asks for the rolling mean "in $O(n)$"; the solution (and its docstring) is `searchsorted` + cumsum, $O(n \log n)$, and says so
Fix the exercise to ask for $O(n \log n)$ and mention the two-pointer sweep as the $O(n)$ variant. STATUS: FIXED: exercises C3 asks for O(n log n) and mentions the two-pointer O(n) sweep as optional; solution unchanged (already said O(n log n) search + O(n) arithmetic).

### 7. [low] "Your turn" item 4 asks to time the conversion "on 10^5 rows of raw" — project.py; `raw` has 52,012 rows
Fix: "on `raw` stacked twice (about 10^5 rows)". STATUS: FIXED: project.py "Your turn" item 4 now says "on raw stacked twice (pd.concat([raw, raw]), about 10^5 rows)".

### 8. [low] Template: two ▶ Try it first prompts have no ■ What just happened — notes.md §3.7 (the timing table) and §3.11 Build 5 (hourly mean by `floor` + `groupby` vs `resample`)
The explanatory paragraphs exist; mark them, and for Build 5 say in one sentence why the two agree here and where they differ (Example 4.6). STATUS: FIXED: notes §3.7 paragraph after the timing table and §3.11 paragraph after Build 5 now start with ■ What just happened; the §3.7 one states the measured ratios (about 200x for iterrows, about 35x for apply, from the stated run), the Build 5 one says why floor+groupby and resample agree here and points to Example 4.6 for where they differ.

### 9. [low] `ewm` parameter mismatch — notes.md §3.11 text says `ewm(halflife=...)`, the block uses `ewm(alpha=alpha, adjust=False)`
State the relation alpha = 1 - 2^(-1/halflife) (verified numerically: the two calls agree with adjust=False). STATUS: FIXED: notes §3.11 now says `ewm(alpha=...)` and gives halflife=h as alpha = 1 - 2^(-1/h).

### 10. [low, no action] Timing lines differ from the stated ones by up to about 30% (notes §3.7, §3.8; solutions B1, B6, C2, C3); the conventions box and the solutions header say only ratios transfer, and every ratio-based claim in the prose survives the re-run (except items 1 and 5 above, which are wrong against the stated numbers themselves). NOT FIXED: not an error.

### 11. [no action] Verified correct (recorded so the next reviewer need not redo them): all 25 notes blocks and all 15 solutions blocks reproduce their stated output byte for byte apart from timings and trailing whitespace; `freq="H"` raises ValueError in 3.0.5; `pd.errors.SettingWithCopyWarning` no longer exists and CoW cannot be disabled; `Series.to_numpy()` is read-only under CoW; `pd.errors.Pandas4Warning` exists; `applymap`, `append`, `fillna(method=)` are gone; `MergeError` is a `ValueError`; `ChainedAssignmentError` is a Warning; resample defaults closed/label left for min/h/D; `rolling("1h")` is closed right; `np.unique` on 200k object strings dominates the NumPy groupby (99 ms vs 0.6 ms scatter-add); the df.plot/ax.plot converter warning and the sharex Period-vs-Timestamp mix-up both reproduce with pandas 3.0.5 / Matplotlib 3.11.2; 55-fold category saving; DST day lengths 24/23/24 and 24/25/24; Example 4.3's sd = sqrt(4 + 1.25) = 2.291; Example 4.6's 73.5 = mean(73.0, 73.5, 74.0); C1's 196 = 49 h x 4; C2's 18,047/20,000 ~ 1 - 1/11; project.py asserts test the quoted numbers to 5e-4 (Example 4.3), 3 decimals (4.4), exact lists (4.1, 4.2, 4.6), and the planted slope with attenuation; cross-references PY-2 Example 4.4, PY-3 Example 4.5, PY-6 §3.4/§3.11, ST-3 (3.3) and §3.11, ST-2 bootstrap all name objects that exist; PY-10, ST-7/8/9, ML-1, F4 are roadmap entries. Attributions (Wickham JSS 59(10) 2014, McKinney 3e, VanderPlas 2e ch. 3) are correct. Template: §0–§10 headings, conventions box, 6 worked examples, A/B/C (5/6/3), §8 specific, §10 matches the brief; no banned phrases; no `inplace`, `iterrows` outside demos, chained assignment outside demos, or deprecated API in project.py; generator seeded; headless run 4.1 s, exit 0, no warnings under -W default.

## Re-verification after fixes

project.py: `MPLBACKEND=Agg python3 -W default project.py` exit 0, 3.3 s, no warnings. All 25 notes.md blocks and all 15 solutions.md blocks re-extracted and re-run: stated outputs reproduce apart from timing lines (ratios within the prose's claims) and trailing whitespace. No code in the module changed; only prose and one comment.

Status: 9 findings (3 medium, 6 low) plus 2 no-action records; 9 fixed, 0 open; script passes.
