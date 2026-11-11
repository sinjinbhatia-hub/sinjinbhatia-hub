# Review: PY-2 Collections, iteration, and the data model

Reviewer pass started 2026-10-25. Module folder: `curriculum/00-python/02-collections-iteration-data-model/`.
Brief: `curriculum/research/briefs/PY-2.md` (has REVIEW FOCUS). Scratch: /tmp/claude-0/py2_*.

Baseline: `MPLBACKEND=Agg python3 project.py` exits 0 in 13.0 s (Python 3.13.16, numpy 2.5.3, matplotlib 3.11.2); every printed
number matches the notes (Part 1 summary, resize points, 10^6-line closed form, 36-point grid, line counts).

Method: (1) every fenced `python`/`pycon` block in notes.md (36) and solutions.md (17) extracted and run by
/tmp/claude-0/py2_harness.py (REPL emulation for bare-expression blocks, doctest for `>>>` blocks, fresh subprocess for scripts),
actual output diffed against the stated output; (2) every closed-form number recomputed independently in /tmp/claude-0/py2_check.py
(sympy/fractions, no module code): sum of squares 333332833333500000, gcd(7919,1000)=1 and 100000 / 94,950,000, the 26/6/4 damping
regimes by exact arithmetic, p*q, p+q, p(2), q(5), p', (pq)', p(q(x)) = 3x^2-4x+2, Knuth's 1/2(1+1/(1-a)) probe counts, 2427.6 and
1012.55, every slice in A1/4.1, A3(b) = [4, 7], hash facts, and the CPython list-growth rule (newsize + newsize>>3 + 6, rounded to 4)
reproducing the quoted resize points [1, 5, 9, 17, 25, 33, 41, 53, 65, 77, 93] — all agree; (3) complexity and hashing claims
checked against CPython facts (list insert(0) O(n), append amortized O(1) with ~12.5% over-allocation, deque O(1) at both ends and
O(n) in the middle, dict/set average O(1), 2/3 load factor, low bits pick the first slot, perturbation mixes in higher bits,
hash(-1) == -2, siphash13 per-process salt, compact-dict holes explaining C2) — all correct and not over-claimed; (4) consistency,
template, cross-references (PY-1 §3.4 isclose and §3.12 Horner; PY-3 `numtools/polynomials.py` "Sparse polynomials from PY-2
section 3.4"; PY-4 notes cite PY-2 §3.4/§3.7; PY-6/PY-9/PY-10 briefs match the §5 claims; CM-6 is small oscillations), and the
REVIEW FOCUS idiom checks (no `open` without `with`, no `for k in d.keys()`, no bare except, no `import *`).

Stated outputs: all 36 notes blocks and all 17 solutions blocks reproduce exactly, except (i) timing/memory numbers, which differ
in value but not in ratio or slope (conventions box says so), and (ii) the B7 solution block, which FAILS (item 1). Exercise and
solution numbering match (A1–A6, B1–B8, C1–C3); every solution is complete and runs. Template: §0–§10 with template headings,
conventions box, ▶/■ pairs around every build-it-yourself sequence, 7 worked examples (4.7 from mechanics), 6 A / 8 B / 3 C
exercises, B7 extends project.py, §8 specific, §10 matches the brief's pairings. History: no attributions beyond "Timsort",
"since Python 3.7", "CPython reserves -1" — all correct.

Format: **ID** | file:location | severity | what is wrong | fix | status

## Findings

**1** | project.py Part 2 (`assert 0.7 < slopes["list: insert(0)+pop(0)"] < 1.3`) and solutions.md B7 (same bound for the queue
variant, in the instruction text and in the `edits` list) | medium | The upper bound 1.3 is not robust: three runs here gave
1.30, 1.22, 1.29 for `insert(0)+pop(0)` (1.30 fails) and 1.37, 1.29, 1.36 for B7's `append+pop(0)`, so the B7 solution's run of the
edited project.py fails its own assert on this machine although the preamble promises timing asserts "pass on any ordinary
machine". Per-n times show the cause: the n=10^3 list (8 KB of arrows) is cache-resident, so the 10^3→10^4 step is ×80 while every
later decade is ×10; the fit over four decades is inflated above 1. The notes say only "grows in proportion to n". | widen the list
front-insertion bounds to (0.7, 1.5) in project.py and in both places in the B7 solution, with a comment; add one paragraph to
notes §3.1 after the timing table saying the fitted exponent often lands at 1.1–1.4 and why, and that "near 1" means linear growth |
FIXED

**2** | solutions.md B8 prose | low | "Measured slopes land a little under the prediction for the list version … and a little above 0
for the dict" contradicts the solution's own printed output (list 1.01, dict −0.05) and the rerun here (0.97, +0.07) | reword:
list near 1; dict near 0 with a sign set by timing noise | FIXED

**3** | notes.md conventions box | trivial | "NumPy appears only in the last cell of `project.py`": NumPy is in Part 9 and three
"Your turn" cells follow it (project.py's own header says Part 9) | "only in Part 9 of `project.py`" | FIXED

**4** | project.py "Your turn" | low | the cell says "Reference solutions follow" but item 2 (lazy `itertools.product` grid) has
none; items 1 and 3 do | add a three-line reference that counts the under-damped points of Part 7's grid lazily and asserts 26 |
FIXED

**5** | project.py Part 8 | low | `assert timings["freq_lists"] > 3 * timings["freq_dict"]` is much weaker than the notes' "far
slower" with 0.1340 s vs 0.0034 s (40×; 40× here too) | tighten to 10× | FIXED

**6** | project.py Part 7 | low | `assert kib_e > 100 * kib_s` against the notes' 87,189 KiB vs 22 KiB (≈4000×; 21 KiB here) |
tighten to 1000× | FIXED

**7** | notes.md §1 (line 26) and `freq_lists` in notes §3.4, solutions B8, project Part 8 | low | `for i in range(len(seen))` with
`seen[i]` where `enumerate` is idiomatic (REVIEW FOCUS item; PY-1 already uses `enumerate`). The §3.4 copy is labelled
"deliberately naive", but the naivety is the linear search, not the indexing | `for i, s in enumerate(seen): … if s == name/w`
in all four copies; line count stays 11 so Part 8's `n_lines == [11, 5, 2]` still holds; §1's comparison count is unchanged |
FIXED

**8** | notes.md §3.4 hashing step 3 | trivial | "(CPython calls this *open addressing*.)" — open addressing is the standard name of
the scheme, not CPython's term | "(The scheme is called *open addressing*; CPython's dicts and sets use it.)" | FIXED

**9** | notes.md Example 4.4 closing paragraph | low | "`groupby` needs sorted input ($O(n\log n)$), and in exchange it streams":
after `sorted(...)` nothing streams; it streams only when the data already arrive sorted | say so | FIXED

**10** | notes.md length | low | 9,147 words including code against the PY addendum's 6,000–8,000 target (ceiling 9,500) | not
changed: under the ceiling, and trimming would mean rewriting sections that are correct | NOT FIXED: within the 9,500 ceiling

**11** | notes.md REPL tracebacks | trivial | tracebacks show `File "<stdin>"` while the conventions box says Python 3.13, whose
interactive REPL prints `<python-input-N>`; a piped/non-interactive 3.13 session still prints `<stdin>` | no change needed |
NOT FIXED: no change needed

(No mathematical, numerical, or complexity errors found; no missing sections; no incomplete solutions; no dubious attributions.)

## Resolution

Items 1–9 fixed: notes.md (6 edits: conventions box, §1 `enumerate`, §3.4 `freq_lists` `enumerate`, open-addressing wording,
new §3.1 paragraph on fitted exponents, Example 4.4 streaming caveat), solutions.md (4 edits: B7 bound 1.5 in the instruction
text and the `edits` list plus one explanatory sentence, B8 prose, B8 `freq_lists` `enumerate`), project.py (5 edits: Part 2 bound
1.5 with comment, Part 7 memory ratio 1000×, Part 8 speed ratio 10×, Part 8 `freq_lists` `enumerate`, "Your turn" lazy-grid
reference asserting 26). Items 10–11 left as is.

Re-run after the fixes: `MPLBACKEND=Agg python3 project.py` exits 0 in 11.9 s, "All checks passed."; printed output identical to
the baseline apart from timings. Harness re-run: notes.md 33/36 blocks match exactly (the 3 others are the timing tables, which
differ in value only); solutions.md: every block runs, the B7 block now exits 0 and prints "All checks passed." (its slopes 1.29 /
0.02 here), and the §1 comparison count (2,448,233) and the `[11, 5, 2]` line counts are unchanged by the `enumerate` edits.

STATUS: review complete — 11 findings (0 high, 1 medium, 10 low/trivial), 9 fixed, 2 not fixed (length within ceiling; `<stdin>`
no change needed); project.py exits 0.
