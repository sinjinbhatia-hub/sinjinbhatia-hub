# Review: PY-3 — Functions as values, closures and decorators, recursion, modules and the standard library

Reviewer pass, adversarial. Module folder `curriculum/00-python/03-functions-modules-stdlib/`.
Method: ran `project.py` (exit 0, 4.4 s, CPython 3.13.16); extracted every fenced `python`/`pycon`/
`console` block from notes.md (52 blocks: 27 scripts, 10 REPL transcripts via doctest, 6 `# file:`
files, 5 shell sessions = 17 commands) and solutions.md (28 blocks: 14 scripts, 1 REPL, 7 files,
3 shell sessions = 9 commands) into scratch dirs and ran them with the conventions-box rules
(scratch dir, or the module folder when `numtools` is mentioned); diffed each against its stated
output; recomputed every quoted number (C(n) = 2F_{n+1} - 1 for n = 20, 25, 30; phi^5; cache_info
hits/misses for fib(30), fib(100), fib(2000); T(n) = n! sum 1/k! for n = 2..7 and n = 20;
binary-search step counts; statistics/datetime/timestamp values; B2 LRU hit counts; B4 depth 997;
B6 mean 8/3; C2 hit rates) in scratch scripts; checked LA-6 (Thm 3.12, section 3.5 cost formula,
Example 4.1 matrix and det 7), PY-1 (eqs 3.9, 3.10, Example 4.6, LEGB) and PY-2 (section 3.4
hashing/sparse polynomial, section 3.10 pathlib/csv/json) cross-references against those notes;
checked template sections, STYLE word bans, outdated idioms (`%` formatting, `os.path`,
class-based decorators: none), the brief's REVIEW FOCUS items.

Result of the output check: every block reproduces its stated output exactly except (i) three
tracebacks whose *path* is shown relative where CPython 3.9+ prints the absolute path (M1),
(ii) `ls numtools` without `__pycache__` (L1), and (iii) timing numbers (Example 4.3, 4.7, section
3.5, B7, C3), which the conventions box declares machine-dependent; all stated orderings and ratios
reproduce (naive fib(25)/fib(20) = 11.3 here vs stated 10-12). The `--help` text, every
`cache_info()`, every error message and every exit status match. The induction for (3.1), the
2n - 1 call count, the binary-search invariant and step bound, the merge-sort bound and the LA-6
cost recurrence all check. No HIGH-severity findings.

## Findings

### M1 — notes.md section 3.8 (lines 785, 829) and solutions.md B5 (line 441) · MEDIUM · stated output
`File "proj/scripts/run.py"`, `File "numtools/cli.py"` and `File "geom/stats.py"` are not what
Python 3.13 prints: since 3.9 `__main__.__file__` is absolute, so the traceback shows
`/<scratch dir>/proj/scripts/run.py` etc. (verified; the line numbers 2, 9, 2 are right, with the
`# file:` comment kept as line 1). The conventions box promises "the real output". Fix: name the
scratch directory (`/tmp/py3`, as PY-1 uses `/tmp/py1`) in the conventions box and in the
solutions header, show those two paths in full, and show the module-folder path with an explicit,
documented abbreviation (`.../03-functions-modules-stdlib/numtools/cli.py`).

**FIXED — conventions box now names the scratch directory `/tmp/py3` and the `.../03-functions-modules-stdlib` abbreviation; the three tracebacks show `/tmp/py3/proj/scripts/run.py`, `.../03-functions-modules-stdlib/numtools/cli.py`, `/tmp/py3/geom/stats.py`; solutions header says the same. Re-run: every block matches up to the documented prefix.**

### M2 — solutions.md B3 (lines 278–319) · MEDIUM · solution does not do what it claims
The exercise asks "can you make that error message clearer?" and the prose answers "Detecting a
callable `lo` lets us write a better message", but `bounded(lo, hi)` has two required parameters,
so `@bounded` (one argument) raises Python's own `TypeError: bounded() missing 1 required positional
argument: 'hi'` before the body runs; the `callable(lo)` guard is dead code for exactly the case
asked about, and the stated output shows Python's message, not the custom one. Verified by running.
Fix: `def bounded(lo, hi=None)` with `if hi is None or callable(lo): raise TypeError(...)`; re-run;
paste the real output; correct the prose.

**FIXED — `def bounded(lo, hi=None)` with `if hi is None or callable(lo): raise TypeError("bounded needs two bounds: write @bounded(lo, hi), not @bounded")`; block re-run, real output pasted (`TypeError: bounded needs two bounds: ...`); prose rewritten to say why the guard needed the `None` default.**

### L1 — notes.md section 3.8 (line 815) · LOW · stated output
`$ ls numtools` shows seven files; once anything has been imported or run, the real listing also
has `__pycache__/` (it is there in this folder, gitignored). The listing is right only for a
package that has never been imported. Fix: one sentence after the block saying so.

**FIXED — one sentence after the listing says `__pycache__/` appears after the first import, is disposable and belongs in `.gitignore`.**

### L2 — solutions.md B2 (line 229) · LOW · wrong claim
"modern dicts remember insertion order, but they cannot cheaply move an existing key to the end":
`d[k] = d.pop(k)` moves a key to the end in O(1) (measured 0.14 us on 10^5 entries). The real
weakness of a plain dict as an LRU is the *oldest* entry: `next(iter(d))` walks over deleted
slots, measured 0.08 us fresh and 82 us after 90,000 pops from the front. Fix wording.

**FIXED — B2 now says a plain dict moves a key to the end in O(1) with `d[k] = d.pop(k)` and that its weakness is `next(iter(d))` walking deleted slots; `OrderedDict`'s linked list keeps both O(1).**

### L3 — solutions.md C1 (line 606) · LOW · wrong location of the failure
"`counted("fib")` would try to wrap a string and fail later, inside `functools.wraps`": `wraps`
tolerates missing attributes (it uses try/except and `getattr(..., {})`), so the failure is at
`f.__name__` with `AttributeError: 'str' object has no attribute '__name__'`, and it happens at
decoration time, not later. Verified. Fix wording.

**FIXED — C1 now says the unchecked `counted("fib")` fails at `f.__name__` with `AttributeError: 'str' object has no attribute '__name__'` at decoration time, and that `functools.wraps` tolerates missing attributes.**

### L4 — notes.md section 3.7 (line 554) · LOW · edge case
"the calls total 2n - 1" holds for n >= 1 (fib(0) is one call, not -1). Add "for n >= 1".

**FIXED — "the calls total $2n - 1$ (for $n \geq 1$)".**

### L5 — numtools/roots.py (line 14) · LOW · cross-reference
`bisect`'s docstring cites "PY-1 eq. 3.9" for the method; PY-1 (3.9) is the iteration-count bound
k = ceil(log2((b - a)/tol)), not the algorithm. Fix: "iteration count: PY-1 eq. (3.9)".

**FIXED — roots.py docstring: "(iteration count: PY-1 eq. 3.9)".**

### L6 — project.py (line 102) · LOW · idiom
`eval(out[0])` parses a subprocess's stdout, in the module whose Example 4.6 says "never pass
untrusted text to `eval`". `ast.literal_eval` is the idiom. Fix.

**FIXED — `import ast`; `ast.literal_eval(out[0])`. project.py exits 0 (3.9 s).**

### L7 — project.py · LOW · PEP 8 line length
Ten lines run 121–135 columns (101, 103, 131, 169, 207, 213, 219, 242, 341, 350) and more exceed
100. Reflowing correct code across the file is not a local fix.

**NOT FIXED: 104 lines of project.py exceed 100 columns; reflowing correct code across the whole file is not a local fix and risks breaking a passing script. Cosmetic.**

### L8 — timing numbers · LOW · note only
Section 3.5 (0.088 s, 0.075 s), Example 4.3 (median us), Example 4.7 (seconds and x-ratios), B7 and
C3 print different numbers on every run; orderings and ratios reproduce and the conventions box
says to trust ratios. No action.

**No action — the conventions box says to trust ratios, not microseconds; all orderings and ratios reproduce.**

## Verification after fixes

- `MPLBACKEND=Agg python3 project.py`: exit 0, 3.9 s, "All checks passed."; two figures written.
- All 46 notes.md blocks and all 19 solutions.md blocks with stated output reproduce it exactly
  (27 + 14 scripts, 10 + 1 REPL transcripts via doctest, 17 + 9 shell commands), up to the documented
  scratch-directory prefix in three tracebacks and the timing lines of L8. The `numtools/__main__.py`
  block is byte-identical to the file. `import numtools` works from a fresh interpreter in the module folder.
- notes.md word count 9,434 (PY addendum ceiling 9,500).
- Not touched: `curriculum/tools/reviews/PY-1.md` shows as modified in git status; that edit belongs
  to a concurrent review of another module, not to this one.

**STATUS: 10 findings (0 high, 2 medium, 8 low); 8 fixed, 1 not fixed (L7, cosmetic line length), 1 no action (L8, timing by convention); project.py exits 0; module consistent.**
