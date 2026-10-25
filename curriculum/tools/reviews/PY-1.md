# Review: PY-1 — Objects, names, numbers, and control flow

Reviewer pass, adversarial. Module folder `curriculum/00-python/01-objects-names-and-control-flow/`.
Method: ran `project.py` (exit 0, 2 s, CPython 3.13.16); extracted every fenced `python`/`pycon`
block from notes.md (32 blocks) and solutions.md (18 blocks) into a scratch dir, ran every script
and diffed against the stated output, ran every REPL transcript through `doctest`; recomputed
every quoted number (ulp values, 0.1 fields, quadratic errors, summation errors, series tails,
H_500000/2, B6 mean tie, 5^22 < 2^53) in a scratch script; checked LA-13 and F1 cross-references
against those modules' notes; checked template sections, STYLE word bans, outdated patterns.

Result of the output check: all 32 notes blocks and all 18 solution blocks reproduce their stated
output exactly, except the two traceback blocks whose *line numbers* are off by one (item M4).
All derivations in §3.3–3.6 and §3.12 check step by step. No HIGH-severity findings.

## Findings

### M1 — notes.md §3.6 (line 532) and §10 (line 1204) · MEDIUM · attribution
"Goldberg 1991, Theorem 8 in the appendix": Theorem 8 (the Kahan summation bound) is in the
"Errors in Summation" part of the paper's "The Details" section, not in an appendix (the Sun
reprint's appendix is Priest's separate note). §10 also says Goldberg "§1-2" holds "the Kahan
theorem of §3.6", which it does not. Fix: cite "Goldberg 1991, Theorem 8" without "appendix";
in §10 point §1–2 at the rounding material and Theorem 8 (in "The Details") at the Kahan bound.

### M2 — notes.md §0 (line 9), §3.5 step 3 (line 440), §5 (line 1156) · MEDIUM · cross-reference
LA-13's cancellation bound is its eq. (3.3); its (3.4) is the *condition number of subtraction*
κ_sub = (|x|+|y|)/|x−y|. LA-13 §3.2 is the scalar bound (same as PY-1's), not "cancellation for
matrices"; the condition number of a linear system is LA-13 §3.3, eq. (3.7). Fix the three
references to name those objects.

### M3 — solutions.md B5 explanation (line ~462) · MEDIUM · wrong reasoning
"(a) works because the tail after a term is about the next term times 1/(k+1)" is wrong: the tail
from the first omitted term onward is about that term itself (later terms shrink by 1/(k+1),
1/((k+1)(k+2)), ...). Measured: first omitted term 1/15! = 7.65e-13, error 8.15e-13. Fix wording.

### M4 — notes.md Example 4.1 (lines 945–958) and solutions.md A4 (lines ~133–140) · MEDIUM · stated output
The traceback line numbers (25/21/15/4 and 10/6) are those of a file saved *without* the
`# file: name.py` first line, but the conventions box tells the reader to save the block as shown,
which gives 26/22/16/5 and 11/7 (verified on 3.13). Fix: paste the tracebacks produced with the
comment line kept.

### L1 — notes.md §3.4 (line 280) · LOW · precision
"a fraction pattern `1001 1001 ...` that repeats and is cut off after 52 bits": the stored field is
rounded, not truncated (last stored nibble is `1010`, project asserts `"1001"*12 + "1010"`), and
that rounding-up is exactly why the stored 0.1 is *above* 0.1. Fix wording.

### L2 — notes.md §3.6 (line 532) · LOW · overstatement
"The method fails if a term can exceed the running total": Goldberg's bound still holds; what is
lost is the compensation for that step (Neumaier's variant keeps it). Soften.

### L3 — notes.md §3.6 (line 532) · LOW · claim vs project
"np.sum ... error growing like log n (PY-6; the project measures it)": the project compares
np.sum at one n = 10^7 only; it does not measure growth. Fix wording.

### L4 — notes.md §3.8 (line 591) · LOW · completeness
Format mini-language summary omits `#` though the example uses `{255:#x}`. Add `[#]`.

### L5 — notes.md Example 4.4 (line 1074) · LOW · wording
"24 integer digits, of which about 9 mean anything": the 9 are the digits of the defined constant;
digits 10–17 are still the correctly rounded double, only the tail is the binary artefact. Reword.

### L6 — notes.md §7 (lines 1166–1172) · LOW · consistency
The project list omits project Part 1 (the floor rule on all sign combinations). Add it.

### L7 — notes.md §3.5/§3.6/§3.12 · LOW · presentation
Build-it-yourself sequences are numbered 1, 2, 3, 6, 4, 5 in reading order (brief's numbering kept);
sequences 3 and 6 have no title. Renumber 6→4, 4→5, 5→6 and title them. No other file refers to
sequences 4–6 (checked with grep).

### L8 — solutions.md header (line 3) · LOW · consistency
Says code lines "may run to 100 columns"; six code lines are 101–108 columns (b1, b2, b3, b5, b6,
c3). Reflow those lines; re-run the blocks.

### L9 — solutions.md C2(c) · LOW · spec
Exercise asks for 100–150 words; the paragraph has 92. Add a sentence.

### L10 — solutions.md B5(c) · LOW · presentation
Prints "error 1.00e-06" (`.2e`) for a value that is 9.99998e-07 and asserted `< 1e-6`; a reader sees
the bound apparently violated. Print more digits; update the stated output.

### L11 — project.py line ~431 · LOW · idiom
`__import__("cmath").exp(...)` inline; import `cmath` at the top like the other modules.

## Checked and found correct (no action)
- All numeric claims: ε, u, ulp(1e3), ulp(1e6), ulp(1e8), ulp(1e16), sqrt(1e16−4), 0.1/0.2/0.3
  Decimal and hex values, 2^55 denominator, tie analysis, naive/stable quadratic roots and errors,
  κ·u = 1.11, naive sum error 1.610e-4, Newton error sequence, bisection counts 20/10/34/40,
  Horner counts, Example 4.4 table, B-series tails, C3 running-total −6.85 and ulp 8.9e-16,
  B6 mean 7.718500000000001, 5^22 < 2^53 < 5^23.
- Derivations: floor rule (3.1), ulp (3.3), rounding model (3.4), isclose (3.5), cancellation bound
  (3.6), rationalized root (3.7) and Vieta, Fast2Sum (3.8), bisection count (3.9), Newton (3.10)
  and quadratic convergence (3.11), Horner (3.12) and n(n+1)/2 count.
- Semantics: small-int cache −5..256, `is` vs `==`, `+=` on lists, UnboundLocalError compile-time
  rule, loop `else`, keyword-only `*`, nan unordered, `nan in [nan]`, `**`/`math.exp` raise vs `*`
  returns inf, round-half-even, `-i` after a crash, 3.12 `sum` is Neumaier.
- Template: §0–§10 with template headings; conventions box; Try-it-first/What-just-happened on
  every sequence; 6 worked examples; A(6)/B(7)/C(3) = 16 exercises, all with complete solutions;
  B7 extends project.py; §8 specific; §10 matches the brief's five sources.
- Code: stdlib + numpy only; headless; exit 0 in 2 s; no `%`-formatting, `str.format`, bare
  `except:`, star imports, or mutable defaults outside the trap examples; asserts test the exact
  quoted values (not weaker bounds); solutions' asserts all pass.
- History: Dekker 1971 (Fast2Sum), Kahan, Neumaier, Goldberg 1991, Batchelder PyCon 2015, PEP 8,
  3.12 sum change — all consistent with my knowledge; only the "appendix" location was wrong (M1).
