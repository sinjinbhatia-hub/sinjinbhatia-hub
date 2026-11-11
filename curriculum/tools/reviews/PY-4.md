# Review: PY-4 Classes and the object model

Reviewer pass started 2026-10-25. Module folder: `curriculum/00-python/04-classes-and-the-object-model/`.
Brief: `curriculum/research/briefs/PY-4.md` (has REVIEW FOCUS). Scratch: /tmp/claude-0/py4/ (harness.py, check.py).

Baseline: `MPLBACKEND=Agg python3 project.py` exits 0 in 2.1 s (Python 3.13.16, numpy 2.5.3, sympy 1.14.0, matplotlib 3.11.2);
every printed number matches the notes (ten-function table, floors, Example 4.5, integrator ratios, slot sizes).

Method: (1) every fenced `python`/`pycon` block in notes.md (35) and solutions.md (16) extracted and run by
/tmp/claude-0/py4/harness.py (notes in one cumulative namespace, with the incremental `class Vec:`/`Rational`/`Matrix` blocks
merged into the existing class as the conventions box says; solutions in one namespace after the preamble's imports), actual
stdout diffed against the stated output; (2) every derivation and hand number recomputed independently in
/tmp/claude-0/py4/check.py with sympy/exact rationals: (3.1), (3.2), (3.3), the product/quotient/chain-rule reading, integer
powers, the 2x2 matrix copy of eps, (e1+e2)^2 = 2 e1 e2; f'(1.3) to 30 digits, u/u'/v/v'/y' and the hand arithmetic
0.27059/7.2361 of Example 4.4, the truncation sizes behind the §1 table (1/2|f''|h = 3.43e-3, 1/6|f'''|h^2 = 1.17e-5); the
pivoted LU of Example 4.5 in exact arithmetic (perm, L, U, Pb, y, x); the Newton iteration counts 5 and 6 and the Kepler root;
the polar Jacobian and det = r; the trapezoid/Simpson error table and its ratios 4/16; A5(a),(b) in Fractions; B4's three f''
values and the mixed partial 2x + cos(xy) - xy sin(xy); B5's determinants; B6's (x-1)(x-2)(x-3); B7's f, f'', f''' at 1 and
h* = sqrt(2 eps f/f''), (3 eps f/f''')^(1/3), E* = sqrt(2 eps f f''), 1.5 eps f/h* (minimised symbolically); C2's two-point
Gauss-Legendre error and int sqrt = 2/3 — all agree; (3) the dispatch rules of §3.3 and Example 4.3 checked against CPython's
`SLOT1BINFULL`/`method_is_overloaded` (same type: no reflected call; subclass overriding `__radd__`, including when the base has
none, goes first, then declines to the left `__add__`) and the data-model reference §3.3.8; the `__eq__`/`__hash__` statements
(implicit `__hash__ = None`, `!=` from `__eq__`, identity fallback, dataclass `eq`/`frozen`/`unsafe_hash` rules) checked by
execution; (4) consistency, template, cross-references (LA-1 §3.1/§3.5/Theorem 3.11; LA-2 Theorem 3.13, (3.3), (3.4), (3.7);
LA-11 shift/Jordan block; LA-14 Jacobian; CM-1 `euler_step`/`rk4_step`/`symplectic_euler_step(a, x, v, h)`; PY-1 §3.5 and
`newton(f, df, x0)` with root 0.7390851332151607; PY-2 §3.4/§3.7/§3.8; PY-3 §3.4 closures; PY-5/6/9/10, ML-2/ML-7 briefs) —
all exist and say what PY-4 says they say; (5) REVIEW FOCUS: mutable class attributes appear only as the named trap (Tally,
Account); inheritance appears as the named anti-example (`CountingDict`), as the compared design (`Integrator`) and as
"subclass so the project's class is untouched" in solutions; brief's allowed imports only; runtime 2.1 s.

Stated outputs: all 35 notes blocks and all 16 solutions blocks reproduce exactly, except (i) the ids in Example 4.1 (conventions
box says they change) and (ii) the timing tables of C1 and C3, which differ in value (ratio 263 vs 279; slope 2.90 vs 2.96) but
not in any conclusion. Exercise and solution numbering match (A1-A6, B1-B7, C1-C3); every solution is complete and runs;
no "similarly"/"left to the reader"/"it can be shown"/"simply"/"obviously" anywhere. Template: §0-§10 with the template
headings, conventions box, 7 worked examples (4.7 from mechanics), 6 A / 7 B / 3 C, B1/B7 extend project.py, §8 specific,
§10 matches the brief. Word count 9,406 (under the PY ceiling of 9,500; above the 6,000-8,000 target, not changed).
History: no attributions beyond Baydin et al. §3.1 (forward mode; correct), Horner, Kepler, Gauss-Legendre, Cauchy-Schwarz,
Lagrange's identity — all correct.

Format: **ID** | file:location | severity | what is wrong | fix | status

## Findings

**1** | solutions.md A6, last sentence of the prose before the code ("A list field could not be hashed even in a frozen class,
so it must be a tuple") and the code (three changes where the exercise asks for two) | medium | False: a field declared with
`compare=False` is also left out of the generated `__hash__`, so `@dataclass(order=True, frozen=True)` with
`tags: list = field(default_factory=list, compare=False)` hashes fine (verified: `hash(V(1, 2))` works). The exercise asks for
"two changes"; the solution makes three and justifies the third with the wrong reason. | give the two-change version
(`frozen=True`, `compare=False`) as the answer with the correct explanation, keep the tuple as a recommendation ("immutable all
the way down, and `default=()` needs no factory"), assert both versions hash | FIXED

**2** | project.py Part 3, the §1-table loop (`0.1 < ratio < 10` forward, `0.01 < ratio < 100` central for all four h) | medium |
Weaker than the notes claim. The entries at h = 1e-2 and 1e-5 are truncation-dominated and deterministic to three digits
(1/2|f''|h = 3.43e-3 and 3.43e-6; 1/6|f'''|h^2 = 1.17e-5; central at 1e-5 is 5.19e-12 and also reproduced); only the 1e-8 and
1e-11 rows are rounding noise. A factor-10/100 window would pass a table that is wrong in every digit. | assert the h = 1e-2 and
1e-5 rows (forward) and the h = 1e-2 row (central) to 1% relative; keep the loose window for the noise rows, with a comment | FIXED

**3** | notes.md §3.7, ■ after Theorem 3.1 ("Where F is not differentiable (`abs` at 0, a branch), you get the derivative of the
branch taken") | low | Implies the module's `Dual` handles `abs` and comparisons. It defines neither `__abs__` nor `__lt__`, so
`abs(Dual(0.0, 1.0))` and `if x > 0` raise `TypeError` (verified) rather than silently differentiating a branch. | say so: an AD
system that defines them returns the derivative of the branch taken; this `Dual` raises `TypeError`, which is the safer
half-working behaviour of §8 | FIXED

**4** | notes.md §3.3 ("Step 4 adds it, with division, length, dot and cross products") | low | There is no step 4; the
build-it-yourself blocks are (1a)-(1d) and the next block is unlabeled | "The next block adds it" | FIXED

**5** | project.py Part 2 (`for p, q, r_, s in rng.integers(-50, 51, size=(1000, 4)): if q == 0 or s == 0: continue`) | low |
Skips pairs with a zero denominator, so 978 of the promised 1000 random pairs are checked (notes §3.6/§7 and the brief say
1000) | draw denominators from 1..50 with a random sign so all 1000 pairs are checked | FIXED

**6** | project.py Part 5 (header "dataclasses, inheritance, protocols ... (notes §3.9-§3.11)") | low | Nothing in Part 5
checks §3.10's protocol/ABC claims: `typing.Protocol` is not imported, `Reading.__array__`, the `Sequence` mixins
(`index`, `count`, `in`, `reversed`), the abstract-instantiation `TypeError`, and `isinstance(Polynomial, Evaluable)` are all
unasserted although the notes print them | add those asserts to Part 5 | FIXED

**7** | notes.md §3.8 `Matrix.lu` and project.py `Matrix.lu` (inner loop `for j in range(k, n)`) | low | The eliminated entry
U[i][k] is computed as U[i][k] - (U[i][k]/U[k][k])*U[k][k], which leaves rounding residue (1.1e-16 measured) below the
diagonal, so the returned U is not exactly upper triangular as LA-2 defines it (Example 4.5 is exact only because its entries
are binary fractions) | start the inner loop at k + 1 and set U[i][k] = 0.0 explicitly, in both files | FIXED

**8** | solutions.md A3 prose ("the trace shows rule 2 of §3.3: the first call is `0 + v`") | low | The first trace line
(`0 + v`) is rule 1 (left operand declines, reflected method runs); rule 2 (subclass asked first) is the second trace line |
reword to name both rules in order | FIXED

**9** | exercises.md A5 ("evaluate with `Dual(2, 1)` by hand: (a) ...; (b) g(x) = x^3 (1/x) at x = 3") | low | (b) is seeded
with `Dual(3, 1)`, not `Dual(2, 1)` | "with the seed `Dual(x, 1)`: (a) ... at x = 2; (b) ... at x = 3" | FIXED

**10** | solutions.md B7 prose (central model `E(h) ~ 1/6|f'''|h^2 + eps|f|/h`) | low | The exercise's own assumption is a
rounding error of about eps|f| in the difference of two function values; the central quotient divides by 2h, so the term is
eps|f|/(2h). The factor 2 moves h* by 2^(1/3) = 1.26 and E* by a similar factor, below the model's accuracy, but the text
should say it is taking the crude bound | add one parenthetical sentence; numbers unchanged | FIXED

**11** | solutions.md C1 | low | The exercise asks for a check "against two forward `Dual` passes and SymPy"; the solution checks
against the Dual passes and the hand formula (y + cos x, x), not SymPy | add a two-line SymPy gradient check and re-run | FIXED

**12** | notes.md §3.1 (0), §3.2 (1a), §3.4 (1c), §3.5 (1d), §3.7 (3), §3.8 (4), §3.9 (5) | low (template) | These ▶ Build it
yourself prompts have no ■ What just happened; the PY addendum says a ▶ that asks for a prediction is followed by a ■ that
explains the gap. Of these only (3) asks for a prediction; (0), (1a) and (4) have the explanatory paragraph unlabeled; (1c),
(1d) and (5) have none. | label the existing paragraphs for (0), (1a), (4); add one short ■ each for (3) and (5). (1c) and (1d)
are left without one on purpose: they ask for no prediction, their outputs are self-explanatory, and notes.md sits at the PY
word ceiling (9,406 before the review; 9,498 after, ceiling 9,500), so every added sentence had to earn its place. | FIXED
(partial by design: (1c), (1d) unchanged)

**13** | solutions.md A2 `FixedMoney`, C1 `Var`, C2 `Counted` | low (REVIEW FOCUS: missing `__repr__`) | Classes offered as
model answers with no `__repr__`, against the module's own rule of thumb (§3.11 "write `__repr__` first") | add one-line
`__repr__`s; outputs unchanged | FIXED

**14** | solutions.md preamble | trivial | C1 and C3 print timings that differ by machine (263 vs 279, 3762x vs 3293x here);
the preamble says only "the outputs shown are the real ones" | add "timings vary by machine; trust the ratios and slopes" | FIXED

**15** | solutions.md B2 last paragraph ("Neither has `<`") | trivial | The exercise asks which comparison operators the
built-in lacks: all four of `<`, `<=`, `>`, `>=` | say so | FIXED

## Verification after fixes

`MPLBACKEND=Agg python3 project.py` exits 0 in 2.4 s (Part 2 now checks all 1000 pairs; Part 3 holds the deterministic rows of
the §1 table to 1%; Part 5 asserts the `Sequence` mixins, the abstract-instantiation `TypeError`, the `Evaluable` protocol and
`Reading.__array__`; `lu()` returns an exactly upper-triangular U, measured 0.0 below the diagonal). Harness re-run: all 35 notes
blocks reproduce (Example 4.1 ids excepted) and all 16 solutions blocks reproduce (C1/C3 timings excepted); the edited A6
(two-change `FrozenVersion` plus `TupleVersion`), A2/C1/C2 (`__repr__`s) and C1 (SymPy gradient) blocks pass their asserts.
notes.md: 9,499 words (ceiling 9,500); 13 ▶, 11 ■ (every ▶ that asks for a prediction is followed by a ■). Files touched:
notes.md, exercises.md, solutions.md, project.py in the module folder, and this file.

STATUS: 15 findings (0 high, 2 medium, 13 low/trivial); 15 fixed, item 12 narrowed on re-verification ((1c), (1d) left as is);
project.py passes; nothing outstanding.
