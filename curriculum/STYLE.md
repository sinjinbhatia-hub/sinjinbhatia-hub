# Author Style Guide

Read this before writing any module. It applies to every course.

## Who you're writing for

One specific student: second-year mechanical engineering and applied math undergraduate, physics concentration. Solid in calculus through multivariable, ODEs, some PDEs. Partial Fourier. Weak linear algebra at the start (rebuilt in course 01). Can write working Python and wants to become genuinely good at it (course 00-python builds that fluency and goes first). Learns best by building the picture, deriving the result, then coding it. Wants bluntness about difficulty, not reassurance.

Write to that one person. "You" throughout.

## Non-negotiables

1. **Intuition before formalism.** Section 2 of every module must let the reader guess the main result before section 3 states it.
2. **Derive, don't declare.** Every central result is derived step by step with numbered steps. "It can be shown" is banned for anything central to the module. If something genuinely can't be derived at this level (e.g., Jordan form existence in full generality), say exactly what is being skipped and where the full proof lives.
3. **Try-it-first prompts** before each key derivation. **What-just-happened** summaries after.
4. **Why it exists.** Section 1 names the concrete problem that forced the idea into existence, with history where it's real and interesting. No invented history.
5. **Honest difficulty.** Section 8 says where people fail and what to do about it. No "with a little practice this becomes natural."
6. **Everything checkable.** Numbers in worked examples are verified by the project code. Derivations are checked symbolically (SymPy) where feasible. If you state a numerical value, it should match the code's output.
7. **Cross-links are concrete.** Name the module code and the specific object ("the Gram matrix $K$ in ML-4 is positive semidefinite by LA-8 Theorem 3.4").

## Voice

- Direct, warm, plain. Short sentences. Explain like a good TA at a whiteboard, not like a textbook preface.
- No hedging filler ("it is interesting to note that"), no cheerleading ("this beautiful result"), no "simply"/"just"/"obviously"/"clearly."
- Use "we" for shared derivation steps ("we expand the left side") and "you" for instructions to the reader.
- Humor only if it clarifies.

## Mathematics

- LaTeX in Markdown: `$...$` inline, `$$...$$` display. Number equations you refer back to with `\tag{3.2}`.
- Define every symbol the first time it appears in a module, even standard ones.
- Prefer concrete dimensions first ($\mathbb{R}^2$, $2\times 2$) then generalize.
- State theorems as **Theorem 3.1 (Name).** with hypotheses explicit. Proofs end with ∎.
- Sign conventions, normalization conventions, and index conventions go in the conventions box at the top *and* are repeated at the point of use.
- When a result has a name, use it (Rank–Nullity, Gauss–Markov, Eckart–Young) so the reader can look it up.

## Worked examples

- At least three per module, increasing difficulty, fully computed, intermediate values shown.
- At least one should be "from" another field in this curriculum (a mechanics example in an LA module, a statistics example in a probability module).

## Exercises and solutions

- A (warm-up): mechanical, 5–10 min. B (core): the module's real content; mixture of compute, prove, explain. C (challenge): one step beyond; may require reading ahead or coding.
- Solutions are complete. Write every step. A reader who is stuck should be unstuck by the solution alone.
- At least one exercise requires modifying or extending `project.py`.

## Code

- Python 3.11+, NumPy, SciPy, Matplotlib, SymPy only (the PY course additionally uses pandas and pytest where its briefs say so). No sklearn/torch unless the module is explicitly about using them (ML-5 onward may use scikit-learn; ML-7 onward may use PyTorch for comparison only, with a NumPy implementation as the primary).
- Implement the algorithm yourself first, then compare to the library.
- `# %%` cells, markdown cells explaining the math, asserts checking theory against computation. Seeded RNG. Headless matplotlib. Figures saved to `figures/`.
- Must run clean with `MPLBACKEND=Agg python3 project.py` in under about 60 seconds and exit 0.
- Docstrings state the mathematical object being computed, with the equation number from the notes if applicable.

## Length

Notes: roughly 3,500–7,000 words depending on module weight (the hours in the roadmap are the guide). Exercises: 12–18 problems. Solutions: as long as they need to be. Project: 150–400 lines.

## Things that get a module sent back

- A derivation replaced by "it can be shown."
- Section 2 that is just section 3 with fewer symbols.
- Worked example numbers that the project code contradicts.
- A solutions file with "similarly" or "left as an exercise."
- Code that doesn't run headless, or imports a library not on the list.
- Invented history or invented attributions.
- Cross-links that name a field but not a module code and object.
