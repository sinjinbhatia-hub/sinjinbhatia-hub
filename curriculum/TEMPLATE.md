# Module Template

Every module in every course uses this structure. Files per module folder `NN-short-name/`:

| File | Purpose |
|---|---|
| `notes.md` | The lesson. Sections 0–9 below, in this order, with these headings. |
| `exercises.md` | Problems only, graded A/B/C. No solutions. |
| `solutions.md` | Full worked solutions to every exercise, step by step. |
| `project.py` | Runnable Python with `# %%` cell markers. Must run headless (`MPLBACKEND=Agg python3 project.py`) in under ~60 s and exit 0. Saves figures to `figures/` next to it. |

## `notes.md` structure

```markdown
# <Course code>-<N>: <Title>

**Time:** ~N hours · **Prerequisites:** <module codes> · **Unlocks:** <module codes>

> **Conventions in this module.** <one-paragraph box: notation, sign conventions, units>

## 0. Where this sits
Two or three sentences. What you must already know. What this module makes possible later (name the modules).

## 1. The problem that forced this idea
History and motivation. What concrete question couldn't be answered without this, who hit it, what they tried first. 2–5 paragraphs. No formalism yet.

## 2. Intuition
The geometric or physical picture. Pictures in words, small examples in R² or with two coins or one spring. ASCII sketches welcome. The reader should be able to *guess* the main theorem by the end of this section.

## 3. Formal development
Definitions, theorems, derivations. Rules:
- Every nontrivial result is derived step by step. Number the steps.
- Before each key derivation, a **▶ Try it first** prompt: one or two sentences telling the reader what to attempt before reading on.
- After each derivation, a **■ What just happened** sentence or two in plain language.
- State what is being assumed. If a proof is omitted for length, say so and point to a source.

## 4. Worked examples
At least three, increasing in difficulty, fully computed. Where numbers appear, they should be checkable by hand or by the project code.

## 5. Why it matters and where it connects
Forward and sideways links: which later modules use this, how the same mathematics appears in other fields in this curriculum. Be concrete ("the inertia tensor in CM-7 is a symmetric matrix; this module's spectral theorem is why principal axes exist").

## 6. Exercises
Pointer to `exercises.md`, with a one-line summary of the A/B/C split and which exercises are essential.

## 7. Coding project
What `project.py` builds, what the reader should modify or extend, and what result they should see. Include the "build it yourself" target: the part of the code the reader should write before looking.

## 8. Where you'll struggle
Honest, specific. Which idea confuses most people, which convention gets inverted, which exercise is harder than it looks. Give the fix, not just the warning.

## 9. Check yourself
8–12 questions answerable in a sentence or a line of math, no solutions. If the reader can't answer them closed-book, they should not move on.

## 10. Pair with
Textbook chapters/sections and video lectures for this module (from the roadmap). One line each on what each source does better than these notes.
```

## `exercises.md` structure

```markdown
# <Code>-<N> Exercises

**A. Warm-up** (do all; 5–10 minutes each)
A1. ...
**B. Core** (do all; these are the module)
B1. ...
**C. Challenge** (do at least one)
C1. ...
```

Target: 4–6 A, 5–8 B, 2–4 C. Mix of computation, proof, and "explain in words". At least one exercise per module should require writing code or extending `project.py`.

## `solutions.md` structure

Mirror the exercise numbering. Every solution is complete: no "left to the reader", no "similarly". Where a proof is required, write the proof. Where computation is required, show intermediate values.

## `project.py` structure

```python
# %% [markdown]
# # <Code>-<N> Project: <title>
# One paragraph: what this builds and what the reader should learn from it.

# %%
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)

# %% [markdown]
# ## Part 1: <name>
# Explain what this part does and the math it implements (cite the notes section).

# %%
def something(...):
    """Docstring stating the math."""
    ...

# ... asserts that check results against theory or against NumPy/SciPy reference implementations

# %% [markdown]
# ## Your turn
# Explicit prompts for what the reader should modify or implement.

# %%
if __name__ == "__main__":
    print("All checks passed.")
```

Rules: deterministic (seeded), no network, no heavy dependencies beyond numpy/scipy/matplotlib/sympy, every claim in the notes that can be checked numerically is checked with an `assert` or a printed comparison. Figures saved, never shown.
