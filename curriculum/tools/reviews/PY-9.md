# Review: PY-9 pandas and working with data

Module: curriculum/00-python/09-pandas-and-data/ (notes.md, exercises.md, solutions.md, project.py)
Brief: curriculum/research/briefs/PY-9.md. Reviewer: independent adversarial pass, 2026-10-31.
Environment: Python 3.13.16, pandas 3.0.5, NumPy 2.5.3, Matplotlib 3.11.2, SciPy 1.18.1, no pyarrow (same as the module's stated run).
Method: ran project.py headless with -W default and timing; extracted and ran every ```python block of notes.md (25 blocks) and solutions.md (15 blocks) in one session each (isolated copy of the folder for solutions) and diffed stated ```text output against real output; recomputed quoted numbers in scratch scripts under /tmp/claude-0/py9; checked semantics claims (loc inclusivity, alignment, merge, resample closed/label, rolling, CoW) against pandas 3.0.5 behaviour and the user guide; checked cross-references against the referenced modules' files; applied the brief's REVIEW FOCUS list.

Format: one item per finding. Severity: high (wrong output/number/claim, failing or weakened assert, missing required section) / medium (misleading or inconsistent) / low (pedagogy, polish). Status: OPEN / FIXED / NOT FIXED: reason.

## Findings

