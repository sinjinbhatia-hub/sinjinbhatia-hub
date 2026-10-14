# Handoff — where the build stopped and how to continue

Paused 2026-10-14 at the user's request, so the next window starts with less context. This file is the only thing a new session needs to read first. After it, read `ROADMAP.md` for the plan and `research/WRITER-BRIEF.md` for the writing rules.

Branch: `claude/brave-knuth-64723m`. All work below is committed and pushed.

## Status (32 of 80 modules have all four files)

Every module listed as having four files has `notes.md`, `exercises.md`, `solutions.md` and `project.py`. At pause time, every `project.py` touched in the last round passed headless (F5, F6, LA-8, LA-9, EM-1, EM-2, ST-1, ST-2).

| State | Modules |
|---|---|
| **Done**: written, independently reviewed, fixes applied | F4, F5, LA-1, LA-2, LA-3, LA-4, LA-13, CM-1, CM-2, CM-3, CM-4, CM-6, CM-7, CM-8 |
| **Reviewed, fixes not applied**: the fix run was stopped and its issue list is lost with the old container | F1, LA-5, LA-6, LA-7, LA-8, LA-9, LA-10, LA-11, LA-12, LA-14 |
| **Written, not yet reviewed** | F2, F3, F6, CM-5, EM-1, EM-2, ST-1, ST-2 |
| **Partial**: only a draft `project.py` exists | EM-3 |
| **Not started** (briefs ready in `research/briefs/`) | EM-4–8, SM-1–7, QM-1–9, ST-3–9, ML-1–10, QC-1–9 (47 modules) |

The second and third rows (18 modules) need one combined review-and-fix pass each. LA-8 and LA-9 may already be partly fixed, because the stopped run edited them. The review re-derives the math anyway, so this is harmless.

## How to continue

The build runs through the Workflow tool with the script at `tools/build_modules_workflow.js`. Pass the script with `scriptPath`, because a script path from an old session will not resolve.

```json
{
  "combined": true,
  "modules": [
    {"code": "EM-3", "dir": "02-physics/b-electromagnetism/03-boundary-value-problems", "brief": "EM-3.md"},
    {"code": "LA-5", "dir": "01-linear-algebra/05-orthogonality-and-least-squares", "brief": "none", "mode": "verify"}
  ]
}
```

- `combined: true` runs one writer, then one agent that reviews and fixes. This is about a third cheaper than the old write → verify → fix pipeline. Use it.
- `mode: "verify"` skips writing and only reviews and fixes an existing module. F1, LA-1–12 and LA-14 have no brief file because they were written from `ROADMAP.md`. Pass `"brief": "none"`; the review prompt treats the brief as optional.
- Folder paths for every module are in the course READMEs and at the top of each brief.

**Concurrency: run 3 workflows at a time, with 3–4 modules each, and no more.** Eight parallel workflows used up a whole usage window in about 15 minutes earlier in this build. After each workflow returns, commit and push before you start the next one. A usage limit kills runs, and whatever is uncommitted at that point is at risk.

Suggested order (each line is roughly one workflow):

1. Verify: F1, F2, F3, F6
2. Verify: LA-5, LA-6, LA-7, LA-8
3. Verify: LA-9, LA-10, LA-11, LA-12, LA-14
4. Verify: CM-5, EM-1, EM-2, ST-1, ST-2
5. Write: EM-3, EM-4, EM-5, EM-6, then EM-7 and EM-8
6. Write: ST-3 to ST-9, three or four per workflow
7. Write: QM-1 to QM-9
8. Write: SM-1 to SM-7
9. Write: ML-1 to ML-10
10. Write: QC-1 to QC-9

The order follows prerequisites: E&M before QM, QM before SM-5 and QC, and ST before ML. Writers can then cite earlier modules by code instead of re-deriving them.

## Known issues to watch

- **Length.** Writers overshoot the notes target. The 6–8k word target with a 9.5k ceiling is in the prompt. Several finished modules run 11–19k words. That is acceptable where the content is dense, but don't let a new module pad.
- **Assert drift.** When a fix changes a number, the number has to change in the notes, the solutions and the assert together. The review prompt says so; spot-check it.
- **Conventions** (fixed; see `STYLE.md`):
  - The complex inner product is conjugate-linear in the first slot.
  - Qubit 0 is the leftmost tensor factor.
  - The QFT uses e^{+2πijk/N}.
  - Fourier uses the Osgood convention in F3 and symmetric 1/√(2π) in the physics modules.
  - Units are SI, with ℏ and k_B written out.

## Final pass, once all 80 modules exist

1. `bash tools/run_all_projects.sh`: every module must print PASS.
2. `python tools/py_to_ipynb.py` builds a `project.ipynb` next to every `project.py`.
3. Update the status column in `README.md`, then delete this file or mark it "complete".
4. Commit and push.
