# Handoff — where the build stands and how to continue

Last updated 2026-11-09 20:45 UTC, at the start of round 17. The branch is `claude/brave-knuth-64723m`. Read this first, then `ROADMAP.md` for the plan and `research/WRITER-BRIEF.md` for the writing rules.

## Status (42 of 80 math modules done; Python course complete, 10 of 10)

"Done" means all four files exist (`notes.md`, `exercises.md`, `solutions.md`, `project.py`), an independent reviewer re-derived the math and recomputed the numbers, the fixes were applied, and the script passes headless. Review notes for each module are in `tools/reviews/<code>.md`. A module with no file there has not been reviewed.

| State | Modules |
|---|---|
| **Done** (42) | F1–F6 · LA-1–14 · CM-1–8 · EM-1–8 · ST-1–6 |
| **Done, Python** (10) | PY-1–10, reviewed with every code block executed |
| **Written, awaiting review** (8) | QM-1–5 · ST-7–9 (Opus; scripts pass; reviews in flight in round 17) |
| **Being written in round 17** (4) | QM-6, QM-7, QM-8, QM-9 |
| **Not started** (26) | SM-1–7 · ML-1–10 · QC-1–9 |

Until a module is in a Done row, treat its derivations and numbers as unverified.

Courses fully done: Python, Foundations, Linear Algebra, Classical Mechanics, E&M.

## What was lost on 2026-11-04, and the fix

The usage limit hit mid-round and the container was reclaimed before the last edits were committed. Lost: the review fixes for QM-1 and QM-2 (their reviews had finished; they are being redone) and partial drafts of QM-6–9. Cause: the autosave script lived in `/tmp` and ran every 10 minutes. It now lives at `tools/autosave.sh`, runs every 3 minutes, and must be started at the beginning of every round:

```bash
cd <repo root>
pip install -q numpy scipy matplotlib sympy pandas pytest numba   # a fresh container has none of these
setsid nohup bash curriculum/tools/autosave.sh > /tmp/autosave.log 2>&1 < /dev/null & disown
```

Also commit and push by hand each time an agent reports. Do not end a turn with uncommitted work.

## How to run a round

Launch one background agent per module with the Agent tool (`subagent_type: general-purpose`). The standing instructions are in two files; the launch message only names the module:

- **Write:** "Read ROOT/tools/AGENT-WRITE.md and follow it exactly, with CODE = QM-6, DIR = 02-physics/d-quantum-mechanics/06-hydrogen-atom, BRIEF = QM-6.md." Pass `model: "opus"` (the user's choice since round 10; Sonnet also works and costs less per module against the usage limit).
- **Review and fix:** "Read ROOT/tools/AGENT-REVIEW.md and follow it exactly, with CODE = …, DIR = …, BRIEF = …." No model override, so the reviewer runs on the session model.

Tell every agent: no git commands, write only inside its module folder (and its review file), save files as soon as each is drafted. `ROOT` is the absolute path of this `curriculum/` folder. Folder paths for every module are on the "Folder:" line of `research/briefs/<code>.md`. Modules without a brief (F1, LA-1–12, LA-14) take BRIEF = none.

`tools/build_modules_workflow.js` holds the same prompts for the Workflow tool (args `{combined: true, writerModel, modules: [{code, dir, brief, mode}]}` with mode `write` or `verify`). Use it when that tool is available; the two AGENT files were generated from it.

A reviewer that is interrupted leaves its findings in `tools/reviews/<code>.md`; the next reviewer for that module reads the file and continues.

## Sizing

There is no usage signal readable from inside the session. Measured: a review costs 2–4M cache-read tokens, an Opus write 12–17M. One 5-hour window has covered about 8 Opus writes plus 3 reviews (round 16) before the limit killed everything in flight. Eight agents at a time is the ceiling on this 4-CPU machine. The user's standing instructions: use the window, but when the limit hits, do not restart on the reset — wait for "continue".

## Build order (prerequisites first)

1. Review QM-1–5 and ST-7–9. Write QM-6–9.
2. Review QM-6–9. Write SM-1–7 (SM-5 needs QM-5).
3. Review SM. Write ML-1–10 (needs ST-1–5, LA-14).
4. Review ML. Write QC-1–9 (needs QM-4, QM-7).
5. Review QC. Final pass.

## Known issues

- **Length.** Almost every module runs over the STYLE.md targets (notes 9–15k words, scripts 500–900 lines). Reviewers decline to trim correct, brief-mandated material. Accepted.
- **Number drift.** ST-5's review found quoted numbers in the notes that no longer matched the script, because the writer changed the data generation after writing the prose. Reviewers of simulation-heavy modules must pin every quoted number with an assert.
- **Signature.** EM-8 uses Griffiths' (−,+,+,+) metric, so ∂_νF^μν = μ₀J^μ. QM-9 and anything relativistic should state its convention.
- **Conventions** are fixed in `STYLE.md` (complex inner product conjugate-linear in the first slot; qubit 0 leftmost; QFT e^{+2πijk/N}; Osgood Fourier convention in F3, symmetric 1/√(2π) in physics; SI with ℏ and k_B explicit).
- pandas 3 is installed (copy-on-write always on); PY-9 is written for it.

## Final pass, once all 80 modules are done

1. `bash tools/run_all_projects.sh` — every module must print PASS.
2. `python3 tools/py_to_ipynb.py` — builds a `project.ipynb` beside every `project.py`.
3. Update the status column in `README.md`; replace this file with a one-line "complete" note.
4. Commit and push.
