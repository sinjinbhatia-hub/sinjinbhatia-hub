# Handoff — where the build stands and how to continue

Last updated 2026-10-24 01:10 UTC, after round 8. The branch is `claude/brave-knuth-64723m`; everything below is committed and pushed. Read this first, then `ROADMAP.md` for the plan and `research/WRITER-BRIEF.md` for the writing rules.

## Status (37 of 80 math modules done; Python course 4 of 10 written)

"Done" means all four files exist (`notes.md`, `exercises.md`, `solutions.md`, `project.py`), an independent reviewer re-derived the math and recomputed the numbers, the fixes were applied, and the script passes headless. Review notes for each module are in `tools/reviews/<code>.md`.

| State | Modules |
|---|---|
| **Done** (37) | F1–F6 · LA-1–14 · CM-1–8 · EM-1–5 · ST-1–4 |
| **Written, awaiting review** (1 + 4) | ST-5 · PY-1, PY-2, PY-3, PY-4 (written by Sonnet; scripts pass; not yet independently checked) |
| **Not started** (41 + 6) | EM-6, EM-7, EM-8 · SM-1–7 · QM-1–9 · ST-6–9 · ML-1–10 · QC-1–9 · PY-5–10 |

**Priority (user's request, 2026-10-23): finish the Python course (`00-python/`, briefs `research/briefs/PY-*.md`) first.** Its briefs end with a REVIEW FOCUS section the reviewer applies (run every code block; a wrong stated output is high severity). The PY writer addendum is at the end of `research/WRITER-BRIEF.md`.

Until a module is in the Done row, treat its derivations and numbers as unverified.

## How to run a round

The build runs through the Workflow tool with `tools/build_modules_workflow.js` (pass it as `scriptPath`). Arguments:

```json
{
  "combined": true,
  "writerModel": "sonnet",
  "modules": [
    {"code": "ST-3", "dir": "03-statistical-modeling/03-linear-regression", "brief": "ST-3.md", "mode": "verify"},
    {"code": "EM-6", "dir": "02-physics/b-electromagnetism/06-electrodynamics-maxwell", "brief": "EM-6.md", "mode": "write"},
    {"code": "EM-7", "dir": "02-physics/b-electromagnetism/07-electromagnetic-waves", "brief": "EM-7.md"}
  ]
}
```

- `mode` omitted: write, then review+fix. `"verify"`: review+fix an existing module. `"write"`: write only, review deferred.
- `writerModel: "sonnet"` puts the writers on Sonnet; reviewers always use the session model (Fable). Keep this. Measured from transcripts: a from-scratch write costs 10–17M cache-read tokens on either model, a review 2–4M, and Sonnet tokens count against the usage limit at a fraction of Fable's rate.
- Modules without a brief file (F1, LA-1–12, LA-14) take `"brief": "none"`. Folder paths for every module are in `research/briefs/*.md` ("Folder:" line) and the course READMEs.
- Reviewers write their findings to `tools/reviews/<code>.md` as they go. If a reviewer is killed, the next one resumes from that file.
- `tokenCap` exists in the script but is not useful: `budget.spent()` does not track the usage limit.

## Sizing a round so it finishes before the usage limit

There is no usage signal readable from inside the session. The limit kills every running agent at once; observed windows died after roughly 20–30M Fable-weighted tokens. A round of **4 Fable reviews + 2–3 Sonnet writes** (about 17M weighted) finished with room to spare in round 7. Run that as two workflows (one verify, one write-only), both at once, then stop and wait for the user before the next round. Do not run more than three workflows at a time. Commit after each workflow returns; an autosave loop (`/tmp/claude-0/autosave.sh`, commits every 10 min) helps but does not survive a worker restart, so restart it with `setsid nohup ... & disown`.

## Build order (prerequisites first)

1. Verify PY-1–4 (Fable). Write PY-5–8 (Sonnet).
2. Verify PY-5–8, ST-5. Write PY-9, PY-10, EM-6.
3. Verify PY-9, PY-10, EM-6. Write EM-7, EM-8, ST-6, ST-7.
3b. Verify those. Write ST-8, ST-9, QM-1.
4. QM-2 through QM-9 (needs EM done; QM-6 reuses EM-3's Legendre machinery).
5. SM-1 through SM-7 (SM-5 needs QM-5).
6. ML-1 through ML-10 (needs ST-1–5).
7. QC-1 through QC-9 (needs QM-4, QM-7).

Writers can cite earlier modules by code rather than re-deriving.

## Known issues

- **Length.** Almost every module runs over the STYLE.md targets (notes 9–15k words, scripts 500–900 lines). Reviewers consistently decline to trim because the material is correct and brief-mandated. Accepted; not worth the tokens to shorten.
- **Equation tags** out of numerical order in a few modules (EM-3, LA-14) because later modules cite them. Cosmetic; leave them.
- **Conventions** are fixed in `STYLE.md` (complex inner product conjugate-linear in the first slot; qubit 0 leftmost; QFT e^{+2πijk/N}; Osgood Fourier convention in F3, symmetric 1/√(2π) in physics; SI with ℏ and k_B explicit). ST-5 uses Bayes factor BF10 while ST-2 used B01; both state it in their conventions box.

## Final pass, once all 80 modules are done

1. `bash tools/run_all_projects.sh` — every module must print PASS.
2. `python3 tools/py_to_ipynb.py` — builds a `project.ipynb` beside every `project.py`.
3. Update the status column in `README.md`; replace this file with a one-line "complete" note.
4. Commit and push.
