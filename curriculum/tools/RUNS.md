# Build workflow runs (for resuming after interruptions)

Script: `curriculum/tools/build_modules_workflow.js`. Resume a run with
`Workflow({scriptPath: <script>, resumeFromRunId: <id>, args: <same args>})`; completed agents replay from cache.

| Run ID | Modules |
|---|---|
| wf_6c57abdb-8b7 | verify+fix: LA-1..LA-12, LA-14, F1 |
| wf_08ca99dc-1c9 | write: F4, F5; verify: F2, F3, F6 |
| wf_fda2e292-c61 | write: LA-13, CM-1, CM-2 |
| wf_476c56eb-6ca | write: CM-3, CM-4, CM-5 |
| wf_5f4c5a5d-78d | write: CM-6, CM-7, CM-8 |
| wf_6814d784-9c6 | write: EM-1..EM-4 |
| wf_9264db32-034 | write: QM-1..QM-4 |
| wf_9f925196-492 | write: SM-1..SM-4 |

## Second round (after the 2026-10-14 04:20 UTC limit reset; combined review+fix mode)

| Run ID | Modules |
|---|---|
| wf_6c57abdb-8b7 | (resumed) verify+fix: LA-1..LA-12, LA-14, F1 |
| wf_6a49dba9-e21 | combined: write F5, LA-13, CM-2, CM-3, CM-5 |
| wf_336d9009-da2 | combined: write CM-6, CM-7, CM-8; review F4, CM-1, CM-4 |
| wf_7b8457b9-8e2 | combined: write EM-1..EM-4 |

## Third round (2026-10-14 14:08 UTC)

| Run ID | Modules |
|---|---|
| wf_6c57abdb-8b7 | (resumed) fixes for LA-5..LA-12, LA-14, F1 |
| wf_6a49dba9-e21 | (resumed, extended) review CM-5; review F2, F3, F6 |
| wf_c814a063-62a | combined: write EM-1..EM-4 |
| wf_84c5b629-4b6 | combined: write ST-1..ST-4 |

## Fourth round (2026-10-15, resumed after pause; combined review+fix mode)

| Run ID | Modules |
|---|---|
| wf_ef976c8b-db7 | verify: F1, F2, F3, F6 |
| wf_707c8802-f23 | verify: LA-5..LA-9 |
| wf_9d716a9f-8af | verify: LA-10, LA-11, LA-12, LA-14, CM-5 |
| wf_41adb0c0-23e | combined: write EM-3..EM-6 |

## Fifth round (2026-10-15 02:25 UTC, after limit reset; reviewers now write tools/reviews/<code>.md)

| Run ID | Modules |
|---|---|
| wf_8522f993-2ee | verify: F3, F6, LA-7, LA-8 |
| wf_29450baf-c1b | verify: LA-9, LA-12, LA-14, CM-5 |
| wf_796f948f-f96 | combined: write EM-3..EM-6 |

## Sixth round (2026-10-16 01:25 UTC, after limit reset)

| Run ID | Modules |
|---|---|
| wf_546d5cfa-8b9 | verify: LA-7, LA-8, LA-14, CM-5 (resumes saved review notes) |
| wf_bd87d67d-5f5 | verify: EM-1, EM-2, EM-3, ST-1, ST-2 |
| wf_e0130d6e-53e | combined: write EM-4..EM-7 (stopped after ~10 min; replaced by Sonnet writers below) |
| wf_bee1ea55-0fb | combined, writerModel=sonnet: write EM-4..EM-8, Fable review |
| wf_f541add1-916 | write-only, writerModel=sonnet: ST-3..ST-6 (review deferred) |

## Seventh round (2026-10-20 02:55 UTC) — sized to finish below the limit (~4 Fable reviews + 3 Sonnet writes)

Cost calibration from transcripts (cache-read tokens): Fable review 2–4M; a from-scratch write 10–17M on either model;
a window dies after ~20–30M Fable tokens, but round 6 survived 20M Fable + 40M Sonnet, so Sonnet counts at a fraction.

| Run ID | Modules |
|---|---|
| wf_29176906-6af | verify (Fable): CM-5, ST-1, ST-2, EM-4 |
| wf_c8e3bc61-b71 | write-only (Fable) EM-5, ST-3 — stopped after 1 min, replaced by Sonnet writers |
| wf_5cea2e5a-bd3 | write-only (Sonnet): EM-5, ST-3, ST-5 |

## Eighth round (2026-10-23 23:36 UTC) — Python course starts; sized to finish below the limit

| Run ID | Modules |
|---|---|
| wf_1b31082c-bf9 | write-only (Sonnet): PY-1, PY-2, PY-3, PY-4 |
| wf_a6929edc-579 | verify (Fable): ST-3, ST-4, EM-5 |
| wf_88f03f6d-842 | verify (Fable): ST-3 (relaunched after a request timeout) |

## Ninth round (2026-10-25 04:17 UTC) — Python course continues; one round per "continue" from the user

| Run ID | Modules |
|---|---|
| wf_b2a58208-cef | verify (Fable): PY-1, PY-2, PY-3, PY-4 |
| wf_140648ac-5f2 | write-only (Sonnet): PY-5, PY-6, PY-7, PY-8 |
