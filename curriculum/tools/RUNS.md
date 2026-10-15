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
